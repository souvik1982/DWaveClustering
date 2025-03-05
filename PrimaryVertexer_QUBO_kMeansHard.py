#import dimod
#from minorminer import find_embedding
#from dwave.system.composites import EmbeddingComposite
#from dwave.system.samplers import DWaveSampler
import numpy
import math
import json
import ROOT
import time

# Quantum Annealing algorithm for Primary Vertexing
# using a modified k-means algorithm
# Updated 11/11/18
# Based on code by Joel Gottleib @ D-Wave
# and calculations of Andrew J Wildridge, Sachin Vaidya and Souvik Das @ Purdue University

# Read the configuration inputFile
# and override the default parameters
nVertices = 2
runMode = "QPU" # CPU (simulated annealing), QPU (quantum adiabatic annealing), or EXACT (all of solution space)
nReads = 100
knob_problemScale = 1.5 # with Lambda scale = 2.0
knob_k = -1 # for g(x) = x
knob_dropout = 0. # Only scaled D_ij > knob_dropout are considered
nSweeps = 1960
configFile = open("../configFile_kMeansHard.txt")
lines = configFile.readlines()
for line in lines:
  words = line.split()
  if words[0] == "nVertices":
    nVertices = int(words[1])
  if words[0] == "runMode":
    runMode = words[1]
  if words[0] == "nReads":
    nReads = int(words[1])
  if words[0] == "knob_problemScale":
    knob_problemScale = float(words[1])
  if words[0] == "knob_k":
    knob_k = float(words[1])
  if words[0] == "knob_dropout":
    knob_dropout = float(words[1])
  if words[0] == "nSweeps":
    nSweeps = float(words[1])
print("Running Configuration")
print(" nVertices = "+str(nVertices))
print(" runMode = "+str(runMode))
print(" nReads = "+str(nReads))
print(" knob_problemScale = "+str(knob_problemScale))
print(" knob_k = "+str(knob_k))
print(" knob_dropout = "+str(knob_dropout))
print(" nSweeps = "+str(nSweeps))

# Read the file containing track zT_i
'''
inputFile = open("inputTracks.txt")
lines = inputFile.readlines()
zT_i = []
for line in lines:
  words = line.split()
  zT_i.append(float(words[2]))
print "Generated tracks at positions = " + str(zT_i)
nTracks = len(zT_i)
'''

# Read the JSON file containing vertex zV_k and zT_i
d_vertextracks = []
inputFile = open("serializedEvents.json")
d_vertextracks = json.load(inputFile)
zT_i = []
zT_unc_i = []
for vertexTracks in d_vertextracks:
  tracks = vertexTracks[1]
  for i in tracks:
    zT_i.append(i[0])
    zT_unc_i.append(i[1])
    print(str(i[0])+" "+str(i[1]))
nTracks = len(zT_i)

def qubit_pik(i, k):
  return nVertices * i + k

# Function to add to QUBO map
def addToQUBO(Q, i, j, value):
  if (i > j):
    Q[(i, j)] = Q.get((i, j), 0) + value
  else:
    Q[(j, i)] = Q.get((j, i), 0) + value

def serializeQUBO(Q, filename):
  serializeQ = []
  for i in Q:
    serializeQ.append((i, Q[i]))
  outputFile = open(filename, "w")
  json.dump(serializeQ, outputFile)

def exponentialSqueeze(x, k):
  if (k>0) :
    # return knob_problemScale*(1.-math.exp(-k * x / knob_problemScale))
    return (1.-math.exp(-k * x))
  else:
    return x

print("Generating QUBO from problem Hamiltonian")

# Initialization of QUBO matrix
Q = {}
# lambert = (nTracks - nVertices)*.01
lambert = 2.0
offset_energy = lambert * nTracks

# raw distances vectors
v_Dij = []
v_gDij = []

# The k-means Hard Hamiltonian

maxDij = 0.
for i in range(nTracks):
  for j in range(nTracks):
    # Dij = (zT_i[i] - zT_i[j])**2 / (zT_unc_i[i]**2 + zT_unc_i[j]**2)
    Dij = math.fabs(zT_i[i] - zT_i[j]) / math.sqrt(zT_unc_i[i]**2. + zT_unc_i[j]**2.)
    if Dij > maxDij :
      maxDij = Dij
print("maxDij = "+str(maxDij))
scaleFactorCoupling = knob_problemScale/maxDij

for k in range(nVertices):
  for i in range(nTracks):
    for j in range(nTracks):
      q_pik = qubit_pik(i, k)
      q_pjk = qubit_pik(j, k)
      # x = scaleFactorCoupling * (zT_i[i] - zT_i[j])**2 / (zT_unc_i[i]**2 + zT_unc_i[j]**2)
      x = scaleFactorCoupling * math.fabs(zT_i[i] - zT_i[j]) / math.sqrt(zT_unc_i[i]**2. + zT_unc_i[j]**2.)
      if x > knob_dropout :
        v_Dij.append(x)
        term0 = exponentialSqueeze(x, knob_k)
        v_gDij.append(term0)
        addToQUBO(Q, q_pik, q_pjk, term0)

for i in range(nTracks):
  for k in range(nVertices):
    q_pik = qubit_pik(i, k)
    term1 = -2. * lambert
    addToQUBO(Q, q_pik, q_pik, term1)
    for l in range(nVertices):
      q_pil = qubit_pik(i, l)
      term2 = lambert
      addToQUBO(Q, q_pik, q_pil, term2)

# Serialize the QUBO
print("Serializing the QUBO")
serializeQUBO(Q, "qubo.json")
print("Serializing the QUBO done")

ROOT.gStyle.SetOptStat(0000)
h_Dij = ROOT.TH1F("h_Dij", "; D(i, j); Entries", len(v_Dij)/10, min(v_Dij), max(v_gDij)); h_Dij.Sumw2()
h_gDij = ROOT.TH1F("h_gDij", "h_gDij", len(v_Dij)/10, min(v_Dij), max(v_gDij)); h_gDij.Sumw2()
h_Dij.SetLineColor(ROOT.kBlue); h_Dij.SetLineWidth(2)
h_gDij.SetLineColor(ROOT.kRed); h_gDij.SetLineWidth(2)
for i in range(0, len(v_Dij)):
  h_Dij.Fill(v_Dij[i])
  h_gDij.Fill(v_gDij[i])
c_Dij = ROOT.TCanvas("c_Dij", "c_Dij", 700, 700)
h_Dij.Draw("hist")
h_Dij.Draw("ep9 same")
h_gDij.Draw("hist same")
h_gDij.Draw("ep same")
leg = ROOT.TLegend(0.5, 0.7, 0.89, 0.89)
leg.SetFillColor(0); leg.SetLineColor(0)
leg.AddEntry(h_Dij, "D(i, j)")
leg.AddEntry(h_gDij, "g(D(i, j), k=5)")
leg.Draw()
c_Dij.SaveAs("c_Dij.png")
c_Dij.SaveAs("c_Dij.root")

g_exponentialSqueeze = ROOT.TGraph()
g_exponentialSqueeze.SetTitle("; D(i, j); g(D(i, j); m = 5)")
g_exponentialSqueeze.GetXaxis().SetRangeUser(0, 1.6)
g_exponentialSqueeze.SetMinimum(0)
for i in range(0, len(v_Dij)):
  g_exponentialSqueeze.SetPoint(i, v_Dij[i]/knob_problemScale, v_gDij[i])
c_exponentialSqueeze = ROOT.TCanvas("c_exponentialSqueeze", "c_exponentialSqueeze", 700, 700)
g_exponentialSqueeze.Draw("A*")
c_exponentialSqueeze.SaveAs("c_exponentialSqueeze.png")
c_exponentialSqueeze.SaveAs("c_exponentialSqueeze.pdf")
quit()

print("Generating QUBO from problem Hamiltonian done.")

'''

# Anneal the QUBO matrix with D-Wave
# Different running pattern depending on:
#   CPU (simulated annealing),
#   QPU (quantum adiabatic annealing), or
#   EXACT (explore all solution space)
if runMode == "QPU":
  print("Processing with D-Wave QPU")
  bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
  #sampler = DWaveSampler()
  sampler = DWaveSampler(solver={'lower_noise': False, 'qpu': True}, token='PURD-87bd7895a858ac9dc635ae9f4d9b9ad6c56602b8')
  embeddingComposite = EmbeddingComposite(sampler)
  #embedding = find_embedding(Q, sampler.edgelist)
  #embeddingComposite = VirtualGraphComposite(sampler, embedding)
  #sampler.properties["annealing_time"] = 10
  #print "EmbeddingComposite properties = "+str(embeddingComposite.properties)
  #print "sampler.properties[annealing_time_range] = "+str(sampler.properties["annealing_time_range"])
  #print "sampler.properties[default_annealing_time] = "+str(sampler.properties["default_annealing_time"])
  # quit()
  response = embeddingComposite.sample(bqm, num_reads=nReads, answer_mode="raw")
  print("response = "+str(response.info["timing"]["qpu_programming_time"])) # Contains timing information
  print("Processing reads done")
elif runMode == "CPU":
  print("Processing with CPU")
  bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
  sampler = dimod.SimulatedAnnealingSampler()
  time_init = time.clock()
  response = sampler.sample(bqm, num_reads=nReads, num_sweeps=int(nSweeps))
  time_stop = time.clock()
  print("Initial time = "+str(time_init))
  print("Final time = "+str(time_stop))
  print("Difference time = "+str(time_stop - time_init))
  # print("response = "+str(response))
  response.info['timing'] = time_stop - time_init
elif runMode == "EXACT":
  print("Processing exact solution")
  bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
  sampler = dimod.ExactSolver()
  response = sampler.sample(bqm)
  print("Processing exact solution done")

# Serialize the response and store to JSON file
# Stored in increasing order of energy
print("Serializing response")
serializableResponse = []
for i_record in response.record:
  energy = i_record.energy
  solution = i_record[0].tolist()
  serializableResponse.append((energy, solution))
outputFile = open("serializableResponse.json", "w")
json.dump(serializableResponse, outputFile)
print("Serializing response done")


for i in range(len(energy_ordered_indices)):
  index = energy_ordered_indices[i]
  energy = response.data_vectors['energy'][index]
  solution = str(response.samples_matrix[index])[1:-1].split()
  if i <= 10:
    print("Energy is " + str(energy) + ", corrected with offset = " + str(energy + offset_energy))
    print("The solution is: " + str(solution))
  serializableResponse.append((energy, solution))
outputFile = open("serializableResponse.json", "w")
json.dump(serializableResponse, outputFile)
print "Serializing response done"
'''
