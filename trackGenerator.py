import numpy
import math
import json
import ROOT

#Sampling through the accept-reject algorithm. h is expected to be a normalized histogram
def sampling_acceptReject(h):
  xMinBin = h.FindFirstBinAbove(0)
  xMaxBin = h.FindLastBinAbove(0)
  xSampleBin = 0
  if xMinBin != -1 and xMaxBin != 1:
    while True:
      xSampleBin = numpy.random.random_integers(xMinBin, xMaxBin)
      #print "xSampleBin = "+str(xSampleBin)
      xRandom = numpy.random.uniform(0, 1)
      if xRandom <= h.GetBinContent(xSampleBin):
        break
  return h.GetBinCenter(xSampleBin)

# track z0 resolution shown in https://twiki.cern.ch/twiki/bin/view/CMSPublic/TrackingPOGPerformance2017MC#Expected_resolutions_on_track_pa
# as a function of eta. Coarsely binned.
def trackz0Resolution_eta(eta):
  eta_fabs = math.fabs(eta)
  dz = 0
  if eta_fabs < 4:
    dz = 0.615
  if eta_fabs < 2.5:
    dz = 0.415
  if eta_fabs < 2:
    dz = 0.275
  if eta_fabs < 1.5:
    dz = 0.175
  if eta_fabs < 1:
    dz = 0.125
  if eta_fabs < 0.5:
    dz = 0.1
  return dz * 0.1 # in cm

# Main program

inFile = ROOT.TFile("/Users/souvik/QuantumComputing/QuantumAnnealing/ThrowTracks/vertexInfo.root")
h_genPV_z = inFile.Get("GenPV_Z")
h_genPV_nTracks = inFile.Get("GenPV_NumTracks")
h_track_pT = inFile.Get("ptSIM")
h_track_eta = inFile.Get("etaSIM")

# Limit genPV_z from -20 to 20
for i in range(1, h_genPV_z.FindBin(-15)):
  h_genPV_z.SetBinContent(i, 0)
for i in range(h_genPV_z.GetNbinsX(), h_genPV_z.FindBin(15)):
  h_genPV_z.SetBinContent(i, 0)
h_genPV_z.GetXaxis().SetRangeUser(-15, 15)

# Zero all track_pT histogram entries with pT < 2 GeV
cut_track_pT = 2 #GeV
for i in range(1, h_track_pT.FindBin(cut_track_pT)):
  h_track_pT.SetBinContent(i, 0)
h_track_pT.GetXaxis().SetRangeUser(2, 60)

h_genPV_z.Scale(1./h_genPV_z.GetSumOfWeights())
h_genPV_nTracks.Scale(1./h_genPV_nTracks.GetSumOfWeights())
h_track_pT.Scale(1./h_track_pT.GetSumOfWeights())
h_track_eta.Scale(1./h_track_eta.GetSumOfWeights())
'''
ROOT.gStyle.SetOptStat(0000)
# Draw histograms we're throwing from, for publication
c_genPV_z = ROOT.TCanvas("c_genPV_z", "c_genPV_z", 700, 700)
h_genPV_z.SetTitle("; Simulated vertex z (cm); Entries")
h_genPV_z.SetLineWidth(2); h_genPV_z.SetLineColor(ROOT.kBlack)
h_genPV_z.Fit("gaus")
h_genPV_z.Draw("")
c_genPV_z.SaveAs("c_genPV_z.png")

c_track_pT = ROOT.TCanvas("c_track_pT", "c_track_pT", 700, 700)
c_track_pT.SetLogy()
h_track_pT.SetLineWidth(2); h_track_pT.SetLineColor(ROOT.kBlack)
h_track_pT.SetTitle("; Simulated track p_{T} (GeV); Entries")
h_track_pT.Draw("hist")
c_track_pT.SaveAs("c_track_pT.png")

c_track_eta = ROOT.TCanvas("c_track_eta", "c_track_eta", 700, 700)
h_track_eta.SetMinimum(0)
h_track_eta.SetTitle("; Simulated track #eta; Entries")
h_track_eta.SetLineWidth(2); h_track_eta.SetLineColor(ROOT.kBlack)
h_track_eta.Draw("hist")
c_track_eta.SaveAs("c_track_eta.png")
'''
nVertices = 4
nTracksMax = 5
d_vertextracks = []
v_vertices = []
for i in range(0, nVertices):

  # Sample position of vertex till a new vertex position is found
  #zVertex = sampling_acceptReject(h_genPV_z)
  zVertex = numpy.random.normal(0, 3.50)
  zVertex_normalized = (zVertex + 15.)/30.

  print "zVertex = "+str(zVertex)+", zVertex_normalized = "+str(zVertex_normalized)
  ''' while True:
    nTracks = int(sampling_acceptReject(h_genPV_nTracks))
    if nTracks == nTracksMax:
      break
  '''
  nTracks = nTracksMax
  print "nTracks = "+str(nTracks)

  v_tracks = []
  for j in range(0, nTracks):
    track_pT = sampling_acceptReject(h_track_pT)
    track_eta = sampling_acceptReject(h_track_eta)
    track_dz0 = trackz0Resolution_eta(track_eta)
    track_dz0_normalized = track_dz0/30.
    track_z0 = numpy.random.normal(zVertex_normalized, track_dz0_normalized)
    print "Track "+str(j)+", pT = "+str(track_pT)+" GeV, eta = "+str(track_eta)+", has track resolution = "+str(track_dz0)+" cm: Assigned position = "+str(track_z0)
    # fill track list
    v_tracks.append((track_z0, track_dz0_normalized))

  d_vertextracks.append((zVertex_normalized, v_tracks))

# Serialize d_vertextracks into a JSON file
outputFile = open("serializedEvents.json", "w")
json.dump(d_vertextracks, outputFile)



''''
# Check sampling through plotting
c_genPV_z = ROOT.TCanvas("c_genPV_z", "c_genPV_z", 700, 700)
h_genPV_z.Draw()
c_genPV_z.SaveAs("c_genPV_z.png")

h_sampled = ROOT.TH1F("h_sampled", "h_sampled", 40, -20, 20)
for i in range(0, 1000):
  x = sampling_acceptReject(h_genPV_z)
  if x != -1:
    h_sampled.Fill(x)

c_sampled = ROOT.TCanvas("c_sampled", "c_sampled", 700, 700)
h_sampled.Draw()
c_sampled.SaveAs("c_sampled.png")
'''
