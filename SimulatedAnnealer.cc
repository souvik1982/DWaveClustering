/* -----------------------------------------------------------------------
Souvik's Simulated Annealer
(C) Purdue University, August 2019
Author: Souvik Das, souvik@purdue.edu
Description: A simulated annealer that works on bit strings and a QUBO.
Created to benchmark D-Wave's quantum annealer.
------------------------------------------------------------------------ */

#include <math.h>
#include <map>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>

typedef std::map<std::pair<unsigned int, unsigned int>, double> Qubo;
typedef std::vector<int> BitState;
typedef std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, double> > > BitCouplings;

bool DEBUG = false;

// Algorithm parameters
double hotFlipProbability = 0.5;
double coldFlipProbability = 0.01;
double nSweeps = 100;
unsigned int nReads = 1000;
double beta_init, beta_final;

// Forward declared functions
void setTemperatureRange(double hotFlipProbability, double coldFlipProbability, double *beta_init, double *beta_final);
clock_t simulatedAnneal(Qubo *qubo, double beta_init, double beta_final, double nSweeps, BitState *finalState, double *finalEnergy);
void printState(BitState *state);
Qubo* readQUBO(std::string filename, unsigned int *nBits);
double energyEvaluate(Qubo *qubo, BitState *state);
BitCouplings* generateBitCouplings(Qubo *qubo);
double energyDifferenceEvaluate(BitCouplings *bitCouplings, BitState *state, unsigned int positionToFlip);

BitCouplings *bitCouplings;

int main()
{
  // Read parameters from the local configFile.txt
  std::ifstream configFile("configFile_kMeansHard.txt");
  std::string line;
  while (getline(configFile, line))
  {
    unsigned int spacePosition = line.find(" ");
    std::string firstWord = line.substr(0, spacePosition);
    std::string secondWord = line.substr(spacePosition);
    if (firstWord == "nReads") nReads = atoi(secondWord.c_str());
    if (firstWord == "nSweeps") nSweeps = atof(secondWord.c_str());
  }
  std::cout<<"=== Welcome to Souvik's Simulated Annealer for QUBO problems === "<<std::endl;
  std::cout<<" Algorithm parameters used"<<std::endl;
  std::cout<<" nReads = "<<nReads<<std::endl;
  std::cout<<" nSweeps = "<<nSweeps<<std::endl;
  std::cout<<" hotFlipProbability = "<<hotFlipProbability<<std::endl;
  std::cout<<" coldFlipProbability = "<<coldFlipProbability<<std::endl;

  // Read the QUBO from a JSON file
  unsigned int nBits;
  Qubo *qubo = readQUBO("qubo.json", &nBits);

  // Initialize BitCouplings object from QUBO to speed up lookup
  bitCouplings = generateBitCouplings(qubo);

  // Find minimum and maximum energy differences to set range of beta
  setTemperatureRange(0.5, 0.01, &beta_init, &beta_final);
  beta_init = 0.1;
  beta_final = 10;
  std::cout<<" beta_init = "<<beta_init<<std::endl;
  std::cout<<" beta_final = "<<beta_final<<std::endl;

  // Ready response output file
  std::ofstream responseFile("serializableResponse.json");
  responseFile<<"[";

  // Iterate over nReads, and solve the problem each time with nSweeps
  BitState *finalState = new BitState(nBits);
  double finalEnergy;
  double time_total=0, timeSq_total=0;
  for (unsigned int i_read=0; i_read < nReads; ++i_read)
  {
    clock_t time = simulatedAnneal(qubo, beta_init, beta_final, nSweeps, finalState, &finalEnergy);
    time_total += time;
    timeSq_total += time*time;
    // write out final state, energy and time into JSON file
    // std::cout<<"read #"<<i_read<<", finalEnergy = "<<finalEnergy<<", time = "<<time<<" microseconds, finalState = ";
    // printState(finalState);
    // Write out energy and state to JSON for interpretation like D-Wave's response
    responseFile<<"["<<finalEnergy<<", [";
    for (unsigned int i=0; i<finalState->size()-1; ++i) responseFile<<finalState->at(i)<<", ";
    responseFile<<finalState->at(finalState->size()-1)<<"]]";
    if (i_read < nReads - 1) responseFile<<", ";
  }
  responseFile<<"]";
  responseFile.close();
  double time_mean = time_total/nReads;
  double time_stddev = sqrt(timeSq_total/nReads - (time_total*time_total)/float(nReads*nReads));
  std::cout<<"Mean time (us) = "<<time_mean<<" +/- "<<time_stddev/sqrt(nReads)<<std::endl;
  std::cout<<"Stddev time (us) = "<<time_stddev<<" +/- "<<time_stddev/sqrt(nReads)<<std::endl;
  return 0;
}

// --- Core algorithm of simulated annealing ---
clock_t simulatedAnneal(Qubo *qubo, double beta_init, double beta_final, double nSweeps, BitState *finalState, double *finalEnergy)
{
  double step = (beta_final - beta_init) / nSweeps;

  std::clock_t time_start = std::clock();

  // Set nBits state to a random vector
  unsigned int nBits = finalState->size();
  for (unsigned int i = 0; i < nBits; ++i)
    finalState->at(i) = rand() % 2;

  // Simulated annealing iterations from beta = beta_init to beta_final in steps of (beta_stop - beta_init)/nSweeps
  for (double beta = beta_init; beta < beta_final; beta+=step)
  {
    // Flip all bits, one by one
    for (unsigned int positionToFlip = 0; positionToFlip < nBits; ++positionToFlip)
    {
      double energyDifference = energyDifferenceEvaluate(bitCouplings, finalState, positionToFlip);
      if (energyDifference < 0)
        (*finalState)[positionToFlip] = 1 - (*finalState)[positionToFlip];
      else if ((float)rand()/RAND_MAX < exp(-beta*energyDifference))
        (*finalState)[positionToFlip] = 1 - (*finalState)[positionToFlip];
    }
    // std::cout<<"state = "; printState(finalState);
  }
  std::clock_t time_stop = std::clock();
  *finalEnergy = energyEvaluate(qubo, finalState);
  return (time_stop - time_start);
}

void printState(BitState *state)
{
  for (unsigned int i=0; i<state->size(); ++i) std::cout<<state->at(i);
  std::cout<<std::endl;
}

double energyEvaluate(Qubo *qubo, BitState *state)
{
  double energy = 0;
  for (Qubo::iterator i_qubo=qubo->begin(); i_qubo!=qubo->end(); ++i_qubo)
  {
    unsigned int i_state = i_qubo->first.first;
    unsigned int j_state = i_qubo->first.second;
    energy += (i_qubo->second) * state->at(i_state) * state->at(j_state);
  }
  return energy;
}

double energyDifferenceEvaluate(BitCouplings *bitCouplings, BitState *state, unsigned int positionToFlip)
{
  std::vector<std::pair<unsigned int, double> > *v_couplings = &(*bitCouplings)[positionToFlip];
  double energyDiff = 0;
  for (unsigned int i=0; i<v_couplings->size(); ++i)
  {
    energyDiff += state->at((*v_couplings)[i].first) * (*v_couplings)[i].second;
  }
  energyDiff *= (1 - 2*(*state)[positionToFlip]);
  return energyDiff;
}

BitCouplings* generateBitCouplings(Qubo *qubo)
{
  BitCouplings *bitCouplings = new BitCouplings();
  for (Qubo::iterator i_qubo=qubo->begin(); i_qubo!=qubo->end(); ++i_qubo)
  {
    unsigned int i_state = i_qubo->first.first;
    unsigned int j_state = i_qubo->first.second;
    double coupling = i_qubo->second;

    if (bitCouplings->find(i_state) == bitCouplings->end()) // Did not find a record for this bit
    {
      std::vector<std::pair<unsigned int, double> > v_coupling;
      v_coupling.push_back(std::make_pair(j_state, coupling));
      (*bitCouplings)[i_state]=v_coupling;
    }
    else
    {
      std::pair<unsigned int, double> p_coupling = std::make_pair(j_state, coupling);
      (*bitCouplings)[i_state].push_back(p_coupling);
    }

    if (j_state != i_state)
    {
      if (bitCouplings->find(j_state) == bitCouplings->end()) // Did not find a record for this bit
      {
        std::vector<std::pair<unsigned int, double> > v_coupling;
        v_coupling.push_back(std::make_pair(i_state, coupling));
        (*bitCouplings)[j_state]=v_coupling;
      }
      else
      {
        std::pair<unsigned int, double> p_coupling = std::make_pair(i_state, coupling);
        (*bitCouplings)[j_state].push_back(p_coupling);
      }
    }
  }
  return bitCouplings;
}

// Set temperature range based on hotFlipProbability and coldFlipProbability
void setTemperatureRange(double hotFlipProbability, double coldFlipProbability, double *beta_init, double *beta_final)
{
  double minCoupling = 999., maxCouplings = 0.;
  for (BitCouplings::iterator i_bitCouplings = bitCouplings->begin(); i_bitCouplings != bitCouplings->end(); ++i_bitCouplings)
  {
    std::vector<std::pair<unsigned int, double> > *v_couplings = &(i_bitCouplings->second);
    double sumCouplings = 0.;
    for (unsigned int i=0; i<v_couplings->size(); ++i)
    {
      double coupling = v_couplings->at(i).second;
      if (fabs(coupling) >0.1 && fabs(coupling) < minCoupling) minCoupling = fabs(coupling);
      sumCouplings += coupling;
    }
    if (fabs(sumCouplings) > maxCouplings) maxCouplings = fabs(sumCouplings);
  }

  // std::cout<<"minCoupling = "<<minCoupling<<std::endl;
  // std::cout<<"maxCouplings = "<<maxCouplings<<std::endl;

  *beta_init = -log(hotFlipProbability)/maxCouplings;
  *beta_final = -log(coldFlipProbability)/minCoupling;
}

// Read the QUBO in JSON format
Qubo* readQUBO(std::string filename, unsigned int *nBits)
{
  Qubo *qubo = new Qubo();
  std::ifstream quboFile(filename.c_str());
  std::string line;
  getline(quboFile, line);

  *nBits = 0;
  unsigned int i=3;
  int state = 0; // 0 = search for first bit, 1 = search for second bit, 2 = search for coupling
  unsigned int firstBit, secondBit;
  double coupling;
  while (i<line.length())
  {
    if (state == 0)
    {
      unsigned int j = line.find(',', i);
      firstBit = atoi(line.substr(i, j-i).c_str());
      if (DEBUG) std::cout<<"firstBit = "<<firstBit<<std::endl;
      if (firstBit > *nBits) *nBits = firstBit;
      i = j+2;
      state = 1;
    }
    if (state == 1)
    {
      unsigned int j = line.find(']', i);
      secondBit = atoi(line.substr(i, j-i).c_str());
      if (DEBUG) std::cout<<"secondBit = "<<secondBit<<std::endl;
      if (secondBit > *nBits) *nBits = secondBit;
      i = j+2;
      state = 2;
    }
    if (state == 2)
    {
      unsigned int j = line.find(']', i);
      coupling = atof(line.substr(i, j-i).c_str());
      if (DEBUG) std::cout<<"coupling = "<<coupling<<std::endl;
      i = j+5;
      state = 0;
      // Fill qubo
      (*qubo)[std::make_pair(firstBit, secondBit)] = coupling;
    }
  }
  ++(*nBits);
  return qubo;
}
