#include "call_libraries.h"  // important header file

// This is added for particle identification (GEN only)

// First, particle ID index
enum particleid{
      Pion,
      Kaon,
      Proton,
      K0s,
      Lambda,
      Xi,
      Omega,
      D0
};

// Particle ID number from PDG: https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
int pid[8] = {
      211,    
      321,
      2212,  
      310,  
      3122,    
      3312,
      3334,  
      421   
};

// Particle ID name
TString pid_str[8] = {
      "Pion",
      "Kaon",
      "Proton",
      "K0s",
      "Lambda",
      "Xi",
      "Omega",
      "D0"  
};