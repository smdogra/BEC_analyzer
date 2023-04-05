//#include "call_libraries.h" // call libraries from ROOT and C++
//#include "function_definition.h" // call for functions
#include  "histogram_definition_bec.h"

/*
Function for random mixing
--> Arguments
ntrkoff_int: range of multiplicity or centrality between first and mixed events
nEvt_to_mix: number of events to mix
ev_ntrkoff: vector with multiplicity or centrality for each event
multiplicity_centrality_bins: multiplicity or centrality bins defined at input_variables.h
vtx_z: vector with z vertex values for each event
vzcut: range of Vz between first and mixed events 
Jet_Vector: vector of vectors with TVector3 (pT, eta, phi) of jets (basically vectors of jets for each event)
Jet_W_Vector: vector of vectors with weights from jets (basically vectors of weights for each jet)
Track_Vector: vector of vectors with TVector3 (pT, eta, phi) of tracks (basically vectors of tracks for each event)
Track_W_Vector: vector of vectors with weights from track (basically vectors of weights for each track)
histo: multidimentional histograms with correlations (DeltaPhi, DeltaEta, TrackBin, MultorCentBin)
trk_pt_bin: track bins defined at input_variables.h
event_weight: event weight vector
histo_jet: multidimentional histograms with jet quantities (Pt, Eta, Phi, MultorCentBin)
histo_trk: multidimentional histograms with track quantities (Pt, Eta, Phi, MultorCentBin)
double_weight: boolean to apply or not double weighting
*/
//void MixEvents_random(int ntrkoff_int, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> multiplicity_centrality_bins, std::vector<double> vtx_z, double vzcut, std::vector<std::vector<TVector3>> Jet_Vector, std::vector<std::vector<double>> Jet_W_Vector, std::vector<std::vector<TVector3>> Track_Vector, std::vector<std::vector<double>> Track_W_Vector, THnSparse *histo, std::vector<double> trk_pt_bin, std::vector<double> event_weight, THnSparse *histo_jet, THnSparse *histo_trk, bool double_weight){
void MixEvents_random(int ntrkoff_int, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> multiplicity_centrality_bins, std::vector<double> vtx_z, double vzcut){
  int aux_n_evts = (int)ev_ntrkoff.size(); // total number of events
  TRandom2 *r = new TRandom2(); // random number producer
  cout << " I am inside of mixing process and the events to be mixed are  :    " <<  aux_n_evts <<endl;
  /*
  // first loop over all events to define the jet event
  for(int nevt_trg = 0; nevt_trg < aux_n_evts; nevt_trg++){
  
  int n_associated = 0; // counter used to find the number to mix 
  int n_associated_check = 0; // counter assure we find the required number of events to mix
  std::vector<int> eventcheck; // vector to make sure we do not have same combinations
  
  std::vector<TVector3> Jet_nevt_trg_vec = Jet_Vector[nevt_trg]; // jet vector for each event
  std::vector<double> Jet_w_nevt_trg_vec = Jet_W_Vector[nevt_trg]; // jet weight vector for each event
  int nMix_nevt_trg = Jet_nevt_trg_vec.size(); // jet vector size
  
  double weight_trgev = event_weight[nevt_trg]; // event weight vector
  
  // while the number of mixed events is less than the required do this loop
  while (n_associated < nEvt_to_mix){ 
  
  // additional check, if we do not have nEvt_to_mix after check all events, increase the values
  n_associated_check = n_associated_check + 1;
  if (n_associated_check == (aux_n_evts - 1)*10){
  cout << "Number of events for mixing is not enough in range: [0," << n_associated_check << "]"<< endl; 
  cout << "Increasing the ranges: vz in 0.2 cm and multiplicity or centrality in 1" << endl; 
  vzcut = vzcut + 0.2;
  ntrkoff_int = ntrkoff_int+1;
  n_associated_check = 0;
  }
  
  int nevt_assoc = r->Integer(aux_n_evts-1); // random events between 0 and aux_n_evts-1
  if(nevt_trg == nevt_assoc) continue; // make sure that we do not choose same events
  if(fabs(ev_ntrkoff[nevt_trg] - ev_ntrkoff[nevt_assoc]) > ntrkoff_int) continue; // multiplicity or centrality matching
  if(fabs(vtx_z[nevt_trg] - vtx_z[nevt_assoc]) > vzcut) continue; // vz matching
  
  // make sure that we do not duplicate the correlation
  bool isduplicated = false; 
  if(eventcheck.size() > 0){for(int i = 0; i < (int)eventcheck.size(); i++){if(nevt_assoc == eventcheck[i]){isduplicated = true; break;}}}
  if(isduplicated == true) continue;
  
  // get weights
  double weight_assev = event_weight[nevt_assoc];
  double f_weight = 1.0;
  if(!double_weight){f_weight =  weight_trgev;}else{f_weight = weight_trgev * weight_assev;} //weighting 
  
  eventcheck.push_back(nevt_assoc); // save all events that pass the requirements
  
  n_associated = n_associated + 1; // if pass the requirements sum 1
  
  std::vector<TVector3> Track_nevt_ass_vec = Track_Vector[nevt_assoc]; // get the vector of tracks
  std::vector<double> Track_w_nevt_ass_vec = Track_W_Vector[nevt_assoc]; // get the vector of tracks      
  int nMix_nevt_ass = (int)Track_nevt_ass_vec.size(); // get the vector of tracks size
  
  // loop and fill correlation histograms
  for(int imix = 0; imix < nMix_nevt_trg; imix++){
  double jet_weight = Jet_w_nevt_trg_vec[imix];
  for(int iimix = 0; iimix < nMix_nevt_ass; iimix++){
  double trkpt = Track_nevt_ass_vec[iimix].Pt();
  // track efficiency correction for reco
  double trk_weight = Track_w_nevt_ass_vec[iimix];
  // find track and multiplicity bins
  int trkbin = (int) find_my_bin(trk_pt_bin,trkpt);
  int multcentbin = (int) find_my_bin(multiplicity_centrality_bins, (float) ev_ntrkoff[nevt_trg]);
              // Fill correlation histograms
              double x4D_jet[4]={Jet_nevt_trg_vec[imix].Pt(),Jet_nevt_trg_vec[imix].Eta(),Jet_nevt_trg_vec[imix].Phi(),(double)multcentbin}; histo_jet->Fill(x4D_jet,jet_weight*f_weight);
              double x4D_trk[4]={Track_nevt_ass_vec[iimix].Pt(),Track_nevt_ass_vec[iimix].Eta(),Track_nevt_ass_vec[iimix].Phi(),(double)multcentbin}; histo_trk->Fill(x4D_trk,trk_weight*f_weight);
              double del_phi = deltaphi2PC(Jet_nevt_trg_vec[imix].Phi(), Track_nevt_ass_vec[iimix].Phi());
              double del_eta = deltaeta(Jet_nevt_trg_vec[imix].Eta(), Track_nevt_ass_vec[iimix].Eta());
              double x4D[4]={del_phi,del_eta,(double)trkbin,(double)multcentbin}; histo->Fill(x4D,jet_weight*trk_weight*f_weight*trkpt);
            }
         } // end of correlation loop

      } // end of while loop (after finding the number of events required)
      eventcheck.clear();
   } // end of all events loop
   */
}

// Function used only to call the mixing in the main code (Arguments the similar as function above)
void call_mix_random(int nEvt_to_mix, int nmix_int, std::vector<int> ev_ntrkoff, std::vector<double> multiplicity_centrality_bins, std::vector<double> vtx_z, double vzcut){
  
  MixEvents_random(nmix_int, nEvt_to_mix, ev_ntrkoff, multiplicity_centrality_bins, vtx_z, vzcut);
   
}
