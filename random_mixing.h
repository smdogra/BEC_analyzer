//#include "call_libraries.h" // call libraries from ROOT and C++
//#include "function_definition.h" // call for functions

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



void MixEvents_random(int ntrkoff_int, int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff_vec, std::vector<double> multiplicity_centrality_bins, std::vector<double> vtx_z, double vzcut, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector_vec, std::vector<std::vector<int>> ev_GoodTrackCharge_vec){
 
 
  int aux_n_evts = (int)ev_ntrkoff_vec.size(); // total number of events

  TRandom2 *r = new TRandom2(); // random number producer
  cout << " I am inside of mixing process and the events to be mixed are  :    " <<  aux_n_evts <<endl;
  /*
  // for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
  for(int nevt1=0; nevt1<5; nevt1++) { // just for checking
        
    std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector_vec[nevt1];
    int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
    std::vector<int> ev_GoodTrackChargeTemp_nevt1_vec = ev_GoodTrackCharge_vec[nevt1];
 
    if(ev_ntrkoff_vec[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff_vec[nevt1])continue;

    int takeAssociated = 0;
    //    for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
    for(int nevt_assoc=nevt1+1; nevt_assoc<3; nevt_assoc++) { // just for checking
      if(ev_ntrkoff_vec[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff_vec[nevt_assoc])continue;

      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector_vec[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      cout << "nMixmult_nevt_assoc = " << nMixmult_nevt_assoc << endl;
      std::vector<int> ev_GoodTrackChargeTemp_nevt_assoc_vec=ev_GoodTrackCharge_vec[nevt_assoc];
      cout << "ev_GoodTrackChargeTemp_nevt_assoc_vec" << ev_GoodTrackCharge_vec[nevt1][nevt_assoc] << endl;

      takeAssociated++;

      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.

      //      DeltazVertex_->Fill((ev_vtx_z_vec)[nevt1] - (ev_vtx_z_vec)[nevt_assoc]);
      for(int imix=0; imix<nMixmult_nevt1; imix++){
	for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
	  //              if(splitcomb(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix])){continue;}
	  Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
	  //Double_t qmixlong = GetQlong(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
	  //              Double_t qmixlongLCMS = GetQlongLCMS(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
	  //              Double_t qmixout = GetQout(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
	  //              Double_t qmixside = GetQside(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
        
	  TLorentzVector psum2 = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
	  Double_t kt=(psum2.Pt())/2.;
        
	  //        Double_t xmix5D_LCMS[5]={qmixlongLCMS,qmixout,qmixside,kt,(Double_t)ev_ntrkoff_vec[nevt1]};
	  //        Double_t xmix3D[3]={qmix,kt,(Double_t)ev_ntrkoff_vec[nevt1]};
 
	  //        Double_t aux_tk1_corr = (Double_t)getTrkCorrWeight(ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Pt(),ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Eta());
	  //        Double_t aux_tk2_corr = (Double_t)getTrkCorrWeight(ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Pt(),ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Eta());
	  //        Double_t aux_tk12_corr= aux_tk1_corr*aux_tk2_corr;
	  //aux_tk12_corr=1.; //just to compare with applying corrections 
	  if(ev_GoodTrackChargeTemp_nevt1_vec[imix]*ev_GoodTrackChargeTemp_nevt_assoc_vec[iimix]>0){

	  //              mixsch_KtVsNch_ntrkoff->Fill(kt,ev_ntrkoff_vec[nevt1],aux_tk12_corr);

	    //        hs_qinv_mixschVsKtVsNch_ntrkoff->Fill(xmix3D,aux_tk12_corr);
	    //        hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS,aux_tk12_corr);

	  }else{

	    //        hs_qinv_mixdchVsKtVsNch_ntrkoff->Fill(xmix3D,aux_tk12_corr);
	    //        hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS,aux_tk12_corr);

	  }

	}
      }
    }//end first for loop 
  }//end MixEvents function
  */
  
}

//Function used only to call the mixing in the main code (Arguments the similar as function above)
//void call_mix_random(int nEvt_to_mix, int nmix_int, std::vector<int> ev_ntrkoff, std::vector<double> multiplicity_centrality_bins, std::vector<double> vtx_z, double vzcut){
void call_mix_random(int ntrkoff_int, int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix, std::vector<int> ev_ntrkoff_vec, std::vector<double> multiplicity_centrality_bins, std::vector<double> vtx_z, double vzcut, std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector_vec, std::vector<std::vector<int>> ev_GoodTrackCharge_vec){
  
  //  MixEvents_random(nmix_int, nEvt_to_mix, ev_ntrkoff, multiplicity_centrality_bins, vtx_z, vzcut);
  MixEvents_random(ntrkoff_int, ntrkoff_min, ntrkoff_max, nEvt_to_mix, ev_ntrkoff_vec, multiplicity_centrality_bins, vtx_z, vzcut, ev_GoodTrackFourVector_vec, ev_GoodTrackCharge_vec);

}
