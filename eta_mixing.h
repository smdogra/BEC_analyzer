
void MixEvents_eta(int ntrkoff_min, int ntrkoff_max, std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vec, std::vector<std::pair<Double_t, std::vector<int>> > ev_GoodTrackCharge_etaMixWeight_vec, std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vec){

///sort vectors of maps for each ntrkoffline range
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > aux_ev_GoodTrackFourVector_etaMixWeight_vec;
std::vector<std::pair<Double_t, std::vector<int>> > aux_ev_GoodTrackCharge_etaMixWeight_vec;
std::vector<std::pair<Double_t,int> > aux_ev_ntrkoff_etaMixWeight_vec;


int aux_n_evts = ev_GoodTrackFourVector_etaMixWeight_vec.size();

for(int ievt=0; ievt<aux_n_evts; ievt++){

   if((ev_ntrkoff_etaMixWeight_vec[ievt]).second<ntrkoff_min || ntrkoff_max<(ev_ntrkoff_etaMixWeight_vec[ievt]).second)continue; //only events in a given range

   aux_ev_GoodTrackFourVector_etaMixWeight_vec.push_back(ev_GoodTrackFourVector_etaMixWeight_vec[ievt]);
   aux_ev_GoodTrackCharge_etaMixWeight_vec.push_back(ev_GoodTrackCharge_etaMixWeight_vec[ievt]);
   aux_ev_ntrkoff_etaMixWeight_vec.push_back(ev_ntrkoff_etaMixWeight_vec[ievt]);

}

std::sort( aux_ev_GoodTrackFourVector_etaMixWeight_vec.begin(), aux_ev_GoodTrackFourVector_etaMixWeight_vec.end(), etaMixSort );
std::sort( aux_ev_GoodTrackCharge_etaMixWeight_vec.begin(), aux_ev_GoodTrackCharge_etaMixWeight_vec.end(), etaMixSort_ch );
std::sort( aux_ev_ntrkoff_etaMixWeight_vec.begin(), aux_ev_ntrkoff_etaMixWeight_vec.end(), etaMixSort_ntrkoff );


int aux_n_evts_inNtrkoffRange = aux_ev_GoodTrackFourVector_etaMixWeight_vec.size();

//auxiliar vectors for mixing
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec;
int nMixmult_nevt1;
std::vector<int> ev_GoodTrackChargeTemp_nevt1_vec;

std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec;
int nMixmult_nevt_assoc;
std::vector<int> ev_GoodTrackChargeTemp_nevt_assoc_vec;

//loop in the vector of maps for mixing - only in a given range of ntrkoffline
for(int nevt=0; nevt+1<aux_n_evts_inNtrkoffRange; nevt+=2) {


   ev_GoodTrackFourVectorTemp_nevt1_vec= (aux_ev_GoodTrackFourVector_etaMixWeight_vec[nevt]).second;
   nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   ev_GoodTrackChargeTemp_nevt1_vec = (aux_ev_GoodTrackCharge_etaMixWeight_vec[nevt]).second;

   ev_GoodTrackFourVectorTemp_nevt_assoc_vec= (aux_ev_GoodTrackFourVector_etaMixWeight_vec[nevt+1]).second;
   nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
   ev_GoodTrackChargeTemp_nevt_assoc_vec= (aux_ev_GoodTrackCharge_etaMixWeight_vec[nevt+1]).second;


   for(int imix=0; imix<nMixmult_nevt1; imix++){
      for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
         if(splitcomb(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix])){continue;}
         Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         Double_t qmixlong = GetQlong(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         Double_t qmixlongLCMS = GetQlongLCMS(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         Double_t qmixout = GetQout(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         Double_t qmixside = GetQside(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);

         TLorentzVector psum2 = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
         Double_t kt=(psum2.Pt())/2.;

         Double_t xmix3D[3]={qmix,kt,(Double_t)(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second};
         Double_t xmix5D_LCMS[5]={qmixlongLCMS,qmixout,qmixside,kt,(Double_t)(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second};

         //Double_t aux_tk1_corr = (Double_t)getTrkCorrWeight(ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Pt(),ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Eta());
         //Double_t aux_tk2_corr = (Double_t)getTrkCorrWeight(ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Pt(),ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Eta());
         //Double_t aux_tk12_corr= aux_tk1_corr*aux_tk2_corr;
	 Double_t aux_tk12_corr=1.;
         if(ev_GoodTrackChargeTemp_nevt1_vec[imix]*ev_GoodTrackChargeTemp_nevt_assoc_vec[iimix]>0){

            hist_mixsch_KtVsNch_ntrkoff->Fill(kt,(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second,aux_tk12_corr);

            hist_qinv_mixschVsKtVsNch_ntrkoff->Fill(xmix3D,aux_tk12_corr);
            hist_qLLCMSqOqSmixschVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS,aux_tk12_corr);

         }else{

            hist_qinv_mixdchVsKtVsNch_ntrkoff->Fill(xmix3D,aux_tk12_corr);
            hist_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS,aux_tk12_corr);

         }
      }
   }



}//end first for loop

//clear vectors for next ntrkoff range
aux_ev_GoodTrackFourVector_etaMixWeight_vec.clear();
aux_ev_GoodTrackCharge_etaMixWeight_vec.clear();
aux_ev_ntrkoff_etaMixWeight_vec.clear();


}//end MixEvents_eta function

void call_mix_eta(int ntrkoff_min, int ntrkoff_max, std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vec, std::vector<std::pair<Double_t, std::vector<int>> > ev_GoodTrackCharge_etaMixWeight_vec, std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vec){

MixEvents_eta(ntrkoff_min, ntrkoff_max, ev_GoodTrackFourVector_etaMixWeight_vec, ev_GoodTrackCharge_etaMixWeight_vec, ev_ntrkoff_etaMixWeight_vec);

}
