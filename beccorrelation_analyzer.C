#include "call_libraries.h"  // call libraries from ROOT and C++
#include "read_tree.h" // read the TChains
#include "vector_definition.h"  // define the vectors for mixing
#include "histogram_definition_bec.h" // define histograms
#include "uiclogo.h" // print UIC jets and start/stop time
#include "function_definition.h"
#include "random_mixing.h" // random mixing function  
#include "eta_mixing.h" // eta mixing function  
/*
Main code to run Jet+Track correlation
Written by Dener Lemos (dener.lemos@cern.ch)
--> Arguments
input_file: text file with a list of root input files: Forest or Skims
ouputfilename: just a text to run on Condor
MCSim: 0 for data and > 0 for MC
pthatmin: pthat min cut for MC only
pthatmax: pthat max cut for MC only
*/

//******************************************************************************************/
//void beccorrelation_analyzer(TString input_file, TString ouputfilename, int MCSim, float pthatmin, float pthatmax){
void beccorrelation_analyzer(TString input_file, TString ouputfilename, int MCSim, bool isEventMix){
  cout<< " BEC code started" << endl;  
  clock_t sec_start, sec_end, sec_start_mix, sec_end_mix;
  sec_start = clock(); // start timing measurement
  
  TDatime* date = new TDatime();
  
  printwelcome(true); // welcome message
  
  print_start(); // start timing print
  
  //-----------------------------------------------------------------------------------------
  
  // var defined for same event
  bool splitcomb(TLorentzVector &vec1,TLorentzVector &vec2);  
  Double_t GetQ( const TLorentzVector &p1, const TLorentzVector &p2 );
  Double_t GetQlong( const TLorentzVector &p1, const TLorentzVector &p2 );
  Double_t GetQlongLCMS( const TLorentzVector &p1, const TLorentzVector &p2 );
  Double_t GetQout( const TLorentzVector &p1, const TLorentzVector &p2 );
  Double_t GetQside( const TLorentzVector &p1, const TLorentzVector &p2 );
  const TLorentzVector InvertPVector( TLorentzVector &vec);
  const Double_t CoulombWpm( const Double_t& q);
  const TLorentzVector InvertXYVector( TLorentzVector &vec);
  const Double_t CoulombW( const Double_t& q);
  Double_t Encode( int w );
  Double_t ComputeEventWeight(std::vector<TLorentzVector> GoodTrackFourVector);
  bool etaMixSort(std::pair<Double_t,std::vector<TLorentzVector>> & a, std::pair<Double_t,std::vector<TLorentzVector>> & b);
  bool etaMixSort_ch(std::pair<Double_t,std::vector<int>> & a, std::pair<Double_t,std::vector<int>> & b);
  bool etaMixSort_ntrkoff(std::pair<Double_t,int> & a, std::pair<Double_t,int> & b);


  // var defined for mixed events

  int mix_procedure_; //which mixing procedure to use. If equal to 1 is random and 2 is eta mix 
  int Multiplicity_; //create histograms in bins of multiplicity according to the trigger.

  //event information                                                                                   
  int evNumber, runNumber, LumiSection;
  double bs_x, bs_y, bs_z;
  int N_vtx;
  double vtx_x, vtx_y, vtx_z, vtx_xError, vtx_yError, vtx_zError;  //primary vertex                     
  std::vector<double> ev_vtx_z_vec;
  std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector_vec;
  std::vector<std::vector<int>> ev_GoodTrackCharge_vec;
  std::vector<int> ev_ntrkoff_vec;
  std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vec;
  std::vector<std::pair<Double_t, std::vector<int>> > ev_GoodTrackCharge_etaMixWeight_vec;
  std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vec;


  //---------------------------------------------------------------------------------------------

  bool is_MC;
  if(MCSim==0){is_MC = false;}else{is_MC = true;}
  
  if(colliding_system!="pPb") do_CM_pPb = false; // Only do center-of-mass for pPb
  
  //print important informations in the output file
  TString data_or_mc;
  if(!is_MC){data_or_mc="Data";}else{data_or_mc="MC";}
  TString simev;
  if(similar_events){simev = "simevs";}else{simev = "";}
  TString ref_sample = "norefsample";
  //  if(do_mixing && !do_rotation){ref_sample = Form("mix%ievsMult%iDVz%.1f%s",N_ev_mix,Mult_or_Cent_range,DVz_range,simev.Data());}else if(!do_mixing && do_rotation){ref_sample = Form("rot%ievs",N_of_rot);}else if(do_mixing && do_rotation){ref_sample = Form("mix%ievsMult%iDVz%.1f%s_rot%ievs",N_ev_mix,Mult_or_Cent_range,DVz_range,simev.Data(),N_of_rot);}
  //  if(!do_pid) particles = "CH";
  // In case of wrong input, printout error message and kill the job
  if(year_of_datataking!=2012 && year_of_datataking!=2016 && year_of_datataking!=2017 && year_of_datataking!=2018){cout << "Data and MC not supported: choose 2012 for pp at 8 TeV, 2016 for pPb at 8.16 TeV, 2017 for pp at 5.02 TeV or XeXe at 5.44 TeV and 2018 for PbPb at 5.02 TeV" << endl; return;}
  if(colliding_system!="pp" && colliding_system!="pPb" && colliding_system!="XeXe" && colliding_system!="PbPb"){cout << "Data and MC not supported: choose pp for proton-proton, pPb for proton-lead, PbPb for lead-lead and XeXe for Xenon-Xenon" << endl; return;}
  if(sNN_energy_GeV!=5020 && sNN_energy_GeV!=5440 && sNN_energy_GeV!=8000 && sNN_energy_GeV!=8160 && sNN_energy_GeV!=13000){cout << "Data and MC not supported: 5020 for pp 2017 or PbPb 2018, 5440 for XeXe, 8000 for pp 2018, 8160 for pPb 2016" << endl; return;}


  //----------------------------------------------------------------------------------------------
  // Track or particle efficiency file
  TFile *fileeff = TFile::Open(Form("aux_files/%s_%i/trk_eff_table/%s",colliding_system.Data(),sNN_energy_GeV,trk_eff_file.Data()));
  cout << endl;
  TH2 *reff2D = nullptr; TH2 *rsec2D = nullptr; TH2 *rfak2D = nullptr; TH2 *rmul2D = nullptr;
  TH3 *reff3D = nullptr; TH3 *rsec3D = nullptr; TH3 *rfak3D = nullptr; TH3 *rmul3D = nullptr;
  fileeff->GetObject("eff", reff2D);	// Absolute Efficiency
  fileeff->GetObject("fake", rfak2D); // Fake Reconstruction Fraction
  fileeff->GetObject("sec", rsec2D);	// Multiple Reconstruction Fraction
  fileeff->GetObject("mult", rmul2D); // Non-Primary Reconstruction Fraction
  vector<TH2*> eff_histos={reff2D, rfak2D, rsec2D, rmul2D};
  vector<TH3*> eff_histos3D={reff3D, rfak3D, rsec3D, rmul3D};

  //-----------------------------------------------------------------------------------------------
  //Print the input in the screen/log 
  //print_input(data_or_mc,fileeff,colliding_system,pthatmin,pthatmax);
  print_input(data_or_mc,fileeff,colliding_system);
  cout << endl;
  
  //------------------------------------------------------------------------------------------------
  // Read the input file(s)
  fstream inputfile;
  inputfile.open(Form("%s",input_file.Data()), ios::in);
  if(!inputfile.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << input_file.Data() << endl;}

  //-------------------------------------------------------------------------------------------------
  // Make a chain and a vector of file names
  std::vector<TString> file_name_vector;
  string file_chain;
  while(getline(inputfile, file_chain)){file_name_vector.push_back(file_chain.c_str());}
  inputfile.close();
  
  //--------------------------------------------------------------------------------------------------
  // Read the trees to be added in the Chain
  TChain *hlt_tree = new TChain("hltanalysis/HltTree");
  TChain *trk_tree = new TChain("ppTrack/trackTree");
  TChain *hea_tree = new TChain("hiEvtAnalyzer/HiTree");
  TChain *gen_tree;
  if(is_MC){gen_tree = new TChain("HiGenParticleAna/hi");}
  TChain *ski_tree = new TChain("skimanalysis/HltTree");
  
  //-----------------------------------------------------------------------------------------------------
  // add all the trees to the chain
  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
    cout << "Adding file " << *listIterator << " to the chains" << endl;
    hlt_tree->Add(*listIterator);
    trk_tree->Add(*listIterator);
    hea_tree->Add(*listIterator);
    ski_tree->Add(*listIterator);
    if(is_MC){gen_tree->Add(*listIterator);}
  }

  //-------------------------------------------------------------------------------------------------------
  // Connect all chains
  hlt_tree->AddFriend(trk_tree);
  hlt_tree->AddFriend(hea_tree);
  hlt_tree->AddFriend(ski_tree);	
  if(is_MC){hlt_tree->AddFriend(gen_tree);}

  //--------------------------------------------------------------------------------------------------------
  // Read the desired branchs in the trees
  read_tree(hlt_tree, is_MC, colliding_system.Data(), sNN_energy_GeV, year_of_datataking, event_filter_str, event_filter_bool); // access the tree informations
  // Use sumw2() to make sure about histogram uncertainties in ROOT
  sw2(); 

  //**********************************************************************************************************/
  int nevents = hlt_tree->GetEntries(); // number of events
  cout << "Total number of events in those files: "<< nevents << endl;
  cout << endl;
  cout << "-------------------------------------------------" << endl;
  
  // Start loop over events
  double nev = (double)nevents;
  nev =1000;
  for(int i = 0; i <nev; i++){
    hlt_tree->GetEntry(i);
    if(i != 0 && (i % 10000) == 0){double alpha = (double)i; cout << " Running -> percentage: " << std::setprecision(3) << ((alpha / nev) * 100) << "%" << endl;}
    Nevents->Fill(0); // filled after each event cut
    // Apply event filters
    //for(int ii = 0; ii < event_filter_bool.size(); ii++) if(event_filter_bool[ii] != 1) continue;
    bool passFilters=true;
    for(int ii = 0; ii < event_filter_bool.size(); ii++){
       if(event_filter_bool[ii] == 0){
         passFilters=false;
         break;
       }	 
    }
    if(!passFilters) continue;
    Nevents->Fill(1);
  

  // Vectors used for objects
    // reco tracks

    std::vector<TVector3> tracks_reco;
    std::vector<int> sube_tracks_reco;
    std::vector<double> track_w_reco;
    // gen tracks

    std::vector<TVector3> tracks_gen;
    std::vector<int> sube_tracks_gen;
    std::vector<double> track_w_gen;

    std::vector<TLorentzVector> GoodTrackFourVector;
    std::vector<int> GoodTrackCharge;
    std::vector<TLorentzVector> GoodTrackFourVector_trkoff;

    //---------------------------------------------------------------------------------------------------
    

    //apply event selections
    //Vz
    if(vertexz <= vz_cut_min || vertexz >= vz_cut_max) continue;
    Nevents->Fill(2);
    //multiplicity or centrality
    int trksize = (int)ntrk;
    int mult;
    
    if(use_centrality){mult = hiBin;}else{mult = get_Ntrkoff(colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkcharge, highpur, trkpterr, trkdcaxy, trkdcaxyerr, trkdcaz, trkdcazerr, trkchi2, trkndof, trknlayer, trknhits, trkalgo, trkmva);}
    if(mult < multiplicity_centrality_bins[0] || mult > multiplicity_centrality_bins[multiplicity_centrality_bins.size()-1])continue; //centrality of multiplicity range	
    int multcentbin = (int) find_my_bin(multiplicity_centrality_bins, (float) mult);
    Nevents->Fill(3);

    // event weight(s), this must be applied in all histograms
    //    double event_weight = get_event_weight(is_MC, use_centrality, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, vertexz, mult, weight, pthat); // get the event weight
    
    double event_weight =1.0;
    // Fill vertex, pthat and multiplicity/centrality histograms
    vzhist->Fill(vertexz);
    vzhist_weighted->Fill(vertexz, event_weight);
    multiplicity->Fill(mult);
    multiplicity_weighted->Fill(mult, event_weight);
    
    // Reconstruction level	(Data and MC)
    // Start loop over reco tracks (trksize is number of reco tracks)
    for (int j = 0; j < trksize; j++){ 
      // Define track/particle kinematics
      float trk_pt = trkpt[j];
      float trk_eta = trketa[j];
      float trk_phi = trkphi[j];
      float trk_charge=  trkcharge[j];

      //for eta-mixing
      if(trk_pt>0.4 && fabs(trk_eta)<2.4 && fabs(trkpterr[j]/trkpt[j])<0.1 && fabs(trkdcaxy[j]/trkdcaxyerr[j])<3 && fabs(trkdcaz[j]/trkdcazerr[j])<3 && highpur[j] ==true){
         TLorentzVector pvector_etamix;
	 pvector_etamix.SetPtEtaPhiM(trk_pt,trk_eta, trk_phi,pi_mass);
         GoodTrackFourVector_trkoff.push_back(pvector_etamix); 	      
      }

      // In pPb case, for the center-of-mass correction if needed
      if(colliding_system=="pPb" && do_CM_pPb){if(is_pgoing){trk_eta = trk_eta - 0.465;}else{trk_eta = -trk_eta - 0.465;}}
      
      // Apply track selection (see read_tree.h to see what each variable means) ___ read this
      if(fabs(trk_eta) > trk_eta_cut) continue;
      if(trkpt[j] <= trk_pt_min_cut) continue;
      if(highpur[j] == false) continue;

      if(fabs(trkpterr[j]/trkpt[j]) >= trk_pt_resolution_cut) continue;
      if(fabs(trkdcaxy[j]/trkdcaxyerr[j]) >= trk_dca_xy_cut) continue;
      if(fabs(trkdcaz[j]/trkdcazerr[j]) >= trk_dca_z_cut) continue;
      if((int) trkpixhits[j]< 1) continue;
      double calomatching = ((pfEcal[j]+pfHcal[j])/cosh(trketa[j]))/trkpt[j];
      if(colliding_system == "PbPb" || colliding_system == "XeXe"){
	if((trkchi2[j]/trkndof[j])/trknlayer[j] >= chi2_ndf_nlayer_cut) continue;
	if(trknhits[j] < nhits) continue;
	if(trkpt[j] > 20.0 && fabs(calomatching) <= calo_matching) continue;
      }
      if(colliding_system=="PbPb" && sNN_energy_GeV==5020 && year_of_datataking==2018){if(trkalgo[j] == 6 && trkmva[j] < 0.98) continue;}
      
      // Track efficiency correction
      double trk_weight = 1.0;
      trk_weight = trk_weight*getTrkCorrWeight(eff_histos[0], eff_histos[1], eff_histos[2], eff_histos[3], trk_pt, trk_eta);

      //cout<< " The weight of the tracks   :   "  <<trk_eta << " " << trk_pt << " " <<  trk_weight << endl;
      
      // Track QA histogram filling
      double x4D_reco_trk[4]={trk_pt,trk_eta,trk_phi,(double) multcentbin}; 
      hist_reco_trk->Fill(x4D_reco_trk);
      hist_reco_trk_corr->Fill(x4D_reco_trk,trk_weight);
      hist_reco_trk_weighted->Fill(x4D_reco_trk,trk_weight*event_weight);
     
      hist_trketa->Fill(trketa[j]); 
      hist_trkpt->Fill(trkpt[j]); 
      hist_trkcharge->Fill(trkcharge[j]); 
      hist_trkpterrbytrkpt->Fill(trkpterr[j]/trkpt[j]);       
      hist_trkdcaxybytrkerrdcaxy->Fill(trkdcaxy[j]/trkdcaxyerr[j]); 
      hist_trkdcazbytrkdcazerr->Fill(trkdcaz[j]/trkdcazerr[j]); 
      //      hist_trkNPixelHit->Fill(trknpixelhit[j]);
      hist_trkchi2bytrkndof->Fill(trkchi2[j]/trkndof[j]); 
      hist_trkchi2bytrkndofbytrknlayer->Fill((trkchi2[j]/trkndof[j])/trknlayer[j]); 
      hist_trknlayer->Fill(trknlayer[j]); 
      hist_trknhits->Fill(trknhits[j]); 
      hist_trkalgo->Fill(trkalgo[j]); 
      hist_trkmva->Fill(trkmva[j]);
     
      //      cout<<  trkchi2[j] <<  " " <<  trkndof[j] << " " <<  trknlayer[j] << " " <<  trknhits[j]<<  endl;
      // Track vector filling


      TLorentzVector pvector;  
      pvector.SetPtEtaPhiM(trk_pt,trk_eta, trk_phi,pi_mass);
      GoodTrackFourVector.push_back(pvector);
      GoodTrackCharge.push_back(trk_charge);

    } // End loop over tracks
   
    //if(Multiplicity_ == 0){if(120<=mult){return;}}                                          
    //if(Multiplicity_ == 0){if(800<mult){return;}} //for extended MinBias                      
    //else if(Multiplicity_ == 1){if(120>mult || 150<=mult){return;}}
    //else if(Multiplicity_ == 2){if(150>mult || 185<=mult){return;}}
    //else if(Multiplicity_ == 3){if(185>mult || 250<=mult){return;}}
    //else if(Multiplicity_ == 4){if(250>mult){return;}}
    //else{}

    //-----------------------------------------------------------------------------------------------------

    if(GoodTrackFourVector.size()<2)continue; //event not used for signal, then do not use for mixing 
    //cout<<" Size of the tracks +++++++    " <<GoodTrackFourVector.size() << " "<< mult <<  endl;
    for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
      //Double_t aux_tk1_corr = 1.0; // No track corrections
      Double_t aux_tk1_corr = (Double_t)getTrkCorrWeight(eff_histos[0], eff_histos[1], eff_histos[2], eff_histos[3],GoodTrackFourVector[itk1].Pt(),GoodTrackFourVector[itk1].Eta());
      Double_t aux_pt = GoodTrackFourVector[itk1].Pt();
      Double_t aux_eta = GoodTrackFourVector[itk1].Eta();
      Double_t aux_phi = GoodTrackFourVector[itk1].Phi();
      if(aux_pt<0.4)cout<< " The weight of the tracks   :   "  <<aux_eta << " " << aux_pt << " " <<  aux_tk1_corr << endl;
      
      hist_ptVsNch_ntrkoff->Fill(aux_pt,mult,   aux_tk1_corr);
      hist_etaVsNch_ntrkoff->Fill(aux_eta,mult, aux_tk1_corr);
      hist_phiVsNch_ntrkoff->Fill(aux_phi,mult, aux_tk1_corr);

      for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
	if(splitcomb(GoodTrackFourVector[itk1],  GoodTrackFourVector[itk2])){continue;}
	
	Float_t d_eta= GoodTrackFourVector[itk1].Eta()-GoodTrackFourVector[itk2].Eta();           
	Float_t d_phi= GoodTrackFourVector[itk1].DeltaPhi(GoodTrackFourVector[itk2]); 
	Double_t costheta= TMath::Abs(GoodTrackFourVector[itk1].Px()*GoodTrackFourVector[itk2].Px() + GoodTrackFourVector[itk1].Py()*GoodTrackFourVector[itk2].Py() + GoodTrackFourVector[itk1].Pz()*GoodTrackFourVector[itk2].Pz())/(GoodTrackFourVector[itk1].P()*GoodTrackFourVector[itk2].P());
	Double_t deltapt1 = TMath::Abs (GoodTrackFourVector[itk1].Pt() - GoodTrackFourVector[itk2].Pt());// added by sunil for for dpt cut                                            
	Double_t q = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);                 //  Qinv                        
	Double_t qo= GetQ(GoodTrackFourVector[itk1],InvertPVector(GoodTrackFourVector[itk2]));   //  Qinv inverted Vector        
	Double_t qs= GetQ(GoodTrackFourVector[itk1],InvertXYVector(GoodTrackFourVector[itk2]));  //  Qinv Rotarted Vector        
       	Double_t qlong =  GetQlong(GoodTrackFourVector[itk1],GoodTrackFourVector[itk2]);      //  3D                          
	Double_t qlongLCMS = GetQlongLCMS(GoodTrackFourVector[itk1],GoodTrackFourVector[itk2]);  //  3D                          
	Double_t qout =  GetQout(GoodTrackFourVector[itk1],GoodTrackFourVector[itk2]);       //  3D                          
	Double_t qside = GetQside(GoodTrackFourVector[itk1],GoodTrackFourVector[itk2]);      //  3D  
	Double_t qlongLCMS_inv = GetQlongLCMS(GoodTrackFourVector[itk1],InvertPVector(GoodTrackFourVector[itk2]));
	Double_t qlongLCMS_rot = GetQlongLCMS(GoodTrackFourVector[itk1],InvertXYVector(GoodTrackFourVector[itk2]));
	Double_t qout_inv = GetQout(GoodTrackFourVector[itk1],InvertPVector(GoodTrackFourVector[itk2]));
	Double_t qout_rot = GetQout(GoodTrackFourVector[itk1],InvertXYVector(GoodTrackFourVector[itk2]));
	Double_t qside_inv = GetQside(GoodTrackFourVector[itk1],InvertPVector(GoodTrackFourVector[itk2]));
	Double_t qside_rot = GetQside(GoodTrackFourVector[itk1],InvertXYVector(GoodTrackFourVector[itk2]));
	
	TLorentzVector psum2 = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
	Double_t kt=(psum2.Pt())/2.;
	TLorentzVector psum2_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector[itk2]);
	Double_t kt_rot = (psum2_rot.Pt())/2.;
	TLorentzVector psum2_inv = GoodTrackFourVector[itk1] + InvertPVector(GoodTrackFourVector[itk2]);
	Double_t kt_inv = (psum2_inv.Pt())/2.;

	Double_t x5D_LCMS[5]={qlongLCMS,qout,qside,kt,(Double_t)mult};
	Double_t x5D_LCMS_inv[5]={qlongLCMS_inv,qout_inv,qside_inv,kt_inv,(Double_t)mult};
	Double_t x5D_LCMS_rot[5]={qlongLCMS_rot,qout_rot,qside_rot,kt_rot,(Double_t)mult};

	Double_t x3D[3]={q,kt,(Double_t)mult};
	Double_t x3D_rot[3]={qs,kt_rot,(Double_t)mult};
	Double_t x3D_inv[3]={qo,kt_inv,(Double_t)mult};


	//Double_t aux_tk2_corr = 1.0; // No track corrections
	Double_t aux_tk2_corr = (Double_t)getTrkCorrWeight(eff_histos[0], eff_histos[1], eff_histos[2], eff_histos[3], GoodTrackFourVector[itk2].Pt(),GoodTrackFourVector[itk2].Eta()); 
	Double_t aux_tk12_corr=aux_tk1_corr*aux_tk2_corr;

	//Cesar: why only for qinv<0.01GeV? 
	if (GoodTrackCharge[itk1]==1 && GoodTrackCharge[itk2]==1   && q < 0.01)hist_dpt_cos_pp->Fill(costheta,deltapt1);
	if (GoodTrackCharge[itk1]==-1 && GoodTrackCharge[itk2]==-1 && q < 0.01)hist_dpt_cos_mm->Fill(costheta,deltapt1);

	if(GoodTrackCharge[itk1]*GoodTrackCharge[itk2]>0){
	  if(q < 0.01)hist_deetadphi->Fill(d_eta, d_phi); //Cesar: why only for qinv<0.01GeV?
	  hist_tk_pairSS_M_->Fill(psum2.M(),aux_tk12_corr);
	  
	  hist_sig_KtVsNch_ntrkoff->Fill(kt,mult,aux_tk12_corr);
	  hist_sig_pairSS_MVsNch_ntrkoff->Fill(psum2.M(),mult,aux_tk12_corr);
	  
	  //for 3D analysis                                                                                                         
	  hist_qschVsKtVsNch_ntrkoff->Fill(x3D,CoulombW(q)*aux_tk12_corr);
	  hist_qschncVsKtVsNch_ntrkoff->Fill(x3D,aux_tk12_corr);
	  hist_qINVschVsKtVsNch_ntrkoff->Fill(x3D_inv,CoulombW(q)*aux_tk12_corr);
	  hist_qINVschncVsKtVsNch_ntrkoff->Fill(x3D_inv,aux_tk12_corr);
	  hist_qROTschVsKtVsNch_ntrkoff->Fill(x3D_rot,CoulombW(q)*aux_tk12_corr);
	  hist_qROTschncVsKtVsNch_ntrkoff->Fill(x3D_rot,aux_tk12_corr);
	  //for 5D analysis                                                                                                         
	  hist_qLLCMSqOqSschVsKtVsNch_ntrkoff->Fill(x5D_LCMS,CoulombW(q)*aux_tk12_corr);
	  hist_qLLCMSqOqSschncVsKtVsNch_ntrkoff->Fill(x5D_LCMS,aux_tk12_corr);
	  hist_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff->Fill(x5D_LCMS_inv,CoulombW(q)*aux_tk12_corr);
	  hist_qLLCMSqOqS_INV_schncVsKtVsNch_ntrkoff->Fill(x5D_LCMS_inv,aux_tk12_corr);
	  hist_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff->Fill(x5D_LCMS_rot,CoulombW(q)*aux_tk12_corr);
	  hist_qLLCMSqOqS_ROT_schncVsKtVsNch_ntrkoff->Fill(x5D_LCMS_rot,aux_tk12_corr);
	
	}else{
	  //now, the same for opposite charge                                                                                       
	  hist_tk_pairOS_M_->Fill(psum2.M(),aux_tk12_corr);
	  hist_pairOS_MVsNch_ntrkoff->Fill(psum2.M(),mult,aux_tk12_corr);
	  //for 3D analysis                                                                                                         
	  hist_qdchVsKtVsNch_ntrkoff->Fill(x3D,CoulombWpm(q)*aux_tk12_corr);
	  hist_qdchncVsKtVsNch_ntrkoff->Fill(x3D,aux_tk12_corr);
	  //for 5D analysis                                                                                                         
	  hist_qLLCMSqOqSdchVsKtVsNch_ntrkoff->Fill(x5D_LCMS,CoulombWpm(q)*aux_tk12_corr);
	  hist_qLLCMSqOqSdchncVsKtVsNch_ntrkoff->Fill(x5D_LCMS,aux_tk12_corr);

	}
      }
    } // End of Loop for tracks
    ev_ntrkoff_vec.push_back(mult);
    ev_GoodTrackFourVector_vec.push_back(GoodTrackFourVector); 
    ev_GoodTrackCharge_vec.push_back(GoodTrackCharge);
    ev_vtx_z_vec.push_back(vertexz);
    ///for eta-mixing
    //build vector of pairs ordered by the etaMixWeight
    Double_t aux_etaMix_w = ComputeEventWeight(GoodTrackFourVector_trkoff);
    //    std::cout<<"aux_etaMix_w : "<<aux_etaMix_w<<std::endl;
    std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight = make_pair(aux_etaMix_w, GoodTrackFourVector);
    ev_GoodTrackFourVector_etaMixWeight_vec.push_back(aux_pair_GoodTrackFourVector_etaMixWeight);
    std::pair<Double_t, std::vector<int>> aux_pair_GoodTrackCharge_etaMixWeight = make_pair(aux_etaMix_w,GoodTrackCharge);
    ev_GoodTrackCharge_etaMixWeight_vec.push_back(aux_pair_GoodTrackCharge_etaMixWeight);
    std::pair<Double_t,int> aux_pair_ntrkoff_etaMixWeight = make_pair(aux_etaMix_w,mult);//IMPORTANT: (Cesar) if you are doing analysis in terms of centrality, is very important you revise this.
    ev_ntrkoff_etaMixWeight_vec.push_back(aux_pair_ntrkoff_etaMixWeight);

  } // End loop over events

  //--------------------------------------------------------------------------------------------------

  int ntrkoff_min = 0; ///just example will need to define several ranges and number of events to mix after
  int ntrkoff_max = 400;
  int nEvt_to_mix = 10;
  
  int MixingID=2;//Cesar: should implement when calling the function ...just temporary here
  if(isEventMix==1 && MixingID==1){
    cout<< " I am starting the random mixing " << endl;
    call_mix_random(ntrkoff_min, ntrkoff_max, nEvt_to_mix, ev_ntrkoff_vec, ev_vtx_z_vec, ev_GoodTrackFourVector_vec, ev_GoodTrackCharge_vec);
  }

  if(isEventMix==1 && MixingID==2){
    cout<< " I am starting the eta mixing " << endl;	  
    call_mix_eta(ntrkoff_min, ntrkoff_max, ev_GoodTrackFourVector_etaMixWeight_vec, ev_GoodTrackCharge_etaMixWeight_vec, ev_ntrkoff_etaMixWeight_vec);
  } 	  
  
  //---------------------------------------------------------------------------------------------------  

  cout <<  " " <<endl;
  cout << "Writing histograms on " << endl;
  cout << " " << endl;

  // Make an output file
  string file_output= Form("%s_%s_%iGeV_%s_%i_%s",colliding_system.Data(), data_or_mc.Data(), sNN_energy_GeV, ref_sample.Data(), date->GetDate(),ouputfilename.Data());
  std::replace(file_output.begin(), file_output.end(), '.', 'p'); // replace . to p
  std::replace(file_output.begin(), file_output.end(), '-', 'N'); // replace - to N for negative
  
  //----------------------------------------------------------------------------------------------------
  
  // Open, write and close the output file

  TFile *MyFile = new TFile(Form("%s.root", file_output.c_str()), "RECREATE");
  if(MyFile->IsOpen()) cout << "output file: " << file_output.c_str() << ".root" << endl;
  MyFile->cd(); 
  // Write in different folders (see histogram_definition.h)___ read this
  MyFile->mkdir("QA_histograms"); 
  MyFile->cd("QA_histograms"); 
  QA_hist();
  
  MyFile->mkdir("SameEventHistograms"); 
  MyFile->cd("SameEventHistograms");
  SameEvents_hist();
  
  MyFile->mkdir("MixedEventHistograms"); 
  MyFile->cd("MixedEventHistograms");
  MixEvents_hist();
  
  MyFile->Close();
  cout << endl;
  cout << "------------------------------------- DONE --------------------------------------" << endl;
  cout << endl;
  
  
  sec_end = clock(); // stop time counting
  cout << "========================================" << endl;
  cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
  cout << "========================================" << endl;
  
  print_stop(); // Print time, date and hour when it stops

}
