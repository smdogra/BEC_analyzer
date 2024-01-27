#include "call_libraries.h"  // call libraries from ROOT and C++
#include "input_variables.h" // call inputs

Double_t pi_mass  = 0.1396;
Double_t cos_cut = 0.99996;
Double_t dpt_cut = 0.04;
Int_t bins5D[5]=   {40, 40, 40, 10, 100  };
Double_t xmin5D[5]={0., 0., 0., 0., 0.  };
Double_t xmax5D[5]={2., 2., 2., 1., 1000.};

Int_t bins3D[3]=   {200,20, 1000  };
Double_t xmin3D[3]={0. ,0.,-0.5 };
Double_t xmax3D[3]={2. ,1.,999.5};

int trkbinsize = (int) trk_pt_bins.size(); // track bins for jet-track correlation
int multbinsize = (int) multiplicity_centrality_bins.size();// multiplicity or centrality bins for jet-track correlation

// -------------------------------- QA plots --------------------------------
// Event quantities
TH1I *Nevents = new TH1I("Nevents", "Nevents", 10, 0, 10);
TH1D *multiplicity = new TH1D("multiplicity", "multiplicity", 500, 0.0, 500.0);
TH1D *multiplicity_weighted = new TH1D("multiplicity_weighted", "multiplicity_weighted", 500, 0.0, 500.0);
TH1D *vzhist = new TH1D("vzhist", "vzhist", 150, -15, 15);
TH1D *vzhist_weighted = new TH1D("vzhist_weighted", "vzhist_weighted", 150, -15, 15);


TH1D *hist_trketa = new TH1D("hist_trketa",  "hist_trketa", 400, -4.0, 4.0);
TH1D *hist_trkpt  = new TH1D("hist_trkpt" ,  "hist_trkpt", 1000, 0, 100);
TH1D *hist_trkcharge= new TH1D("hist_trkcharge", "hist_trkcharge", 3, -1.5, 1.5);
TH1D *hist_trkpterrbytrkpt= new TH1D("hist_trkpterrbytrkpt", "hist_trkpterrbytrkpt", 1000, 0, 10);
TH1D *hist_trkNPixelHit= new TH1D("hist_trkNPixelHit", "hist_trkNPixelHit",10, -0.5, 9.5);
TH1D *hist_trkdcaxybytrkerrdcaxy = new TH1D("hist_trkdcaxybytrkerrdcaxy" , "hist_trkdcaxybytrkerrdcaxy", 100, -5, 5);
TH1D *hist_trkdcazbytrkdcazerr = new TH1D("hist_trkdcazbytrkdcazerr" , "hist_trkdcazbytrkdcazerr",100, -5, 5);
TH1D *hist_trkchi2bytrkndof = new TH1D("hist_trkchi2bytrkndof" , "hist_trkchi2bytrkndof",1000,0,100);
TH1D *hist_trkchi2bytrkndofbytrknlayer = new TH1D("hist_trkchi2bytrkndofbytrknlayer" , "hist_trkchi2bytrkndofbytrknlayer", 1000 ,0,100);
TH1D *hist_trknlayer = new TH1D("hist_trknlayer" , "hist_trknlayer",50,-0.5,49.5);
TH1D *hist_trknhits = new TH1D("hist_trknhits" , "hist_trknhits",50,-0.5,49.5);
TH1D *hist_trkalgo = new TH1D("hist_trkalgo" , "hist_trkalgo",100,-0.5,10);
TH1D *hist_trkmva = new TH1D("hist_trkmva" , "hist_trkmva",100,-0.5,99.5);

TH2F *hist_ptVsNch_ntrkoff = new TH2F("hist_ptVsNch_ntrkoff", "hist_ptVsNch_ntrkoff",1000, 0,100, 400, -0.5, 399.5);
TH2F *hist_etaVsNch_ntrkoff = new TH2F("hist_etaVsNch_ntrkoff", "hist_etaVsNch_ntrkoff",400, -4.0, 4.0, 400, -0.5, 399.5);
TH2F *hist_phiVsNch_ntrkoff = new TH2F("hist_phiVsNch_ntrkoff", "hist_phiVsNch_ntrkoff",400, -4.0, 4.0, 400, -0.5, 399.5);

TH2F *hist_dpt_cos_pp = new TH2F("hist_dpt_cos_pp", "hist_dpt_cos_pp",1000, 0.99910, 1.0001, 200, -0.5, 0.5);
TH2F *hist_dpt_cos_mm = new TH2F("hist_dpt_cos_mm", "hist_dpt_cos_mm",1000, 0.99910, 1.0001, 200, -0.5, 0.5);
TH2F *hist_dpt_cos_pm = new TH2F("hist_dpt_cos_pm", "hist_dpt_cos_pm",1000, 0.99910, 1.0001, 200, -0.5, 0.5);


TH2F *hist_deetadphi = new TH2F("hist_deetadphi", "hist_deetadphi", 1000, -6, 6, 1000, -4.0, 4.0);
TH1F *hist_tk_pairSS_M_ = new TH1F("hist_tk_pairSS_M_", "Invariant mass same-sign tracks", 1000, 0, 1);
 
TH2F *hist_sig_KtVsNch_ntrkoff= new TH2F("sig_KtVsNch_ntrkoff","sig_KtVsNch_ntrkoff",200,0.,20.,1000,-0.5,999.5);
TH2F *hist_sig_pairSS_MVsNch_ntrkoff= new TH2F("sig_pairSS_MVsNch_ntrkoff","sig_pairSS_MVsNch_ntrkoff",200,0,2,1000,-0.5,999.5);


TH1F *hist_DeltazVertex_=new TH1F("DeltazVertex_","DeltazVertex_", 100, -2.5, 2.5); //between two events after mixing                                        
TH2F *hist_mixsch_KtVsNch_ntrkoff=new TH2F("mixsch_KtVsNch_ntrkoff","mixsch_KtVsNch_ntrkoff",200,0.,20.,1000,-0.5,999.5);
THnSparse *hist_qinv_mixschVsKtVsNch_ntrkoff=new THnSparseD("hs_qinv_mixschVsKtVsNch_ntrkoff","hs_qinv_mixschVsKtVsNch_ntrkoff",3, bins3D,  xmin3D, xmax3D); //MixingRef Same-sign (NocoulombCorr)                    
THnSparse *hist_qinv_mixdchVsKtVsNch_ntrkoff=new THnSparseD("hs_qinv_mixdchVsKtVsNch_ntrkoff","hs_qinv_mixdchVsKtVsNch_ntrkoff",3, bins3D,  xmin3D, xmax3D); //MixingRef Opp-sign (NoCoulombCorr)    
THnSparse *hist_qLLCMSqOqSmixschVsKtVsNch_ntrkoff=new THnSparseD("hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff","hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D); //MixingRef Same-sign (NocoulombCorr)              
THnSparse *hist_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff=new THnSparseD("hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff","hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D); //MixingRef Opp-sign (NoCoulombCorr)


THnSparseD *hist_qschVsKtVsNch_ntrkoff =new THnSparseD("hist_qschVsKtVsNch_ntrkoff","hist_qschVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);

THnSparseD *hist_qschncVsKtVsNch_ntrkoff=new THnSparseD("hist_qschncVsKtVsNch_ntrkoff","hist_qschncVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
THnSparseD *hist_qINVschVsKtVsNch_ntrkoff=new THnSparseD("hist_qINVschVsKtVsNch_ntrkoff","hist_qINVschVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
THnSparseD *hist_qINVschncVsKtVsNch_ntrkoff=new THnSparseD("hist_qINVschncVsKtVsNch_ntrkoff","hist_qINVschncVsKtVsNch_ntrkoff",3,bins3D,xmin3D);
THnSparseD *hist_qROTschVsKtVsNch_ntrkoff=new THnSparseD("hist_qROTschVsKtVsNch_ntrkoff","hist_qROTschVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
THnSparseD *hist_qROTschncVsKtVsNch_ntrkoff=new THnSparseD("hist_qROTschncVsKtVsNch_ntrkoff","hist_qROTschncVsKtVsNch_ntrkoff" ,3,bins3D,xmin3D,xmax3D);                                                                                                                   
THnSparseD *hist_qLLCMSqOqSschVsKtVsNch_ntrkoff=new THnSparseD("hist_qLLCMSqOqSschVsKtVsNch_ntrkoff","hist_qLLCMSqOqSschVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
THnSparseD *hist_qLLCMSqOqSschncVsKtVsNch_ntrkoff=new THnSparseD("hist_qLLCMSqOqSschncVsKtVsNch_ntrkoff","hist_qLLCMSqOqSschncVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
THnSparseD *hist_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff=new THnSparseD("hist_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff","hist_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
THnSparseD *hist_qLLCMSqOqS_INV_schncVsKtVsNch_ntrkoff=new THnSparseD("hist_qLLCMSqOqS_INV_schncVsKtVsNch_ntrkoff","hist_qLLCMSqOqS_INV_schncVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
THnSparseD *hist_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff=new THnSparseD("hist_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff","hist_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
THnSparseD *hist_qLLCMSqOqS_ROT_schncVsKtVsNch_ntrkoff=new THnSparseD("hist_qLLCMSqOqS_ROT_schncVsKtVsNch_ntrkoff","hist_qLLCMSqOqS_ROT_schncVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
TH1F *hist_tk_pairOS_M_=new TH1F("hist_tk_pairOS_M_", "Invariant mass opposite-sign tracks", 1000, 0, 1);
TH2F *hist_pairOS_MVsNch_ntrkoff=new TH2F("hist_pairOS_MVsNch_ntrkoff","hist_pairOS_MVsNch_ntrkoff",200,0,2,1000,-0.5,999.5);
THnSparseD *hist_qdchVsKtVsNch_ntrkoff=new THnSparseD("hist_qdchVsKtVsNch_ntrkoff","hist_qdchVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
THnSparseD *hist_qdchncVsKtVsNch_ntrkoff=new THnSparseD("hist_qdchncVsKtVsNch_ntrkoff","hist_qdchncVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
THnSparseD *hist_qLLCMSqOqSdchVsKtVsNch_ntrkoff=new THnSparseD("hist_qLLCMSqOqSdchVsKtVsNch_ntrkoff","hist_qLLCMSqOqSdchVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
THnSparseD *hist_qLLCMSqOqSdchncVsKtVsNch_ntrkoff=new THnSparseD("hist_qLLCMSqOqSdchncVsKtVsNch_ntrkoff","hist_qLLCMSqOqSdchncVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D); 


// Track/Particle histograms
int    bins4D_trk[4]   =   { 500   ,  50  ,   64           , multbinsize-1};
double xmin4D_trk[4]   =   { 0.0   , -2.5 ,   -TMath::Pi() , 0};
double xmax4D_trk[4]   =   { 20.0  ,  2.5 ,   TMath::Pi()  , (double) multbinsize-1};

// --> Reco
THnSparseD *hist_reco_trk = new THnSparseD("hist_reco_trk", "hist_reco_trk", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_reco_trk_corr = new THnSparseD("hist_reco_trk_corr", "hist_reco_trk_corr", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_reco_trk_weighted = new THnSparseD("hist_reco_trk_weighted", "hist_reco_trk_weighted", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);

// --> Gen
THnSparseD *hist_gen_trk = new THnSparseD("hist_gen_trk", "hist_gen_trk", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_gen_trk_weighted = new THnSparseD("hist_gen_trk_weighted", "hist_gen_trk_weighted", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);


/*
bool splitcomb(TLorentzVector &vec1,TLorentzVector &vec2){
  bool issplit=false;
  Double_t cosa = TMath::Abs(vec1.Px()*vec2.Px() + vec1.Py()*vec2.Py() + vec1.Pz()*vec2.Pz())/(vec1.P()*vec2.P());
  Double_t deltapt = TMath::Abs(vec1.Pt() - vec2.Pt());
  if ( (cosa >cos_cut) && (deltapt < dpt_cut)) { issplit = true;}
  //std::cout << "cosa: " << cosa << " dpt: " << deltapt << " is split: " << issplit << std::endl;                               
  return issplit;
}

Double_t GetQ(const TLorentzVector &p1, const TLorentzVector &p2){
  TLorentzVector Sum4V = p1+p2;
  Double_t q = Sum4V.Mag2() - 4*pi_mass*pi_mass;
  //  std::cout<<(  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  ) <<std::endl;                                                     
  return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}
Double_t GetQlong(const TLorentzVector &p1, const TLorentzVector &p2){
  TLorentzVector Diff4V = p1-p2;
  Double_t qlong = fabs(Diff4V.Pz());
  return qlong;
}

Double_t GetQlongLCMS(const TLorentzVector &p1, const TLorentzVector &p2){
  Double_t num = 2*( (p1.Pz())*(p2.E()) - (p2.Pz())*(p1.E()) );
  Double_t den = TMath::Sqrt( (p1.E()+p2.E())*(p1.E()+p2.E()) - (p1.Pz()+p2.Pz())*(p1.Pz()+p2.Pz()) );
  Double_t qlongLCMS = 0.0;
  if(den != 0) qlongLCMS = fabs(num/den);
  return qlongLCMS;
}

Double_t GetQout(const TLorentzVector &p1, const TLorentzVector &p2){
  TVector3 qT;
  qT.SetXYZ(p1.Px()-p2.Px(),p1.Py()-p2.Py(),0.0);
  TVector3 kT;
  kT.SetXYZ( (p1.Px()+p2.Px())/2.0 , (p1.Py()+p2.Py())/2.0 ,0.0);
  TVector3 qout;
  qout = qT.Dot(kT.Unit())*kT.Unit();

  Double_t absValue = qout.Mag();
  ///if (qT.Dot(kT.Unit()) < 0)signed_absValue=-1.*qout.Mag();                                                                   

  return absValue;


}

Double_t GetQside(const TLorentzVector &p1, const TLorentzVector &p2){

  TVector3 qT;
  qT.SetXYZ(p1.Px()-p2.Px(),p1.Py()-p2.Py(),0.0);
  TVector3 kT;
  kT.SetXYZ( (p1.Px()+p2.Px())/2.0 , (p1.Py()+p2.Py())/2.0 ,0.0);
  TVector3 qout;
  qout = qT.Dot(kT.Unit())*kT.Unit();
  TVector3 qsid;
  qsid = qT - qout;

  Double_t absValue = qsid.Mag();

  return absValue;

}
                                   

const TLorentzVector InvertPVector( TLorentzVector &vec){
  TLorentzVector ovec = vec;
  ovec.SetPx(-vec.Px());
  ovec.SetPy(-vec.Py());
  ovec.SetPz(-vec.Pz());
  return ovec;
}

const TLorentzVector InvertXYVector( TLorentzVector &vec){
  TLorentzVector ovec = vec;
  ovec.SetX(-vec.X());
  ovec.SetY(-vec.Y());
  return ovec;
}

//return the weight factor due to Coloumb repulsion [Gamow] same charge                                                         
const Double_t CoulombW(const Double_t& q){
  const Double_t alpha=1./137.;
  Double_t x=2*TMath::Pi()*(alpha*pi_mass/q);
  //return (TMath::Exp(x)-1.)/x; // OLD MATTIA's DEFINITION                                                                     

  //Double_t ws=scf*((exp(arg)-1)/arg-1)+1; // PAOLO's DEFINITION                                                               
  Double_t weight = 1;//0.85; // TEMPORARY SET TO 0.85 * GAMOW FACTOR                                                           
  //Double_t weight = 1.15; //for syst. +15%                                                                                    
  //Double_t weight = 0.85; //for syst. -15%                                                                                    
  return weight*( (TMath::Exp(x)-1.)/x -1 ) + 1;
 }

//return the weight factor due to Coloumb attraction [Gamow] opposite charge                                                    
const  Double_t CoulombWpm(const Double_t& q){
  const Double_t alpha=1./137.;
  Double_t x=2*TMath::Pi()*(alpha*pi_mass/q);
  // return (1.-TMath::Exp(-x))/x; // OLD MATTIA's DEFINITION                                                                   

  // Double_t wd=scf*((1-exp(-arg))/arg-1)+1; // PAOLO's DEFINITION                                                             
  Double_t weight = 1;//0.85; // TEMPORARY SET TO 0.85 * GAMOW FACTOR                                                           
  //Double_t weight = 1.15; //for syst. +15%                                                                                    
  //Double_t weight = 0.85; //for syst. -15%                                                                                    
  return weight*( (1.-TMath::Exp(-x))/x -1 ) + 1;
}
*/






// Evaluate uncertainties correctly at ROOT
void sw2(){
	Nevents->Sumw2();
	hist_trketa->Sumw2();
	hist_trkpt->Sumw2();
	hist_trkcharge->Sumw2();
	hist_trkpterrbytrkpt->Sumw2();

	hist_trkdcaxybytrkerrdcaxy->Sumw2();
	hist_trkdcazbytrkdcazerr->Sumw2();
       	hist_trkchi2bytrkndof->Sumw2();
	hist_trkchi2bytrkndofbytrknlayer->Sumw2();
	hist_trknlayer->Sumw2();
	hist_trknhits->Sumw2();
	hist_trkalgo->Sumw2();
	hist_trkmva->Sumw2();
	
}

// write QA histograms
/*
--> Arguments
isMC: true for MC and false for Data
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void QA_hist(){
  Nevents->Write();
  multiplicity->Write();
  multiplicity_weighted->Write();
  vzhist->Write();
  vzhist_weighted->Write();
  
  hist_trketa->Write();               
  hist_trkpt->Write();                                                                                                           
  hist_trkcharge->Write();                                                                                                       
  hist_trkpterrbytrkpt->Write();      
  hist_trkNPixelHit->Write();
  hist_trkdcaxybytrkerrdcaxy->Write();
  hist_trkdcazbytrkdcazerr->Write();
  hist_trkchi2bytrkndof->Write();
  hist_trkchi2bytrkndofbytrknlayer->Write();
  hist_trknlayer->Write();
  hist_trknhits->Write();
  hist_trkalgo->Write();
  hist_trkmva->Write();

  hist_ptVsNch_ntrkoff -> Write();
  hist_etaVsNch_ntrkoff -> Write();
  hist_phiVsNch_ntrkoff -> Write();

  hist_dpt_cos_pp -> Write();
  hist_dpt_cos_mm -> Write();
  hist_dpt_cos_pm -> Write();

  hist_deetadphi -> Write();
  hist_tk_pairSS_M_ -> Write();

  hist_sig_KtVsNch_ntrkoff-> Write();
  hist_sig_pairSS_MVsNch_ntrkoff-> Write();
  hist_tk_pairOS_M_-> Write();
  hist_pairOS_MVsNch_ntrkoff-> Write();

}

void SameEvents_hist(){

  hist_qschVsKtVsNch_ntrkoff-> Write();
  hist_qschncVsKtVsNch_ntrkoff-> Write();
  hist_qINVschVsKtVsNch_ntrkoff-> Write();
  hist_qINVschncVsKtVsNch_ntrkoff-> Write();
  hist_qROTschVsKtVsNch_ntrkoff-> Write();
  hist_qROTschncVsKtVsNch_ntrkoff-> Write();                                                                                              

  hist_qLLCMSqOqSschVsKtVsNch_ntrkoff-> Write();
  hist_qLLCMSqOqSschncVsKtVsNch_ntrkoff-> Write();
  hist_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff-> Write();
  hist_qLLCMSqOqS_INV_schncVsKtVsNch_ntrkoff-> Write();
  hist_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff-> Write();
  hist_qLLCMSqOqS_ROT_schncVsKtVsNch_ntrkoff-> Write();

  hist_qdchVsKtVsNch_ntrkoff-> Write();
  hist_qdchncVsKtVsNch_ntrkoff-> Write();
  hist_qLLCMSqOqSdchVsKtVsNch_ntrkoff-> Write();
  hist_qLLCMSqOqSdchncVsKtVsNch_ntrkoff-> Write();
}

void MixEvents_hist(){
  hist_DeltazVertex_->Write();                                      
  hist_mixsch_KtVsNch_ntrkoff->Write();
  hist_qinv_mixschVsKtVsNch_ntrkoff ->Write();                                                             
  hist_qinv_mixdchVsKtVsNch_ntrkoff->Write();
  hist_qLLCMSqOqSmixschVsKtVsNch_ntrkoff->Write();
  hist_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff->Write();
}
