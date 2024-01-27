#include "call_libraries.h"  // call libraries from ROOT and C++
#include "trk_efficiency_correction.h" // track efficiency correction
#include "weights.h" // weights applied

/*
Find Ntrk offline -> updated for all systems (and easy to update for future systems)
The Ntrk offline is a definition with specific cuts (we should not change it). The track systematics must be applied using the input_variables.h!
--> Arguments
col_sys: colliding system
col_energy:  colliding energy
yearofdatataking: year of data-taking
size: track collection size per event
eta: track eta
pt: track pT
charge: track charge
hp: track high purity workflow
pterr: track pT uncertainty
dcaxy: track DCA in the transverse plane
dcaxyerr: track DCA in the transverse plane uncertainty
dcaz: track DCA in the longitudinal plane
dcazerr: track DCA in the longitudinal plane uncertainty
chi2: track chi2 of reconstruction
ndof: track number of degrees of freedom reconstruction
nlayer: track number of layers with measurements
nhits: track number of hits with measurements
algo: track MVA algorith step
mva: track MVA algorith value [-1,1]
*/
int get_Ntrkoff(TString col_sys, int col_energy, int yearofdatataking, int size, float *eta, float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr, float* chi2, int* ndof, int* nlayer, int* nhits, int* algo, float* mva){
	float Ntrk_off = 0;
	for(int ii=0; ii<size; ii++){ 
		if(fabs(eta[ii]) > 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] == false) continue;
		if(fabs(pterr[ii]/pt[ii]) >= 0.1) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 3.0) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= 3.0) continue;
		double calomatching = ((pfEcal[ii]+pfHcal[ii])/cosh(eta[ii]))/pt[ii];
		if(col_sys=="pPb" && col_energy==8160 && yearofdatataking==2016){ if(pt[ii] <= 0.4) continue;}
		if(col_sys=="pp" && col_energy==5020 && yearofdatataking==2017) if(pt[ii] <= 0.5) continue; 
		if(col_sys=="pp" && col_energy==13000 && yearofdatataking==2017) if(pt[ii] <= 0.5) continue; 
		if(col_sys=="XeXe" && col_energy==5440 && yearofdatataking==2017){
			if(pt[ii] <= 0.5) continue; 
			if((chi2[ii]/ndof[ii])/nlayer[ii] >= 0.15) continue;
		 	if(nhits[ii] < 11) continue;
		 	if(pt[ii] > 20.0){if(calomatching<=0.5)continue;} //is this applicable in pp or pPb?
		}
		if(col_sys=="PbPb" && col_energy==5020 && yearofdatataking==2018){
			if(pt[ii] <= 0.5) continue; 
			if((chi2[ii]/ndof[ii])/nlayer[ii] >= 0.18) continue;
		 	if(nhits[ii] < 11) continue;
		 	if(pt[ii] > 20.0){if(calomatching<=0.5)continue;} //is this applicable in pp or pPb?
			if(algo[ii]==6 && mva[ii]<0.98) continue;
		}
		Ntrk_off=Ntrk_off+1;
	}
	return Ntrk_off;
}

/*
Calculate jet Aj asymmetry
--> Arguments
pt_leading: leading jet pT
pt_subleading: subleading jet pT
*/
float asymmetry(float pt_leading, float pt_subleading){
	float Avariable = (pt_leading - pt_subleading)/(pt_leading + pt_subleading);
	return Avariable;
}

/*
Calculate jet Xj asymmetry
--> Arguments
pt_leading: leading jet pT
pt_subleading: subleading jet pT
*/
float xjvar(float pt_leading, float pt_subleading){
	float XJvariable = pt_subleading/pt_leading;
	return XJvariable;
}

/*
Calculate Delta Eta
--> Arguments
eta1: eta of first object
eta2: eta of second object
*/
float deltaeta(float eta1, float eta2){
	float deltaEta = ( eta1 - eta2 );
	return deltaEta;
}

/*
Calculate Delta Phi
--> Arguments
phi1: eta of first object
phi2: eta of second object
*/
float deltaphi(float phi1, float phi2){
	float deltaPhi = ( phi1 - phi2 );
	return deltaPhi;
}

/*
Calculate Delta Phi in the range [-pi/2 , 3/2 pi]
--> Arguments
phi1: eta of first object
phi2: eta of second object
*/
float deltaphi2PC(float phi1, float phi2){     
	float deltaPhi = (phi1 - phi2);
	if( deltaPhi > TMath::Pi() ) deltaPhi = deltaPhi - 2.*TMath::Pi();
	if( deltaPhi <= -TMath::Pi()/2.) deltaPhi = deltaPhi + 2.*TMath::Pi();
	return deltaPhi;
}

/*
Calculate delta R (distance)
--> Arguments
eta1: eta of first object
phi1: eta of first object
eta2: eta of second object
phi2: eta of second object
*/
float deltaR(float eta1, float phi1, float eta2, float phi2){
	float deltaR = sqrt(pow(deltaeta(eta1,eta2),2) + pow(deltaphi(phi1,phi2),2));
	return deltaR;
}

/*
Find the leading and subleading jets, return jet leading and subleading pt, eta and phi
--> Arguments
pt: jet pT
eta: jet Eta
phi: jet Phi
leadpt: leading jet pT
leadeta: leading jet Eta
leadphi: leading jet Phi
sublpt: subleading jet pT
subleta: subleading jet Eta
sublphi: subleading jet Phi
*/
void find_leading_subleading(float pt, float eta, float phi, float &leadpt, float &leadeta, float &leadphi, float &sublpt, float &subleta, float &sublphi){
	if( pt > leadpt ) {
    	sublpt = leadpt;
        leadpt = pt;
        leadeta = eta;
        leadphi = phi;
    } else if( sublpt < pt) {
    	sublpt = pt;
        subleta = eta;
        sublphi = phi;
    }
}

/*
Find bin dynamically 
--> Arguments
quant_vec: vector with binning
quant: variable
*/
int find_my_bin(std::vector<double> quant_vec, double quant){
	int bin = -999;
	for(int ii = 0; ii < quant_vec.size()-1; ii++) {if(quant >= quant_vec[ii] && quant < quant_vec[ii+1]){bin = ii;} }
    return bin;
}

/*
Measure the correlation between objects
--> Arguments
jets: vector with jet informations
jet_w: vector with jet weight informations
tracks: vector with track informations
trk_w: vector with track weight informations
histo_corr: multidimentional histogram for correlations {Delta Phi, Delta Eta, track pT bin, multiplicity or centrality}
histjet: multidimentional histogram for jets in correlations {pT, Eta, Phi}
histtrk: multidimentional histogram for tracks in correlations {pT, Eta, Phi}
event_weight: event weight vector for each event
mult: multiplicity or centrality vector for each event
do_rotation: true means apply/fill rotation method, otherwise use false
N_rot: number of rotation (only use if do_rotation is true)
histo_rot: histogram using rotation method (only use if do_rotation is true)
sube_trk: vector with sube track (MC embedded samples) sube == 0 means PYTHIA embedded tracks while sube > 0 means the other MC (HYDJET, EPOS, ...)
histo_corr_subeg0: if sube > 0 save in this histogram
*/
void correlation(std::vector<TVector3> jets, std::vector<double> jets_w, std::vector<TVector3> tracks, std::vector<double> tracks_w, THnSparse* histo_corr, THnSparse* histjet, THnSparse* histtrk, float event_weight, int mult, bool do_rotation, int N_rot, THnSparse* histo_rot, std::vector<int> sube_trk, THnSparse* histo_corr_subeg0){
	// get correlation histograms
	for (int a = 0; a < jets.size(); a++){ // start loop over jets
        double jet_weight = jets_w[a];
		for (int b = 0; b < tracks.size(); b++){ // start loop over tracks
			double trkpt = tracks[b].Pt();
			double trketa = tracks[b].Eta();
			int subetrk = sube_trk[b];
            // track efficiency correction for reco
            double trk_weight = tracks_w[b];
            // Find track and multiplicity bins
			int trkbin = (int) find_my_bin(trk_pt_bins,trkpt);
			int multcentbin = (int) find_my_bin(multiplicity_centrality_bins, (float) mult);
			// Fill jet and track quantities
			double x4D_jet[4]={jets[a].Pt(),jets[a].Eta(),jets[a].Phi(), (double)multcentbin}; histjet->Fill(x4D_jet,jet_weight*event_weight);
			double x4D_trk[4]={tracks[b].Pt(),tracks[b].Eta(),tracks[b].Phi(), (double)multcentbin}; histtrk->Fill(x4D_trk,trk_weight*event_weight);
			// Fill correlation histograms
			double del_phi = deltaphi2PC(jets[a].Phi(), tracks[b].Phi());
			double del_eta = deltaeta(jets[a].Eta(), tracks[b].Eta());
			double x4D[4]={del_phi,del_eta,(double)trkbin,(double)multcentbin}; 
			if(subetrk==0){histo_corr->Fill(x4D,jet_weight*trk_weight*event_weight*trkpt);}else{histo_corr_subeg0->Fill(x4D,jet_weight*trk_weight*event_weight*trkpt);}
		}
		// get rotation histograms 
		if(do_rotation){
			for(int c = 0; c < N_rot; c++){
				TRandom2 *r = new TRandom2(); // make a random number
				float alpha = r->Uniform(2.0*TMath::Pi()); // make a random number between 0 and 2pi
				float newphi = jets[a].Phi() + alpha; // add to the jet phi and reorder the axis
				if (newphi > 3.0*TMath::Pi()/2.) newphi -= 2.0*TMath::Pi();
				if (newphi < -TMath::Pi()/2.) newphi += 2.0*TMath::Pi();
				float neweta = -1.0*jets[a].Eta(); // invert the jet eta
				for (int d = 0; d < tracks.size(); d++){ // start loop over tracks using the new jet eta and phi (similar as above)
					double trkpt = tracks[d].Pt();
					double trketa = tracks[d].Eta();
           			// track efficiency correction for reco
           			double trk_weight = tracks_w[d];
            		// Find track and multiplicity bins
					int trkbin = (int) find_my_bin(trk_pt_bins,trkpt);
					int multcentbin = (int) find_my_bin(multiplicity_centrality_bins, (float)mult);
					// Fill correlation histograms     
					double del_phi_rot = deltaphi2PC(newphi, tracks[d].Phi());
					double del_eta_rot = deltaeta(neweta, tracks[d].Eta());
					double x4D_rot[4]={del_phi_rot,del_eta_rot,(double)trkbin,(double)multcentbin}; histo_rot->Fill(x4D_rot,jet_weight*trk_weight*event_weight*trkpt);
				}
			}
		}//end rotation
	}//end jet track loop
}

/*
Function to fill vectors to be used during the mixing (the ones with &)
--> Arguments
similar_events: true for using tracks only if event has one jet within the jet cut or falt for all tracks
nev: histogram with number of events stored to be used in the mixing
jets: vector with jet informations
jet_weight: vector with jet weight informations
tracks: vector with track informations
trk_weight: vector with track weight informations
mult: multiplicity or centrality
vertexz: Z vertex in centimeters
weight: event weight
ev_jet_vector: vector to be used in the mixing with jet information for each event
ev_jet_weight_vector: vector to be used in the mixing with jet weight information for each event
ev_track_vector: vector to be used in the mixing with track information for each event
ev_trk_weight_vector: vector to be used in the mixing with track weight information for each event
multvec: vector to be used in the mixing with event multiplicity or centrality information
vzvec: vector to be used in the mixing with event Z vertex position information
weightvec: vector to be used in the mixing with event weight information
*/






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
  /*Double_t kT_norm = TMath::Hypot((p1.Px()+p2.Px())/2.0,(p1.Py()+p2.Py())/2.0);                                                
   Double_t qT_dot_kT = (p1.Px()-p2.Px())*((p1.Px()+p2.Px())/2) + (p1.Py()-p2.Py())*((p1.Py()+p2.Py())/2);                        
                                                                                                                                  
   Double_t qout = 0.0;                                                                                                           
   if(kT_norm != 0) qout = qT_dot_kT/kT_norm;                                                                                     
   return qout;                                                                                                                   
  */
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

  /*Double_t qT_x = p1.Px()-p2.Px();                                                                                             
   Double_t qT_y = p1.Py()-p2.Py();                                                                                               
                                                                                                                                  
   Double_t kT_norm = TMath::Hypot((p1.Px()+p2.Px())/2.0,(p1.Py()+p2.Py())/2.0);                                                  
   Double_t qT_dot_kT = (p1.Px()-p2.Px())*((p1.Px()+p2.Px())/2) + (p1.Py()-p2.Py())*((p1.Py()+p2.Py())/2);                        
   Double_t kT_norm_x = ((p1.Px()+p2.Px())/2.0)/kT_norm;                                                                          
   Double_t qOut_x = (qT_dot_kT*kT_norm_x)/kT_norm;                                                                               
   Double_t kT_norm_y = ((p1.Py()+p2.Py())/2.0)/kT_norm;                                                                          
   Double_t qOut_y = (qT_dot_kT*kT_norm_y)/kT_norm;                                                                               
                                                                                                                                  
   Double_t qSide_x = qT_x - qOut_x;                                                                                              
   Double_t qSide_y = qT_y - qOut_y;                                                                                              
   Double_t qSide = TMath::Hypot(qSide_x,qSide_y);                                                                                
                                                                                                                                  
   if (qSide_y<0)qSide=-1.*qSide;                                                                                                 
                                                                                                                                  
   return qSide;                                                                                                                  
  */

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

/*
double getTrkCorrWeight(double pT, double eta){

  double eff = reff2D->GetBinContent(
				     reff2D->GetXaxis()->FindBin(eta),
				     reff2D->GetYaxis()->FindBin(pT) );
  if(eff >= 0.9999 || eff <= 0.0001) eff = 1;

  double sec = rsec2D->GetBinContent(
				     rsec2D->GetXaxis()->FindBin(eta),
				     rsec2D->GetYaxis()->FindBin(pT));
  if( sec >= 0.9999 || sec <= 0.0001) sec = 0;

  double fak = rfak2D->GetBinContent(
				     rfak2D->GetXaxis()->FindBin(eta),
				     rfak2D->GetYaxis()->FindBin(pT));
  if( fak >= 0.9999 || fak <= 0.0001) fak = 0;

  double mul = rmul2D->GetBinContent(
				     rmul2D->GetXaxis()->FindBin(eta),
				     rmul2D->GetYaxis()->FindBin(pT));
  if( mul >= 0.9999 || mul <= 0.0001) mul = 0;


  return (1. - fak ) * ( 1. - sec ) / eff  / (1. + mul );
  //return (1. - fak ) / eff;                                                                                                     
  //return 1. / eff;                                                                                                              

}
*/

Double_t Encode( int w )
{
        //int Range[]={3,5,8,10,13,16,20,25,30,200} ;
        int Range[]={3,5,8,10,13,16,20,25,30,200,10000} ;//added a protection for high-multiplicity 
        int i(0), j(0) ;
        while( w >= Range[i++] ) j++ ;
        //henc->Fill(w, (Double_t) j) ;
        return (Long64_t) j ;
}

Double_t ComputeEventWeight(std::vector<TLorentzVector> GoodTrackFourVector) // compute eta of event
{
        int nFwd(0), nBkw(0), nCtr(0) ;
        for( unsigned int i=0 ; i< GoodTrackFourVector.size() ; i++ )
    {
                Double_t eta = GoodTrackFourVector[i].Eta() ;
                if( eta < -0.8 )
        { nBkw++ ;}
                else if( eta > 0.8 )
        { nFwd++ ;}
                else
        { nCtr++ ;}
    }


        Double_t ReturnValue = (100*Encode( nFwd ) + 10 * Encode( nCtr ) + Encode( nBkw )) ;
        return ReturnValue ;
}


bool etaMixSort(std::pair<Double_t,std::vector<TLorentzVector>> & a, std::pair<Double_t,std::vector<TLorentzVector>> & b){

   return a.first < b.first ? true : false ;

}

bool etaMixSort_ch(std::pair<Double_t,std::vector<int>> & a, std::pair<Double_t,std::vector<int>> & b){

   return a.first < b.first ? true : false ;

}

bool etaMixSort_ntrkoff(std::pair<Double_t,int> & a, std::pair<Double_t,int> & b){

   return a.first < b.first ? true : false ;

}
