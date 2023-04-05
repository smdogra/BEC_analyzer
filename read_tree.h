#include "call_libraries.h"  // call libraries from ROOT and C++

// declare variables

// event quantities
float vertexz; // event z vertex
int hiBin; // event centrality (used if use_centrality = true in input_variables.h)

// trigger quantities
int jet_trigger_bit; // jet HLT path trigger used for analysis (jet_trigger variable in input_variables.h)

// reco jets
//int nref;         // number of jets
//float jteta[9999]; // jet eta
//float jtphi[9999]; // jet phi
//float rawpt[9999]; // jet pT without JEC
//float trackMax[9999]; // track maximum pT in a jet
// reco tracks
int ntrk;                 // number of track
float trkpt[9999];       // track pT
float trketa[9999];      // track eta
float trkphi[9999];      // track phi
float trkpterr[9999];    // track pT error (uncertainty)
float trkdcaxy[9999];    // track dxy impact parameter (transverse distance between primary vertex and collision - distance of closest approuch - DCA)
float trkdcaz[9999];     // track dz impact parameter (longitudinal distance between primary vertex and collision - distance of closest approuch - DCA)
float trkdcaxyerr[9999]; // track dxy error (uncertainty)
float trkdcazerr[9999];  // track dxy error (uncertainty)
float trkchi2[9999];     // track reconstruction chi2 of the fitting
float pfEcal[9999];      // particle flow energy deposit in ECAL
float pfHcal[9999];      // particle flow energy deposit in HCAL
float trkmva[9999];      // track mva for each step
int trkalgo[9999];       // track algorithm/step
int trkndof[9999];       // track number of degrees of freedom in the fitting 
int trkcharge[9999];     // track charge
int trknhits[9999];      // number of hits in the tracker
int trknlayer[9999];     // number of layers with measurement in the tracker
bool highpur[9999];      // tracker steps MVA selection
//unsigned char trkpixhits[9999];  // 
UChar_t trkpixhits[9999];  //
// events quantities from gen
float weight; // event weight --> pthat weight
float pthat;  // pthat (initial parton pT)
/*-------
// gen jets
int ngen;             // number of gen jets
float gen_jtpt[9999];  // gen jet pT
float gen_jteta[9999]; // gen jet eta
float gen_jtphi[9999]; // gen jet phi

// matched jets
float refpt[9999]; // jet pT matched with Gen pT
float refeta[9999]; // jet eta matched with Gen eta
float refphi[9999]; // jet phi matched with Gen phi
int refparton_flavor[9999]; // jet phi matched with Gen phi
int refparton_flavorForB[9999]; // jet phi matched with Gen phi

// gen tracks
std::vector<float> *gen_trkpt = 0;  // gen particle pT
std::vector<float> *gen_trketa = 0; // gen particle eta
std::vector<float> *gen_trkphi = 0; // gen particle phi
std::vector<int> *gen_trkchg = 0;   // gen particle charge
std::vector<int> *gen_trkpid = 0;   // gen particle pid
std::vector<int> *gen_trksube = 0;   // gen particle pid
-------*/
//All variables listed above are readed in the function bellow
/*
Function to read the Forest/Skim tree
Arguments ->  transfer quantities from trees to our variables
tree: input TChain from jet_analyzer.C file
all the arguments bellow are defined in input_variables.h
is_MC: true -> MC; false -> Data
use_WTA: true -> use WTA (winner-takes-all); false -> use E-Scheme
jet_trigger: string with trigger name
colliding_system: pp, pPb, PbPb, XeXe, ... (future)
colliding_energy: colliding energy in GeV -> 5020, 5440, 8160, 13000, ... 
year_of_datataking: year of data taking
event_filterstr: string of event filters
event_filters: integer (0 or 1) from event filters
*/
void read_tree(TChain *tree, bool is_MC, TString colliding_system, int colliding_energy, int year_of_datataking, std::vector<TString> event_filterstr, std::vector<int> event_filters){

    tree->SetBranchStatus("*", 0); // disable all branches - this is important while reading big files

    // enable branches of interest -> see definition of each variables above

    // event quantities
    //tree->SetBranchStatus(Form("%s",jet_trigger.Data()), 1);
    tree->SetBranchStatus("vz", 1);
    if(colliding_system=="PbPb" || colliding_system=="XeXe") tree->SetBranchStatus("hiBin", 1); //centrality only for PbPb and XeXe
    for(int i = 0; i < event_filterstr.size(); i++) tree->SetBranchStatus(Form("%s",event_filterstr[i].Data()), 1); //event filters
    
    //tree->SetBranchAddress(Form("%s",jet_trigger.Data()), &jet_trigger_bit);
    tree->SetBranchAddress("vz", &vertexz);
    if(colliding_system=="PbPb" || colliding_system=="XeXe") tree->SetBranchAddress("hiBin", &hiBin); //centrality only for PbPb and XeXe
    for(int i = 0; i < event_filters.size(); i++) tree->SetBranchAddress(Form("%s",event_filterstr[i].Data()),&event_filters[i]);

    // track quantities
    tree->SetBranchStatus("nTrk", 1);
    tree->SetBranchStatus("trkPt", 1);
    tree->SetBranchStatus("trkEta", 1);
    tree->SetBranchStatus("trkPhi", 1);
    tree->SetBranchStatus("trkPtError", 1);
    tree->SetBranchStatus("trkDxy1", 1);
    tree->SetBranchStatus("trkDxyError1", 1);
    tree->SetBranchStatus("trkDz1", 1);
    tree->SetBranchStatus("trkDzError1", 1);
    tree->SetBranchStatus("trkNPixelHit", 1);

    if(colliding_system=="PbPb" || colliding_system=="XeXe"){
      tree->SetBranchStatus("trkChi2", 1);
      tree->SetBranchStatus("trkNdof", 1);
      tree->SetBranchStatus("trkNHit", 1);
      tree->SetBranchStatus("trkNlayer", 1);
    }
    tree->SetBranchStatus("trkCharge", 1);
    tree->SetBranchStatus("highPurity", 1);
    tree->SetBranchStatus("pfEcal", 1);
    tree->SetBranchStatus("pfHcal", 1);
    
    tree->SetBranchAddress("nTrk", &ntrk);
    tree->SetBranchAddress("trkPt", &trkpt);
    tree->SetBranchAddress("trkEta", &trketa);
    tree->SetBranchAddress("trkPhi", &trkphi);
    tree->SetBranchAddress("trkPtError", &trkpterr);
    tree->SetBranchAddress("trkNPixelHit", &trkpixhits);

    tree->SetBranchAddress("trkDxy1", &trkdcaxy);
    tree->SetBranchAddress("trkDxyError1", &trkdcaxyerr);
    tree->SetBranchAddress("trkDz1", &trkdcaz);
    tree->SetBranchAddress("trkDzError1", &trkdcazerr);
    tree->SetBranchAddress("trkCharge", &trkcharge);
    if(colliding_system=="PbPb" || colliding_system=="XeXe"){
      tree->SetBranchAddress("trkChi2", &trkchi2);
      tree->SetBranchAddress("trkNdof", &trkndof);
      tree->SetBranchAddress("trkNHit", &trknhits);
      tree->SetBranchAddress("trkNlayer", &trknlayer);
    }
    tree->SetBranchAddress("highPurity", &highpur);
    tree->SetBranchAddress("pfEcal", &pfEcal);
    tree->SetBranchAddress("pfHcal", &pfHcal);
    
    if(colliding_system=="PbPb" && colliding_energy==5020 && year_of_datataking==2018){ //special for 2018 PbPb MC
      tree->SetBranchStatus("trkMVA", 1);
      tree->SetBranchStatus("trkAlgo", 1);
      tree->SetBranchAddress("trkMVA", &trkmva);
      tree->SetBranchAddress("trkAlgo", &trkalgo);
    }
    /*--Sunil
    // gen particle quantities
    if(is_MC){
        tree->SetBranchStatus("pt", 1);
        tree->SetBranchStatus("eta", 1);
        tree->SetBranchStatus("phi", 1);
        tree->SetBranchStatus("chg", 1);
        tree->SetBranchStatus("pdg", 1);
        tree->SetBranchStatus("sube", 1);

        tree->SetBranchAddress("pt", &gen_trkpt);
        tree->SetBranchAddress("eta", &gen_trketa);
        tree->SetBranchAddress("phi", &gen_trkphi);
        tree->SetBranchAddress("chg", &gen_trkchg);
        tree->SetBranchAddress("pdg", &gen_trkpid);
        tree->SetBranchAddress("sube", &gen_trksube);
    }
    --*/
}
