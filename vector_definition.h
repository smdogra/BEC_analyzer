#include "call_libraries.h"  // call libraries from ROOT and C++

//Define the vectors used in the mixing events

// reco jet  + reco track
std::vector<std::vector<TVector3>> ev_track_vector_reco_reco;
std::vector<std::vector<TVector3>> ev_jet_vector_reco_reco;
std::vector<int> multvec_reco_reco;
std::vector<double> vzvec_reco_reco;
std::vector<double> weights_reco_reco;
std::vector<std::vector<double>> jet_weights_reco_reco;
std::vector<std::vector<double>> trk_weights_reco_reco;

// reco leading jet  + reco track
std::vector<std::vector<TVector3>> ev_track_vector_leadjet_reco_reco;
std::vector<std::vector<TVector3>> ev_jet_vector_leadjet_reco_reco;
std::vector<int> multvec_leadjet_reco_reco;
std::vector<double> vzvec_leadjet_reco_reco;
std::vector<double> weights_leadjet_reco_reco;
std::vector<std::vector<double>> jet_weights_leadjet_reco_reco;
std::vector<std::vector<double>> trk_weights_leadjet_reco_reco;

// reco subleading jet  + reco track
std::vector<std::vector<TVector3>> ev_track_vector_subleadjet_reco_reco;
std::vector<std::vector<TVector3>> ev_jet_vector_subleadjet_reco_reco;
std::vector<int> multvec_subleadjet_reco_reco;
std::vector<double> vzvec_subleadjet_reco_reco;
std::vector<double> weights_subleadjet_reco_reco;
std::vector<std::vector<double>> jet_weights_subleadjet_reco_reco;
std::vector<std::vector<double>> trk_weights_subleadjet_reco_reco;

// reco jet  + gen track
std::vector<std::vector<TVector3>> ev_track_vector_reco_gen;
std::vector<std::vector<TVector3>> ev_jet_vector_reco_gen;
std::vector<int> multvec_reco_gen;
std::vector<double> vzvec_reco_gen;
std::vector<double> weights_reco_gen;
std::vector<std::vector<double>> jet_weights_reco_gen;
std::vector<std::vector<double>> trk_weights_reco_gen;

// reco leading jet  + gen track
std::vector<std::vector<TVector3>> ev_track_vector_leadjet_reco_gen;
std::vector<std::vector<TVector3>> ev_jet_vector_leadjet_reco_gen;
std::vector<int> multvec_leadjet_reco_gen;
std::vector<double> vzvec_leadjet_reco_gen;
std::vector<double> weights_leadjet_reco_gen;
std::vector<std::vector<double>> jet_weights_leadjet_reco_gen;
std::vector<std::vector<double>> trk_weights_leadjet_reco_gen;

// reco subleading jet  + gen track
std::vector<std::vector<TVector3>> ev_track_vector_subleadjet_reco_gen;
std::vector<std::vector<TVector3>> ev_jet_vector_subleadjet_reco_gen;
std::vector<int> multvec_subleadjet_reco_gen;
std::vector<double> vzvec_subleadjet_reco_gen;
std::vector<double> weights_subleadjet_reco_gen;
std::vector<std::vector<double>> jet_weights_subleadjet_reco_gen;
std::vector<std::vector<double>> trk_weights_subleadjet_reco_gen;

// gen jet  + reco track
std::vector<std::vector<TVector3>> ev_track_vector_gen_reco;
std::vector<std::vector<TVector3>> ev_jet_vector_gen_reco;
std::vector<int> multvec_gen_reco;
std::vector<double> vzvec_gen_reco;
std::vector<double> weights_gen_reco;
std::vector<std::vector<double>> jet_weights_gen_reco;
std::vector<std::vector<double>> trk_weights_gen_reco;

// gen leading jet  + reco track
std::vector<std::vector<TVector3>> ev_track_vector_leadjet_gen_reco;
std::vector<std::vector<TVector3>> ev_jet_vector_leadjet_gen_reco;
std::vector<int> multvec_leadjet_gen_reco;
std::vector<double> vzvec_leadjet_gen_reco;
std::vector<double> weights_leadjet_gen_reco;
std::vector<std::vector<double>> jet_weights_leadjet_gen_reco;
std::vector<std::vector<double>> trk_weights_leadjet_gen_reco;

// gen subleading jet  + reco track
std::vector<std::vector<TVector3>> ev_track_vector_subleadjet_gen_reco;
std::vector<std::vector<TVector3>> ev_jet_vector_subleadjet_gen_reco;
std::vector<int> multvec_subleadjet_gen_reco;
std::vector<double> vzvec_subleadjet_gen_reco;
std::vector<double> weights_subleadjet_gen_reco;
std::vector<std::vector<double>> jet_weights_subleadjet_gen_reco;
std::vector<std::vector<double>> trk_weights_subleadjet_gen_reco;

// gen jet  + gen track
std::vector<std::vector<TVector3>> ev_track_vector_gen_gen;
std::vector<std::vector<TVector3>> ev_jet_vector_gen_gen;
std::vector<int> multvec_gen_gen;
std::vector<double> vzvec_gen_gen;
std::vector<double> weights_gen_gen;
std::vector<std::vector<double>> jet_weights_gen_gen;
std::vector<std::vector<double>> trk_weights_gen_gen;

// gen leading jet  + reco gen
std::vector<std::vector<TVector3>> ev_track_vector_leadjet_gen_gen;
std::vector<std::vector<TVector3>> ev_jet_vector_leadjet_gen_gen;
std::vector<int> multvec_leadjet_gen_gen;
std::vector<double> vzvec_leadjet_gen_gen;
std::vector<double> weights_leadjet_gen_gen;
std::vector<std::vector<double>> jet_weights_leadjet_gen_gen;
std::vector<std::vector<double>> trk_weights_leadjet_gen_gen;

// gen subleading jet  + reco track
std::vector<std::vector<TVector3>> ev_track_vector_subleadjet_gen_gen;
std::vector<std::vector<TVector3>> ev_jet_vector_subleadjet_gen_gen;
std::vector<int> multvec_subleadjet_gen_gen;
std::vector<double> vzvec_subleadjet_gen_gen;
std::vector<double> weights_subleadjet_gen_gen;
std::vector<std::vector<double>> jet_weights_subleadjet_gen_gen;
std::vector<std::vector<double>> trk_weights_subleadjet_gen_gen;
