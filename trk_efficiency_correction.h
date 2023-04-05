#include "call_libraries.h"  // call libraries from ROOT and C++

// Obtain weight from tracking efficiency correction - default used by CMS
// it returns a efficiency correction to be applied in histograms as -> Fill(variable,this_weight)

bool checkBounds(double pt, double eta){
  if( TMath::Abs(eta) > 2.4 ){return false;}
  if( pt < 0 || pt > 500 ){return false;}
  return true;
}

// inputs are 2D histograms: reff2D for efficiency, rfak2D for fakes, rsec2D for secondary (decays), rmul2D for multiple reconstruction
// pT and eta are the transverse momentum and pseudorapidity of the track (considering a 2D histogram where X is eta axis and Y pT axis)
double getTrkCorrWeight(TH2 *reff2D, TH2 *rfak2D, TH2 *rsec2D, TH2 *rmul2D, double pT, double eta){
  if( !checkBounds(pT, eta) ) return 0;
  double factor = 1.0;
  double eff = reff2D->GetBinContent(
                  reff2D->GetXaxis()->FindBin(eta),
                  reff2D->GetYaxis()->FindBin(pT) );
  if(eff >= 0.9999 || eff <= 0.0001) eff = 1;

  double fak = rfak2D->GetBinContent(
              rfak2D->GetXaxis()->FindBin(eta),
              rfak2D->GetYaxis()->FindBin(pT));
  if( fak >= 0.9999 || fak <= 0.0001) fak = 0;

  double sec = rsec2D->GetBinContent(
              rsec2D->GetXaxis()->FindBin(eta),
              rsec2D->GetYaxis()->FindBin(pT));
  if( sec >= 0.9999 || sec <= 0.0001) sec = 0;

  double mul = rmul2D->GetBinContent(
              rmul2D->GetXaxis()->FindBin(eta),
              rmul2D->GetYaxis()->FindBin(pT));
  if( mul >= 0.9999 || mul <= 0.0001) mul = 0;

  // uncomment as needed
  // return (1. - fak ) * ( 1. - sec ) / eff  / (1. + mul ); //complete correction
  return (1. - fak ) / eff; //only fake and eff
  // return 1. / eff; //only efficiency

}
