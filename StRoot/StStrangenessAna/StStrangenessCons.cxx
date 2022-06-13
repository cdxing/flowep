#include "StStrangenessCons.h"
#include "TMath.h"

ClassImp(Strangeness)

//-----------------------------------------------------------------
TString Strangeness::v0_tree[4]  = {"XuPhiMesonEvent","LambdaEvent","LambdaEvent","K0sEvent"}; // TODO: 2: LambdaEvent->antiLambdaEvent
TString Strangeness::v0_branch[4] = {"phi_flow_branch","Lambda_flow_branch","Lambda_flow_branch","K0s_flow_branch"}; // TODO: 2: Lambda_flow_branch->antiLambda_flow_branch
TString Strangeness::Energy[5] = {"200GeV","54GeV","39GeV","7GeV","15GeV"};  // shaowei
TString Strangeness::Partype[4] = {"Phi","Lambda","AntiLambda","K0s"};
Int_t   Strangeness::mList_Delta = 10;  // shaowei change from 20

//--------------------------------------------------------------------
// Cuts 
Float_t Strangeness::mDcaEPMax[5] = {3.0,1.0,1.0,1.0,1.0}; // for event plane reconstruction: 3.0 for 200GeV, 1.0 for BES
Int_t   Strangeness::mTrackMin = 2;
Float_t Strangeness::mShiftOrder2[5] = {2.0, 4.0, 6.0, 8.0, 10.0};
Float_t Strangeness::mShiftOrder3[5] = {3.0, 6.0, 9.0, 12.0, 15.0};
Int_t   Strangeness::mHitsFitTPCMin = 15;
Int_t   Strangeness::mHitsMaxTPCMin = 0;
Float_t Strangeness::mHitsRatioTPCMin = 0.52;
Float_t Strangeness::mPrimPtMin[5] = {0.15,0.2,0.2,0.2,0.2}; // for event plane reconstruction and for pion, kaon, proton: 0.15 for 200 GeV, 0.2 for BES
Float_t Strangeness::mPrimPtMax = 2.0;
Float_t Strangeness::mPrimPtWeight = 2.0;
Float_t Strangeness::mPrimMomMax = 10.0;
Float_t Strangeness::mEtaMax = 1.0;
Float_t Strangeness::mEta_Gap[4] = {0.05,0.10,0.20,0.50};
Int_t   Strangeness::mEtaGap_total = 4;

Float_t Strangeness::mBeamEnergy[5] = {200.0,54.0,39.0,27.0,14.5};   // shaowei
Int_t   Strangeness::mBeamYear[4] = {2011,2010,2011,2014};

//--------------------------------------------------------------------
// Histogram
// pt bin
//                                       0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,10 ,21 ,22
Float_t Strangeness::pt_low_phi[23] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2};
Float_t Strangeness::pt_up_phi[23]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6};

// Centrality bin
Int_t Strangeness::cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
Int_t Strangeness::cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
TString Strangeness::Centrality_01[4] = {"0080","0010","1040","4080"};
TString Strangeness::Centrality_23[4] = {"0070","0010","1040","4070"};

// phi-Psi bin
Double_t Strangeness::phi_Psi2_low[7] = {0.0,TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0};
Double_t Strangeness::phi_Psi2_up[7]  = {TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0,7.0*TMath::Pi()/14.0};
Double_t Strangeness::phi_Psi3_low[7] = {0.0,TMath::Pi()/21.0,2.0*TMath::Pi()/21.0,3.0*TMath::Pi()/21.0,4.0*TMath::Pi()/21.0,5.0*TMath::Pi()/21.0,6.0*TMath::Pi()/21.0};
Double_t Strangeness::phi_Psi3_up[7]  = {TMath::Pi()/21.0,2.0*TMath::Pi()/21.0,3.0*TMath::Pi()/21.0,4.0*TMath::Pi()/21.0,5.0*TMath::Pi()/21.0,6.0*TMath::Pi()/21.0,7.0*TMath::Pi()/21.0};

Double_t Strangeness::Psi2_low[3] = {-3.0*TMath::Pi()/2.0,-1.0*TMath::Pi()/2.0,1.0*TMath::Pi()/2.0};
Double_t Strangeness::Psi2_up[3]  = {-1.0*TMath::Pi()/2.0, 1.0*TMath::Pi()/2.0,3.0*TMath::Pi()/2.0};
Double_t Strangeness::Psi3_low[5] = {-4.0*TMath::Pi()/3.0,-3.0*TMath::Pi()/3.0,-1.0*TMath::Pi()/3.0,1.0*TMath::Pi()/3.0,3.0*TMath::Pi()/3.0};
Double_t Strangeness::Psi3_up[5]  = {-3.0*TMath::Pi()/3.0,-1.0*TMath::Pi()/3.0, 1.0*TMath::Pi()/3.0,3.0*TMath::Pi()/3.0,4.0*TMath::Pi()/3.0};

Int_t Strangeness::pt_total_phi = 23;

Int_t Strangeness::Centrality_total = 4;    // shaowei
Int_t Strangeness::Centrality_start = 0;
Int_t Strangeness::Centrality_stop  = 4;    // shaowei

Int_t Strangeness::EtaGap_total = 4;
Int_t Strangeness::EtaGap_start = 0;
Int_t Strangeness::EtaGap_stop  = 4;

Int_t Strangeness::Phi_Psi_total = 7;

Float_t Strangeness::InvMass_low[4] = {0.98,1.06,1.06,0.4};
Float_t Strangeness::InvMass_high[4] = {1.05,1.20,1.20,0.4};
