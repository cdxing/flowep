#ifndef StStrangenessCons_h
#define StStrangenessCons_h

#include "TString.h"

class Strangeness
{
  public:
    static TString v0_tree[4]; // XuPhiMesonEvent, Lambda_flow_tree, ,antiLambda_flow_tree, K0s_flow_tree
    static TString v0_branch[4]; // Events, Lambda_flow_branch, antiLambda_flow_branch, K0s_flow_branch
    static TString Energy[5]; // 0 for 200GeV, 2 for 39GeV, 3 for 27GeV  // shaowei  1 for 54GeV
    static TString Partype[4]; // 0 for Phi, 1 for Lambda, 2 for antiLambda, 3 for K0s 
    static Int_t   mList_Delta; // number of list in one job
    //------------------------------------------------------------------
    // Cuts
    static Float_t mDcaEPMax[5]; // 0: 200 GeV, 1: 39 GeV, 2: 27 GeV, 3: 19.6 GeV, 4: 62.4 GeV
    static Int_t   mTrackMin; // 2
    static Float_t mShiftOrder2[5];
    static Float_t mShiftOrder3[5];
    static Int_t   mHitsFitTPCMin; // 15
    static Int_t   mHitsMaxTPCMin; // 0
    static Float_t mHitsRatioTPCMin; // 0.51
    static Float_t mPrimPtMin[5]; // 0: 0.15(200 GeV), 1-4: 0.2 (BES)
    static Float_t mPrimPtMax; // 2.0
    static Float_t mPrimPtWeight; // 2.0
    static Float_t mPrimMomMax;  // 10.0
    static Float_t mEtaMax; // 1.0
    static Float_t mEta_Gap[4]; // 0 = 0.05, 1 = 0.1, 2 = 0.2, 3 = 0.5
    static Int_t   mEtaGap_total; // 4
    static Float_t mBeamEnergy[5]; // 200GeV, 39GeV, 15GeV, 27GeV // shaowei
    static Int_t   mBeamYear[4]; 
    //--------------------------------------------------------------------
    // Histogram
    // pt bin
    static Float_t pt_low_phi[23];
    static Float_t pt_up_phi[23];

    // ycm bin
    static Float_t ycm_low_phi[16];
    static Float_t ycm_up_phi[16];

    // Centrality bin
    static Int_t   cent_low[4];   // shaowei
    static Int_t   cent_up[4];
    static TString Centrality_01[4]; // {"0080","0010","1040","4080"};
    static TString Centrality_23[4]; // {"0070","0010","1040","4070"};

    // phi-Psi bin
    static Double_t phi_Psi1_low[7];
    static Double_t phi_Psi1_up[7];

    static Double_t phi_Psi2_low[7];
    static Double_t phi_Psi2_up[7];
    static Double_t phi_Psi3_low[7];
    static Double_t phi_Psi3_up[7];

    static Double_t Psi1_low[3];
    static Double_t Psi1_up[3];
    static Double_t Psi2_low[3];
    static Double_t Psi2_up[3];
    static Double_t Psi3_low[5];
    static Double_t Psi3_up[5];

    static Int_t pt_total_phi;
    static Int_t ycm_total_phi;

    static Int_t Centrality_total;
    static Int_t Centrality_start;
    static Int_t Centrality_stop;

    static Int_t EtaGap_total;
    static Int_t EtaGap_start;
    static Int_t EtaGap_stop;

    static Int_t Phi_Psi_total;
    static Float_t InvMass_low[4];
    static Float_t InvMass_high[4];

  ClassDef(Strangeness,1)
};
#endif
