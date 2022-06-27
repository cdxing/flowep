#ifndef StStrangenessHistoManger_h
#define StStrangenessHistoManger_h

#include "StMessMgr.h"

class TH1F;
class TH2F;
class TH3F;

class StStrangenessHistoManger
{
  public:
    StStrangenessHistoManger();
    ~StStrangenessHistoManger();

    void Init(Int_t X_flag, Int_t mode);
    void Fill(Float_t ycm, Int_t Cent9, Int_t eta_gap, Float_t phi_psi1, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t Mass2, Double_t reweight);
    void FillPsiShiftFull_wt(Float_t Psi);
    void Fill_EP_QA_East(Float_t Psi2_east);
    void Fill_EP_QA_West(Float_t Psi2_west);
    void Fill_EPs_QA(Float_t Psi2_east, Float_t Psi2_west);

    void Fill_sub(Float_t pt, Int_t Cent9, Int_t eta_gap, Float_t phi_psi1, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t Mass2, Double_t reweight);
    void Write();
    void FillAcc(Float_t, Float_t, Float_t);

  private:
    // flow analysis
    // 0 = pt bin
    // 1 = centrality: 0 = 0-80%(0-70%), 1 = 0-10%, 2 = 10-40%, 3 = 40-80%(40-70%)
    // 2 = eta_gap
    // 3 = phi - Psi
    TH1F *h_mMass1_EP[16][4][4][7];   // dchen

    TH1F *h_mMass2_EP[23][4][4][7];   // shaowei
    TH1F *h_mMass3_EP[23][4][4][7];  // shaowei
    // subtract k0s
    TH1F *h_mMass2_EP_sub[23][4][4][7];
    TH1F *h_mMass3_EP_sub[23][4][4][7];

    TH3F *h_pt_y_mass_se;
    TH3F *h_pt_y_mass_me;

    // raw pt spectra
    // 0 = pt bin
    // 1 = centrality: 0 = 0-80%(0-70%), 1 = 0-10%, 2 = 10-40%, 3 = 40-80%(40-70%)
    // 2 = eta_gap
    TH1F *h_mMass_Spec[16][4][4];
    //TH1F *h_mMass_Spec[23][4][4];
    TH1F *h_mMass_Spec_sub[23][4][4];

    // event plane resolution correction
    // 0 = centrality
    // 1 = eta_gap
    TH1F *h_mMass_Yields[9][4];
    // subtract k0s
    TH1F *h_mMass_Yields_sub[9][4];
    // Event Plane QA
    //TH1F *h_psi1_epd_ABCD_shifted_wt_sub[2];
    TH1F *h_psi1_epd_ABCD_shifted_wt_full;
    //TH1F *h_psi1_epd_ABCD_shifted_wt_sub[_numSubEvents];
    TH2F *h_psi2_tpc_AB_shifted_subs;

    TH1F *h_psi2_tpc_AB_shifted_sub[2];
    //TH1F *h_psi2_tpc_AB_shifted_sub[_numSubEvents];

  ClassDef(StStrangenessHistoManger,1)
};
#endif
