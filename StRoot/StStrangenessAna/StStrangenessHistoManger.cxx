#include "StStrangenessHistoManger.h"
#include "StStrangenessCons.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

ClassImp(StStrangenessHistoManger)

    //-------------------------------------------------------------
StStrangenessHistoManger::StStrangenessHistoManger()
{
}

StStrangenessHistoManger::~StStrangenessHistoManger()
{
}
//-------------------------------------------------------------

void StStrangenessHistoManger::Init(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
{
    // flow analysis
    for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt bin | TODO: increase pt_bin to 8 GeV/c
    {
        for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
            {
                for(Int_t m = 0; m < Strangeness::Phi_Psi_total; m++) // phi-psi bin
                {
                    TString Mode[2] = {"SE","ME"};
                    TString HistName;
                    HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_%s_%s",i,j,l,m,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
                    h_mMass2_EP[i][j][l][m] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
                    h_mMass2_EP[i][j][l][m]->Sumw2();
                    HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_%s_%s",i,j,l,m,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
                    h_mMass3_EP[i][j][l][m] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
                    h_mMass3_EP[i][j][l][m]->Sumw2();

                    // subtract K0s
                }
            }
        }
    }

    // raw pt spectra | TODO: use finer pt_bin
    for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt bin
    {
        for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
            {
                TString Mode[2] = {"SE","ME"};
                TString HistName; 
                HistName = Form("Spec_pt_%d_Centrality_%d_EtaGap_%d_%s_%s",i,j,l,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
                h_mMass_Spec[i][j][l] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
                h_mMass_Spec[i][j][l]->Sumw2();

                // subtract K0s
            }
        }
    }

    // Yields
    for(Int_t j = 0; j < 9; j++) // centrality bin
    {
        for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
        {
            TString Mode[2] = {"SE","ME"};
            TString HistName;
            HistName = Form("Yields_Centrality_%d_EtaGap_%d_%s_%s",j,l,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
            h_mMass_Yields[j][l] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
            h_mMass_Yields[j][l]->Sumw2();

            // subtract K0s
        }
    }

    h_pt_y_mass_se = new TH3F("h_pt_y_mass_se", "", 15, 0.0, 3.0, 40, -2.0, 2.0, 100, 0.98, 1.05);
    h_pt_y_mass_me = new TH3F("h_pt_y_mass_me", "", 15, 0.0, 3.0, 40, -2.0, 2.0, 100, 0.98, 1.05);
    //for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
    for (int sub = 0; sub < 2; sub++) // event plane Psi histograms
    {
    	h_psi1_epd_ABCD_shifted_wt_sub[sub] = new TH1F(Form("h_psi1_epd_ABCD_shifted_wt_sub_%d",sub + 1), Form("#Psi_{1}^{EPD} %d w/ eta weighting distribution (shifted)",sub + 1 ) ,1024,-7.0,7.0);
	h_psi2_tpc_AB_shifted_sub[sub] = new TH1F(Form("h_psi2_tpc_AB_shifted_sub_%d",sub + 1), Form("#Psi_{2}^{TPC} %d distribution (shifted)", sub + 1) ,1024,-7.0,7.0);
    }
    h_psi2_tpc_AB_shifted_subs = new TH2F(Form("h_psi2_tpc_AB_shifted_subs"), "#Psi_{2} distribution (shifted) #Psi_{2}^{East}  vs. #Psi_{2}^{West} " ,35,-3.5,3.5,35,-3.5,3.5);


}
//-------------------------------------------------------------
void StStrangenessHistoManger::Fill_EP_QA_East(Float_t Psi2_east)
{
	h_psi2_tpc_AB_shifted_sub[0]->Fill(Psi2_east);
}
void StStrangenessHistoManger::Fill_EP_QA_West(Float_t Psi2_west)
{
	h_psi2_tpc_AB_shifted_sub[1]->Fill(Psi2_west);
}
void StStrangenessHistoManger::Fill_EPs_QA(Float_t Psi2_east, Float_t Psi2_west)
{
	h_psi2_tpc_AB_shifted_subs->Fill(Psi2_east,Psi2_west);
}

void StStrangenessHistoManger::Fill(Float_t pt, Int_t Cent9, Int_t eta_gap, Float_t phi_psi2, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t InvMass, Double_t reweight)
{
    if(Res2 > 0.0)
    {
        for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt_bin
        {
            if(pt > Strangeness::pt_low_phi[i] && pt < Strangeness::pt_up_phi[i])
            {
                for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
                {
                    if(Cent9 >= Strangeness::cent_low[j] && Cent9 <= Strangeness::cent_up[j])
                    {
                        for(Int_t psi_bin = 0; psi_bin < 3; psi_bin++)
                        {
                            if(phi_psi2 >= Strangeness::Psi2_low[psi_bin] && phi_psi2 < Strangeness::Psi2_up[psi_bin])
                            {
                                Float_t phi_psi2_final = phi_psi2 - (psi_bin-1)*2.0*TMath::Pi()/2.0;
                                for(Int_t m = 0; m < Strangeness::Phi_Psi_total; m++) // phi-psi2 bin
                                {
                                    if(TMath::Abs(phi_psi2_final) >= Strangeness::phi_Psi2_low[m] && TMath::Abs(phi_psi2_final) < Strangeness::phi_Psi2_up[m])
                                    {
                                        // flow
                                        h_mMass2_EP[i][j][eta_gap][m]->Fill(InvMass,(reweight/Res2));
                                        // raw pt spectra
                                        h_mMass_Spec[i][j][eta_gap]->Fill(InvMass,reweight);
                                        //		    cout << "m = " << m << endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if(Res3 > 0.0)
    {
        for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt_bin
        {
            if(pt > Strangeness::pt_low_phi[i] && pt < Strangeness::pt_up_phi[i])
            {
                for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
                {
                    if(Cent9 >= Strangeness::cent_low[j] && Cent9 <= Strangeness::cent_up[j])
                    {
                        for(Int_t psi_bin = 0; psi_bin < 5; psi_bin++)
                        {
                            if(phi_psi3 >= Strangeness::Psi3_low[psi_bin] && phi_psi3 < Strangeness::Psi3_up[psi_bin])
                            {
                                Float_t phi_psi3_final = phi_psi3 - (psi_bin-2)*2.0*TMath::Pi()/3.0;
                                for(Int_t m = 0; m < Strangeness::Phi_Psi_total; m++) // phi-psi3 bin
                                {
                                    if(TMath::Abs(phi_psi3_final) >= Strangeness::phi_Psi3_low[m] && TMath::Abs(phi_psi3_final) < Strangeness::phi_Psi3_up[m])
                                    {
                                        // flow
                                        h_mMass3_EP[i][j][eta_gap][m]->Fill(InvMass,(reweight/Res3));
                                        //		    cout << "phi_psi3_bin = " << m << endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    h_mMass_Yields[Cent9][eta_gap]->Fill(InvMass,reweight);
}

void StStrangenessHistoManger::FillAcc(Float_t pt, Float_t rapidity, Float_t InvMass)
{
    h_pt_y_mass_se->Fill(pt, rapidity, InvMass);
    h_pt_y_mass_me->Fill(pt, rapidity, InvMass);
}

//-------------------------------------------------------------
void StStrangenessHistoManger::Write()
{
    // flow
    for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt bin
    {
        for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
            {
                for(Int_t m = 0; m < Strangeness::Phi_Psi_total; m ++) // phi-psi bin
                {
                    h_mMass2_EP[i][j][l][m]->Write();
                    h_mMass3_EP[i][j][l][m]->Write();
                }
            }
        }
    }

    // Yields
    for(Int_t j = 0; j < 9; j++) // centrality bin
    {
        for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
        {
            h_mMass_Yields[j][l]->Write();
        }
    }

    // raw pt spectra
    for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt bin
    {
        for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
        {
            for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
            {
                h_mMass_Spec[i][j][l]->Write();
            }
        }
    }
    h_pt_y_mass_se->Write();
    h_pt_y_mass_me->Write();
    h_psi2_tpc_AB_shifted_sub[0]->Write();
    h_psi2_tpc_AB_shifted_sub[1]->Write();
    h_psi2_tpc_AB_shifted_subs->Write();
}
