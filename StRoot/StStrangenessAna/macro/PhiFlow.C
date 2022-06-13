#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLegend.h"

Double_t PolyBreitWigner(Double_t *x_val, Double_t *par) 
{
  Double_t x = x_val[0];
  Double_t m0 = par[0];
  Double_t Gamma = par[1];
  Double_t Norm = par[2];

  Double_t denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  Double_t BW = Norm*Gamma/denom;

  Double_t Poly = par[3] + par[4]*x;

  Double_t y = BW + Poly;

  return y;
}

Double_t Poly(Double_t *x_val, Double_t *par) 
{
  Double_t x = x_val[0];
  Double_t y = par[0] + par[1]*x;

  return y;
}

Double_t BreitWigner(Double_t *x_val, Double_t *par) 
{
  Double_t x = x_val[0];
  Double_t m0 = par[0];
  Double_t Gamma = par[1];
  Double_t Norm = par[2];

  Double_t denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  Double_t BW = Norm*Gamma/denom;

  return BW;
}

Double_t flow_2(Double_t *x_val, Double_t *par)
{
  Double_t x, y;
  Double_t Ampl, v2;
  x = x_val[0];
  Ampl = par[0];
  v2 = par[1];

  y = Ampl*(1.0 + 2.0*v2*TMath::Cos(2.0*x));

  return y;
}

Double_t flow_3(Double_t *x_val, Double_t *par)
{
  Double_t x, y;
  Double_t Ampl, v3;
  x = x_val[0];
  Ampl = par[0];
  v3 = par[1];

  y = Ampl*(1.0 + 2.0*v3*TMath::Cos(3.0*x));

  return y;
}

Double_t Gaussion(Double_t *x_val, Double_t *par)
{
  Double_t x, y;
  x = x_val[0];

  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t norm = par[2];

  y = norm*TMath::Exp(-1.0*(x-mu)*(x-mu)/(2.0*sigma*sigma))/(sigma*TMath::Sqrt(2*TMath::Pi()));

  return y;
}

TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1)
{
  // coordinates in NDC!
  // plots the string label in position x and y in NDC coordinates
  // size is the text size
  // color is the text color

  //    if(x<0||y<0)
  //    {   // defaults
  //      x=gPad->GetLeftMargin()*1.15;
  //      y=(1-gPad->GetTopMargin())*1.04;
  //    }
  TLatex* text=new TLatex(x,y,label);
  text->SetTextFont(font);
  text->SetTextSize(size);
  if(NDC == 1) text->SetNDC();
  text->SetTextColor(color);
  text->SetTextAngle(angle);
  text->Draw();
  return text;
}

void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
  TLine* Zero_line = new TLine();
  Zero_line -> SetX1(x1_val);
  Zero_line -> SetX2(x2_val);
  Zero_line -> SetY1(y1_val);
  Zero_line -> SetY2(y2_val);
  Zero_line -> SetLineWidth(LineWidth);
  Zero_line -> SetLineStyle(LineStyle);
  Zero_line -> SetLineColor(Line_Col);
  Zero_line -> Draw();
  //delete Zero_line;
}

static const TString Energy[3] = {"200GeV","20GeV","15GeV"};
static const Int_t pt_total_phi = 16;
static const Int_t Centrality_total = 1;
static const Int_t Centrality_start = 0;
static const Int_t Centrality_stop = 1;
static const Int_t EtaGap_total = 1;
static const Int_t EtaGap_start = 0;
static const Int_t EtaGap_stop = 1;
static const Int_t Phi_Psi_total = 7;
static const Float_t nSigmaPhi = 2.0;
static const Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};
static const Int_t cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
static const Int_t cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
static const Float_t BW_Start = 0.994;
static const Float_t BW_Stop  = 1.050;

/*
// pt bin
//                                           0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15
static const Int_t pt_total_New_phi = 15;
static const Float_t pt_low_phi[pt_total_phi] = {0.2,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
static const Float_t pt_up_phi[pt_total_phi]  = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};
static const Int_t pt_new_bin_start[pt_total_New_phi] = {0,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
static const Int_t pt_new_bin_stop[pt_total_New_phi]  = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
// pt bin
// //                                       0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,10 ,21 ,22
// Float_t Strangeness::pt_low_phi[23] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2};
// Float_t Strangeness::pt_up_phi[23]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6};
*/

/*
// pt bin
//                                           0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15
static const Int_t pt_total_New_phi = 6;
static Float_t pt_low_phi[pt_total_phi] = {0.2,0.8,1.2,1.6,2.0,2.8};
static Float_t pt_up_phi[pt_total_phi]  = {0.8,1.2,1.6,2.0,2.8,4.2};
static Int_t pt_new_bin_start[pt_total_New_phi] = {0,3,5,7, 9,13};
static Int_t pt_new_bin_stop[pt_total_New_phi]  = {2,4,6,8,12,15};
*/


static const Int_t pt_total_New_phi = 3;
static Float_t pt_low_phi[pt_total_phi] = {0.2,1.0,1.8};
static Float_t pt_up_phi[pt_total_phi]  = {1.0,1.8,2.6};
static Int_t pt_new_bin_start[pt_total_New_phi] = {0,4, 8};
static Int_t pt_new_bin_stop[pt_total_New_phi]  = {3,7,11};

/*
static const Int_t pt_total_New_phi = 8;
static Float_t pt_low_phi[pt_total_New_phi] = {0.4,0.6,1.0,1.2,1.6,2.0,2.4,2.8};
static Float_t pt_up_phi[pt_total_New_phi]  = {0.6,1.0,1.2,1.6,2.0,2.4,2.8,4.2};
static Int_t pt_new_bin_start[pt_total_New_phi] = {1,2,4,5,7, 9,11,13};
static Int_t pt_new_bin_stop[pt_total_New_phi]  = {1,3,4,6,8,10,12,16};
*/
/*
static const Int_t pt_total_New_phi = 3;
static Float_t pt_low_phi[pt_total_New_phi] = {0.2,0.8,1.6};
static Float_t pt_up_phi[pt_total_New_phi]  = {0.8,1.6,3.4};
static Int_t pt_new_bin_start[pt_total_New_phi] = {0,3,7};
static Int_t pt_new_bin_stop[pt_total_New_phi]  = {2,6,14};
*/

void PhiFlow(const Int_t energy = 2)
{
//  gStyle->SetOptFit(1111);
  gStyle->SetTitleX(0.55);
  // open input file for same event and mixed event
  TString inputfile_SE = Form("/project/projectdirs/star/sss/OutPut/AuAu15GeV/Phi/flow_Phi/Yields_SE_%s.root",Energy[energy].Data());
  TFile *File_SE = TFile::Open(inputfile_SE.Data());
  TString inputfile_ME = Form("/project/projectdirs/star/sss/OutPut/AuAu15GeV/Phi/flow_Phi/Yields_ME_%s.root",Energy[energy].Data());
  TFile *File_ME = TFile::Open(inputfile_ME.Data());

  // read histogram for same event and mixed event
  TH1F *h_mMass_Phi2_SE[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TH1F *h_mMass_Phi3_SE[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];

  TH1F *h_mMass_Phi2_ME[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TH1F *h_mMass_Phi3_ME[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];

  for(Int_t i = 0; i < pt_total_phi; i++) // pt bin
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  TString HistName;

	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_Phi_SE",i,j,l,m);
	  h_mMass_Phi2_SE[i][j][l][m] = (TH1F*)File_SE->Get(HistName.Data())->Clone();
	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_Phi_SE",i,j,l,m);
	  h_mMass_Phi3_SE[i][j][l][m] = (TH1F*)File_SE->Get(HistName.Data())->Clone();

	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_Phi_ME",i,j,l,m);
	  h_mMass_Phi2_ME[i][j][l][m] = (TH1F*)File_ME->Get(HistName.Data())->Clone();
	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_Phi_ME",i,j,l,m);
	  h_mMass_Phi3_ME[i][j][l][m] = (TH1F*)File_ME->Get(HistName.Data())->Clone();
	}
      }
    }
  }
  /*
  File_SE->Close();
  File_ME->Close();
  */

  // calculate scaling factor and subtract mixed event from same event
  TH1F *h_mMass_Phi2_SM[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TH1F *h_mMass_Phi3_SM[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  for(Int_t i = 0; i < pt_total_phi; i++) // pt bin
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  Int_t bin2_SE_start = h_mMass_Phi2_SE[i][j][l][m]->FindBin(1.04);
	  Int_t bin2_SE_stop  = h_mMass_Phi2_SE[i][j][l][m]->FindBin(1.05);
	  Float_t Inte2_SE = h_mMass_Phi2_SE[i][j][l][m]->Integral(bin2_SE_start,bin2_SE_stop);
	  h_mMass_Phi2_SM[i][j][l][m] = (TH1F*)h_mMass_Phi2_SE[i][j][l][m]->Clone();

	  Int_t bin3_SE_start = h_mMass_Phi3_SE[i][j][l][m]->FindBin(1.04);
	  Int_t bin3_SE_stop  = h_mMass_Phi3_SE[i][j][l][m]->FindBin(1.05);
	  Float_t Inte3_SE = h_mMass_Phi3_SE[i][j][l][m]->Integral(bin3_SE_start,bin3_SE_stop);
	  h_mMass_Phi3_SM[i][j][l][m] = (TH1F*)h_mMass_Phi3_SE[i][j][l][m]->Clone();

	  Int_t bin2_ME_start = h_mMass_Phi2_ME[i][j][l][m]->FindBin(1.04);
	  Int_t bin2_ME_stop  = h_mMass_Phi2_ME[i][j][l][m]->FindBin(1.05);
	  Float_t Inte2_ME = h_mMass_Phi2_ME[i][j][l][m]->Integral(bin2_ME_start,bin2_ME_stop);
	  h_mMass_Phi2_ME[i][j][l][m]->Scale(Inte2_SE/Inte2_ME);
	  h_mMass_Phi2_SM[i][j][l][m]->Add(h_mMass_Phi2_ME[i][j][l][m],-1.0);

	  Int_t bin3_ME_start = h_mMass_Phi3_ME[i][j][l][m]->FindBin(1.04);
	  Int_t bin3_ME_stop  = h_mMass_Phi3_ME[i][j][l][m]->FindBin(1.05);
	  Float_t Inte3_ME = h_mMass_Phi3_ME[i][j][l][m]->Integral(bin3_ME_start,bin3_ME_stop);
	  h_mMass_Phi3_ME[i][j][l][m]->Scale(Inte3_SE/Inte3_ME);
	  h_mMass_Phi3_SM[i][j][l][m]->Add(h_mMass_Phi3_ME[i][j][l][m],-1.0);
	}
      }
    }
  }

  // merged original pT bins to get new pT bins
  TH1F *h_mMass_Phi2_SM_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TH1F *h_mMass_Phi3_SM_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];

  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
      {
	for(Int_t i = 0; i < pt_total_New_phi; i++)
	{
	  for(Int_t pt_bin = pt_new_bin_start[i]; pt_bin <= pt_new_bin_stop[i]; pt_bin++)
	  {
	    if(pt_bin == pt_new_bin_start[i])
	    {
	      h_mMass_Phi2_SM_New_Bin[i][j][l][m] = (TH1F*)h_mMass_Phi2_SM[pt_bin][j][l][m]->Clone();
	      h_mMass_Phi3_SM_New_Bin[i][j][l][m] = (TH1F*)h_mMass_Phi3_SM[pt_bin][j][l][m]->Clone();
	    }
	    else
	    {
	      h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Add(h_mMass_Phi2_SM[pt_bin][j][l][m],1.0);
	      h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Add(h_mMass_Phi3_SM[pt_bin][j][l][m],1.0);
	    }
	  }
	  if(i >= 6)
	  {
	    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Rebin(2);
	    h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Rebin(2);
	  }
	}
      }
    }
  }

  /*
  // Draw new bins
  TCanvas *c_Phi2_New_Bin[Centrality_total][EtaGap_total][Phi_Psi_total];
  TCanvas *c_Phi3_New_Bin[Centrality_total][EtaGap_total][Phi_Psi_total];

  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
      {
	TString CanName;

	CanName = Form("c_Phi2_SM_New_Bin_Centrality_%d_EtaGap_%d_phi_Psi_%d",j,l,m);
	c_Phi2_New_Bin[j][l][m] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
	c_Phi2_New_Bin[j][l][m]->Divide(3,3);
	for(Int_t i = 0; i < 8; i++)
	{
	  c_Phi2_New_Bin[j][l][m]->cd(i+1);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetTitle("");
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetStats(0);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->DrawCopy("PE");
	}

	CanName = Form("c_Phi3_SM_New_Bin_Centrality_%d_EtaGap_%d_phi_Psi_%d",j,l,m);
	c_Phi3_New_Bin[j][l][m] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
	c_Phi3_New_Bin[j][l][m]->Divide(3,3);
	for(Int_t i = 0; i < 8; i++)
	{
	  c_Phi3_New_Bin[j][l][m]->cd(i+1);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetTitle("");
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetStats(0);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->DrawCopy("PE");
	}
      }
    }
  }
  */

  // Poly+BW Fit for merged phi-psi bin
  TH1F *h_mMass_Phi2_SM_total[pt_total_New_phi][Centrality_total][EtaGap_total];
  TH1F *h_mMass_Phi3_SM_total[pt_total_New_phi][Centrality_total][EtaGap_total];

  TF1  *f_PolyBW_Phi2_total[pt_total_New_phi][Centrality_total][EtaGap_total];
  TF1  *f_PolyBW_Phi3_total[pt_total_New_phi][Centrality_total][EtaGap_total];

  Float_t ParFit2_total[pt_total_New_phi][Centrality_total][EtaGap_total][5];
  Float_t ParFit3_total[pt_total_New_phi][Centrality_total][EtaGap_total][5];

  for(Int_t i = 0; i < pt_total_New_phi; i++) // pt bin
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  TString HistName;
	  if(m == 0)
	  {
	    HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_Total_2nd_Phi_SE",i,j,l);
	    h_mMass_Phi2_SM_total[i][j][l] = (TH1F*)h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Clone(HistName.Data());
	    HistName = "f_"+HistName;
	    f_PolyBW_Phi2_total[i][j][l] = new TF1(HistName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
	    for(Int_t n_par = 0; n_par < 5; n_par++)
	    {
	      f_PolyBW_Phi2_total[i][j][l]->ReleaseParameter(n_par);
	    }
	    f_PolyBW_Phi2_total[i][j][l]->SetParameter(0,1.019);
	    f_PolyBW_Phi2_total[i][j][l]->SetParLimits(0,1.014,1.024);
	    f_PolyBW_Phi2_total[i][j][l]->SetParameter(1,0.0055);
	    f_PolyBW_Phi2_total[i][j][l]->SetParameter(2,10000);
	    f_PolyBW_Phi2_total[i][j][l]->SetParameter(3,-6000);
	    f_PolyBW_Phi2_total[i][j][l]->SetParameter(4,0.5);
	    f_PolyBW_Phi2_total[i][j][l]->SetRange(BW_Start,BW_Stop);

	    HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_Total_3rd_Phi_SE",i,j,l);
	    h_mMass_Phi3_SM_total[i][j][l] = (TH1F*)h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Clone(HistName.Data());
	    HistName = "f_"+HistName;
	    f_PolyBW_Phi3_total[i][j][l] = new TF1(HistName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
	    for(Int_t n_par = 0; n_par < 4; n_par++)
	    {
	      f_PolyBW_Phi3_total[i][j][l]->ReleaseParameter(n_par);
	    }
	    f_PolyBW_Phi3_total[i][j][l]->SetParameter(0,1.019);
	    f_PolyBW_Phi3_total[i][j][l]->SetParLimits(0,1.014,1.024);
	    f_PolyBW_Phi3_total[i][j][l]->SetParameter(1,0.0055);
	    f_PolyBW_Phi3_total[i][j][l]->SetParameter(2,10000);
	    f_PolyBW_Phi3_total[i][j][l]->SetParameter(3,-6000);
	    f_PolyBW_Phi3_total[i][j][l]->SetParameter(4,0.5);
	    f_PolyBW_Phi3_total[i][j][l]->SetRange(BW_Start,BW_Stop);
	  }
	  else
	  {
	    h_mMass_Phi2_SM_total[i][j][l]->Add(h_mMass_Phi2_SM_New_Bin[i][j][l][m],1.0);
	    h_mMass_Phi3_SM_total[i][j][l]->Add(h_mMass_Phi3_SM_New_Bin[i][j][l][m],1.0);
	  }
	}
	h_mMass_Phi2_SM_total[i][j][l]->Fit(f_PolyBW_Phi2_total[i][j][l],"NQR");
	h_mMass_Phi3_SM_total[i][j][l]->Fit(f_PolyBW_Phi3_total[i][j][l],"NQR");
	for(Int_t n_par = 0; n_par < 5; n_par++)
	{
	  ParFit2_total[i][j][l][n_par] = f_PolyBW_Phi2_total[i][j][l]->GetParameter(n_par);
	  ParFit3_total[i][j][l][n_par] = f_PolyBW_Phi3_total[i][j][l]->GetParameter(n_par);
	}
      }
    }
  }

  /*
  // Draw the merged histogram and fit line
  TCanvas *c_mMass_Phi2_SM_total[Centrality_total][EtaGap_total];
  TCanvas *c_mMass_Phi3_SM_total[Centrality_total][EtaGap_total];

  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString CanName;

      CanName = Form("c_Phi2_SM_Total_Centrality_%d_EtaGap_%d",j,l);
      c_mMass_Phi2_SM_total[j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
      c_mMass_Phi2_SM_total[j][l]->Divide(3,3);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	c_mMass_Phi2_SM_total[j][l]->cd(i+1);
	h_mMass_Phi2_SM_total[i][j][l]->Draw("PE");
	f_PolyBW_Phi2_total[i][j][l]->SetLineColor(2);
	f_PolyBW_Phi2_total[i][j][l]->Draw("l same");
      }

      CanName = Form("c_Phi3_SM_Total_Centrality_%d_EtaGap_%d",j,l);
      c_mMass_Phi3_SM_total[j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
      c_mMass_Phi3_SM_total[j][l]->Divide(3,3);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	c_mMass_Phi3_SM_total[j][l]->cd(i+1);
	h_mMass_Phi3_SM_total[i][j][l]->Draw("PE");
	f_PolyBW_Phi3_total[i][j][l]->SetLineColor(2);
	f_PolyBW_Phi3_total[i][j][l]->Draw("l same");
      }
    }
  }
  */

  // Poly+BW Fit for phi-psi bin
  TF1 *f_PolyBW_Phi2[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TF1 *f_PolyBW_Phi3[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  for(Int_t i = 0; i < pt_total_New_phi; i++) // pt bin
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  TString HistName;

	  HistName = Form("f_pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_Phi_SE",i,j,l,m);
	  f_PolyBW_Phi2[i][j][l][m] = new TF1(HistName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
	  for(Int_t n_par = 0; n_par < 5; n_par++)
	  {
	    f_PolyBW_Phi2[i][j][l][m]->ReleaseParameter(n_par);
	  }
	  f_PolyBW_Phi2[i][j][l][m]->FixParameter(0,ParFit2_total[i][j][l][0]);
	  f_PolyBW_Phi2[i][j][l][m]->FixParameter(1,ParFit2_total[i][j][l][1]);
	  f_PolyBW_Phi2[i][j][l][m]->SetParameter(2,ParFit2_total[i][j][l][2]/7.0);
	  f_PolyBW_Phi2[i][j][l][m]->SetParameter(3,ParFit2_total[i][j][l][3]);
	  f_PolyBW_Phi2[i][j][l][m]->SetParameter(4,ParFit2_total[i][j][l][4]);
	  f_PolyBW_Phi2[i][j][l][m]->SetRange(BW_Start,BW_Stop);

	  HistName = Form("f_pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_Phi_SE",i,j,l,m);
	  f_PolyBW_Phi3[i][j][l][m] = new TF1(HistName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
	  for(Int_t n_par = 0; n_par < 5; n_par++)
	  {
	    f_PolyBW_Phi3[i][j][l][m]->ReleaseParameter(n_par);
	  }
	  f_PolyBW_Phi3[i][j][l][m]->FixParameter(0,ParFit3_total[i][j][l][0]);
	  f_PolyBW_Phi3[i][j][l][m]->FixParameter(1,ParFit3_total[i][j][l][1]);
	  f_PolyBW_Phi3[i][j][l][m]->SetParameter(2,ParFit3_total[i][j][l][2]/7.0);
	  f_PolyBW_Phi3[i][j][l][m]->SetParameter(3,ParFit3_total[i][j][l][3]);
	  f_PolyBW_Phi3[i][j][l][m]->SetParameter(4,ParFit3_total[i][j][l][4]);
	  f_PolyBW_Phi3[i][j][l][m]->SetRange(BW_Start,BW_Stop);

	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Fit(f_PolyBW_Phi2[i][j][l][m],"NQR");
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Fit(f_PolyBW_Phi3[i][j][l][m],"NQR");
	}
      }
    }
  }

  // subtract linear background from signal
  TF1 *f_Poly_Phi2[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TF1 *f_Poly_Phi3[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  TString FuncName;

	  FuncName = Form("f_poly2_pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d",i,j,l,m);
	  f_Poly_Phi2[i][j][l][m] = new TF1(FuncName.Data(),Poly,BW_Start,BW_Stop,2);
	  f_Poly_Phi2[i][j][l][m]->FixParameter(0,f_PolyBW_Phi2[i][j][l][m]->GetParameter(3));
	  f_Poly_Phi2[i][j][l][m]->FixParameter(1,f_PolyBW_Phi2[i][j][l][m]->GetParameter(4));
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Add(f_Poly_Phi2[i][j][l][m],-1.0);

	  FuncName = Form("f_poly3_pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d",i,j,l,m);
	  f_Poly_Phi3[i][j][l][m] = new TF1(FuncName.Data(),Poly,BW_Start,BW_Stop,2);
	  f_Poly_Phi3[i][j][l][m]->FixParameter(0,f_PolyBW_Phi3[i][j][l][m]->GetParameter(3));
	  f_Poly_Phi3[i][j][l][m]->FixParameter(1,f_PolyBW_Phi3[i][j][l][m]->GetParameter(4));
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Add(f_Poly_Phi3[i][j][l][m],-1.0);
	}
      }
    }
  }

  // do gaussian fit and extract width of Phi from merged phi-Psi distribution
  // do breit wigner fit to get fit parameter for counting
  TH1F *h_mMass_Phi2_Sig_total[pt_total_New_phi][Centrality_total][EtaGap_total];
  TH1F *h_mMass_Phi3_Sig_total[pt_total_New_phi][Centrality_total][EtaGap_total];

  TF1 *f_gauss_2[pt_total_New_phi][Centrality_total][EtaGap_total];
  TF1 *f_gauss_3[pt_total_New_phi][Centrality_total][EtaGap_total];

  Float_t gaussSig = 0.60;
  Float_t ParGaus_2[pt_total_New_phi][Centrality_total][EtaGap_total][2];
  Float_t ParGaus_3[pt_total_New_phi][Centrality_total][EtaGap_total][2];

  TF1 *f_bw_2[pt_total_New_phi][Centrality_total][EtaGap_total];
  TF1 *f_bw_3[pt_total_New_phi][Centrality_total][EtaGap_total];

  Float_t ParBW_2[pt_total_New_phi][Centrality_total][EtaGap_total][2];
  Float_t ParBW_3[pt_total_New_phi][Centrality_total][EtaGap_total][2];
  for(Int_t i = 0; i < pt_total_New_phi; i++) // pt bin
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	TString HistName;
	TString FuncName_g;
	TString FuncName_bw;

	HistName = Form("Sig_pt_%d_Centality_%d_EtaGap_%d_2nd",i,j,l);
	FuncName_g = Form("f_gauss_pt_%d_Centrality_%d_EtaGap_%d_2nd",i,j,l);
	FuncName_bw = Form("f_bw_pt_%d_Centrality_%d_EtaGap_%d_2nd",i,j,l);
	for(Int_t m = 0; m < Phi_Psi_total; m++)
	{
	  if(m == 0)
	  {
	    h_mMass_Phi2_Sig_total[i][j][l] = (TH1F*)h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Clone(HistName.Data());
	    h_mMass_Phi2_Sig_total[i][j][l]->SetTitle("");
	    f_gauss_2[i][j][l] = new TF1(FuncName_g.Data(),Gaussion,BW_Start,BW_Stop,3);
	    f_gauss_2[i][j][l]->SetParameter(0,ParFit2_total[i][j][l][0]);
//	    f_gauss_2[i][j][l]->SetParameter(1,ParFit2_total[i][j][l][1]);
	    f_gauss_2[i][j][l]->SetParameter(1,0.07);
	    f_gauss_2[i][j][l]->SetParameter(2,1000);
	    f_gauss_2[i][j][l]->SetRange(ParFit2_total[i][j][l][0]-gaussSig*ParFit2_total[i][j][l][1],ParFit2_total[i][j][l][0]+gaussSig*ParFit2_total[i][j][l][1]);

	    f_bw_2[i][j][l] = new TF1(FuncName_bw.Data(),BreitWigner,BW_Start,BW_Stop,3);
	    f_bw_2[i][j][l]->SetParameter(0,ParFit2_total[i][j][l][0]);
	    f_bw_2[i][j][l]->SetParameter(1,ParFit2_total[i][j][l][1]);
	    f_bw_2[i][j][l]->SetParameter(2,1000);
	    f_bw_2[i][j][l]->SetRange(BW_Start,BW_Stop);
	  }
	  else
	  {
	    h_mMass_Phi2_Sig_total[i][j][l]->Add(h_mMass_Phi2_SM_New_Bin[i][j][l][m],1.0);
	  }
	}
	h_mMass_Phi2_Sig_total[i][j][l]->Fit(f_gauss_2[i][j][l],"MQNR");
	ParGaus_2[i][j][l][0] = f_gauss_2[i][j][l]->GetParameter(0);
	ParGaus_2[i][j][l][1] = f_gauss_2[i][j][l]->GetParameter(1);

	h_mMass_Phi2_Sig_total[i][j][l]->Fit(f_bw_2[i][j][l],"MQNR");
	ParBW_2[i][j][l][0] = f_bw_2[i][j][l]->GetParameter(0);
	ParBW_2[i][j][l][1] = f_bw_2[i][j][l]->GetParameter(1);

	HistName = Form("Sig_pt_%d_Centality_%d_EtaGap_%d_3rd",i,j,l);
	FuncName_g = Form("f_gauss_pt_%d_Centrality_%d_EtaGap_%d_3rd",i,j,l);
	FuncName_bw = Form("f_bw_pt_%d_Centrality_%d_EtaGap_%d_3rd",i,j,l);
	for(Int_t m = 0; m < Phi_Psi_total; m++)
	{
	  if(m == 0)
	  {
	    h_mMass_Phi3_Sig_total[i][j][l] = (TH1F*)h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Clone(HistName.Data());
	    h_mMass_Phi3_Sig_total[i][j][l]->SetTitle("");
	    f_gauss_3[i][j][l] = new TF1(FuncName_g.Data(),Gaussion,BW_Start,BW_Stop,3);
	    f_gauss_3[i][j][l]->SetParameter(0,ParFit3_total[i][j][l][0]);
//	    f_gauss_3[i][j][l]->SetParameter(1,ParFit3_total[i][j][l][1]);
	    f_gauss_3[i][j][l]->SetParameter(1,0.07);
	    f_gauss_3[i][j][l]->SetParameter(2,1000);
	    f_gauss_3[i][j][l]->SetRange(ParFit3_total[i][j][l][0]-gaussSig*ParFit3_total[i][j][l][1],ParFit3_total[i][j][l][0]+gaussSig*ParFit3_total[i][j][l][1]);

	    f_bw_3[i][j][l] = new TF1(FuncName_bw.Data(),BreitWigner,BW_Start,BW_Stop,3);
	    f_bw_3[i][j][l]->SetParameter(0,ParFit3_total[i][j][l][0]);
	    f_bw_3[i][j][l]->SetParameter(1,ParFit3_total[i][j][l][1]);
	    f_bw_3[i][j][l]->SetParameter(2,1000);
	    f_bw_3[i][j][l]->SetRange(BW_Start,BW_Stop);
	  }
	  else
	  {
	    h_mMass_Phi3_Sig_total[i][j][l]->Add(h_mMass_Phi3_SM_New_Bin[i][j][l][m],1.0);
	  }
	}
	h_mMass_Phi3_Sig_total[i][j][l]->Fit(f_gauss_3[i][j][l],"MQNR");
	ParGaus_3[i][j][l][0] = f_gauss_3[i][j][l]->GetParameter(0);
	ParGaus_3[i][j][l][1] = f_gauss_3[i][j][l]->GetParameter(1);

	h_mMass_Phi3_Sig_total[i][j][l]->Fit(f_bw_3[i][j][l],"MQNR");
	ParBW_3[i][j][l][0] = f_bw_3[i][j][l]->GetParameter(0);
	ParBW_3[i][j][l][1] = f_bw_3[i][j][l]->GetParameter(1);
      }
    }
  }

  // Draw InvMass distribution for all pt bin with Gaussian fit
  TCanvas *c2_phi_Psi[pt_total_New_phi][Centrality_total][EtaGap_total];
  TCanvas *c3_phi_Psi[pt_total_New_phi][Centrality_total][EtaGap_total];
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	TString CanName;

	CanName = Form("c2_phi_Psi_pt_%d_Centrality_%d_EtaGap_%d",i,j,l);
	c2_phi_Psi[i][j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
	c2_phi_Psi[i][j][l]->Divide(3,3);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  c2_phi_Psi[i][j][l]->cd(m+1);
	  c2_phi_Psi[i][j][l]->cd(m+1)->SetLeftMargin(0.2);
	  c2_phi_Psi[i][j][l]->cd(m+1)->SetBottomMargin(0.2);
	  c2_phi_Psi[i][j][l]->cd(m+1)->SetTicks(1,1);
	  c2_phi_Psi[i][j][l]->cd(m+1)->SetGrid(0,0);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetStats(0);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetTitle("");
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetNdivisions(505,"X");
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetNdivisions(505,"Y");
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetXaxis()->SetLabelSize(0.05);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetYaxis()->SetLabelSize(0.05);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetYaxis()->SetTitle("Counts");
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetXaxis()->SetTitleSize(0.06);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetYaxis()->SetTitleSize(0.06);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetYaxis()->SetTitleOffset(1.7);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetXaxis()->CenterTitle();
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetYaxis()->CenterTitle();
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Draw("pE");
	  Float_t x1 = ParGaus_2[i][j][l][0] - nSigmaPhi*ParGaus_2[i][j][l][1];
	  Float_t x2 = ParGaus_2[i][j][l][0] + nSigmaPhi*ParGaus_2[i][j][l][1];
	  Float_t y = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinContent(h_mMass_Phi2_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_2[i][j][l][0]));
	  PlotLine(x1,x1,0,y,4,2,2);
	  PlotLine(x2,x2,0,y,4,2,2);
	  PlotLine(0.98,1.05,0,0,1,2,2);
	}
	c2_phi_Psi[i][j][l]->cd(8);
	c2_phi_Psi[i][j][l]->cd(8)->SetLeftMargin(0.2);
	c2_phi_Psi[i][j][l]->cd(8)->SetBottomMargin(0.2);
	c2_phi_Psi[i][j][l]->cd(8)->SetTicks(1,1);
	c2_phi_Psi[i][j][l]->cd(8)->SetGrid(0,0);
	h_mMass_Phi2_Sig_total[i][j][l]->SetStats(0);
	h_mMass_Phi2_Sig_total[i][j][l]->SetTitle("Merged Inv. Mass distribution");
	h_mMass_Phi2_Sig_total[i][j][l]->SetNdivisions(505,"X");
	h_mMass_Phi2_Sig_total[i][j][l]->SetNdivisions(505,"Y");
	h_mMass_Phi2_Sig_total[i][j][l]->GetXaxis()->SetLabelSize(0.05);
	h_mMass_Phi2_Sig_total[i][j][l]->GetYaxis()->SetLabelSize(0.05);
	h_mMass_Phi2_Sig_total[i][j][l]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
	h_mMass_Phi2_Sig_total[i][j][l]->GetYaxis()->SetTitle("Counts");
	h_mMass_Phi2_Sig_total[i][j][l]->GetXaxis()->SetTitleSize(0.06);
	h_mMass_Phi2_Sig_total[i][j][l]->GetYaxis()->SetTitleSize(0.06);
	h_mMass_Phi2_Sig_total[i][j][l]->GetYaxis()->SetTitleOffset(1.7);
	h_mMass_Phi2_Sig_total[i][j][l]->GetXaxis()->CenterTitle();
	h_mMass_Phi2_Sig_total[i][j][l]->GetYaxis()->CenterTitle();
	h_mMass_Phi2_Sig_total[i][j][l]->Draw("pE");
	f_gauss_2[i][j][l]->SetLineColor(2);
	f_gauss_2[i][j][l]->Draw("l same");
	PlotLine(0.98,1.05,0,0,1,2,2);

	CanName = Form("c3_phi_Psi_pt_%d_Centrality_%d_EtaGap_%d",i,j,l);
	c3_phi_Psi[i][j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
	c3_phi_Psi[i][j][l]->Divide(3,3);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  c3_phi_Psi[i][j][l]->cd(m+1);
	  c3_phi_Psi[i][j][l]->cd(m+1)->SetLeftMargin(0.2);
	  c3_phi_Psi[i][j][l]->cd(m+1)->SetBottomMargin(0.2);
	  c3_phi_Psi[i][j][l]->cd(m+1)->SetTicks(1,1);
	  c3_phi_Psi[i][j][l]->cd(m+1)->SetGrid(0,0);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetStats(0);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetTitle("");
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetNdivisions(505,"X");
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetNdivisions(505,"Y");
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetXaxis()->SetLabelSize(0.05);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetYaxis()->SetLabelSize(0.05);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetYaxis()->SetTitle("Counts");
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetXaxis()->SetTitleSize(0.06);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetYaxis()->SetTitleSize(0.06);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetYaxis()->SetTitleOffset(1.7);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetXaxis()->CenterTitle();
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetYaxis()->CenterTitle();
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Draw("pE");
	  Float_t x1 = ParGaus_3[i][j][l][0] - nSigmaPhi*ParGaus_3[i][j][l][1];
	  Float_t x2 = ParGaus_3[i][j][l][0] + nSigmaPhi*ParGaus_3[i][j][l][1];
	  Float_t y = h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetBinContent(h_mMass_Phi3_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_3[i][j][l][0]));
	  PlotLine(x1,x1,0,y,4,2,2);
	  PlotLine(x2,x2,0,y,4,2,2);
	  PlotLine(0.98,1.05,0,0,1,2,2);
	}
	c3_phi_Psi[i][j][l]->cd(8);
	c3_phi_Psi[i][j][l]->cd(8)->SetLeftMargin(0.2);
	c3_phi_Psi[i][j][l]->cd(8)->SetBottomMargin(0.2);
	c3_phi_Psi[i][j][l]->cd(8)->SetTicks(1,1);
	c3_phi_Psi[i][j][l]->cd(8)->SetGrid(0,0);
	h_mMass_Phi3_Sig_total[i][j][l]->SetStats(0);
	h_mMass_Phi3_Sig_total[i][j][l]->SetTitle("Merged Inv. Mass distribution");
	h_mMass_Phi3_Sig_total[i][j][l]->SetNdivisions(505,"X");
	h_mMass_Phi3_Sig_total[i][j][l]->SetNdivisions(505,"Y");
	h_mMass_Phi3_Sig_total[i][j][l]->GetXaxis()->SetLabelSize(0.05);
	h_mMass_Phi3_Sig_total[i][j][l]->GetYaxis()->SetLabelSize(0.05);
	h_mMass_Phi3_Sig_total[i][j][l]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
	h_mMass_Phi3_Sig_total[i][j][l]->GetYaxis()->SetTitle("Counts");
	h_mMass_Phi3_Sig_total[i][j][l]->GetXaxis()->SetTitleSize(0.06);
	h_mMass_Phi3_Sig_total[i][j][l]->GetYaxis()->SetTitleSize(0.06);
	h_mMass_Phi3_Sig_total[i][j][l]->GetYaxis()->SetTitleOffset(1.7);
	h_mMass_Phi3_Sig_total[i][j][l]->GetXaxis()->CenterTitle();
	h_mMass_Phi3_Sig_total[i][j][l]->GetYaxis()->CenterTitle();
	h_mMass_Phi3_Sig_total[i][j][l]->Draw("pE");
	f_gauss_3[i][j][l]->SetLineColor(2);
	f_gauss_3[i][j][l]->Draw("l same");
	PlotLine(0.98,1.05,0,0,1,2,2);
      }
    }
  }

  // calculate total counts and errors for each phi-Psi bin by bin counting
  TH1F *h_Counts2[pt_total_New_phi][Centrality_total][EtaGap_total];
  TH1F *h_Counts3[pt_total_New_phi][Centrality_total][EtaGap_total];
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	TString HistName;

	HistName = Form("Counts2_pt_%d_Centrality_%d_EtaGap_%d",i,j,l);
	h_Counts2[i][j][l] = new TH1F(HistName.Data(),HistName.Data(),7,0,PI_max[0]);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  Float_t counts = 0.0;
	  Float_t errors = 0.0;
	  Float_t bin_center = PI_max[0]/14.0+m*PI_max[0]/7.0;
	  Int_t bin_start = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_2[i][j][l][0]-nSigmaPhi*ParGaus_2[i][j][l][1]);
          Int_t bin_stop  = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_2[i][j][l][0]+nSigmaPhi*ParGaus_2[i][j][l][1]);
	  for(Int_t bin = bin_start; bin <= bin_stop; bin++)
	  {
	    counts += h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinContent(bin);
	    errors += h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinError(bin)*h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinError(bin);
	  }
	  h_Counts2[i][j][l]->SetBinContent(h_Counts2[i][j][l]->FindBin(bin_center),counts);
	  h_Counts2[i][j][l]->SetBinError(h_Counts2[i][j][l]->FindBin(bin_center),TMath::Sqrt(errors));
	}

	HistName = Form("Counts3_pt_%d_Centrality_%d_EtaGap_%d",i,j,l);
	h_Counts3[i][j][l] = new TH1F(HistName.Data(),HistName.Data(),7,0,PI_max[1]);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  Float_t counts = 0.0;
	  Float_t errors = 0.0;
	  Float_t bin_center = PI_max[1]/14.0+m*PI_max[1]/7.0;
	  Int_t bin_start = h_mMass_Phi3_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_3[i][j][l][0]-nSigmaPhi*ParGaus_3[i][j][l][1]);
	  Int_t bin_stop  = h_mMass_Phi3_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_3[i][j][l][0]+nSigmaPhi*ParGaus_3[i][j][l][1]);
	  for(Int_t bin = bin_start; bin <= bin_stop; bin++)
	  {
	    counts += h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetBinContent(bin);
	    errors += h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetBinError(bin)*h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetBinError(bin);
	  }
	  h_Counts3[i][j][l]->SetBinContent(h_Counts3[i][j][l]->FindBin(bin_center),counts);
	  h_Counts3[i][j][l]->SetBinError(h_Counts3[i][j][l]->FindBin(bin_center),TMath::Sqrt(errors));
	}
      }
    }
  }

  // Draw phi-Psi distribution with counting
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	c2_phi_Psi[i][j][l]->cd(9);
	c2_phi_Psi[i][j][l]->cd(9)->SetLeftMargin(0.2);
	c2_phi_Psi[i][j][l]->cd(9)->SetBottomMargin(0.2);
	c2_phi_Psi[i][j][l]->cd(9)->SetTicks(1,1);
	c2_phi_Psi[i][j][l]->cd(9)->SetGrid(0,0);
	h_Counts2[i][j][l]->SetStats(0);
	TString Title = Form("%1.1f < p_{T} < %1.1f GeV/c",pt_low_phi[i],pt_up_phi[i]);
	h_Counts2[i][j][l]->SetTitle(Title.Data());
	h_Counts2[i][j][l]->SetNdivisions(505,"X");
	h_Counts2[i][j][l]->SetNdivisions(505,"Y");
	h_Counts2[i][j][l]->GetXaxis()->SetLabelSize(0.05);
	h_Counts2[i][j][l]->GetYaxis()->SetLabelSize(0.05);
	h_Counts2[i][j][l]->GetXaxis()->SetTitle("#phi-#Psi_{2}");
	h_Counts2[i][j][l]->GetYaxis()->SetTitle("Counts");
	h_Counts2[i][j][l]->GetXaxis()->SetTitleSize(0.06);
	h_Counts2[i][j][l]->GetYaxis()->SetTitleSize(0.06);
	h_Counts2[i][j][l]->GetYaxis()->SetTitleOffset(1.7);
	h_Counts2[i][j][l]->GetXaxis()->CenterTitle();
	h_Counts2[i][j][l]->GetYaxis()->CenterTitle();
	Int_t bin_low = h_Counts2[i][j][l]->FindBin(PI_max[0]/14.0+6*PI_max[0]/7.0);
	Int_t bin_up  = h_Counts2[i][j][l]->FindBin(PI_max[0]/14.0+0*PI_max[0]/7.0);
	h_Counts2[i][j][l]->GetYaxis()->SetRangeUser(h_Counts2[i][j][l]->GetBinContent(bin_low)*0.8,h_Counts2[i][j][l]->GetBinContent(bin_up)*1.2);
	h_Counts2[i][j][l]->Draw("pE");
	TLegend *leg2 = new TLegend(0.5,0.7,0.88,0.8);
	leg2->SetFillColor(0);
	leg2->SetBorderSize(0);
	leg2->AddEntry(h_Counts2[i][j][l],"Counting","p");
	leg2->Draw("same");

	c3_phi_Psi[i][j][l]->cd(9);
	c3_phi_Psi[i][j][l]->cd(9)->SetLeftMargin(0.2);
	c3_phi_Psi[i][j][l]->cd(9)->SetBottomMargin(0.2);
	c3_phi_Psi[i][j][l]->cd(9)->SetTicks(1,1);
	c3_phi_Psi[i][j][l]->cd(9)->SetGrid(0,0);
	h_Counts3[i][j][l]->SetStats(0);
	h_Counts3[i][j][l]->SetTitle(Title.Data());
	h_Counts3[i][j][l]->SetNdivisions(505,"X");
	h_Counts3[i][j][l]->SetNdivisions(505,"Y");
	h_Counts3[i][j][l]->GetXaxis()->SetLabelSize(0.05);
	h_Counts3[i][j][l]->GetYaxis()->SetLabelSize(0.05);
	h_Counts3[i][j][l]->GetXaxis()->SetTitle("#phi-#Psi_{3}");
	h_Counts3[i][j][l]->GetYaxis()->SetTitle("Counts");
	h_Counts3[i][j][l]->GetXaxis()->SetTitleSize(0.06);
	h_Counts3[i][j][l]->GetYaxis()->SetTitleSize(0.06);
	h_Counts3[i][j][l]->GetYaxis()->SetTitleOffset(1.7);
	h_Counts3[i][j][l]->GetXaxis()->CenterTitle();
	h_Counts3[i][j][l]->GetYaxis()->CenterTitle();
	bin_low = h_Counts3[i][j][l]->FindBin(PI_max[1]/14.0+6*PI_max[1]/7.0);
	bin_up  = h_Counts3[i][j][l]->FindBin(PI_max[1]/14.0+0*PI_max[1]/7.0);
	h_Counts3[i][j][l]->GetYaxis()->SetRangeUser(h_Counts3[i][j][l]->GetBinContent(bin_low)*0.8,h_Counts3[i][j][l]->GetBinContent(bin_up)*1.2);
	h_Counts3[i][j][l]->Draw("pE");
	TLegend *leg3 = new TLegend(0.5,0.7,0.88,0.8);
	leg3->SetFillColor(0);
	leg3->SetBorderSize(0);
	leg3->AddEntry(h_Counts3[i][j][l],"Counting","p");
	leg3->Draw("same");
      }
    }
  }

  // Breit Wigner fit for phi-Psi bin 
  TF1 *f_bw_phi_2[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TF1 *f_bw_phi_3[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];

  // calculate total counts and errors for each phi-Psi bin by BW Integration
  TH1F *h_Counts2_bw[pt_total_New_phi][Centrality_total][EtaGap_total];
  TH1F *h_Counts3_bw[pt_total_New_phi][Centrality_total][EtaGap_total];
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	TString HistName;

	HistName = Form("Counts2_pt_%d_Centrality_%d_EtaGap_%d_bw",i,j,l);
	h_Counts2_bw[i][j][l] = new TH1F(HistName.Data(),HistName.Data(),7,0,PI_max[0]);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  TString FuncName_bw = Form("f_bw_pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd",i,j,l,m);
	  f_bw_phi_2[i][j][l][m] = new TF1(FuncName_bw.Data(),BreitWigner,BW_Start,BW_Stop,3);;
	  f_bw_phi_2[i][j][l][m]->FixParameter(0,ParBW_2[i][j][l][0]);
	  f_bw_phi_2[i][j][l][m]->FixParameter(1,ParBW_2[i][j][l][1]);
	  f_bw_phi_2[i][j][l][m]->SetParameter(2,1000);
	  f_bw_phi_2[i][j][l][m]->SetRange(BW_Start,BW_Stop);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Fit(f_bw_phi_2[i][j][l][m],"NMQR");
	  Float_t bin_width = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinWidth(1);
	  Float_t Inte_start = ParGaus_2[i][j][l][0]-nSigmaPhi*ParGaus_2[i][j][l][1]-0.5*bin_width;
          Float_t Inte_stop  = ParGaus_2[i][j][l][0]+nSigmaPhi*ParGaus_2[i][j][l][1]+0.5*bin_width;
	  Float_t counts = f_bw_phi_2[i][j][l][m]->Integral(Inte_start,Inte_stop)/bin_width;
	  Float_t errors = f_bw_phi_2[i][j][l][m]->IntegralError(Inte_start,Inte_stop)/bin_width;
	  Float_t bin_center = PI_max[0]/14.0+m*PI_max[0]/7.0;
	  h_Counts2_bw[i][j][l]->SetBinContent(h_Counts2_bw[i][j][l]->FindBin(bin_center),counts);
	  h_Counts2_bw[i][j][l]->SetBinError(h_Counts2_bw[i][j][l]->FindBin(bin_center),errors);
	}

	HistName = Form("Counts3_pt_%d_Centrality_%d_EtaGap_%d_bw",i,j,l);
	h_Counts3_bw[i][j][l] = new TH1F(HistName.Data(),HistName.Data(),7,0,PI_max[1]);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  TString FuncName_bw = Form("f_bw_pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd",i,j,l,m);
	  f_bw_phi_3[i][j][l][m] = new TF1(FuncName_bw.Data(),BreitWigner,BW_Start,BW_Stop,3);;
	  f_bw_phi_3[i][j][l][m]->FixParameter(0,ParBW_3[i][j][l][0]);
	  f_bw_phi_3[i][j][l][m]->FixParameter(1,ParBW_3[i][j][l][1]);
	  f_bw_phi_3[i][j][l][m]->SetParameter(2,1000);
	  f_bw_phi_3[i][j][l][m]->SetRange(BW_Start,BW_Stop);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Fit(f_bw_phi_3[i][j][l][m],"NMQR");
	  Float_t bin_width = h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetBinWidth(1);
	  Float_t Inte_start = ParGaus_3[i][j][l][0]-nSigmaPhi*ParGaus_3[i][j][l][1]-0.5*bin_width;
          Float_t Inte_stop  = ParGaus_3[i][j][l][0]+nSigmaPhi*ParGaus_3[i][j][l][1]+0.5*bin_width;
	  Float_t counts = f_bw_phi_3[i][j][l][m]->Integral(Inte_start,Inte_stop)/bin_width;
	  Float_t errors = f_bw_phi_3[i][j][l][m]->IntegralError(Inte_start,Inte_stop)/bin_width;
	  Float_t bin_center = PI_max[1]/14.0+m*PI_max[1]/7.0;
	  h_Counts3_bw[i][j][l]->SetBinContent(h_Counts3_bw[i][j][l]->FindBin(bin_center),counts);
	  h_Counts3_bw[i][j][l]->SetBinError(h_Counts3_bw[i][j][l]->FindBin(bin_center),errors);
	}
      }
    }
  }

  // Draw Breit Wigner fit for InvMass distribution of all pt bin
  TCanvas *c2_phi_Psi_bw[pt_total_New_phi][Centrality_total][EtaGap_total];
  TCanvas *c3_phi_Psi_bw[pt_total_New_phi][Centrality_total][EtaGap_total];
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	TString CanName;

	CanName = Form("c2_phi_Psi_pt_%d_Centrality_%d_EtaGap_%d_bw",i,j,l);
	c2_phi_Psi_bw[i][j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
	c2_phi_Psi_bw[i][j][l]->Divide(3,3);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  c2_phi_Psi_bw[i][j][l]->cd(m+1);
	  c2_phi_Psi_bw[i][j][l]->cd(m+1)->SetLeftMargin(0.2);
	  c2_phi_Psi_bw[i][j][l]->cd(m+1)->SetBottomMargin(0.2);
	  c2_phi_Psi_bw[i][j][l]->cd(m+1)->SetTicks(1,1);
	  c2_phi_Psi_bw[i][j][l]->cd(m+1)->SetGrid(0,0);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Draw("pE");
	  f_bw_phi_2[i][j][l][m]->SetLineColor(2);
	  f_bw_phi_2[i][j][l][m]->Draw("l Same");
	  Float_t x1 = ParGaus_2[i][j][l][0] - nSigmaPhi*ParGaus_2[i][j][l][1];
	  Float_t x2 = ParGaus_2[i][j][l][0] + nSigmaPhi*ParGaus_2[i][j][l][1];
	  Float_t y = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinContent(h_mMass_Phi2_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_2[i][j][l][0]));
	  PlotLine(x1,x1,0,y,4,2,2);
	  PlotLine(x2,x2,0,y,4,2,2);
	  PlotLine(0.98,1.05,0,0,1,2,2);
	}
	c2_phi_Psi_bw[i][j][l]->cd(8);
	c2_phi_Psi_bw[i][j][l]->cd(8)->SetLeftMargin(0.2);
	c2_phi_Psi_bw[i][j][l]->cd(8)->SetBottomMargin(0.2);
	c2_phi_Psi_bw[i][j][l]->cd(8)->SetTicks(1,1);
	c2_phi_Psi_bw[i][j][l]->cd(8)->SetGrid(0,0);
	h_mMass_Phi2_Sig_total[i][j][l]->Draw("pE");
	f_bw_2[i][j][l]->SetLineColor(2);
	f_bw_2[i][j][l]->Draw("l same");
	PlotLine(0.98,1.05,0,0,1,2,2);

	CanName = Form("c3_phi_Psi_pt_%d_Centrality_%d_EtaGap_%d_bw",i,j,l);
	c3_phi_Psi_bw[i][j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
	c3_phi_Psi_bw[i][j][l]->Divide(3,3);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  c3_phi_Psi_bw[i][j][l]->cd(m+1);
	  c3_phi_Psi_bw[i][j][l]->cd(m+1)->SetLeftMargin(0.2);
	  c3_phi_Psi_bw[i][j][l]->cd(m+1)->SetBottomMargin(0.2);
	  c3_phi_Psi_bw[i][j][l]->cd(m+1)->SetTicks(1,1);
	  c3_phi_Psi_bw[i][j][l]->cd(m+1)->SetGrid(0,0);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Draw("pE");
	  f_bw_phi_3[i][j][l][m]->SetLineColor(2);
	  f_bw_phi_3[i][j][l][m]->Draw("l Same");
	  Float_t x1 = ParGaus_3[i][j][l][0] - nSigmaPhi*ParGaus_3[i][j][l][1];
	  Float_t x2 = ParGaus_3[i][j][l][0] + nSigmaPhi*ParGaus_3[i][j][l][1];
	  Float_t y = h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetBinContent(h_mMass_Phi3_SM_New_Bin[i][j][l][m]->FindBin(ParGaus_3[i][j][l][0]));
	  PlotLine(x1,x1,0,y,4,2,2);
	  PlotLine(x2,x2,0,y,4,2,2);
	  PlotLine(0.98,1.05,0,0,1,2,2);
	}
	c3_phi_Psi_bw[i][j][l]->cd(8);
	c3_phi_Psi_bw[i][j][l]->cd(8)->SetLeftMargin(0.2);
	c3_phi_Psi_bw[i][j][l]->cd(8)->SetBottomMargin(0.2);
	c3_phi_Psi_bw[i][j][l]->cd(8)->SetTicks(1,1);
	c3_phi_Psi_bw[i][j][l]->cd(8)->SetGrid(0,0);
	h_mMass_Phi3_Sig_total[i][j][l]->Draw("pE");
	f_bw_3[i][j][l]->SetLineColor(2);
	f_bw_3[i][j][l]->Draw("l same");
	PlotLine(0.98,1.05,0,0,1,2,2);
      }
    }
  }

  // Draw phi-Psi distribution with Breit Wignar Fit 
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	c2_phi_Psi_bw[i][j][l]->cd(9);
	c2_phi_Psi_bw[i][j][l]->cd(9)->SetLeftMargin(0.2);
	c2_phi_Psi_bw[i][j][l]->cd(9)->SetBottomMargin(0.2);
	c2_phi_Psi_bw[i][j][l]->cd(9)->SetTicks(1,1);
	c2_phi_Psi_bw[i][j][l]->cd(9)->SetGrid(0,0);
	h_Counts2_bw[i][j][l]->SetStats(0);
	TString Title = Form("%1.1f < p_{T} < %1.1f GeV/c",pt_low_phi[i],pt_up_phi[i]);
	h_Counts2_bw[i][j][l]->SetTitle(Title.Data());
	h_Counts2_bw[i][j][l]->SetNdivisions(505,"X");
	h_Counts2_bw[i][j][l]->SetNdivisions(505,"Y");
	h_Counts2_bw[i][j][l]->GetXaxis()->SetLabelSize(0.05);
	h_Counts2_bw[i][j][l]->GetYaxis()->SetLabelSize(0.05);
	h_Counts2_bw[i][j][l]->GetXaxis()->SetTitle("#phi-#Psi_{2}");
	h_Counts2_bw[i][j][l]->GetYaxis()->SetTitle("Counts");
	h_Counts2_bw[i][j][l]->GetXaxis()->SetTitleSize(0.06);
	h_Counts2_bw[i][j][l]->GetYaxis()->SetTitleSize(0.06);
	h_Counts2_bw[i][j][l]->GetYaxis()->SetTitleOffset(1.7);
	h_Counts2_bw[i][j][l]->GetXaxis()->CenterTitle();
	h_Counts2_bw[i][j][l]->GetYaxis()->CenterTitle();
	Int_t bin_low = h_Counts2_bw[i][j][l]->FindBin(PI_max[0]/14.0+6*PI_max[0]/7.0);
	Int_t bin_up  = h_Counts2_bw[i][j][l]->FindBin(PI_max[0]/14.0+0*PI_max[0]/7.0);
	h_Counts2_bw[i][j][l]->GetYaxis()->SetRangeUser(h_Counts2_bw[i][j][l]->GetBinContent(bin_low)*0.8,h_Counts2_bw[i][j][l]->GetBinContent(bin_up)*1.2);
	h_Counts2_bw[i][j][l]->Draw("pE");
	TLegend *leg2 = new TLegend(0.5,0.7,0.88,0.8);
	leg2->SetFillColor(0);
	leg2->SetBorderSize(0);
	leg2->AddEntry(h_Counts2_bw[i][j][l],"Breit Wigner","p");
	leg2->Draw("same");

	c3_phi_Psi_bw[i][j][l]->cd(9);
	c3_phi_Psi_bw[i][j][l]->cd(9)->SetLeftMargin(0.2);
	c3_phi_Psi_bw[i][j][l]->cd(9)->SetBottomMargin(0.2);
	c3_phi_Psi_bw[i][j][l]->cd(9)->SetTicks(1,1);
	c3_phi_Psi_bw[i][j][l]->cd(9)->SetGrid(0,0);
	h_Counts3_bw[i][j][l]->SetStats(0);
	h_Counts3_bw[i][j][l]->SetTitle(Title.Data());
	h_Counts3_bw[i][j][l]->SetNdivisions(505,"X");
	h_Counts3_bw[i][j][l]->SetNdivisions(505,"Y");
	h_Counts3_bw[i][j][l]->GetXaxis()->SetLabelSize(0.05);
	h_Counts3_bw[i][j][l]->GetYaxis()->SetLabelSize(0.05);
	h_Counts3_bw[i][j][l]->GetXaxis()->SetTitle("#phi-#Psi_{3}");
	h_Counts3_bw[i][j][l]->GetYaxis()->SetTitle("Counts");
	h_Counts3_bw[i][j][l]->GetXaxis()->SetTitleSize(0.06);
	h_Counts3_bw[i][j][l]->GetYaxis()->SetTitleSize(0.06);
	h_Counts3_bw[i][j][l]->GetYaxis()->SetTitleOffset(1.7);
	h_Counts3_bw[i][j][l]->GetXaxis()->CenterTitle();
	h_Counts3_bw[i][j][l]->GetYaxis()->CenterTitle();
	bin_low = h_Counts3_bw[i][j][l]->FindBin(PI_max[1]/14.0+6*PI_max[1]/7.0);
	bin_up  = h_Counts3_bw[i][j][l]->FindBin(PI_max[1]/14.0+0*PI_max[1]/7.0);
	h_Counts3_bw[i][j][l]->GetYaxis()->SetRangeUser(h_Counts3_bw[i][j][l]->GetBinContent(bin_low)*0.8,h_Counts3_bw[i][j][l]->GetBinContent(bin_up)*1.2);
	h_Counts3_bw[i][j][l]->Draw("pE");
	TLegend *leg3 = new TLegend(0.5,0.7,0.88,0.8);
	leg3->SetFillColor(0);
	leg3->SetBorderSize(0);
	leg3->AddEntry(h_Counts3_bw[i][j][l],"Breit Wigner","p");
	leg3->Draw("same");
      }
    }
  }

  // do cos fit to extract raw v2 and v3
  TF1 *f_phi2[pt_total_New_phi][Centrality_total][EtaGap_total];
  TF1 *f_phi3[pt_total_New_phi][Centrality_total][EtaGap_total];

  TF1 *f_phi2_bw[pt_total_New_phi][Centrality_total][EtaGap_total];
  TF1 *f_phi3_bw[pt_total_New_phi][Centrality_total][EtaGap_total];
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	TString Flow_phi;

	Flow_phi = Form("flow_pt_%d_Centrality_%d_etagap_%d_2nd_phi",i,j,l);
	f_phi2[i][j][l] = new TF1(Flow_phi.Data(),flow_2,0.0,PI_max[0],2);
	f_phi2[i][j][l]->SetParameter(0,2.0);
	f_phi2[i][j][l]->SetParameter(1,1.0);
	h_Counts2[i][j][l]->Fit(f_phi2[i][j][l],"NQMI");

	Flow_phi = Form("flow_pt_%d_Centrality_%d_etagap_%d_3rd_phi",i,j,l);
	f_phi3[i][j][l] = new TF1(Flow_phi.Data(),flow_3,0.0,PI_max[1],2);
	f_phi3[i][j][l]->SetParameter(0,2.0);
	f_phi3[i][j][l]->SetParameter(1,1.0);
	h_Counts3[i][j][l]->Fit(f_phi3[i][j][l],"NQMI");

	Flow_phi = Form("flow_pt_%d_Centrality_%d_etagap_%d_2nd_phi_bw",i,j,l);
	f_phi2_bw[i][j][l] = new TF1(Flow_phi.Data(),flow_2,0.0,PI_max[0],2);
	f_phi2_bw[i][j][l]->SetParameter(0,2.0);
	f_phi2_bw[i][j][l]->SetParameter(1,1.0);
	h_Counts2_bw[i][j][l]->Fit(f_phi2_bw[i][j][l],"NQMI");

	Flow_phi = Form("flow_pt_%d_Centrality_%d_etagap_%d_3rd_phi_bw",i,j,l);
	f_phi3_bw[i][j][l] = new TF1(Flow_phi.Data(),flow_3,0.0,PI_max[1],2);
	f_phi3_bw[i][j][l]->SetParameter(0,2.0);
	f_phi3_bw[i][j][l]->SetParameter(1,1.0);
	h_Counts3_bw[i][j][l]->Fit(f_phi3_bw[i][j][l],"NQMI");
      }
    }
  }

  // Draw the fit line
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	c2_phi_Psi[i][j][l]->cd(9);
	f_phi2[i][j][l]->SetLineColor(2);
	f_phi2[i][j][l]->Draw("l Same");

	c3_phi_Psi[i][j][l]->cd(9);
	f_phi3[i][j][l]->SetLineColor(2);
	f_phi3[i][j][l]->Draw("l Same");

	c2_phi_Psi_bw[i][j][l]->cd(9);
	f_phi2_bw[i][j][l]->SetLineColor(2);
	f_phi2_bw[i][j][l]->Draw("l Same");

	c3_phi_Psi_bw[i][j][l]->cd(9);
	f_phi3_bw[i][j][l]->SetLineColor(2);
	f_phi3_bw[i][j][l]->Draw("l Same");
      }
    }
  }

  // Draw phi-Psi distribution from Counting and Breit Wignar Fit in the same plot with fit line
  TCanvas *c_phi_Psi2[Centrality_total][EtaGap_total];
  TCanvas *c_phi_Psi3[Centrality_total][EtaGap_total];
  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString CanName;

      CanName = Form("c_phi_Psi2_Centrality_%d_EtaGap_%d",j,l);
      c_phi_Psi2[j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1600,1600);
      c_phi_Psi2[j][l]->Divide(4,4);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	c_phi_Psi2[j][l]->cd(i+1);
	c_phi_Psi2[j][l]->cd(i+1)->SetLeftMargin(0.2);
	c_phi_Psi2[j][l]->cd(i+1)->SetBottomMargin(0.2);
	c_phi_Psi2[j][l]->cd(i+1)->SetTicks(1,1);
	c_phi_Psi2[j][l]->cd(i+1)->SetGrid(0,0);

	h_Counts2[i][j][l]->SetMarkerStyle(20);
	h_Counts2[i][j][l]->SetMarkerSize(1.0);
	h_Counts2[i][j][l]->SetMarkerColor(1);
	h_Counts2[i][j][l]->Draw("pE");
	TLegend *leg2 = new TLegend(0.5,0.6,0.88,0.8);
	leg2->SetFillColor(0);
	leg2->SetBorderSize(0);
	leg2->AddEntry(h_Counts2[i][j][l],"Counting","p");
	f_phi2[i][j][l]->SetLineColor(1);
	f_phi2[i][j][l]->SetLineStyle(1);
	f_phi2[i][j][l]->Draw("l Same");

	h_Counts2_bw[i][j][l]->SetMarkerStyle(24);
	h_Counts2_bw[i][j][l]->SetMarkerSize(1.0);
	h_Counts2_bw[i][j][l]->SetMarkerColor(2);
	h_Counts2_bw[i][j][l]->Draw("pE same");
	leg2->AddEntry(h_Counts2_bw[i][j][l],"Breit Wigner","p");
	f_phi2_bw[i][j][l]->SetLineColor(2);
	f_phi2_bw[i][j][l]->SetLineStyle(2);
	f_phi2_bw[i][j][l]->Draw("l Same");
	leg2->Draw("same");
      }

      CanName = Form("c_phi_Psi3_Centrality_%d_EtaGap_%d",j,l);
      c_phi_Psi3[j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1600,1600);
      c_phi_Psi3[j][l]->Divide(4,4);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	c_phi_Psi3[j][l]->cd(i+1);
	c_phi_Psi3[j][l]->cd(i+1)->SetLeftMargin(0.2);
	c_phi_Psi3[j][l]->cd(i+1)->SetBottomMargin(0.2);
	c_phi_Psi3[j][l]->cd(i+1)->SetTicks(1,1);
	c_phi_Psi3[j][l]->cd(i+1)->SetGrid(0,0);

	h_Counts3[i][j][l]->SetMarkerStyle(20);
	h_Counts3[i][j][l]->SetMarkerSize(1.0);
	h_Counts3[i][j][l]->SetMarkerColor(1);
	h_Counts3[i][j][l]->Draw("pE");
	TLegend *leg3 = new TLegend(0.5,0.6,0.88,0.8);
	leg3->SetFillColor(0);
	leg3->SetBorderSize(0);
	leg3->AddEntry(h_Counts3[i][j][l],"Counting","p");
	f_phi3[i][j][l]->SetLineColor(1);
	f_phi3[i][j][l]->SetLineStyle(1);
	f_phi3[i][j][l]->Draw("l Same");

	h_Counts3_bw[i][j][l]->SetMarkerStyle(24);
	h_Counts3_bw[i][j][l]->SetMarkerSize(1.0);
	h_Counts3_bw[i][j][l]->SetMarkerColor(2);
	h_Counts3_bw[i][j][l]->Draw("pE same");
	leg3->AddEntry(h_Counts3_bw[i][j][l],"Breit Wigner","p");
	f_phi3_bw[i][j][l]->SetLineColor(2);
	f_phi3_bw[i][j][l]->SetLineStyle(2);
	f_phi3_bw[i][j][l]->Draw("l Same");
	leg3->Draw("same");
      }
    }
  }

  // fill raw v2 and v3 into histogram for counting and Breit Wignar
  TH1F *h_flow_2[Centrality_total][EtaGap_total];
  TH1F *h_flow_3[Centrality_total][EtaGap_total];

  TH1F *h_flow_2_bw[Centrality_total][EtaGap_total];
  TH1F *h_flow_3_bw[Centrality_total][EtaGap_total];
  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString HistName;

      HistName = Form("flow_2_Centrality_%d_EtaGap_%d_phi",j,l);
      h_flow_2[j][l] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,3.6);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	Int_t bin_center = h_flow_2[j][l]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
	Float_t bin_content = f_phi2[i][j][l]->GetParameter(1);
	Float_t bin_error = f_phi2[i][j][l]->GetParError(1);
	h_flow_2[j][l]->SetBinContent(bin_center,bin_content);
	h_flow_2[j][l]->SetBinError(bin_center,bin_error);
      }

      HistName = Form("flow_3_Centrality_%d_EtaGap_%d_phi",j,l);
      h_flow_3[j][l] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,3.6);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	Int_t bin_center = h_flow_3[j][l]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
	Float_t bin_content = f_phi3[i][j][l]->GetParameter(1);
	Float_t bin_error = f_phi3[i][j][l]->GetParError(1);
	h_flow_3[j][l]->SetBinContent(bin_center,bin_content);
	h_flow_3[j][l]->SetBinError(bin_center,bin_error);
      }

      HistName = Form("flow_2_Centrality_%d_EtaGap_%d_phi_bw",j,l);
      h_flow_2_bw[j][l] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,3.6);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	Int_t bin_center = h_flow_2_bw[j][l]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
	Float_t bin_content = f_phi2_bw[i][j][l]->GetParameter(1);
	Float_t bin_error = f_phi2_bw[i][j][l]->GetParError(1);
	h_flow_2_bw[j][l]->SetBinContent(bin_center,bin_content);
	h_flow_2_bw[j][l]->SetBinError(bin_center,bin_error);
      }

      HistName = Form("flow_3_Centrality_%d_EtaGap_%d_phi_bw",j,l);
      h_flow_3_bw[j][l] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,3.6);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	Int_t bin_center = h_flow_3_bw[j][l]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
	Float_t bin_content = f_phi3_bw[i][j][l]->GetParameter(1);
	Float_t bin_error = f_phi3_bw[i][j][l]->GetParError(1);
	h_flow_3_bw[j][l]->SetBinContent(bin_center,bin_content);
	h_flow_3_bw[j][l]->SetBinError(bin_center,bin_error);
      }
    }
  }

  // Draw raw v2 and v3 into phi-Psi distribution
  TCanvas *c_raw_v2[Centrality_total][EtaGap_total];
  TCanvas *c_raw_v3[Centrality_total][EtaGap_total];
  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString CanName;
      CanName = Form("c_raw_v2_Centrality_%d_EtaGap_%d",j,l);
      c_raw_v2[j][l] = new TCanvas(CanName.Data(),CanName.Data(),10,10,800,800);
      c_raw_v2[j][l]->cd();
      c_raw_v2[j][l]->cd()->SetLeftMargin(0.2);
      c_raw_v2[j][l]->cd()->SetBottomMargin(0.2);
      c_raw_v2[j][l]->cd()->SetTicks(1,1);
      c_raw_v2[j][l]->cd()->SetGrid(0,0);

      h_flow_2[j][l]->SetStats(0);
      h_flow_2[j][l]->SetMarkerStyle(20);
      h_flow_2[j][l]->SetMarkerSize(1.0);
      h_flow_2[j][l]->SetMarkerColor(1);
      h_flow_2[j][l]->DrawCopy("pE");

      h_flow_2_bw[j][l]->SetStats(0);
      h_flow_2_bw[j][l]->SetMarkerStyle(24);
      h_flow_2_bw[j][l]->SetMarkerSize(1.0);
      h_flow_2_bw[j][l]->SetMarkerColor(2);
      h_flow_2_bw[j][l]->DrawCopy("pE same");

      CanName = Form("c_raw_v3_Centrality_%d_EtaGap_%d",j,l);
      c_raw_v3[j][l] = new TCanvas(CanName.Data(),CanName.Data(),10,10,800,800);
      c_raw_v3[j][l]->cd();
      c_raw_v3[j][l]->cd()->SetLeftMargin(0.2);
      c_raw_v3[j][l]->cd()->SetBottomMargin(0.2);
      c_raw_v3[j][l]->cd()->SetTicks(1,1);
      c_raw_v3[j][l]->cd()->SetGrid(0,0);

      h_flow_3[j][l]->SetStats(0);
      h_flow_3[j][l]->SetMarkerStyle(20);
      h_flow_3[j][l]->SetMarkerSize(1.0);
      h_flow_3[j][l]->SetMarkerColor(1);
      h_flow_3[j][l]->DrawCopy("pE");

      h_flow_3_bw[j][l]->SetStats(0);
      h_flow_3_bw[j][l]->SetMarkerStyle(24);
      h_flow_3_bw[j][l]->SetMarkerSize(1.0);
      h_flow_3_bw[j][l]->SetMarkerColor(2);
      h_flow_3_bw[j][l]->DrawCopy("pE same");
    }
  }

  //--------------------------------------------------------------------------------------------------
  // calculate the final resolution correction
  // read the file from Same Event and Mixed Event
  TH1F *h_Yields_Phi_SE[9][EtaGap_total];
  TH1F *h_Yields_Phi_ME[9][EtaGap_total];
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString HistName;

      HistName = Form("Yields_Centrality_%d_EtaGap_%d_Phi_SE",j,l);
      h_Yields_Phi_SE[j][l] = (TH1F*)File_SE->Get(HistName.Data())->Clone();

      HistName = Form("Yields_Centrality_%d_EtaGap_%d_Phi_ME",j,l);
      h_Yields_Phi_ME[j][l] = (TH1F*)File_ME->Get(HistName.Data())->Clone();
    }
  }

  TH1F *h_Yields_Phi_SM[9][EtaGap_total];
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      Int_t bin_SE_start = h_Yields_Phi_SE[j][l]->FindBin(1.04);
      Int_t bin_SE_stop  = h_Yields_Phi_SE[j][l]->FindBin(1.05);
      Float_t Inte_SE = h_Yields_Phi_SE[j][l]->Integral(bin_SE_start,bin_SE_stop);
      h_Yields_Phi_SM[j][l] = (TH1F*)h_Yields_Phi_SE[j][l]->Clone();

      Int_t bin_ME_start = h_Yields_Phi_ME[j][l]->FindBin(1.04);
      Int_t bin_ME_stop  = h_Yields_Phi_ME[j][l]->FindBin(1.05);
      Float_t Inte_ME = h_Yields_Phi_ME[j][l]->Integral(bin_ME_start,bin_ME_stop);
      h_Yields_Phi_ME[j][l]->Scale(Inte_SE/Inte_ME);
      h_Yields_Phi_SM[j][l]->Add(h_Yields_Phi_ME[j][l],-1.0);
    }
  }

  TF1  *f_PolyBW_Phi_Yields[9][EtaGap_total];

  Float_t ParFit_Yields[9][EtaGap_total][5];

  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString HistName;

      HistName = Form("f_Centrality_%d_EtaGap_%d_Yields_Phi_SE",j,l);
      f_PolyBW_Phi_Yields[j][l] = new TF1(HistName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
      for(Int_t n_par = 0; n_par < 5; n_par++)
      {
	f_PolyBW_Phi_Yields[j][l]->ReleaseParameter(n_par);
      }
      f_PolyBW_Phi_Yields[j][l]->SetParameter(0,1.019);
      f_PolyBW_Phi_Yields[j][l]->SetParLimits(0,1.014,1.024);
      f_PolyBW_Phi_Yields[j][l]->SetParameter(1,0.0055);
      f_PolyBW_Phi_Yields[j][l]->SetParameter(2,10000);
      f_PolyBW_Phi_Yields[j][l]->SetParameter(3,-6000);
      f_PolyBW_Phi_Yields[j][l]->SetParameter(4,0.5);
      f_PolyBW_Phi_Yields[j][l]->SetRange(BW_Start,BW_Stop);

      h_Yields_Phi_SM[j][l]->Fit(f_PolyBW_Phi_Yields[j][l],"NQR");
      for(Int_t n_par = 0; n_par < 5; n_par++)
      {
	ParFit_Yields[j][l][n_par] = f_PolyBW_Phi_Yields[j][l]->GetParameter(n_par);
      }
    }
  }

  TF1 *f_Poly_Phi_Yields[9][EtaGap_total];

  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString FuncName;

      FuncName = Form("f_poly_Centrality_%d_EtaGap_%d_Yields",j,l);
      f_Poly_Phi_Yields[j][l] = new TF1(FuncName.Data(),Poly,BW_Start,BW_Stop,2);
      f_Poly_Phi_Yields[j][l]->FixParameter(0,ParFit_Yields[j][l][3]);
      f_Poly_Phi_Yields[j][l]->FixParameter(1,ParFit_Yields[j][l][4]);
      h_Yields_Phi_SM[j][l]->Add(f_Poly_Phi_Yields[j][l],-1.0);
    }
  }

  // calculate total counts and for each centrality bin and eta bin
  Float_t Counts_Yields[9][EtaGap_total];
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      Counts_Yields[j][l] = 0.0;
      Int_t bin_start = h_Yields_Phi_SM[j][l]->FindBin(ParFit_Yields[j][l][0]-2.5*ParFit_Yields[j][l][1]);
      Int_t bin_stop  = h_Yields_Phi_SM[j][l]->FindBin(ParFit_Yields[j][l][0]+2.5*ParFit_Yields[j][l][1]);
      for(Int_t bin = bin_start; bin <= bin_stop; bin++)
      {
	Counts_Yields[j][l] += h_Yields_Phi_SM[j][l]->GetBinContent(bin);
      }
    }
  }

  // get resolution
  TString input_res = Form("/project/projectdirs/star/sss/OutPut/AuAu15GeV/Resolution/file_%s_Resolution.root",Energy[energy].Data());
  TFile *input = TFile::Open(input_res.Data());
  TProfile *p_res2[EtaGap_total];
  TProfile *p_res3[EtaGap_total];
  Double_t mean_res_2_phi[Centrality_total][EtaGap_total];
  Double_t mean_res_3_phi[Centrality_total][EtaGap_total];
  for(Int_t l = EtaGap_start; l < EtaGap_stop; l++)
  {
    TString ProName;
    ProName = Form("Res2_EtaGap_%d_EP",l);
    p_res2[l] = (TProfile*)input->Get(ProName.Data());
    ProName = Form("Res3_EtaGap_%d_EP",l);
    p_res3[l] = (TProfile*)input->Get(ProName.Data());
  }

  for(Int_t j = Centrality_start; j < Centrality_stop; j++)
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      mean_res_2_phi[j][l] = 0.0;
      mean_res_3_phi[j][l] = 0.0;

      Double_t res_2[9];
      Float_t phi2_total = 0.0;
      for(Int_t cent = cent_low[j]; cent <= cent_up[j]; cent++)
      {
	res_2[cent] = TMath::Sqrt(p_res2[l]->GetBinContent(cent+1));
	phi2_total += Counts_Yields[cent][l];
      }
      for(Int_t cent = cent_low[j]; cent <= cent_up[j]; cent++)
      {
	mean_res_2_phi[j][l] += Counts_Yields[cent][l]/(res_2[cent]*phi2_total);
      }
      cout <<  "centality_bin = " << j << ", eta_bin = " << l << ", mean_res_2_phi = " << mean_res_2_phi[j][l] << endl;

      Float_t phi3_total = 0.0;
      Double_t res_3[9];
      for(Int_t i = 0; i < 9; i++)
      {
	res_3[i] = -999.9;
      }
      for(Int_t cent = cent_low[j]; cent <= cent_up[j]; cent++)
      {
	if(p_res3[l]->GetBinContent(cent+1) > 0)
	{
	  res_3[cent] = TMath::Sqrt(p_res3[l]->GetBinContent(cent+1));
	  phi3_total += Counts_Yields[cent][l];
	}
      }
      for(Int_t cent = cent_low[j]; cent <= cent_up[j]; cent++)
      {
	mean_res_3_phi[j][l] += Counts_Yields[cent][l]/(res_3[cent]*phi3_total);
      }
      cout <<  "centality_bin = " << j << ", eta_bin = " << l << ", mean_res_3_phi = " << mean_res_3_phi[j][l] << endl;
    }
  }

  // scale raw v2 and v3 with resoltuion correction
  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
	h_flow_2[j][l]->Scale(mean_res_2_phi[j][l]);
	h_flow_3[j][l]->Scale(mean_res_3_phi[j][l]);
	h_flow_2_bw[j][l]->Scale(mean_res_2_phi[j][l]);
	h_flow_3_bw[j][l]->Scale(mean_res_3_phi[j][l]);
    }
  }

  // save v2 and v3 into TGraphAsymmErrors and into file
  TGraphAsymmErrors *g_flow_2[Centrality_total][EtaGap_total];
  TGraphAsymmErrors *g_flow_3[Centrality_total][EtaGap_total];
  TGraphAsymmErrors *g_flow_2_bw[Centrality_total][EtaGap_total];
  TGraphAsymmErrors *g_flow_3_bw[Centrality_total][EtaGap_total];
  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString g_Name_2 = Form("g_Centrality_%d_EtaGap_%d_2nd",j,l);
      g_flow_2[j][l] = new TGraphAsymmErrors();
      g_flow_2[j][l]->SetName(g_Name_2.Data());
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	Int_t bin_center = h_flow_2[j][l]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
	Float_t bin_content = h_flow_2[j][l]->GetBinContent(bin_center);
	Float_t bin_error   = h_flow_2[j][l]->GetBinError(bin_center);
	g_flow_2[j][l]->SetPoint(i,(pt_low_phi[i]+pt_up_phi[i])/2.0,bin_content);
	g_flow_2[j][l]->SetPointError(i,0.0,0.0,bin_error,bin_error);
      }

      TString g_Name_3 = Form("g_Centrality_%d_EtaGap_%d_3rd",j,l);
      g_flow_3[j][l] = new TGraphAsymmErrors();
      g_flow_3[j][l]->SetName(g_Name_3.Data());
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	Int_t bin_center = h_flow_3[j][l]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
	Float_t bin_content = h_flow_3[j][l]->GetBinContent(bin_center);
	Float_t bin_error   = h_flow_3[j][l]->GetBinError(bin_center);
	g_flow_3[j][l]->SetPoint(i,(pt_low_phi[i]+pt_up_phi[i])/2.0,bin_content);
	g_flow_3[j][l]->SetPointError(i,0.0,0.0,bin_error,bin_error);
      }

      TString g_Name_2_bw = Form("g_Centrality_%d_EtaGap_%d_2nd_bw",j,l);
      g_flow_2_bw[j][l] = new TGraphAsymmErrors();
      g_flow_2_bw[j][l]->SetName(g_Name_2_bw.Data());
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	Int_t bin_center = h_flow_2_bw[j][l]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
	Float_t bin_content = h_flow_2_bw[j][l]->GetBinContent(bin_center);
	Float_t bin_error   = h_flow_2_bw[j][l]->GetBinError(bin_center);
	g_flow_2_bw[j][l]->SetPoint(i,(pt_low_phi[i]+pt_up_phi[i])/2.0+0.05,bin_content);
	g_flow_2_bw[j][l]->SetPointError(i,0.0,0.0,bin_error,bin_error);
      }

      TString g_Name_3_bw = Form("g_Centrality_%d_EtaGap_%d_3rd_bw",j,l);
      g_flow_3_bw[j][l] = new TGraphAsymmErrors();
      g_flow_3_bw[j][l]->SetName(g_Name_3_bw.Data());
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	Int_t bin_center = h_flow_3_bw[j][l]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
	Float_t bin_content = h_flow_3_bw[j][l]->GetBinContent(bin_center);
	Float_t bin_error   = h_flow_3_bw[j][l]->GetBinError(bin_center);
	g_flow_3_bw[j][l]->SetPoint(i,(pt_low_phi[i]+pt_up_phi[i])/2.0+0.05,bin_content);
	g_flow_3_bw[j][l]->SetPointError(i,0.0,0.0,bin_error,bin_error);
      }
    }
  }

  TFile *File_OutPut = new TFile("/project/projectdirs/star/sss/OutPut/AuAu15GeV/Phi/flow_Phi/flow/Flow_Phi.root","RECREATE");
  File_OutPut->cd();
  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      g_flow_2[j][l]->Write();
      g_flow_3[j][l]->Write();
      g_flow_2_bw[j][l]->Write();
      g_flow_3_bw[j][l]->Write();
    }
  }
  File_OutPut->Close();
//  File_SE->Close();
//  File_ME->Close();
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	TString CanName;

	CanName = Form("./figures_phi/c2_phi_Psi_pt_%d_Centrality_%d_EtaGap_%d.gif",i,j,l);
	c2_phi_Psi[i][j][l]->SaveAs(CanName.Data());
	CanName = Form("./figures_phi/c3_phi_Psi_pt_%d_Centrality_%d_EtaGap_%d.gif",i,j,l);
	c3_phi_Psi[i][j][l]->SaveAs(CanName.Data());
	CanName = Form("./figures_phi/c2_phi_Psi_pt_%d_Centrality_%d_EtaGap_%d_bw.gif",i,j,l);
	c2_phi_Psi_bw[i][j][l]->SaveAs(CanName.Data());
	CanName = Form("./figures_phi/c3_phi_Psi_pt_%d_Centrality_%d_EtaGap_%d_bw.gif",i,j,l);
	c3_phi_Psi_bw[i][j][l]->SaveAs(CanName.Data());
      }
    }
  }
  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString CanName;

      CanName = Form("./figures_phi/c_phi_Psi2_Centrality_%d_EtaGap_%d.gif",j,l);
      c_phi_Psi2[j][l]->SaveAs(CanName.Data());
      CanName = Form("./figures_phi/c_phi_Psi3_Centrality_%d_EtaGap_%d.gif",j,l);
      c_phi_Psi3[j][l]->SaveAs(CanName.Data());
    }
  }

  // draw poster fit figure
  TCanvas *c_play = new TCanvas("c_play","c_play",1400,10,900,900);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetTicks(1,1);
  c_play->SetGrid(0,0);
  c_play->SetFillColor(0);
  c_play->SetFillStyle(4000);
  c_play->SetFrameFillColor(18);
//  h_Counts3[10][0][0]->SetTitle("");
  h_Counts3[2][0][0]->SetStats(0);
  h_Counts3[2][0][0]->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  h_Counts3[2][0][0]->GetYaxis()->SetTitle("dN/d(#phi-#Psi_{3})");
  h_Counts3[2][0][0]->GetXaxis()->SetTitleSize(0.06);
  h_Counts3[2][0][0]->GetYaxis()->SetTitleSize(0.06);
  h_Counts3[2][0][0]->GetXaxis()->SetTitleOffset(0.9);
  h_Counts3[2][0][0]->GetXaxis()->CenterTitle();
  h_Counts3[2][0][0]->GetYaxis()->CenterTitle();
  h_Counts3[2][0][0]->GetXaxis()->SetTitleColor(0);
  h_Counts3[2][0][0]->GetYaxis()->SetTitleColor(0);
  h_Counts3[2][0][0]->SetNdivisions(505,"X");
  h_Counts3[2][0][0]->SetNdivisions(505,"Y");
  h_Counts3[2][0][0]->GetXaxis()->SetLabelSize(0.04);
  h_Counts3[2][0][0]->GetYaxis()->SetLabelSize(0.04);
  h_Counts3[2][0][0]->GetXaxis()->SetLabelColor(0);
  h_Counts3[2][0][0]->GetYaxis()->SetLabelColor(0);
  h_Counts3[2][0][0]->Draw("pE");
  TFile *file_counts = new TFile("/project/projectdirs/star/sss/OutPut/AuAu15GeV/Phi/flow_Phi/flow/Counts.root","RECREATE");
  file_counts->cd();
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	h_Counts2[i][j][l]->Write();
	h_Counts3[i][j][l]->Write();
      }
    }
  }
  file_counts->Close();
}
