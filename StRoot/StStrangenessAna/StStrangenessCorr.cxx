#include "StStrangenessCorr.h"
#include "StStrangenessCons.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TMath.h"
#include "TFile.h"
#include "StMessMgr.h"

ClassImp(StStrangenessCorr)

TString StStrangenessCorr::mVStr[2] = {"pos","neg"};
TString StStrangenessCorr::mOrder[2] = {"2nd","3rd"};
TString StStrangenessCorr::mMethod[2] = {"EP","SP"};
//--------------------------------------------------------
StStrangenessCorr::StStrangenessCorr()
{
}

StStrangenessCorr::~StStrangenessCorr()
{
}
//--------------------------------------------------------
// ReCenter Correction
void StStrangenessCorr::InitReCenterCorrection(Int_t mEnergy)
{
  //TString InPutFile_ReCenter = Form("/star/u/slan/pwg/fastoffline/7p7gev/recenterpar/file_%s_ReCenterPar.root",Strangeness::Energy[mEnergy].Data());
  TString InPutFile_ReCenter = Form("/star/data01/pwg/dchen/Ana/19p6GeV/shift/19gev_recenterpar.root");

  mInPutFile_ReCenter = TFile::Open(InPutFile_ReCenter.Data());
}

// get ReCenter Parameter
TVector2 StStrangenessCorr::getReCenterPar_East(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap, Int_t method)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  //TString ProName_x = Form("qx_%s_Vertex_%s_EtaGap_%d_East_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());
  //TString ProName_y = Form("qy_%s_Vertex_%s_EtaGap_%d_East_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());
  TString ProName_x = Form("p_mq2x_tpc_AB_sub_1");
  TString ProName_y = Form("p_mq2y_tpc_AB_sub_1");

  TProfile2D *p_x = (TProfile2D*)mInPutFile_ReCenter->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile_ReCenter->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

TVector2 StStrangenessCorr::getReCenterPar_West(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap, Int_t method)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  //TString ProName_x = Form("qx_%s_Vertex_%s_EtaGap_%d_West_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());
  //TString ProName_y = Form("qy_%s_Vertex_%s_EtaGap_%d_West_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());
  TString ProName_x = Form("p_mq2x_tpc_AB_sub_1");
  TString ProName_y = Form("p_mq2y_tpc_AB_sub_1");

  TProfile2D *p_x = (TProfile2D*)mInPutFile_ReCenter->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile_ReCenter->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

// calculate qVector of Track
TVector2 StStrangenessCorr::calq2Vector(TLorentzVector lTrack)
{
  const Float_t phi = lTrack.Phi();
  TVector2 q2Vector(0.0,0.0);

  const Float_t q2x = TMath::Cos(2.0*phi);
  const Float_t q2y = TMath::Sin(2.0*phi);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StStrangenessCorr::calq3Vector(TLorentzVector lTrack)
{
  const Float_t phi = lTrack.Phi();
  TVector2 q3Vector(0.0,0.0);

  const Float_t q3x = TMath::Cos(3.0*phi);
  const Float_t q3y = TMath::Sin(3.0*phi);
  q3Vector.Set(q3x,q3y);

  return q3Vector;
}

Float_t StStrangenessCorr::getWeight(TLorentzVector lTrack)
{
  Float_t pt = lTrack.Perp();
  Float_t w;
  if(pt <= Strangeness::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > Strangeness::mPrimPtWeight)
  {
    w = Strangeness::mPrimPtWeight;
  }

  return w;
}
//--------------------------------------------------------
// Shift Correction
void StStrangenessCorr::InitShiftCorrection(Int_t mEnergy)
{
  //TString InPutFile_Shift = Form("/star/u/slan/pwg/fastoffline/7p7gev/recenter/production/out_Shift/shift.root");
  TString InPutFile_Shift = Form("/star/data01/pwg/dchen/Ana/19p6GeV/shift/19gev_shiftpar.root");


  mInPutFile_Shift = TFile::Open(InPutFile_Shift.Data());

  //TString InPutFile_Res = Form("/star/u/slan/pwg/fastoffline/7p7gev/resolution/out_Resolution/Resolution.root");
  //TString InPutFile_Res = Form("/star/u/slan/pwg/fastoffline/7p7gev/resolution/out_Resolution/Resolution.root");
  TString InPutFile_Res = Form("/star/data01/pwg/dchen/Ana/19p6GeV/shift/19gev_shifted.root");
  mInPutFile_Res = TFile::Open(InPutFile_Res.Data());
}

bool StStrangenessCorr::passTrackNumCut(Int_t NumTrackEast, Int_t NumTrackWest)
{
  if(!(NumTrackEast > Strangeness::mTrackMin && NumTrackWest > Strangeness::mTrackMin))
  {
    return kFALSE;
  }

  return kTRUE;
}

Float_t StStrangenessCorr::AngleShift(Float_t Psi_raw, Float_t order)
{
  Float_t Psi_Corr = Psi_raw;
  if(Psi_raw > TMath::Pi()/order)
  {
    Psi_Corr = Psi_raw - 2.0*TMath::Pi()/order;
  }
  if(Psi_raw < -1.0*TMath::Pi()/order)
  {
    Psi_Corr = Psi_raw + 2.0*TMath::Pi()/order;
  }

  return Psi_Corr;
}

// calculate EP angle after Shift Correction EPD
// 1st
Float_t StStrangenessCorr::calShiftAngle1EPDfull_EP(TVector2 Q1Vector_EPD, Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(Q1Vector_EPD.Y(),Q1Vector_EPD.X());
  Float_t mean_sin[20], mean_cos[20];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  //for(Int_t k = 0; k < 5; k++) // Shift Order loop
  for(Int_t k=0; k<20; k++)
  {
    TString ProName_cos, ProName_sin;
    TProfile3D *p_cos, *p_sin;

    //ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    ProName_cos = Form("EPshiftpar_epd_ABCD_wt_full_cos");

    p_cos = (TProfile3D*)mInPutFile_Shift->Get(ProName_cos.Data());
    //mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));
    mean_cos[k] = p_cos->GetBinContent(Cent9+1, k+1, runIndex+1);
    	//cout << "SE test k  "<< k <<  endl;

    //ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    ProName_sin = Form("EPshiftpar_epd_ABCD_wt_full_sin");
    	//cout << "p_cos =" << p_cos <<  endl;
    	//cout << "SE test kk  "<< k <<  endl;
    p_sin = (TProfile3D*)mInPutFile_Shift->Get(ProName_sin.Data());
    	//cout << "SE test kkk  "<< k <<  endl;
    //mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));
    	//cout << "p_sin =" << p_sin <<  endl;
    	//cout << "p_sin->GetBinContent(Cent9+1, k+1, runIndex+1)" << p_sin->GetBinContent(Cent9+1, k+1, runIndex+1) <<  endl;
    mean_sin[k] = p_sin->GetBinContent(Cent9+1, k+1, runIndex+1);
    	//cout << "SE test kkkk  "<< k <<  endl;
               // implement shift correction
    float tmp = (float)(1*(k+1));
    float sinAve = p_sin->GetBinContent(Cent9+1, k+1, runIndex+1);
    float cosAve = p_cos->GetBinContent(Cent9+1, k+1, runIndex+1);
    delta_Psi +=
    2.0*(cosAve*TMath::Sin(tmp* Psi_ReCenter) - sinAve*TMath::Cos(tmp* Psi_ReCenter))/tmp;

    //delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(Strangeness::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(Strangeness::mShiftOrder2[k]*Psi_ReCenter));
    	//cout << "SE test kkkkk  "<< k <<  endl;
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,1.0);

  return Psi_Shift;
}

// calculate EP angle after Shift Correction
// 2nd
Float_t StStrangenessCorr::calShiftAngle2East_EP(TVector2 Q2Vector_East, Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(Q2Vector_East.Y(),Q2Vector_East.X())/2.0;
  Float_t mean_sin[20], mean_cos[20];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  //for(Int_t k = 0; k < 5; k++) // Shift Order loop
  for(Int_t k=0; k<20; k++)
  {
    TString ProName_cos, ProName_sin;
    TProfile3D *p_cos, *p_sin;

    //ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    /*ProName_cos = Form("EPshiftpar_tpc_AB_sub_1_cos");
	    
    p_cos = (TProfile3D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    //ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    ProName_sin = Form("EPshiftpar_tpc_AB_sub_1_sin");
    p_sin = (TProfile3D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(Strangeness::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(Strangeness::mShiftOrder2[k]*Psi_ReCenter));*/

    //ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    ProName_cos = Form("EPshiftpar_tpc_AB_sub_1_cos");

    p_cos = (TProfile3D*)mInPutFile_Shift->Get(ProName_cos.Data());
    //mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));
    mean_cos[k] = p_cos->GetBinContent(Cent9+1, k+1, runIndex+1);
    	//cout << "SE test k  "<< k <<  endl;

    //ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    ProName_sin = Form("EPshiftpar_tpc_AB_sub_1_sin"); //cout << "p_cos =" << p_cos <<  endl; //cout << "SE test kk  "<< k <<  endl;
    p_sin = (TProfile3D*)mInPutFile_Shift->Get(ProName_sin.Data());
    	//cout << "SE test kkk  "<< k <<  endl;
    //mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));
    	//cout << "p_sin =" << p_sin <<  endl;
    	//cout << "p_sin->GetBinContent(Cent9+1, k+1, runIndex+1)" << p_sin->GetBinContent(Cent9+1, k+1, runIndex+1) <<  endl;
    mean_sin[k] = p_sin->GetBinContent(Cent9+1, k+1, runIndex+1);
    	//cout << "SE test kkkk  "<< k <<  endl;

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(Strangeness::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(Strangeness::mShiftOrder2[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,2.0);

  return Psi_Shift;
}

Float_t StStrangenessCorr::calShiftAngle2West_EP(TVector2 Q2Vector_West, Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(Q2Vector_West.Y(),Q2Vector_West.X())/2.0;
  Float_t mean_sin[20], mean_cos[20];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  //for(Int_t k = 0; k < 5; k++) // Shift Order loop
  for(Int_t k=0; k<20; k++)
  {
    TString ProName_cos, ProName_sin;
    TProfile3D *p_cos, *p_sin;

    //ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    ProName_cos = Form("EPshiftpar_tpc_AB_sub_2_cos");

    p_cos = (TProfile3D*)mInPutFile_Shift->Get(ProName_cos.Data());
    //mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));
    mean_cos[k] = p_cos->GetBinContent(Cent9+1, k+1, runIndex+1);
    	//cout << "SE test k  "<< k <<  endl;

    //ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    ProName_sin = Form("EPshiftpar_tpc_AB_sub_2_sin");
    	//cout << "p_cos =" << p_cos <<  endl;
    	//cout << "SE test kk  "<< k <<  endl;
    p_sin = (TProfile3D*)mInPutFile_Shift->Get(ProName_sin.Data());
    	//cout << "SE test kkk  "<< k <<  endl;
    //mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));
    	//cout << "p_sin =" << p_sin <<  endl;
    	//cout << "p_sin->GetBinContent(Cent9+1, k+1, runIndex+1)" << p_sin->GetBinContent(Cent9+1, k+1, runIndex+1) <<  endl;
    mean_sin[k] = p_sin->GetBinContent(Cent9+1, k+1, runIndex+1);
    	//cout << "SE test kkkk  "<< k <<  endl;

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(Strangeness::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(Strangeness::mShiftOrder2[k]*Psi_ReCenter));
    	//cout << "SE test kkkkk  "<< k <<  endl;
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,2.0);

  return Psi_Shift;
}

// 3rd 
Float_t StStrangenessCorr::calShiftAngle3East_EP(TVector2 Q3Vector_East, Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(Q3Vector_East.Y(),Q3Vector_East.X())/3.0;
  Float_t mean_sin[20], mean_cos[20];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  //for(Int_t k = 0; k < 5; k++) // Shift Order loop
  for(Int_t k=0; k<20; k++)
  {
    TString ProName_cos, ProName_sin;
    TProfile3D *p_cos, *p_sin;

    //ProName_cos = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    ProName_cos = Form("EPshiftpar_tpc_AB_sub_1_cos");
    p_cos = (TProfile3D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    //ProName_sin = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    ProName_sin = Form("EPshiftpar_tpc_AB_sub_1_sin");
    p_sin = (TProfile3D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(Strangeness::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(Strangeness::mShiftOrder3[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,3.0);

  return Psi_Shift;
}

Float_t StStrangenessCorr::calShiftAngle3West_EP(TVector2 Q3Vector_West, Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(Q3Vector_West.Y(),Q3Vector_West.X())/3.0;
  Float_t mean_sin[20], mean_cos[20];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  //for(Int_t k = 0; k < 5; k++) // Shift Order loop
  for(Int_t k=0; k<20; k++)
  {
    TString ProName_cos, ProName_sin;
    TProfile3D *p_cos, *p_sin;

    //ProName_cos = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    ProName_cos = Form("EPshiftpar_tpc_AB_sub_2_cos");
    p_cos = (TProfile3D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    //ProName_sin = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    ProName_sin = Form("EPshiftpar_tpc_AB_sub_2_sin");
    p_sin = (TProfile3D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(Strangeness::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(Strangeness::mShiftOrder3[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw,3.0);

  return Psi_Shift;
}
//--------------------------------------------------------
// Resolution Correction
Float_t StStrangenessCorr::getResolution1_EP(Int_t Cent9, Int_t eta_gap)
{
  //TString ProName = Form("Res2_EtaGap_%d_EP",eta_gap);
  TString ProName = Form("p_r1_epd_ABCD_wt_sub0_1");
  TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
  Float_t Res_raw = p_res2->GetBinContent(p_res2->FindBin(Cent9));
  if(Res_raw <= 0)
  {
    return -999.9;
  }
  else
  {
    Float_t Res = TMath::Sqrt(Res_raw) * TMath::Sqrt(2); // full event plane resolution
    return Res;
  }
}

Float_t StStrangenessCorr::getResolution2_EP(Int_t Cent9, Int_t eta_gap)
{
  //TString ProName = Form("Res2_EtaGap_%d_EP",eta_gap);
  TString ProName = Form("p_r2_tpc_AB_sub0_1");
  TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
  Float_t Res_raw = p_res2->GetBinContent(p_res2->FindBin(Cent9));
  if(Res_raw <= 0)
  {
    return -999.9;
  }
  else
  {
    Float_t Res = TMath::Sqrt(Res_raw);
    return Res;
  }
}

Float_t StStrangenessCorr::getResolution3_EP(Int_t Cent9, Int_t eta_gap)
{
  //TString ProName = Form("Res3_EtaGap_%d_EP",eta_gap);
  TString ProName = Form("p_r2_tpc_AB_sub0_1");
  TProfile *p_res3 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
  Float_t Res_raw = p_res3->GetBinContent(p_res3->FindBin(Cent9));
  if(Res_raw <= 0)
  {
    return -999.9;
  }
  else
  {
    Float_t Res = TMath::Sqrt(Res_raw);
    return Res;
  }
}
