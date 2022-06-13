#include "StStrangenessCut.h"
#include "StStrangenessCons.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StMessMgr.h"

ClassImp(StStrangenessCut)

//---------------------------------------------------------------------------------

StStrangenessCut::StStrangenessCut(Int_t energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StStrangenessCut::~StStrangenessCut()
{
}

//---------------------------------------------------------------------------------
bool StStrangenessCut::passTrackEP(TLorentzVector lTrack, Float_t dca)
{
  // dca cut for event plane reconstruction
  //if(fabs(dca) > Strangeness::mDcaEPMax[mEnergy])
    if(dca > 1.0)   //  should be always be less than 1
    {
        return kFALSE;
    }

    // pt cut 0.2 - 2.0 GeV/c
  Float_t pt = lTrack.Perp();
  if(!(pt > Strangeness::mPrimPtMin[mEnergy] && pt < Strangeness::mPrimPtMax))
  {
    return kFALSE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------
bool StStrangenessCut::passTrackEtaEast(TLorentzVector lTrack, Int_t j) // neg || j = different eta_gap
{
  Float_t eta = lTrack.Eta();
  
  // eta cut
  // eta_gap between two sub event plane is 2*mEta_Gap[i]
  if(!(eta > -1.0*Strangeness::mEtaMax && eta < -1.0*Strangeness::mEta_Gap[j]))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StStrangenessCut::passTrackEtaWest(TLorentzVector lTrack, Int_t j) // pos || j = different eta_gap
{
  Float_t eta = lTrack.Eta();

  // eta cut
  // eta_gap between two sub event plane is 2*mEta_Gap[i]
  if(!(eta > Strangeness::mEta_Gap[j] && eta < Strangeness::mEtaMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StStrangenessCut::passPhiEtaEast(TLorentzVector lTrack) // neg
{
  Float_t eta = lTrack.Eta();
  
  // eta cut
  // eta_gap between phi meson and sub event plane is mEta_Gap[i]
  if(!(eta > -1.0*Strangeness::mEtaMax && eta < 0.0))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StStrangenessCut::passPhiEtaWest(TLorentzVector lTrack) // pos
{
  Float_t eta = lTrack.Eta();

  // eta cut
  // eta_gap between two sub event plane is 2*mEta_Gap[i]
  if(!(eta > 0.0 && eta < Strangeness::mEtaMax))
  {
    return kFALSE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------
