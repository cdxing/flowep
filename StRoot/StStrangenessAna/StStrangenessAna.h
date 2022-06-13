#ifndef StStrangenessAna_h
#define StStrangenessAna_h

#include "TObject.h"
#include "TString.h"

class StRefMultCorr;
class TFile;
class TChain;
class TTree;
class StAlexPhiMesonEvent;
class StAlexPhiMesonTrack;
class StStrangenessCorr;
class StStrangenessCut;
class StStrangenessHistoManger;
class StRunIdEventsDb;
class StV0Event;
class StV0Track;

class StStrangenessAna : public TObject
{
  public:
    //StStrangenessAna(Int_t energy, Int_t X_flag, Int_t List, Long64_t start_event, Long64_t stop_event, Int_t mode); // energy: 0 for 200GeV, 1 for 39GeV | X_flag: 0 for Same Event, 1 for Mixed Event | List: number of list to use | mode: 0 for phi, 1 for Lambda, 2 for anti-Lambda, 3 for K0s
    StStrangenessAna(Int_t energy, Int_t X_flag,const Char_t *inputFile, char* jobid, Long64_t start_event, Long64_t stop_event, Int_t mode); // energy: 0 for 200GeV, 1 for 39GeV | X_flag: 0 for Same Event, 1 for Mixed Event | List: number of list to use | mode: 0 for phi, 1 for Lambda, 2 for anti-Lambda, 3 for K0s
    ~StStrangenessAna();

    void setInputDir(const TString inputdir);
    void setOutputfile(const TString outputfile);
    void setSEList(const TString iSEList);
    void setMEList(const TString iMEList);
    void setStopEvent_SE(const Long64_t StopEvent_SE);
    void setStartEvent_SE(const Long64_t StartEvent_SE);
    void setStopEvent_ME(const Long64_t StopEvent_ME);
    void setStartEvent_ME(const Long64_t StartEvent_ME);

    void Init();
    void InitSE();
    void InitME();
    void Make();
    void MakePhiSE();
    void MakePhiME();
    void Finish();

    int GetRunIndex(int runId);

  private:
    TString mInputdir;
    TString mInputFile;
    TString mOutputfile;
    TString mJobId;
    TString mSEList;
    TString mMEList;
    Long64_t mStopEvent_SE;
    Long64_t mStartEvent_SE;
    Long64_t mStopEvent_ME;
    Long64_t mStartEvent_ME;

    Long64_t mStart_Event;
    Long64_t mStop_Event;

    TFile *mFile_OutPut;
    //TChain *mInPut_SE;
    TTree *mInPut_SE;
    //TChain *mInPut_ME;
    TTree *mInPut_ME;
    Int_t mEnergy;
    Int_t mX_flag; // 0 for Same Event, 1 for Mixed Event
    Int_t mList;
    Int_t mMode; // 0 for phi, 1 for Lambda, 2 for anti-Lambda, 3 for K0s
    StAlexPhiMesonEvent *mXuPhiMeson_event;
    StAlexPhiMesonTrack *mXuPhiMeson_track;
    StStrangenessCorr *mStrangenessCorr;
    StStrangenessCut *mStrangenessCut;
    StStrangenessHistoManger *mStrangenessHistoManger;
    StRunIdEventsDb *mRunIdEventsDb;

    static StRefMultCorr *mRefMultCorr;
    static Int_t mSE_input_flag;
    static Int_t mME_input_flag;
    static char* XUV0_EVENT_TREE;
    static char* XUV0_EVENT_BRANCH;

  ClassDef(StStrangenessAna,1)
};
#endif
