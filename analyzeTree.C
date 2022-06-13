#ifdef __CINT__ 
#define ROOT_GuiTypes 0 
#endif
#include <TSystem>
#include <TFile>

//void Analysis(const Int_t energy = 3, const Int_t X_flag = 0, const Int_t List = 1, const Long64_t start_event = 1000, const Long64_t stop_event = 10024, const Int_t mode = 0)

void analyzeTree(const Int_t energy = 3, const Int_t X_flag = 0, const Char_t *inputFile="test.list", char *jobid = "test", const Long64_t start_event = 1000, const Long64_t stop_event = 10024, const Int_t mode = 0)
{
	cout << "test 0"<< endl;
  // energy: 0 for 200GeV, 1 for 39GeV, 2 for 27 GeV, 3 for 15 GeV
  // X_flag: 0 for Same Event, 1 for Mixed Event
  // List: different number for different TTree list
  // mode: 0 for phi meson, 1 for Lambda, 2 for anti-Lambda, 3 for K0s
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  // root4star -b -q PhiFlow.C\(3,0,0,1000,1024,1\)

  gSystem->Load("StAlexPhiMesonEvent");
  //gSystem->Load("StV0Event");
  gSystem->Load("StStrangenessAna");
  gSystem->Load("StRunIdEventsDb");
  gSystem->Load("StRefMultCorr");

  cout << "All libraries are loaded!!!!" << endl;

  //**************************************************************************************

  cout << "Start to Read Trees!" << endl;

  //StStrangenessAna *mStrangenessAna = new StStrangenessAna(energy,X_flag,List,start_event,stop_event,mode);
  StStrangenessAna *mStrangenessAna = new StStrangenessAna(energy,X_flag,inputFile,jobid,start_event,stop_event,mode);
  mStrangenessAna->Init();
  mStrangenessAna->Make();
  mStrangenessAna->Finish();

//  cout << "End of the Calculation!!" << endl;
}
