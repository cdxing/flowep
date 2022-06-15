#include "StStrangenessAna.h"
#include "StStrangenessCons.h"
#include "StStrangenessCut.h"
#include "StStrangenessCorr.h"
#include "StStrangenessHistoManger.h"
#include "StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.h"
//#include "StRoot/StRefMultCorr/StRefMultCorr.h"
//#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StThreeVectorF.hh"
#include "StMessMgr.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include <fstream>
#include "StRoot/Run/run.h"

ClassImp(StStrangenessAna)

    //StRefMultCorr* StStrangenessAna::mRefMultCorr = NULL;   // shaowei
    Int_t StStrangenessAna::mSE_input_flag = 1;
    Int_t StStrangenessAna::mME_input_flag = 1;
    char* StStrangenessAna::XUV0_EVENT_TREE = NULL;
    char* StStrangenessAna::XUV0_EVENT_BRANCH = NULL;
    //----------------------------------------------------
//StStrangenessAna::StStrangenessAna(Int_t energy, Int_t X_flag, Int_t List, Long64_t start_event, Long64_t stop_event, Int_t mode)
StStrangenessAna::StStrangenessAna(Int_t energy,Int_t X_flag,const Char_t *inputFile, char* jobid, Long64_t start_event, Long64_t stop_event, Int_t mode)
{
    mEnergy = energy;
    mX_flag = X_flag;
    //mList = List;
    mInputFile =  inputFile;
    mJobId  = jobid;
    mStart_Event = start_event;
    mStop_Event = stop_event;
    mMode = mode;
    mStrangenessCorr = new StStrangenessCorr();
    mStrangenessCut = new StStrangenessCut(mEnergy);
    mStrangenessHistoManger = new StStrangenessHistoManger();
}

StStrangenessAna::~StStrangenessAna()
{
}
//----------------------------------------------------
// set Input/Output
void StStrangenessAna::setInputDir(const TString inputdir)
{
    mInputdir = inputdir.Copy();
    cout << "Input directory was set to: " << mInputdir.Data() << endl;
}
void StStrangenessAna::setOutputfile(const TString outputfile)
{
    mOutputfile = outputfile.Copy();
    cout << "Output file was set to: " << mOutputfile.Data() << endl;
}
void StStrangenessAna::setSEList(const TString iSEList)
{
    mSEList = iSEList.Copy();
    cout << "Same event list was set to: " << mSEList.Data() << endl;
}
void StStrangenessAna::setMEList(const TString iMEList)
{
    mMEList = iMEList.Copy();
    cout << "Mixed event list was set to: " << mMEList.Data() << endl;
}
void StStrangenessAna::setStopEvent_SE(const Long64_t StopEvent_SE)
{
    mStopEvent_SE = StopEvent_SE;
    cout << "nStopEvent_SE = " << mStopEvent_SE << endl;
}
void StStrangenessAna::setStartEvent_SE(const Long64_t StartEvent_SE)
{
    mStartEvent_SE = StartEvent_SE;
    cout << "nStartEvent_SE = " << mStartEvent_SE << endl;
}
void StStrangenessAna::setStopEvent_ME(const Long64_t StopEvent_ME)
{
    mStopEvent_ME = StopEvent_ME;
    cout << "nStopEvent_ME = " << mStopEvent_ME << endl;
}
void StStrangenessAna::setStartEvent_ME(const Long64_t StartEvent_ME)
{
    mStartEvent_ME = StartEvent_ME;
    cout << "nStartEvent_ME = " << mStartEvent_ME << endl;
}
//----------------------------------------------------
// initial functions
void StStrangenessAna::Init()
{

    mStrangenessCorr->InitReCenterCorrection(mEnergy);
    mStrangenessCorr->InitShiftCorrection(mEnergy);
    mStrangenessHistoManger->Init(mX_flag,mMode);

    //TString inputdir = Form("/star/u/slan/pwg/fastoffline/7p7gev/phitree/out_Phi/");
    //TString inputdir = Form("/star/data01/pwg/dchen/Ana/19p6GeV/parentparticle/production/");
    //TString inputdir = Form("/star/data01/pwg/dchen/Ana/19p6GeV/parentparticle_me/production/");
    TString inputdir = Form("/star/u/dchen/ana/19gev_2019/flowep/");
    setInputDir(inputdir);

    //const Int_t list_start = Strangeness::mList_Delta*mList + 1; // start list
    //const Int_t list_stop  = Strangeness::mList_Delta*(mList+1); // stop list

    if(mX_flag == 0)
    {
        //TString SEList = Form("/star/u/slan/pwg/fastoffline/7p7gev/phitree/List/Split_SE_%s_%d_%d.list",Strangeness::Energy[mEnergy].Data(),list_start,list_stop);
        //TString SEList = Form("/star/u/dchen/ana/19gev_2019/parentparticle/phitreelist/Split_SE_19GeV_%d_%d.list",list_start,list_stop);
        TString SEList = Form("/star/u/dchen/ana/19gev_2019/parentparticle/phitreelist/SElist");
        setSEList(SEList);
        cout << SEList << endl;

        //TString outputfile = Form("/star/u/slan/pwg/fastoffline/7p7gev/phiflow/flow_%s_SE/Yields_SE_%s_%d_%d.root",Strangeness::Partype[mMode].Data(),Strangeness::Energy[mEnergy].Data(),list_start,list_stop);
        TString outputfile = Form("/star/data01/pwg/dchen/Ana/19p6GeV/flowep/flow_%s_SE/Yields_SE_19GeV_%s.root",Strangeness::Partype[mMode].Data(),mJobId.Data());
        setOutputfile(outputfile);

        setStartEvent_SE(Long64_t(mStart_Event));
        setStopEvent_SE(Long64_t(mStop_Event));

        InitSE();
    }

    if(mX_flag == 1)
    {
        //TString MEList = Form("/star/u/slan/pwg/fastoffline/7p7gev/phitree/List/Split_ME_%s_%d_%d.list",Strangeness::Energy[mEnergy].Data(),list_start,list_stop);
        //TString MEList = Form("/star/u/dchen/ana/19gev_2019/parentparticle/phitreelist/Split_ME_19GeV_%d_%d.list",list_start,list_stop);
        TString MEList = Form("/star/u/dchen/ana/19gev_2019/parentparticle/phitreelist/MElist");
        setMEList(MEList);

        //TString outputfile = Form("/star/u/slan/pwg/fastoffline/7p7gev/phiflow/flow_%s_SE/Yields_SE_19GeV_%d_%d.root",Strangeness::Partype[mMode].Data(),list_start,list_stop);
        TString outputfile = Form("/star/data01/pwg/dchen/Ana/19p6GeV/flowep/flow_%s_ME/Yields_ME_19GeV_%s.root",Strangeness::Partype[mMode].Data(),mJobId.Data());
        //TString outputfile = Form("/star/u/slan/pwg/fastoffline/7p7gev/phiflow/flow_%s_ME/Yields_ME_%s_%d_%d.root",Strangeness::Partype[mMode].Data(),Strangeness::Energy[mEnergy].Data(),list_start,list_stop);
        setOutputfile(outputfile);

        setStartEvent_ME(Long64_t(mStart_Event));
        setStopEvent_ME(Long64_t(mStop_Event));

        InitME();
    }
}

// Initialize Same Event
void StStrangenessAna::InitSE()
{
    TString Notification = Form("Initializing parameters and input/output for %s Same Event",Strangeness::Partype[mMode].Data());
    cout << Notification.Data() << endl;
    mFile_OutPut = new TFile(mOutputfile.Data(),"RECREATE");

    XUV0_EVENT_TREE       = (char*)Strangeness::v0_tree[mMode].Data();
    XUV0_EVENT_BRANCH     = (char*)Strangeness::v0_branch[mMode].Data();

    //----------------------------------------------------------------------------------------------------
    // Same event input
    if (gSystem->AccessPathName(mInputFile)) { std::cout << "Error reading input file!" << std::endl;}
     TFile *inputFile = TFile::Open(mInputFile);
       if (!inputFile) { std::cout << "Input file could not be opened properly!" << std::endl;  }
    Long64_t entries_save = 0;
    if(inputFile)
    {
	  mInPut_SE = (TTree*)inputFile->Get(XUV0_EVENT_TREE);  

        Long64_t file_entries = mInPut_SE->GetEntries();
        cout << "File added to data chain: " << mInputFile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
        entries_save = file_entries;
    }
       if (!inputFile) { std::cout << "Input file could not be opened properly!" << std::endl;  }
    /*if (!mSEList.IsNull())   // if input file is ok
    {
        cout << "Open same event file list " << mSEList << endl;
        ifstream in(mSEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            mInPut_SE  = new TChain( XUV0_EVENT_TREE, XUV0_EVENT_TREE );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile;
                    addfile = str;
                    addfile = mInputdir+addfile;
                    mInPut_SE->AddFile(addfile.Data(),-1, XUV0_EVENT_TREE );
                    Long64_t file_entries = mInPut_SE->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
        }
        else
        {
            cout << "WARNING: SE file input is problemtic" << endl;
            mSE_input_flag = 0;
        }
    }*/

    // Set the input tree
    if (mSE_input_flag == 1 && !mInPut_SE->GetBranch( XUV0_EVENT_BRANCH ))
    {
        cerr << "ERROR: Could not find branch '"
            << XUV0_EVENT_BRANCH << "'in tree!" << endl;
    }

    if(mMode == 0) mXuPhiMeson_event = new StAlexPhiMesonEvent();

    if(mSE_input_flag == 1)
    {
        if(mMode == 0) mInPut_SE->SetBranchAddress( XUV0_EVENT_BRANCH, &mXuPhiMeson_event );

        Int_t num_events_SE = mInPut_SE->GetEntriesFast();
        cout << "Number of events in file(s) = " << num_events_SE << endl;
        if(mStartEvent_SE > num_events_SE) mStartEvent_SE = num_events_SE;
        if(mStopEvent_SE  > num_events_SE) mStopEvent_SE  = num_events_SE;

        cout << "New nStartEvent_SE = " << mStartEvent_SE << ", new nStopEvent_SE = " << mStopEvent_SE << endl;
    }
}

// Initialize Mixed Event
void StStrangenessAna::InitME()
{
    TString Notification = Form("Initializing parameters and input/output for %s Mixed Event",Strangeness::Partype[mMode].Data());
    cout << Notification.Data() << endl;
    cout << "test 0 " << endl;

    mFile_OutPut = new TFile(mOutputfile.Data(),"RECREATE");

    XUV0_EVENT_TREE       = (char*)Strangeness::v0_tree[mMode].Data();
    XUV0_EVENT_BRANCH     = (char*)Strangeness::v0_branch[mMode].Data();
    cout << "test 1 " << endl;

    //----------------------------------------------------------------------------------------------------
    // Mixed event input
    if (gSystem->AccessPathName(mInputFile)) { std::cout << "Error reading input file!" << std::endl;}
     TFile *inputFile = TFile::Open(mInputFile);
       if (!inputFile) { std::cout << "Input file could not be opened properly!" << std::endl;  }
    Long64_t entries_save = 0;
    cout << "test 2 " << endl;
    if(inputFile)
    {
    cout << "test 3 " << endl;
        //TString addfile;
    cout << "test 3.1 " << endl;
        //addfile = mInputFile;
    cout << "test 3.2 " << endl;
        //addfile = mInputdir+addfile;
    cout << "test 3.3 " << endl;
    //cout << "addfile.Data(): "<< addfile.Data() << endl;
        //mInPut_ME->AddFile(addfile.Data(),-1, XUV0_EVENT_TREE );
	  mInPut_ME = (TTree*)inputFile->Get(XUV0_EVENT_TREE);  

    cout << "test 3.4 " << endl;
        Long64_t file_entries = mInPut_ME->GetEntries();
    cout << "test 3.5 " << endl;
        cout << "File added to data chain: " << mInputFile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
    cout << "test 3.6 " << endl;
        entries_save = file_entries;
    cout << "test 3.7 " << endl;
    }
    /*if (!mMEList.IsNull())   // if input file is ok
    {
        cout << "Open mixed event file list " << mMEList << endl;
        ifstream in(mMEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            mInPut_ME  = new TChain( XUV0_EVENT_TREE, XUV0_EVENT_TREE );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile;
                    addfile = str;
                    addfile = mInputdir+addfile;
    			cout << "addfile.Data(): "<< addfile.Data() << endl;
                    mInPut_ME->AddFile(addfile.Data(),-1, XUV0_EVENT_TREE );
                    Long64_t file_entries = mInPut_ME->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
        }
        else
        {
            cout << "WARNING: ME file input is problemtic" << endl;
            mME_input_flag = 0;
        }
    }*/

    // Set the input tree
    if (mME_input_flag == 1 && !mInPut_ME->GetBranch( XUV0_EVENT_BRANCH ))
    {
        cerr << "ERROR: Could not find branch '"
            << XUV0_EVENT_BRANCH << "'in tree!" << endl;
    }
    cout << "test 3 " << endl;

    if(mMode == 0) mXuPhiMeson_event = new StAlexPhiMesonEvent();
    cout << "test 4 " << endl;

    if(mME_input_flag == 1)
    {
        if(mMode == 0) mInPut_ME->SetBranchAddress( XUV0_EVENT_BRANCH, &mXuPhiMeson_event );

        Int_t num_events_ME = mInPut_ME->GetEntriesFast();
        cout << "Number of events in file(s) = " << num_events_ME << endl;
        if(mStartEvent_ME > num_events_ME) mStartEvent_ME = num_events_ME;
        if(mStopEvent_ME  > num_events_ME) mStopEvent_ME  = num_events_ME;

        cout << "New nStartEvent_ME = " << mStartEvent_ME << ", new nStopEvent_ME = " << mStopEvent_ME << endl;
    }
    cout << "test 5 " << endl;
}
//----------------------------------------------------
void StStrangenessAna::Make()
{
    if(mX_flag == 0)
    {
        if(mMode == 0) MakePhiSE();
    cout << "test 5.0 " << endl;
    }

    if(mX_flag == 1)
    {
        if(mMode == 0) MakePhiME();
    cout << "test 5.1 " << endl;
    }
}


// loop phi meson Same Event
void StStrangenessAna::MakePhiSE()
{
    Long64_t start_event_use = mStartEvent_SE;
    Long64_t stop_event_use  = mStopEvent_SE;

    mInPut_SE->SetBranchAddress( XUV0_EVENT_BRANCH, &mXuPhiMeson_event );
    mInPut_SE->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

    // Initialise Event Head
    TVector3 PrimaryVertex(0.0,0.0,0.0);
    Int_t          RunId = 0;
    //Int_t          RefMult = 0;
    Int_t          Centrality = 0;
    Float_t        Reweight = 0;
    //Int_t          N_prim = 0;
    //Int_t          N_non_prim = 0;
    //Int_t          N_Tof_match = 0;
    //Float_t        ZDCx = 0.0; 
    //Float_t        BBCx = 0.0; 
    //Float_t        VzVpd = 0.0;
    Int_t          NumTrackUsed = 0;
    // ---------------------------------------QVector---------------------------------------------
    TVector2       Q1Full[4];
    TVector2       Q2East[4];
    TVector2       Q2West[4];
    TVector2       Q3East[4];
    TVector2       Q3West[4];
    // -----------------------------------Number of Tracks----------------------------------------
    Int_t          NumTrackEast[4];
    Int_t          NumTrackWest[4];
    for(Int_t j = 0; j < 4; j++)
    {
        Q1Full[j].Set(0.0,0.0);
        Q2East[j].Set(0.0,0.0);
        Q2West[j].Set(0.0,0.0);
        Q3East[j].Set(0.0,0.0);
        Q3West[j].Set(0.0,0.0);
        NumTrackEast[j] = 0;
        NumTrackWest[j] = 0;
    }

    for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
    {
        if (!mInPut_SE->GetEntry( counter )) // take the event -> information is stored in event
            break;  // end of data chunk
    	//cout << "SE test 0 " <<  endl;

        // get Event Header
        PrimaryVertex     = mXuPhiMeson_event->getPrimaryVertex();
        RunId             = mXuPhiMeson_event->getRunId();
        //RefMult           = mXuPhiMeson_event->getRefMult();
        Centrality        = mXuPhiMeson_event->getCentrality();
        Reweight          = mXuPhiMeson_event->getReweight();
        //N_prim            = mXuPhiMeson_event->getN_prim();
        //N_non_prim        = mXuPhiMeson_event->getN_non_prim();
        //N_Tof_match       = mXuPhiMeson_event->getN_Tof_match();
        //ZDCx              = mXuPhiMeson_event->getZDCx(); 
        //BBCx              = mXuPhiMeson_event->getBBCx(); 
        //VzVpd             = mXuPhiMeson_event->getVzVpd();
        NumTrackUsed      = mXuPhiMeson_event->getNumTracks();

        for(Int_t j = 0; j < 4; j++)
        {
            Q1Full[j]       = mXuPhiMeson_event->getQ1Full(j);
            Q2East[j]       = mXuPhiMeson_event->getQ2East(j);
            Q2West[j]       = mXuPhiMeson_event->getQ2West(j);
            Q3East[j]       = mXuPhiMeson_event->getQ3East(j);
            Q3West[j]       = mXuPhiMeson_event->getQ3West(j);
            NumTrackEast[j] = mXuPhiMeson_event->getNumTrackEast(j);
            NumTrackWest[j] = mXuPhiMeson_event->getNumTrackWest(j);
        }
        Float_t Psi1_Full = -999.9;
        Float_t Psi2_East_ltrack = -999.9, Psi2_West_ltrack = -999.9;

        // Initialise Track 
        Float_t m2A = -999.9;
        Float_t m2B = -999.9;
        Float_t nsA = -999.9;
        Float_t nsB = -999.9;
        Float_t dcaA = -999.9;
        Float_t dcaB = -999.9;
        TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
        TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
        Int_t flagA = -1;
        Int_t flagB = -1;


        //take the Vz cut |Vz|<40
        //if(fabs(PrimaryVertex.Z()) > 40.0)  continue;

        // vz sign
        Int_t vz_sign;
        if(PrimaryVertex.Z() > 0.0)
        {
            vz_sign = 0;
        }
        else
        {
            vz_sign = 1;
        }


        //float vx=PrimaryVertex.X();
        //float vy=PrimaryVertex.Y();
        //float vz=PrimaryVertex.Z();

        //if( (vx < 1.e-5 && vx > -1.e-5) &&
        //        (vy < 1.e-5 && vy > -1.e-5) &&
        //        (vz < 1.e-5 && vz > -1.e-5)  )
        //    continue;


        // Centrality
        //mRefMultCorr->init(RunId);
        //mRefMultCorr->initEvent(RefMult,PrimaryVertex.z());
        const Int_t cent9 = Centrality;
        const Double_t reweight = Reweight;
        if(cent9 > 8 || cent9 < 0) continue;
    	//cout << "SE test 1 " <<  endl;

        // runIndex
        // cout << runIndex << endl;
        const int runIndex = GetRunIndex(RunId);

        if (counter != 0  &&  counter % 1000 == 0)
            cout << "." << flush;
        if (counter != 0  &&  counter % 10000 == 0)
        {
            if((stop_event_use-start_event_use) > 0)
            {
                Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
                cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (strangeness_flow) " << flush;
            }
        }
        TVector2 Q1EPD_Full = Q1Full[0];
        //Psi1_Full = TMath::ATan2(Q1EPD_Full.Y(),Q1EPD_Full.X());
        Psi1_Full = mStrangenessCorr->calShiftAngle1EPDfull_EP(Q1EPD_Full,runIndex,cent9,vz_sign,0);
        mStrangenessHistoManger->FillPsiShiftFull_wt(Psi1_Full);

        // get Track Information
        for(Int_t j = 0; j < Strangeness::mEtaGap_total; j++)
        {
    	//cout << "SE test 2 " <<  endl;
            if(mStrangenessCorr->passTrackNumCut(NumTrackEast[j],NumTrackWest[j]))
            { // passTrackNumCut
    	//cout << "SE test 3 " <<  endl;
                for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
                {
    	//cout << "SE test 4 " <<  endl;
                    mXuPhiMeson_track = mXuPhiMeson_event->getTrack(nTracks);
                    m2A = mXuPhiMeson_track->getMass2A();
                    m2B = mXuPhiMeson_track->getMass2B();
                    nsA = mXuPhiMeson_track->getNSigKaonA();
                    nsB = mXuPhiMeson_track->getNSigKaonB();
                    dcaA = mXuPhiMeson_track->getDcaA();
                    dcaB = mXuPhiMeson_track->getDcaB();
                    lTrackA = mXuPhiMeson_track->getTrackA();
                    lTrackB = mXuPhiMeson_track->getTrackB();
                    flagA = mXuPhiMeson_track->getFlagA();
                    flagB = mXuPhiMeson_track->getFlagB();

                    Float_t pA = lTrackA.P();
                    Float_t pB = lTrackB.P();
                    TLorentzVector lTrack = lTrackA + lTrackB;
                    Float_t pt_lTrack = lTrack.Perp();

                    // apply additional PID cut to increase significance
                    if(
                            ((fabs(pA) <= 0.65 && m2A < -10) || (m2A > 0 && ((fabs(pA) < 1.5 && m2A > 0.16 && m2A < 0.36) || (fabs(pA) >= 1.5 && m2A > 0.125 && m2A < 0.36)) )) &&
                            ((fabs(pB) <= 0.65 && m2B < -10) || (m2B > 0 && ((fabs(pB) < 1.5 && m2B > 0.16 && m2B < 0.36) || (fabs(pB) >= 1.5 && m2B > 0.125 && m2B < 0.36)) )) &&
                            ((pt_lTrack) < 0.8 || ((pt_lTrack) >= 0.8 && ( (m2A > 0.16 && m2A < 0.36) || (m2B > 0.16 && m2B < 0.36)))) &&
                            (
                             ((m2A < -10 && nsA < 2.5 && nsA > -(1.5)) || (m2A > 0.16 && m2A < 0.36)) &&
                             ((m2B < -10 && nsB < 2.5 && nsB > -(1.5)) || (m2B > 0.16 && m2B < 0.36))
                            )
                      )
                    {
                        Float_t phi_lTrack = lTrack.Phi();
                        Float_t InvMass_lTrack = lTrack.M();
                        Float_t y_lTrack = lTrack.Rapidity();
                        mStrangenessHistoManger->FillAcc(pt_lTrack, y_lTrack, InvMass_lTrack);
    	//cout << "SE test 5 " <<  endl;


                        if(mStrangenessCut->passPhiEtaEast(lTrack)) // neg eta(east)
                        { // Below is West Only
                            TVector2 Q2Vector = Q2West[j];
                            TVector2 Q3Vector = Q3West[j];
    	//cout << "SE test 6 " <<  endl;
                            // subtract auto-correlation from pos eta(west) event plane
                            if(flagA == 0 && mStrangenessCut->passTrackEP(lTrackA,dcaA) && mStrangenessCut->passTrackEtaWest(lTrackA,j)) // trackA
                            {
    	//cout << "SE test 6.0 " <<  endl;
                                Float_t  w = mStrangenessCorr->getWeight(lTrackA);

                                TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lTrackA);
                                TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
                                Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

                                TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lTrackA);
                                TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
                                Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
                            }
                            if(flagB == 0 && mStrangenessCut->passTrackEP(lTrackB,dcaB) && mStrangenessCut->passTrackEtaWest(lTrackB,j)) // trackB
                            {
    	//cout << "SE test 6.1 " <<  endl;
                                Float_t  w = mStrangenessCorr->getWeight(lTrackB);

                                TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lTrackB);
                                TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
                                Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

                                TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lTrackB);
                                TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
                                Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
                            }
    	//cout << "SE test 6.2 " <<  endl;
                            Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
    	//cout << "SE test 6.3 " <<  endl;
                            Float_t Psi2_west = mStrangenessCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign,j);
			    Psi2_West_ltrack = Psi2_west ; 
    	//cout << "SE test 6.4 " <<  endl;
                            Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
    	//cout << "SE test 6.5 " <<  endl;
                            Float_t Psi3_west = mStrangenessCorr->calShiftAngle3West_EP(Q3Vector,runIndex,cent9,vz_sign,j);
    	//cout << "SE test 6.6 " <<  endl;
                            Float_t phi_Psi2 = phi_lTrack - Psi2_west;
                            Float_t phi_Psi3 = phi_lTrack - Psi3_west;

    	//cout << "SE test 6.7 " <<  endl;
                            mStrangenessHistoManger->Fill(pt_lTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lTrack,reweight);
    	//cout << "SE test 6.8 " <<  endl;
                        }

                        if(mStrangenessCut->passPhiEtaWest(lTrack)) // pos eta
                        { // Below is East Only

                            TVector2 Q2Vector = Q2East[j];
                            TVector2 Q3Vector = Q3East[j];
                            // subtract auto-correlation from pos eta(west) event plane
                            if(flagA == 0 && mStrangenessCut->passTrackEP(lTrackA,dcaA) && mStrangenessCut->passTrackEtaEast(lTrackA,j)) // trackA
                            {
                                Float_t  w = mStrangenessCorr->getWeight(lTrackA);

                                TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lTrackA);
                                TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
                                Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

                                TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lTrackA);
                                TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
                                Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
                            }
                            if(flagB == 0 && mStrangenessCut->passTrackEP(lTrackB,dcaB) && mStrangenessCut->passTrackEtaEast(lTrackB,j)) // trackB
                            {
                                Float_t  w = mStrangenessCorr->getWeight(lTrackB);

                                TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lTrackB);
                                TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
                                Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

                                TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lTrackB);
                                TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
                                Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
                            }
                            Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
                            Float_t Psi2_east = mStrangenessCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign,j);
			    Psi2_East_ltrack = Psi2_east ; 
                            Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
                            Float_t Psi3_east = mStrangenessCorr->calShiftAngle3East_EP(Q3Vector,runIndex,cent9,vz_sign,j);
                            Float_t phi_Psi2 = phi_lTrack - Psi2_east;
                            Float_t phi_Psi3 = phi_lTrack - Psi3_east;

                            mStrangenessHistoManger->Fill(pt_lTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lTrack,reweight);
                        }
                    }// Below is East Only
                }// apply additional PID cut to increase significance
            }// loop over all tracks of the actual event
        }// passTrackNumCut
        if(Psi2_East_ltrack > -999. && Psi2_West_ltrack > -999.){
        	//cout << "SE  psi2 east: " << Psi2_East_ltrack << " psi2 west: " << Psi2_West_ltrack << endl;
                mStrangenessHistoManger->Fill_EP_QA_West(Psi2_West_ltrack);
                mStrangenessHistoManger->Fill_EP_QA_East(Psi2_East_ltrack);
        	mStrangenessHistoManger->Fill_EPs_QA(Psi2_East_ltrack, Psi2_West_ltrack);
        }
    }// get Track Information

    cout << "." << flush;
    cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
    cout << endl;
}

// loop phi meson Mixed Event
void StStrangenessAna::MakePhiME()
{
    Long64_t start_event_use = mStartEvent_ME;
    Long64_t stop_event_use  = mStopEvent_ME;
    cout << "test ME 0 " << endl;

    mInPut_ME->SetBranchAddress( XUV0_EVENT_BRANCH, &mXuPhiMeson_event );
    mInPut_ME->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry
    cout << "test ME 1 " << endl;

    // Initialise Event Head
    TVector3 PrimaryVertex(0.0,0.0,0.0);
    Int_t          RunId = 0;
    //Int_t          RefMult = 0;
    Int_t          Centrality = 0;
    Float_t        Reweight = 0;
    //Int_t          N_prim = 0;
    //Int_t          N_non_prim = 0;
    //Int_t          N_Tof_match = 0;
    //Float_t        ZDCx = 0.0; 
    //Float_t        BBCx = 0.0; 
    //Float_t        VzVpd = 0.0;
    Int_t          NumTrackUsed = 0;
    // ---------------------------------------QVector---------------------------------------------
    TVector2       Q1Full[4];
    TVector2       Q2East[4];
    TVector2       Q2West[4];
    TVector2       Q3East[4];
    TVector2       Q3West[4];
    // -----------------------------------Number of Tracks----------------------------------------
    Int_t          NumTrackEast[4];
    Int_t          NumTrackWest[4];
    for(Int_t j = 0; j < 4; j++)
    {
        Q1Full[j].Set(0.0,0.0);
        Q2East[j].Set(0.0,0.0);
        Q2West[j].Set(0.0,0.0);
        Q3East[j].Set(0.0,0.0);
        Q3West[j].Set(0.0,0.0);
        NumTrackEast[j] = 0;
        NumTrackWest[j] = 0;
    }
    cout << "test ME 2 " << endl;

    for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
    {
        if (!mInPut_ME->GetEntry( counter )) // take the event -> information is stored in event
            break;  // end of data chunk

        // get Event Header
        PrimaryVertex     = mXuPhiMeson_event->getPrimaryVertex();
        RunId             = mXuPhiMeson_event->getRunId();
        //RefMult           = mXuPhiMeson_event->getRefMult();
        Centrality        = mXuPhiMeson_event->getCentrality();
        Reweight          = mXuPhiMeson_event->getReweight();
        //N_prim            = mXuPhiMeson_event->getN_prim();
        //N_non_prim        = mXuPhiMeson_event->getN_non_prim();
        //N_Tof_match       = mXuPhiMeson_event->getN_Tof_match();
        //ZDCx              = mXuPhiMeson_event->getZDCx(); 
        //BBCx              = mXuPhiMeson_event->getBBCx(); 
        //VzVpd             = mXuPhiMeson_event->getVzVpd();
        NumTrackUsed      = mXuPhiMeson_event->getNumTracks();

        for(Int_t j = 0; j < 4; j++)
        {
            Q1Full[j]       = mXuPhiMeson_event->getQ1Full(j);
            Q2East[j]       = mXuPhiMeson_event->getQ2East(j);
            Q2West[j]       = mXuPhiMeson_event->getQ2West(j);
            Q3East[j]       = mXuPhiMeson_event->getQ3East(j);
            Q3West[j]       = mXuPhiMeson_event->getQ3West(j);
            NumTrackEast[j] = mXuPhiMeson_event->getNumTrackEast(j);
            NumTrackWest[j] = mXuPhiMeson_event->getNumTrackWest(j);
        }
        Float_t Psi1_Full = -999.9;
        Float_t Psi2_East_ltrack = -999.9, Psi2_West_ltrack = -999.9;

        // Initialise Track 
        Float_t m2A = -999.9;
        Float_t m2B = -999.9;
        Float_t nsA = -999.9;
        Float_t nsB = -999.9;
        Float_t dcaA = -999.9;
        Float_t dcaB = -999.9;
        TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
        TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
        Int_t flagA = -1;
        Int_t flagB = -1;

        //cout << "vz = " << PrimaryVertex.Z() << endl;
        //if(fabs(PrimaryVertex.Z()) > 40.0)  continue;
        //cout << "vz corr = " << PrimaryVertex.Z() << endl;

        // vz sign
        Int_t vz_sign;
        if(PrimaryVertex.Z() > 0.0)
        {
            vz_sign = 0;
        }
        else
        {
            vz_sign = 1;
        }

        //float vx=PrimaryVertex.X();
        //float vy=PrimaryVertex.Y();
        //float vz=PrimaryVertex.Z();

        //if( (vx < 1.e-5 && vx > -1.e-5) &&
        //        (vy < 1.e-5 && vy > -1.e-5) &&
        //        (vz < 1.e-5 && vz > -1.e-5)  )
        //    continue;


        // Centrality
        //mRefMultCorr->init(RunId);
        //mRefMultCorr->initEvent(RefMult,PrimaryVertex.z());
        const Int_t cent9 = Centrality;
        const Double_t reweight = Reweight;
        if(cent9 > 8 || cent9 < 0) continue;

        // runIndex
        //mRunIdEventsDb = StRunIdEventsDb::Instance(Strangeness::mBeamEnergy[mEnergy],Strangeness::mBeamYear[mEnergy]);
        //const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(RunId); // expensive
        // cout << runIndex << endl;
        const int runIndex = GetRunIndex(RunId);

            if (counter != 0  &&  counter % 1000 == 0)
                cout << "." << flush;
        if (counter != 0  &&  counter % 10000 == 0)
        {
            if((stop_event_use-start_event_use) > 0)
            {
                Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
                cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (strangeness_flow) " << flush;
            }
        }
        TVector2 Q1EPD_Full = Q1Full[0];
        //Psi1_Full = TMath::ATan2(Q1EPD_Full.Y(),Q1EPD_Full.X());
        Psi1_Full = mStrangenessCorr->calShiftAngle1EPDfull_EP(Q1EPD_Full,runIndex,cent9,vz_sign,0);
        mStrangenessHistoManger->FillPsiShiftFull_wt(Psi1_Full);

        // get Track Information
        for(Int_t j = 0; j < Strangeness::mEtaGap_total; j++)
        {
            if(mStrangenessCorr->passTrackNumCut(NumTrackEast[j],NumTrackWest[j]))
            {
                for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
                {
                    mXuPhiMeson_track = mXuPhiMeson_event->getTrack(nTracks);
                    m2A = mXuPhiMeson_track->getMass2A();
                    m2B = mXuPhiMeson_track->getMass2B();
                    nsA = mXuPhiMeson_track->getNSigKaonA();
                    nsB = mXuPhiMeson_track->getNSigKaonB();
                    dcaA = mXuPhiMeson_track->getDcaA();
                    dcaB = mXuPhiMeson_track->getDcaB();
                    lTrackA = mXuPhiMeson_track->getTrackA();
                    lTrackB = mXuPhiMeson_track->getTrackB();
                    Float_t pA = lTrackA.P();
                    Float_t pB = lTrackB.P();
                    flagA = mXuPhiMeson_track->getFlagA();
                    flagB = mXuPhiMeson_track->getFlagB();
                    TLorentzVector lTrack = lTrackA + lTrackB;
                    Float_t pt_lTrack = lTrack.Perp();

                    // apply additional PID cut to increase significance 
                    if(
                            ((fabs(pA) <= 0.65 && m2A < -10) || (m2A > 0 && ((fabs(pA) < 1.5 && m2A > 0.16 && m2A < 0.36) || (fabs(pA) >= 1.5 && m2A > 0.125 && m2A < 0.36)) )) &&
                            ((fabs(pB) <= 0.65 && m2B < -10) || (m2B > 0 && ((fabs(pB) < 1.5 && m2B > 0.16 && m2B < 0.36) || (fabs(pB) >= 1.5 && m2B > 0.125 && m2B < 0.36)) )) &&
                            ((pt_lTrack) < 0.8 || ((pt_lTrack) >= 0.8 && ( (m2A > 0.16 && m2A < 0.36) || (m2B > 0.16 && m2B < 0.36)))) &&
                            (
                             ((m2A < -10 && nsA < 2.5 && nsA > -(1.5)) || (m2A > 0.16 && m2A < 0.36)) &&
                             ((m2B < -10 && nsB < 2.5 && nsB > -(1.5)) || (m2B > 0.16 && m2B < 0.36))
                            )
                      )
                    {
                        Float_t phi_lTrack = lTrack.Phi();
                        Float_t InvMass_lTrack = lTrack.M();
                        Float_t y_lTrack = lTrack.Rapidity();
                        mStrangenessHistoManger->FillAcc(pt_lTrack, y_lTrack, InvMass_lTrack);

                        if(mStrangenessCut->passPhiEtaEast(lTrack)) // neg eta(east)
                        { // Below is West Only
                            TVector2 Q2Vector = Q2West[j];
                            TVector2 Q3Vector = Q3West[j];
                            // subtract auto-correlation from pos eta(west) event plane
                            if(flagA == 0 && mStrangenessCut->passTrackEP(lTrackA,dcaA) && mStrangenessCut->passTrackEtaWest(lTrackA,j)) // trackA
                            {
                                Float_t  w = mStrangenessCorr->getWeight(lTrackA);

                                TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lTrackA);
                                TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
                                Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

                                TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lTrackA);
                                TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
                                Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
                            }
                            if(flagB == 0 && mStrangenessCut->passTrackEP(lTrackB,dcaB) && mStrangenessCut->passTrackEtaWest(lTrackB,j)) // trackB
                            {
                                Float_t  w = mStrangenessCorr->getWeight(lTrackB);

                                TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lTrackB);
                                TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
                                Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

                                TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lTrackB);
                                TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
                                Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
                            }
                            Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
                            Float_t Psi2_west = mStrangenessCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign,j);
			    Psi2_West_ltrack = Psi2_west ; 
                            Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
                            Float_t Psi3_west = mStrangenessCorr->calShiftAngle3West_EP(Q3Vector,runIndex,cent9,vz_sign,j);
                            Float_t phi_Psi2 = phi_lTrack - Psi2_west;
                            Float_t phi_Psi3 = phi_lTrack - Psi3_west;

                            mStrangenessHistoManger->Fill(pt_lTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lTrack,reweight);
                        }

                        if(mStrangenessCut->passPhiEtaWest(lTrack)) // pos eta
                        { // Below is East Only
                            TVector2 Q2Vector = Q2East[j];
                            TVector2 Q3Vector = Q3East[j];
                            // subtract auto-correlation from neg eta(east) event plane
                            if(flagA == 0 && mStrangenessCut->passTrackEP(lTrackA,dcaA) && mStrangenessCut->passTrackEtaEast(lTrackA,j)) // trackA
                            {
                                Float_t  w = mStrangenessCorr->getWeight(lTrackA);

                                TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lTrackA);
                                TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
                                Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

                                TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lTrackA);
                                TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
                                Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
                            }
                            if(flagB == 0 && mStrangenessCut->passTrackEP(lTrackB,dcaB) && mStrangenessCut->passTrackEtaEast(lTrackB,j)) // trackB
                            {
                                Float_t  w = mStrangenessCorr->getWeight(lTrackB);

                                TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lTrackB);
                                TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
                                Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

                                TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lTrackB);
                                TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
                                Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
                            }
                            Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
                            Float_t Psi2_east = mStrangenessCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign,j);
			    Psi2_East_ltrack = Psi2_east ; 
                            Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
                            Float_t Psi3_east = mStrangenessCorr->calShiftAngle3East_EP(Q3Vector,runIndex,cent9,vz_sign,j);
                            Float_t phi_Psi2 = phi_lTrack - Psi2_east;
                            Float_t phi_Psi3 = phi_lTrack - Psi3_east;

                            mStrangenessHistoManger->Fill(pt_lTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lTrack,reweight);
                        }
                    }
                }
            }
        }
        if(Psi2_East_ltrack > -999. && Psi2_West_ltrack > -999.){
        	//cout << "ME psi2 east: " << Psi2_East_ltrack << " psi2 west: " << Psi2_West_ltrack << endl;
                mStrangenessHistoManger->Fill_EP_QA_West(Psi2_West_ltrack);
                mStrangenessHistoManger->Fill_EP_QA_East(Psi2_East_ltrack);
        	mStrangenessHistoManger->Fill_EPs_QA(Psi2_East_ltrack, Psi2_West_ltrack);
        }
    }

    cout << "." << flush;
    cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
    cout << endl;
}

// loop Lambda Same Event
// loop Lambda Mixed Event
//-------------------------------------------------------------------
void StStrangenessAna::Finish()
{
    mFile_OutPut->cd();
    mStrangenessHistoManger->Write();
    mFile_OutPut->Close();
}

int StStrangenessAna::GetRunIndex(int runID)
{
    int runIndex=-999;
    for(int i=0; i<2704; i++)
    {
        if(runID==numbers[i])
        {
            runIndex=i;
        }
    }
    if(runIndex == -999) cout << "Run numbers are not found!!!!" << endl;
    return runIndex;
}

