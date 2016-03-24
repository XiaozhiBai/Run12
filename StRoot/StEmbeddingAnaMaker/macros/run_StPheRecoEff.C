void run_StPheRecoEff(const char* infiles="@test_file.list", const char* outfile) 
{
    //Check STAR Library. Please set SL_version to the original star library used 
    // in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
    //string SL_version = "SL12d";
    string SL_version = "SL12d_embed"; // DEBUG
    string env_SL = getenv("STAR");
    if (env_SL.find(SL_version) == string::npos) {
        cout << "Environment Star Library does not match the requested library in run_st_etree.C. Exiting..." << endl;
        return;
    }
    // Load shared libraries
    Load();

    // Create chain
    StChain* chain = new StChain;

    // FILES
    StFile* files = new StFile();
    files->AddFile(infiles);

    // I/O maker
    // StIOMaker* ioMaker = new StIOMaker;
    // ioMaker->SetFile(file);
    // ioMaker->SetIOMode("r");
    StIOMaker* ioMaker = new StIOMaker("ioMaker","r",files);
    // ioMaker->SetFile(file);
    // ioMaker->SetIOMode("r");
    ioMaker->SetBranch("*", 0, "0");
    //ioMaker->SetBranch("McEventBranch",0,"r"); 
    ioMaker->SetBranch("geantBranch", 0, "r");
    ioMaker->SetBranch("eventBranch", 0, "r");


    StMcEventMaker *mcEventMaker = new StMcEventMaker();
    mcEventMaker->doPrintEventInfo = false;
    mcEventMaker->doPrintMemoryInfo = false;

    StAssociationMaker* assoc = new StAssociationMaker;
    assoc->useInTracker();

    //.. see example in CVS: StRoot/macros/mudst/exampleEmc.C
    // Need St_db_Maker for Emc calibration
    St_db_Maker* dbMk = new St_db_Maker("StarDb", "MySQL:StarDb", "MySQL:StarDb", "$STAR/StarDb");
    //dbMk->SetMaxEntryTime(20100301,0);
    //dbMk->SetDateTime(20080101,000001);

    StEmcSimulatorMaker* emcSim = new StEmcSimulatorMaker(); //use this instead to "redo" converstion from geant->adc
    // Makers for clusterfinding
    StPreEclMaker *pre_ecl = new StPreEclMaker();
    pre_ecl->setPrint(kFALSE);
    StEpcMaker *epc = new StEpcMaker();
    epc->setPrint(kFALSE);

    int col_energy = 200; //.. option: 200, 62, 39
    // std::vector<unsigned int> triggers;
    // triggers.push_back(370542); //
    // triggers.push_back(370522); //

    // Monte Carlo event maker
    StPheRecoEff* analysis = new StPheRecoEff(outfile,1,col_energy);//,triggers);

    // Initialize chain
    chain->Init();
    chain->EventLoop(1000000);
    chain->Finish();
    
    //delete chain;
}

void Load() {
    gROOT->Macro("loadMuDst.C");
    gROOT->Macro("LoadLogger.C");
    gSystem->Load("StMcEvent");
    gSystem->Load("StMcEventMaker");
    gSystem->Load("StEEmcUtil");
    gSystem->Load("StAssociationMaker");
    gSystem->Load("StDbLib.so");
    gSystem->Load("StDbBroker.so");
    gSystem->Load("libglobal_Tables.so");
    gSystem->Load("St_db_Maker.so");
    gSystem->Load("StDetectorDbMaker");
    gSystem->Load("StTpcDb");
    gSystem->Load("StDbUtilities");
    gSystem->Load("StMcEvent");
    gSystem->Load("StMcEventMaker");
    gSystem->Load("StAssociationMaker");
    gSystem->Load("StDaqLib");
    gSystem->Load("StEmcRawMaker");
    gSystem->Load("StEmcADCtoEMaker");
    gSystem->Load("StPreEclMaker");
    gSystem->Load("StEpcMaker");
    gSystem->Load("StEmcSimulatorMaker");
    gSystem->Load("StEmcUtil");
    gSystem->Load("StEEmcUtil");
    gSystem->Load("StEEmcDbMaker");
    gSystem->Load("StEmcTriggerMaker");
    //  gSystem->Load("StMcAnalysisMaker");
    gSystem->Load("StEmbeddingAnaMaker");
}
