
void AnalysisTaskRunXi1820(double mom=2.7, int nevts=0,  TString pre = ""){
	TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);

	//Output File
	if (pre==""){
		TString Path = "/home/ikp1/puetz/panda/myscripts/simChain/";
		TString outPath = Path;
		TString OutputFile = outPath + "analysis_output.root";

		//Input simulation Files
		TString inPIDFile = Path + "pid_complete.root";
		TString inParFile = Path + "simparams.root";
	}
	else{
		TString Path = pre;
		TString outPath = Path + "_test_";
		TString OutputFile = Path + "_analysis_output.root";

		//Input simulation Files
		TString inPIDFile = Path + "_pid_complete.root";
		TString inParFile = Path + "_simparams.root";
	}

	TString PIDParFile = TString( gSystem->Getenv("VMCWORKDIR")) + "/macro/params/all.par";

	//Initialization
	FairLogger::GetLogger()->SetLogToFile(kFALSE);
	FairRunAna* RunAna = new FairRunAna();
	FairRuntimeDb* rtdb = RunAna->GetRuntimeDb();
	RunAna->SetInputFile(inPIDFile);

	//setup parameter database
	FairParRootFileIo* parIo = new FairParRootFileIo();
	parIo->open(inParFile);
	FairParAsciiFileIo* parIoPID = new FairParAsciiFileIo();
	parIoPID->open(PIDParFile.Data(),"in");

	rtdb->setFirstInput(parIo);
	rtdb->setSecondInput(parIoPID);
	rtdb->setOutput(parIo);


	RunAna->SetOutputFile(OutputFile);

	// *** HERE OUR TASK GOES!
	AnalysisTaskXi1820 *anaTask = new AnalysisTaskXi1820();
	anaTask->SetOutPutDir(outPath);
	anaTask->SetNEvents(nevts);
	anaTask->SetMom(mom);
	RunAna->AddTask(anaTask);

	RunAna->Init();
	RunAna->Run(0.,1.);

	exit(0);
}
