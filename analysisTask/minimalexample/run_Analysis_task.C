
void run_Analysis_task(int nevts=0){

	TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);

	//Output File
	TString Path ="/private/puetz/mysimulations/test/boxgenerator/2pi_vtx_000/";
	TString outPath =Path;

	TString OutputFile = outPath + "analysis_output_test.root";

	//Input simulation Files
	TString inPIDFile = Path + "pid_complete.root";
	TString inParFile = Path + "simparams.root";
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
	AnalysisTask *anaTask = new AnalysisTask();
	RunAna->AddTask(anaTask);

	RunAna->Init();
	RunAna->Run(0.,1.);

}
