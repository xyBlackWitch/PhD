//bool checkfile(TString fn)
//{
//	bool fileok=true;
//	TFile fff(fn); 
//	if (fff.IsZombie()) fileok=false;
//	TTree *t=(TTree*)fff.Get("cbmsim");
//	if (t==0x0) fileok=false;
//	
//	if (!fileok) cout <<"Skipping broken file '"<<fn<<"'"<<endl;
//	return fileok;
//}

void prod_ana_Xi1820(TString outpre="M9999", int nevts=0, double mom=4.6)
{
	TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);


 	TString OutFile1   = TString::Format("%s_ana.root",outpre.Data()); 
 	TString OutFile2   = TString::Format("%s_",outpre.Data()); 
	TString inParFile = TString::Format("%s_par.root",outpre.Data());
	
  	FairRunAna *fRun= new FairRunAna();
  
   	//bool firstfile=true;

  	// *** Add pid files
  	//for (int i=id;i<=to;++i)
  	//{
	
	TString fname = TString::Format("%s_pid.root",outpre.Data());

	//	if ( checkfile(fname) )
	//	{
	//		if (firstfile) 
				fRun->SetInputFile(fname);
	//		else 
	//			fRun->AddFile(fname);
	//		
	//		firstfile=false;
	//	}
  	//}
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

	RunAna->AddFriend(RecoFile);
	RunAna->SetOutputFile(OutputFile);


	// *** HERE OUR TASK GOES!
	AnalysisTaskXi1820 *anaTask = new AnalysisTaskXi1820();
	anaTask->SetOutPutDir(outPath);
	anaTask->SetNEvents(nevts);
	anaTask->SetMom(mom);
	RunAna->AddTask(anaTask);

	RunAna->Init();
	RunAna->Run(0.,1.);

	
}

















