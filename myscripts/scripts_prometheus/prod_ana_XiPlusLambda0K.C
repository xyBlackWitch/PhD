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

void prod_ana_XiPlusLambda0K(TString outpre="M9999", int nevts=0, double mom=4.6)
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
  	
	// *** PID table with selection thresholds; can be modified by the user
	TString pidParFile = TString(gSystem->Getenv("VMCWORKDIR"))+"/macro/params/all.par";
	
	// *** initialization
	FairLogger::GetLogger()->SetLogToFile(kFALSE);
	FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
	
	// *** setup parameter database 	
	
	FairParRootFileIo* parIO = new FairParRootFileIo();
	parIO->open(inParFile);
	FairParAsciiFileIo* parIOPid = new FairParAsciiFileIo();
	parIOPid->open(pidParFile.Data(),"in");
	
	rtdb->setFirstInput(parIO);
	rtdb->setSecondInput(parIOPid);
	rtdb->setOutput(parIO);  
	
	
	fRun->SetOutputFile(OutFile1);
	
	//---------------------Create and Set the Field(s)---------- 
//  	PndMultiField *fField= new PndMultiField("FULL");
//  	fRun->SetField(fField);
	
	//RhoCalculationTools::ForceConstantBz(20.0);

	// ***
	// *** HERE YOUR ANALYSIS CODE GOES!
	// ***
	AnalysisTaskXiPlusLambda0K *anaTask = new AnalysisTaskXiPlusLambda0K();
	anaTask->SetOutPutDir(OutFile2);
	anaTask->SetNEvents(nevts);
	anaTask->SetMom(mom);
	fRun->AddTask(anaTask);

	fRun->Init();
	fRun->Run(0.,1.);

	
}

















