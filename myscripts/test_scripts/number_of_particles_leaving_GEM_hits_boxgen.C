
void GemHits(RhoCandidate *c){

	int tag = 0;

	PndPidCandidate * pidCand = (PndPidCandidate*)c->GetRecoCandidate();


	if(pidCand){
		int gemHits = pidCand->GetGemHits();

		if(gemHits>0) tag=1;
		else tag=0;

	}

	return tag;
}

void number_of_particles_leaving_GEM_hits_boxgen(TString pre="", int nevts=0, double mom=4.1){
	  TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);

	  TStopwatch timer;

	  if (pre==""){
		  //Output File
		  TString OutputFile = "test_analysis_output.root";
		  TString outPath = "";
		  //Input simulation Files
		 TString inPIDFile = "pid_complete.root";
		 TString inParFile = "simparams.root";
	  }
	  else {
		  //Output File
		  TString outPath = pre + "_";
		  TString OutputFile = pre + "_test_analysis_output.root";

		  //Input simulation Files
		  TString inPIDFile = pre + "_pid_complete.root";
		  TString inParFile = pre + "_simparams.root";
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
	  RunAna->Init();


	  /*************************************************************************
	   *  Create new ntuple and fill them with information
	   ************************************************************************/

	  //*** create tuples
	  RhoTuple * ntpPiMinus = new RhoTuple("ntpPiMinus", "PiMinus info");
	  RhoTuple * ntpPiPlus = new RhoTuple("ntpPiPlus", "PiPlus info");
	  RhoTuple * ntpKaonMinus = new RhoTuple("ntpKaonMinus", "KaonMinus info");
	  RhoTuple * ntpKaonPlus = new RhoTuple("ntpKaonPlus", "KaonPlus info");
	  RhoTuple * ntpProton = new RhoTuple("ntpProton", "Proton info");
	  RhoTuple * ntpAntiProton = new RhoTuple("ntpAntiProton", "Antiproton info");


	  //Create output file
	  TFile *out = TFile::Open(outPath+"test_output_ana.root","RECREATE");

	  // data reader Object
	  PndAnalysis* theAnalysis = new PndAnalysis();
	  if (nevts==0) nevts = theAnalysis->GetEntries();


	  //RhoCandLists for analysis
	  RhoCandList piplus, piminus, proton, antiproton, kaonminus, kaonplus;

	  RhoCandidate * dummyCand = new RhoCandidate(); //dummy candidate for empty candidate usage


	  double p_m0 = TDatabasePDG::Instance()->GetParticle("proton")->Mass();


	  TLorentzVector ini (0,0, mom, sqrt(p_m0*p_m0+ mom*mom)+p_m0);
	  TVector3 beamBoost = ini.BoostVector();

	  PndRhoTupleQA qa(theAnalysis, mom);

	  int evt=-1;
	  while (theAnalysis->GetEvent() && ++evt<nevts){

	    if ((evt%100)==0) cout << "evt "<< evt <<endl;


	    TString PidSelection = "PidAlgoIdealCharged";//"PidAlgoMvd;PidAlgoStt;PidAlgoDrc"; to change from ideal PID to realistic PID uncomment this!


	    //***Selection with no PID info
	    theAnalysis->FillList(piminus, "PionBestMinus", PidSelection);
	    theAnalysis->FillList(piplus, "PionBestPlus", PidSelection);
	    theAnalysis->FillList(kaonminus, "KaonBestMinus", PidSelection);
	    theAnalysis->FillList(kaonplus, "KaonBestPlus", PidSelection);
	    theAnalysis->FillList(proton, "ProtonBestPlus", PidSelection);
	    theAnalysis->FillList(antiproton, "ProtonBestMinus", PidSelection);


	    //Get piminus information
	    ntpPiMinus->Column("ev",     (Float_t) evt);

	    for (int j=0; j<piminus.GetLength(); ++j){


			//information about the mother and MCTruth Candidate

			TLorentzVector l;
			float costheta = -999.;

			RhoCandidate * truth = piminus[j]->GetMcTruth();
			RhoCandidate * mother;
			if (truth)  mother = truth->TheMother();

			int moth = (mother==0x0) ? 88888 : mother->PdgCode();
			ntpPiMinus->Column("Mother", (Int_t) moth);

			bool truthmatch = theAnalysis->McTruthMatch(piminus[j]);
			ntpPiMinus->Column("MCTruthMatch", (bool) truthmatch);

			int gemhit = GemHits(piminus[j]);

			int count = 0;
			if (moth==88888 && gemhit==1 && truthmatch==1) count=1;

			ntpPiMinus->Column("GemHit", (int) count, 0);


	    }
	    ntpPiMinus->DumpData();

		//Get PiPlus information

	    ntpPiPlus->Column("ev", (int) evt);

	    for (int j=0; j<piplus.GetLength(); ++j){


			//information about the mother and MCTruth Candidate

			TLorentzVector l;
			float costheta = -999.;

			RhoCandidate * truth = piplus[j]->GetMcTruth();
			RhoCandidate * mother;
			if (truth)  mother = truth->TheMother();

			int moth = (mother==0x0) ? 88888 : mother->PdgCode();
			ntpPiPlus->Column("Mother", (Int_t) moth);

			bool truthmatch = theAnalysis->McTruthMatch(piplus[j]);
			ntpPiPlus->Column("MCTruthMatch", (bool) truthmatch);

			int gemhit = GemHits(piplus[j]);

			int count = 0;
			if (moth==88888 && gemhit==1 && truthmatch==1) count=1;

			ntpPiPlus->Column("GemHit", (int) count, 0);


	    }

	    ntpPiPlus->DumpData();

	    ntpKaonMinus->Column("ev", (int) evt);

	    for (int j=0; j<kaonminus.GetLength(); ++j){


			//information about the mother and MCTruth Candidate

			TLorentzVector l;
			float costheta = -999.;

			RhoCandidate * truth = kaonminus[j]->GetMcTruth();
			RhoCandidate * mother;
			if (truth)  mother = truth->TheMother();

			int moth = (mother==0x0) ? 88888 : mother->PdgCode();
			ntpKaonMinus->Column("Mother", (Int_t) moth);

			bool truthmatch = theAnalysis->McTruthMatch(kaonminus[j]);
			ntpKaonMinus->Column("MCTruthMatch", (bool) truthmatch);

			int gemhit = GemHits(kaonminus[j]);

			int count = 0;
			if (moth==88888 && gemhit==1 && truthmatch==1) count=1;

			ntpKaonMinus->Column("GemHit", (int) count, 0);


	    }

	    ntpKaonMinus->DumpData();


	    ntpKaonPlus->Column("ev", (int) evt);

	    for (int j=0; j<kaonplus.GetLength(); ++j){


			//information about the mother and MCTruth Candidate

			TLorentzVector l;
			float costheta = -999.;

			RhoCandidate * truth = kaonplus[j]->GetMcTruth();
			RhoCandidate * mother;
			if (truth)  mother = truth->TheMother();

			int moth = (mother==0x0) ? 88888 : mother->PdgCode();
			ntpKaonPlus->Column("Mother", (Int_t) moth);

			bool truthmatch = theAnalysis->McTruthMatch(kaonplus[j]);
			ntpKaonPlus->Column("MCTruthMatch", (bool) truthmatch);

			int gemhit = GemHits(kaonplus[j]);

			int count = 0;
			if (moth==88888 && gemhit==1 && truthmatch==1) count=1;

			ntpKaonPlus->Column("GemHit", (int) count, 0);


	    }

	    ntpKaonPlus->DumpData();

//		Get Proton information
	    ntpProton->Column("ev", (int) evt);

	    for (int j=0; j<proton.GetLength(); ++j){


			//information about the mother and MCTruth Candidate

			TLorentzVector l;
			float costheta = -999.;

			RhoCandidate * truth = proton[j]->GetMcTruth();
			RhoCandidate * mother;
			if (truth)  mother = truth->TheMother();

			int moth = (mother==0x0) ? 88888 : mother->PdgCode();
			ntpProton->Column("Mother", (Int_t) moth);

			bool truthmatch = theAnalysis->McTruthMatch(proton[j]);
			ntpProton->Column("MCTruthMatch", (bool) truthmatch);

			int gemhit = GemHits(proton[j]);

			int count = 0;
			if (moth==88888 && gemhit==1 && truthmatch==1) count=1;

			ntpProton->Column("GemHit", (int) count, 0);


	    }

	    ntpProton->DumpData();

//			Get Antiproton
	    ntpAntiProton->Column("ev", (int) evt);

	    for (int j=0; j<antiproton.GetLength(); ++j){


			//information about the mother and MCTruth Candidate

			TLorentzVector l;
			float costheta = -999.;

			RhoCandidate * truth = antiproton[j]->GetMcTruth();
			RhoCandidate * mother;
			if (truth)  mother = truth->TheMother();

			int moth = (mother==0x0) ? 88888 : mother->PdgCode();
			ntpAntiProton->Column("Mother", (Int_t) moth);

			bool truthmatch = theAnalysis->McTruthMatch(antiproton[j]);
			ntpAntiProton->Column("MCTruthMatch", (bool) truthmatch);

			int gemhit = GemHits(antiproton[j]);

			int count = 0;
			if (moth==88888 && gemhit==1 && truthmatch==1) count=1;

			ntpAntiProton->Column("GemHit", (int) count, 0);


	    }

	    ntpAntiProton->DumpData();

	  }


	  //Write output
	  out->cd();

	  ntpPiMinus ->GetInternalTree()->Write();
	  ntpPiPlus->GetInternalTree()->Write();
	  ntpKaonMinus ->GetInternalTree()->Write();
	  ntpKaonPlus->GetInternalTree()->Write();
	  ntpProton->GetInternalTree()->Write();
	  ntpAntiProton->GetInternalTree()->Write();


	  out->Save();

	  timer.Stop();
	  Double_t rtime = timer.RealTime();
	  Double_t ctime = timer.CpuTime();

	  cout<<"Macro finisched successfully."<<endl;
	  cout<<"Realtime: "<<rtime<<" s, CPU time: "<<ctime<<" s"<<endl;
	  cout<<endl;


	  exit(0);

	}

