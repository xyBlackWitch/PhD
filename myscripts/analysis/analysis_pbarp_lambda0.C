class PndAnalysis;
class PndAnaPidSelector;
class RhoCandList;
class RhoTuple;

//#include <map>

/******************************************
 * Methods needed in the analysis task
 *******************************************/

void numberOfHitsInSubdetector(TString pre, RhoCandidate *c, RhoTuple *n){

	/* This method saves the number of Hits in the MVD, STT and GEM detector
	 * into the RhoTuple.
	 */

	PndPidCandidate *pidCand = (PndPidCandidate*)c->GetRecoCandidate();

	if(pidCand){
		n->Column(pre + "MvdHits", (Int_t) pidCand->GetMvdHits(), 0);
		n->Column(pre + "SttHits", (Int_t) pidCand->GetSttHits(), 0);
		n->Column(pre + "GemHits", (Int_t) pidCand->GetGemHits(), 0);

	}
	else{
		n->Column(pre + "MvdHits", (Int_t) -999, 0);
		n->Column(pre + "SttHits", (Int_t) -999, 0);
		n->Column(pre + "GemHits", (Int_t) -999, 0);
	}
}

void tagNHits(TString pre, RhoCandidate *c, RhoTuple *n){

	/**@brief Tag the particle with different integers
	 * @details Tag the particle with different integers:
	 * 0: if there is less than 4 hits in any inner tracking detector
	 * 1: sttHits>3 or mvdHits>3 or gemHit>3
	 */

	int tag=0;

	PndPidCandidate * pidCand = (PndPidCandidate*)c->GetRecoCandidate();

	if(pidCand){
		int mvdHits = pidCand->GetMvdHits();
		int sttHits = pidCand->GetSttHits();
		int gemHits = pidCand->GetGemHits();

		if(mvdHits>3 || sttHits>3 || gemHits>3) tag=1;


	}

	n->Column(pre + "HitTag", (Int_t) tag, 0);
}

int tagHits(RhoCandidate *c){

	/**@brief Tag the particle with different integers
	 * @details Tag the particle with different integers:
	 * 0: if there is no hit in the detector
	 * 1: sttHits>3 or mvdHits>3 or gemHit>3
	 */
	int tag = 0;

	PndPidCandidate * pidCand = (PndPidCandidate*)c->GetRecoCandidate();

	if(pidCand){
		int mvdHits = pidCand->GetMvdHits();
		int sttHits = pidCand->GetSttHits();
		int gemHits = pidCand->GetGemHits();

		if(mvdHits>3 || sttHits>3 || gemHits>3) tag=1;
	//		cout << "mvd: " << mvdHits << " stt: " << sttHits << " gem: " << gemHits << endl;

	}

	return tag;
}

void qaVtxDiff(TString pre, RhoCandidate * c, RhoTuple * n){

	  if(n==0) return;
	  if(c==0) return;

	  RhoCandidate * mct = c->GetMcTruth();

	  if(mct){
		  TVector3 v = c->DecayVtx();
		  TVector3 mcv = mct->Daughter(0)->Pos();
		  TVector3 vdiff = v-mcv;
		  TMatrixD cov7 = c->Cov7();

		  n->Column(pre + "diffvx", (Float_t) vdiff.X(), 0.0f );
		  n->Column(pre + "diffvy", (Float_t) vdiff.Y(), 0.0f );
		  n->Column(pre + "diffvz", (Float_t) vdiff.Z(), 0.0f );

		  n->Column(pre + "pullvx", (Float_t) (vdiff.X()/TMath::Sqrt(cov7(0,0))), 0.0f);
		  n->Column(pre + "pullvy", (Float_t) (vdiff.Y()/TMath::Sqrt(cov7(1,1))), 0.0f);
		  n->Column(pre + "pullvz", (Float_t) (vdiff.Z()/TMath::Sqrt(cov7(2,2))), 0.0f);

	  }
	  else{
		  n->Column(pre + "diffvx", (Float_t) -999.0, 0.0f );
		  n->Column(pre + "diffvy", (Float_t) -999.0, 0.0f );
		  n->Column(pre + "diffvz", (Float_t) -999.0, 0.0f );

		  n->Column(pre + "pullvx", (Float_t) -999.0, 0.0f);
		  n->Column(pre + "pullvy", (Float_t) -999.0, 0.0f);
		  n->Column(pre + "pullvz", (Float_t) -999.0, 0.0f);

	  }
}

void qaMomRes(TString pre, RhoCandidate * c, RhoTuple * n){

	  if(n==0 || c==0) return;

	  RhoCandidate * mct = c->GetMcTruth();
	  float momres = -999.0;

	  if(mct){
		  float p = c->P();
		  float mcp = mct->P();

		  momres = (p-mcp)/mcp;
	  }

	  n->Column(pre + "mom_res", (Float_t) momres, 0.0f);

}

//std::map<int,int> VertexQaIndex(RhoCandList* candList, float probLimit=0.01){
//	  /** @brief  give back the order of the best chi2
//	   * @details give back the order of the best chi2!  1 means best, 2: second best (same with negative valuesfor bad chi2 )
//	   */
//
//	  std::map<double, int> chi2_good, chi2_bad;
//
//	  for (int j=0; j<candList->GetLength(); ++j){
//
//		  PndKinVtxFitter vtxfitter(candList->Get(j));
//		  vtxfitter.Fit();
//
//		  bool failedchi2 = TMath::IsNaN(vtxfitter.GetChi2());
//		  bool failedprob = TMath::IsNaN(vtxfitter.GetProb());
//
//		  if(!failedchi2 && !failedprob){
//
//			  if (vtxfitter.GetProb() > probLimit){ //Prob > 0.01
//				  chi2_good[vtxfitter.GetChi2()]=j;
//			  }
//			  else{ //Prob <= 0.01
//				  chi2_bad[vtxfitter.GetChi2()]=j;
//			  }
//
//		  }
//	  }
//
//	  std::map<double, int>::iterator is_good, is_bad;
//	  std::map<int, int> indexBestFit;
//	  int running = 0;
//
//	  for (is_good = chi2_good.begin(); is_good != chi2_good.end(); is_good++, running++){
//		   indexBestFit[is_good->second] = running + 1;
//	  }
//
//	  running =0;
//
//	  for (is_bad = chi2_bad.begin(); is_bad != chi2_bad.end(); is_bad++, running++){
//		  indexBestFit[is_bad->second] = - (running + 1);
//	  }
//
//
//	  return indexBestFit;
//}

//std::map<int,int> MassFitQaIndex(RhoCandList* candList, float m0, float probLimit=0.01){
//	  /** @brief  give back the order of the best chi2 for MassFit
//	   * @details give back the order of the best chi2 for the MassFit!  1 means best, 2: second best (analoge for bad chi2 with negative values)
//	   */
//
//	  if(m0==0) std::cout << "Mass is missing for mass fit" << std::endl;
//
//	  std::map<double, int> chi2_good, chi2_bad;
//
//	  for (int i=0; i<candList->GetLength(); i++){
//
//		  PndKinFitter massfitter(candList->Get(i));
//		  massfitter.AddMassConstraint(m0);
//		  massfitter.Fit();
//
//		  bool failedchi2 = TMath::IsNaN(massfitter.GetChi2());
//		  bool failedprob = TMath::IsNaN(massfitter.GetProb());
//
//		  if(!failedchi2 && !failedprob){
//
//			  if (massfitter.GetProb() > probLimit){
//				  chi2_good[massfitter.GetChi2()]=i;
//			  }
//			  else{
//				  chi2_bad[massfitter.GetChi2()]=i;
//			  }
//		  }
//	  }
//
//	  std::map<doubl//e,int>::iterator is_good, is_bad;
//	  std::map<int,int> bestMassFit;
//
//	  int run =0;
//
//	  for (is_good = chi2_good.begin(); is_good != chi2_good.end(); is_good++, run++){
//		  bestMassFit[is_good->second] = run + 1;
//	  }
//
//	  run = 0;
//
//	  for (is_bad = chi2_bad.begin(); is_bad != chi2_bad.end(); is_bad++, run++){
//		  bestMassFit[is_bad->second] = - (run + 1);
//	  }
//
//
//	  return bestMassFit;
//}




/**********************************************
 * The analysis part starts
 **********************************************/

void analysis_pbarp_lambda0(double mom=1.712, int nevts=0, TString pre=""){
  
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
  RhoTuple * ntpMC = new RhoTuple("ntpMC", "MCTruth info");
  RhoTuple * ntpPiMinus = new RhoTuple("ntpPiMinus", "PiMinus info");
  RhoTuple * ntpPiPlus = new RhoTuple("ntpPiPlus", "PiPlus info");
  RhoTuple * ntpProton = new RhoTuple("ntpProton", "Proton info");
  RhoTuple * ntpAntiProton = new RhoTuple("ntpAntiProton", "Antiproton info");
  RhoTuple * ntpLambda0 = new RhoTuple("ntpLambda0", "Lambda0 info");
  RhoTuple * ntpAntiLambda0 = new RhoTuple("ntpAntiLambda0", "AntiLambda0 info");
  RhoTuple * ntpSys = new RhoTuple("ntpSys", "pbar p system info");

  //Create output file 
  TFile *out = TFile::Open(outPath+"test_output_ana.root","RECREATE");

  // data reader Object
  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts = theAnalysis->GetEntries();
  

  //RhoCandLists for analysis
  RhoCandList piplus, piminus, lambda0, antiLambda0, proton, antiProton, sys, Lambda0Fit, AntiLambda0Fit;
  RhoCandList mclist, all;

  RhoCandidate * dummyCand = new RhoCandidate(); //dummy candidate for empty candidate usage

  //***Mass selector for Lambda0
  double m0_lambda0= TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
  cout<<"Mass of Lambda: "<<m0_lambda0<<endl;
  RhoMassParticleSelector * lambdaMassSelector = new RhoMassParticleSelector("lambda0", m0_lambda0, 0.3);
 
  //*** get information about pbar p system
  double m0_pbarpsystem = TDatabasePDG::Instance()->GetParticle("pbarpSystem")->Mass();
  cout << "Mass pbar p system: "<< m0_pbarpsystem << endl;
  double p_m0 = TDatabasePDG::Instance()->GetParticle("proton")->Mass();

  
  TLorentzVector ini (0,0, mom, sqrt(p_m0*p_m0+ mom*mom)+p_m0);
  TVector3 beamBoost = ini.BoostVector();
  
  PndRhoTupleQA qa(theAnalysis, mom);
 
  int evt=-1;
  while (theAnalysis->GetEvent() && ++evt<nevts){

    if ((evt%100)==0) cout << "evt "<< evt <<endl;
    cout << "Event number: " << evt << endl;

    //***get MC list and store info
    theAnalysis->FillList(mclist, "McTruth");
    qa.qaMcList("", mclist, ntpMC);
    ntpMC->DumpData();
   
		
	//if you want to print the hole MCTree uncomment the following
    /*
    for (int j=0;j<mclist.GetLength();++j)
    {
      RhoCandidate *mcmother = mclist[j]->TheMother();        // mother of mc particle         
      int muid = (mcmother==0x0) ? -1 : mcmother->GetTrackNumber(); // track ID of mother, if existing 
        
      cout << "Track "<< mclist[j]->GetTrackNumber()<<" (PDG:"<<mclist[j]->PdgCode() <<") has mother "<<muid;
      if (mclist[j]->NDaughters()>0) cout <<" and daughter(s) ";
	 for (k=0;k<mclist[j]->NDaughters();++k) cout <<mclist[j]->Daughter(k)->GetTrackNumber()<<"  ";
	cout<<endl;        
    }*/


    //***Setup event shape object

    TString PidSelection = "PidAlgoIdealCharged";//"PidAlgoMvd;PidAlgoStt;PidAlgoDrc"; to change from ideal PID to realistic PID uncomment this!

    theAnalysis->FillList(all, "All", PidSelection);
    PndEventShape evsh(all, ini, 0.05, 0.1);
    
    //***Selection with no PID info
    theAnalysis->FillList(piminus, "PionBestMinus", PidSelection);
    theAnalysis->FillList(piplus, "PionBestPlus", PidSelection);
    theAnalysis->FillList(proton, "ProtonBestPlus", PidSelection);
    theAnalysis->FillList(antiProton, "ProtonBestMinus", PidSelection);


    //Get piminus information
    for (int j=0; j<piminus.GetLength(); ++j){

		//general info about event
		ntpPiMinus->Column("ev",     (Float_t) evt);
		ntpPiMinus->Column("cand",    (Float_t) j);
		ntpPiMinus->Column("ncand",   (Float_t) piminus.GetLength());


		//info about 4-vector
		qa.qaP4("PiMinus_",piminus[j]->P4(), ntpPiMinus);
		ntpPiMinus->Column("PiMinus_CosTheta", (Float_t) piminus[j]->GetMomentum().CosTheta());
		qa.qaCand("PiMinus_", piminus[j], ntpPiMinus);

		//get number of Hits in subdetector (only needed if ideal PR is used!!!)
		numberOfHitsInSubdetector("PiMinus_", piminus[j], ntpPiMinus);
		tagNHits("PiMinus_", piminus[j], ntpPiMinus);

		//information about the mother and MCTruth Candidate

		TLorentzVector l;
		float costheta = -999.;

		RhoCandidate * truth = piminus[j]->GetMcTruth();
		RhoCandidate * mother;
		if (truth)  mother = truth->TheMother();

		int moth = (mother==0x0) ? 88888 : mother->PdgCode();
		ntpPiMinus->Column("Mother", (Int_t) moth);

		ntpPiMinus->Column("MCTruthMatch", (bool) theAnalysis->McTruthMatch(piminus[j]));

		if (truth){
			l = truth->P4();
			costheta = truth->GetMomentum().CosTheta();
			qa.qaCand("PiMiuns_MC_", truth, ntpPiMinus);
		}
		else{
			qa.qaCand("PiMiuns_MC_", dummyCand, ntpPiMinus);
		}

		ntpPiMinus->Column("PiMinus_MC_CosTheta", (Float_t) costheta);
		qa.qaP4("PiMinus_MC_",  l, ntpPiMinus);


		ntpPiMinus->DumpData();
    }


	//Get PiPlus information
    for (int j=0; j<piplus.GetLength(); ++j){

		//general info about event
		ntpPiPlus->Column("ev",     (Float_t) evt);
		ntpPiPlus->Column("cand",    (Float_t) j);
		ntpPiPlus->Column("ncand",   (Float_t) piplus.GetLength());


		//info about 4-vector
		qa.qaP4("PiPlus_",piplus[j]->P4(), ntpPiPlus);
		ntpPiPlus->Column("PiPlus_CosTheta", (Float_t) piplus[j]->GetMomentum().CosTheta());
		qa.qaCand("PiPlus_", piplus[j], ntpPiPlus);

		//get number of Hits in subdetector (only needed if ideal PR is used!!!)
		numberOfHitsInSubdetector("PiPlus_", piplus[j], ntpPiPlus);
		tagNHits("PiPlus_", piplus[j], ntpPiPlus);

		//information about the mother and MCTruth Candidate

		TLorentzVector l;
		float costheta = -999.;

		RhoCandidate * truth = piplus[j]->GetMcTruth();
		RhoCandidate * mother;
		if (truth) mother = truth->TheMother();

		int moth = (mother==0x0) ? 88888 : mother->PdgCode();
		ntpPiPlus->Column("Mother", (Int_t) moth);

		ntpPiPlus->Column("MCTruthMatch", (bool) theAnalysis->McTruthMatch(piplus[j]));

		if (truth){
			l = truth->P4();
			costheta = truth->GetMomentum().CosTheta();
			qa.qaCand("PiMiuns_MC_", truth, ntpPiPlus);
		}
		else{
			qa.qaCand("PiMiuns_MC_", dummyCand, ntpPiPlus);
		}

		ntpPiPlus->Column("PiPlus_MC_CosTheta", (Float_t) costheta);
		qa.qaP4("PiPlus_MC_",  l, ntpPiPlus);


		ntpPiPlus->DumpData();

    }

	//Get Proton information
    for (int j=0; j<proton.GetLength(); j++){

		//general info about event
		ntpProton->Column("ev",     (Float_t) evt);
		ntpProton->Column("cand",    (Float_t) j);
		ntpProton->Column("ncand",   (Float_t) proton.GetLength());


		//info about 4-vector
		qa.qaP4("Proton_",proton[j]->P4(), ntpProton);
		ntpProton->Column("Proton_CosTheta", (Float_t) proton[j]->GetMomentum().CosTheta());
		qa.qaCand("Proton_", proton[j], ntpProton);

		//get number of Hits in subdetector (only needed if ideal PR is used!!!)
		numberOfHitsInSubdetector("Proton_", proton[j], ntpProton);
		tagNHits("Proton_", proton[j], ntpProton);

		//information about the mother and MCTruth Candidate

		TLorentzVector l;
		float costheta = -999.;

		RhoCandidate * truth = proton[j]->GetMcTruth();
		RhoCandidate * mother;
		if (truth) mother = truth->TheMother();

		int moth = (mother==0x0) ? 88888 : mother->PdgCode();
		ntpProton->Column("Mother", (Int_t) moth);

		ntpProton->Column("MCTruthMatch", (bool) theAnalysis->McTruthMatch(proton[j]));

		if (truth){
			l = truth->P4();
			costheta = truth->GetMomentum().CosTheta();
			qa.qaCand("PiMiuns_MC_", truth, ntpProton);
		}
		else{
			qa.qaCand("PiMiuns_MC_", dummyCand, ntpProton);
		}

		ntpProton->Column("Proton_MC_CosTheta", (Float_t) costheta);
		qa.qaP4("Proton_MC_",  l, ntpProton);


		ntpProton->DumpData();
    }
    

		//Get Antiproton
    for (int j=0; j<antiProton.GetLength(); j++){


		//general info about event
		ntpAntiProton->Column("ev",     (Float_t) evt);
		ntpAntiProton->Column("cand",    (Float_t) j);
		ntpAntiProton->Column("ncand",   (Float_t) antiProton.GetLength());


		//info about 4-vector
		qa.qaP4("AntiProton_",antiProton[j]->P4(), ntpAntiProton);
		ntpAntiProton->Column("AntiProton_CosTheta", (Float_t) antiProton[j]->GetMomentum().CosTheta());
		qa.qaCand("AntiProton_", antiProton[j], ntpAntiProton);

		//get number of Hits in subdetector (only needed if ideal PR is used!!!)
		numberOfHitsInSubdetector("AntiProton_", antiProton[j], ntpAntiProton);
		tagNHits("AntiProton_", antiProton[j], ntpAntiProton);

		//information about the mother and MCTruth Candidate

		TLorentzVector l;
		float costheta = -999.;

		RhoCandidate * truth = antiProton[j]->GetMcTruth();
		RhoCandidate * mother;
		if (truth) mother = truth->TheMother();

		int moth = (mother==0x0) ? 88888 : mother->PdgCode();
		ntpAntiProton->Column("Mother", (Int_t) moth);

		ntpAntiProton->Column("MCTruthMatch", (bool) theAnalysis->McTruthMatch(antiProton[j]));

		if (truth){
			l = truth->P4();
			costheta = truth->GetMomentum().CosTheta();
			qa.qaCand("PiMiuns_MC_", truth, ntpAntiProton);
		}
		else{
			qa.qaCand("PiMiuns_MC_", dummyCand, ntpAntiProton);
		}

		ntpAntiProton->Column("AntiProton_MC_CosTheta", (Float_t) costheta);
		qa.qaP4("AntiProton_MC_",  l, ntpAntiProton);


		ntpAntiProton->DumpData();
    }
        
			
    //***Lambda0 -> PiMinus + Proton
    lambda0.Combine(piminus,proton);
	lambda0.Select(lambdaMassSelector);
    lambda0.SetType(3122);

//    //*****sort the candidates according to their fit information
//    std::map<int,int> bestVtxFit, bestMassFitLambda0;
//    bestVtxFit = VertexQaIndex(&lambda0);
//    bestMassFitLambda0 = MassFitQaIndex(&lambda0, m0_lambda0);



    for (int j=0; j<lambda0.GetLength(); ++j){

		//general info about event
		ntpLambda0->Column("ev",     (Float_t) evt);
		ntpLambda0->Column("cand",    (Float_t) j);
		ntpLambda0->Column("ncand",   (Float_t) lambda0.GetLength());
		ntpLambda0->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(lambda0[j]));

		//inforamtion about mother
		RhoCandidate * truth = lambda0[j]->GetMcTruth();
		RhoCandidate * mother;
	 	if (truth) mother = truth->TheMother();
	   	int moth = (mother==0x0) ? 88888 : mother->PdgCode();

	 	ntpLambda0->Column("Mother", (Float_t) moth);

		//information about candidate
		qa.qaP4("Lambda0_", lambda0[j]->P4(), ntpLambda0);
		qa.qaCand("Lambda0_", lambda0[j], ntpLambda0);
		qa.qaComp("Lambda0_", lambda0[j], ntpLambda0);

		//get number of Hits in subdetector (only needed if ideal PR is used!!!)
		numberOfHitsInSubdetector("Lambda0_", lambda0[j], ntpLambda0);

		//check if daughter particles leave more than 3 hit in any inner tracking detector (only needed if ideal PR is used!!!)
		int tag = 0;
		int ndau = lambda0[j]->NDaughters();
		int dtag[2]={0,0};


		for(int dau=0; dau<ndau; dau++){
		  RhoCandidate * daughter = lambda0[j]->Daughter(dau);
		  dtag[dau] = tagHits(daughter);
		}

		if(dtag[0]==1 && dtag[1]==1) tag=1;

		ntpLambda0->Column("Lambda0_HitTag", (Int_t) tag);



		// do vertex fit

		PndKinVtxFitter vertexfitterLambda0 (lambda0[j]);
		vertexfitterLambda0.Fit();
		RhoCandidate * lambda0Fit = lambda0[j]->GetFit();

		//Get information about fit
		qa.qaFitter("VtxFit_", &vertexfitterLambda0, ntpLambda0);
//		ntpLambda0->Column("VtxFit_HowGood", (Int_t)  bestVtxFit[j]);
		qa.qaVtx("VtxFit_", lambda0Fit , ntpLambda0);
		qa.qaCand("VtxFit_", lambda0Fit, ntpLambda0);
		qa.qaP4("VtxFit_", lambda0Fit->P4(), ntpLambda0);
		qa.qaComp("VtxFit_", lambda0Fit, ntpLambda0);


		// differenz to MCTruth
		qa.qaMcDiff("VtxFit_", lambda0Fit, ntpLambda0);
		qaVtxDiff("VtxFit_", lambda0Fit, ntpLambda0);
		qaMomRes("VtxFit_", lambda0Fit, ntpLambda0);


		// do mass fit
		PndKinFitter massFitterLambda0(lambda0Fit);
		massFitterLambda0.AddMassConstraint(m0_lambda0);
		massFitterLambda0.Fit();

		RhoCandidate * lambda0Fit_mass = lambda0Fit->GetFit();
		qa.qaFitter("MassFit_", &massFitterLambda0, ntpLambda0);

//		ntpLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitLambda0[j]);

		//Get MCTruth inforamtion
		TLorentzVector l;

		if(0x0 != truth){
			l = truth->P4();
			qa.qaVtx("McTruth_", truth, ntpLambda0);
			qa.qaComp("McTruth_", truth, ntpLambda0);
		}
		else{
			qa.qaVtx("McTruth_", dummyCand, ntpLambda0);
			qa.qaComp("McTruth_", dummyCand, ntpLambda0);
		}

		qa.qaP4("McTruth_", l, ntpLambda0);



//		if (bestVtxFit[j]==1 && bestMassFitLambda0[j]>0 && tag==1){ //*** use only bestChi2Cand
		if(vertexfitterLambda0.GetProb()>0.01 && massFitterLambda0.GetProb()>0.01 /* && tag==1 */){ //*** use only candidates with prob>0.01 in both fitter
		  Lambda0Fit.Append(lambda0Fit);
		}

	  //***information of boosted particle
	  lambda0Fit->Boost(-beamBoost);
	  qa.qaComp("boost_", lambda0Fit, ntpLambda0);

	  ntpLambda0->DumpData();
  }



    //***AntiLambda0 -> PiPlus + AntiProton
    antiLambda0.Combine(piplus,antiProton);
	antiLambda0.Select(lambdaMassSelector);
    antiLambda0.SetType(-3122);

//    //*****sort the candidates according to their fit information
//    std::map<int,int> bestVtxFitAntiLambda0, bestMassFitAntiLambda0;
//    bestVtxFitAntiLambda0 = VertexQaIndex(&antiLambda0);
//    bestMassFitAntiLambda0 = MassFitQaIndex(&antiLambda0, m0_lambda0);



    for (int j=0; j<antiLambda0.GetLength(); ++j){

		//general info about event
		ntpAntiLambda0->Column("ev",     (Float_t) evt);
		ntpAntiLambda0->Column("cand",    (Float_t) j);
		ntpAntiLambda0->Column("ncand",   (Float_t) lambda0.GetLength());
		ntpAntiLambda0->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(antiLambda0[j]));

		//inforamtion about mother
		RhoCandidate * mother;
		RhoCandidate * truth = antiLambda0[j]->GetMcTruth();
		if (truth) mother = truth->TheMother();
		int moth = (mother==0x0) ? 88888 : mother->PdgCode();

		ntpAntiLambda0->Column("Mother", (Float_t) moth);

		//information about candidate
		qa.qaP4("Lambda0_", antiLambda0[j]->P4(), ntpAntiLambda0);
		qa.qaCand("Lambda0_", antiLambda0[j], ntpAntiLambda0);
		qa.qaComp("Lambda0_", antiLambda0[j], ntpAntiLambda0);

		//get number of Hits in subdetector (only needed if ideal PR is used!!!)
		numberOfHitsInSubdetector("Lambda0_", antiLambda0[j], ntpAntiLambda0);

		//check if daughter particles leave more than 3 hit in any inner tracking detector (only needed if ideal PR is used!!!)
		int tag = 0;
		int ndau = antiLambda0[j]->NDaughters();
		int dtag[2]={0,0};


		for(int dau=0; dau<ndau; dau++){
		  RhoCandidate * daughter = antiLambda0[j]->Daughter(dau);
		  dtag[dau] = tagHits(daughter);
		}

		if(dtag[0]==1 && dtag[1]==1) tag=1;

		ntpAntiLambda0->Column("Lambda0_HitTag", (Int_t) tag);



		// do vertex fit

		PndKinVtxFitter vertexfitterAntiLambda0 (antiLambda0[j]);
		vertexfitterAntiLambda0.Fit();
		RhoCandidate * antiLambda0Fit = antiLambda0[j]->GetFit();

		//Get information about fit
		qa.qaFitter("VtxFit_", &vertexfitterAntiLambda0, ntpAntiLambda0);
//		ntpAntiLambda0->Column("VtxFit_HowGood", (Int_t)  bestVtxFitAntiLambda0[j]);
		qa.qaVtx("VtxFit_", antiLambda0Fit , ntpAntiLambda0);
		qa.qaCand("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
		qa.qaP4("VtxFit_", antiLambda0Fit->P4(), ntpAntiLambda0);
		qa.qaComp("VtxFit_", antiLambda0Fit, ntpAntiLambda0);


		// differenz to MCTruth
		qa.qaMcDiff("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
		qaVtxDiff("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
		qaMomRes("VtxFit_", antiLambda0Fit, ntpAntiLambda0);


		// do mass fit
		PndKinFitter massFitterAntiLambda0(antiLambda0Fit);
		massFitterAntiLambda0.AddMassConstraint(m0_lambda0);
		massFitterAntiLambda0.Fit();

		RhoCandidate * lambda0Fit_mass = antiLambda0Fit->GetFit();
		qa.qaFitter("MassFit_", &massFitterAntiLambda0, ntpAntiLambda0);

//		ntpAntiLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitAntiLambda0[j]);


		//Get MCTruth information
		TLorentzVector l;

		if(0x0 != truth){
			l = truth->P4();
			qa.qaVtx("McTruth_", truth, ntpAntiLambda0);
			qa.qaComp("McTruth_", truth, ntpAntiLambda0);
		}
		else{
			qa.qaVtx("McTruth_", dummyCand, ntpAntiLambda0);
			qa.qaComp("McTruth_", dummyCand, ntpAntiLambda0);
		}

		qa.qaP4("McTruth_", l, ntpAntiLambda0);



//		if (bestVtxFitAntiLambda0[j]==1 && bestMassFitAntiLambda0[j]>0 && tag==1){ //*** use only bestChi2Cand
		if (vertexfitterAntiLambda0.GetProb()>0.01 &&  massFitterAntiLambda0.GetProb()>0.01 /* && tag==1 */){ //*** use only candidate with prob>0.01 in both fitter
			AntiLambda0Fit.Append(antiLambda0Fit);
		}


	  //***information of boosted particle
	  antiLambda0Fit->Boost(-beamBoost);
	  qa.qaComp("boost_", antiLambda0Fit, ntpAntiLambda0);

	  ntpAntiLambda0->DumpData();
    }



    //*** pbar + p -> Lambda0 + AntiLambda0

	sys.Combine(Lambda0Fit, AntiLambda0Fit);
	sys.SetType(88888);


	for (int j=0; j<sys.GetLength(); ++j){

		//general information about event
		ntpSys->Column("ev",     (Float_t) evt);
		ntpSys->Column("cand",    (Float_t) j);
		ntpSys->Column("ncand",   (Float_t) sys.GetLength());
		ntpSys->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(sys[j]));

		qa.qaP4("", sys[j]->P4(), ntpSys);
		qa.qaComp("", sys[j], ntpSys);
		qa.qaPoca("", sys[j], ntpSys);
		qa.qaEventShapeShort("es_", &evsh, ntpAntiLambda0);


		//do vertex fit

		PndKinVtxFitter vertexFitter_cc (sys[j]);
		vertexFitter_cc.Fit();
		RhoCandidate * ccFit = sys[j]->GetFit();


		//***do 4c fit
		PndKinFitter cc_Fitter4c (sys[j]);
		cc_Fitter4c.Add4MomConstraint(ini);
		cc_Fitter4c.Fit();

		RhoCandidate * ccFit = sys[j]->GetFit();

		// store info of 4c Fit
		ntpSys->Column("f4c_Chi2_cc", (Float_t) cc_Fitter4c.GetChi2());
		ntpSys->Column("f4c_NDF_cc", (Float_t) cc_Fitter4c.GetNdf());
		ntpSys->Column("f4c_Prob_cc", (Float_t) cc_Fitter4c.GetProb());

		qa.qaComp("f4c_", ccFit, ntpSys);

		//difference to MC Truth
		qa.qaMcDiff("f4c_", ccFit, ntpSys);
		qaVtxDiff("f4c_", ccFit, ntpSys);
		qaMomRes("f4c_", ccFit, ntpSys);


		ntpSys->DumpData();
	}
	 Lambda0Fit.Cleanup();
	 AntiLambda0Fit.Cleanup();
  }


  //Write output
  out->cd();

  ntpMC -> GetInternalTree()->Write();
  ntpPiMinus ->GetInternalTree()->Write();
  ntpPiPlus->GetInternalTree()->Write();
  ntpProton->GetInternalTree()->Write();
  ntpAntiProton->GetInternalTree()->Write();
  ntpLambda0->GetInternalTree()->Write();
  ntpAntiLambda0->GetInternalTree()->Write();
  ntpSys->GetInternalTree()->Write();

  out->Save();

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();

  cout<<"Macro finisched successfully."<<endl;
  cout<<"Realtime: "<<rtime<<" s, CPU time: "<<ctime<<" s"<<endl;
  cout<<endl;


  exit(0);
 
}
