class PndAnalysis;
class PndAnaPidSelector;
class RhoCandList;
class RhoTuple;

#include "common_jenny.cpp"


void analysis_pbarp_lambda0(TString pre="", int nevts=0, double mom = 1.7){
  
  TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);
  
  TStopwatch timer; 


  //Output File
  TString OutputFile = pre + "_ana.root";

  //Input simulation Files
  TString inPIDFile = pre + "_pid.root";
  TString inParFile = pre + "_par.root";
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

  //*** create tuples
  RhoTuple * ntpMC = new RhoTuple("ntpMC", "MCTruth info");
  RhoTuple * ntpPiMinus = new RhoTuple("ntpPiMinus", "PiMinus info");
  RhoTuple * ntpPiPlus = new RhoTuple("ntpPiPlus", "PiPlus info");
  RhoTuple * ntpProton = new RhoTuple("ntpProton", "Proton info");
  RhoTuple * ntpAntiProton = new RhoTuple("ntpAntiProton", "Antiproton info");
  RhoTuple * ntpLambda0 = new RhoTuple("ntpLambda0", "Lambda0 info");
  RhoTuple * ntpAntiLambda0 = new RhoTuple("ntpAntiLambda0", "AntiLambda0 info");
  RhoTuple * ntpCrossCheck = new RhoTuple("ntpCrossCheck", "CrossCheck info");

  //Create output file 
  TFile *out = TFile::Open(pre+"_output_ana.root","RECREATE");

  // data reader Object
  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts = theAnalysis->GetEntries();
  
  //RhoCandLists for analysis
  RhoCandList piplus, piminus, lambda0, antiLambda0, proton, antiProton, crossCheck, Lambda0Fit, AntiLambda0Fit;
  RhoCandList mclist, all;

  RhoCandidate * dummyCand = new RhoCandidate();

  //***Mass selector
  double m0_lambda0= TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
  cout<<"Mass of Lambda0: "<<m0_lambda0<<endl;
	RhoMassParticleSelector * lambdaMassSelector = new RhoMassParticleSelector("lambda0", m0_lambda0, 0.3);   
 
	double m0_pbarpsystem = TDatabasePDG::Instance()->GetParticle("pbarpSystem")->Mass();
  
  double pbarmom = mom;
  double p_m0 = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
  TLorentzVector ini (0,0, pbarmom, sqrt(p_m0*p_m0+ pbarmom*pbarmom)+p_m0);
  TVector3 beamBoost = ini.BoostVector();
  
  PndRhoTupleQA qa(theAnalysis, pbarmom);
 
	int evt=0;
  while (theAnalysis->GetEvent() && ++evt<nevts){

    if ((evt%100)==0) cout << "evt "<< evt <<endl;
    cout << "Event number: " << evt << endl;

    //***get MC list and store info
    theAnalysis->FillList(mclist, "McTruth");
    qa.qaMcList("", mclist, ntpMC);
    ntpMC->DumpData();
   
		
		//if you want to print out the hole MCTree uncomment the following
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

    TString PidSelection = "PidAlgoIdealCharged";//"PidAlgoMvd;PidAlgoStt;PidAlgoDrc";

    theAnalysis->FillList(all, "All", PidSelection);
    PndEventShape evsh(all, ini, 0.05, 0.1);
    
    //***Selection with no PID info
    theAnalysis->FillList(piminus, "PionAllMinus", PidSelection);
    theAnalysis->FillList(piplus, "PionAllPlus", PidSelection);
    theAnalysis->FillList(proton, "ProtonAllPlus", PidSelection);
    theAnalysis->FillList(antiProton, "ProtonAllMinus", PidSelection);

    //Get PiMinusanalysis/pbarp_lambda0_antilambda0/10000_events/idealtracking
    for (int j=0; j<piminus.GetLength(); ++j){
			
			//general info about event
      ntpPiMinus->Column("ev",     (Float_t) evt);
      ntpPiMinus->Column("cand",    (Float_t) j);
      ntpPiMinus->Column("ncand",   (Float_t) piminus.GetLength());
      ntpPiMinus->Column("MCTruthMatch", (bool) theAnalysis->McTruthMatch(piminus[j]));
      
      //info about 4-vector
      qa.qaP4("PiMinus_",piminus[j]->P4(), ntpPiMinus);
	  qa.qaCand("PiMinus_", piminus[j], ntpPiMinus);
	  ntpPiMinus->Column("PiMinus_CosTheta", (Float_t) piminus[j]->GetMomentum().CosTheta());


	  RhoCandidate * mother = piminus[j]->TheMother();
      TLorentzVector l;

      RhoCandidate * truth = piminus[j]->GetMcTruth();
      
      double costheta = -999.;

      if (truth != 0x0){
    	  l=truth->P4();
    	  costheta = truth->GetMomentum().CosTheta();
      }

      int moth = (mother==0x0) ? 88888 : mother->PdgCode();
      ntpPiMinus->Column("Mother", (Int_t) moth);
      qa.qaP4("PiMinus_MC_", l, ntpPiMinus);
      ntpPiMinus->Column("PiMinus_MC_CosTheta", (Double_t) costheta);
			
      ntpPiMinus->DumpData();
    }

		//Get PiPlus
    for (int j=0; j<piplus.GetLength(); ++j){

      //general info about event
      ntpPiPlus->Column("ev",     (Float_t) evt);
      ntpPiPlus->Column("cand",    (Float_t) j);
      ntpPiPlus->Column("ncand",   (Float_t) piplus.GetLength());
      ntpPiPlus->Column("MCTruthMatch", (bool) theAnalysis->McTruthMatch(piplus[j]));

      //info about 4-vector
      qa.qaP4("PiPlus_", piplus[j]->P4(), ntpPiPlus);
      ntpPiPlus->Column("PiPlus_CosTheta", (Float_t) piplus[j]->GetMomentum().CosTheta());
      qa.qaCand("PiPlus_", piplus[j], ntpPiPlus);

      RhoCandidate * mother = piplus[j]->TheMother();
      TLorentzVector l;

      RhoCandidate * truth = piplus[j]->GetMcTruth();

      double costheta = -999.;

      if (truth != 0x0){
    	  l=truth->P4();
    	  costheta = truth->GetMomentum().CosTheta();
      }

      int moth = (mother==0x0) ? 88888 : mother->PdgCode();
      ntpPiPlus->Column("Mother", (Int_t) moth);
      qa.qaP4("PiPlus_MC_", l, ntpPiPlus);
      ntpPiPlus->Column("PiPlus_MC_CosTheta", (Double_t) costheta);

			
      ntpPiPlus->DumpData();

    }

		//Get Proton
	for (int j=0; j<proton.GetLength(); j++){

		//general info about event
		ntpProton->Column("ev",     (Float_t) evt);
		ntpProton->Column("cand",    (Float_t) j);
		ntpProton->Column("ncand",   (Float_t) proton.GetLength());
		ntpProton->Column("MCTruthMatch", (bool) theAnalysis->McTruthMatch(proton[j]));

		//info about 4-vector
		qa.qaP4("Proton_", proton[j]->P4(), ntpProton);
		qa.qaCand("Proton_", proton[j], ntpProton);
		ntpProton->Column("Proton_CosTheta", (Float_t) proton[j]->GetMomentum().CosTheta());

		RhoCandidate * mother = proton[j]->TheMother();
		TLorentzVector l;

		RhoCandidate * truth = proton[j]->GetMcTruth();

		double costheta = -999.;

		if (truth != 0x0){
			l=truth->P4();
			costheta = truth->GetMomentum().CosTheta();
		}

		int moth = (mother==0x0) ? 88888 : mother->PdgCode();
		ntpProton->Column("Mother", (Int_t) moth);
		qa.qaP4("Proton_MC_", l, ntpProton);
		ntpProton->Column("Proton_MC_CosTheta", (Double_t) costheta);
		
		ntpProton->DumpData();
}
    

		//Get Antiproton
	for (int j=0; j<antiProton.GetLength(); j++){

		//general info about event
		ntpAntiProton->Column("ev",     (Float_t) evt);
		ntpAntiProton->Column("cand",    (Float_t) j);
		ntpAntiProton->Column("ncand",   (Float_t) antiProton.GetLength());
		ntpAntiProton->Column("MCTruthMatch", (bool) theAnalysis->McTruthMatch(antiProton[j]));

		//info about 4-vector
		qa.qaP4("AntiProton_", antiProton[j]->P4(), ntpAntiProton);
		qa.qaCand("AntiProton_", antiProton[j], ntpAntiProton);
		ntpAntiProton->Column("AntiProton_CosTheta", (Float_t) antiProton[j]->GetMomentum().CosTheta());

		RhoCandidate * mother = antiProton[j]->TheMother();
		TLorentzVector l;

		RhoCandidate * truth = antiProton[j]->GetMcTruth();

		double costheta = -999.;

		if (truth != 0x0){
			l=truth->P4();
			costheta = truth->GetMomentum().CosTheta();
		}
		
		int moth = (mother==0x0) ? 88888 : mother->PdgCode();
		ntpAntiProton->Column("Mother", (Int_t) moth);
		qa.qaP4("AntiProton_MC_", l, ntpAntiProton);
		ntpAntiProton->Column("AntiProton_MC_CosTheta", (Double_t) costheta);
			
		ntpAntiProton->DumpData();
}
        
			
    //***Lambda0 -> PiMinus + Proton
    lambda0.Combine(piminus,proton);
	lambda0.Select(lambdaMassSelector);
    lambda0.SetType(3122);

//    std::map<int,int> bestVtxFit, bestMassFitLambda0;
//    bestVtxFit = jenny::VertexQaIndex(&lambda0);
//    bestMassFitLambda0 = jenny::MassFitQaIndex(&lambda0, m0_lambda0);



    for (int j=0; j<lambda0.GetLength(); ++j){

      
      //general info about event
      ntpLambda0->Column("ev",     (Float_t) evt);
      ntpLambda0->Column("cand",    (Float_t) j);
      ntpLambda0->Column("ncand",   (Float_t) lambda0.GetLength());
      ntpLambda0->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(lambda0[j]));


			RhoCandidate * mother = lambda0[j]->TheMother();
  	 	int moth = (mother==0x0) ? 88888 : mother->PdgCode();

			ntpLambda0->Column("Mother", (Float_t) moth);
      
      qa.qaP4("Lambda0_", lambda0[j]->P4(), ntpLambda0);
      qa.qaCand("Lambda0_", lambda0[j], ntpLambda0);
      qa.qaComp("Lambda0_", lambda0[j], ntpLambda0);
//      qa.qaEventShapeShort("es_", &evsh, ntpLambda0);
			


      // do vertex fit

      PndKinVtxFitter vertexfitterLambda0 (lambda0[j]);
      vertexfitterLambda0.Fit();
      RhoCandidate * lambda0Fit = lambda0[j]->GetFit();

      qa.qaFitter("VtxFit_", &vertexfitterLambda0, ntpLambda0);
//      ntpLambda0->Column("VtxFit_HowGood", (Int_t)  bestVtxFit[j]);
	  qa.qaVtx("VtxFit_", lambda0Fit , ntpLambda0);
	  qa.qaCand("VtxFit_", lambda0Fit, ntpLambda0);
	  qa.qaP4("VtxFit_", lambda0Fit->P4(), ntpLambda0);
	  qa.qaComp("VtxFit_", lambda0Fit, ntpLambda0);

      // differenz to MCTruth
      qa.qaMcDiff("VtxFit_", lambda0Fit, ntpLambda0);
      jenny::qaVtxDiff("VtxFit_", lambda0Fit, ntpLambda0);
      jenny::qaMomRes("VtxFit_", lambda0Fit, ntpLambda0);





      // do mass fit
      PndKinFitter massFitterLambda0(lambda0Fit);
      massFitterLambda0.AddMassConstraint(m0_lambda0);
      massFitterLambda0.Fit();

      RhoCandidate * lambda0Fit_mass = lambda0Fit->GetFit();

      qa.qaFitter("MassFit_", &massFitterLambda0, ntpLambda0);
//      ntpLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitLambda0[j]);
      qa.qaMcDiff("MassFit_", lambda0Fit_mass, ntpLambda0);
      jenny::qaVtxDiff("VtxFit_", lambda0Fit, ntpLambda0);
      jenny::qaMomRes("VtxFit_", lambda0Fit, ntpLambda0);


      // store only best fitted candidate
//      if( bestVtxFit[j]==1){// && bestMassFitLambda0[j]>0){
    	  Lambda0Fit.Append(lambda0Fit);
//      }

      RhoCandidate * truth = lambda0[j]->GetMcTruth();
      TLorentzVector l;



	    if(0x0 != truth){
	    	l = truth->P4();
	    	qa.qaVtx("McTruth_", truth, ntpLambda0);
	    	dl = truth->Daughter(0)->Pos();
	    }
	    else{
	    	qa.qaVtx("McTruth_", dummyCand, ntpLambda0);
	    }


      qa.qaP4("McTruth_", l, ntpLambda0);




	  //***information of boosted particle
	  lambda0Fit->Boost(-beamBoost);
	  qa.qaComp("boost_", lambda0Fit, ntpLambda0);
    		
	  ntpLambda0->DumpData();
    }



    //***AntiLambda0 -> PiMinus + Proton
    antiLambda0.Combine(piplus,antiProton);
	antiLambda0.Select(lambdaMassSelector);
    antiLambda0.SetType(-3122);

//    std::map<int,int> bestVtxFitAntiLambda0, bestMassFitAntiLambda0;
//    bestVtxFitAntiLambda0 = jenny::VertexQaIndex(&antiLambda0);
//    bestMassFitAntiLambda0 = jenny::MassFitQaIndex(&antiLambda0, m0_lambda0);




//    float bestchi2=9999;

    for (int j=0; j<antiLambda0.GetLength(); ++j){

      //general info about event
      ntpAntiLambda0->Column("ev",     (Float_t) evt);
      ntpAntiLambda0->Column("cand",    (Float_t) j);
      ntpAntiLambda0->Column("ncand",   (Float_t) antiLambda0.GetLength());
      ntpAntiLambda0->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(antiLambda0[j]));

			RhoCandidate * mother = antiLambda0[j]->TheMother();
			int moth = (mother==0x0) ? 88888 : mother->PdgCode();
			ntpAntiLambda0->Column("Mother", (Float_t) moth);

      qa.qaP4("AntiLambda0_", antiLambda0[j]->P4(), ntpAntiLambda0);
      qa.qaComp("AntiLambda0_", antiLambda0[j], ntpAntiLambda0);
      qa.qaEventShapeShort("es_", &evsh, ntpAntiLambda0);




      // do vertex fit


      PndKinVtxFitter vertexfitterAntiLambda0 (antiLambda0[j]);
      vertexfitterAntiLambda0.Fit();
      RhoCandidate * antiLambda0Fit = antiLambda0[j]->GetFit();


      // store info of vertex fit


      qa.qaFitter("VtxFit_", &vertexfitterAntiLambda0, ntpAntiLambda0);
//      ntpAntiLambda0->Column("VtxFit_HowGood", (Int_t) bestVtxFitAntiLambda0[j]);

	  qa.qaVtx("VtxFit_", antiLambda0Fit, ntpAntiLambda0);


	  // difference to MCTruth
      qa.qaMcDiff("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
      jenny::qaVtxDiff("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
      jenny::qaMomRes("VtxFit_", antiLambda0Fit, ntpAntiLambda0);



//      if(prob>0.01 && chi2<bestchi2){
//    	  bestchi2 = chi2;
//    	  AntiLambda0Fit.InsertAt(0,antiLambda0Fit);
//      }


      // do mass fit
      PndKinFitter massFitterAntiLambda0(antiLambda0Fit);
      massFitterAntiLambda0.AddMassConstraint(m0_lambda0);
      massFitterAntiLambda0.Fit();

      RhoCandidate * antiLambda0Fit_mass = antiLambda0Fit->GetFit();

      qa.qaFitter("MassFit_", &massFitterAntiLambda0, ntpAntiLambda0);
//      ntpAntiLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitAntiLambda0[j]);
      qa.qaMcDiff("MassFit_", antiLambda0Fit_mass, ntpAntiLambda0);


//      if(bestVtxFitAntiLambda0[j]==1){ //&& bestMassFitAntiLambda0[j]>0){
    	  AntiLambda0Fit.Append(antiLambda0Fit);
//      }

      RhoCandidate * truth = antiLambda0[j]->GetMcTruth();
      TLorentzVector l;
      if(0x0 != truth){
				l = truth->P4();
				qa.qaVtx("truth_", truth, ntpAntiLambda0);
      }
      qa.qaP4("truth_", l, ntpAntiLambda0);

			//***information of boosted particle
			antiLambda0Fit->Boost(-beamBoost);
			qa.qaComp("boost_", antiLambda0Fit, ntpAntiLambda0);


      ntpAntiLambda0->DumpData();
    }



		//*** Cross check: pbar + p -> Lambda0 + AntiLambda0

		crossCheck.Combine(Lambda0Fit, AntiLambda0Fit);
		crossCheck.SetType(88888);


		for (int j=0; j<crossCheck.GetLength(); ++j){

			//general information about event
			ntpCrossCheck->Column("ev",     (Float_t) evt);
			ntpCrossCheck->Column("cand",    (Float_t) j);
			ntpCrossCheck->Column("ncand",   (Float_t) crossCheck.GetLength());
			ntpCrossCheck->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(crossCheck[j]));

			qa.qaP4("", crossCheck[j]->P4();, ntpCrossCheck);
			qa.qaComp("", crossCheck[j], ntpCrossCheck);
			qa.qaPoca("", crossCheck[j], ntpCrossCheck);
			qa.qaEventShapeShort("es_", &evsh, ntpAntiLambda0);


			//do vertex fit



			PndKinVtxFitter vertexFitter_cc (crossCheck[j]);
			vertexFitter_cc.Fit();
			RhoCandidate * ccFit = crossCheck[j]->GetFit();

			//store info of vertex fit
			qa.qaFitter("VtxFit_", &vertexFitter_cc, ntpCrossCheck);
			qa.qaVtx("VtxFit_", ccFit, ntpCrossCheck);
			qa.qaMcDiff("VtxFit_", ccFit, ntpCrossCheck);
			jenny::qaVtxDiff("VtxFit_", ccFit, ntpCrossCheck);
			jenny::qaMomRes("VtxFit_", ccFit, ntpCrossCheck);



			float mc_mass_l0 = 0.;

			RhoCandidate * truth = ccFit->GetMcTruth();
		    TLorentzVector l;


		    if(0x0 != truth){
		    	l = truth->P4();
		    	qa.qaVtx("McTruth_", truth, ntpCrossCheck);
		    }
		    else{
		    	qa.qaVtx("McTruth_", dummyCand, ntpCrossCheck);
		    }

  		    qa.qaP4("McTruth_", l, ntpCrossCheck);


			//***do 4c fit
			PndKinFitter cc_Fitter4c (crossCheck[j]);
			cc_Fitter4c.Add4MomConstraint(ini);
			cc_Fitter4c.Fit();

			RhoCandidate * ccFit = crossCheck[j]->GetFit();

			// store info of 4c Fit
			ntpCrossCheck->Column("f4c_Chi2_cc", (Float_t) cc_Fitter4c.GetChi2());
			ntpCrossCheck->Column("f4c_NDF_cc", (Float_t) cc_Fitter4c.GetNdf());
			ntpCrossCheck->Column("f4c_Prob_cc", (Float_t) cc_Fitter4c.GetProb());

			qa.qaComp("f4c_", ccFit, ntpCrossCheck);

			//difference to MC Truth
			qa.qaMcDiff("f4c_", ccFit, ntpCrossCheck);
			jenny::qaVtxDiff("f4c_", ccFit, ntpCrossCheck);
			jenny::qaMomRes("f4c_", ccFit, ntpCrossCheck);


			//do kalman vertex fit

			PndKalmanVtxFitter kalmanfitter (crossCheck[j]);
			kalmanfitter.Fit();
			RhoCandidate * kalmanFit = crossCheck[j]->GetFit();

			qa.qaFitter("KalmanFit_", &kalmanfitter, ntpCrossCheck);
			qa.qaVtx("KalmanFit_", kalmanFit, ntpCrossCheck);

			jenny::qaVtxDiff("KalmanFit_", kalmanFit, ntpCrossCheck);
			jenny::qaMomRes("KalmanFit_", kalmanFit, ntpCrossCheck);

			//do mass fit
			PndKinFitter massFitter_cc (ccFit);
			massFitter_cc.AddMassConstraint(m0_pbarpsystem);
			massFitter_cc.Fit();

			RhoCandidate * massFit_cc = ccFit->GetFit();

			//store fit results
			ntpCrossCheck->Column("fmass_Chi2_cc", (Float_t) massFitter_cc.GetChi2());
			ntpCrossCheck->Column("fmass_NDF_cc", (Float_t) massFitter_cc.GetNdf());
			ntpCrossCheck->Column("fmass_Prob_cc", (Float_t) massFitter_cc.GetProb());

			qa.qaMcDiff("fmassMCDiff_", massFit_cc, ntpCrossCheck);

			ntpCrossCheck->DumpData();
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
  ntpCrossCheck->GetInternalTree()->Write();

  out->Save();
  
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();

  cout<<"Macro finisched successfully."<<endl;
  cout<<"Realtime: "<<rtime<<" s, CPU time: "<<ctime<<" s"<<endl;
  cout<<endl;

 
  exit(0);
 
}
