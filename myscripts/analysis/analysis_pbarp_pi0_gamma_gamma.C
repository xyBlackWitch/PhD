class PndAnalysis;
class PndAnaPidSelector;
class RhoCandList;
class RhoTuple;

#include "../common_jenny.cpp"


void analysis_pbarp_pi0_gamma_gamma(int nevts=0, float mom=3){
  
  TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);
  
  TStopwatch timer; 

  //Output File
  TString Path ="/home/ikp1/puetz/panda/mysimulations/test/boxgenerator/pi0_gamma_gamma/";
  TString outPath ="";//"/home/ikp1/puetz/panda/mysimulations/test/boxgenerator/pi0_gamma_gamma/";
 
  TString OutputFile = outPath + "analysis_output.root";
  
  //Input simulation Files
  TString inPIDFile = Path + "Pion_pid_complete.root";
  TString inParFile = Path + "Pion_simparams.root";
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
  RhoTuple * ntpGamma = new RhoTuple("ntpGamma", "Gamma info");
  RhoTuple * ntpPiNull = new RhoTuple("ntpPiNull", "PiNull info");

  //Create output file for histograms
  TFile *out = TFile::Open(outPath+"output_ana.root","RECREATE");


  // data reader Object
  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts = theAnalysis->GetEntries();
  
  //RhoCandLists for analysis
  RhoCandList pinull, gamma;
  RhoCandList mclist, all;

  //***Mass selector
  double m0_pi0 = TDatabasePDG::Instance()->GetParticle("pi0")->Mass();
  cout<<"Mass pbar pi0 system: "<<m0_pi0<<endl;
  
  Pi0MassSelector = new RhoMassParticleSelector("pi0", m0_pi0, 0.3);

  int evt = 0;
  //*** lorentz vector of the initial particle
  
  
  PndRhoTupleQA qa(theAnalysis, mom);
  //int nscattered = 0;

  RhoCandidate * dummyCand = new RhoCandidate();

  while (theAnalysis->GetEvent() && evt++<nevts){
    if ((evt%100)==0) cout << "evt "<< evt <<endl;

    //***get MC list and store info
    theAnalysis->FillList(mclist, "McTruth");
    qa.qaMcList("", mclist, ntpMC);
    ntpMC->DumpData();

    //***Setup event shape object

	TString PidSelection = "PidAlgoIdealNeutral";//"PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoDisc";

    theAnalysis->FillList(all, "All", PidSelection);
    
    //***Selection with no PID info
    theAnalysis->FillList(gamma, "Neutral", PidSelection);

    //Momentum histograms
	
    for (int k=0; k<gamma.GetLength(); ++k){
		ntpGamma->Column("ev",     (Float_t) evt);
		ntpGamma->Column("cand",    (Float_t) k);
		ntpGamma->Column("ncand",   (Float_t) gamma.GetLength());
		ntpGamma->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(gamma[k]));

		qa.qaP4("gamma_", gamma[k]->P4(), ntpGamma);
		qa.qaCand("gamma_", gamma[k], ntpGamma);

		jenny::numberOfHitsInSubdetector("gamma_", gamma[k], ntpGamma);
		jenny::tagNHits("gamma_", gamma[k], ntpGamma);

		RhoCandidate * mother_k =0;
		RhoCandidate * truth = gamma[k]->GetMcTruth();


		TLorentzVector l;
		float costheta = -999.;
		if(truth!=0x0){
		  l=truth->P4();
		  mother_k =  truth->TheMother();
		  costheta = truth->GetMomentum().CosTheta();
		  qa.qaCand("gamma_MC_", gamma[k]->GetMcTruth(), ntpGamma);
		}
		else{
		  qa.qaCand("gamma_MC_", dummyCand, ntpGamma);
		}

		int moth_pip = (mother_k==0x0)? 88888 : mother_k->PdgCode();
		ntpGamma->Column("Mother", (Float_t) moth_pip);

		qa.qaP4("gamma_MC_", l, ntpGamma);
		ntpGamma->Column("gamma_MC_CosTheta", (Float_t) costheta);

		ntpGamma->DumpData();
	}

    //***Pi0 -> gamma + gamma
    pinull.Combine(gamma, gamma);
    pinull.Select(Pi0MassSelector);
    pinull.SetType(111);

    std::map<int,int> bestVtxFitpinull, bestMassFitpinull;

    bestVtxFitpinull = jenny::VertexQaIndex(&pinull);
    bestMassFitpinull = jenny::MassFitQaIndex(&pinull, m0_pi0);


    for (int j=0; j<pinull.GetLength(); ++j){


      //general info about event
      ntpPiNull->Column("ev",     (Float_t) evt);
      ntpPiNull->Column("cand",    (Float_t) j);
      ntpPiNull->Column("ncand",   (Float_t) pinull.GetLength());
      ntpPiNull->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(pinull[j]));
      ntpPiNull->Column("pinull_Pdg", (Float_t) pinull[j]->PdgCode());

	  RhoCandidate * mother = pinull[j]->TheMother();
  	  int moth = (mother==0x0) ? 88888 : mother->PdgCode();

	  ntpPiNull->Column("Mother", (Float_t) moth);

      qa.qaP4("pinull_", pinull[j]->P4(), ntpPiNull);
      qa.qaComp("pinull_", pinull[j], ntpPiNull);

      int tag = 0;
      int ndau = pinull[j]->NDaughters();
      int dtag[2]={0,0};


      for(int dau=0; dau<ndau; dau++){
    	  RhoCandidate * daughter = pinull[j]->Daughter(dau);
    	  dtag[dau] = jenny::tagHits(daughter);
      }

      if(dtag[0]==1 && dtag[1]==1) tag=1;


      ntpPiNull->Column("pinull_HitTag", (Int_t) tag);


      // do vertex fit
      PndKinVtxFitter vertexfitterpinull (pinull[j]);
      vertexfitterpinull.Fit();
      RhoCandidate * pinullFit = pinull[j]->GetFit();


      // store info of vertex fit
      qa.qaFitter("VtxFit_", &vertexfitterpinull, ntpPiNull);
      ntpPiNull->Column("VtxFit_HowGood", (Int_t) bestVtxFitpinull[j]);
      qa.qaVtx("VtxFit_", pinullFit, ntpPiNull);
      qa.qaComp("VtxFit_", pinullFit, ntpPiNull);
      qa.qaP4("VtxFit_", pinullFit->P4(), ntpPiNull);

      // differenz to MCTruth
       qa.qaMcDiff("VtxFit_", pinullFit, ntpPiNull);
       jenny::qaVtxDiff("VtxFit_", pinullFit, ntpPiNull);
       jenny::qaMomRes("VtxFit_", pinullFit, ntpPiNull);


	   // do Kalman vertex fit
	   PndKalmanVtxFitter kalmanfitterpinull (pinull[j]);
	   kalmanfitterpinull.Fit();
	   RhoCandidate * pinullKalmanFit = pinull[j]->GetFit();


	   // store info of vertex fit
	   qa.qaFitter("KalmanFit_", &kalmanfitterpinull, ntpPiNull);
	   //ntpPiNull->Column("KalmanFit_HowGood", (Int_t) bestKalmanFitpinull[j]);
	   qa.qaVtx("KalmanFit_", pinullKalmanFit, ntpPiNull);
	   qa.qaComp("KalmanFit_", pinullKalmanFit, ntpPiNull);
	   qa.qaP4("KalmanFit_", pinullKalmanFit->P4(), ntpPiNull);

	   // differenz to MCTruth
		qa.qaMcDiff("KalmanFit_", pinullKalmanFit, ntpPiNull);
		jenny::qaVtxDiff("KalmanFit_", pinullKalmanFit, ntpPiNull);
		jenny::qaMomRes("KalmanFit_", pinullKalmanFit, ntpPiNull);

      // do mass fit
      PndKinFitter massFitterpinull(pinullFit);
      massFitterpinull.AddMassConstraint(m0_pi0);
      massFitterpinull.Fit();

      RhoCandidate * pinullFit_mass = pinullFit->GetFit();
      qa.qaFitter("MassFit_", &massFitterpinull, ntpPiNull);

      ntpPiNull->Column("MassFit_HowGood", (Int_t) bestMassFitpinull[j]);


      RhoCandidate * truth = pinull[j]->GetMcTruth();
      RhoCandidate * truthDaughter = pinull[j]->Daughter(0)->GetMcTruth();
      TLorentzVector l;
      TVector3 dl;

	    if(0x0 != truth){
	    	l = truth->P4();
	    	qa.qaVtx("McTruth_", truth, ntpPiNull);
	    	dl = truth->Daughter(0)->Pos();
	    }
	    else{
	    	qa.qaVtx("McTruth_", dummyCand, ntpPiNull);
	    }

      qa.qaP4("McTruth_", l, ntpPiNull);

      ntpPiNull->DumpData();


   }
  }


  //Write output
  out->cd();
  ntpMC -> GetInternalTree()->Write();
  ntpGamma ->GetInternalTree()->Write();
  ntpPiNull->GetInternalTree()->Write();

  out->Save();
  
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout<<endl<<endl;
  cout<<"Macro finisched successfully."<<endl;
  cout<<"Realtime: "<<rtime<<" s, CPU time: "<<ctime<<" s"<<endl;
  cout<<endl;

 
  exit(0);
 
}
