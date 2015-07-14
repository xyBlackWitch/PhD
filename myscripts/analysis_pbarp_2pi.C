class PndAnalysis;
class PndAnaPidSelector;
class RhoCandList;
class RhoTuple;


void analysis_pbarp_2pi(int nevts=0){
  
  TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);
  
  TStopwatch timer; 

  //Output File
  TString Path ="";// "/private/puetz/fairsoft_mar15/pandaroot/mysimulations/test/tracking/test_2pi_100_events_second_run/bugfix/";
  TString outPath ="";//"/private/puetz/fairsoft_mar15/pandaroot/mysimulations/test/tracking/test_2pi_100_events_second_run/bugfix/";
 
  TString OutputFile = outPath + "analysis_output.root";
  
  //Input simulation Files
  TString inPIDFile = Path + "pidideal_complete.root";
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
  RunAna->Init();

  //*** create tuples
  RhoTuple * ntpMC = new RhoTuple("ntpMC", "MCTruth info");
  RhoTuple * ntpPiMinus = new RhoTuple("ntpPiMinus", "PiMinus info");
  RhoTuple * ntpPiPlus = new RhoTuple("ntpPiPlus", "PiPlus info");
  RhoTuple * ntpBeam = new RhoTuple("ntpBeam", "Beam info");

  //Create output file for histograms
  TFile *out = TFile::Open(outPath+"output_ana.root","RECREATE");


  // data reader Object
  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts = theAnalysis->GetEntries();
  
  //RhoCandLists for analysis
  RhoCandList charged, piplus, piminus, beam;
  RhoCandList mclist, all;

  //***Mass selector
  double m0_beam = TDatabasePDG::Instance()->GetParticle("pbarpSystem")->Mass();
  cout<<"Mass pbar p system: "<<m0_beam<<endl;
  
  int evt = 0;
  //*** lorentz vector of the initial particle
  double pbarmom = 6.231552;
  double p_m0 = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
  TLorentzVector ini (0,0, pbarmom, sqrt(p_m0*p_m0+ pbarmom*pbarmom)+p_m0);
  
  
  PndRhoTupleQA qa(theAnalysis, pbarmom);
  //int nscattered = 0;

  while (theAnalysis->GetEvent() && evt++<nevts){
    if ((evt%100)==0) cout << "evt "<< evt <<endl;

    //***get MC list and store info
    theAnalysis->FillList(mclist, "McTruth");
    qa.qaMcList("", mclist, ntpMC);
    ntpMC->DumpData();

 
    
    

    //***Setup event shape object

	TString PidSelection = "";//"PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoDisc";	

    theAnalysis->FillList(all, "All", PidSelection);
    PndEventShape evsh(all, ini, 0.05, 0.1);
    
    //***Selection with no PID info
    theAnalysis->FillList(piminus, "PionAllMinus", PidSelection);
    theAnalysis->FillList(piplus, "PionAllPlus", PidSelection);

    //Momentum histograms
	
	for (int j=0; j<piminus.GetLength(); ++j){
			  
	  //info about about 4-vector
      qa.qaP4("PiMinus_",piminus[j]->P4(), ntpPiMinus);
      qa.qaCand("PiMinus_", piminus[j], ntpPiMinus);
      qa.qaVtx("PiMinus_", piminus[j], ntpPiMinus);
      qa.qaPoca("PiMinus_", piminus[j], ntpPiMinus);

	  RhoCandidate * mother = piminus[j]->GetMcTruth()->TheMother();
      int moth = (mother==0x0) ? 88888 : mother->PdgCode(); 			
			
			ntpPiMinus->Column("PiMinus_CosTheta", (Float_t) piminus[j]->GetMomentum().CosTheta());
			ntpPiMinus->Column("Mother", (Float_t) moth);
      ntpPiMinus->Column("MCTruthMatch", (bool) theAnalysis->McTruthMatch(piminus[j]));
			
			qa.qaP4("PiMinus_MC_", piminus[j]->GetMcTruth()->P4(), ntpPiMinus);
			ntpPiMinus->Column("PiMinus_MC_CosTheta", (Float_t) piminus[j]->GetMcTruth()->GetMomentum().CosTheta());
	
      ntpPiMinus->DumpData();


    }

    for (int j=0; j<piplus.GetLength(); ++j){
      qa.qaP4("PiPlus_", piplus[j]->P4(), ntpPiPlus);
      qa.qaCand("PiPlus_", piplus[j], ntpPiPlus);
      qa.qaVtx("PiPlus_", piplus[j], ntpPiPlus);
      qa.qaPoca("PiPlus_", piplus[j], ntpPiPlus);
	  
			RhoCandidate * mother = piplus[j]->GetMcTruth()->TheMother();
			int moth = (mother==0x0)? 88888 : mother->PdgCode() ;

 
			ntpPiPlus->Column("PiPlus_CosTheta", (Float_t) piplus[j]->GetMomentum().CosTheta());
			ntpPiPlus->Column("Mother", (Float_t) moth);
      ntpPiPlus->Column("MCTruthMatch", (bool) theAnalysis->McTruthMatch(piplus[j]));
			
			qa.qaP4("PiPlus_MC_", piplus[j]->GetMcTruth()->P4(), ntpPiPlus);
			ntpPiPlus->Column("PiPlus_MC_CosTheta", (Float_t) piplus[j]->GetMcTruth()->GetMomentum().CosTheta());


			//if(piplus[j]->GetMcTruth()->Pz()>-0.1 && piplus[j]->GetMcTruth()->Pz()<3 && piplus[j]->GetMcTruth()->Pt()<0.5 && theAnalysis->McTruthMatch(piplus[j])==1){
      //	nscattered+=1;
			/*RhoCandidate * mcmother = piplus[j]->GetMcTruth()->TheMother(); 
			int mother_pdg = (mcmother==0x0)? 88888 : mcmother->PdgCode() ;
			bool mctruthmatch = theAnalysis->McTruthMatch(piplus[j]);
			//if (mother_pdg==88888&&mctruthmatch){
				out << "Event: "<< evt <<" PiPlus Track "<< piplus[j]->GetTrackNumber()<<" (PDG:"<<piplus[j]->PdgCode() <<") has mother "<<mother_pdg;
				if (piplus[j]->NDaughters()>0) cout <<" and daughter(s) ";
				for (k=0;k<piplus[j]->NDaughters();++k) cout <<piplus[j]->Daughter(k)->GetTrackNumber()<<"  ";
				cout<<endl;
			//}*/
      ntpPiPlus->DumpData();


    }

    //***pbar p -> PiPlus + PiMinus
    beam.Combine(piminus,piplus);
    beam.SetType(113);
    
    for (int j=0; j<beam.GetLength(); ++j){
      
      //general info about event
      ntpBeam->Column("ev",     (Float_t) evt);
      ntpBeam->Column("cand",    (Float_t) j);
      ntpBeam->Column("ncand",   (Float_t) beam.GetLength());
      ntpBeam->Column("McTruth", (bool) theAnalysis->McTruthMatch(beam[j]));


      //beamPz->Fill(beam[j]->P());
      //qa.qaP4("Beam_", beam[j]->P4(), ntpBeam);
      qa.qaComp("Beam_", beam[j], ntpBeam);
      qa.qaEventShapeShort("es_", &evsh, ntpBeam);

      // do vertex fit
      PndKinVtxFitter vertexfitter (beam[j]);
      vertexfitter.Fit();
      RhoCandidate * beamFit = beam[j]->GetFit();

      // store info of vertex fit
      ntpBeam->Column("fvtx_Chi2", (float) vertexfitter.GetChi2());
      ntpBeam->Column("fvtx_NDF", (float) vertexfitter.GetNdf());
      ntpBeam->Column("fvtx_Prob", (float) vertexfitter.GetProb());  

      qa.qaMcDiff("fvtxMcDiff_", beamFit, ntpBeam);


      // do mass fit
      PndKinFitter massFitter(beamFit);
      massFitter.AddMassConstraint(m0_beam);
      massFitter.Fit();

      RhoCandidate * beamFit_mass = beamFit->GetFit();
      ntpBeam->Column("fMass_Chi2", (float) massFitter.GetChi2());
      ntpBeam->Column("fMass_NDF", (float) massFitter.GetNdf());
      ntpBeam->Column("fMass_Prob", (float) massFitter.GetProb()); 

      qa.qaMcDiff("Beam_", beamFit_mass, ntpBeam);

      RhoCandidate * truth = beam[j]->GetMcTruth();
      TLorentzVector l;
      if(0x0 != truth){
			l = truth->P4();
			qa.qaVtx("truth_", truth, ntpBeam);
      }
      qa.qaP4("truth_", l, ntpBeam);
      

      ntpBeam->DumpData();
    }
    
 
    
    //*** 4C Fit
    for (int j=0; j<beam.GetLength(); ++j){
      PndKinFitter fitter4c (beam[j]);
      fitter4c.Add4MomConstraint(ini);
      fitter4c.Fit();
      
      ntpBeam->Column("f4c_Chi2", (float) fitter4c.GetChi2());
      ntpBeam->Column("f4c_prob", (float) fitter4c.GetProb());
      ntpBeam->Column("f4c_NDF", (float) fitter4c.GetNdf());

      RhoCandidate * beam4cFit = beam[j]->GetFit();
      qa.qaMcDiff("f4c_", beam4cFit, ntpBeam);
      ntpBeam->DumpData();
      
    }
   }
	//cout<<"Scattered: "<<nscattered<<endl;

  //Write output
  out->cd();
  ntpMC -> GetInternalTree()->Write();
  ntpPiMinus ->GetInternalTree()->Write();
  ntpPiPlus->GetInternalTree()->Write();
  ntpBeam->GetInternalTree()->Write();

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
