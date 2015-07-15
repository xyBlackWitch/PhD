class PndAnalysis;
class PndAnaPidSelector;
class RhoCandList;
class RhoTuple;

#include "common_jenny.cpp"


enum pidNumbers {
	kPip = 211, kPim = -211,
	kPp = 2212, kaPm = -2212,
	kl0 = 3122, kal0 = -3122,
	kXim = 3312, kaXip = -3312
};

void analysis_pbarp_Xi(int nevts=0){
  
  TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);
  
  TStopwatch timer; 



  //Output File
  TString Path = "/private/puetz/mysimulations/test/boxgenerator/Xi/run2/old/";
  TString outPath = Path;
  TString OutputFile = outPath + "analysis_output_oldidealtracking.root";
  
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
  RunAna->Init();

  //*** create tuples
  RhoTuple * ntpMC = new RhoTuple("ntpMC", "MCTruth info");
  RhoTuple * ntpPiMinus = new RhoTuple("ntpPiMinus", "PiMinus info");
  RhoTuple * ntpPiPlus = new RhoTuple("ntpPiPlus", "PiPlus info");
  RhoTuple * ntpProton = new RhoTuple("ntpProton", "Proton info");
  RhoTuple * ntpAntiProton = new RhoTuple("ntpAntiProton", "Antiproton info");
  RhoTuple * ntpLambda0 = new RhoTuple("ntpLambda0", "Lambda0 info");
  RhoTuple * ntpAntiLambda0 = new RhoTuple("ntpAntiLambda0", "AntiLambda0 info");
  RhoTuple * ntpXiMinus = new RhoTuple("ntpXiMinus", "XiMinus info");
  RhoTuple * ntpXiPlus = new RhoTuple("ntpXiPlus", "XiPlus info");
  RhoTuple * ntpXiSys = new RhoTuple("ntpXiSys", "XiMinus XiPlus system info");

  //Create output file 
  TFile *out = TFile::Open(outPath+"output_ana_oldidealtracking.root","RECREATE");

  // data reader Object
  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts = theAnalysis->GetEntries();
  
  //RhoCandLists for analysis
  RhoCandList piplus, piminus, lambda0, antiLambda0, proton, antiProton, xiplus, ximinus, xiSys;
//  RhoCandList BestCandLambda0, BestCandAntiLambda0, NotCombinedPiMinus, CombinedPiMinus, CombinedPiPlus, NotCombinedPiPlus;
  RhoCandList Lambda0Fit, AntiLambda0Fit, XiMinusFit, XiPlusFit;
  RhoCandList mclist, all;

  //***Mass selector
  double m0_lambda0= TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
  cout<<"Mass of Lambda0: "<<m0_lambda0<<endl;
  RhoMassParticleSelector * lambdaMassSelector = new RhoMassParticleSelector("lambda0", m0_lambda0, 0.3);
 
  double m0_Xi = TDatabasePDG::Instance()->GetParticle("Xi-")->Mass();
  cout<<"Mass of Xi-: "<<m0_Xi<<endl;
  RhoMassParticleSelector * xiMassSelector = new RhoMassParticleSelector("Xi-", m0_Xi, 0.3);

  double m0_pbarpsystem = TDatabasePDG::Instance()->GetParticle("pbarpSystem")->Mass();

  
  double pbarmom = 2.7;
  double p_m0 = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
  TLorentzVector ini (0,0, pbarmom, sqrt(p_m0*p_m0+ pbarmom*pbarmom)+p_m0);
  TVector3 beamBoost = ini.BoostVector();
  
  PndRhoTupleQA qa(theAnalysis, pbarmom);
 
  int evt=-1;
  int index=0;

  while (theAnalysis->GetEvent() && ++evt<nevts){

    if ((evt%100)==0) cout << "evt "<< evt <<endl;
//    cout << "Running event " << evt << endl;
		
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
		
    TString PidSelection = "PidAlgoIdealCharged";//"PidAlgoMvd;PidAlgoStt;PidAlgoDrc";

    theAnalysis->FillList(all, "All", PidSelection);
    PndEventShape evsh(all, ini, 0.05, 0.1);
    
    //***Selection with no PID info
    theAnalysis->FillList(piminus, "PionAllMinus", PidSelection);
//    theAnalysis->FillList(NotCombinedPiMinus, "PionAllMinus", PidSelection);
    theAnalysis->FillList(piplus, "PionAllPlus", PidSelection);
    theAnalysis->FillList(proton, "ProtonAllPlus", PidSelection);
    theAnalysis->FillList(antiProton, "ProtonAllMinus", PidSelection);

//
//    for (int pip=0; pip<piplus.GetLength(); ++pip){
//        ntpPiPlus->Column("ev",     (Float_t) evt);
//        ntpPiPlus->Column("cand",    (Float_t) pip);
//        ntpPiPlus->Column("ncand",   (Float_t) piplus.GetLength());
//        ntpPiPlus->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(piplus[pip]));
//
//        qa.qaP4("PiPlus_", piplus[pip]->P4(), ntpPiPlus);
//
//        RhoCandidate * mother_pip = piplus[pip]->GetMcTruth()->TheMother();
//        int moth_pip = (0x0==mother_pip)? 88888 : mother_pip->PdgCode();
//
//        ntpPiPlus->Column("Mother", (Float_t) moth_pip);
//        ntpPiPlus->Column("PiPlus_CosTheta", (Float_t) piplus[pip]->GetMomentum().CosTheta());
//        ntpPiPlus->Column("PiPlus_Pdg", (Float_t) piplus[pip]->PdgCode());
//
//        qa.qaP4("PiPlus_MC_", piplus[pip]->GetMcTruth()->P4(), ntpPiPlus);
//        ntpPiPlus->Column("PiPlus_MC_CosTheta", (Float_t) piplus[pip]->GetMcTruth()->GetMomentum().CosTheta());
//
//        ntpPiPlus->DumpData();
//    }
//
//    for (int pim=0; pim<piminus.GetLength(); ++pim){
//        ntpPiMinus->Column("ev",     (Float_t) evt);
//        ntpPiMinus->Column("cand",    (Float_t) pim);
//        ntpPiMinus->Column("ncand",   (Float_t) piminus.GetLength());
//        ntpPiMinus->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(piminus[pim]));
//
//        qa.qaP4("piminus_", piminus[pim]->P4(), ntpPiMinus);
//
//        RhoCandidate * mother_pim = piminus[pim]->GetMcTruth()->TheMother();
//        int moth_pim = (0x0==mother_pim)? 88888 : mother_pim->PdgCode();
//
//        ntpPiMinus->Column("Mother", (Float_t) moth_pim);
//        ntpPiMinus->Column("PiMinus_CosTheta", (Float_t) piminus[pim]->GetMomentum().CosTheta());
//        ntpPiMinus->Column("PiMinus_Pdg", (Float_t) piminus[pim]->PdgCode());
//
//        qa.qaP4("piminus_MC_", piminus[pim]->GetMcTruth()->P4(), ntpPiMinus);
//        ntpPiMinus->Column("piminus_MC_CosTheta", (Float_t) piminus[pim]->GetMcTruth()->GetMomentum().CosTheta());
//
//        ntpPiMinus->DumpData();
//    }
//
//    for (int prot=0; prot<proton.GetLength(); ++prot){
//        ntpProton->Column("ev",     (Float_t) evt);
//        ntpProton->Column("cand",    (Float_t) prot);
//        ntpProton->Column("ncand",   (Float_t) proton.GetLength());
//        ntpProton->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(proton[prot]));
//
//        qa.qaP4("proton_", proton[prot]->P4(), ntpProton);
//
//        RhoCandidate * mother_prot = proton[prot]->GetMcTruth()->TheMother();
//        int moth_prot = (0x0==mother_prot)? 88888 : mother_prot->PdgCode();
//
//        ntpProton->Column("Mother", (Float_t) moth_prot);
//        ntpProton->Column("proton_CosTheta", (Float_t) proton[prot]->GetMomentum().CosTheta());
//        ntpProton->Column("proton_Pdg", (Float_t) proton[prot]->PdgCode());
//
//        qa.qaP4("proton_MC_", proton[prot]->GetMcTruth()->P4(), ntpProton);
//        ntpProton->Column("proton_MC_CosTheta", (Float_t) proton[prot]->GetMcTruth()->GetMomentum().CosTheta());
//
//        ntpProton->DumpData();
//    }
//
//    for (int aProt=0; aProt<antiProton.GetLength(); ++aProt){
//        ntpAntiProton->Column("ev",     (Float_t) evt);
//        ntpAntiProton->Column("cand",    (Float_t) aProt);
//        ntpAntiProton->Column("ncand",   (Float_t) antiProton.GetLength());
//        ntpAntiProton->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(antiProton[aProt]));
//
//        qa.qaP4("antiProton_", antiProton[aProt]->P4(), ntpAntiProton);
//
//        RhoCandidate * mother_aProt = antiProton[aProt]->GetMcTruth()->TheMother();
//        int moth_aProt = (0x0==mother_aProt)? 88888 : mother_aProt->PdgCode();
//
//        ntpAntiProton->Column("Mother", (Float_t) moth_aProt);
//        ntpAntiProton->Column("antiProton_CosTheta", (Float_t) antiProton[aProt]->GetMomentum().CosTheta());
//        ntpAntiProton->Column("antiProton_Pdg", (Float_t) antiProton[aProt]->PdgCode());
//
//        qa.qaP4("antiProton_MC_", antiProton[aProt]->GetMcTruth()->P4(), ntpAntiProton);
//        ntpAntiProton->Column("antiProton_MC_CosTheta", (Float_t) antiProton[aProt]->GetMcTruth()->GetMomentum().CosTheta());
//
//        ntpAntiProton->DumpData();
//    }





    //***Lambda0 -> PiMinus + Proton

    lambda0.Combine(piminus,proton);
	lambda0.Select(lambdaMassSelector);
    lambda0.SetType(kl0);
//
//    std::map<int,int> bestVtxFitLambda0, bestMassFitLambda0;
//
//    bestVtxFitLambda0 = jenny::VertexQaIndex(&lambda0);
//    bestMassFitLambda0 = jenny::MassFitQaIndex(&lambda0, m0_lambda0);


    for (int j=0; j<lambda0.GetLength(); ++j){


      //general info about event
      ntpLambda0->Column("ev",     (Float_t) evt);
      ntpLambda0->Column("cand",    (Float_t) j);
      ntpLambda0->Column("ncand",   (Float_t) lambda0.GetLength());
      ntpLambda0->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(lambda0[j]));
      ntpLambda0->Column("Lambda0_Pdg", (Float_t) lambda0[j]->PdgCode());

	  RhoCandidate * mother = lambda0[j]->TheMother();
  	  int moth = (mother==0x0) ? 88888 : mother->PdgCode();

	  ntpLambda0->Column("Mother", (Float_t) moth);

      qa.qaP4("Lambda0_", lambda0[j]->P4(), ntpLambda0);
      qa.qaComp("Lambda0_", lambda0[j], ntpLambda0);



      // do vertex fit
      PndKinVtxFitter vertexfitterLambda0 (lambda0[j]);
	  vertexfitterLambda0.Fit();
      RhoCandidate * lambda0Fit = lambda0[j]->GetFit();


      // store info of vertex fit
      qa.qaFitter("VtxFit_", &vertexfitterLambda0, ntpLambda0);
//      ntpLambda0->Column("VtxFit_HowGood", (Int_t) bestVtxFitLambda0[j]);
      qa.qaVtx("VtxFit_", lambda0Fit, ntpLambda0);

      // differenz to MCTruth
       qa.qaMcDiff("fvtxMcDiff_", lambda0Fit, ntpLambda0);


//      jenny::CombinedList(lambda0[j], &CombinedPiMinus, -211);


//      for (int daughter=0; daughter<lambda0[j]->NDaughters(); daughter++){
//    	  RhoCandidate * daughterCand = lambda0[j]->Daughter(daughter);
//    	  if (daughterCand->PdgCode()==-211) CombinedPiMinus.Append(daughterCand);
//      }
//
//      CombinedPiMinus.RemoveClones();


      // do mass fit
      PndKinFitter massFitterLambda0(lambda0Fit);
      massFitterLambda0.AddMassConstraint(m0_lambda0);
      massFitterLambda0.Fit();

      RhoCandidate * lambda0Fit_mass = lambda0Fit->GetFit();
      ntpLambda0->Column("fMass_Chi2", (float) massFitterLambda0.GetChi2());
      ntpLambda0->Column("fMass_NDF", (float) massFitterLambda0.GetNdf());
      ntpLambda0->Column("fMass_Prob", (float) massFitterLambda0.GetProb());

      qa.qaMcDiff("Lambda0_", lambda0Fit_mass, ntpLambda0);

      Lambda0Fit.Append(lambda0Fit);

      RhoCandidate * truth = lambda0[j]->GetMcTruth();
      TLorentzVector l;

      if(0x0 != truth){
				l = truth->P4();
				qa.qaVtx("MCTruth_", truth, ntpLambda0);
      }
      qa.qaP4("MCTruth_", l, ntpLambda0);


//			//***information of boosted particle
//			lambda0Fit->Boost(-beamBoost);
//			qa.qaComp("boost_", lambda0Fit, ntpLambda0);

		  ntpLambda0->DumpData();

   }

//    for (int pim=0; pim<CombinedPiMinus.GetLength(); pim++){
//        	  NotCombinedPiMinus.Remove(CombinedPiMinus[pim]);
//          }
//

//     CombinedPiMinus.Cleanup();

//    //***AntiLambda0 -> PiPlus + AntiProton
//    antiLambda0.Combine(piplus,antiProton);
//	antiLambda0.Select(lambdaMassSelector);
//    antiLambda0.SetType(kal0);
//
//    for (int j=0; j<antiLambda0.GetLength(); ++j){
//
//      //general info about event
//      ntpAntiLambda0->Column("ev",     (Float_t) evt);
//      ntpAntiLambda0->Column("cand",    (Float_t) j);
//      ntpAntiLambda0->Column("ncand",   (Float_t) antiLambda0.GetLength());
//      ntpAntiLambda0->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(antiLambda0[j]));
//      ntpAntiLambda0->Column("AntiLambda0_Pdg", (Float_t) antiLambda0[j]->PdgCode());
//
//	  RhoCandidate * mother = antiLambda0[j]->TheMother();
//	  int moth = (mother==0x0) ? 88888 : mother->PdgCode();
//	  ntpAntiLambda0->Column("Mother", (Float_t) moth);
//
//      qa.qaP4("AntiLambda0_", antiLambda0[j]->P4(), ntpAntiLambda0);
//      qa.qaComp("AntiLambda0_", antiLambda0[j], ntpAntiLambda0);
//      //qa.qaEventShapeShort("es_", &evsh, ntpAntiLambda0);
//
//
//
//      // do vertex fit
//      PndKinVtxFitter vertexfitterAntiLambda0 (antiLambda0[j]);
//      vertexfitterAntiLambda0.Fit();
//      RhoCandidate * antiLambda0Fit = antiLambda0[j]->GetFit();
//
//      AntiLambda0Fit.Append(antiLambda0Fit);
//
//      // store info of vertex fit
//      ntpAntiLambda0->Column("fvtx_Chi2", (float) vertexfitterAntiLambda0.GetChi2());
//      ntpAntiLambda0->Column("fvtx_NDF", (float) vertexfitterAntiLambda0.GetNdf());
//      ntpAntiLambda0->Column("fvtx_Prob", (float) vertexfitterAntiLambda0.GetProb());
//
//
//			// difference to MCTruth
//      qa.qaMcDiff("fvtxMcDiff_", antiLambda0Fit, ntpAntiLambda0);
//
//
//      float chi2 = vertexfitterAntiLambda0.GetChi2();
//      float prob = vertexfitterAntiLambda0.GetProb();
//
////      if (chi2 < bestchi && prob>0.01){
////		  bestchi = chi2;
////		  bestchicand=j;
////		  BestCandAntiLambda0.Add(antiLambda0Fit);//InsertAt(0, antiLambda0Fit);
////      }
////
////      for (int daughter=0; daughter<antiLambda0[j]->NDaughters(); daughter++){
////      	  RhoCandidate * daughterCand = antiLambda0[j]->Daughter(daughter);
////      	  if (daughterCand->PdgCode()==-211) CombinedPiPlus.Append(daughterCand);
////       }
////
////       CombinedPiPlus.RemoveClones();
//
//
//
//      // do mass fit
//      PndKinFitter massFitterAntiLambda0(antiLambda0Fit);
//      massFitterAntiLambda0.AddMassConstraint(m0_lambda0);
//      massFitterAntiLambda0.Fit();
//
//      RhoCandidate * antiLambda0Fit_mass = antiLambda0Fit->GetFit();
//      ntpAntiLambda0->Column("fMass_Chi2", (float) massFitterAntiLambda0.GetChi2());
//      ntpAntiLambda0->Column("fMass_NDF", (float) massFitterAntiLambda0.GetNdf());
//      ntpAntiLambda0->Column("fMass_Prob", (float) massFitterAntiLambda0.GetProb());
//
//      qa.qaMcDiff("AntiLambda0_", antiLambda0Fit_mass, ntpAntiLambda0);
//
//      RhoCandidate * truth = antiLambda0[j]->GetMcTruth();
//
//      TLorentzVector l;
//      if(0x0 != truth){
//				l = truth->P4();
//				qa.qaVtx("MCTruth_", truth, ntpAntiLambda0);
//      }
//      qa.qaP4("MCTruth_", l, ntpAntiLambda0);
//
////			//***information of boosted particle
////			antiLambda0Fit->Boost(-beamBoost);
////	 		qa.qaComp("boost_", antiLambda0Fit, ntpAntiLambda0);
//
//      ntpAntiLambda0->DumpData();
//    }
//
//
//

    //*** Xi- -> Lambda0 + Pi-
	ximinus.Combine(Lambda0Fit, piminus);//BestCandLambda0,NotCombinedPiMinus);
	ximinus.Select(xiMassSelector);
	ximinus.SetType(kXim);

    for (int j=0; j<ximinus.GetLength(); ++j){

      //general info about event
      ntpXiMinus->Column("ev",     (Float_t) evt);
      ntpXiMinus->Column("cand",    (Float_t) j);
      ntpXiMinus->Column("ncand",   (Float_t) ximinus.GetLength());
      ntpXiMinus->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(ximinus[j]));

	  RhoCandidate * mother = ximinus[j]->TheMother();
	  int moth = (mother==0x0) ? 88888 : mother->PdgCode();
	  ntpXiMinus->Column("Mother", (Float_t) moth);

      qa.qaP4("XiMinus_", ximinus[j]->P4(), ntpXiMinus);
      qa.qaComp("XiMinus_", ximinus[j], ntpXiMinus);
      qa.qaPoca("XiMinus_", ximinus[j], ntpXiMinus);
      //qa.qaEventShapeShort("es_", &evsh, ntpXiMinus);



      // do vertex-fit



      PndKinVtxFitter vertexfitterXiMinus (ximinus[j]);
      vertexfitterXiMinus.Fit();
      RhoCandidate * ximinusFit = ximinus[j]->GetFit();

      XiMinusFit.Append(ximinusFit);



      // store info of vertex-fit
      ntpXiMinus->Column("FitVertex_Chi2", (float) vertexfitterXiMinus.GetChi2());
      ntpXiMinus->Column("FitVertex_NDF", (float) vertexfitterXiMinus.GetNdf());
      ntpXiMinus->Column("FitVertex_Prob", (float) vertexfitterXiMinus.GetProb());

//      if (vertexfitterXiMinus.GetChi2()<-1e-5){
//          		  cout << "Chi2 = " << vertexfitterXiMinus.GetChi2() <<" Probability = " << vertexfitterXiMinus.GetProb() << endl;
//            	  	  counter+=1;
//      }

      //qa.qaFitter("Fitvtx_", vertexfitterXiMinus, ntpXiMinus);
      qa.qaVtx("XiMinusFit_", ximinusFit, ntpXiMinus);
      qa.qaCand("XiMinusFit_", ximinusFit, ntpXiMinus);

	  // difference to MCTruth
      qa.qaMcDiff("fvtxMcDiff_", ximinusFit, ntpXiMinus);


      // do mass fit
      PndKinFitter massFitterXiMinus(ximinusFit);
      massFitterXiMinus.AddMassConstraint(m0_lambda0);
      massFitterXiMinus.Fit();

      RhoCandidate * ximinusFit_mass = ximinusFit->GetFit();
      ntpXiMinus->Column("fMass_Chi2", (float) massFitterXiMinus.GetChi2());
      ntpXiMinus->Column("fMass_NDF", (float) massFitterXiMinus.GetNdf());
      ntpXiMinus->Column("fMass_Prob", (float) massFitterXiMinus.GetProb());

      qa.qaMcDiff("MasFit_", ximinusFit_mass, ntpXiMinus);



      RhoCandidate * truth = ximinus[j]->GetMcTruth();
      TLorentzVector l;
      TVector3 decayvtx;
      if(0x0 != truth){
				l = truth->P4();
      }

      RhoCandidate * daughterTruth = ximinus[j]->Daughter(0)->GetMcTruth();
      TVector3 dl;

      if (0x0 != daughterTruth) dl = daughterTruth->Pos();

      jenny::qaP3("MCTruth_", dl, ntpXiMinus);
      qa.qaP4("MCTruth_", l, ntpXiMinus);

//			//***information of boosted particle
//			ximinusFit->Boost(-beamBoost);
//	 		qa.qaComp("boost_", ximinusFit, ntpXiMinus);
//
      ntpXiMinus->DumpData();
	}
    Lambda0Fit.Cleanup();
//    BestCandLambda0.Cleanup();
//    NotCombinedPiMinus.Cleanup();
//
//
//
//	//*** Xi+ -> AntiLambda0 + Pi+
//		xiplus.Combine(AntiLambda0Fit,piplus);
//		xiplus.Select(xiMassSelector);
//		xiplus.SetType(kaXip);
//
//    for (int j=0; j<xiplus.GetLength(); ++j){
//
//      //general info about event
//      ntpXiPlus->Column("ev",     (Float_t) evt);
//      ntpXiPlus->Column("cand",    (Float_t) j);
//      ntpXiPlus->Column("ncand",   (Float_t) xiplus.GetLength());
//      ntpXiPlus->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(xiplus[j]));
//
//      RhoCandidate * mother = xiplus[j]->TheMother();
//      int moth = (mother==0x0) ? 88888 : mother->PdgCode();
//      ntpXiPlus->Column("Mother", (Float_t) moth);
//
//      qa.qaP4("xiplus_", xiplus[j]->P4(), ntpXiPlus);
//      //qa.qaComp("xiplus_", xiplus[j], ntpXiPlus);
//
//
//
//
//      // do vertex-fit
//      PndKinVtxFitter vertexfitterxiplus (xiplus[j]);
//      vertexfitterxiplus.Fit();
//      RhoCandidate * xiplusFit = xiplus[j]->GetFit();
//
//      XiPlusFit.Append(xiplusFit);
//
//      // store info of vertex-fit
//      ntpXiPlus->Column("FitVertex_Chi2", (float) vertexfitterxiplus.GetChi2());
//      ntpXiPlus->Column("FitVertex_NDF", (float) vertexfitterxiplus.GetNdf());
//      ntpXiPlus->Column("FitVertex_Prob", (float) vertexfitterxiplus.GetProb());
//
//
//	  // difference to MCTruth
//      qa.qaMcDiff("fvtxMcDiff_", xiplusFit, ntpXiPlus);
//
//
//      // do mass fit
//      PndKinFitter massFitterxiplus(xiplusFit);
//      massFitterxiplus.AddMassConstraint(m0_lambda0);
//      massFitterxiplus.Fit();
//
//      RhoCandidate * xiplusFit_mass = xiplusFit->GetFit();
//      ntpXiPlus->Column("fMass_Chi2", (float) massFitterxiplus.GetChi2());
//      ntpXiPlus->Column("fMass_NDF", (float) massFitterxiplus.GetNdf());
//      ntpXiPlus->Column("fMass_Prob", (float) massFitterxiplus.GetProb());
//
//      qa.qaMcDiff("xiplus_", xiplusFit_mass, ntpXiPlus);
//
//      RhoCandidate * truth = xiplus[j]->GetMcTruth();
//      TLorentzVector l;
//      if(0x0 != truth){
//				l = truth->P4();
//				qa.qaVtx("MCTruth_", truth, ntpXiPlus);
//      }
//      qa.qaP4("MCTruth_", l, ntpXiPlus);
//
////			//***information of boosted particle
////			xiplusFit->Boost(-beamBoost);
////	 		qa.qaComp("boost_", xiplusFit, ntpXiPlus);
//
//      ntpXiPlus->DumpData();
//	 }
//    AntiLambda0Fit.Cleanup();
////    BestCandAntiLambda0.Cleanup();

//
//
//    //******* Xi+ Xi- System*****************************
//
//    xiSys.Combine(XiPlusFit, XiMinusFit);
//    xiSys.SetType(88888);
//
//    for (int syscand=0; syscand<xiSys.GetLength(); ++syscand){
//
//        ntpXiSys->Column("ev",     (Float_t) evt);
//        ntpXiSys->Column("cand",    (Float_t) j);
//        ntpXiSys->Column("ncand",   (Float_t) ximinus.GetLength());
//        ntpXiSys->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(xiSys[syscand]));
//
//  	  RhoCandidate * mother = xiSys[syscand]->TheMother();
//  	  int moth = (mother==0x0) ? 88888 : mother->PdgCode();
//  	  ntpXiSys->Column("Mother", (Float_t) moth);
//
//        qa.qaP4("XiSys__", xiSys[syscand]->P4(), ntpXiSys);
//        qa.qaComp("XiSys__", xiSys[syscand], ntpXiSys);
//        qa.qaPoca("XiSys__", xiSys[syscand], ntpXiSys);
//
//        RhoCandidate *  truth = xiSys[syscand]->GetMcTruth();
//        if(truth != 0x0) qa.qaComp("MC_", truth, ntpXiSys);
//
//        //4C-Fitter
//
//        PndKinFitter fitter4c (xiSys[syscand]);
//        fitter4c.Add4MomConstraint(ini);
//        fitter4c.Fit();
//
//        RhoCandidate * xiSysFit4c = xiSys[syscand]->GetFit();
//
//        ntpXiSys->Column("Fit4c_Chi2", (float) fitter4c.GetChi2());
//        ntpXiSys->Column("Fit4c_Prob", (float) fitter4c.GetProb());
//        ntpXiSys->Column("Fit4c_NDF", (float) fitter4c.GetNdf());
//
//
//        qa.qaComp("XiSys4cFit_", xiSysFit4c, ntpXiSys);
//
//
//
//        ntpXiSys->DumpData();
//
//
//    }
    XiMinusFit.Cleanup();
    XiPlusFit.Cleanup();
  }



   

  //Write output
  out->cd();

  ntpMC -> GetInternalTree()->Write();
//  ntpPiMinus ->GetInternalTree()->Write();
//  ntpPiPlus->GetInternalTree()->Write();
//  ntpProton->GetInternalTree()->Write();
//  ntpAntiProton->GetInternalTree()->Write();
  ntpLambda0->GetInternalTree()->Write();
//  ntpAntiLambda0->GetInternalTree()->Write();
  ntpXiMinus->GetInternalTree()->Write();
//  ntpXiPlus->GetInternalTree()->Write();
//  ntpXiSys->GetInternalTree()->Write();

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
