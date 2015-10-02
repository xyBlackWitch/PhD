/**
* @file analysis_pbarp_xi_1820.C
* @mainpage analysis_pbarp_xi_1820.C Analysis macro for the reaction pbar p -> Xi+ Xi(1820)-
*
* @author Jennifer Puetz (jennifer.puetz@fz-juelich.de)
* @date 2015
* @brief analysis macro
* @details This file holds the analysis of the reaction
* reaction pbar p -> Xi+ Xi(1820)-
* 					  |   |
*  					  |   -> Lambda0 + K-
*  					  |			|
*  					  |			-> p + Pi-	
* 					   -> AntiLambda0 + Pi+
* 					   		|
* 					   		-> pbar + Pi+
*
*/

class PndAnalysis;
class PndAnaPidSelector;
class RhoCandList;
class RhoTuple;

#include "../common_jenny.cpp"


enum pidNumbers {
	kPip = 211, kPim = -211,
	kPp = 2212, kaPm = -2212,
	kl0 = 3122, kal0 = -3122,
	kKm = -321,
	kXim = 23314, kaXip = -3312
};



void analysis_pbarp_Xi_1820(int nevts=0, double mom=4.6){
  
  TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);
  
  TStopwatch timer; 


  //Output File
  TString Path = "/home/ikp1/puetz/panda/myscripts/simChain/";
  TString outPath = Path;
  TString OutputFile = outPath + "analysis_output.root";
  
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
  RhoTuple * ntpKaonMinus = new RhoTuple("ntpKaonMinus", "KaonMinus info");
  RhoTuple * ntpLambda0 = new RhoTuple("ntpLambda0", "Lambda0 info");
  RhoTuple * ntpAntiLambda0 = new RhoTuple("ntpAntiLambda0", "AntiLambda0 info");
  RhoTuple * ntpXiMinus = new RhoTuple("ntpXi1820Minus", "XiMinus info");
  RhoTuple * ntpXiPlus = new RhoTuple("ntpXiPlus", "XiPlus info");
  RhoTuple * ntpXiSys = new RhoTuple("ntpXiSys", "XiMinus XiPlus system info");

  //Create output file 
  TFile *out = TFile::Open(outPath+"output_ana.root","RECREATE");

  // data reader Object
  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts = theAnalysis->GetEntries();
  
  //RhoCandLists for analysis
  RhoCandList piplus, piminus, lambda0, antiLambda0, kaonMinus, proton, antiProton, xiplus, ximinus, xiSys;
  RhoCandList CombinedPiPlus, NotCombinedPiPlus;
  RhoCandList Lambda0Fit, AntiLambda0Fit, XiMinusFit, XiPlusFit;
  RhoCandList mclist, all;

  //Dummy RhoCandidate
  RhoCandidate * dummyCand = new RhoCandidate();


  //***Mass selector
  double m0_lambda0= TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
  cout<<"Mass of Lambda0: "<<m0_lambda0<<endl;
  RhoMassParticleSelector * lambdaMassSelector = new RhoMassParticleSelector("lambda0", m0_lambda0, 0.3);

  double m0_Xi = TDatabasePDG::Instance()->GetParticle("Xi-")->Mass();
  cout<<"Mass of Xi-: "<<m0_Xi<<endl;
  RhoMassParticleSelector * xiMassSelector = new RhoMassParticleSelector("Xi-", m0_Xi, 0.3);
 
  double m0_Xi1820 = 1.823;//TDatabasePDG::Instance()->GetParticle("Xi(1820)-")->Mass();
  cout<<"Mass of Xi(1820)-: "<<m0_Xi1820<<endl;
  RhoMassParticleSelector * xi1820MassSelector = new RhoMassParticleSelector("Xi(1820)-", m0_Xi1820, 0.3);

  double m0_pbarpsystem = TDatabasePDG::Instance()->GetParticle("pbarpSystem")->Mass();

  
  double pbarmom = mom;
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
    theAnalysis->FillList(piminus, "PionBestMinus", PidSelection);
    theAnalysis->FillList(NotCombinedPiPlus, "PionBestPlus", PidSelection);
    theAnalysis->FillList(piplus, "PionBestPlus", PidSelection);
    theAnalysis->FillList(proton, "ProtonBestPlus", PidSelection);
    theAnalysis->FillList(antiProton, "ProtonBestMinus", PidSelection);
    theAnalysis->FillList(kaonMinus, "KaonBestMinus", PidSelection);



	
    for (int pip=0; pip<piplus.GetLength(); ++pip){
    	
        ntpPiPlus->Column("ev",     (Float_t) evt);
        ntpPiPlus->Column("cand",    (Float_t) pip);
        ntpPiPlus->Column("ncand",   (Float_t) piplus.GetLength());
        ntpPiPlus->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(piplus[pip]));


        qa.qaP4("PiPlus_", piplus[pip]->P4(), ntpPiPlus);
        qa.qaCand("PiPlus_", piplus[pip], ntpPiPlus);

        jenny::numberOfHitsInSubdetector("PiPlus_", piplus[pip], ntpPiPlus);
        jenny::tagNHits("PiPlus_", piplus[pip], ntpPiPlus);

        RhoCandidate * truth = piplus[pip]->GetMcTruth();
        RhoCandidate * mother=0;


        TLorentzVector l;
        float costht = -999.;
        
        if(truth!=0x0){
        	l=truth->P4();
        	mother = truth->TheMother();
        	costht= truth->GetMomentum().CosTheta();
        }
                
        int moth_pip = (0x0==mother)? 88888 : mother->PdgCode();


        ntpPiPlus->Column("Mother", (Float_t) moth_pip);
        ntpPiPlus->Column("PiPlus_CosTheta", (Float_t) piplus[pip]->GetMomentum().CosTheta());

        qa.qaP4("PiPlus_MC_", l, ntpPiPlus);
        qa.qaCand("PiPlus_MC_", truth, ntpPiPlus);
        ntpPiPlus->Column("PiPlus_MC_CosTheta", (Float_t) costht);

        
        ntpPiPlus->DumpData();
    }

	
    for (int pim=0; pim<piminus.GetLength(); ++pim){
    	
        ntpPiMinus->Column("ev",     (Float_t) evt);
        ntpPiMinus->Column("cand",    (Float_t) pim);
        ntpPiMinus->Column("ncand",   (Float_t) piminus.GetLength());
        ntpPiMinus->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(piminus[pim]));

        qa.qaP4("piminus_", piminus[pim]->P4(), ntpPiMinus);
        qa.qaCand("piminus_", piminus[pim], ntpPiMinus);

        jenny::numberOfHitsInSubdetector("piminus_", piminus[pim], ntpPiMinus);
        jenny::tagNHits("piminus_", piminus[pim], ntpPiMinus);

        RhoCandidate * truth = piminus[pim]->GetMcTruth();
        RhoCandidate * mother_pim=0;
        TLorentzVector l;
        float costht = -999.;
        
        if(truth!=0x0){
        	l=truth->P4();
        	costht= truth->GetMomentum().CosTheta();
        	mother_pim=truth->TheMother();
        }
                
        int moth_pim = (mother_pim==0x0) ? 88888 : mother_pim->PdgCode();

        ntpPiMinus->Column("Mother", (Float_t) moth_pim);
        ntpPiMinus->Column("piminus_CosTheta", (Float_t) piminus[pim]->GetMomentum().CosTheta());

        qa.qaP4("piminus_MC_", l, ntpPiMinus);
        qa.qaCand("piminus_MC_", truth, ntpPiMinus);
        ntpPiMinus->Column("piminus_MC_CosTheta", (Float_t) costht);
        
        ntpPiMinus->DumpData();
    }
	
    for (int prot=0; prot<proton.GetLength(); ++prot){
        ntpProton->Column("ev",     (Float_t) evt);
        ntpProton->Column("cand",    (Float_t) prot);
        ntpProton->Column("ncand",   (Float_t) proton.GetLength());
        ntpProton->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(proton[prot]));

        qa.qaP4("proton_", proton[prot]->P4(), ntpProton);
        qa.qaCand("proton_", proton[prot], ntpProton);

        jenny::numberOfHitsInSubdetector("proton_", proton[prot], ntpProton);
        jenny::tagNHits("proton_", proton[prot], ntpProton);

		RhoCandidate * truth = proton[prot]->GetMcTruth();
		RhoCandidate * mother_prot=0;
		TLorentzVector l;
		float costht = -999.;
		
		if(truth!=0x0){
			l=truth->P4();
			costht= truth->GetMomentum().CosTheta();
			mother_prot= truth->TheMother();
		}
			 
		int moth_prot = (0x0==mother_prot)? 88888 : mother_prot->PdgCode();
		
		ntpProton->Column("Mother", (Float_t) moth_prot);
		ntpProton->Column("proton_CosTheta", (Float_t) proton[prot]->GetMomentum().CosTheta());
		
		qa.qaP4("proton_MC_", l, ntpProton);
		qa.qaCand("proton_MC_", truth, ntpProton);
		ntpProton->Column("proton_MC_CosTheta", (Float_t) costht);
         

        ntpProton->DumpData();
    }
	
    for (int aProt=0; aProt<antiProton.GetLength(); ++aProt){
        ntpAntiProton->Column("ev",     (Float_t) evt);
        ntpAntiProton->Column("cand",    (Float_t) aProt);
        ntpAntiProton->Column("ncand",   (Float_t) antiProton.GetLength());
        ntpAntiProton->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(antiProton[aProt]));

        qa.qaP4("antiProton_", antiProton[aProt]->P4(), ntpAntiProton);
        qa.qaCand("antiProton_", antiProton[aProt], ntpAntiProton);

        jenny::numberOfHitsInSubdetector("antiProton_", antiProton[aProt], ntpAntiProton);
        jenny::tagNHits("antiProton_", antiProton[aProt], ntpAntiProton);

		RhoCandidate * truth = antiProton[aProt]->GetMcTruth();
		RhoCandidate * mother_aProt=0;
		TLorentzVector l;
		float costht = -999.;
		
		if(truth!=0x0){
			l=truth->P4();
			costht= truth->GetMomentum().CosTheta();
			mother_aProt=truth->TheMother();
		}
			 
		int moth_aProt = (0x0==mother_aProt)? 88888 : mother_aProt->PdgCode();
		
		ntpAntiProton->Column("Mother", (Float_t) moth_aProt);
		ntpAntiProton->Column("antiProton_CosTheta", (Float_t) antiProton[aProt]->GetMomentum().CosTheta());
		
		qa.qaP4("antiProton_MC_", l, ntpAntiProton);
		qa.qaCand("antiProton_MC_", truth, ntpAntiProton);
		ntpAntiProton->Column("antiProton_MC_CosTheta", (Float_t) costht);
         

        ntpAntiProton->DumpData();
    }
	
    for (int kMin=0; kMin<kaonMinus.GetLength(); ++kMin){
         ntpKaonMinus->Column("ev",     (Float_t) evt);
         ntpKaonMinus->Column("cand",    (Float_t) kMin);
         ntpKaonMinus->Column("ncand",   (Float_t) kaonMinus.GetLength());
         ntpKaonMinus->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(kaonMinus[kMin]));

         qa.qaP4("kaonMinus_", kaonMinus[kMin]->P4(), ntpKaonMinus);
         qa.qaCand("kaonMinus_", kaonMinus[kMin], ntpKaonMinus);

         jenny::numberOfHitsInSubdetector("kaonMinus_", kaonMinus[kMin], ntpKaonMinus);
         jenny::tagNHits("kaonMinus_", kaonMinus[kMin], ntpKaonMinus);

 		RhoCandidate * truth = kaonMinus[kMin]->GetMcTruth();
 		RhoCandidate * mother_kMin=0;
 		TLorentzVector l;
 		float costht = -999.;
 		
 		if(truth!=0x0){
 			l=truth->P4();
 			costht= truth->GetMomentum().CosTheta();
 			mother_kMin=truth->TheMother();
 		}
 			 
 		int moth_kMin = (0x0==mother_kMin)? 88888 : mother_kMin->PdgCode();
 		
 		ntpKaonMinus->Column("Mother", (Float_t) moth_kMin);
 		ntpKaonMinus->Column("kaonMinus_CosTheta", (Float_t) kaonMinus[kMin]->GetMomentum().CosTheta());
 		
 		qa.qaP4("kaonMinus_MC_", l, ntpKaonMinus);
 		qa.qaCand("kaonMinus_MC_", truth, ntpKaonMinus);
 		ntpKaonMinus->Column("kaonMinus_MC_CosTheta", (Float_t) costht);
          

        ntpKaonMinus->DumpData();
     }
    


    //***Lambda0 -> PiMinus + Proton
	
    lambda0.Combine(piminus,proton);
	lambda0.Select(lambdaMassSelector);
    lambda0.SetType(kl0);

    std::map<int,int> bestVtxFitLambda0, bestMassFitLambda0;

    bestVtxFitLambda0 = jenny::VertexQaIndex(&lambda0);
    bestMassFitLambda0 = jenny::MassFitQaIndex(&lambda0, m0_lambda0);


    for (int j=0; j<lambda0.GetLength(); ++j){

    	
      //general info about event
      ntpLambda0->Column("ev",     (Float_t) evt);
      ntpLambda0->Column("cand",    (Float_t) j);
      ntpLambda0->Column("ncand",   (Float_t) lambda0.GetLength());
      ntpLambda0->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(lambda0[j]));
      ntpLambda0->Column("Lambda0_Pdg", (Float_t) lambda0[j]->PdgCode());

      RhoCandidate * truth = lambda0[j]->GetMcTruth();
	  RhoCandidate * mother_lamb=0;


	  if (truth)  mother_lamb=truth->TheMother();
  	  int moth = (mother_lamb==0x0) ? 88888 : mother_lamb->PdgCode();


	  ntpLambda0->Column("Mother", (Float_t) moth);

      qa.qaP4("Lambda0_", lambda0[j]->P4(), ntpLambda0);
      qa.qaComp("Lambda0_", lambda0[j], ntpLambda0);

      int tag = 0;
      int ndau = lambda0[j]->NDaughters();
      int dtag[2]={0,0};

  	
      for(int dau=0; dau<ndau; dau++){
    	  RhoCandidate * daughter = lambda0[j]->Daughter(dau);
    	  dtag[dau] = jenny::tagHits(daughter);
      }

      if(dtag[0]==1 && dtag[1]==1) tag=1;


      ntpLambda0->Column("Lambda0_HitTag", (Int_t) tag);

  	
      // do vertex fit
      PndKinVtxFitter vertexfitterLambda0 (lambda0[j]);
	  vertexfitterLambda0.Fit();
      RhoCandidate * lambda0Fit = lambda0[j]->GetFit();


      // store info of vertex fit
      qa.qaFitter("VtxFit_", &vertexfitterLambda0, ntpLambda0);
      ntpLambda0->Column("VtxFit_HowGood", (Int_t) bestVtxFitLambda0[j]);
      qa.qaVtx("VtxFit_", lambda0Fit, ntpLambda0);
      qa.qaCand("VtxFit_", lambda0Fit, ntpLambda0);
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
      ntpLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitLambda0[j]);

  	

      TLorentzVector l;
      
	    if(0x0 != truth){
	    	l = truth->P4();
	    	qa.qaVtx("McTruth_", truth, ntpLambda0);
	    	qa.qaCand("McTruth_", truth, ntpLambda0);
	    	
	    }
	    else{
	    
	    	qa.qaVtx("McTruth_", dummyCand, ntpLambda0);
	    	qa.qaCand("McTruth_", dummyCand, ntpLambda0);

	    }

      qa.qaP4("McTruth_", l, ntpLambda0);


      //*** use for Xi only bestChi2Cand

      if (bestVtxFitLambda0[j]==1 && bestMassFitLambda0[j]>0 && tag==1){
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
    antiLambda0.SetType(kal0);

    std::map<int,int> bestVtxFitAntiLambda0, bestMassFitAntiLambda0;

    bestVtxFitAntiLambda0 = jenny::VertexQaIndex(&antiLambda0);
    bestMassFitAntiLambda0 = jenny::MassFitQaIndex(&antiLambda0, m0_lambda0);


    for (int j=0; j<antiLambda0.GetLength(); ++j){

      //general info about event
      ntpAntiLambda0->Column("ev",     (Float_t) evt);
      ntpAntiLambda0->Column("cand",    (Float_t) j);
      ntpAntiLambda0->Column("ncand",   (Float_t) antiLambda0.GetLength());
      ntpAntiLambda0->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(antiLambda0[j]));
      ntpAntiLambda0->Column("AntiLambda0_Pdg", (Float_t) antiLambda0[j]->PdgCode());

      RhoCandidate * truth = antiLambda0[j]->GetMcTruth();
	  RhoCandidate * mother_ala =0;

	  if(truth) mother_ala = truth->TheMother();

	  int moth = (mother_ala==0x0) ? 88888 : mother_ala->PdgCode();
	  ntpAntiLambda0->Column("Mother", (Float_t) moth);

      qa.qaP4("AntiLambda0_", antiLambda0[j]->P4(), ntpAntiLambda0);
      qa.qaCand("AntiLambda0_", antiLambda0[j], ntpAntiLambda0);
      qa.qaComp("AntiLambda0_", antiLambda0[j], ntpAntiLambda0);
      

      int tag = 0;
      int ndau = antiLambda0[j]->NDaughters();
      int dtag[2]={0,0};


      for(int dau=0; dau<ndau; dau++){
    	  	 RhoCandidate * daughter = antiLambda0[j]->Daughter(dau);
    	  	 dtag[dau] = jenny::tagHits(daughter);
      }

      if(dtag[0]==1 && dtag[1]==1) tag=1;


      ntpAntiLambda0->Column("AntiLambda0_HitTag", (Int_t) tag);



      // do vertex fit
      PndKinVtxFitter vertexfitterAntiLambda0 (antiLambda0[j]);
      vertexfitterAntiLambda0.Fit();
      RhoCandidate * antiLambda0Fit = antiLambda0[j]->GetFit();



      // store info of vertex fit
      qa.qaFitter("VtxFit_", &vertexfitterAntiLambda0, ntpAntiLambda0);
      ntpAntiLambda0->Column("VtxFit_HowGood", (Int_t) bestVtxFitAntiLambda0[j]);
      qa.qaVtx("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
      qa.qaComp("VtxFit_", antiLambda0Fit, ntpAntiLambda0);

      qa.qaMcDiff("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
      jenny::qaVtxDiff("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
      jenny::qaMomRes("VtxFit_", antiLambda0Fit, ntpAntiLambda0);


      // do mass fit
      PndKinFitter massFitterAntiLambda0(antiLambda0Fit);
      massFitterAntiLambda0.AddMassConstraint(m0_lambda0);
      massFitterAntiLambda0.Fit();

      RhoCandidate * antiLambda0Fit_mass = antiLambda0Fit->GetFit();
      qa.qaFitter("MassFit_", &massFitterAntiLambda0, ntpAntiLambda0);
      ntpAntiLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitAntiLambda0[j]);




      TLorentzVector l;
      if(0x0 != truth){
			l = truth->P4();
			qa.qaVtx("MCTruth_", truth, ntpAntiLambda0);
			qa.qaCand("MCTruth_", truth, ntpAntiLambda0);

      }
      else{
    	  	qa.qaVtx("McTruth_", dummyCand, ntpAntiLambda0);
    	  	qa.qaCand("McTruth_", dummyCand, ntpAntiLambda0);

      }

      qa.qaP4("MCTruth_", l, ntpAntiLambda0);


      //***information of boosted particle
      antiLambda0Fit->Boost(-beamBoost);
      qa.qaComp("boost_", antiLambda0Fit, ntpAntiLambda0);



      if(bestVtxFitAntiLambda0[j]==1 && bestMassFitAntiLambda0[j]>0 && tag==1){
		  AntiLambda0Fit.Append(antiLambda0Fit);
		  jenny::CombinedList(antiLambda0Fit, &CombinedPiPlus, 211);
      }

      ntpAntiLambda0->DumpData();
    }

    jenny::GetNotCombinedList(CombinedPiPlus, &NotCombinedPiPlus);
    CombinedPiPlus.Cleanup();





    //*** Xi(1820)- -> Lambda0 + K-
	ximinus.Combine(Lambda0Fit, kaonMinus);
	ximinus.Select(xi1820MassSelector);
	ximinus.SetType(kXim);

	std::map<int,int> BestVtxFitXiMinus, BestMassFitXiMinus;

	BestVtxFitXiMinus = jenny::VertexQaIndex(&ximinus);
	BestMassFitXiMinus = jenny::MassFitQaIndex(&ximinus, m0_Xi1820);


    for (int j=0; j<ximinus.GetLength(); ++j){

      //general info about event
      ntpXiMinus->Column("ev",     (Float_t) evt);
      ntpXiMinus->Column("cand",    (Float_t) j);
      ntpXiMinus->Column("ncand",   (Float_t) ximinus.GetLength());
      ntpXiMinus->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(ximinus[j]));
      ntpXiMinus->Column("XiMinus_Pdg", (Float_t) ximinus[j]->PdgCode());

      RhoCandidate * truth = ximinus[j]->GetMcTruth();
	  RhoCandidate * mother_xim = 0;

	  if(truth) mother_xim = truth->TheMother();

	  int moth = (mother_xim==0x0) ? 88888 : mother_xim->PdgCode();
	  ntpXiMinus->Column("Mother", (Float_t) moth);

      qa.qaP4("XiMinus_", ximinus[j]->P4(), ntpXiMinus);
      qa.qaCand("XiMinus_", ximinus[j], ntpXiMinus);
      qa.qaComp("XiMinus_", ximinus[j], ntpXiMinus);
      qa.qaPoca("XiMinus_", ximinus[j], ntpXiMinus);



      // do vertex-fit

      PndKinVtxFitter vertexfitterXiMinus (ximinus[j]);
      vertexfitterXiMinus.Fit();
      RhoCandidate * ximinusFit = ximinus[j]->GetFit();


      // store info of vertex-fit

      qa.qaFitter("VtxFit_", &vertexfitterXiMinus, ntpXiMinus);
      ntpXiMinus->Column("VtxFit_HowGood", (Int_t) BestVtxFitXiMinus[j]);

      qa.qaVtx("VtxFit_", ximinusFit, ntpXiMinus);
      qa.qaCand("VtxFit_", ximinusFit, ntpXiMinus);
      qa.qaComp("VtxFit_", ximinusFit, ntpXiMinus);


	  // difference to MCTruth
      qa.qaMcDiff("VtxFit_", ximinusFit, ntpXiMinus);
      jenny::qaVtxDiff("VtxFit_", ximinusFit, ntpXiMinus);
      jenny::qaMomRes("VtxFit_", ximinusFit, ntpXiMinus);


      // do mass fit
      PndKinFitter massFitterXiMinus(ximinusFit);
      massFitterXiMinus.AddMassConstraint(m0_lambda0);
      massFitterXiMinus.Fit();

      RhoCandidate * ximinusFit_mass = ximinusFit->GetFit();
      qa.qaFitter("MassFit_", &massFitterXiMinus, ntpXiMinus);
      ntpXiMinus->Column("MassFit_HowGood", (Int_t) BestMassFitXiMinus[j]);

      qa.qaMcDiff("MassFit_", ximinusFit_mass, ntpXiMinus);
      jenny::qaVtxDiff("MassFit_", ximinusFit, ntpXiMinus);



      TLorentzVector l;

      if(0x0 != truth){
				l = truth->P4();
				qa.qaVtx("MCTruth_", truth, ntpXiMinus);
				qa.qaCand("MCTruth_", truth, ntpXiMinus);

      }
      else{
    	  qa.qaVtx("MCTruth_", dummyCand, ntpXiMinus);
    	  qa.qaCand("MCTruth_", dummyCand, ntpXiMinus);

      }

      qa.qaP4("MCTruth_", l, ntpXiMinus);


      if (BestVtxFitXiMinus[j]==1 && BestMassFitXiMinus[j]>0){
    	  XiMinusFit.Append(ximinusFit);
      }


      //***information of boosted particle
      ximinusFit->Boost(-beamBoost);
      qa.qaComp("boost_", ximinusFit, ntpXiMinus);

      ntpXiMinus->DumpData();


	}
    Lambda0Fit.Cleanup();



	//*** Xi+ -> AntiLambda0 + Pi+
	xiplus.Combine(AntiLambda0Fit,NotCombinedPiPlus);
	xiplus.Select(xiMassSelector);
	xiplus.SetType(kaXip);

	std::map<int,int> BestVtxFitXiPlus, BestMassFitXiPlus;

	BestVtxFitXiPlus = jenny::VertexQaIndex(&xiplus);
	BestMassFitXiPlus = jenny::MassFitQaIndex(&xiplus, m0_Xi);

    for (int j=0; j<xiplus.GetLength(); ++j){

      //general info about event
      ntpXiPlus->Column("ev",     (Float_t) evt);
      ntpXiPlus->Column("cand",    (Float_t) j);
      ntpXiPlus->Column("ncand",   (Float_t) xiplus.GetLength());
      ntpXiPlus->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(xiplus[j]));

      RhoCandidate * truth = xiplus[j]->GetMcTruth();
      RhoCandidate * mother =0;

      if(truth) mother = truth->TheMother();

      int moth = (mother==0x0) ? 88888 : mother->PdgCode();
      ntpXiPlus->Column("Mother", (Float_t) moth);

      qa.qaP4("Xiplus_", xiplus[j]->P4(), ntpXiPlus);
      qa.qaCand("Xiplus_", xiplus[j], ntpXiPlus);
      qa.qaComp("Xiplus_", xiplus[j], ntpXiPlus);



      //******** do vertex-fit
      PndKinVtxFitter vertexfitterxiplus (xiplus[j]);
      vertexfitterxiplus.Fit();
      RhoCandidate * xiplusFit = xiplus[j]->GetFit();


      // store info of vertex-fit
      qa.qaFitter("VtxFit_", &vertexfitterxiplus, ntpXiPlus);
      ntpXiPlus->Column("VtxFit_HowGood", (Int_t) BestVtxFitXiPlus[j]);
      qa.qaVtx("VtxFit_", xiplusFit, ntpXiPlus);
      qa.qaCand("VtxFit_", xiplusFit, ntpXiPlus);
      qa.qaComp("VtxFit_", xiplusFit, ntpXiPlus);

	  // difference to MCTruth
      qa.qaMcDiff("VtxFit_", xiplusFit, ntpXiPlus);
      jenny::qaVtxDiff("VtxFit_", xiplusFit, ntpXiPlus);
      jenny::qaMomRes("VtxFit_", xiplusFit, ntpXiPlus);


      //****** do mass fit
      PndKinFitter massFitterxiplus(xiplusFit);
      massFitterxiplus.AddMassConstraint(m0_lambda0);
      massFitterxiplus.Fit();

      RhoCandidate * xiplusFit_mass = xiplusFit->GetFit();
      qa.qaFitter("MassFit_", &massFitterxiplus, ntpXiPlus);
      ntpXiPlus->Column("MassFit_HowGood", (float) BestMassFitXiPlus[j]);
      qa.qaVtx("MassFit_", xiplusFit_mass, ntpXiPlus);

      qa.qaMcDiff("MassFit_", xiplusFit_mass, ntpXiPlus);
      jenny::qaVtxDiff("MassFit_", xiplusFit_mass, ntpXiPlus);
      jenny::qaMomRes("MassFit_", xiplusFit_mass, ntpXiPlus);


      TLorentzVector l;
      if(0x0 != truth){
				l = truth->P4();
				qa.qaVtx("MCTruth_", truth, ntpXiPlus);
				qa.qaCand("MCTruth_", truth, ntpXiPlus);

      }
      else{
    	  qa.qaVtx("MCTruth_", dummyCand, ntpXiPlus);
    	  qa.qaCand("MCTruth_", dummyCand, ntpXiPlus);

      }

      qa.qaP4("MCTruth_", l, ntpXiPlus);


      if(BestVtxFitXiPlus[j]==1 && BestMassFitXiPlus[j]>0){
    	  XiPlusFit.Append(xiplusFit);
      }

      //***information of boosted particle
      xiplusFit->Boost(-beamBoost);
      qa.qaComp("boost_", xiplusFit, ntpXiPlus);

      ntpXiPlus->DumpData();
	 }

    AntiLambda0Fit.Cleanup();
    NotCombinedPiPlus.Cleanup();



    //******* Xi+ Xi(1820)- System*****************************

    xiSys.Combine(XiPlusFit, XiMinusFit);
    xiSys.SetType(88888);
    
	std::map<int,int> BestVtxFitXiSys, Best4CFitXiSys;

	BestVtxFitXiSys = jenny::VertexQaIndex(&xiSys);
	Best4CFitXiSys = jenny::FourConstraintFitQaIndex(&xiSys, mom);

    for (int syscand=0; syscand<xiSys.GetLength(); ++syscand){
    	
    	cout << "XiSys" << endl;
    	
    	cout << "inforamtion"<< endl;
    	
		ntpXiSys->Column("ev",     (Float_t) evt);
		ntpXiSys->Column("cand",    (Float_t) j);
		ntpXiSys->Column("ncand",   (Float_t) ximinus.GetLength());
		ntpXiSys->Column("McTruthMatch", (bool) theAnalysis->McTruthMatch(xiSys[syscand]));


		qa.qaP4("XiSys_", xiSys[syscand]->P4(), ntpXiSys);
		qa.qaCand("XiSys_", xiSys[syscand], ntpXiSys);
		qa.qaComp("XiSys_", xiSys[syscand], ntpXiSys);
		qa.qaPoca("XiSys_", xiSys[syscand], ntpXiSys);


		RhoCandidate *  truth = xiSys[syscand]->GetMcTruth();
		TLorentzVector l;

		if (truth != 0x0){
			qa.qaCand("McTruth_", truth, ntpXiSys);
			qa.qaVtx("McTruth_", truth, ntpXiSys);
			l = truth->P4();
		}
		else{
			qa.qaCand("McTruth_", dummyCand, ntpXiSys);
			qa.qaVtx("McTruth_", dummyCand, ntpXiSys);

		}
		qa.qaP4("McTruth_", l, ntpXiSys);

		//Vertex Fitter
		


		PndKinVtxFitter vertexfitter (xiSys[syscand]);
		vertexfitter.Fit();

		RhoCandidate * xiSysVtxFit = xiSys[syscand]->GetFit();

		qa.qaFitter("VtxFit_", &vertexfitter, ntpXiSys);
		ntpXiSys->Colum("VtxFit_HowGood", (Int_t) BestVtxFitXiSys[syscand]);
		
		qa.qaCand("VtxFit_", xiSysVtxFit, ntpXiSys);
		qa.qaComp("VtxFit_", xiSysVtxFit, ntpXiSys);
		qa.qaVtx("VtxFit_", xiSysVtxFit, ntpXiSys);
		

		//4C-Fitter
		


		PndKinFitter fitter4c (xiSys[syscand]);
		fitter4c.Add4MomConstraint(ini);
		fitter4c.Fit();

		RhoCandidate * xiSysFit4c = xiSys[syscand]->GetFit();

		qa.qaFitter("4CFit_", &fitter4c, ntpXiSys);
		ntpXiSys->Colum("4CFit_HowGood", (Int_t) Best4CFitXiSys[syscand]);
		qa.qaCand("4cFit_", xiSysFit4c, ntpXiSys);
		qa.qaComp("4cFit_", xiSysFit4c, ntpXiSys);
		qa.qaVtx("4CFit_", xiSysFit4c, ntpXiSys);


		ntpXiSys->DumpData();


    }
    XiMinusFit.Cleanup();
    XiPlusFit.Cleanup();
  }



   

  //Write output
  out->cd();

  ntpMC -> GetInternalTree()->Write();
  ntpPiMinus ->GetInternalTree()->Write();
  ntpPiPlus->GetInternalTree()->Write();
  ntpProton->GetInternalTree()->Write();
  ntpAntiProton->GetInternalTree()->Write();
  ntpKaonMinus->GetInternalTree()->Write();
  ntpLambda0->GetInternalTree()->Write();
  ntpAntiLambda0->GetInternalTree()->Write();
  ntpXiMinus->GetInternalTree()->Write();
  ntpXiPlus->GetInternalTree()->Write();
  ntpXiSys->GetInternalTree()->Write();

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
