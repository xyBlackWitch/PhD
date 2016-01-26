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

using std::cout;
using std::endl;


enum pidNumbers {
	kPip = 211, kPim = -211,
	kPp = 2212, kaPm = -2212,
	kl0 = 3122, kal0 = -3122,
	kKm = -321,
	kXim = 23314, kaXip = -3312
};

void CombinedList(RhoCandidate* cand, RhoCandList* combinedList, int pdg){
	  /**
	   * @brief: gives back a list of already combined particles
	   * @details: The function creates a list of already combined particles for the analysis
	   */
	  for (int daughter=0; daughter<cand->NDaughters(); daughter++){
		  RhoCandidate * daughterCand = cand->Daughter(daughter);
		  if (daughterCand->PdgCode()==pdg){
			  combinedList->Append(daughterCand);
		  }

	  }

	  combinedList->RemoveClones();


}

void GetNotCombinedList(RhoCandList combinedList, RhoCandList * candList){
	  for (int j=0; j<combinedList.GetLength(); j++){
		  RhoCandidate * combinedCand = combinedList[j];
		  candList->Remove(combinedCand);
	  }
}


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
	 * 0: if there is no hit in the detector
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



std::map<int,int> VertexQaIndex(RhoCandList* candList, float probLimit=0.01){
	  /** @brief  give back the order of the best chi2
	   * @details give back the order of the best chi2!  1 means best, 2: second best (same with negative valuesfor bad chi2 )
	   */

	  std::map<double, int> chi2_good, chi2_bad;

	  for (int j=0; j<candList->GetLength(); ++j){

		  PndKinVtxFitter vtxfitter(candList->Get(j));
		  vtxfitter.Fit();

		  bool failedchi2 = TMath::IsNaN(vtxfitter.GetChi2());
		  bool failedprob = TMath::IsNaN(vtxfitter.GetProb());

		  if(!failedchi2 && !failedprob){

			  if (vtxfitter.GetProb() > probLimit){ //Prob > 0.01
				  chi2_good[vtxfitter.GetChi2()]=j;
			  }
			  else{ //Prob <= 0.01
				  chi2_bad[vtxfitter.GetChi2()]=j;
			  }

		  }
	  }

	  std::map<double, int>::iterator is_good, is_bad;
	  std::map<int, int> indexBestFit;
	  int running = 0;

	  for (is_good = chi2_good.begin(); is_good != chi2_good.end(); is_good++, running++){
		   indexBestFit[is_good->second] = running + 1;
	  }

	  running =0;

	  for (is_bad = chi2_bad.begin(); is_bad != chi2_bad.end(); is_bad++, running++){
		  indexBestFit[is_bad->second] = - (running + 1);
	  }


	  return indexBestFit;
}

std::map<int,int> MassFitQaIndex(RhoCandList* candList, float m0, float probLimit=0.01){
	  /** @brief  give back the order of the best chi2 for MassFit
	   * @details give back the order of the best chi2 for the MassFit!  1 means best, 2: second best (analoge for bad chi2 with negative values)
	   */

	  if(m0==0) std::cout << "Mass is missing for mass fit" << std::endl;

	  std::map<double, int> chi2_good, chi2_bad;

	  for (int i=0; i<candList->GetLength(); i++){

		  PndKinFitter massfitter(candList->Get(i));
		  massfitter.AddMassConstraint(m0);
		  massfitter.Fit();

		  bool failedchi2 = TMath::IsNaN(massfitter.GetChi2());
		  bool failedprob = TMath::IsNaN(massfitter.GetProb());

		  if(!failedchi2 && !failedprob){

			  if (massfitter.GetProb() > probLimit){
				  chi2_good[massfitter.GetChi2()]=i;
			  }
			  else{
				  chi2_bad[massfitter.GetChi2()]=i;
			  }
		  }
	  }

	  std::map<double,int>::iterator is_good, is_bad;
	  std::map<int,int> bestMassFit;

	  int run =0;

	  for (is_good = chi2_good.begin(); is_good != chi2_good.end(); is_good++, run++){
		  bestMassFit[is_good->second] = run + 1;
	  }

	  run = 0;

	  for (is_bad = chi2_bad.begin(); is_bad != chi2_bad.end(); is_bad++, run++){
		  bestMassFit[is_bad->second] = - (run + 1);
	  }


	  return bestMassFit;
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


void analysis_Xi1820(int nevts=0, double mom=4.6){

	  TStopwatch timer;
	
	  TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);

	  TString OutputFile = "analysis_output.root";

	  //Input simulation Files
	  TString inPIDFile = "pid_complete.root";
	  TString inParFile = "simparams.root";

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
	  ntpMc = new RhoTuple("ntpMC", "MCTruth info");
	  ntpPiMinus = new RhoTuple("ntpPiMinus", "PiMinus info");
	  ntpPiPlus = new RhoTuple("ntpPiPlus", "PiPlus info");
	  ntpProton = new RhoTuple("ntpProton", "Proton info");
	  ntpAntiProton = new RhoTuple("ntpAntiProton", "Antiproton info");
	  ntpKaonMinus = new RhoTuple("ntpKaonMinus", "KaonMinus info");
	  ntpLambda0 = new RhoTuple("ntpLambda0", "Lambda0 info");
	  ntpAntiLambda0 = new RhoTuple("ntpAntiLambda0", "AntiLambda0 info");
	  ntpXiMinus1820 = new RhoTuple("ntpXiMinus1820", "XiMinus info");
	  ntpXiPlus = new RhoTuple("ntpXiPlus", "XiPlus info");
	  ntpXiSys = new RhoTuple("ntpXiSys", "XiMinus XiPlus system info");

	  //Create output file for histograms

	  TFile * output = TFile::Open("output_ana.root","RECREATE");


	  // data reader Object
	  PndAnalysis * Analysis = new PndAnalysis();


	  //***Mass selector
	  double m0_lambda0= TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
	  cout<<"Mass of Lambda0: "<<m0_lambda0<<endl;
	  RhoMassParticleSelector * lambdaMassSelector = new RhoMassParticleSelector("lambda0", m0_lambda0, 0.3);

	  double m0_Xi = TDatabasePDG::Instance()->GetParticle("Xi-")->Mass();
	  cout<<"Mass of Xi-: "<<m0_Xi<<endl;
	  xiMassSelector = new RhoMassParticleSelector("Xi-", m0_Xi, 0.3);

	  double m0_Xi1820 = 1.823;
	  cout<<"Mass of Xi(1820)-: "<<m0_Xi1820<<endl;
	  xi1820MassSelector = new RhoMassParticleSelector("Xi(1820)-", m0_Xi1820, 0.3);

	  double m0_beam = TDatabasePDG::Instance()->GetParticle("pbarpSystem")->Mass();
	  cout<<"Mass pbar p system: "<<m0_beam<<endl;


	  //*** lorentz vector of the initial particle
	  cout << "BeamMomentum: " << mom << endl;
	  double p_m0 = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
	  TLorentzVector fini (0,0, mom, sqrt(p_m0*p_m0+ mom*mom)+p_m0);

	  if(nevts==0) nevts = Analysis->GetEntries();


		 //RhoCandLists for analysis
		RhoCandList piplus, piminus, lambda0, antiLambda0, proton, antiProton, kaonminus, xiplus, ximinus, xiSys;
		RhoCandList CombinedPiPlus, NotCombinedPiPlus; //NotCombinedPiMinus, CombinedPiMinus;
		RhoCandList Lambda0Fit, AntiLambda0Fit, XiMinusFit, XiPlusFit;
		RhoCandList mclist, all;

		//Dummy RhoCandidate
		RhoCandidate * dummyCand = new RhoCandidate();

		//****************************Analysis*************************************************
		TVector3 beamBoost = fini.BoostVector();
		PndRhoTupleQA qa(Analysis, mom);
		int fEvtCount = -1;


		while (Analysis->GetEvent() && ++fEvtCount<nevts){
			if ((fEvtCount%100)==0) cout << "evt "<< fEvtCount <<endl;

			//***get MC list and store info
			Analysis->FillList(mclist, "McTruth");
			qa.qaMcList("", mclist, ntpMc);
			ntpMc->DumpData();


			//***Setup event shape object

		    TString PidSelection = "PidAlgoIdealCharged";//"PidAlgoMvd;PidAlgoStt;PidAlgoDrc";

			Analysis->FillList(all, "All", PidSelection);
			PndEventShape evsh(all, fini, 0.05, 0.1);

			//***Selection
		    Analysis->FillList(piminus, "PionBestMinus", PidSelection);
		    Analysis->FillList(NotCombinedPiPlus, "PionBestPlus", PidSelection);
		    Analysis->FillList(piplus, "PionBestPlus", PidSelection);
		    Analysis->FillList(proton, "ProtonBestPlus", PidSelection);
		    Analysis->FillList(antiProton, "ProtonBestMinus", PidSelection);
		    Analysis->FillList(kaonminus, "KaonBestMinus", PidSelection);


		    for (int pip=0; pip<piplus.GetLength(); ++pip){
		        ntpPiPlus->Column("ev",     (Float_t) fEvtCount);
		        ntpPiPlus->Column("cand",    (Float_t) pip);
		        ntpPiPlus->Column("ncand",   (Float_t) piplus.GetLength());
		        ntpPiPlus->Column("McTruthMatch", (bool) Analysis->McTruthMatch(piplus[pip]));

		        qa.qaP4("piplus_", piplus[pip]->P4(), ntpPiPlus);
		        qa.qaCand("piplus_", piplus[pip], ntpPiPlus);

		        numberOfHitsInSubdetector("piplus_", piplus[pip], ntpPiPlus);
		        tagNHits("piplus_", piplus[pip], ntpPiPlus);

				RhoCandidate * mother_pip;
				RhoCandidate * truth = piplus[pip]->GetMcTruth();
				if (truth) mother_pip = truth->TheMother();

				int moth_pip;
				if (mother_pip==0x0){
				  moth_pip = 88888;
				}
				else
				  moth_pip = mother_pip->PdgCode();

				ntpPiPlus->Column("Mother", (Float_t) moth_pip);

				TLorentzVector l;
				float costheta = -999.;
				if(truth!=0x0){
				  l=truth->P4();
				  costheta = truth->GetMomentum().CosTheta();
				  qa.qaCand("piplus_MC_", piplus[pip]->GetMcTruth(), ntpPiPlus);
				}
				else{
				  qa.qaCand("piplus_MC_", dummyCand, ntpPiPlus);
				}

				qa.qaP4("piplus_MC_", l, ntpPiPlus);
				ntpPiPlus->Column("piplus_MC_CosTheta", (Float_t) costheta);

		        ntpPiPlus->DumpData();
		    }

		    for (int pim=0; pim<piminus.GetLength(); ++pim){
		        ntpPiMinus->Column("ev",     (Float_t) fEvtCount);
		        ntpPiMinus->Column("cand",    (Float_t) pim);
		        ntpPiMinus->Column("ncand",   (Float_t) piminus.GetLength());
		        ntpPiMinus->Column("McTruthMatch", (bool) Analysis->McTruthMatch(piminus[pim]));

		        qa.qaP4("piminus_", piminus[pim]->P4(), ntpPiMinus);
		        qa.qaCand("piminus_", piminus[pim], ntpPiMinus);

		        numberOfHitsInSubdetector("PiMinus_", piminus[pim], ntpPiMinus);
		        tagNHits("piminus_", piminus[pim], ntpPiMinus);

				RhoCandidate * mother_pim;
				RhoCandidate * truth = piminus[pim]->GetMcTruth();
				if (truth) mother_pim = truth->TheMother();

				int moth_pim;
				if (mother_pim==0x0){
				  moth_pim = 88888;
				}
				else
				  moth_pim = mother_pim->PdgCode();

				ntpPiMinus->Column("Mother", (Float_t) moth_pim);

				TLorentzVector l;
				float costheta = -999.;
				if(truth!=0x0){
				  l=truth->P4();
				  costheta = truth->GetMomentum().CosTheta();
				  qa.qaCand("piminus_MC_", piminus[pim]->GetMcTruth(), ntpPiMinus);
				}
				else{
				  qa.qaCand("piminus_MC_", dummyCand, ntpPiMinus);
				}

				qa.qaP4("piminus_MC_", l, ntpPiMinus);
				ntpPiMinus->Column("piminus_MC_CosTheta", (Float_t) costheta);

		        ntpPiMinus->DumpData();
		    }

		    for (int prot=0; prot<proton.GetLength(); ++prot){
		        ntpProton->Column("ev",     (Float_t) fEvtCount);
		        ntpProton->Column("cand",    (Float_t) prot);
		        ntpProton->Column("ncand",   (Float_t) proton.GetLength());
		        ntpProton->Column("McTruthMatch", (bool) Analysis->McTruthMatch(proton[prot]));

		        qa.qaP4("proton_", proton[prot]->P4(), ntpProton);
		        qa.qaCand("proton_", proton[prot], ntpProton);

		        numberOfHitsInSubdetector("proton_", proton[prot], ntpProton);
		        tagNHits("proton_", proton[prot], ntpProton);

				RhoCandidate * truth = proton[prot]->GetMcTruth();

				RhoCandidate * mother_prot;
				if (truth) mother_prot = truth->TheMother();
				int moth_prot;
				if (mother_prot==0x0){
				  moth_prot = 88888;
				}
				else
				  moth_prot = mother_prot->PdgCode();

				ntpProton->Column("MC_Mother_PDG", (Float_t) moth_prot);

				TLorentzVector l;
				float costheta = -999.;
				if(truth!=0x0){
				  l=truth->P4();
				  costheta = truth->GetMomentum().CosTheta();
				  qa.qaCand("proton_MC_", proton[prot]->GetMcTruth(), ntpProton);
				}
				else{
				  qa.qaCand("proton_MC_", dummyCand, ntpProton);
				}

				qa.qaP4("proton_MC_", l, ntpProton);
				ntpProton->Column("proton_MC_CosTheta", (Float_t) costheta);

		        ntpProton->DumpData();
		    }

		    for (int aProt=0; aProt<antiProton.GetLength(); ++aProt){
		        ntpAntiProton->Column("ev",     (Float_t) fEvtCount);
		        ntpAntiProton->Column("cand",    (Float_t) aProt);
		        ntpAntiProton->Column("ncand",   (Float_t) antiProton.GetLength());
		        ntpAntiProton->Column("McTruthMatch", (bool) Analysis->McTruthMatch(antiProton[aProt]));

		        qa.qaP4("AntiProton_", antiProton[aProt]->P4(), ntpAntiProton);
		        qa.qaCand("AntiProton_", antiProton[aProt], ntpAntiProton);

		        numberOfHitsInSubdetector("AntiProton_", antiProton[aProt], ntpAntiProton);
		        tagNHits("AntiProton_", antiProton[aProt], ntpAntiProton);

				RhoCandidate * mother_aprot;
				RhoCandidate * truth = antiProton[aProt]->GetMcTruth();
				if (truth) mother_aprot = truth->TheMother();

				int moth_aprot;
				if (mother_aprot==0x0){
				  moth_aprot = 88888;
				}
				else
				  moth_aprot = mother_aprot->PdgCode();

				ntpAntiProton->Column("Mother", (Float_t) moth_aprot);


				TLorentzVector l;
				float costheta = -999.;
				if(truth!=0x0){
				  l=truth->P4();
				  costheta = truth->GetMomentum().CosTheta();
				  qa.qaCand("AntiProton_MC_", antiProton[aProt]->GetMcTruth(), ntpAntiProton);
				}
				else{
				  qa.qaCand("AntiProton_MC_", dummyCand, ntpAntiProton);
				}

				qa.qaP4("AntiProton_MC_", l, ntpAntiProton);
				ntpAntiProton->Column("AntiProton_MC_CosTheta", (Float_t) costheta);

		        ntpAntiProton->DumpData();
		    }


		    for (int k=0; k<kaonminus.GetLength(); ++k){
				ntpKaonMinus->Column("ev",     (Float_t) fEvtCount);
				ntpKaonMinus->Column("cand",    (Float_t) k);
				ntpKaonMinus->Column("ncand",   (Float_t) kaonminus.GetLength());
				ntpKaonMinus->Column("McTruthMatch", (bool) Analysis->McTruthMatch(kaonminus[k]));

				qa.qaP4("kaonminus_", kaonminus[k]->P4(), ntpKaonMinus);
				qa.qaCand("kaonminus_", kaonminus[k], ntpKaonMinus);

				numberOfHitsInSubdetector("kaonminus_", kaonminus[k], ntpKaonMinus);
				tagNHits("kaonminus_", kaonminus[k], ntpKaonMinus);

				RhoCandidate * mother_k =0;
				RhoCandidate * truth = kaonminus[k]->GetMcTruth();


				TLorentzVector l;
				float costheta = -999.;
				if(truth!=0x0){
				  l=truth->P4();
				  mother_k =  truth->TheMother();
				  costheta = truth->GetMomentum().CosTheta();
				  qa.qaCand("kaonminus_MC_", kaonminus[k]->GetMcTruth(), ntpKaonMinus);
				}
				else{
				  qa.qaCand("kaonminus_MC_", dummyCand, ntpKaonMinus);
				}

				int moth_pip = (mother_k==0x0)? 88888 : mother_k->PdgCode();
				ntpKaonMinus->Column("Mother", (Float_t) moth_pip);

				qa.qaP4("kaonminus_MC_", l, ntpKaonMinus);
				ntpKaonMinus->Column("kaonminus_MC_CosTheta", (Float_t) costheta);

				ntpKaonMinus->DumpData();
			}


		    //***Lambda0 -> PiMinus + Proton

		     lambda0.Combine(piminus,proton);
		     lambda0.Select(lambdaMassSelector);
		     lambda0.SetType(kl0);

		     std::map<int,int> bestVtxFitLambda0, bestMassFitLambda0;

		     bestVtxFitLambda0 = VertexQaIndex(&lambda0);
		     bestMassFitLambda0 = MassFitQaIndex(&lambda0, m0_lambda0);


		     for (int j=0; j<lambda0.GetLength(); ++j){


		       //general info about event
		       ntpLambda0->Column("ev",     (Float_t) fEvtCount);
		       ntpLambda0->Column("cand",    (Float_t) j);
		       ntpLambda0->Column("ncand",   (Float_t) lambda0.GetLength());
		       ntpLambda0->Column("McTruthMatch", (bool) Analysis->McTruthMatch(lambda0[j]));
		       ntpLambda0->Column("Lambda0_Pdg", (Float_t) lambda0[j]->PdgCode());





		       qa.qaP4("Lambda0_", lambda0[j]->P4(), ntpLambda0);
		       qa.qaComp("Lambda0_", lambda0[j], ntpLambda0);

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

		       RhoCandidate * mother = lambda0Fit->TheMother();
		       int moth = (mother==0x0) ? 88888 : mother->PdgCode();
		       ntpLambda0->Column("Mother", (Float_t) moth);


		       // store info of vertex fit
		       qa.qaFitter("VtxFit_", &vertexfitterLambda0, ntpLambda0);
		       ntpLambda0->Column("VtxFit_HowGood", (Int_t) bestVtxFitLambda0[j]);
		       qa.qaVtx("VtxFit_", lambda0Fit, ntpLambda0);
		       qa.qaComp("VtxFit_", lambda0Fit, ntpLambda0);
		       qa.qaP4("VtxFit_", lambda0Fit->P4(), ntpLambda0);

		       // differenz to MCTruth
		        qa.qaMcDiff("VtxFit_", lambda0Fit, ntpLambda0);
		        qaVtxDiff("VtxFit_", lambda0Fit, ntpLambda0);
		        qaMomRes("VtxFit_", lambda0Fit, ntpLambda0);


			   // do Kalman vertex fit
			   PndKalmanVtxFitter kalmanfitterLambda0 (lambda0[j]);
			   kalmanfitterLambda0.Fit();
			   RhoCandidate * lambda0KalmanFit = lambda0[j]->GetFit();


			   // store info of vertex fit
			   qa.qaFitter("KalmanFit_", &kalmanfitterLambda0, ntpLambda0);
			   //ntpLambda0->Column("KalmanFit_HowGood", (Int_t) bestKalmanFitLambda0[j]);
			   qa.qaVtx("KalmanFit_", lambda0KalmanFit, ntpLambda0);
			   qa.qaComp("KalmanFit_", lambda0KalmanFit, ntpLambda0);
			   qa.qaP4("KalmanFit_", lambda0KalmanFit->P4(), ntpLambda0);

			   // differenz to MCTruth
				qa.qaMcDiff("KalmanFit_", lambda0KalmanFit, ntpLambda0);
				qaVtxDiff("KalmanFit_", lambda0KalmanFit, ntpLambda0);
				qaMomRes("KalmanFit_", lambda0KalmanFit, ntpLambda0);

		       // do mass fit
		       PndKinFitter massFitterLambda0(lambda0Fit);
		       massFitterLambda0.AddMassConstraint(m0_lambda0);
		       massFitterLambda0.Fit();

		       RhoCandidate * lambda0Fit_mass = lambda0Fit->GetFit();
		       qa.qaFitter("MassFit_", &massFitterLambda0, ntpLambda0);

		       ntpLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitLambda0[j]);


		       RhoCandidate * truth = lambda0[j]->GetMcTruth();
		       RhoCandidate * truthDaughter = lambda0[j]->Daughter(0)->GetMcTruth();
		       TLorentzVector l;
		       TVector3 dl;

		 	    if(0x0 != truth){
		 	    	l = truth->P4();
		 	    	qa.qaVtx("McTruth_", truth, ntpLambda0);
		 	    	dl = truth->Daughter(0)->Pos();
		 	    }
		 	    else{
		 	    	qa.qaVtx("McTruth_", dummyCand, ntpLambda0);
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

		     std::map<int,int> bestVtxFitantiLambda0, bestMassFitantiLambda0;

		     bestVtxFitantiLambda0 = VertexQaIndex(&antiLambda0);
		     bestMassFitantiLambda0 = MassFitQaIndex(&antiLambda0, m0_lambda0);


		     for (int j=0; j<antiLambda0.GetLength(); ++j){


		       //general info about event
		       ntpAntiLambda0->Column("ev",     (Float_t) fEvtCount);
		       ntpAntiLambda0->Column("cand",    (Float_t) j);
		       ntpAntiLambda0->Column("ncand",   (Float_t) antiLambda0.GetLength());
		       ntpAntiLambda0->Column("McTruthMatch", (bool) Analysis->McTruthMatch(antiLambda0[j]));
		       ntpAntiLambda0->Column("antiLambda0_Pdg", (Float_t) antiLambda0[j]->PdgCode());



		       qa.qaP4("antiLambda0_", antiLambda0[j]->P4(), ntpAntiLambda0);
		       qa.qaComp("antiLambda0_", antiLambda0[j], ntpAntiLambda0);

		       int tag = 0;
		       int ndau = antiLambda0[j]->NDaughters();
		       int dtag[2]={0,0};


		       for(int dau=0; dau<ndau; dau++){
		     	  RhoCandidate * daughter = antiLambda0[j]->Daughter(dau);
		     	  dtag[dau] = tagHits(daughter);
		       }

		       if(dtag[0]==1 && dtag[1]==1) tag=1;


		       ntpAntiLambda0->Column("antiLambda0_HitTag", (Int_t) tag);


		       // do vertex fit
		       PndKinVtxFitter vertexfitterantiLambda0 (antiLambda0[j]);
		       vertexfitterantiLambda0.Fit();
		       RhoCandidate * antiLambda0Fit = antiLambda0[j]->GetFit();

		       RhoCandidate * mother = antiLambda0Fit->TheMother();
		       int moth = (mother==0x0) ? 88888 : mother->PdgCode();

		       ntpAntiLambda0->Column("Mother", (Float_t) moth);

		       // store info of vertex fit
		       qa.qaFitter("VtxFit_", &vertexfitterantiLambda0, ntpAntiLambda0);
		       ntpAntiLambda0->Column("VtxFit_HowGood", (Int_t) bestVtxFitantiLambda0[j]);
		       qa.qaVtx("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
		       qa.qaComp("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
		       qa.qaP4("VtxFit_", antiLambda0Fit->P4(), ntpAntiLambda0);

		       // differenz to MCTruth
		        qa.qaMcDiff("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
		        qaVtxDiff("VtxFit_", antiLambda0Fit, ntpAntiLambda0);
		        qaMomRes("VtxFit_", antiLambda0Fit, ntpAntiLambda0);


			   // do Kalman vertex fit
			   PndKalmanVtxFitter kalmanfitterantiLambda0 (antiLambda0[j]);
			   kalmanfitterantiLambda0.Fit();
			   RhoCandidate * antiLambda0KalmanFit = antiLambda0[j]->GetFit();


			   // store info of vertex fit
			   qa.qaFitter("KalmanFit_", &kalmanfitterantiLambda0, ntpAntiLambda0);
			   //ntpAntiLambda0->Column("KalmanFit_HowGood", (Int_t) bestKalmanFitantiLambda0[j]);
			   qa.qaVtx("KalmanFit_", antiLambda0KalmanFit, ntpAntiLambda0);
			   qa.qaComp("KalmanFit_", antiLambda0KalmanFit, ntpAntiLambda0);
			   qa.qaP4("KalmanFit_", antiLambda0KalmanFit->P4(), ntpAntiLambda0);

			   // differenz to MCTruth
				qa.qaMcDiff("KalmanFit_", antiLambda0KalmanFit, ntpAntiLambda0);
				qaVtxDiff("KalmanFit_", antiLambda0KalmanFit, ntpAntiLambda0);
				qaMomRes("KalmanFit_", antiLambda0KalmanFit, ntpAntiLambda0);

		       // do mass fit
		       PndKinFitter massFitterantiLambda0(antiLambda0Fit);
		       massFitterantiLambda0.AddMassConstraint(m0_lambda0);
		       massFitterantiLambda0.Fit();

		       RhoCandidate * antiLambda0Fit_mass = antiLambda0Fit->GetFit();
		       qa.qaFitter("MassFit_", &massFitterantiLambda0, ntpAntiLambda0);

		       ntpAntiLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitantiLambda0[j]);


		       RhoCandidate * truth = antiLambda0[j]->GetMcTruth();
		       RhoCandidate * truthDaughter = antiLambda0[j]->Daughter(0)->GetMcTruth();
		       TLorentzVector l;
		       TVector3 dl;

		 	    if(0x0 != truth){
		 	    	l = truth->P4();
		 	    	qa.qaVtx("McTruth_", truth, ntpAntiLambda0);
		 	    	dl = truth->Daughter(0)->Pos();
		 	    }
		 	    else{
		 	    	qa.qaVtx("McTruth_", dummyCand, ntpAntiLambda0);
		 	    }

		       qa.qaP4("McTruth_", l, ntpAntiLambda0);


		       //*** use for Xi only bestChi2Cand

		       if (bestVtxFitantiLambda0[j]==1 && bestMassFitantiLambda0[j]>0 && tag==1){
		 		  AntiLambda0Fit.Append(antiLambda0Fit);
		 		  CombinedList(antiLambda0Fit, &CombinedPiPlus, 211);
		       }


		       //***information of boosted particle
		       antiLambda0Fit->Boost(-beamBoost);
		       qa.qaComp("boost_", antiLambda0Fit, ntpAntiLambda0);

		       ntpAntiLambda0->DumpData();


		      }

		      GetNotCombinedList(CombinedPiPlus, &NotCombinedPiPlus);
		      CombinedPiPlus.Cleanup();



		      //*** Xi(1820)- -> Lambda0 + K-
		   	ximinus.Combine(Lambda0Fit, kaonminus);
		   	ximinus.Select(xi1820MassSelector);
		   	ximinus.SetType(kXim);

		   	std::map<int,int> BestVtxFitXiMinus, BestMassFitXiMinus;

		   	BestVtxFitXiMinus = VertexQaIndex(&ximinus);
		   	BestMassFitXiMinus = MassFitQaIndex(&ximinus, m0_Xi1820);


			for (int j=0; j<ximinus.GetLength(); ++j){

				//general info about event
				ntpXiMinus1820->Column("ev",     (Float_t) fEvtCount);
				ntpXiMinus1820->Column("cand",    (Float_t) j);
				ntpXiMinus1820->Column("ncand",   (Float_t) ximinus.GetLength());
				ntpXiMinus1820->Column("McTruthMatch", (bool) Analysis->McTruthMatch(ximinus[j]));
				ntpXiMinus1820->Column("XiMinus_Pdg", (Float_t) ximinus[j]->PdgCode());




				qa.qaP4("XiMinus_", ximinus[j]->P4(), ntpXiMinus1820);
				qa.qaComp("XiMinus_", ximinus[j], ntpXiMinus1820);
				qa.qaPoca("XiMinus_", ximinus[j], ntpXiMinus1820);



				// do vertex-fit

				PndKinVtxFitter vertexfitterXiMinus (ximinus[j]);
				vertexfitterXiMinus.Fit();
				RhoCandidate * ximinusFit = ximinus[j]->GetFit();

				RhoCandidate * mother = ximinusFit->TheMother();

				int moth = (mother==0x0) ? 88888 : mother->PdgCode();
				ntpXiMinus1820->Column("Mother", (Float_t) moth);

				// store info of vertex-fit

				qa.qaFitter("VtxFit_", &vertexfitterXiMinus, ntpXiMinus1820);
				ntpXiMinus1820->Column("VtxFit_HowGood", (Int_t) BestVtxFitXiMinus[j]);

				qa.qaVtx("VtxFit_", ximinusFit, ntpXiMinus1820);
				qa.qaP4("VtxFit_", ximinusFit->P4(), ntpXiMinus1820);
				qa.qaComp("VtxFit_", ximinusFit, ntpXiMinus1820);


				// difference to MCTruth
				qa.qaMcDiff("VtxFit_", ximinusFit, ntpXiMinus1820);
				qaVtxDiff("VtxFit_", ximinusFit, ntpXiMinus1820);
				qaMomRes("VtxFit_", ximinusFit, ntpXiMinus1820);



				// do mass fit
				PndKinFitter massFitterXiMinus(ximinusFit);
				massFitterXiMinus.AddMassConstraint(m0_Xi1820);
				massFitterXiMinus.Fit();

				RhoCandidate * ximinusFit_mass = ximinusFit->GetFit();
				qa.qaFitter("MassFit_", &massFitterXiMinus, ntpXiMinus1820);
				ntpXiMinus1820->Column("MassFit_HowGood", (Int_t) BestMassFitXiMinus[j]);

				qa.qaMcDiff("MassFit_", ximinusFit_mass, ntpXiMinus1820);
				qaVtxDiff("MassFit_", ximinusFit, ntpXiMinus1820);

				RhoCandidate * truth = ximinus[j]->GetMcTruth();
				TLorentzVector l;

				if(0x0 != truth){
					l = truth->P4();
					qa.qaVtx("MCTruth_", truth, ntpXiMinus1820);
				}
				else{
					qa.qaVtx("MCTruth_", dummyCand, ntpXiMinus1820);
				}

				qa.qaP4("MCTruth_", l, ntpXiMinus1820);


				if (BestVtxFitXiMinus[j]==1){
					XiMinusFit.Append(ximinusFit);
				}


				//***information of boosted particle
				ximinusFit->Boost(-beamBoost);
				qa.qaComp("boost_", ximinusFit, ntpXiMinus1820);

				ntpXiMinus1820->DumpData();


			}
			Lambda0Fit.Cleanup();




			//*** Xi+ -> AntiLambda0 + Pi+
		   	xiplus.Combine(AntiLambda0Fit, NotCombinedPiPlus);
		   	xiplus.Select(xiMassSelector);
		   	xiplus.SetType(kaXip);

		   	std::map<int,int> BestVtxFitxiplus, BestMassFitxiplus;

		   	BestVtxFitxiplus = VertexQaIndex(&xiplus);
		   	BestMassFitxiplus = MassFitQaIndex(&xiplus, m0_Xi);


			for (int j=0; j<xiplus.GetLength(); ++j){

				//general info about event
				ntpXiPlus->Column("ev",     (Float_t) fEvtCount);
				ntpXiPlus->Column("cand",    (Float_t) j);
				ntpXiPlus->Column("ncand",   (Float_t) xiplus.GetLength());
				ntpXiPlus->Column("McTruthMatch", (bool) Analysis->McTruthMatch(xiplus[j]));
				ntpXiPlus->Column("xiplus_Pdg", (Float_t) xiplus[j]->PdgCode());



				qa.qaP4("xiplus_", xiplus[j]->P4(), ntpXiPlus);
				qa.qaComp("xiplus_", xiplus[j], ntpXiPlus);
				qa.qaPoca("xiplus_", xiplus[j], ntpXiPlus);



				// do vertex-fit

				PndKinVtxFitter vertexfitterxiplus (xiplus[j]);
				vertexfitterxiplus.Fit();
				RhoCandidate * xiplusFit = xiplus[j]->GetFit();


				RhoCandidate * mother = xiplusFit->TheMother();
				int moth = (mother==0x0) ? 88888 : mother->PdgCode();
				ntpXiPlus->Column("Mother", (Float_t) moth);

				// store info of vertex-fit

				qa.qaFitter("VtxFit_", &vertexfitterxiplus, ntpXiPlus);
				ntpXiPlus->Column("VtxFit_HowGood", (Int_t) BestVtxFitxiplus[j]);

				qa.qaVtx("VtxFit_", xiplusFit, ntpXiPlus);
				qa.qaP4("VtxFit_", xiplusFit->P4(), ntpXiPlus);



				// difference to MCTruth
				qa.qaMcDiff("VtxFit_", xiplusFit, ntpXiPlus);
				qaVtxDiff("VtxFit_", xiplusFit, ntpXiPlus);
				qaMomRes("VtxFit_", xiplusFit, ntpXiPlus);


				// do kalman vertex-fit

				PndKalmanVtxFitter kalmanfitterxiplus (xiplus[j]);
				kalmanfitterxiplus.Fit();
				RhoCandidate * xiplusKalmanFit = xiplus[j]->GetFit();


				// store info of vertex-fit

				qa.qaFitter("KalmanFit_", &kalmanfitterxiplus, ntpXiPlus);

				qa.qaVtx("KalmanFit_", xiplusKalmanFit, ntpXiPlus);
				qa.qaP4("KalmanFit_", xiplusKalmanFit->P4(), ntpXiPlus);



				// difference to MCTruth
				qa.qaMcDiff("KalmanFit_", xiplusKalmanFit, ntpXiPlus);
				qaVtxDiff("KalmanFit_", xiplusKalmanFit, ntpXiPlus);
				qaMomRes("KalmanFit_", xiplusKalmanFit, ntpXiPlus);



				// do mass fit
				PndKinFitter massFitterxiplus(xiplusFit);
				massFitterxiplus.AddMassConstraint(m0_Xi);
				massFitterxiplus.Fit();

				RhoCandidate * xiplusFit_mass = xiplusFit->GetFit();
				qa.qaFitter("MassFit_", &massFitterxiplus, ntpXiPlus);
				ntpXiPlus->Column("MassFit_HowGood", (Int_t) BestMassFitxiplus[j]);

				qa.qaMcDiff("MassFit_", xiplusFit_mass, ntpXiPlus);
				qaVtxDiff("MassFit_", xiplusFit, ntpXiPlus);


				RhoCandidate * truth = xiplus[j]->GetMcTruth();
				TLorentzVector l;

				if(0x0 != truth){
					l = truth->P4();
					qa.qaVtx("MCTruth_", truth, ntpXiPlus);
				}
				else{
					qa.qaVtx("MCTruth_", dummyCand, ntpXiPlus);
				}

				qa.qaP4("MCTruth_", l, ntpXiPlus);


				if (BestVtxFitxiplus[j]==1 && BestMassFitxiplus[j]>0){
					XiPlusFit.Append(xiplusFit);
				}


	//			***information of boosted particle
				xiplusFit->Boost(-beamBoost);
				qa.qaComp("boost_", xiplusFit, ntpXiPlus);

				ntpXiPlus->DumpData();
			 }

		    AntiLambda0Fit.Cleanup();
		    NotCombinedPiPlus.Cleanup();





	//	    ******* Xi+ Xi- System*****************************

		     xiSys.Combine(XiPlusFit, XiMinusFit);
		     xiSys.SetType(88888);

		     for (int syscand=0; syscand<xiSys.GetLength(); ++syscand){

		 		ntpXiSys->Column("ev",     (Float_t) fEvtCount);
		 		ntpXiSys->Column("cand",    (Float_t) syscand);
		 		ntpXiSys->Column("ncand",   (Float_t) xiSys.GetLength());
		 		ntpXiSys->Column("McTruthMatch", (bool) Analysis->McTruthMatch(xiSys[syscand]));


		 		qa.qaP4("XiSys_", xiSys[syscand]->P4(), ntpXiSys);
		 		qa.qaComp("XiSys_", xiSys[syscand], ntpXiSys);
		  		qa.qaPoca("XiSys_", xiSys[syscand], ntpXiSys);


		 		RhoCandidate *  truth = xiSys[syscand]->GetMcTruth();
		 		TLorentzVector l;
		 		TLorentzVector l0;
		 		TLorentzVector l1;
		 		TLorentzVector l0l0;
		 		TLorentzVector l0l1;
		 		TLorentzVector l1l0;
				TLorentzVector l1l1;

		 		if (truth != 0x0){
		 			qa.qaVtx("McTruth_", truth, ntpXiSys);
		 			l = truth->P4();

		 			RhoCandidate * d0 = truth->Daughter(0);
		 			RhoCandidate * d1 = truth->Daughter(1);

		 			l1 = d1->P4();
		 			l0 = d0->P4();

		 			int d0_pdg =d0->PdgCode();
		 			int d1_pdg =d1->PdgCode();

		 			RhoCandidate * d0d0 = d0->Daughter(0);
		 			RhoCandidate * d0d1 = d0->Daughter(1);
		 			RhoCandidate * d1d0 = d1->Daughter(0);
					RhoCandidate * d1d1 = d1->Daughter(1);

		 			l0l0 = d0d0->P4();
		 			l0l1 = d0d1->P4();
		 			l1l0 = d1d0->P4();
					l1l1 = d1d1->P4();

		 			int d0d0_pdg = d0d0->PdgCode();
		 			int d0d1_pdg = d0d1->PdgCode();
		 			int d1d0_pdg = d1d0->PdgCode();
		 			int d1d1_pdg = d1d1->PdgCode();

			 		ntpXiSys->Column("McTruth_d0_pdg",   (Int_t) d0_pdg);
			 		ntpXiSys->Column("McTruth_d0d0_pdg",   (Int_t) d0d0_pdg);
					ntpXiSys->Column("McTruth_d0d1_pdg",   (Int_t) d0d1_pdg);
			 		ntpXiSys->Column("McTruth_d1_pdg",   (Int_t) d1_pdg);
			 		ntpXiSys->Column("McTruth_d0d0_pdg",   (Int_t) d1d0_pdg);
					ntpXiSys->Column("McTruth_d0d1_pdg",   (Int_t) d1d1_pdg);



		 		}
		 		else{
		 			qa.qaVtx("McTruth_", dummyCand, ntpXiSys);
			 		ntpXiSys->Column("McTruth_d0_pdg",   0);
			 		ntpXiSys->Column("McTruth_d0d0_pdg",   0);
					ntpXiSys->Column("McTruth_d0d1_pdg",   0);
			 		ntpXiSys->Column("McTruth_d1_pdg",   0);
			 		ntpXiSys->Column("McTruth_d0d0_pdg",   0);
					ntpXiSys->Column("McTruth_d0d1_pdg",   0);
		 		}
		 		qa.qaP4("McTruth_", l, ntpXiSys);
		 		qa.qaP4("McTruth_d0_", l0, ntpXiSys);
		 		qa.qaP4("McTruth_d0d0_", l0l0, ntpXiSys);
				qa.qaP4("McTruth_d0d1_", l0l1, ntpXiSys);
				qa.qaP4("McTruth_d1_", l1, ntpXiSys);
		 		qa.qaP4("McTruth_d1d0_", l1l0, ntpXiSys);
				qa.qaP4("McTruth_d1d1_", l1l1, ntpXiSys);


		 		//4C-Fitter

		 		PndKinFitter fitter4c (xiSys[syscand]);
		 		fitter4c.Add4MomConstraint(fini);
		 		fitter4c.Fit();

		 		RhoCandidate * xiSysFit4c = xiSys[syscand]->GetFit();

		 		qa.qaFitter("4CFit_", &fitter4c, ntpXiSys);
		 		qa.qaComp("4cFit_", xiSysFit4c, ntpXiSys);
		 		qa.qaVtx("4CFit_", xiSysFit4c, ntpXiSys);
		 		qa.qaMcDiff("4CFit_", xiSysFit4c, ntpXiSys);
		 		qaVtxDiff("4CFit_", xiSysFit4c, ntpXiSys);
		 		qaMomRes("4CFit_", xiSysFit4c, ntpXiSys);


		 		ntpXiSys->DumpData();


		     }
		     XiMinusFit.Cleanup();
		     XiPlusFit.Cleanup();

	   }

		//Write output
		output->cd();
		ntpMc -> GetInternalTree()->Write();
		ntpPiMinus ->GetInternalTree()->Write();
		ntpPiPlus->GetInternalTree()->Write();
		ntpProton->GetInternalTree()->Write();
		ntpAntiProton->GetInternalTree()->Write();
		ntpKaonMinus->GetInternalTree()->Write();
		ntpLambda0->GetInternalTree()->Write();
		ntpAntiLambda0->GetInternalTree()->Write();
		ntpXiMinus1820->GetInternalTree()->Write();
		ntpXiPlus->GetInternalTree()->Write();
		ntpXiSys->GetInternalTree()->Write();


		output->Save();

		timer.Stop();
		Double_t rtime = timer.RealTime();
		Double_t ctime = timer.CpuTime();
		cout<<endl<<endl;
		cout<<"Macro finisched successfully."<<endl;
		cout<<"Realtime: "<<rtime<<" s, CPU time: "<<ctime<<" s"<<endl;

		exit(0);

}
