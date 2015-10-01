//class PndAnalysis;
//class PndAnaPidSelector;
//class RhoCandList;
//class RhoTuple;

#include "AnalysisTaskXi1820.h"

// C++ headers
#include <string>
#include <iostream>
#include <map>

//Fair headers
#include "FairLogger.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairParAsciiFileIo.h"

//ROOT headers
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TVector3.h"

//Rho headers
#include "RhoTuple.h"
#include "PndRhoTupleQA.h"
#include "RhoCandList.h"
#include "RhoCandidate.h"
#include "RhoMassParticleSelector.h"

//Analysis headers
#include "PndAnalysis.h"
#include "PndEventShape.h"
#include "PndKinVtxFitter.h"
#include "PndKalmanVtxFitter.h"
#include "PndKinFitter.h"
#include "PndPidCandidate.h"


using std::cout;
using std::endl;

ClassImp(AnalysisTaskXi1820)

//---------------Default constructor------------
AnalysisTaskXi1820::AnalysisTaskXi1820():
	FairTask("Panda Tutorial Analysis Task") {
	}
//------------------------------------------

//--------------- destructor ----------------
AnalysisTaskXi1820::~AnalysisTaskXi1820(){}
//-----------------------------------------


enum pidNumbers {
	kPip = 211, kPim = -211,
	kPp = 2212, kaPm = -2212,
	kl0 = 3122, kal0 = -3122,
	kKm = -321,
	kXim = 23314, kaXip = -3312
};

void AnalysisTaskXi1820::CombinedList(RhoCandidate* cand, RhoCandList* combinedList, int pdg){
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

void AnalysisTaskXi1820::GetNotCombinedList(RhoCandList combinedList, RhoCandList * candList){
	  for (int j=0; j<combinedList.GetLength(); j++){
		  RhoCandidate * combinedCand = combinedList[j];
		  candList->Remove(combinedCand);
	  }
}


void AnalysisTaskXi1820::numberOfHitsInSubdetector(TString pre, RhoCandidate *c, RhoTuple *n){

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

void AnalysisTaskXi1820::tagNHits(TString pre, RhoCandidate *c, RhoTuple *n){

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

int AnalysisTaskXi1820::tagHits(RhoCandidate *c){

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



std::map<int,int> AnalysisTaskXi1820::VertexQaIndex(RhoCandList* candList, float probLimit=0.01){
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

std::map<int,int> AnalysisTaskXi1820::MassFitQaIndex(RhoCandList* candList, float m0, float probLimit=0.01){
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

void AnalysisTaskXi1820::qaVtxDiff(TString pre, RhoCandidate * c, RhoTuple * n){

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

void AnalysisTaskXi1820::qaMomRes(TString pre, RhoCandidate * c, RhoTuple * n){

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

InitStatus AnalysisTaskXi1820::Init(){




  //*** create tuples
  fntpMc = new RhoTuple("ntpMC", "MCTruth info");
  fntpPiMinus = new RhoTuple("ntpPiMinus", "PiMinus info");
  fntpPiPlus = new RhoTuple("ntpPiPlus", "PiPlus info");
  fntpProton = new RhoTuple("ntpProton", "Proton info");
  fntpAntiProton = new RhoTuple("ntpAntiProton", "Antiproton info");
  fntpKaonMinus = new RhoTuple("ntpKaonMinus", "KaonMinus info");
  fntpLambda0 = new RhoTuple("ntpLambda0", "Lambda0 info");
  fntpAntiLambda0 = new RhoTuple("ntpAntiLambda0", "AntiLambda0 info");
  fntpXiMinus1820 = new RhoTuple("ntpXiMinus1820", "XiMinus info");
  fntpXiPlus = new RhoTuple("ntpXiPlus", "XiPlus info");
  fntpXiSys = new RhoTuple("ntpXiSys", "XiMinus XiPlus system info");

  //Create output file for histograms
  TString outpath;
  outpath.Append(fOutPath);
  foutput = TFile::Open( outpath + "output_ana.root","RECREATE");


  // data reader Object
  fAnalysis = new PndAnalysis();
  

  //***Mass selector
  fm0_lambda0= TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
  cout<<"Mass of Lambda0: "<<fm0_lambda0<<endl;
  lambdaMassSelector = new RhoMassParticleSelector("lambda0", fm0_lambda0, 0.3);

  fm0_Xi = TDatabasePDG::Instance()->GetParticle("Xi-")->Mass();
  cout<<"Mass of Xi-: "<<fm0_Xi<<endl;
  xiMassSelector = new RhoMassParticleSelector("Xi-", fm0_Xi, 0.3);

  fm0_Xi1820 = 1.823;
  cout<<"Mass of Xi(1820)-: "<<fm0_Xi1820<<endl;
  xi1820MassSelector = new RhoMassParticleSelector("Xi(1820)-", fm0_Xi1820, 0.3);

  fm0_beam = TDatabasePDG::Instance()->GetParticle("pbarpSystem")->Mass();
  cout<<"Mass pbar p system: "<<fm0_beam<<endl;
  

  //*** lorentz vector of the initial particle
  cout << "BeamMomentum: " << fmom << endl;
  double p_m0 = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
  fini.SetXYZT(0,0, fmom, sqrt(p_m0*p_m0+ fmom*fmom)+p_m0);
  
  if(fnevts==0) fnevts = fAnalysis->GetEntries();


  return kSUCCESS;

}

void AnalysisTaskXi1820::Exec(Option_t* op)
{
	 //RhoCandLists for analysis
	RhoCandList piplus, piminus, lambda0, antiLambda0, proton, antiProton, kaonminus, xiplus, ximinus, xiSys;
	RhoCandList CombinedPiPlus, NotCombinedPiPlus; //NotCombinedPiMinus, CombinedPiMinus;
	RhoCandList Lambda0Fit, AntiLambda0Fit, XiMinusFit, XiPlusFit;
	RhoCandList mclist, all;

	//Dummy RhoCandidate
	RhoCandidate * dummyCand = new RhoCandidate();

	//****************************Analysis*************************************************
	TVector3 beamBoost = fini.BoostVector();
	PndRhoTupleQA qa(fAnalysis, fmom);
	fEvtCount = -1;


	while (fAnalysis->GetEvent() && ++fEvtCount<fnevts){
		if ((fEvtCount%100)==0) cout << "evt "<< fEvtCount <<endl;

		//***get MC list and store info
		fAnalysis->FillList(mclist, "McTruth");
		qa.qaMcList("", mclist, fntpMc);
		fntpMc->DumpData();


		//***Setup event shape object

	    TString PidSelection = "PidAlgoIdealCharged";//"PidAlgoMvd;PidAlgoStt;PidAlgoDrc";

		fAnalysis->FillList(all, "All", PidSelection);
		PndEventShape evsh(all, fini, 0.05, 0.1);

		//***Selection
	    fAnalysis->FillList(piminus, "PionAllMinus", PidSelection);
	    fAnalysis->FillList(NotCombinedPiPlus, "PionAllPlus", PidSelection);
	    fAnalysis->FillList(piplus, "PionAllPlus", PidSelection);
	    fAnalysis->FillList(proton, "ProtonAllPlus", PidSelection);
	    fAnalysis->FillList(antiProton, "ProtonAllMinus", PidSelection);
	    fAnalysis->FillList(kaonminus, "KaonAllMinus", PidSelection);


	    for (int pip=0; pip<piplus.GetLength(); ++pip){
	          fntpPiPlus->Column("ev",     (Float_t) fEvtCount);
	          fntpPiPlus->Column("cand",    (Float_t) pip);
	          fntpPiPlus->Column("ncand",   (Float_t) piplus.GetLength());
	          fntpPiPlus->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(piplus[pip]));

	          qa.qaP4("PiPlus_", piplus[pip]->P4(), fntpPiPlus);
	          qa.qaCand("PiPlus_", piplus[pip], fntpPiPlus);

	          numberOfHitsInSubdetector("PiPlus_", piplus[pip], fntpPiPlus);
	          tagNHits("PiPlus_", piplus[pip], fntpPiPlus);

	          fntpPiPlus->Column("PiPlus_CosTheta", (Float_t) piplus[pip]->GetMomentum().CosTheta());

	          RhoCandidate * mother_pip =0;
	          RhoCandidate * truth = piplus[pip]->GetMcTruth();


	          TLorentzVector l;
	          float costheta = -999.;
	          if(truth!=0x0){
	        	  l=truth->P4();
	        	  mother_pip=truth->TheMother();
	        	  costheta = truth->GetMomentum().CosTheta();
		          qa.qaCand("PiPlus_MC_", piplus[pip]->GetMcTruth(), fntpPiPlus);
	          }
	          else{
	        	  qa.qaCand("PiPlus_MC_", dummyCand, fntpPiPlus);
	          }

	          int moth_pip = (mother_pip==0x0)? 88888 : mother_pip->PdgCode();

			  fntpPiPlus->Column("Mother", (Float_t) moth_pip);


	          qa.qaP4("PiPlus_MC_", l, fntpPiPlus);
	          fntpPiPlus->Column("PiPlus_MC_CosTheta", (Float_t) costheta);


	          fntpPiPlus->DumpData();
	    }

	    for (int pim=0; pim<piminus.GetLength(); ++pim){
	        fntpPiMinus->Column("ev",     (Float_t) fEvtCount);
	        fntpPiMinus->Column("cand",    (Float_t) pim);
	        fntpPiMinus->Column("ncand",   (Float_t) piminus.GetLength());
	        fntpPiMinus->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(piminus[pim]));

	        qa.qaP4("piminus_", piminus[pim]->P4(), fntpPiMinus);
	        qa.qaCand("piminus_", piminus[pim], fntpPiMinus);

	        numberOfHitsInSubdetector("PiMinus_", piminus[pim], fntpPiMinus);
	        tagNHits("piminus_", piminus[pim], fntpPiMinus);

			RhoCandidate * mother_pim =0;
			RhoCandidate * truth = piminus[pim]->GetMcTruth();

			TLorentzVector l;
			float costheta = -999.;
			if(truth!=0x0){
			  l=truth->P4();
			  mother_pim = truth->TheMother();
			  costheta = truth->GetMomentum().CosTheta();
			  qa.qaCand("piminus_MC_", piminus[pim]->GetMcTruth(), fntpPiMinus);
			}
			else{
			  qa.qaCand("piminus_MC_", dummyCand, fntpPiMinus);
			}

			int moth_pim = (mother_pim==0x0)? 88888 : mother_pim->PdgCode();
			fntpPiMinus->Column("Mother", (Float_t) moth_pim);

			qa.qaP4("piminus_MC_", l, fntpPiMinus);
			fntpPiMinus->Column("piminus_MC_CosTheta", (Float_t) costheta);

	        fntpPiMinus->DumpData();
	    }

	    for (int prot=0; prot<proton.GetLength(); ++prot){
	        fntpProton->Column("ev",     (Float_t) fEvtCount);
	        fntpProton->Column("cand",    (Float_t) prot);
	        fntpProton->Column("ncand",   (Float_t) proton.GetLength());
	        fntpProton->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(proton[prot]));

	        qa.qaP4("proton_", proton[prot]->P4(), fntpProton);
	        qa.qaCand("proton_", proton[prot], fntpProton);

	        numberOfHitsInSubdetector("proton_", proton[prot], fntpProton);
	        tagNHits("proton_", proton[prot], fntpProton);

			RhoCandidate * mother_prot = 0;
			RhoCandidate * truth = proton[prot]->GetMcTruth();


			TLorentzVector l;
			float costheta = -999.;
			if(truth!=0x0){
			  l=truth->P4();
			  mother_prot = truth->TheMother();
			  costheta = truth->GetMomentum().CosTheta();
			  qa.qaCand("proton_MC_", proton[prot]->GetMcTruth(), fntpProton);
			}
			else{
			  qa.qaCand("proton_MC_", dummyCand, fntpProton);
			}

			int moth_pip = (mother_prot==0x0)? 88888 : mother_prot->PdgCode();
			fntpProton->Column("Mother", (Float_t) moth_pip);

			qa.qaP4("proton_MC_", l, fntpProton);
			fntpProton->Column("proton_MC_CosTheta", (Float_t) costheta);

	        fntpProton->DumpData();
	    }

	    for (int aProt=0; aProt<antiProton.GetLength(); ++aProt){
	        fntpAntiProton->Column("ev",     (Float_t) fEvtCount);
	        fntpAntiProton->Column("cand",    (Float_t) aProt);
	        fntpAntiProton->Column("ncand",   (Float_t) antiProton.GetLength());
	        fntpAntiProton->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(antiProton[aProt]));

	        qa.qaP4("AntiProton_", antiProton[aProt]->P4(), fntpAntiProton);
	        qa.qaCand("AntiProton_", antiProton[aProt], fntpAntiProton);

	        numberOfHitsInSubdetector("AntiProton_", antiProton[aProt], fntpAntiProton);
	        tagNHits("AntiProton_", antiProton[aProt], fntpAntiProton);

			RhoCandidate * mother_aprot = 0;
			RhoCandidate * truth = antiProton[aProt]->GetMcTruth();

			TLorentzVector l;
			float costheta = -999.;
			if(truth!=0x0){
			  l=truth->P4();
			  mother_aprot= truth->TheMother();
			  costheta = truth->GetMomentum().CosTheta();
			  qa.qaCand("AntiProton_MC_", antiProton[aProt]->GetMcTruth(), fntpAntiProton);
			}
			else{
			  qa.qaCand("AntiProton_MC_", dummyCand, fntpAntiProton);
			}

			int moth_pip = (mother_aprot==0x0)? 88888 : mother_aprot->PdgCode();
			fntpAntiProton->Column("Mother", (Float_t) moth_pip);

			qa.qaP4("AntiProton_MC_", l, fntpAntiProton);
			fntpAntiProton->Column("AntiProton_MC_CosTheta", (Float_t) costheta);

	        fntpAntiProton->DumpData();
	    }

	    for (int k=0; k<kaonminus.GetLength(); ++k){
			fntpKaonMinus->Column("ev",     (Float_t) fEvtCount);
			fntpKaonMinus->Column("cand",    (Float_t) k);
			fntpKaonMinus->Column("ncand",   (Float_t) kaonminus.GetLength());
			fntpKaonMinus->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(kaonminus[k]));

			qa.qaP4("kaonminus_", kaonminus[k]->P4(), fntpKaonMinus);
			qa.qaCand("kaonminus_", kaonminus[k], fntpKaonMinus);

			numberOfHitsInSubdetector("kaonminus_", kaonminus[k], fntpKaonMinus);
			tagNHits("kaonminus_", kaonminus[k], fntpKaonMinus);

			RhoCandidate * mother_k =0;
			RhoCandidate * truth = kaonminus[k]->GetMcTruth();


			TLorentzVector l;
			float costheta = -999.;
			if(truth!=0x0){
			  l=truth->P4();
			  mother_k =  truth->TheMother();
			  costheta = truth->GetMomentum().CosTheta();
			  qa.qaCand("kaonminus_MC_", kaonminus[k]->GetMcTruth(), fntpKaonMinus);
			}
			else{
			  qa.qaCand("kaonminus_MC_", dummyCand, fntpKaonMinus);
			}

			int moth_pip = (mother_k==0x0)? 88888 : mother_k->PdgCode();
			fntpKaonMinus->Column("Mother", (Float_t) moth_pip);

			qa.qaP4("kaonminus_MC_", l, fntpKaonMinus);
			fntpKaonMinus->Column("kaonminus_MC_CosTheta", (Float_t) costheta);

			fntpKaonMinus->DumpData();
		}


	    //***Lambda0 -> PiMinus + Proton

	     lambda0.Combine(piminus,proton);
	     lambda0.Select(lambdaMassSelector);
	     lambda0.SetType(kl0);

	     std::map<int,int> bestVtxFitLambda0, bestMassFitLambda0;

	     bestVtxFitLambda0 = VertexQaIndex(&lambda0);
	     bestMassFitLambda0 = MassFitQaIndex(&lambda0, fm0_lambda0);


	     for (int j=0; j<lambda0.GetLength(); ++j){


	       //general info about event
	       fntpLambda0->Column("ev",     (Float_t) fEvtCount);
	       fntpLambda0->Column("cand",    (Float_t) j);
	       fntpLambda0->Column("ncand",   (Float_t) lambda0.GetLength());
	       fntpLambda0->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(lambda0[j]));
	       fntpLambda0->Column("Lambda0_Pdg", (Float_t) lambda0[j]->PdgCode());

	 	  RhoCandidate * mother = 0;
	 	  RhoCandidate * truth = lambda0[j]->GetMcTruth();

	 	  if(truth) mother =  truth->TheMother();
	 	  int moth = (mother==0x0) ? 88888 : mother->PdgCode();

	 	  fntpLambda0->Column("Mother", (Float_t) moth);

	       qa.qaP4("Lambda0_", lambda0[j]->P4(), fntpLambda0);
	       qa.qaComp("Lambda0_", lambda0[j], fntpLambda0);

	       int tag = 0;
	       int ndau = lambda0[j]->NDaughters();
	       int dtag[2]={0,0};


	       for(int dau=0; dau<ndau; dau++){
	     	  RhoCandidate * daughter = lambda0[j]->Daughter(dau);
	     	  dtag[dau] = tagHits(daughter);
	       }

	       if(dtag[0]==1 && dtag[1]==1) tag=1;


	       fntpLambda0->Column("Lambda0_HitTag", (Int_t) tag);


	       // do vertex fit
	       PndKinVtxFitter vertexfitterLambda0 (lambda0[j]);
	       vertexfitterLambda0.Fit();
	       RhoCandidate * lambda0Fit = lambda0[j]->GetFit();


	       // store info of vertex fit
	       qa.qaFitter("VtxFit_", &vertexfitterLambda0, fntpLambda0);
	       fntpLambda0->Column("VtxFit_HowGood", (Int_t) bestVtxFitLambda0[j]);
	       qa.qaVtx("VtxFit_", lambda0Fit, fntpLambda0);
	       qa.qaComp("VtxFit_", lambda0Fit, fntpLambda0);
	       qa.qaP4("VtxFit_", lambda0Fit->P4(), fntpLambda0);

	       // differenz to MCTruth
	        qa.qaMcDiff("VtxFit_", lambda0Fit, fntpLambda0);
	        qaVtxDiff("VtxFit_", lambda0Fit, fntpLambda0);
	        qaMomRes("VtxFit_", lambda0Fit, fntpLambda0);


//		   // do Kalman vertex fit
//		   PndKalmanVtxFitter kalmanfitterLambda0 (lambda0[j]);
//		   kalmanfitterLambda0.Fit();
//		   RhoCandidate * lambda0KalmanFit = lambda0[j]->GetFit();
//
//
//		   // store info of vertex fit
//		   qa.qaFitter("KalmanFit_", &kalmanfitterLambda0, fntpLambda0);
//		   //fntpLambda0->Column("KalmanFit_HowGood", (Int_t) bestKalmanFitLambda0[j]);
//		   qa.qaVtx("KalmanFit_", lambda0KalmanFit, fntpLambda0);
//		   qa.qaComp("KalmanFit_", lambda0KalmanFit, fntpLambda0);
//		   qa.qaP4("KalmanFit_", lambda0KalmanFit->P4(), fntpLambda0);
//
//		   // differenz to MCTruth
//			qa.qaMcDiff("KalmanFit_", lambda0KalmanFit, fntpLambda0);
//			qaVtxDiff("KalmanFit_", lambda0KalmanFit, fntpLambda0);
//			qaMomRes("KalmanFit_", lambda0KalmanFit, fntpLambda0);

	       // do mass fit
	       PndKinFitter massFitterLambda0(lambda0Fit);
	       massFitterLambda0.AddMassConstraint(fm0_lambda0);
	       massFitterLambda0.Fit();

	       RhoCandidate * lambda0Fit_mass = lambda0Fit->GetFit();
	       qa.qaFitter("MassFit_", &massFitterLambda0, fntpLambda0);

	       fntpLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitLambda0[j]);



	       TLorentzVector l;


	 	    if(0x0 != truth){
	 	    	l = truth->P4();
	 	    	qa.qaVtx("McTruth_", truth, fntpLambda0);
	 	    }
	 	    else{
	 	    	qa.qaVtx("McTruth_", dummyCand, fntpLambda0);
	 	    }

	       qa.qaP4("McTruth_", l, fntpLambda0);


	       //*** use for Xi only bestChi2Cand which satisfied the mass fit

	       if (bestVtxFitLambda0[j]==1 && bestMassFitLambda0[j]>0 && tag==1){
	 		  Lambda0Fit.Append(lambda0Fit);
//	 		  CombinedList(lambda0Fit, &CombinedPiMinus, -211);
	       }


	       //***information of boosted particle
	       lambda0Fit->Boost(-beamBoost);
	       qa.qaComp("boost_", lambda0Fit, fntpLambda0);

	       fntpLambda0->DumpData();


	    }

//	     GetNotCombinedList(CombinedPiMinus, &NotCombinedPiMinus);
//	     CombinedPiMinus.Cleanup();



	     //***AntiLambda0 -> PiPlus + AntiProton
	      antiLambda0.Combine(piplus,antiProton);
	      antiLambda0.Select(lambdaMassSelector);
	      antiLambda0.SetType(kal0);

	      std::map<int,int> bestVtxFitAntiLambda0, bestMassFitAntiLambda0;

	      bestVtxFitAntiLambda0 = VertexQaIndex(&antiLambda0);
	      bestMassFitAntiLambda0 = MassFitQaIndex(&antiLambda0, fm0_lambda0);


	      for (int j=0; j<antiLambda0.GetLength(); ++j){

	        //general info about event
	        fntpAntiLambda0->Column("ev",     (Float_t) fEvtCount);
	        fntpAntiLambda0->Column("cand",    (Float_t) j);
	        fntpAntiLambda0->Column("ncand",   (Float_t) antiLambda0.GetLength());
	        fntpAntiLambda0->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(antiLambda0[j]));
	        fntpAntiLambda0->Column("AntiLambda0_Pdg", (Float_t) antiLambda0[j]->PdgCode());

			RhoCandidate * mother = 0;
			RhoCandidate * truth = antiLambda0[j]->GetMcTruth();

			if(truth) mother = truth->TheMother();

			int moth = (mother==0x0) ? 88888 : mother->PdgCode();
			fntpAntiLambda0->Column("Mother", (Float_t) moth);

	        qa.qaP4("AntiLambda0_", antiLambda0[j]->P4(), fntpAntiLambda0);
	        qa.qaComp("AntiLambda0_", antiLambda0[j], fntpAntiLambda0);


	        int tag = 0;
	        int ndau = antiLambda0[j]->NDaughters();
	        int dtag[2]={0,0};


	        for(int dau=0; dau<ndau; dau++){
	      	  	 RhoCandidate * daughter = antiLambda0[j]->Daughter(dau);
	      	  	 dtag[dau] = tagHits(daughter);
	        }

	        if(dtag[0]==1 && dtag[1]==1) tag=1;


	        fntpAntiLambda0->Column("AntiLambda0_HitTag", (Int_t) tag);



	        // do vertex fit
	        PndKinVtxFitter vertexfitterAntiLambda0 (antiLambda0[j]);
	        vertexfitterAntiLambda0.Fit();
	        RhoCandidate * antiLambda0Fit = antiLambda0[j]->GetFit();



	        // store info of vertex fit
	        qa.qaFitter("VtxFit_", &vertexfitterAntiLambda0, fntpAntiLambda0);
	        qa.qaVtx("VtxFit_", antiLambda0Fit, fntpAntiLambda0);
	        qa.qaP4("VtxFit_", antiLambda0Fit->P4(), fntpAntiLambda0);
	        qa.qaComp("VtxFit_",antiLambda0Fit, fntpAntiLambda0);

	        fntpAntiLambda0->Column("VtxFit_HowGood", (Int_t) bestVtxFitAntiLambda0[j]);

	        qa.qaMcDiff("VtxFit_", antiLambda0Fit, fntpAntiLambda0);
	        qaVtxDiff("VtxFit_", antiLambda0Fit, fntpAntiLambda0);
	        qaMomRes("VtxFit_", antiLambda0Fit, fntpAntiLambda0);


//	        // do kalman vertex fit
//	        PndKalmanVtxFitter kalmanfitterAntiLambda0 (antiLambda0[j]);
//	        kalmanfitterAntiLambda0.Fit();
//	        RhoCandidate * antiLambda0KalmanFit = antiLambda0[j]->GetFit();
//
//
//
//	        // store info of kalman fit
//	        qa.qaFitter("KalmanFit_", &kalmanfitterAntiLambda0, fntpAntiLambda0);
//	        qa.qaVtx("KalmanFit_", antiLambda0KalmanFit, fntpAntiLambda0);
//
//	        //fntpAntiLambda0->Column("KalmanFit_HowGood", (Int_t) bestKalmanFitAntiLambda0[j]);
//
//	        qa.qaMcDiff("KalmanFit_", antiLambda0KalmanFit, fntpAntiLambda0);
//	        qaVtxDiff("KalmanFit_", antiLambda0KalmanFit, fntpAntiLambda0);
//	        qaMomRes("KalmanFit_", antiLambda0KalmanFit, fntpAntiLambda0);



	        // do mass fit
	        PndKinFitter massFitterAntiLambda0(antiLambda0Fit);
	        massFitterAntiLambda0.AddMassConstraint(fm0_lambda0);
	        massFitterAntiLambda0.Fit();

	        RhoCandidate * antiLambda0Fit_mass = antiLambda0Fit->GetFit();
	        qa.qaFitter("MassFit_", &massFitterAntiLambda0, fntpAntiLambda0);
	        fntpAntiLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitAntiLambda0[j]);




	        TLorentzVector l;
	        if(0x0 != truth){
				l = truth->P4();
				qa.qaVtx("MCTruth_", truth, fntpAntiLambda0);
	        }
	        else{
	        	qa.qaVtx("MCTruth_", dummyCand, fntpAntiLambda0);
	        }

	        qa.qaP4("MCTruth_", l, fntpAntiLambda0);


	        //***information of boosted particle
	        antiLambda0Fit->Boost(-beamBoost);
	        qa.qaComp("boost_", antiLambda0Fit, fntpAntiLambda0);



	        if(bestVtxFitAntiLambda0[j]==1 && bestMassFitAntiLambda0[j]>0 && tag==1){
	  		  AntiLambda0Fit.Append(antiLambda0Fit);
	  		  CombinedList(antiLambda0Fit, &CombinedPiPlus, 211);
	        }

	        fntpAntiLambda0->DumpData();
	      }

	      GetNotCombinedList(CombinedPiPlus, &NotCombinedPiPlus);
	      CombinedPiPlus.Cleanup();



	      //*** Xi(1820)- -> Lambda0 + K-
	   	ximinus.Combine(Lambda0Fit, kaonminus);
	   	ximinus.Select(xi1820MassSelector);
	   	ximinus.SetType(kXim);

	   	std::map<int,int> BestVtxFitXiMinus, BestMassFitXiMinus;

	   	BestVtxFitXiMinus = VertexQaIndex(&ximinus);
	   	BestMassFitXiMinus = MassFitQaIndex(&ximinus, fm0_Xi1820);


		for (int j=0; j<ximinus.GetLength(); ++j){

			//general info about event
			fntpXiMinus1820->Column("ev",     (Float_t) fEvtCount);
			fntpXiMinus1820->Column("cand",    (Float_t) j);
			fntpXiMinus1820->Column("ncand",   (Float_t) ximinus.GetLength());
			fntpXiMinus1820->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(ximinus[j]));
			fntpXiMinus1820->Column("XiMinus_Pdg", (Float_t) ximinus[j]->PdgCode());


			RhoCandidate * mother = ximinus[j]->TheMother();

			int moth = (mother==0x0) ? 88888 : mother->PdgCode();
			fntpXiMinus1820->Column("Mother", (Float_t) moth);

			qa.qaP4("XiMinus_", ximinus[j]->P4(), fntpXiMinus1820);
			qa.qaComp("XiMinus_", ximinus[j], fntpXiMinus1820);
			qa.qaPoca("XiMinus_", ximinus[j], fntpXiMinus1820);



			// do vertex-fit

			PndKinVtxFitter vertexfitterXiMinus (ximinus[j]);
			vertexfitterXiMinus.Fit();
			RhoCandidate * ximinusFit = ximinus[j]->GetFit();


			// store info of vertex-fit

			qa.qaFitter("VtxFit_", &vertexfitterXiMinus, fntpXiMinus1820);
			fntpXiMinus1820->Column("VtxFit_HowGood", (Int_t) BestVtxFitXiMinus[j]);

			qa.qaVtx("VtxFit_", ximinusFit, fntpXiMinus1820);
			qa.qaP4("VtxFit_", ximinusFit->P4(), fntpXiMinus1820);
			qa.qaComp("VtxFit_", ximinusFit, fntpXiMinus1820);


			// difference to MCTruth
			qa.qaMcDiff("VtxFit_", ximinusFit, fntpXiMinus1820);
			qaVtxDiff("VtxFit_", ximinusFit, fntpXiMinus1820);
			qaMomRes("VtxFit_", ximinusFit, fntpXiMinus1820);


//			// do kalman vertex-fit
//
//			PndKalmanVtxFitter kalmanfitterXiMinus (ximinus[j]);
//			kalmanfitterXiMinus.Fit();
//			RhoCandidate * ximinusKalmanFit = ximinus[j]->GetFit();
//
//
//			// store info of Kalman-fit
//
//			qa.qaFitter("KalmanFit_", &kalmanfitterXiMinus, fntpXiMinus1820);
//			//fntpXiMinus1820->Column("KalmanFit_HowGood", (Int_t) BestKalmanFitXiMinus[j]);
//
//			qa.qaVtx("KalmanFit_", ximinusKalmanFit, fntpXiMinus1820);
//			qa.qaP4("KalmanFit_", ximinusKalmanFit->P4(), fntpXiMinus1820);
////			qa.qaCand("KalmanFit_", ximinusKalmanFit, fntpXiMinus1820);
//
//
//			// difference to MCTruth
//			qa.qaMcDiff("KalmanFit_", ximinusKalmanFit, fntpXiMinus1820);
//			qaVtxDiff("KalmanFit_", ximinusKalmanFit, fntpXiMinus1820);
//			qaMomRes("KalmanFit_", ximinusKalmanFit, fntpXiMinus1820);



			// do mass fit
			PndKinFitter massFitterXiMinus(ximinusFit);
			massFitterXiMinus.AddMassConstraint(fm0_Xi1820);
			massFitterXiMinus.Fit();

			RhoCandidate * ximinusFit_mass = ximinusFit->GetFit();
			qa.qaFitter("MassFit_", &massFitterXiMinus, fntpXiMinus1820);
			fntpXiMinus1820->Column("MassFit_HowGood", (Int_t) BestMassFitXiMinus[j]);

			qa.qaMcDiff("MassFit_", ximinusFit_mass, fntpXiMinus1820);
			qaVtxDiff("MassFit_", ximinusFit, fntpXiMinus1820);

			RhoCandidate * truth = ximinus[j]->GetMcTruth();
			TLorentzVector l;

			if(0x0 != truth){
				l = truth->P4();
				qa.qaVtx("MCTruth_", truth, fntpXiMinus1820);
			}
			else{
				qa.qaVtx("MCTruth_", dummyCand, fntpXiMinus1820);
			}

			qa.qaP4("MCTruth_", l, fntpXiMinus1820);


			if (BestVtxFitXiMinus[j]==1 && BestMassFitXiMinus[j]>0){
				XiMinusFit.Append(ximinusFit);
			}


			//***information of boosted particle
			ximinusFit->Boost(-beamBoost);
			qa.qaComp("boost_", ximinusFit, fntpXiMinus1820);

			fntpXiMinus1820->DumpData();


		}
		Lambda0Fit.Cleanup();
//		NotCombinedPiMinus.Cleanup();



		//*** Xi+ -> AntiLambda0 + Pi+
		xiplus.Combine(AntiLambda0Fit,NotCombinedPiPlus);
		xiplus.Select(xiMassSelector);
		xiplus.SetType(kaXip);

		std::map<int,int> BestVtxFitXiPlus, BestMassFitXiPlus;

		BestVtxFitXiPlus = VertexQaIndex(&xiplus);
		BestMassFitXiPlus = MassFitQaIndex(&xiplus, fm0_Xi);

	    for (int j=0; j<xiplus.GetLength(); ++j){

	      //general info about event
	      fntpXiPlus->Column("ev",     (Float_t) fEvtCount);
	      fntpXiPlus->Column("cand",    (Float_t) j);
	      fntpXiPlus->Column("ncand",   (Float_t) xiplus.GetLength());
	      fntpXiPlus->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(xiplus[j]));
	      fntpXiPlus->Column("XiPlus_Pdg", (Float_t) xiplus[j]->PdgCode());

	      RhoCandidate * mother = xiplus[j]->TheMother();
	      int moth = (mother==0x0) ? 88888 : mother->PdgCode();
	      fntpXiPlus->Column("Mother", (Float_t) moth);

	      qa.qaP4("Xiplus_", xiplus[j]->P4(), fntpXiPlus);
	      qa.qaComp("Xiplus_", xiplus[j], fntpXiPlus);



	      //******** do vertex-fit
	      PndKinVtxFitter vertexfitterxiplus (xiplus[j]);
	      vertexfitterxiplus.Fit();
	      RhoCandidate * xiplusFit = xiplus[j]->GetFit();


	      // store info of vertex-fit
	      qa.qaFitter("VtxFit_", &vertexfitterxiplus, fntpXiPlus);
	      fntpXiPlus->Column("VtxFit_HowGood", (Int_t) BestVtxFitXiPlus[j]);
	      qa.qaVtx("VtxFit_", xiplusFit, fntpXiPlus);
	      qa.qaP4("VtxFit_", xiplusFit->P4(), fntpXiPlus);
	      qa.qaComp("VtxFit_", xiplusFit, fntpXiPlus);

		  // difference to MCTruth
	      qa.qaMcDiff("VtxFit_", xiplusFit, fntpXiPlus);
	      qaVtxDiff("VtxFit_", xiplusFit, fntpXiPlus);
	      qaMomRes("VtxFit_", xiplusFit, fntpXiPlus);


//	      //******** do kalman vertex-fit
//	      PndKalmanVtxFitter kalmanfitterxiplus (xiplus[j]);
//	      kalmanfitterxiplus.Fit();
//	      RhoCandidate * xiplusKalmanFit = xiplus[j]->GetFit();
//
//
//	      // store info of kalman fit
//	      qa.qaFitter("KalmanFit_", &kalmanfitterxiplus, fntpXiPlus);
//	      //fntpXiPlus->Column("KalmanFit_HowGood", (Int_t) BestKalmanFitXiPlus[j]);
//	      qa.qaVtx("KalmanFit_", xiplusKalmanFit, fntpXiPlus);
//	      qa.qaP4("KalmanFit_", xiplusKalmanFit->P4(), fntpXiPlus);
//
//		  // difference to MCTruth
//	      qa.qaMcDiff("KalmanFit_", xiplusKalmanFit, fntpXiPlus);
//	      qaVtxDiff("KalmanFit_", xiplusKalmanFit, fntpXiPlus);
//	      qaMomRes("KalmanFit_", xiplusKalmanFit, fntpXiPlus);



	      //****** do mass fit
	      PndKinFitter massFitterxiplus(xiplusFit);
	      massFitterxiplus.AddMassConstraint(fm0_Xi);
	      massFitterxiplus.Fit();

	      RhoCandidate * xiplusFit_mass = xiplusFit->GetFit();
	      qa.qaFitter("MassFit_", &massFitterxiplus, fntpXiPlus);
	      fntpXiPlus->Column("MassFit_HowGood", (float) BestMassFitXiPlus[j]);
	      qa.qaVtx("MassFit_", xiplusFit_mass, fntpXiPlus);
	      qa.qaP4("MassFit_", xiplusFit_mass->P4(), fntpXiPlus);

	      qa.qaMcDiff("MassFit_", xiplusFit_mass, fntpXiPlus);
	      qaVtxDiff("MassFit_", xiplusFit_mass, fntpXiPlus);
	      qaMomRes("MassFit_", xiplusFit_mass, fntpXiPlus);

	      RhoCandidate * truth = xiplus[j]->GetMcTruth();
	      TLorentzVector l;
	      if(0x0 != truth){
					l = truth->P4();
					qa.qaVtx("MCTruth_", truth, fntpXiPlus);
	      }
	      else{
	    	  qa.qaVtx("MCTruth_", dummyCand, fntpXiPlus);
	      }

	      qa.qaP4("MCTruth_", l, fntpXiPlus);


	      if(BestVtxFitXiPlus[j]==1 && BestMassFitXiPlus[j]>0){
	    	  XiPlusFit.Append(xiplusFit);
	      }

	      //***information of boosted particle
	      xiplusFit->Boost(-beamBoost);
	      qa.qaComp("boost_", xiplusFit, fntpXiPlus);

	      fntpXiPlus->DumpData();
		 }

	    AntiLambda0Fit.Cleanup();
	    NotCombinedPiPlus.Cleanup();




//	    ******* Xi+ Xi- System*****************************

	     xiSys.Combine(XiPlusFit, XiMinusFit);
	     xiSys.SetType(88888);

	     for (int syscand=0; syscand<xiSys.GetLength(); ++syscand){

	 		fntpXiSys->Column("ev",     (Float_t) fEvtCount);
	 		fntpXiSys->Column("cand",    (Float_t) syscand);
	 		fntpXiSys->Column("ncand",   (Float_t) ximinus.GetLength());
	 		fntpXiSys->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(xiSys[syscand]));


	 		RhoCandidate * mother = xiSys[syscand]->TheMother();
	 		int moth = (mother==0x0) ? 88888 : mother->PdgCode();
	 		fntpXiSys->Column("Mother", (Float_t) moth);


	 		qa.qaP4("XiSys_", xiSys[syscand]->P4(), fntpXiSys);
	 		qa.qaComp("XiSys_", xiSys[syscand], fntpXiSys);
	  		qa.qaPoca("XiSys_", xiSys[syscand], fntpXiSys);


	 		RhoCandidate *  truth = xiSys[syscand]->GetMcTruth();
	 		TLorentzVector l;

	 		if (truth != 0x0){
	 			qa.qaVtx("McTruth_", truth, fntpXiSys);
	 			l = truth->P4();
	 		}
	 		else{
	 			qa.qaVtx("McTruth_", dummyCand, fntpXiSys);
	 		}
	 		qa.qaP4("McTruth_", l, fntpXiSys);

	 		//4C-Fitter

	 		PndKinFitter fitter4c (xiSys[syscand]);
	 		fitter4c.Add4MomConstraint(fini);
	 		fitter4c.Fit();

	 		RhoCandidate * xiSysFit4c = xiSys[syscand]->GetFit();

	 		qa.qaFitter("4CFit_", &fitter4c, fntpXiSys);
	 		qa.qaComp("4cFit_", xiSysFit4c, fntpXiSys);
	 		qa.qaVtx("4CFit_", xiSysFit4c, fntpXiSys);
	 		qa.qaMcDiff("4CFit_", xiSysFit4c, fntpXiSys);
	 		qaVtxDiff("4CFit_", xiSysFit4c, fntpXiSys);
	 		qaMomRes("4CFit_", xiSysFit4c, fntpXiSys);


	 		fntpXiSys->DumpData();


	     }
	     XiMinusFit.Cleanup();
	     XiPlusFit.Cleanup();

   }
}

void AnalysisTaskXi1820::Finish()
{

	//Write output
	foutput->cd();
	fntpMc -> GetInternalTree()->Write();
	fntpPiMinus ->GetInternalTree()->Write();
	fntpPiPlus->GetInternalTree()->Write();
	fntpProton->GetInternalTree()->Write();
	fntpAntiProton->GetInternalTree()->Write();
	fntpKaonMinus->GetInternalTree()->Write();
	fntpLambda0->GetInternalTree()->Write();
	fntpAntiLambda0->GetInternalTree()->Write();
	fntpXiMinus1820->GetInternalTree()->Write();
	fntpXiPlus->GetInternalTree()->Write();
	fntpXiSys->GetInternalTree()->Write();


	foutput->Save();

	ftimer.Stop();
	Double_t rtime = ftimer.RealTime();
	Double_t ctime = ftimer.CpuTime();
	cout<<endl<<endl;
	cout<<"Macro finisched successfully."<<endl;
	cout<<"Realtime: "<<rtime<<" s, CPU time: "<<ctime<<" s"<<endl;
}

