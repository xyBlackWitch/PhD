//class PndAnalysis;
//class PndAnaPidSelector;
//class RhoCandList;
//class RhoTuple;

#include "AnalysisTaskLambda0.h"

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
//#include "PndTrack.h"


using std::cout;
using std::endl;

ClassImp(AnalysisTaskLambda0)

//---------------Default constructor------------
AnalysisTaskLambda0::AnalysisTaskLambda0():
	FairTask("Panda Tutorial Analysis Task") {
	}
//------------------------------------------

//--------------- destructor ----------------
AnalysisTaskLambda0::~AnalysisTaskLambda0(){}
//-----------------------------------------


enum pidNumbers {
	kPip = 211, kPim = -211,
	kPp = 2212, kaPm = -2212,
	kl0 = 3122, kal0 = -3122,
	kKm = -321,
	kXim = 23314, kaXip = -3312
};

void AnalysisTaskLambda0::CombinedList(RhoCandidate* cand, RhoCandList* combinedList, int pdg){
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

void AnalysisTaskLambda0::GetNotCombinedList(RhoCandList combinedList, RhoCandList * candList){
	  for (int j=0; j<combinedList.GetLength(); j++){
		  RhoCandidate * combinedCand = combinedList[j];
		  candList->Remove(combinedCand);
	  }
}


void AnalysisTaskLambda0::numberOfHitsInSubdetector(TString pre, RhoCandidate *c, RhoTuple *n){

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

//void AnalysisTaskLambda0::tagNHits(TString pre, RhoCandidate *c, RhoTuple *n){
//
//	/**@brief Tag the particle with different integers
//	 * @details Tag the particle with different integers:
//	 * 0: if there is no hit in the detector
//	 * 1: sttHits>3 or mvdHits>3 or gemHit>3
//	 */
//
//	int tag=0;
//
//	PndPidCandidate * pidCand = (PndPidCandidate*)c->GetRecoCandidate();
//
//	if(pidCand){
//		int mvdHits = pidCand->GetMvdHits();
//		int sttHits = pidCand->GetSttHits();
//		int gemHits = pidCand->GetGemHits();
//
//		if(mvdHits>3 || sttHits>3 || gemHits>3) tag=1;
//
//
//	}
//
//	n->Column(pre + "HitTag", (Int_t) tag, 0);
//}

void AnalysisTaskLambda0::tagNHits(TString pre, RhoCandidate *c, RhoTuple *n){

	/**@brief Tag the particle with different integers
	 * @details Tag the particle with different integers:
	 * 0: if there is no hit in the detector
	 * 1: sttHits>3 or mvdHits>3 or gemHit>3
	 */

	int tag=0;

	PndPidCandidate * pidCand = (PndPidCandidate*)c->GetRecoCandidate();

	int branch = trackBranch(c);

	if(pidCand){
		int mvdHits = pidCand->GetMvdHits();
		int sttHits = pidCand->GetSttHits();
		int gemHits = pidCand->GetGemHits();


		if(mvdHits>3 || gemHits>3) tag=1;
		else if (sttHits>3 && branch==48) tag=1;
		else tag=0;
	}

	n->Column(pre + "HitTag", (Int_t) tag, 0);
}

int AnalysisTaskLambda0::tagHits(RhoCandidate *c){

	/**@brief Tag the particle with different integers
	 * @details Tag the particle with different integers:
	 * 0: if there is no hit in the detector
	 * 1: sttHits>3 or mvdHits>3 or gemHit>3
	 */
	int tag = 0;

	PndPidCandidate * pidCand = (PndPidCandidate*)c->GetRecoCandidate();

	int branch = trackBranch(c);

	if(pidCand){
		int mvdHits = pidCand->GetMvdHits();
		int sttHits = pidCand->GetSttHits();
		int gemHits = pidCand->GetGemHits();

		if(mvdHits>3 || gemHits>3) tag=1;
		else if (sttHits>3 && branch==48) tag=1;
		else tag=0;

	}

	return tag;
}

int AnalysisTaskLambda0::trackBranch(RhoCandidate *c){

	int branch=0;

	PndPidCandidate * pid = (PndPidCandidate*)c->GetRecoCandidate();
	if(pid){
		branch = pid->GetTrackBranch();
	}

	return branch;
}

void AnalysisTaskLambda0::TagTrackBranch(RhoCandidate *d0, RhoCandidate *d1, RhoTuple *n){
	/* @brief check if daughter particles cause no hit in the FTS
	 * @details check if daughter particles cause no hit in the FTS. 0 means cause a hit in FTS, 1 means cause no hin in FTS
	 */

	int tagbranch=0;

	PndPidCandidate * pidd0 = (PndPidCandidate*) d0->GetRecoCandidate();
	PndPidCandidate * pidd1 = (PndPidCandidate*) d1->GetRecoCandidate();

	if(pidd0 && pidd1){
		int branchd0 = pidd0->GetTrackBranch();
		int branchd1 = pidd1->GetTrackBranch();

		if(branchd0==48 & branchd1==48){
				tagbranch=1;
		}
	}
	 n->Column("NoFTSHit", (Int_t) tagbranch, -999);
}

void AnalysisTaskLambda0::TagTrackBranch(RhoCandidate *d0, RhoCandidate *d1, RhoCandidate *d2, RhoTuple *n){
	/* @brief check if daughter particles cause no hit in the FTS
	 * @details check if daughter particles cause no hit in the FTS. 0 means cause a hit in FTS, 1 means cause no hin in FTS
	 */

	int tagbranch=0;

	PndPidCandidate * pidd0 = (PndPidCandidate*) d0->GetRecoCandidate();
	PndPidCandidate * pidd1 = (PndPidCandidate*) d1->GetRecoCandidate();
	PndPidCandidate * pidd2 = (PndPidCandidate*) d2->GetRecoCandidate();

	if(pidd0 && pidd1 && pidd2){
		int branchd0 = pidd0->GetTrackBranch();
		int branchd1 = pidd1->GetTrackBranch();
		int branchd2 = pidd2->GetTrackBranch();

		if(branchd0==48 && branchd1==48 && branchd2==48){
				tagbranch=1;
		}
	}
	 n->Column("NoFTSHit", (Int_t) tagbranch, -999);
}


std::map<int,int> AnalysisTaskLambda0::VertexQaIndex(RhoCandList* candList, float probLimit=0.01){
	  /** @brief  give back the order of the best chi2
	   * @details give back the order of the best chi2!  1 means best, 2: second best (same with negative values for bad chi2 )
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

std::map<int,int> AnalysisTaskLambda0::MassFitQaIndex(RhoCandList* candList, float m0, float probLimit=0.01){
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

void AnalysisTaskLambda0::qaVtxDiff(TString pre, RhoCandidate * c, RhoTuple * n){

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

void AnalysisTaskLambda0::qaMomRes(TString pre, RhoCandidate * c, RhoTuple * n){

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

InitStatus AnalysisTaskLambda0::Init(){




  //*** create tuples
  fntpMC = new RhoTuple("ntpMC", "MCTruth info");
  fntpPiMinus = new RhoTuple("ntpPiMinus", "PiMinus info");
  fntpPiPlus = new RhoTuple("ntpPiPlus", "PiPlus info");
  fntpProton = new RhoTuple("ntpProton", "Proton info");
  fntpAntiProton = new RhoTuple("ntpAntiProton", "Antiproton info");
  fntpLambda0 = new RhoTuple("ntpLambda0", "Lambda0 info");
  fntpAntiLambda0 = new RhoTuple("ntpAntiLambda0", "AntiLambda0 info");
  fntpSys = new RhoTuple("ntpXiSys", "XiMinus XiPlus system info");

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


  fm0_beam = TDatabasePDG::Instance()->GetParticle("pbarpSystem")->Mass();
  cout<<"Mass pbar p system: "<<fm0_beam<<endl;


  //*** lorentz vector of the initial particle
  cout << "BeamMomentum: " << fmom << endl;
  double p_m0 = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
  fini.SetXYZT(0,0, fmom, sqrt(p_m0*p_m0+ fmom*fmom)+p_m0);

  if(fnevts==0) fnevts = fAnalysis->GetEntries();


  return kSUCCESS;

}

void AnalysisTaskLambda0::Exec(Option_t* op)
{
	 //RhoCandLists for analysis
	RhoCandList piplus, piminus, lambda0, antiLambda0, proton, antiProton, sys;
	RhoCandList CombinedPiPlus, NotCombinedPiPlus;
	RhoCandList Lambda0Fit, AntiLambda0Fit;
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
		qa.qaMcList("", mclist, fntpMC);
		fntpMC->DumpData();


		//***Setup event shape object

	    TString PidSelection = "PidAlgoMvd;PidAlgoStt;PidAlgoDrc";//;PidAlgoDisc;PidAlgoEmcBayes";//"PidAlgoIdealCharged";//

		fAnalysis->FillList(all, "All", PidSelection);
		PndEventShape evsh(all, fini, 0.05, 0.1);

		//***Selection
	    fAnalysis->FillList(piminus, "PionAllMinus", PidSelection);
	    fAnalysis->FillList(NotCombinedPiPlus, "PionAllPlus", PidSelection);
	    fAnalysis->FillList(piplus, "PionAllPlus", PidSelection);
	    fAnalysis->FillList(proton, "ProtonAllPlus", PidSelection);
	    fAnalysis->FillList(antiProton, "ProtonAllMinus", PidSelection);


	    for (int pip=0; pip<piplus.GetLength(); ++pip){
	        fntpPiPlus->Column("ev",     (Float_t) fEvtCount);
	        fntpPiPlus->Column("cand",    (Float_t) pip);
	        fntpPiPlus->Column("ncand",   (Float_t) piplus.GetLength());
	        fntpPiPlus->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(piplus[pip]));

	        fntpPiPlus->Column("TrackBranch", (Int_t) trackBranch(piplus[pip]));

	        qa.qaP4("piplus_", piplus[pip]->P4(), fntpPiPlus);
	        qa.qaCand("piplus_", piplus[pip], fntpPiPlus);

	        numberOfHitsInSubdetector("piplus_", piplus[pip], fntpPiPlus);
	        tagNHits("piplus_", piplus[pip], fntpPiPlus);

			RhoCandidate * mother_pip;
			RhoCandidate * truth = piplus[pip]->GetMcTruth();
			if (truth) mother_pip = truth->TheMother();

			int moth_pip;
			if (mother_pip==0x0){
			  moth_pip = 88888;
			}
			else
			  moth_pip = mother_pip->PdgCode();

			fntpPiPlus->Column("Mother", (Float_t) moth_pip);

			qa.qaMcDiff("piplus_", piplus[pip], fntpPiPlus);

			TLorentzVector l;
			float costheta = -999.;
			if(truth!=0x0){
			  l=truth->P4();
			  costheta = truth->GetMomentum().CosTheta();
			  qa.qaCand("piplus_MC_", piplus[pip]->GetMcTruth(), fntpPiPlus);
			}
			else{
			  qa.qaCand("piplus_MC_", dummyCand, fntpPiPlus);
			}

			qa.qaP4("piplus_MC_", l, fntpPiPlus);
			fntpPiPlus->Column("piplus_MC_CosTheta", (Float_t) costheta);

	        fntpPiPlus->DumpData();
	    }

	    for (int pim=0; pim<piminus.GetLength(); ++pim){
	        fntpPiMinus->Column("ev",     (Float_t) fEvtCount);
	        fntpPiMinus->Column("cand",    (Float_t) pim);
	        fntpPiMinus->Column("ncand",   (Float_t) piminus.GetLength());
	        fntpPiMinus->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(piminus[pim]));

	        int branch = trackBranch(piminus[pim]);
	        fntpPiMinus->Column("TrackBranch", (Int_t) branch);


	        qa.qaP4("piminus_", piminus[pim]->P4(), fntpPiMinus);
	        qa.qaCand("piminus_", piminus[pim], fntpPiMinus);

	        numberOfHitsInSubdetector("PiMinus_", piminus[pim], fntpPiMinus);
	        tagNHits("piminus_", piminus[pim], fntpPiMinus);

			RhoCandidate * mother_pim;
			RhoCandidate * truth = piminus[pim]->GetMcTruth();
			if (truth) mother_pim = truth->TheMother();

			int moth_pim;
			if (mother_pim==0x0){
			  moth_pim = 88888;
			}
			else
			  moth_pim = mother_pim->PdgCode();

			fntpPiMinus->Column("Mother", (Float_t) moth_pim);

			qa.qaMcDiff("piminus_", piminus[pim], fntpPiMinus);

			TLorentzVector l;
			float costheta = -999.;
			if(truth!=0x0){
			  l=truth->P4();
			  costheta = truth->GetMomentum().CosTheta();
			  qa.qaCand("piminus_MC_", piminus[pim]->GetMcTruth(), fntpPiMinus);
			}
			else{
			  qa.qaCand("piminus_MC_", dummyCand, fntpPiMinus);
			}

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

			RhoCandidate * truth = proton[prot]->GetMcTruth();

			RhoCandidate * mother_prot;
			if (truth) mother_prot = truth->TheMother();
			int moth_prot;
			if (mother_prot==0x0){
			  moth_prot = 88888;
			}
			else
			  moth_prot = mother_prot->PdgCode();

			fntpProton->Column("Mother", (Float_t) moth_prot);

			TLorentzVector l;
			float costheta = -999.;
			if(truth!=0x0){
			  l=truth->P4();
			  costheta = truth->GetMomentum().CosTheta();
			  qa.qaCand("proton_MC_", proton[prot]->GetMcTruth(), fntpProton);
			}
			else{
			  qa.qaCand("proton_MC_", dummyCand, fntpProton);
			}

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

			RhoCandidate * mother_aprot;
			RhoCandidate * truth = antiProton[aProt]->GetMcTruth();
			if (truth) mother_aprot = truth->TheMother();

			int moth_aprot;
			if (mother_aprot==0x0){
			  moth_aprot = 88888;
			}
			else
			  moth_aprot = mother_aprot->PdgCode();

			fntpAntiProton->Column("Mother", (Float_t) moth_aprot);


			TLorentzVector l;
			float costheta = -999.;
			if(truth!=0x0){
			  l=truth->P4();
			  costheta = truth->GetMomentum().CosTheta();
			  qa.qaCand("AntiProton_MC_", antiProton[aProt]->GetMcTruth(), fntpAntiProton);
			}
			else{
			  qa.qaCand("AntiProton_MC_", dummyCand, fntpAntiProton);
			}

			qa.qaP4("AntiProton_MC_", l, fntpAntiProton);
			fntpAntiProton->Column("AntiProton_MC_CosTheta", (Float_t) costheta);

	        fntpAntiProton->DumpData();
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

	       TagTrackBranch(lambda0[j]->Daughter(0), lambda0[j]->Daughter(1), fntpLambda0);


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


	       fntpLambda0->Column("HitTag", (Int_t) tag);



	       qa.qaMcDiff("Lambda0_", lambda0[j], fntpLambda0);
		   qaVtxDiff("Lambda0_", lambda0[j], fntpLambda0);
		   qa.qaMcDiff("Lambda0d0_", lambda0[j]->Daughter(0), fntpLambda0);
		   qa.qaMcDiff("Lambda0d1_", lambda0[j]->Daughter(1), fntpLambda0);

//		   fAnalysis->ResetDaughters(lambda0[j]);


	       // do vertex fit
	       PndKinVtxFitter vertexfitterLambda0 (lambda0[j]);
	       vertexfitterLambda0.Fit();
	       RhoCandidate * lambda0Fit = lambda0[j]->GetFit();

	       RhoCandidate * mother = lambda0Fit->TheMother();
	       int moth = (mother==0x0) ? 88888 : mother->PdgCode();
	       fntpLambda0->Column("Mother", (Float_t) moth);


	       // store info of vertex fit
	       qa.qaFitter("VtxFit_", &vertexfitterLambda0, fntpLambda0);
	       fntpLambda0->Column("VtxFit_HowGood", (Int_t) bestVtxFitLambda0[j]);
	       qa.qaVtx("VtxFit_", lambda0Fit, fntpLambda0);
	       qa.qaComp("VtxFit_", lambda0Fit, fntpLambda0);
	       qa.qaP4("VtxFit_", lambda0Fit->P4(), fntpLambda0);
	       qa.qaMcDiff("VtxFitd0_", lambda0Fit->Daughter(0), fntpLambda0);
	       qa.qaMcDiff("VtxFitd1_", lambda0Fit->Daughter(1), fntpLambda0);



	       // differenz to MCTruth
	        qa.qaMcDiff("VtxFit_", lambda0Fit, fntpLambda0);
	        qaVtxDiff("VtxFit_", lambda0Fit, fntpLambda0);
	        qaMomRes("VtxFit_", lambda0Fit, fntpLambda0);


		   // do Kalman vertex fit
		   PndKalmanVtxFitter kalmanfitterLambda0 (lambda0[j]);
		   kalmanfitterLambda0.Fit();
		   RhoCandidate * lambda0KalmanFit = lambda0[j]->GetFit();


		   // store info of vertex fit
		   qa.qaFitter("KalmanFit_", &kalmanfitterLambda0, fntpLambda0);
		   //fntpLambda0->Column("KalmanFit_HowGood", (Int_t) bestKalmanFitLambda0[j]);
		   qa.qaVtx("KalmanFit_", lambda0KalmanFit, fntpLambda0);
		   qa.qaComp("KalmanFit_", lambda0KalmanFit, fntpLambda0);
		   qa.qaP4("KalmanFit_", lambda0KalmanFit->P4(), fntpLambda0);

		   // differenz to MCTruth
			qa.qaMcDiff("KalmanFit_", lambda0KalmanFit, fntpLambda0);
			qaVtxDiff("KalmanFit_", lambda0KalmanFit, fntpLambda0);
			qaMomRes("KalmanFit_", lambda0KalmanFit, fntpLambda0);

	       // do mass fit
	       PndKinFitter massFitterLambda0(lambda0Fit);
	       massFitterLambda0.AddMassConstraint(fm0_lambda0);
	       massFitterLambda0.Fit();

	       RhoCandidate * lambda0Fit_mass = lambda0Fit->GetFit();
	       qa.qaFitter("MassFit_", &massFitterLambda0, fntpLambda0);

	       fntpLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitLambda0[j]);


	       RhoCandidate * truth = lambda0[j]->GetMcTruth();
	       RhoCandidate * truthDaughter = lambda0[j]->Daughter(0)->GetMcTruth();
	       TLorentzVector l;
	       TVector3 dl;

	 	    if(0x0 != truth){
	 	    	l = truth->P4();
	 	    	qa.qaVtx("McTruth_", truth, fntpLambda0);
	 	    	dl = truth->Daughter(0)->Pos();
	 	    }
	 	    else{
	 	    	qa.qaVtx("McTruth_", dummyCand, fntpLambda0);
	 	    }

	       qa.qaP4("McTruth_", l, fntpLambda0);


	       //*** use for Xi only bestChi2Cand

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

	     std::map<int,int> bestVtxFitantiLambda0, bestMassFitantiLambda0;

	     bestVtxFitantiLambda0 = VertexQaIndex(&antiLambda0);
	     bestMassFitantiLambda0 = MassFitQaIndex(&antiLambda0, fm0_lambda0);


	     for (int j=0; j<antiLambda0.GetLength(); ++j){


	       //general info about event
	       fntpAntiLambda0->Column("ev",     (Float_t) fEvtCount);
	       fntpAntiLambda0->Column("cand",    (Float_t) j);
	       fntpAntiLambda0->Column("ncand",   (Float_t) antiLambda0.GetLength());
	       fntpAntiLambda0->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(antiLambda0[j]));
	       fntpAntiLambda0->Column("antiLambda0_Pdg", (Float_t) antiLambda0[j]->PdgCode());



	       qa.qaP4("antiLambda0_", antiLambda0[j]->P4(), fntpAntiLambda0);
	       qa.qaComp("antiLambda0_", antiLambda0[j], fntpAntiLambda0);

	       int tag = 0;
	       int ndau = antiLambda0[j]->NDaughters();
	       int dtag[2]={0,0};


	       for(int dau=0; dau<ndau; dau++){
	     	  RhoCandidate * daughter = antiLambda0[j]->Daughter(dau);
	     	  dtag[dau] = tagHits(daughter);
	       }

	       if(dtag[0]==1 && dtag[1]==1) tag=1;


	       fntpAntiLambda0->Column("HitTag", (Int_t) tag);

	       qa.qaMcDiff("antiLambda0_", antiLambda0[j], fntpAntiLambda0);
		   qaVtxDiff("antiLambda0_", antiLambda0[j], fntpAntiLambda0);


	       // do vertex fit
	       PndKinVtxFitter vertexfitterantiLambda0 (antiLambda0[j]);
	       vertexfitterantiLambda0.Fit();
	       RhoCandidate * antiLambda0Fit = antiLambda0[j]->GetFit();

	       RhoCandidate * mother = antiLambda0Fit->TheMother();
	       int moth = (mother==0x0) ? 88888 : mother->PdgCode();

	       fntpAntiLambda0->Column("Mother", (Float_t) moth);

	       // store info of vertex fit
	       qa.qaFitter("VtxFit_", &vertexfitterantiLambda0, fntpAntiLambda0);
	       fntpAntiLambda0->Column("VtxFit_HowGood", (Int_t) bestVtxFitantiLambda0[j]);
	       qa.qaVtx("VtxFit_", antiLambda0Fit, fntpAntiLambda0);
	       qa.qaComp("VtxFit_", antiLambda0Fit, fntpAntiLambda0);
	       qa.qaP4("VtxFit_", antiLambda0Fit->P4(), fntpAntiLambda0);

	       // differenz to MCTruth
	        qa.qaMcDiff("VtxFit_", antiLambda0Fit, fntpAntiLambda0);
	        qaVtxDiff("VtxFit_", antiLambda0Fit, fntpAntiLambda0);
	        qaMomRes("VtxFit_", antiLambda0Fit, fntpAntiLambda0);



		   // do Kalman vertex fit
		   PndKalmanVtxFitter kalmanfitterantiLambda0 (antiLambda0[j]);
		   kalmanfitterantiLambda0.Fit();
		   RhoCandidate * antiLambda0KalmanFit = antiLambda0[j]->GetFit();


		   // store info of vertex fit
		   qa.qaFitter("KalmanFit_", &kalmanfitterantiLambda0, fntpAntiLambda0);
		   //fntpAntiLambda0->Column("KalmanFit_HowGood", (Int_t) bestKalmanFitantiLambda0[j]);
		   qa.qaVtx("KalmanFit_", antiLambda0KalmanFit, fntpAntiLambda0);
		   qa.qaComp("KalmanFit_", antiLambda0KalmanFit, fntpAntiLambda0);
		   qa.qaP4("KalmanFit_", antiLambda0KalmanFit->P4(), fntpAntiLambda0);

		   // differenz to MCTruth
			qa.qaMcDiff("KalmanFit_", antiLambda0KalmanFit, fntpAntiLambda0);
			qaVtxDiff("KalmanFit_", antiLambda0KalmanFit, fntpAntiLambda0);
			qaMomRes("KalmanFit_", antiLambda0KalmanFit, fntpAntiLambda0);

	       // do mass fit
	       PndKinFitter massFitterantiLambda0(antiLambda0Fit);
	       massFitterantiLambda0.AddMassConstraint(fm0_lambda0);
	       massFitterantiLambda0.Fit();

	       RhoCandidate * antiLambda0Fit_mass = antiLambda0Fit->GetFit();
	       qa.qaFitter("MassFit_", &massFitterantiLambda0, fntpAntiLambda0);

	       fntpAntiLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitantiLambda0[j]);


	       RhoCandidate * truth = antiLambda0[j]->GetMcTruth();
	       RhoCandidate * truthDaughter = antiLambda0[j]->Daughter(0)->GetMcTruth();
	       TLorentzVector l;
	       TVector3 dl;

	 	    if(0x0 != truth){
	 	    	l = truth->P4();
	 	    	qa.qaVtx("McTruth_", truth, fntpAntiLambda0);
	 	    	dl = truth->Daughter(0)->Pos();
	 	    }
	 	    else{
	 	    	qa.qaVtx("McTruth_", dummyCand, fntpAntiLambda0);
	 	    }

	       qa.qaP4("McTruth_", l, fntpAntiLambda0);


	       //*** use for Xi only bestChi2Cand

	       if (bestVtxFitantiLambda0[j]==1 && bestMassFitantiLambda0[j]>0 && tag==1){
	 		  AntiLambda0Fit.Append(antiLambda0Fit);
	 		  CombinedList(antiLambda0Fit, &CombinedPiPlus, 211);
	       }


	       //***information of boosted particle
	       antiLambda0Fit->Boost(-beamBoost);
	       qa.qaComp("boost_", antiLambda0Fit, fntpAntiLambda0);

	       fntpAntiLambda0->DumpData();


	      }

	      GetNotCombinedList(CombinedPiPlus, &NotCombinedPiPlus);
	      CombinedPiPlus.Cleanup();






	    //*** pbar + p -> Lambda0 + AntiLambda0

	 	sys.Combine(Lambda0Fit, AntiLambda0Fit);
	 	sys.SetType(88888);


	 	for (int j=0; j<sys.GetLength(); ++j){

	 		//general information about event
	 		fntpSys->Column("ev",     (Float_t) fEvtCount);
	 		fntpSys->Column("cand",    (Float_t) j);
	 		fntpSys->Column("ncand",   (Float_t) sys.GetLength());
	 		fntpSys->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(sys[j]));

	 		qa.qaP4("", sys[j]->P4(), fntpSys);
	 		qa.qaComp("", sys[j], fntpSys);
	 		qa.qaPoca("", sys[j], fntpSys);
	 		qa.qaEventShapeShort("es_", &evsh, fntpAntiLambda0);


	 		//***do 4c fit
	 		PndKinFitter cc_Fitter4c (sys[j]);
	 		cc_Fitter4c.Add4MomConstraint(fini);
	 		cc_Fitter4c.Fit();

	 		RhoCandidate * ccFit = sys[j]->GetFit();

	 		// store info of 4c Fit
	 		fntpSys->Column("f4c_Chi2_cc", (Float_t) cc_Fitter4c.GetChi2());
	 		fntpSys->Column("f4c_NDF_cc", (Float_t) cc_Fitter4c.GetNdf());
	 		fntpSys->Column("f4c_Prob_cc", (Float_t) cc_Fitter4c.GetProb());

	 		qa.qaComp("f4c_", ccFit, fntpSys);

	 		//difference to MC Truth
	 		qa.qaMcDiff("f4c_", ccFit, fntpSys);
	 		qaVtxDiff("f4c_", ccFit, fntpSys);
	 		qaMomRes("f4c_", ccFit, fntpSys);


	 		fntpSys->DumpData();
	 	}
	 	 Lambda0Fit.Cleanup();
	 	 AntiLambda0Fit.Cleanup();

   }
}

void AnalysisTaskLambda0::Finish()
{

	//Write output
	foutput->cd();
	fntpMC -> GetInternalTree()->Write();
	fntpPiMinus ->GetInternalTree()->Write();
	fntpPiPlus->GetInternalTree()->Write();
	fntpProton->GetInternalTree()->Write();
	fntpAntiProton->GetInternalTree()->Write();
	fntpLambda0->GetInternalTree()->Write();
	fntpAntiLambda0->GetInternalTree()->Write();
	fntpSys->GetInternalTree()->Write();


	foutput->Save();

	ftimer.Stop();
	Double_t rtime = ftimer.RealTime();
	Double_t ctime = ftimer.CpuTime();
	cout<<endl<<endl;
	cout<<"Macro finisched successfully."<<endl;
	cout<<"Realtime: "<<rtime<<" s, CPU time: "<<ctime<<" s"<<endl;
	cout << "AnalysisTaskLambda0 finished" << endl;
}

