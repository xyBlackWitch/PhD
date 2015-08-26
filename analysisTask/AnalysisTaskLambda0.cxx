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
	kXim = 3312, kaXip = -3312
};


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

void AnalysisTaskLambda0::tagNHits(TString pre, RhoCandidate *c, RhoTuple *n){

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

int AnalysisTaskLambda0::tagHits(RhoCandidate *c){

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



std::map<int,int> AnalysisTaskLambda0::VertexQaIndex(RhoCandList* candList, float probLimit=0.01){
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
	fntpMC = new RhoTuple("fntpMC", "MCTruth info");
	fntpPiMinus = new RhoTuple("fntpPiMinus", "PiMinus info");
	fntpPiPlus = new RhoTuple("fntpPiPlus", "PiPlus info");
	fntpProton = new RhoTuple("fntpProton", "Proton info");
	fntpAntiProton = new RhoTuple("fntpAntiProton", "Antiproton info");
	fntpLambda0 = new RhoTuple("fntpLambda0", "Lambda0 info");
	fntpAntiLambda0 = new RhoTuple("fntpAntiLambda0", "AntiLambda0 info");
	fntpCrossCheck = new RhoTuple("fntpCrossCheck", "CrossCheck info");

	//Create output file
	TString outpath;
	outpath.Append(fOutPath);
	foutput = TFile::Open(outpath+"output_ana.root","RECREATE");

	// data reader Object
	fAnalysis = new PndAnalysis();
	if (fnevts==0) fnevts = fAnalysis->GetEntries();

	//***Mass selector
	fm0_lambda0 = TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
	cout<<"Mass of Lambda0: "<<fm0_lambda0<<endl;
	lambdaMassSelector = new RhoMassParticleSelector("lambda0", fm0_lambda0, 0.3);

	fm0_beam= TDatabasePDG::Instance()->GetParticle("pbarpSystem")->Mass();

	double p_m0 = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
	fini.SetXYZT(0,0, fmom, sqrt(p_m0*p_m0+ fmom*fmom)+p_m0);



	return kSUCCESS;

}

void AnalysisTaskLambda0::Exec(Option_t* op)
{
	TVector3 beamBoost = fini.BoostVector();

	//RhoCandLists for analysis
	RhoCandList piplus, piminus, lambda0, antiLambda0, proton, antiProton, crossCheck, Lambda0Fit, AntiLambda0Fit;
	RhoCandList mclist, all;

	RhoCandidate * dummyCand = new RhoCandidate();

	PndRhoTupleQA qa(fAnalysis, fmom);

	int evt=-1;
	while (fAnalysis->GetEvent() && ++evt<fnevts){

		if ((evt%100)==0) cout << "evt "<< evt <<endl;

		//***get MC list and store info
		fAnalysis->FillList(mclist, "McTruth");
		qa.qaMcList("", mclist, fntpMC);
		fntpMC->DumpData();


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

		fAnalysis->FillList(all, "All", PidSelection);
		PndEventShape evsh(all, fini, 0.05, 0.1);

		//***Selection with no PID info
		fAnalysis->FillList(piminus, "PionAllMinus", PidSelection);
		fAnalysis->FillList(piplus, "PionAllPlus", PidSelection);
		fAnalysis->FillList(proton, "ProtonAllPlus", PidSelection);
		fAnalysis->FillList(antiProton, "ProtonAllMinus", PidSelection);

	    //Get piminus
	    for (int j=0; j<piminus.GetLength(); ++j){

			//general info about event
			fntpPiMinus->Column("ev",     (Float_t) evt);
			fntpPiMinus->Column("cand",    (Float_t) j);
			fntpPiMinus->Column("ncand",   (Float_t) piminus.GetLength());


			//info about 4-vector
			qa.qaP4("PiMinus_",piminus[j]->P4(), fntpPiMinus);
			qa.qaP4("PiMinus_MC_",  piminus[j]->GetMcTruth()->P4(), fntpPiMinus);
			qa.qaCand("PiMinus_", piminus[j], fntpPiMinus);


			RhoCandidate * truth = piminus[j]->GetMcTruth();
			RhoCandidate * mother;
			float costht = -999.;

			if (truth!=0x0){
				mother = truth->TheMother();
				costht = truth->GetMomentum().CosTheta();
			}

			int moth = (mother==0x0) ? 88888 : mother->PdgCode();

			fntpPiMinus->Column("MCTruthMatch", (bool) fAnalysis->McTruthMatch(piminus[j]));
			fntpPiMinus->Column("Mother", (Int_t) moth);
			fntpPiMinus->Column("PiMinus_CosTheta", (Float_t) piminus[j]->GetMomentum().CosTheta());
			fntpPiMinus->Column("PiMinus_MC_CosTheta", (Float_t) costht);


			fntpPiMinus->DumpData();
	    }

		//Get PiPlus
		for (int j=0; j<piplus.GetLength(); ++j){
			qa.qaP4("PiPlus_", piplus[j]->P4(), fntpPiPlus);
			qa.qaCand("PiPlus_", piplus[j], fntpPiPlus);

			//general info about event
			fntpPiPlus->Column("ev",     (Float_t) evt);
			fntpPiPlus->Column("cand",    (Float_t) j);
			fntpPiPlus->Column("ncand",   (Float_t) piplus.GetLength());


			//info about 4-vector
			qa.qaP4("PiPlus_", piplus[j]->P4(), fntpPiPlus);
			qa.qaP4("PiPlus_MC_", piplus[j]->GetMcTruth()->P4(), fntpPiPlus);

			RhoCandidate * truth = piplus[j]->GetMcTruth();
			RhoCandidate * mother;
			float tht = -999.;

			if (truth!=0x0){
				mother=truth->TheMother();
				tht = truth->GetMomentum().CosTheta();
			}

			int moth = (mother==0x0) ? 88888 : mother->PdgCode();

			fntpPiPlus->Column("MCTruthMatch", (bool) fAnalysis->McTruthMatch(piplus[j]));
			fntpPiPlus->Column("Mother", (Int_t) moth);
			fntpPiPlus->Column("PiPlus_CosTheta", (Float_t) piplus[j]->GetMomentum().CosTheta());
			fntpPiPlus->Column("PiPlus_MC_CosTheta", (Float_t) tht);

			fntpPiPlus->DumpData();

		}

		//Get Proton
		for (int j=0; j<proton.GetLength(); j++){

			//general info about event
			fntpProton->Column("ev",     (Float_t) evt);
			fntpProton->Column("cand",    (Float_t) j);
			fntpProton->Column("ncand",   (Float_t) proton.GetLength());



			//info about 4-vector
			qa.qaP4("Proton_", proton[j]->P4(), fntpProton);
			qa.qaP4("Proton_MC_", proton[j]->GetMcTruth()->P4(), fntpProton);
			qa.qaCand("Proton_", proton[j], fntpProton);
			qa.qaVtx("Proton_", proton[j], fntpProton);

			RhoCandidate * truth  = proton[j]->GetMcTruth();
			RhoCandidate * mother;
			float tht = -999.;

			if (truth != 0x0){
				mother = truth->TheMother();
				tht = truth->GetMomentum().CosTheta();
			}

			int moth = (mother==0x0) ? 88888 : mother->PdgCode();

			fntpProton->Column("MCTruthMatch", (bool) fAnalysis->McTruthMatch(proton[j]));
			fntpProton->Column("Mother", (Int_t) moth);
			fntpProton->Column("Proton_CosTheta", (Float_t) proton[j]->GetMomentum().CosTheta());
			fntpProton->Column("Proton_MC_CosTheta", (Float_t) tht);



			fntpProton->DumpData();
		}

		//Get Antiproton
		for (int j=0; j<antiProton.GetLength(); j++){

			//general info about event
			fntpAntiProton->Column("ev",     (Float_t) evt);
			fntpAntiProton->Column("cand",    (Float_t) j);
			fntpAntiProton->Column("ncand",   (Float_t) antiProton.GetLength());

			//info about 4-vector
			qa.qaP4("AntiProton_", antiProton[j]->P4(), fntpAntiProton);
			qa.qaP4("AntiProton_MC_", antiProton[j]->GetMcTruth()->P4(), fntpAntiProton);
			qa.qaCand("AntiProton_", antiProton[j], fntpAntiProton);


			RhoCandidate * mother = antiProton[j]->GetMcTruth()->TheMother();
			int moth = (mother==0x0) ? 88888 : mother->PdgCode();

			fntpAntiProton->Column("MCTruthMatch", (bool) fAnalysis->McTruthMatch(antiProton[j]));
			fntpAntiProton->Column("Mother", (Int_t) moth);
			fntpAntiProton->Column("AntiProton_CosTheta", (Float_t) antiProton[j]->GetMomentum().CosTheta());
			fntpAntiProton->Column("AntiProton_MC_CosTheta", (Float_t) antiProton[j]->GetMcTruth()->GetMomentum().CosTheta());

			fntpAntiProton->DumpData();
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
			fntpLambda0->Column("ev",     (Float_t) evt);
			fntpLambda0->Column("cand",    (Float_t) j);
			fntpLambda0->Column("ncand",   (Float_t) lambda0.GetLength());
			fntpLambda0->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(lambda0[j]));


			RhoCandidate * mother = lambda0[j]->TheMother();
			int moth = (mother==0x0) ? 88888 : mother->PdgCode();

			fntpLambda0->Column("Mother", (Float_t) moth);

			qa.qaP4("Lambda0_", lambda0[j]->P4(), fntpLambda0);
			qa.qaCand("Lambda0_", lambda0[j], fntpLambda0);
			qa.qaComp("Lambda0_", lambda0[j], fntpLambda0);
			qa.qaMcDiff("Lambda0_", lambda0[j], fntpLambda0);
			//      qa.qaEventShapeShort("es_", &evsh, fntpLambda0);



			// do vertex fit

			PndKinVtxFitter vertexfitterLambda0 (lambda0[j]);
			vertexfitterLambda0.Fit();
			RhoCandidate * lambda0Fit = lambda0[j]->GetFit();

			qa.qaFitter("VtxFit_", &vertexfitterLambda0, fntpLambda0);
			//      fntpLambda0->Column("VtxFit_HowGood", (Int_t)  bestVtxFit[j]);
			qa.qaVtx("VtxFit_", lambda0Fit , fntpLambda0);
			qa.qaCand("VtxFit_", lambda0Fit, fntpLambda0);
			qa.qaP4("VtxFit_", lambda0Fit->P4(), fntpLambda0);
			qa.qaComp("VtxFit_", lambda0Fit, fntpLambda0);

			// differenz to MCTruth
			qa.qaMcDiff("VtxFit_", lambda0Fit, fntpLambda0);
			qaVtxDiff("VtxFit_", lambda0Fit, fntpLambda0);
			qaMomRes("VtxFit_", lambda0Fit, fntpLambda0);


			PndKalmanVtxFitter kalmanFitterLambda0 (lambda0[j]);
			kalmanFitterLambda0.Fit();
			RhoCandidate * kalmanFitLambda0 = lambda0[j]->GetFit();

			qa.qaFitter("KalmanFit_", &kalmanFitterLambda0, fntpLambda0);
			//      fntpLambda0->Column("KalmanFit_HowGood", (Int_t)  bestKalmanFitLambda0[j]);
			qa.qaVtx("KalmanFit_", kalmanFitLambda0 , fntpLambda0);
			qa.qaCand("KalmanFit_", kalmanFitLambda0, fntpLambda0);
			qa.qaP4("KalmanFit_", kalmanFitLambda0->P4(), fntpLambda0);
			qa.qaComp("KalmanFit_", kalmanFitLambda0, fntpLambda0);

			// differenz to MCTruth
			qa.qaMcDiff("KalmanFit_", kalmanFitLambda0, fntpLambda0);
			qaVtxDiff("KalmanFit_", kalmanFitLambda0, fntpLambda0);
			qaMomRes("KalmanFit_", kalmanFitLambda0, fntpLambda0);



			// do mass fit
			PndKinFitter massFitterLambda0(lambda0Fit);
			massFitterLambda0.AddMassConstraint(fm0_lambda0);
			massFitterLambda0.Fit();

			RhoCandidate * lambda0Fit_mass = lambda0Fit->GetFit();

			qa.qaFitter("MassFit_", &massFitterLambda0, fntpLambda0);
			//      fntpLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitLambda0[j]);
			qa.qaMcDiff("MassFit_", lambda0Fit_mass, fntpLambda0);
			qaVtxDiff("VtxFit_", lambda0Fit, fntpLambda0);
			qaMomRes("VtxFit_", lambda0Fit, fntpLambda0);


			// store only best fitted candidate
			//      if( bestVtxFit[j]==1){// && bestMassFitLambda0[j]>0){
			Lambda0Fit.Append(lambda0Fit);
			//      }

			RhoCandidate * truth = lambda0[j]->GetMcTruth();
			TLorentzVector l;



			if(0x0 != truth){
				l = truth->P4();
				qa.qaVtx("McTruth_", truth, fntpLambda0);
	//			dl = truth->Daughter(0)->Pos();
			}
			else{
				qa.qaVtx("McTruth_", dummyCand, fntpLambda0);
			}


			qa.qaP4("McTruth_", l, fntpLambda0);




			//***information of boosted particle
			lambda0Fit->Boost(-beamBoost);
			qa.qaComp("boost_", lambda0Fit, fntpLambda0);

			fntpLambda0->DumpData();
		}

	    //***AntiLambda0 -> PiMinus + Proton
	    antiLambda0.Combine(piplus,antiProton);
		antiLambda0.Select(lambdaMassSelector);
	    antiLambda0.SetType(-3122);

	//    std::map<int,int> bestVtxFitAntiLambda0, bestMassFitAntiLambda0;
	//    bestVtxFitAntiLambda0 = jenny::VertexQaIndex(&antiLambda0);
	//    bestMassFitAntiLambda0 = jenny::MassFitQaIndex(&antiLambda0, m0_lambda0);


	    for (int j=0; j<antiLambda0.GetLength(); ++j){

			//general info about event
			fntpAntiLambda0->Column("ev",     (Float_t) evt);
			fntpAntiLambda0->Column("cand",    (Float_t) j);
			fntpAntiLambda0->Column("ncand",   (Float_t) antiLambda0.GetLength());
			fntpAntiLambda0->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(antiLambda0[j]));

			RhoCandidate * mother = antiLambda0[j]->TheMother();
			int moth = (mother==0x0) ? 88888 : mother->PdgCode();
			fntpAntiLambda0->Column("Mother", (Float_t) moth);

			qa.qaP4("AntiLambda0_", antiLambda0[j]->P4(), fntpAntiLambda0);
			qa.qaComp("AntiLambda0_", antiLambda0[j], fntpAntiLambda0);
			qa.qaEventShapeShort("es_", &evsh, fntpAntiLambda0);




			// do vertex fit


			PndKinVtxFitter vertexfitterAntiLambda0 (antiLambda0[j]);
			vertexfitterAntiLambda0.Fit();
			RhoCandidate * antiLambda0Fit = antiLambda0[j]->GetFit();


			// store info of vertex fit


			qa.qaFitter("VtxFit_", &vertexfitterAntiLambda0, fntpAntiLambda0);
			//      fntpAntiLambda0->Column("VtxFit_HowGood", (Int_t) bestVtxFitAntiLambda0[j]);

			qa.qaVtx("VtxFit_", antiLambda0Fit, fntpAntiLambda0);


			// difference to MCTruth
			qa.qaMcDiff("VtxFit_", antiLambda0Fit, fntpAntiLambda0);
			qaVtxDiff("VtxFit_", antiLambda0Fit, fntpAntiLambda0);
			qaMomRes("VtxFit_", antiLambda0Fit, fntpAntiLambda0);



			//      if(prob>0.01 && chi2<bestchi2){
			//    	  bestchi2 = chi2;
			//    	  AntiLambda0Fit.InsertAt(0,antiLambda0Fit);
			//      }


			// do mass fit
			PndKinFitter massFitterAntiLambda0(antiLambda0Fit);
			massFitterAntiLambda0.AddMassConstraint(fm0_lambda0);
			massFitterAntiLambda0.Fit();

			RhoCandidate * antiLambda0Fit_mass = antiLambda0Fit->GetFit();

			qa.qaFitter("MassFit_", &massFitterAntiLambda0, fntpAntiLambda0);
			//      fntpAntiLambda0->Column("MassFit_HowGood", (Int_t) bestMassFitAntiLambda0[j]);
			qa.qaMcDiff("MassFit_", antiLambda0Fit_mass, fntpAntiLambda0);


			//      if(bestVtxFitAntiLambda0[j]==1){ //&& bestMassFitAntiLambda0[j]>0){
			AntiLambda0Fit.Append(antiLambda0Fit);
			//      }

			RhoCandidate * truth = antiLambda0[j]->GetMcTruth();
			TLorentzVector l;
			if(0x0 != truth){
				l = truth->P4();
				qa.qaVtx("truth_", truth, fntpAntiLambda0);
			}
			qa.qaP4("truth_", l, fntpAntiLambda0);

			//***information of boosted particle
			antiLambda0Fit->Boost(-beamBoost);
			qa.qaComp("boost_", antiLambda0Fit, fntpAntiLambda0);


			fntpAntiLambda0->DumpData();
	    }

		//*** Cross check: pbar + p -> Lambda0 + AntiLambda0

		crossCheck.Combine(Lambda0Fit, AntiLambda0Fit);
		crossCheck.SetType(88888);


		for (int j=0; j<crossCheck.GetLength(); ++j){

			//general information about event
			fntpCrossCheck->Column("ev",     (Float_t) evt);
			fntpCrossCheck->Column("cand",    (Float_t) j);
			fntpCrossCheck->Column("ncand",   (Float_t) crossCheck.GetLength());
			fntpCrossCheck->Column("McTruthMatch", (bool) fAnalysis->McTruthMatch(crossCheck[j]));

			qa.qaP4("", crossCheck[j]->P4(), fntpCrossCheck);
			qa.qaComp("", crossCheck[j], fntpCrossCheck);
			qa.qaPoca("", crossCheck[j], fntpCrossCheck);
			qa.qaEventShapeShort("es_", &evsh, fntpAntiLambda0);


			//do vertex fit



			PndKinVtxFitter vertexFitter_cc (crossCheck[j]);
			vertexFitter_cc.Fit();
			RhoCandidate * ccFit = crossCheck[j]->GetFit();

			//store info of vertex fit
			qa.qaFitter("VtxFit_", &vertexFitter_cc, fntpCrossCheck);
			qa.qaP4("VtxFit_", ccFit->P4(), fntpCrossCheck);
			qa.qaVtx("VtxFit_", ccFit, fntpCrossCheck);
			qa.qaMcDiff("VtxFit_", ccFit, fntpCrossCheck);
			qaVtxDiff("VtxFit_", ccFit, fntpCrossCheck);
			qaMomRes("VtxFit_", ccFit, fntpCrossCheck);



			float mc_mass_l0 = 0.;

			RhoCandidate * truth = ccFit->GetMcTruth();
		    TLorentzVector l;


		    if(0x0 != truth){
		    	l = truth->P4();
		    	qa.qaVtx("McTruth_", truth, fntpCrossCheck);
		    }
		    else{
		    	qa.qaVtx("McTruth_", dummyCand, fntpCrossCheck);
		    }

  		    qa.qaP4("McTruth_", l, fntpCrossCheck);


			//***do 4c fit
			PndKinFitter cc_Fitter4c (crossCheck[j]);
			cc_Fitter4c.Add4MomConstraint(fini);
			cc_Fitter4c.Fit();

			RhoCandidate * ccFit4C = crossCheck[j]->GetFit();

			// store info of 4c Fit
			fntpCrossCheck->Column("f4c_Chi2_cc", (Float_t) cc_Fitter4c.GetChi2());
			fntpCrossCheck->Column("f4c_NDF_cc", (Float_t) cc_Fitter4c.GetNdf());
			fntpCrossCheck->Column("f4c_Prob_cc", (Float_t) cc_Fitter4c.GetProb());

			qa.qaComp("f4c_", ccFit4C, fntpCrossCheck);
			qa.qaP4("f4C_", ccFit4C->P4(), fntpCrossCheck);
			qa.qaVtx("f4C_", ccFit4C, fntpCrossCheck);

			//difference to MC Truth
			qa.qaMcDiff("f4c_", ccFit4C, fntpCrossCheck);
			qaVtxDiff("f4c_", ccFit4C, fntpCrossCheck);
			qaMomRes("f4c_", ccFit4C, fntpCrossCheck);


			//do kalman vertex fit

			PndKalmanVtxFitter kalmanfitter (crossCheck[j]);
			kalmanfitter.Fit();
			RhoCandidate * kalmanFit = crossCheck[j]->GetFit();

			qa.qaFitter("KalmanFit_", &kalmanfitter, fntpCrossCheck);
			qa.qaVtx("KalmanFit_", kalmanFit, fntpCrossCheck);
			qa.qaP4("KalmanFit_", kalmanFit->P4(), fntpCrossCheck);

			qaVtxDiff("KalmanFit_", kalmanFit, fntpCrossCheck);
			qaMomRes("KalmanFit_", kalmanFit, fntpCrossCheck);

			//do mass fit
			PndKinFitter massFitter_cc (ccFit);
			massFitter_cc.AddMassConstraint(fm0_beam);
			massFitter_cc.Fit();

			RhoCandidate * massFit_cc = ccFit->GetFit();

			//store fit results
			fntpCrossCheck->Column("fmass_Chi2_cc", (Float_t) massFitter_cc.GetChi2());
			fntpCrossCheck->Column("fmass_NDF_cc", (Float_t) massFitter_cc.GetNdf());
			fntpCrossCheck->Column("fmass_Prob_cc", (Float_t) massFitter_cc.GetProb());

			qa.qaMcDiff("fmassMCDiff_", massFit_cc, fntpCrossCheck);

			fntpCrossCheck->DumpData();
		}
	 Lambda0Fit.Cleanup();
	 AntiLambda0Fit.Cleanup();
	}

}

void AnalysisTaskLambda0::Finish()
{
	  foutput->cd();

	  fntpMC -> GetInternalTree()->Write();
	  fntpPiMinus ->GetInternalTree()->Write();
	  fntpPiPlus->GetInternalTree()->Write();
	  fntpProton->GetInternalTree()->Write();
	  fntpAntiProton->GetInternalTree()->Write();
	  fntpLambda0->GetInternalTree()->Write();
	  fntpAntiLambda0->GetInternalTree()->Write();
	  fntpCrossCheck->GetInternalTree()->Write();

	  foutput->Save();

	  ftimer.Stop();
	  Double_t rtime = ftimer.RealTime();
	  Double_t ctime = ftimer.CpuTime();

	  cout<<"Macro finisched successfully."<<endl;
	  cout<<"Realtime: "<<rtime<<" s, CPU time: "<<ctime<<" s"<<endl;
	  cout<<endl;
}

