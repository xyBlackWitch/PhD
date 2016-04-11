#ifndef AnalysisTaskXi1690_H
#define AnalysisTaskXi1690_H 1

#include "FairTask.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include <map>



class PndAnalysis;
class RhoFitters;
class RhoTools;
class RhoCandList;
class RhoCandidate;
class RhoTuple;
class RhoMassParticleSelector;

class AnalysisTaskXi1690 : public FairTask
{

	public:

		//** default constructor
		AnalysisTaskXi1690();

		//** destructor
		~AnalysisTaskXi1690();

		virtual InitStatus Init();

		virtual void Exec(Option_t* opt);

		virtual void Finish();

		void SetOutPutDir(const char* fname=""){
			fOutPath=fname;
		}

		void SetNEvents(double nevents=0){
			fnevts = nevents;
		}
		void SetMom(double mom=4.1){
			fmom=mom;
		}



	protected:


	private:

//		enum pidNumbers;

		void CombinedList(RhoCandidate *cand, RhoCandList *combinedList, int pdg);
		void GetNotCombinedList(RhoCandList combinedList, RhoCandList * candList);

		void numberOfHitsInSubdetector(TString pre, RhoCandidate *c, RhoTuple *n);
		void tagNHits(TString pre, RhoCandidate *c, RhoTuple *n);
		int tagHits(RhoCandidate *c);

		std::map<int,int> VertexQaIndex(RhoCandList* candList, float probLimit);
		std::map<int,int> MassFitQaIndex(RhoCandList* candList, float m0, float probLimit);

		void qaVtxDiff(TString pre, RhoCandidate * c, RhoTuple * n);
		void qaMomRes(TString pre, RhoCandidate * c, RhoTuple * n);

		TStopwatch ftimer;

		int fEvtCount;
		double fnevts;

		double fmom;
		double fm0_lambda0;
		double fm0_Xi;
		double fm0_Xi1690;
		double fm0_beam;

		const char* fOutPath;

		RhoTuple * fntpMc;
		RhoTuple * fntpPiMinus;
		RhoTuple * fntpPiPlus;
		RhoTuple * fntpProton;
		RhoTuple * fntpAntiProton;
		RhoTuple * fntpKaonMinus;
		RhoTuple * fntpKaonPlus;
		RhoTuple * fntpLambda0;
		RhoTuple * fntpAntiLambda0;
		RhoTuple * fntpXiMinus1690;
		RhoTuple * fntpXiPlus;
		RhoTuple * fntpXiSys;

		TFile* foutput;

		PndAnalysis* fAnalysis;

		RhoMassParticleSelector* lambdaMassSelector;
		RhoMassParticleSelector* xiMassSelector;
		RhoMassParticleSelector* Xi1690MassSelector;

		TLorentzVector fini;

		ClassDef(AnalysisTaskXi1690,1);



};
#endif
