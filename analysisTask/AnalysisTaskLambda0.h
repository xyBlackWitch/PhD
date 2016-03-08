#ifndef AnalysisTaskLambda0_H
#define AnalysisTaskLambda0_H 1

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

class AnalysisTaskLambda0 : public FairTask
{

	public:

		//** default constructor
		AnalysisTaskLambda0();

		//** destructor
		~AnalysisTaskLambda0();

		virtual InitStatus Init();

		virtual void Exec(Option_t* opt);

		virtual void Finish();

		void SetOutPutDir(const char* fname=""){
			fOutPath=fname;
		}

		void SetOutputName(const char* outname="output_ana.root"){
			foutname=outname;
		}

		void SetNEvents(double nevents=0){
			fnevts = nevents;
		}
		void SetMom(double mom=1.7){
			fmom=mom;
		}



	protected:


	private:

//		enum pidNumbers;


		void numberOfHitsInSubdetector(TString pre, RhoCandidate *c, RhoTuple *n);
		void tagNHits(TString pre, RhoCandidate *c, RhoTuple *n);
		int tagHits(RhoCandidate *c);

		std::map<int,int> VertexQaIndex(RhoCandList* candList, float probLimit);
		std::map<int,int> MassFitQaIndex(RhoCandList* candList, float m0, float probLimit);

		void qaVtxDiff(TString pre, RhoCandidate * c, RhoTuple * n);
		void qaMomRes(TString pre, RhoCandidate * c, RhoTuple * n);

		TStopwatch ftimer;

		double fnevts;
		int fEvtCount;

		double fmom;
		double fm0_lambda0;
		double fm0_beam;

		const char* fOutPath;
		const char* foutname;

		RhoTuple * fntpMC;
		RhoTuple * fntpPiMinus;
		RhoTuple * fntpPiPlus;
		RhoTuple * fntpProton;
		RhoTuple * fntpAntiProton;
		RhoTuple * fntpLambda0;
		RhoTuple * fntpAntiLambda0;
		RhoTuple * fntpSys;

		TFile* foutput;

		PndAnalysis* fAnalysis;

		RhoMassParticleSelector* lambdaMassSelector;

		TLorentzVector fini;

		ClassDef(AnalysisTaskLambda0,1);



};
#endif
