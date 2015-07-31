#ifndef AnalysisTask_H
#define AnalysisTask_H 1

#include "FairTask.h"
#include "TFile.h"
#include "TLorentzVector.h"

class PndAnalysis;
class RhoCandList;
class RhoTuple;

class AnalysisTask : public FairTask
{

	public:

		//** default constructor
		AnalysisTask();

		//** destructor
		~AnalysisTask();

		virtual InitStatus Init();

		virtual void Exec(Option_t* opt);

		virtual void Finish();

	protected:


	private:

		int fEvtCount;
		double fmom;

		RhoTuple *fntpMc;

		TFile* foutput;

		PndAnalysis* fAnalysis;

		TLorentzVector fini;

		ClassDef(AnalysisTask,1);



};
#endif
