//class PndAnalysis;
//class PndAnaPidSelector;
//class RhoCandList;
//class RhoTuple;

#include "AnalysisTask.h"

// C++ headers
#include <string>
#include <iostream>

//ROOT headers
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

//Rho headers
#include "RhoTuple.h"
#include "PndRhoTupleQA.h"
#include "RhoCandList.h"


//Analysis headers
#include "PndAnalysis.h"

using std::cout;
using std::endl;


ClassImp(AnalysisTask)

//---------------Default constructor------------
AnalysisTask::AnalysisTask():
	FairTask("Panda Tutorial Analysis Task") {
	}
//------------------------------------------

//--------------- destructor ----------------
AnalysisTask::~AnalysisTask(){}
//-----------------------------------------


InitStatus AnalysisTask::Init(){
  

  //*** create tuples
  fntpMc = new RhoTuple("ntpMC", "MCTruth info");


  //Create output file for histograms
  foutput = TFile::Open(outPath+"output_ana_test.root","RECREATE");

  // data reader Object
  fAnalysis = new PndAnalysis();

  //***Mass selector
  double m0_beam = TDatabasePDG::Instance()->GetParticle("pbarpSystem")->Mass();
  cout<<"Mass pbar p system: "<<m0_beam<<endl;

  //*** lorentz vector of the initial particle
  fmom = 6.231552;
  double p_m0 = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
  fini.SetXYZT(0,0, fmom, sqrt(p_m0*p_m0+ fmom*fmom)+p_m0);
  
  return kSUCCESS;
}

void AnalysisTask::Exec(Option_t* op)
{
	RhoCandList mclist, all;
	PndRhoTupleQA qa(fAnalysis, fmom);
	fEvtCount = 0;
	int nevts = fAnalysis->GetEntries();

	while (fAnalysis->GetEvent() && fEvtCount++<nevts){
    if ((fEvtCount%100)==0) cout << "evt "<< fEvtCount <<endl;

    //***get MC list and store info
    fAnalysis->FillList(mclist, "McTruth");
    qa.qaMcList("", mclist, fntpMc);
    fntpMc->DumpData();
   }
}

void AnalysisTask::Finish()
{

  //Write output
  foutput->cd();
  fntpMc -> GetInternalTree()->Write();

  foutput->Save();

}

