/**
* @file common.cpp
* @mainpage common.cpp Helper Functions
* @author Jennifer Puetz (jennifer.puetz@gmx.de)
* @date 2015
* @brief little functions for ROOT
* @details This file holds all the common functions I wrote while working on stuff.
* Although they are intended to be used from ROOT macros, they included the `#include`s to let it compile as well.
*
* If you want to use the methods as well from ROOT's command line, make sure to include this file in your `~/.rootrc` file with
* ~~~
* Rint.Load: ~/Documents/Coding/PhysicsAnalysis/common.cpp
* ~~~
*/


#include "TH1.h"
#include "TString.h"
#include "TCanvas.h"
#include <map>

//#include "RhoTuple.h"
//#include "RhoCandidate.h"

namespace jenny{

	void CreateDrawAndSaveHistogram(TH1* &histo, TString outputdir, TString outputname, bool saveoutput, bool close){

		/** @brief  saves Histogramm as *.root and *.png and if wanted closes the histograms at the end
		*	@details This mehtod create a histogramm and save it as root and png file. If you choose close, the canvas is closed after the histogram was saved
		*/

		TString name = TString(histo->GetName());
		TString title = TString(histo->GetTitle());

		TCanvas * canvas = new TCanvas("c_"+name, title, 0,0,800,500);
		//gStyle->SetOptStat(11);
		histo->Draw();


		if (saveoutput){
			canvas->Print(outputdir + "root-files/" + outputname + ".root");
			canvas->Print(outputdir + "png-files/" + outputname + ".png");
		}

		if (close) canvas->Close();


	}

	void CreateDrawAndSaveHistogramWithFit(TH1* &histo, TString outputdir, TString outputname, bool saveoutput, bool close, bool autorange = true, double innerRange = 0.1 , double outerRange=1){

			/** @brief  saves Histogramm as *.root and *.png and if wanted closes the histograms at the end
			*	@details This mehtod create a histogramm and save it as root and png file. If you choose close, the canvas is closed after the histogram was saved
			*/

			TString name = TString(histo->GetName());
			TString title = TString(histo->GetTitle());

			TCanvas * canvas = new TCanvas("c_"+name, title, 0,0,800,500);
			//gStyle->SetOptStat(11);
			histo->Draw();

			TF1 * fit = doubleGaussFit(histo, autorange, innerRange, outerRange);
			fit->Draw("SAME");


			if (saveoutput){
				canvas->Print(outputdir + "root-files/" + outputname + ".root");
				canvas->Print(outputdir + "png-files/" + outputname + ".png");
			}

			if (close) canvas->Close();


		}

	void CreateDrawAndSaveHistogram(TH2* &histo, TString outputdir, TString outputname, bool saveoutput, bool close){

		/** @brief  saves Histogramm as *.root and *.png and if wanted closes the histograms at the end
		*	@details This mehtod create a histogramm and save it as root and png file. If you choose close, the canvas is closed after the histogram was saved
		*/

		TString name = TString(histo->GetName());
		TString title = TString(histo->GetTitle());

		TCanvas * canvas = new TCanvas("c_"+name, title, 0,0,800,500);
		//gStyle->SetOptStat(0);
		histo->Draw("COLZ");


		if (saveoutput){
			canvas->Print(outputdir + "root-files/" + outputname + ".root");
			canvas->Print(outputdir + "png-files/" + outputname + ".png");
		}

		if (close) canvas->Close();


	}

	  std::map<int,int> VertexQaIndex(RhoCandList* candList, float probLimit=0.01){
		  /** @brief  give back the order of the best chi2
		   * @details give back the order of the best chi2!  1 means best, 2: second best (analoge for bad chi2 with negative values)
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

		  int running =0;

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

		  int run = 0;

		  for (is_bad = chi2_bad.begin(); is_bad != chi2_bad.end(); is_bad++, run++){
			  bestMassFit[is_bad->second] = - (run + 1);
		  }


		  return bestMassFit;
	  }

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
		  return combinedList;

	  }

	  TF1 * doubleGaussFit(TH1 * hist, bool autorange, double inner, double outer){
		  /**
		   * @brief performes a double gaussian Fit
		   */
		  Double_t parameters[6] = {0,0,0,0,0,0};
		  gStyle->SetOptFit(1);

		  double center = hist->GetMean();

		  if (autorange){
			  double inner = fabs(center - hist->GetXaxis()->GetXmax()) /10;
			  double outer = fabs(center - hist->GetXaxis()->GetXmax()) /10*8;
		  }

		  TF1 * fit1 = new TF1("fit1", "gaus", center-inner, center+inner );

		  hist->Fit(fit1, "Q0R");
		  fit1->GetParameters(&parameters[0]);

		  TF1 * fit2 = new TF1("fit2", "gaus", center-outer, center+outer );

		  hist->Fit(fit2, "Q0R");
		  fit1->GetParameters(&parameters[3]);

		  TF1 * fitproper = new TF1("fitproper", "gaus(0)+gaus(3)", center-outer, center+outer);
		  fitproper->SetParameters(parameters);
		  fitproper->SetParName(0, "Const.(inner)");
		  fitproper->SetParName(1, "Mean (inner)");
		  fitproper->SetParName(2, "Sigma (inner)");
		  fitproper->SetParName(3, "Const.(outer)");
		  fitproper->SetParName(4, "Mean (outer)");
		  fitproper->SetParName(5, "Sigma (outer)");

		  hist->Fit(fitproper, "Q0R");

		  fitproper->SetLineColor(hist->GetLineColor());
		  fitproper->SetLineWidth(hist->GetLineWidth());
		  fitproper->SetLineStyle(2);

		  return fitproper;

	  }

	  void qaP3(TString pre, TVector3 v, RhoTuple *n, bool skip=false){

		  if (n==0) return;

		  if(!skip){
			  n->Column(pre+"decayvx", (Float_t) v.x(), 0.0f);
			  n->Column(pre+"decayvy", (Float_t) v.y(), 0.0f);
			  n->Column(pre+"decayvz", (Float_t) v.z(), 0.0f);
		  }
		  else{
			  n->Column(pre+"decayvx", (Float_t) -999, 0.0f);
			  n->Column(pre+"decayvy", (Float_t) -999, 0.0f);
			  n->Column(pre+"decayvz", (Float_t) -999, 0.0f);
		  }

	  }


	  void numberOfHitsInSubdetector(TString pre="", RhoCandidate *c, RhoTuple *n){

		/* This void method saves the number of Hits in the MVD, STT and GEM detector
		 * into the RhoTuple.
		 * Useful only for final state particle!
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


	  void tagNHits(TString pre= "", RhoCandidate *c, RhoTuple *n){

	  	/**@brief Tag the particle with different integers
	  	 * @details Tag the particle with different integers:
	  	 * 0: if there is no hit in the detector
	  	 * 1: sttHits>3 or mvdHits>3 or gemHit>3
	  	 */

	  	bool comp = c->IsComposite();
	  	int tag=0;

	  	if (!comp){

	  		PndPidCandidate * pidCand = (PndPidCandidate*)c->GetRecoCandidate();

	  		if(pidCand){
	  			int mvdHits = pidCand->GetMvdHits();
	  			int sttHits = pidCand->GetSttHits();
	  			int gemHits = pidCand->GetGemHits();

	  			if(mvdHits>3 || sttHits>3 || gemHits>3) tag=1;


	  		}

	  	}
	  	else{
	  		int dtag[2] = {0,0};
	  		int ndau = c->NDaughters();

	  		for(int dau=0; dau<ndau; dau++){
	  			RhoCandidate * daughter = c->Daughter(dau);
	  			dtag[dau] = tagNHits(daughter);
	  		}
	  		if(dtag[0]==1 && dtag[1]==1) tag=1;

	  	}

	  	n->Column(pre + "HitTag", (Int_t) tag, 0);
	  }



	  int tagNHits(RhoCandidate *c){

		/**@brief Tag the particle with different integers
		 * @details Tag the particle with different integers:
		 * 0: if there is no hit in the detector
		 * 1: sttHits>3 or mvdHits>3 or gemHit>3
		 */
	  	int tag = 0;
	  	bool comp = c->IsComposite();


	  	if(!comp){

	  		PndPidCandidate * pidCand = (PndPidCandidate*)c->GetRecoCandidate();

	  		if(pidCand){
	  			int mvdHits = pidCand->GetMvdHits();
	  			int sttHits = pidCand->GetSttHits();
	  			int gemHits = pidCand->GetGemHits();

	  			if(mvdHits>3 || sttHits>3 || gemHits>3) tag=1;
	  		//		cout << "mvd: " << mvdHits << " stt: " << sttHits << " gem: " << gemHits << endl;

	  		}
	  	}
	  	else{
	  		int dtag[2] = {0,0};
	  		int ndau = c->NDaughters();

	  		for(int dau=0; dau<ndau; dau++){
	  			RhoCandidate * daughter = c->Daughter(dau);
	  			dtag[dau] = tagNHits(daughter);
	  		}
	  		if(dtag[0]==1 && dtag[1]==1) tag=1;
	  	}

	  	return tag;
	  }

}
