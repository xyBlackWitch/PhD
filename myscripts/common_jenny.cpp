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
#include "TMath.h"
#include <map>
#include "TDatabasePDG"
#include "common_andi.cpp"
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/PandaSmartLabel.C"

//#include "RhoTuple.h"
//#include "RhoCandidate.h"

namespace jenny{

	void CreateDrawAndSaveHistogram(TH1* &histo, TString outputdir, TString outputname, bool saveoutput, bool close, bool prelim=false){

		/** @brief  saves Histogramm as *.root and *.png and if wanted closes the histograms at the end
		*	@details This mehtod create a histogramm and save it as root and png file. If you choose close, the canvas is closed after the histogram was saved
		*/

		TString name = TString(histo->GetName());
		TString title = TString(histo->GetTitle());

		TCanvas * canvas = new TCanvas("c_"+name, title, 0,0,1500,1000);
		gStyle->SetOptStat(1111);
		gStyle->SetOptFit(0);
		gStyle->SetStatY(0.85);
		gStyle->SetStatX(0.85);
		gStyle->SetStatW(0.15);
		gStyle->SetStatH(0.15);


		histo->GetXaxis()->SetLabelSize(0.045);
		histo->GetXaxis()->SetTitleSize(0.05);
		histo->GetXaxis()->SetTitleOffset(0.90);
		histo->GetYaxis()->SetLabelSize(0.045);
		histo->GetYaxis()->SetTitleSize(0.05);
		histo->GetYaxis()->SetTitleOffset(0.80);
		histo->Draw();

		if(prelim) PandaSmartLabel();


		if (saveoutput){
			canvas->Print(outputdir + "root-files/" + outputname + ".root");
			canvas->Print(outputdir + "png-files/" + outputname + ".png");
			canvas->Print(outputdir + "pdf-files/" + outputname + ".pdf");
		}

		if (close) canvas->Close();


	}

	void CreateDrawAndSaveHistogram(TH2* &histo, TString outputdir, TString outputname, bool saveoutput, bool close, bool prelim=false){

			/** @brief  saves Histogramm as *.root and *.png and if wanted closes the histograms at the end
			*	@details This mehtod create a histogramm and save it as root and png file. If you choose close, the canvas is closed after the histogram was saved
			*/

			TString name = TString(histo->GetName());
			TString title = TString(histo->GetTitle());

			TCanvas * canvas = new TCanvas("c_"+name, title, 0,0,1500,1000);
			gStyle->SetOptStat(0);
			gStyle->SetStatY(0.85);
			gStyle->SetStatX(0.85);
			gStyle->SetStatW(0.15);
			gStyle->SetStatH(0.15);

			histo->GetXaxis()->SetLabelSize(0.045);
			histo->GetXaxis()->SetTitleSize(0.05);
			histo->GetXaxis()->SetTitleOffset(0.90);
			histo->GetYaxis()->SetLabelSize(0.045);
			histo->GetYaxis()->SetTitleSize(0.05);
			histo->GetYaxis()->SetTitleOffset(0.80);

			histo->Draw("COLZ");
			if(prelim) PandaSmartLabel();


			if (saveoutput){
				canvas->Print(outputdir + "root-files/" + outputname + ".root");
				canvas->Print(outputdir + "png-files/" + outputname + ".png");
				canvas->Print(outputdir + "pdf-files/" + outputname + ".pdf");
			}

			if (close) canvas->Close();


		}

	void CreateDrawAndSaveHistogramWithFit(TH1* &histo, TString outputdir, TString outputname, bool saveoutput, bool close, bool autorange = true, double innerRange = 0.1 , double outerRange=1, bool excludeCenter=false, int fittype=0){

			/** @brief  saves Histogramm as *.root and *.png and if wanted closes the histograms at the end
			*	@details This mehtod create a histogramm and save it as root and png file. If you choose close, the canvas is closed after the histogram was saved
			*/

			Double_t parameters[6] = {0.,0.,0.,0.,0.,0.};
			TString name = TString(histo->GetName());
			TString title = TString(histo->GetTitle());

			TCanvas * canvas = new TCanvas("c_"+name, title, 0,0,1500,1000);
			gStyle->SetOptStat(1000000001);
			gStyle->SetOptFit(111);
			gStyle->SetFitFormat("5.8g");
			gStyle->SetStatY(0.9);
			gStyle->SetStatX(0.9);
			gStyle->SetStatW(0.35);
			gStyle->SetStatH(0.45);

			histo->GetXaxis()->SetLabelSize(0.045);
			histo->GetXaxis()->SetTitleSize(0.05);
			histo->GetXaxis()->SetTitleOffset(0.90);
			histo->GetYaxis()->SetLabelSize(0.045);
			histo->GetYaxis()->SetTitleSize(0.05);
			histo->GetYaxis()->SetTitleOffset(0.80);


			histo->Draw();

			TF1 * fit;
//
			if(excludeCenter){
				fit = andi::doubleGaussFitExcludeCenter(histo, false, innerRange, outerRange);
			}
			else if (fittype==1){
				fit = andi::gaussFit(histo);
			}
			else if (fittype==2){
				fit = doubleGaussFit(histo, autorange, innerRange, outerRange);
			}
			else{
				std::cout << "Type of fit is not defined"<< std::endl;
			}
			fit->SetLineColor(kRed);
			fit->SetLineStyle(7);
			fit->SetLineWidth(3);
			fit->Draw("SAME");

			if (saveoutput){
				canvas->Print(outputdir + "root-files/" + outputname + ".root");
				canvas->Print(outputdir + "png-files/" + outputname + ".png");
				canvas->Print(outputdir + "pdf-files/" + outputname + ".pdf");
			}

			if (close) canvas->Close();


		}

	void CreateDrawAndSaveHistogramFit(TH1* &histo, TString outputdir, TString outputname, bool saveoutput, bool close, bool autorange = true, double innerRange = 0.1 , double outerRange=1, bool excludeCenter=false){
		return CreateDrawAndSaveHistogramWithFit(histo, outputdir, outputname, saveoutput, close, autorange, innerRange, outerRange,excludeCenter,1);
	}

	void CreateDrawAndSaveHistogramDoulbeFit(TH1* &histo, TString outputdir, TString outputname, bool saveoutput, bool close, bool autorange = true, double innerRange=0.1 , double outerRange=1, bool excludeCenter=false){
		return CreateDrawAndSaveHistogramWithFit(histo, outputdir, outputname, saveoutput, close, autorange, innerRange, outerRange,excludeCenter,2);
	}

	void GetFitParameterDoubleFit(TH1* &histo, bool autorange = true, double innerRange=0.1 , double outerRange=1, bool excludeCenter=false){

		Double_t parameter[6] = {0,0,0,0,0,0};
		TF1 * fit;

		if(excludeCenter){
			fit = andi::doubleGaussFitExcludeCenter(histo, false, innerRange, outerRange);
		}
		else{
			fit = doubleGaussFit(histo, autorange, innerRange, outerRange);
		}

		fit->GetParameters(&parameter[0]);
		//return parameter;
	}

	void GetFitParameterErrorDoubleFit(TH1* &histo, bool autorange = true, double innerRange=0.1 , double outerRange=1, bool excludeCenter=false){

			Double_t errors[6] = {0,0,0,0,0,0};
			TF1 * fit;

			if(excludeCenter){
				fit = andi::doubleGaussFitExcludeCenter(histo, false, innerRange, outerRange);
			}
			else{
				fit = doubleGaussFit(histo, autorange, innerRange, outerRange);
			}

			for (int i=0; i<6; i++){

				errors[i]=fit->GetParError(i);
			}


			//return parameter;
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

		  int running =0;

		  for (is_bad = chi2_bad.begin(); is_bad != chi2_bad.end(); is_bad++, running++){
			  indexBestFit[is_bad->second] = - (running + 1);
		  }


		  return indexBestFit;
	  }

	  void CreateDrawAndSaveNHistograms(TH1* &h1, TH1* &h2, TString leg1="", TString leg2="", TString outputdir, TString outputname, bool saveoutput, bool close, bool prelim=false){
	  // should be added to common_jenny.cpp

	  	/** @brief  saves 2 histograms in same file as *.root and *.png and if wanted closes the canvas at the end
	  	*	@details This mehtod create 2 histograms and save it as root and png file. If you choose close, the canvas is closed after the histograms are saved
	  	*/

	  	TString name = TString(h1->GetName());
	  	TString title = TString(h1->GetTitle());

	  	h1->SetLineColor(kBlue);
	  	h2->SetLineColor(kRed);


	  	TLegend * legend = new TLegend(0.7,0.62,0.86,0.795, "");
	  	legend->AddEntry(h1, leg1, "l");
	  	legend->AddEntry(h2, leg2, "l");


	  	TCanvas * canvas = new TCanvas("c_"+name, title, 0,0,1500,1000);
	  	gStyle->SetOptStat(0);

	  	h1->Draw();
	  	h2->Draw("SAME");
	  	legend->Draw();

	  	if(prelim) PandaSmartLabel("Lprel");


	  	if (saveoutput){
	  		canvas->Print(outputdir + "root-files/" + outputname + ".root");
	  		canvas->Print(outputdir + "png-files/" + outputname + ".png");
			canvas->Print(outputdir + "pdf-files/" + outputname + ".pdf");
	  	}

	  	if (close) canvas->Close();

	  }

	  void CreateDrawAndSaveNHistograms(TH1* &h1, TH1* &h2, TH1* &h3, TString leg1="", TString leg2="", TString leg3="", TString outputdir, TString outputname, bool saveoutput, bool close,bool prelim=false){
	  // should be added to common_jenny.cpp

	  	/** @brief  saves 3 histograms in same file as *.root and *.png and if wanted closes the canvas at the end
	  	*	@details This mehtod create 3 histograms and save it as root and png file. If you choose close, the canvas is closed after the histograms are saved
	  	*/

	  	TString name = TString(h1->GetName());
	  	TString title = TString(h1->GetTitle());

	  	h1->SetLineColor(kBlue);
	  	h2->SetLineColor(kRed);
	  	h3->SetLineColor(kBlack);

	  	TLegend * legend = new TLegend(0.7,0.62,0.86,0.795, "");
	  	legend->AddEntry(h1, leg1, "l");
	  	legend->AddEntry(h2, leg2, "l");
	  	legend->AddEntry(h3, leg3, "l");


	  	TCanvas * canvas = new TCanvas("c_"+name, title, 0,0,1500,1000);
	  	gStyle->SetOptStat(0);
	  	h1->Draw();
	  	h2->Draw("SAME");
	  	h3->Draw("SAME");
	  	legend->Draw();
	  	if(prelim) PandaSmartLabel();


	  	if (saveoutput){
	  		canvas->Print(outputdir + "root-files/" + outputname + ".root");
	  		canvas->Print(outputdir + "png-files/" + outputname + ".png");
			canvas->Print(outputdir + "pdf-files/" + outputname + ".pdf");
	  	}

	  	if (close) canvas->Close();

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


	  std::map<int,int> FourConstraintFitQaIndex(RhoCandList* candList, double mom=0, float probLimit=0.01){
	 		  /** @brief  give back the order of the best chi2 for MassFit
	 		   * @details give back the order of the best chi2 for the MassFit!  1 means best, 2: second best (analoge for bad chi2 with negative values)
	 		   */

	 		  if(mom==0) std::cout << "Initial 4 mometum is missing for 4C fit" << std::endl;

	 		  double p_m0 = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
	 		  TLorentzVector ini (0,0, mom, sqrt(p_m0*p_m0+ mom*mom)+p_m0);
	 		  
	 		  std::map<double, int> chi2_good, chi2_bad;

	 		  for (int i=0; i<candList->GetLength(); i++){

	 			  PndKinFitter fitter4C(candList->Get(i));
	 			  fitter4c.Add4MomConstraint(ini);
	 			  fitter4C.Fit();

	 			  bool failedchi2 = TMath::IsNaN(fitter4C.GetChi2());
	 			  bool failedprob = TMath::IsNaN(fitter4C.GetProb());

	 			  if(!failedchi2 && !failedprob){

	 				  if (fitter4C.GetProb() > probLimit){
	 					  chi2_good[fitter4C.GetChi2()]=i;
	 				  }
	 				  else{
	 					  chi2_bad[fitter4C.GetChi2()]=i;
	 				  }
	 			  }
	 		  }

	 		  std::map<double,int>::iterator is_good, is_bad;
	 		  std::map<int,int> best4CFit;

	 		  int run =0;

	 		  for (is_good = chi2_good.begin(); is_good != chi2_good.end(); is_good++, run++){
	 			  best4CFit[is_good->second] = run + 1;
	 		  }

	 		  int run = 0;

	 		  for (is_bad = chi2_bad.begin(); is_bad != chi2_bad.end(); is_bad++, run++){
	 			  best4CFit[is_bad->second] = - (run + 1);
	 		  }


	 		  return best4CFit;
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
		  //return combinedList;

	  }

	  void GetNotCombinedList(RhoCandList combinedList, RhoCandList * candList){
		  for (int j=0; j<combinedList.GetLength(); j++){
			  RhoCandidate * combinedCand = combinedList[j];
			  candList->Remove(combinedCand);
		  }
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


	  void qaVtxDiff(TString pre="", RhoCandidate * c, RhoTuple * n){

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

	  void qaMomRes(TString pre="", RhoCandidate * c, RhoTuple * n){

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


	  void numberOfHitsInSubdetector(TString pre="", RhoCandidate *c, RhoTuple *n){

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


	  void tagNHits(TString pre= "", RhoCandidate *c, RhoTuple *n){

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

}
