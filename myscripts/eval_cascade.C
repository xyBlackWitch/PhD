/**
 * @file eval_cascade.C
 * @mainpage eval_cascade.C Get Data from analysis file and creates and saves different histograms
 *
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Get Data from analysis file and create and save different histograms
 * @details This files get the data from the analysis of
 * pbar p -> Xi+ Xi-
 * 			 |	  |-> Lambda0 + Pi-
 * 			 |			|-> p + Pi-
 * 			 |
 * 			 |->  AntiLambda0 + Pi+
 * 			 		|-> pbar + Pi+
 * 	 Different useful histogramms are created and saved.
 */

class RhoTuple;
class RhoCandidate;

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "common_jenny.cpp"


void CreateDrawAndSaveNHistograms(TH1* &h1, TH1* &h2, TString leg1="", TString leg2="", TString outputdir, TString outputname, bool saveoutput, bool close){
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


	TCanvas * canvas = new TCanvas("c_"+name, title, 0,0,800,500);
	gStyle->SetOptStat(0);
	h1->Draw();
	h2->Draw("SAME");
	legend->Draw();


	if (saveoutput){
		canvas->Print(outputdir + "root-files/" + outputname + ".root");
		canvas->Print(outputdir + "png-files/" + outputname + ".png");
	}

	if (close) canvas->Close();

}

void CreateDrawAndSaveNHistograms(TH1* &h1, TH1* &h2, TH1* &h3, TString leg1="", TString leg2="", TString leg3="", TString outputdir, TString outputname, bool saveoutput, bool close){
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


	TCanvas * canvas = new TCanvas("c_"+name, title, 0,0,800,500);
	gStyle->SetOptStat(0);
	h1->Draw();
	h2->Draw("SAME");
	h3->Draw("SAME");
	legend->Draw();


	if (saveoutput){
		canvas->Print(outputdir + "root-files/" + outputname + ".root");
		canvas->Print(outputdir + "png-files/" + outputname + ".png");
	}

	if (close) canvas->Close();

}



void eval_cascade(bool save=true, bool close=false){

	//Input file
	TString inPath = "/private/puetz/mysimulations/analysis/pbarp_Xiplus_Ximinus/idealtracking/1000_events/momentum_2.7/run2/";
	TFile * inFile = new TFile(inPath+"output_ana.root", "READ");

	//Output file
	TString outPath = "";

	//Get NTuple
	TTree * ntpMC = (TTree*) inFile->Get("ntpMC");
	TTree * ntpPiMinus = (TTree*) inFile->Get("ntpPiMinus");
	TTree * ntpPiPlus = (TTree*) inFile->Get("ntpPiPlus");
	TTree * ntpProton = (TTree*) inFile->Get("ntpProton");
	TTree * ntpAntiProton= (TTree*) inFile->Get("ntpAntiProton");
	TTree * ntpLambda0 = (TTree*) inFile->Get("ntpLambda0");
	TTree * ntpAntiLambda0= (TTree*) inFile->Get("ntpAntiLambda0");
	TTree * ntpXiPlus = (TTree*) inFile->Get("ntpXiPlus");
	TTree * ntpXiMinus = (TTree*) inFile->Get("ntpXiMinus");


	//************************** Projections for PiMinus from Xi- decay ********************************************

	//Create histograms
	TH1D * h_piminus_tht = new TH1D("h_piminus_tht", "#Theta distribution for #pi^{-}(#Xi^{-}); #theta/rad; counts" ,100,-1,1);
	ntpPiMinus->Project("h_piminus_tht", "piminus_tht", "McTruthMatch==1 && Mother==3312");

	TH1D * h_piminus_pz = new TH1D("h_piminus_pz", "distribution of longitudinal momentum for #pi^{-}(#Xi^{-}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiMinus->Project("h_piminus_pz", "piminus_pz", "McTruthMatch==1 && Mother==3312");

	TH1D * h_piminus_npart = new TH1D("h_piminus_npart", "number of reconstructed #pi^{-}(#Xi^{-}) per Event; reco. #pi^{-} per Event; counts", 7,1,7);
	ntpPiMinus->Project("h_piminus_npart", "ncand", "McTruthMatch==1 && Mother==3312");

	TH2D * h_piminus_pt_vs_pz = new TH2D("h_piminus_pt_vs_pz", "p_{t} vs p_{z} for #pi^{-}(#Xi^{-}); p_{z}/Gev/c; p_{t}/Gev/c", 100,-0.1,0.7,100,0,0.2);
	ntpPiMinus->Project("h_piminus_pt_vs_pz", "piminus_pt:piminus_pz", "McTruthMatch==1 && Mother==3312");

	TH2D * h_piminus_costht_vs_pz = new TH2D("h_piminus_costht_vs_pz", "cos(#theta) vs p_{z} for #pi^{-}(#Xi^{-}); p_{z}/Gev/c; cos(#theta)", 100,-0.1,0.7,100,-1.,1.);
	ntpPiMinus->Project("h_piminus_costht_vs_pz", "PiMinus_CosTheta:piminus_pz", "McTruthMatch==1 && Mother==3312");


	//historgrams for efficiency
	TH1D * h_piminus_MC_tht = new TH1D("h_piminus_MC_tht", "cos(#Theta) distribution for #pi^{-}_{MC}(#Xi^{-}); cos(#theta); counts" ,100,-1,1);
	ntpPiMinus->Project("h_piminus_MC_tht", "cos(piminus_MC_tht)", "McTruthMatch==1 && Mother==3312");

	TH1D * h_piminus_MC_pz = new TH1D("h_piminus_MC_pz", "distribution of longitudinal momentum for #pi^{-}_{MC}(#Xi^{-}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiMinus->Project("h_piminus_MC_pz", "piminus_MC_pz", "McTruthMatch==1 && Mother==3312");

	TH1D * h_MC_tht = new TH1D("h_MC_tht", "#Theta_{MC} distribution for #pi^{-}(#Xi^{-}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_MC_tht", "cos(tht)", "moth==1 && pdg==-211");

	TH1D * h_MC_pz = new TH1D("h_MC_pz", "distribution of longitudinal momentum (MC) for #pi^{-}(#Xi^{-}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_MC_pz", "pz", "moth==1 && pdg==-211");

	TH1D * h_eff_pz = new TH1D("h_eff_pz", "reconstruction efficiency for #pi^{-} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_eff_pz->Divide(h_piminus_MC_pz, h_MC_pz);
	h_eff_pz->SetStats(0);
	h_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_eff_tht = new TH1D("h_eff_tht", "reconstruction efficiency for #pi^{-} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_eff_tht->Divide(h_piminus_MC_tht, h_MC_tht);
	h_eff_tht->SetStats(0);
	h_eff_tht->GetYaxis()->SetRangeUser(0,1.1);




	//	save histograms
//	jenny::CreateDrawAndSaveHistogram(h_piminus_MC_tht, outPath, "h_piminus_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminus_pz, outPath, "h_piminus_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminus_npart, outPath, "h_piminus_npart", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminus_pt_vs_pz, outPath, "h_piminus_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminus_costht_vs_pz, outPath, "h_piminus_costht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_eff_tht, outPath, "h_eff_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_eff_pz, outPath, "h_eff_pz", save, close);



	//************************** Projections for piplus from Xi- decay ********************************************

	//Create histograms
	TH1D * h_piplus_tht = new TH1D("h_piplus_tht", "cos(#Theta) distribution for #pi^{+}(#Xi^{+}); cos(#theta); counts" ,100,-1,1);
	ntpPiPlus->Project("h_piplus_tht", "PiPlus_tht", "McTruthMatch==1 && Mother==-3312");

	TH1D * h_piplus_pz = new TH1D("h_piplus_pz", "distribution of longitudinal momentum for #pi^{+}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiPlus->Project("h_piplus_pz", "PiPlus_pz", "McTruthMatch==1 && Mother==-3312");

	TH1D * h_piplus_npart = new TH1D("h_piplus_npart", "number of reconstructed #pi^{+}(#Xi^{+}) per Event; reco. #pi^{+} per Event; counts", 7,1,7);
	ntpPiPlus->Project("h_piplus_npart", "ncand", "McTruthMatch==1 && Mother==-3312");

	TH2D * h_piplus_pt_vs_pz = new TH2D("h_piplus_pt_vs_pz", "p_{t} vs p_{z} for #pi^{+}(#Xi^{+}); p_{z}/Gev/c; p_{t}/Gev/c", 100,-0.1,0.7,100,0,0.2);
	ntpPiPlus->Project("h_piplus_pt_vs_pz", "PiPlus_pt:PiPlus_pz", "McTruthMatch==1 && Mother==-3312");

	TH2D * h_piplus_costht_vs_pz = new TH2D("h_piplus_costht_vs_pz", "cos(#theta) vs p_{z} for #pi^{+}(#Xi^{+}); p_{z}/Gev/c; cos(#theta)", 100,-0.1,0.7,100,-1.,1.);
	ntpPiPlus->Project("h_piplus_costht_vs_pz", "PiPlus_CosTheta:PiPlus_pz", "McTruthMatch==1 && Mother==-3312");


	//historgrams for efficiency
	TH1D * h_piplus_McTruth_tht = new TH1D("h_piplus_McTruth_tht", "#Theta distribution for #pi^{+}_{MC}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpPiPlus->Project("h_piplus_McTruth_tht", "cos(PiPlus_MC_tht)", "McTruthMatch==1 && Mother==-3312");

	TH1D * h_piplus_McTruth_pz = new TH1D("h_piplus_McTruth_pz", "distribution of longitudinal momentum for #pi^{+}_{MC}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiPlus->Project("h_piplus_McTruth_pz", "PiPlus_MC_pz", "McTruthMatch==1 && Mother==-3312");

	TH1D * h_piplus_MC_tht = new TH1D("h_piplus_MC_tht", "#Theta_{MC} distribution for #pi^{+}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_piplus_MC_tht", "cos(tht)", "moth==2 && pdg==211");

	TH1D * h_piplus_MC_pz = new TH1D("h_piplus_MC_pz", "distribution of longitudinal momentum (MC) for #pi^{+}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_piplus_MC_pz", "pz", "moth==2 && pdg==211");

	TH1D * h_piplus_eff_pz = new TH1D("h_piplus_eff_pz", "reconstruction efficiency for #pi^{+} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_piplus_eff_pz->Divide(h_piplus_McTruth_pz, h_piplus_MC_pz);
	h_piplus_eff_pz->SetStats(0);
	h_piplus_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_piplus_eff_tht = new TH1D("h_piplus_eff_tht", "reconstruction efficiency for #pi^{+} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_piplus_eff_tht->Divide(h_piplus_McTruth_tht, h_piplus_MC_tht);
	h_piplus_eff_tht->SetStats(0);
	h_piplus_eff_tht->GetYaxis()->SetRangeUser(0,1.1);




	//save histograms
//	jenny::CreateDrawAndSaveHistogram(h_piplus_MC_tht, outPath, "h_piplus_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_pz, outPath, "h_piplus_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_npart, outPath, "h_piplus_npart", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_pt_vs_pz, outPath, "h_piplus_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_costht_vs_pz, outPath, "h_piplus_costht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_eff_tht, outPath, "h_piplus_eff_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_eff_pz, outPath, "h_piplus_eff_pz", save, close);


	//*********************** plots for comparison of pi+ and pi- *********************************************

//	CreateDrawAndSaveNHistograms(h_piminus_npart, h_piplus_npart, "#pi^{-}","#pi^{+}",  outPath, "h_piplus_piminus_npart", save, close);



	//************************ projections for lambda0 *******************************************************

	TH1D * h_lambda0_tht = new TH1D("h_lambda0_tht", "#theta distribution for #Lambda^{0}; #theta/rad; counts", 100,0,1);
	ntpLambda0->Project("h_lambda0_tht", "Lambda0_tht", "McTruthMatch==1");

	TH1D * h_lambda0_pz = new TH1D("h_lambda0_pz", "distribution of longitudinal momentum for #Lambda^{0}; p_{z}/GeV/c; counts", 100,0,2.7);
	ntpLambda0->Project("h_lambda0_pz", "Lambda0_pz", "McTruthMatch==1");


	TString cut1 = "VtxFit_prob>0.01";
	TString cut2 = cut1 + "&& fMass_Prob>0.01";


	//mass resolution with different cuts
	TH1D * h_lambda0_mass = new TH1D("h_lambda0_mass", "mass resolution of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpLambda0->Project("h_lambda0_mass", "(Lambda0_m-MCTruth_m)/1.115683", "McTruthMatch==1");

	TH1D * h_lambda0_mass_cut1 = new TH1D("h_lambda0_mass_cut1", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpLambda0->Project("h_lambda0_mass_cut1", "(Lambda0_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut1);

	TH1D * h_lambda0_mass_cut2 = new TH1D("h_lambda0_mass_cut2", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpLambda0->Project("h_lambda0_mass_cut2", "(Lambda0_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut2);


	//momentum resolution for different cuts

	TH1D * h_lambda0_p_diff = new TH1D("h_lambda0_p_diff", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpLambda0->Project("h_lambda0_p_diff", "Lambda0_p-MCTruth_p", "McTruthMatch==1");

	TH1D * h_lambda0_p_cut1 = new TH1D("h_lambda0_p_cut1", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpLambda0->Project("h_lambda0_p_cut1", "Lambda0_p-MCTruth_p", "McTruthMatch==1&&"+cut1);

	TH1D * h_lambda0_p_cut2 = new TH1D("h_lambda0_p_cut2", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpLambda0->Project("h_lambda0_p_cut2", "Lambda0_p-MCTruth_p", "McTruthMatch==1&&"+cut2);

//	TH1D * h_lambda0_p = new TH1D("h_lambda0_p", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,0,3);
//	ntpLambda0->Project("h_lambda0_p", "Lambda0_p", "McTruthMatch==1");
//
//
//	TH1D * h_lambda0_p_mc = new TH1D("h_lambda0_p_mc", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,0,3);
//	ntpLambda0->Project("h_lambda0_p_mc", "MCTruth_p", "McTruthMatch==1");

//	TH1D * h_lambda0_p_diff2 = new TH1D("h_lambda0_p_diff2", "test", 2,0,2);
//
//
//	float p_reco, p_mc;
//	UChar_t truthmatch;
//
//	ntpLambda0->SetBranchAddress("MCTruth_p", &p_mc);
//	ntpLambda0->SetBranchAddress("Lambda0_p", &p_reco);
//	ntpLambda0->SetBranchAddress("McTruthMatch", &truthmatch);
//
//
//	int nentries = ntpLambda0->GetEntriesFast();
//
//	for (int i=-1; i<nentries; ++i){
//		float truth = (float) truthmatch;
//		ntpLambda0->GetEntry(i);
////		cout << truth << endl;
////		cout << p_reco << endl;
////		if (truth==1)
//		h_lambda0_p_diff2->Fill(truth);
//	}


	//Vertex resolution

	TH1D * h_lambda0_vtx_resx = new TH1D("h_lambda0_vtx_resx", "vertex resolution for #Lambda^{0} (x coordinate); #Delta x[cm]; counts", 30,-3,3);
	ntpLambda0->Project("h_lambda0_vtx_resx", "VtxFit_vx-MCTruth_vx", "McTruthMatch==1");

	TH1D * h_lambda0_vtx_resy = new TH1D("h_lambda0_vtx_resy", "vertex resolution for #Lambda^{0} (y coordinate); #Delta y[cm]; counts", 30,-3,3);
	ntpLambda0->Project("h_lambda0_vtx_resy", "VtxFit_vy-MCTruth_vy", "McTruthMatch==1");

	TH1D * h_lambda0_vtx_resz = new TH1D("h_lambda0_vtx_resz", "vertex resolution for #Lambda^{0} (z coordinate); #Delta z[cm]; counts", 30,-3,3);
	ntpLambda0->Project("h_lambda0_vtx_resz", "VtxFit_vz-MCTruth_vz", "McTruthMatch==1");


	//2d projections

	TH2D * h_lambda0_pt_vs_pz = new TH2D("h_lambda0_pt_vs_pz", "p_{t} vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; p_{t}/GeV/c", 100,0,2.7, 100,0,0.4);
	ntpLambda0->Project("h_lambda0_pt_vs_pz", "Lambda0_pt:Lambda0_pz", "McTruthMatch==1");

	TH2D * h_lambda0_tht_vs_pz = new TH2D("h_lambda0_tht_vs_pz", "cos(#theta) vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; cos(#theta)", 100,0,2.7, 100,-1,1);
	ntpLambda0->Project("h_lambda0_tht_vs_pz", "cos(Lambda0_tht):Lambda0_pz", "McTruthMatch==1");

	//reconstruction efficiency

	TH1D * h_piplus_McTruth_tht = new TH1D("h_piplus_McTruth_tht", "#Theta distribution for #pi^{+}_{MC}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpPiPlus->Project("h_piplus_McTruth_tht", "cos(PiPlus_MC_tht)", "McTruthMatch==1 && Mother==-3312");

	TH1D * h_piplus_McTruth_pz = new TH1D("h_piplus_McTruth_pz", "distribution of longitudinal momentum for #pi^{+}_{MC}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiPlus->Project("h_piplus_McTruth_pz", "PiPlus_MC_pz", "McTruthMatch==1 && Mother==-3312");

	TH1D * h_piplus_MC_tht = new TH1D("h_piplus_MC_tht", "#Theta_{MC} distribution for #pi^{+}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_piplus_MC_tht", "cos(tht)", "moth==2 && pdg==211");

	TH1D * h_piplus_MC_pz = new TH1D("h_piplus_MC_pz", "distribution of longitudinal momentum (MC) for #pi^{+}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_piplus_MC_pz", "pz", "moth==2 && pdg==211");

	TH1D * h_piplus_eff_pz = new TH1D("h_piplus_eff_pz", "reconstruction efficiency for #pi^{+} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_piplus_eff_pz->Divide(h_piplus_McTruth_pz, h_piplus_MC_pz);
	h_piplus_eff_pz->SetStats(0);
	h_piplus_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_piplus_eff_tht = new TH1D("h_piplus_eff_tht", "reconstruction efficiency for #pi^{+} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_piplus_eff_tht->Divide(h_piplus_McTruth_tht, h_piplus_MC_tht);
	h_piplus_eff_tht->SetStats(0);
	h_piplus_eff_tht->GetYaxis()->SetRangeUser(0,1.1);


	// save histograms
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_tht, outPath, "h_lambda0_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_pz, outPath, "h_lambda0_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_mass, outPath, "h_lambda0_massres", save, close);
//	CreateDrawAndSaveNHistograms(h_lambda0_mass, h_lambda0_mass_cut1,h_lambda0_mass_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_lambda0_mass_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_p_diff2, outPath, "h_lambda0_p_diff", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_p_mc, outPath, "h_lambda0", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_pt_vs_pz, outPath, "h_lambda0_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_tht_vs_pz, outPath, "h_lambda0_tht_vs_pz", save, close);
//	CreateDrawAndSaveNHistograms(h_lambda0_p_diff, h_lambda0_p_cut1,h_lambda0_p_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_lambda0_pres_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_vtx_resx, outPath, "h_lambda0_vtx_resx", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_vtx_resy, outPath, "h_lambda0_vtx_resy", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_vtx_resz, outPath, "h_lambda0_vtx_resz", save, close);


}
