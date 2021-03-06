/**
 * @file number_of_events_table.C
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Creates a table with the number of evts before and after cut for each particle.
 * @details Methode Creates a table with the number of evts before and after cut for each particle.. Data is comming from the analysis which is stored in a file called output_ana.root
 * Basically ROOT.
 */

class RhoTuple;

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "../../myscripts/common_jenny.cpp"

#include <fstream>
#include <iostream>

void number_of_events_table_DPM(TString inFile=""){

	std::ifstream filelist(inFile);

	if (!filelist.is_open()){
		cerr << "failed in open." << endl;
		exit();
	}


	//*** In file
//	TChain ntpPiMinus("ntpPiMinus");
//	TFile * input = new TFile(chain, "READ");

	TChain ntpMC ("ntpMC");
	TChain ntpPiMinus ("ntpPiMinus");
	TChain ntpPiPlus ("ntpPiPlus");
	TChain ntpKaonMinus ("ntpKaonMinus");
	TChain ntpProton  ("ntpProton");
	TChain ntpAntiProton  ("ntpAntiProton");
	TChain ntpLambda0  ("ntpLambda0");
	TChain ntpAntiLambda0  ("ntpAntiLambda0");
	TChain ntpXiPlus  ("ntpXiPlus");
	TChain ntpXiMinus1820  ("ntpXiMinus1820");
	TChain ntpXiSys  ("ntpXiSys");

	char line[BUFSIZ];

	int counter=0;


	while (!filelist.eof()){

		filelist.getline(line, sizeof(line));
		if (strcmp(line, "") == 0) continue;


		ntpMC.Add(line);
		ntpPiMinus.Add(line);
		ntpPiPlus.Add(line);
		ntpKaonMinus.Add(line);
		ntpProton.Add(line);
		ntpAntiProton.Add(line);
		ntpLambda0.Add(line);
		ntpAntiLambda0.Add(line);
		ntpXiPlus.Add(line);
		ntpXiMinus1820.Add(line);
		ntpXiSys.Add(line);

		counter+=1;
	}

	double nevents_mc = counter * 2000; //ntpMC.GetEntriesFast();
	TString cuts = "  VtxFit_HowGood==1 && MassFit_HowGood>0";
	TString VtxCut = "   VtxFit_HowGood==1";
	TString cut4c = " 4CFit_prob>0.01";

	cout << nevents_mc << endl;

	cout << "particle|   #evts (uncut)|    #evts (ratio in %)|   MC ratio in %|   dp/p in %" << endl;


	//**** PiMinus
	TH1D * h_piminus_tht_uncut = new TH1D("h_piminus_tht_uncut", "h_piminus_tht", 100, 0,10);
	ntpPiMinus.Project("h_piminus_tht_uncut", "piminus_tht", " Mother==3122");
	double piminus_uncut =  h_piminus_tht_uncut->GetEntries();

	TH1D * h_piminus_tht = new TH1D("h_piminus_tht", "h_piminus_tht", 100, 0,10);
	ntpPiMinus.Project("h_piminus_tht", "piminus_tht", " piminus_HitTag && Mother==3122");
	int piminus =  h_piminus_tht->GetEntries();

	TH1D * h_piminus_dp = new TH1D("h_piminus_dp", "h_piminus_dp", 250, -0.1,0.1);
	ntpPiMinus.Project("h_piminus_dp", "(piminus_p-piminus_MC_p)/piminus_MC_p", "  piminus_HitTag && Mother==3122");

	//Double_t param[6] = jenny::GetFitParameterDoubleFit(h_piminus_dp, false, 0.02,0.1, true);

	double ratio_piminus_cut = (piminus_uncut==0)? 0:  piminus/piminus_uncut;
	double ratio_piminus_mc = piminus/nevents_mc;

	cout << "PiMinus|   " <<  piminus_uncut << "|   " <<  piminus << "(" << ratio_piminus_cut*100 << ")|   " << ratio_piminus_mc*100 << endl;// "|   " << param[2]*100 << endl;


	//**** PiPlus (AntiLambda0)
	TH1D * h_piplus_tht_uncut = new TH1D("h_piplus_tht_uncut", "h_piplus_tht", 100, 0,10);
	ntpPiPlus.Project("h_piplus_tht_uncut", "piplus_tht", "  Mother==-3122");
	double piplus_uncut =  h_piplus_tht_uncut->GetEntries();

	TH1D * h_piplus_tht = new TH1D("h_piplus_tht", "h_piplus_tht", 100, 0,10);
	ntpPiPlus.Project("h_piplus_tht", "piplus_tht", "  piplus_HitTag && Mother==-3122");
	int piplus =  h_piplus_tht->GetEntries();

	double ratio_piplus_cut = piplus/piplus_uncut;
	double ratio_piplus_mc = piplus/nevents_mc;

	TH1D * h_piplus_dp = new TH1D("h_piplus_dp", "h_piplus_dp", 250, -0.1,0.1);
	ntpPiPlus.Project("h_piplus_dp", "(piplus_p-piplus_MC_p)/piplus_MC_p", "  piplus_HitTag && Mother==-3122");

	//Double_t parampip[6] = jenny::GetFitParameterDoubleFit(h_piplus_dp, false, 0.02,0.1, true);

	cout << "PiPlus(AL0)|   " <<  piplus_uncut << "|   " <<  piplus << "(" << ratio_piplus_cut*100 << ")|   " << ratio_piplus_mc*100 << endl;//"|   " << parampip[2]*100 << endl;



	//**** piplus (Xi+)
	TH1D * h_piplus2_tht_uncut = new TH1D("h_piplus2_tht_uncut", "h_piplus2_tht", 100, 0,10);
	ntpPiPlus.Project("h_piplus2_tht_uncut", "piplus_tht", "  Mother==-3312");
	double piplus2_uncut =  h_piplus2_tht_uncut->GetEntries();


	TH1D * h_piplus2_tht = new TH1D("h_piplus2_tht", "h_piplus2_tht", 100, 0,10);
	ntpPiPlus.Project("h_piplus2_tht", "piplus_tht", "  piplus_HitTag && Mother==-3312");
	int piplus2 =  h_piplus2_tht->GetEntries();

	TH1D * h_piplus2_dp = new TH1D("h_piplus2_dp", "h_piplus_dp", 250, -0.1,0.1);
	ntpPiPlus.Project("h_piplus2_dp", "(piplus_p-piplus_MC_p)/piplus_MC_p", "  piplus_HitTag && Mother==-3312");

	//Double_t parampip2[6] = jenny::GetFitParameterDoubleFit(h_piplus2_dp, false, 0.02,0.1, true);

	if(piplus2_uncut==0) double ratio_piplus2_cut = 0;
	else{
		double ratio_piplus2_cut = piplus2/piplus2_uncut;
	}
	double ratio_piplus2_mc = piplus2/nevents_mc;
	

	cout << "PiPlus(Xi+)|" <<  piplus2_uncut << "|" <<  piplus2 << "(" << ratio_piplus2_cut*100 << ")|" << ratio_piplus2_mc*100 << endl;//"|   " << parampip2[2]*100 << endl;


	//**** kaonMinus
	TH1D * h_kaonMinus_tht_uncut = new TH1D("h_kaonMinus_tht_uncut", "h_kaonMinus_tht", 100, 0,10);
	ntpKaonMinus.Project("h_kaonMinus_tht_uncut", "kaonminus_tht", "  Mother==23314");
	double kaonMinus_uncut =  h_kaonMinus_tht_uncut->GetEntries();

	TH1D * h_kaonMinus_tht = new TH1D("h_kaonMinus_tht", "h_kaonMinus_tht", 100, 0,10);
	ntpKaonMinus.Project("h_kaonMinus_tht", "kaonminus_tht", "  kaonminus_HitTag && Mother==23314");
	int kaonMinus =  h_kaonMinus_tht->GetEntries();

	TH1D * h_kaonminus_dp = new TH1D("h_kaonminus_dp", "h_kaonminus_dp", 250, -0.1,0.1);
	ntpKaonMinus.Project("h_kaonminus_dp", "(kaonminus_p-kaonminus_MC_p)/kaonminus_MC_p", " kaonminus_HitTag && Mother==23314");

	//Double_t paramk[6] = jenny::GetFitParameterDoubleFit(h_kaonminus_dp, false, 0.02,0.1, true);

	if(kaonMinus_uncut==0) double ratio_kaonMinus_cut = 0;
	else{double ratio_kaonMinus_cut = kaonMinus/kaonMinus_uncut;
	}
	double ratio_kaonMinus_mc = kaonMinus/nevents_mc;

	cout << "kaonMinus|   " <<  kaonMinus_uncut << "|   " <<  kaonMinus << "(" << ratio_kaonMinus_cut*100 << ")|   " << ratio_kaonMinus_mc*100 << endl;//"|   " << paramk[2]*100 << endl;


	//**** Proton
	TH1D * h_proton_tht_uncut = new TH1D("h_proton_tht_uncut", "h_proton_tht", 100, 0,10);
	ntpProton.Project("h_proton_tht_uncut", "proton_tht", "  Mother==3122");
	double proton_uncut =  h_proton_tht_uncut->GetEntries();

	TH1D * h_proton_tht = new TH1D("h_proton_tht", "h_proton_tht", 100, 0,10);
	ntpProton.Project("h_proton_tht", "proton_tht", "  proton_HitTag && Mother==3122");
	int proton =  h_proton_tht->GetEntries();

	TH1D * h_proton_dp = new TH1D("h_proton_dp", "h_proton_dp", 250, -0.1,0.1);
	ntpProton.Project("h_proton_dp", "(proton_p-proton_MC_p)/proton_MC_p", "  proton_HitTag && Mother==3122");

	//Double_t paramProt[6] = jenny::GetFitParameterDoubleFit(h_proton_dp, false, 0.02,0.1, true);
	if(proton_uncut==0) double ratio_proton_cut = 0;
	else{double ratio_proton_cut = proton/proton_uncut;
	}
	double ratio_proton_mc = proton/nevents_mc;

	cout << "proton|   " <<  proton_uncut << "|   " <<  proton << "(" << ratio_proton_cut*100 << ")|   " << ratio_proton_mc*100 <<  endl;//"|   " << paramProt[2]*100 << endl;




	//**** AntiProton
	TH1D * h_AntiProton_tht_uncut = new TH1D("h_AntiProton_tht_uncut", "h_AntiProton_tht", 100, 0,10);
	ntpAntiProton.Project("h_AntiProton_tht_uncut", "AntiProton_tht", "");
	double AntiProton_uncut =  h_AntiProton_tht_uncut->GetEntries();

	TH1D * h_AntiProton_tht = new TH1D("h_AntiProton_tht", "h_AntiProton_tht", 100, 0,10);
	ntpAntiProton.Project("h_AntiProton_tht", "AntiProton_tht", "  AntiProton_HitTag");
	int AntiProton =  h_AntiProton_tht->GetEntries();

	if(AntiProton_uncut==0) double ratio_AntiProton_cut = 0;
	else{double ratio_AntiProton_cut = AntiProton/AntiProton_uncut;
	}

	double ratio_AntiProton_mc = AntiProton/nevents_mc;

	cout << "AntiProton|   " <<  AntiProton_uncut << "|   " <<  AntiProton << "(" << ratio_AntiProton_cut*100 << ")|   " << ratio_AntiProton_mc*100 << endl;//"|   " << paramProt[2]*100 << endl;



	//**** lambda0
	TH1D * h_Lambda0_tht_uncut = new TH1D("h_Lambda0_tht_uncut", "h_Lambda0_tht", 100, 0,10);
	ntpLambda0.Project("h_Lambda0_tht_uncut", "Lambda0_tht", "");
	double Lambda0_uncut =  h_Lambda0_tht_uncut->GetEntries();

	TH1D * h_Lambda0_tht = new TH1D("h_Lambda0_tht", "h_Lambda0_tht", 100, 0,10);
	ntpLambda0.Project("h_Lambda0_tht", "Lambda0_tht", "HitTag && "+cuts);
	int lambda0 =  h_Lambda0_tht->GetEntries();
	TH1D * h_Lambda0_dp = new TH1D("h_Lambda0_dp", "h_Lambda0_dp", 250, -0.1,0.1);
	ntpLambda0.Project("h_Lambda0_dp", "(Lambda0_p-McTruth_p)/McTruth_p", "HitTag && "+cuts );

	//Double_t paraml0[6] = jenny::GetFitParameterDoubleFit(h_Lambda0_dp, false, 0.02,0.1, true);

	double ratio_Lambda0_cut = lambda0/Lambda0_uncut;
	double ratio_Lambda0_mc = lambda0/nevents_mc;

	cout << "lambda0|   " <<  Lambda0_uncut << "|   " <<  lambda0 << "(" << ratio_Lambda0_cut*100 << ")|   " << ratio_Lambda0_mc*100 << endl;//"|   " << paraml0[2]*100 << endl;


	//**** AntiLambda0
	TH1D * h_antiLambda0_tht_uncut = new TH1D("h_antiLambda0_tht_uncut", "h_antiLambda0_tht", 100, 0,10);
	ntpAntiLambda0.Project("h_antiLambda0_tht_uncut", "antiLambda0_tht", "");
	double antiLambda0_uncut =  h_antiLambda0_tht_uncut->GetEntries();

	TH1D * h_antiLambda0_tht = new TH1D("h_antiLambda0_tht", "h_antiLambda0_tht", 100, 0,10);
	ntpAntiLambda0.Project("h_antiLambda0_tht", "antiLambda0_tht", "HitTag && "+cuts);
	int AntiLambda0 =  h_antiLambda0_tht->GetEntries();

	TH1D * h_antiLambda0_dp = new TH1D("h_antiLambda0_dp", "h_antiLambda0_dp", 250, -0.1,0.1);
	ntpAntiLambda0.Project("h_antiLambda0_dp", "(antiLambda0_p-McTruth_p)/McTruth_p", "HitTag && "+cuts);

	//Double_t paramAL0[6] = jenny::GetFitParameterDoubleFit(h_antiLambda0_dp, false, 0.02,0.1, true);

	double ratio_antiLambda0_cut = AntiLambda0/antiLambda0_uncut;
	double ratio_antiLambda0_mc = AntiLambda0/nevents_mc;

	cout << "AntiLambda0|   " <<  antiLambda0_uncut << "|   " <<  AntiLambda0 << "(" << ratio_antiLambda0_cut*100 << ")|   " << ratio_antiLambda0_mc*100 << endl;//"|   " << paramAL0[2]*100 << endl;



	//**** XiPlus
	TH1D * h_xiplus_tht_uncut = new TH1D("h_xiplus_tht_uncut", "h_xiplus_tht", 100, 0,10);
	ntpXiPlus.Project("h_xiplus_tht_uncut", "xiplus_tht", "");
	double xiplus_uncut =  h_xiplus_tht_uncut->GetEntries();

	TH1D * h_xiplus_tht = new TH1D("h_xiplus_tht", "h_xiplus_tht", 100, 0,10);
	ntpXiPlus.Project("h_xiplus_tht", "xiplus_tht", cuts);
	int XiPlus =  h_xiplus_tht->GetEntries();

	TH1D * h_xiplus_dp = new TH1D("h_xiplus_dp", "h_xiplus_dp", 250, -0.1,0.1);
	ntpXiPlus.Project("h_xiplus_dp", "(xiplus_p-MCTruth_p)/MCTruth_p", cuts);

	//Double_t paramxip[6] = jenny::GetFitParameterDoubleFit(h_xiplus_dp, false, 0.02,0.1, true);

	if(xiplus_uncut==0) double ratio_xiplus_cut = 0;
	else{double ratio_xiplus_cut = XiPlus/xiplus_uncut;
	}

	double ratio_XiPlus_McTruth = XiPlus/nevents_mc;

	cout << "XiPlus|   " <<  xiplus_uncut << "|   " <<  XiPlus << "(" << ratio_xiplus_cut*100 << ")|   " << ratio_XiPlus_McTruth*100 << endl;//"|   " << paramxip[2]*100 << endl;

	//**** XiMinus1820
	TH1D * h_XiMinus_tht_uncut = new TH1D("h_XiMinus_tht_uncut", "h_XiMinus_tht", 100, 0,10);
	ntpXiMinus1820.Project("h_XiMinus_tht_uncut", "XiMinus_tht", "");
	double XiMinus_uncut =  h_XiMinus_tht_uncut->GetEntries();

	TH1D * h_XiMinus_tht = new TH1D("h_XiMinus_tht", "h_XiMinus_tht", 100, 0,10);
	ntpXiMinus1820.Project("h_XiMinus_tht", "XiMinus_tht", VtxCut);
	int XiMinus1820 =  h_XiMinus_tht->GetEntries();

	TH1D * h_xiMinus_dp = new TH1D("h_xiMinus_dp", "h_xiMinus_dp", 250, -0.1,0.1);
	ntpXiMinus1820.Project("h_xiMinus_dp", "VtxFit_p-MCTruth_p", VtxCut);

	//	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_xiMinus_dp, " ", " ", false, false, false, 0.02 , 0.1, true);
	//Double_t paramxim[6] = jenny::GetFitParameterDoubleFit(h_xiMinus_dp, false, 0.02,0.1, true);
	
	if(XiMinus_uncut==0) double ratio_XiMinus_cut = 0;
	else{double ratio_XiMinus_cut = XiMinus1820/XiMinus_uncut;
	}

	double ratio_XiMinus_mc = XiMinus1820/nevents_mc;

	cout << "XiMinus1820|   " <<  XiMinus_uncut << "|   " <<  XiMinus1820 << "(" << ratio_XiMinus_cut*100 << ")|   " << ratio_XiMinus_mc*100 << endl;// "|   " << paramxim[2]*100 << endl;


	//**** XiSys
	TH1D * h_XiSys_tht_uncut = new TH1D("h_XiSys_tht_uncut", "h_XiSys_tht", 100, 0,10);
	ntpXiSys.Project("h_XiSys_tht_uncut", "XiSys_tht", "");
	double XiSys_uncut =  h_XiSys_tht_uncut->GetEntries();

	TH1D * h_XiSys_tht = new TH1D("h_XiSys_tht", "h_XiSys_tht", 100, 0,10);
	ntpXiSys.Project("h_XiSys_tht", "XiSys_tht", cut4c);
	int XiSys =  h_XiSys_tht->GetEntries();

	TH1D * h_XiSys_dp = new TH1D("h_XiSys_dp", "h_XiSys_dp", 250, -0.1,0.1);
	ntpXiSys.Project("h_XiSys_dp", "(XiSys_p-McTruth_p)/McTruth_p", cut4c);

	//Double_t paramxisys[6] = jenny::GetFitParameterDoubleFit(h_XiSys_dp, false, 0.02,0.1, true);

	if(XiSys_uncut==0) double ratio_XiSys_cut = 0;
	else{double ratio_XiSys_cut = XiSys/XiSys_uncut;
	}

	double ratio_XiSys_mc = XiSys/nevents_mc;

	cout << "XiSys|   " <<  XiSys_uncut << "|   " <<  XiSys << "(" << ratio_XiSys_cut*100 << ")|   " << ratio_XiSys_mc*100 << endl;//"|   " << paramxisys[2]*100 << endl;


//	TCanvas *c = new TCanvas("c","c", 0,0, 800,500);
//
//	TH1D * reco = new TH1D("reco", "fraction of reconstructed final state particles; particle type; reco efficiency in %", 5,0,5);
//	reco->Fill("#pi^{-}", ratio_piminus_mc*100);
//	reco->Fill("#pi^{+}", (ratio_piplus_mc+ratio_piplus2_mc)/2*100);
//	reco->Fill("K^{-}", ratio_kaonMinus_mc*100);
//	reco->Fill("p", ratio_proton_mc*100);
//	reco->Fill("#bar{p}", ratio_AntiProton_mc*100);
//
//
//	gStyle->SetOptStat(0);
//	reco->GetYaxis()->SetRangeUser(0,100);
//	reco->Draw();

}
