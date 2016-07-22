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
#include "/home/ikp1/puetz/panda/myscripts/common_jenny.cpp"


void number_of_events_table_cc(TString inFile=""){

	//*** In file

	TFile * input = new TFile(inFile, "READ");

	TTree * ntpMC = (TTree*) input->Get("ntpMC");
	TTree * ntpPiMinus = (TTree*) input->Get("ntpPiMinus");
	TTree * ntpPiPlus = (TTree*) input->Get("ntpPiPlus");
	TTree * ntpkaonplus = (TTree*) input->Get("ntpKaonPlus");
	TTree * ntpProton = (TTree*) input->Get("ntpProton");
	TTree * ntpAntiProton = (TTree*) input->Get("ntpAntiProton");
	TTree * ntpLambda0 = (TTree*) input->Get("ntpLambda0");
	TTree * ntpAntiLambda0 = (TTree*) input->Get("ntpAntiLambda0");
	TTree * ntpXiPlus1820 = (TTree*) input->Get("ntpXiPlus1820");
	TTree * ntpXiMinus = (TTree*) input->Get("ntpXiMinus");
	TTree * ntpXiSys = (TTree*) input->Get("ntpXiSys");


	double nevents_mc = ntpMC->GetEntriesFast();
	TString cuts = " McTruthMatch && VtxFit_HowGood==1 && MassFit_prob>0.01";
	TString VtxCut = " McTruthMatch && VtxFit_HowGood==1 & HitTag==1";
	TString cut4c = "McTruthMatch && 4CFit_prob>0.01";


	cout << "particle|   #evts (uncut)|    #evts (ratio in %)|   MC ratio in %|   dp/p in %" << endl;


	//**** PiMinus(Lambda0)
	TH1D * h_piminus_tht_uncut = new TH1D("h_piminus_tht_uncut", "h_piminus_tht", 100, 0,10);
	ntpPiMinus->Project("h_piminus_tht_uncut", "piminus_tht", "McTruthMatch && Mother==3122");
	double piminus_uncut =  h_piminus_tht_uncut->GetEntries();

	TH1D * h_piminus_tht = new TH1D("h_piminus_tht", "h_piminus_tht", 100, 0,10);
	ntpPiMinus->Project("h_piminus_tht", "piminus_tht", "McTruthMatch && piminus_HitTag && Mother==3122");
	int piminus =  h_piminus_tht->GetEntries();

	TH1D * h_piminus_dp = new TH1D("h_piminus_dp", "h_piminus_dp", 250, -0.1,0.1);
	ntpPiMinus->Project("h_piminus_dp", "(piminus_p-piminus_MC_p)/piminus_MC_p", "McTruthMatch && piminus_HitTag && Mother==3122");

	Double_t param[6] = jenny::GetFitParameterDoubleFit(h_piminus_dp, false, 0.02,0.1, true);
	//jenny::CreateDrawAndSaveHistogramDoulbeFit(h_piminus_dp, "","", false, false, false, 0.02,0.1, true);


	double ratio_piminus_cut = piminus/piminus_uncut;
	double ratio_piminus_mc = piminus/nevents_mc;

	cout << "PiMinus(L0)|   " <<  piminus_uncut << "|   " <<  piminus << "(" << ratio_piminus_cut*100 << ")|   " << ratio_piminus_mc*100 << endl;//<< "|   " << param[2]*100 << endl;



	//**** piminus (Xi)
	TH1D * h_piminus2_tht_uncut2 = new TH1D("h_piminus2_tht_uncut2", "h_piminus2_tht", 100, 0,10);
	ntpPiMinus->Project("h_piminus2_tht_uncut2", "piminus_tht", "McTruthMatch && Mother==3312");
	double piminus_uncut2 =  h_piminus2_tht_uncut2->GetEntries();

	TH1D * h_piminus2_tht = new TH1D("h_piminus2_tht", "h_piminus2_tht", 100, 0,10);
	ntpPiMinus->Project("h_piminus2_tht", "piminus_tht", "McTruthMatch && piminus_HitTag && Mother==3312");
	int piminus2 =  h_piminus2_tht->GetEntries();

	TH1D * h_piminus2_dp = new TH1D("h_piminus2_dp", "h_piminus2_dp", 250, -0.1,0.1);
	ntpPiMinus->Project("h_piminus2_dp", "(piminus_p-piminus_MC_p)/piminus_MC_p", "McTruthMatch && piminus_HitTag && Mother==3312");

	Double_t parampip2[6] = jenny::GetFitParameterDoubleFit(h_piminus2_dp, false, 0.02,0.1, true);
	//jenny::CreateDrawAndSaveHistogramDoulbeFit(h_piminus2_dp, "","", false, false, false, 0.02,0.1, true);

	double ratio_piminus_cut2 = piminus2/piminus_uncut2;
	double ratio_piminus_mc2 = piminus2/nevents_mc;

	cout << "PiMinus(Xi)|   " <<  piminus_uncut2 << "|   " <<  piminus2 << "(" << ratio_piminus_cut2*100 << ")|   " << ratio_piminus_mc2*100<< endl;// << "|   " << parampip2[2]*100 << endl;


	//**** PiPlus (AntiLambda0)
	TH1D * h_piplus_tht_uncut = new TH1D("h_piplus_tht_uncut", "h_piplus_tht", 100, 0,10);
	ntpPiPlus->Project("h_piplus_tht_uncut", "piplus_tht", "McTruthMatch && Mother==-3122");
	double piplus_uncut =  h_piplus_tht_uncut->GetEntries();

	TH1D * h_piplus_tht = new TH1D("h_piplus_tht", "h_piplus_tht", 100, 0,10);
	ntpPiPlus->Project("h_piplus_tht", "piplus_tht", "McTruthMatch && piplus_HitTag && Mother==-3122");
	int piplus =  h_piplus_tht->GetEntries();

	TH1D * h_piplus_dp = new TH1D("h_piplus_dp", "h_piplus_dp", 250, -0.1,0.1);
	ntpPiPlus->Project("h_piplus_dp", "(piplus_p-piplus_MC_p)/piplus_MC_p", "McTruthMatch && piplus_HitTag && Mother==-3122");

	Double_t parampip[6] = jenny::GetFitParameterDoubleFit(h_piplus_dp, false, 0.02,0.1, true);
	//jenny::CreateDrawAndSaveHistogramDoulbeFit(h_piplus_dp, "","", false, false, false, 0.02,0.1, true);

	double ratio_piplus_cut = piplus/piplus_uncut;
	double ratio_piplus_mc = piplus/nevents_mc;

	cout << "PiPlus(AL0)|   " <<  piplus_uncut << "|   " <<  piplus << "(" << ratio_piplus_cut*100 << ")|   " << ratio_piplus_mc*100 << endl;//<< "|   " << parampip[2]*100 << endl;




	//**** kaonplus
	TH1D * h_kaonplus_tht_uncut = new TH1D("h_kaonplus_tht_uncut", "h_kaonplus_tht", 100, 0,10);
	ntpkaonplus->Project("h_kaonplus_tht_uncut", "kaonplus_tht", "McTruthMatch && Mother==-23314");
	double kaonplus_uncut =  h_kaonplus_tht_uncut->GetEntries();

	TH1D * h_kaonplus_tht = new TH1D("h_kaonplus_tht", "h_kaonplus_tht", 100, 0,10);
	ntpkaonplus->Project("h_kaonplus_tht", "kaonplus_tht", "McTruthMatch && kaonplus_HitTag && Mother==-23314");
	int kaonplus =  h_kaonplus_tht->GetEntries();

	TH1D * h_kaonplus_dp = new TH1D("h_kaonplus_dp", "h_kaonplus_dp", 250, -0.1,0.1);
	ntpkaonplus->Project("h_kaonplus_dp", "(kaonplus_p-kaonplus_MC_p)/kaonplus_MC_p", "McTruthMatch && kaonplus_HitTag && Mother==-23314");

	Double_t paramk[6] = jenny::GetFitParameterDoubleFit(h_kaonplus_dp, false, 0.02,0.1, true);
	//jenny::CreateDrawAndSaveHistogramDoulbeFit(h_kaonplus_dp, "","", false, false, false, 0.02,0.1, true);

	double ratio_kaonplus_cut = kaonplus/kaonplus_uncut;
	double ratio_kaonplus_mc = kaonplus/nevents_mc;

	cout << "kaonplus|   " <<  kaonplus_uncut << "|   " <<  kaonplus << "(" << ratio_kaonplus_cut*100 << ")|   " << ratio_kaonplus_mc*100 << endl;//<< "|   " << paramk[2]*100 << endl;



	//**** Proton
	TH1D * h_proton_tht_uncut = new TH1D("h_proton_tht_uncut", "h_proton_tht", 100, 0,10);
	ntpProton->Project("h_proton_tht_uncut", "proton_tht", "McTruthMatch && Mother==3122");
	double proton_uncut =  h_proton_tht_uncut->GetEntries();

	TH1D * h_proton_tht = new TH1D("h_proton_tht", "h_proton_tht", 100, 0,10);
	ntpProton->Project("h_proton_tht", "proton_tht", "McTruthMatch && proton_HitTag && Mother==3122");
	int proton =  h_proton_tht->GetEntries();

	TH1D * h_proton_dp = new TH1D("h_proton_dp", "h_proton_dp", 250, -0.1,0.1);
	ntpProton->Project("h_proton_dp", "(proton_p-proton_MC_p)/proton_MC_p", "McTruthMatch && proton_HitTag && Mother==3122");

	Double_t paramProt[6] = jenny::GetFitParameterDoubleFit(h_proton_dp, false, 0.02,0.1, true);
	//jenny::CreateDrawAndSaveHistogramDoulbeFit(h_proton_dp, "","", false, false, false, 0.02,0.1, true);

	double ratio_proton_cut = proton/proton_uncut;
	double ratio_proton_mc = proton/nevents_mc;

	cout << "proton|   " <<  proton_uncut << "|   " <<  proton << "(" << ratio_proton_cut*100 << ")|   " << ratio_proton_mc*100 << endl;//<< "|   " << paramProt[2]*100 << endl;




	//**** AntiProton
	TH1D * h_AntiProton_tht_uncut = new TH1D("h_AntiProton_tht_uncut", "h_AntiProton_tht", 100, 0,10);
	ntpAntiProton->Project("h_AntiProton_tht_uncut", "AntiProton_tht", "McTruthMatch");
	double AntiProton_uncut =  h_AntiProton_tht_uncut->GetEntries();

	TH1D * h_AntiProton_tht = new TH1D("h_AntiProton_tht", "h_AntiProton_tht", 100, 0,10);
	ntpAntiProton->Project("h_AntiProton_tht", "AntiProton_tht", "McTruthMatch && AntiProton_HitTag");
	int AntiProton =  h_AntiProton_tht->GetEntries();

	TH1D * h_AntiProton_dp = new TH1D("h_AntiProton_dp", "h_AntiProton_dp", 250, -0.1,0.1);
	ntpAntiProton->Project("h_AntiProton_dp", "(AntiProton_p-AntiProton_MC_p)/AntiProton_MC_p", "McTruthMatch && AntiProton_HitTag");

	Double_t paramAProt[6] = jenny::GetFitParameterDoubleFit(h_AntiProton_dp, false, 0.02,0.1, true);
	//jenny::CreateDrawAndSaveHistogramDoulbeFit(h_AntiProton_dp, "","", false, false, false, 0.02,0.1, true);

	double ratio_AntiProton_cut = AntiProton/AntiProton_uncut;
	double ratio_AntiProton_mc = AntiProton/nevents_mc;

	cout << "AntiProton|   " <<  AntiProton_uncut << "|   " <<  AntiProton << "(" << ratio_AntiProton_cut*100 << ")|   " << ratio_AntiProton_mc*100 << endl;//<< "|   " << paramAProt[2]*100 << endl;



	//**** lambda0
	TH1D * h_Lambda0_tht_uncut = new TH1D("h_Lambda0_tht_uncut", "h_Lambda0_tht", 100, 0,10);
	ntpLambda0->Project("h_Lambda0_tht_uncut", "Lambda0_tht", "McTruthMatch & HitTag==1");
	double Lambda0_uncut =  h_Lambda0_tht_uncut->GetEntries();

	TH1D * h_Lambda0_tht = new TH1D("h_Lambda0_tht", "h_Lambda0_tht", 100, 0,10);
	ntpLambda0->Project("h_Lambda0_tht", "Lambda0_tht", "HitTag && "+cuts);
	int lambda0 =  h_Lambda0_tht->GetEntries();

	TH1D * h_Lambda0_dp = new TH1D("h_Lambda0_dp", "h_Lambda0_dp", 250, -0.1,0.1);
	ntpLambda0->Project("h_Lambda0_dp", "(Lambda0_p-McTruth_p)/McTruth_p", "HitTag && "+cuts );

	Double_t paraml0[6] = jenny::GetFitParameterDoubleFit(h_Lambda0_dp, false, 0.02,0.1, true);
	//jenny::CreateDrawAndSaveHistogramDoulbeFit(h_Lambda0_dp, "","", false, false, false, 0.02,0.1, true);

	double ratio_Lambda0_cut = lambda0/Lambda0_uncut;
	double ratio_Lambda0_mc = lambda0/nevents_mc;

	cout << "lambda0|   " <<  Lambda0_uncut << "|   " <<  lambda0 << "(" << ratio_Lambda0_cut*100 << ")|   " << ratio_Lambda0_mc*100 << endl;//<< "|   " << paraml0[2]*100 << endl;



	//**** AntiLambda0
	TH1D * h_antiLambda0_tht_uncut = new TH1D("h_antiLambda0_tht_uncut", "h_antiLambda0_tht", 100, 0,10);
	ntpAntiLambda0->Project("h_antiLambda0_tht_uncut", "antiLambda0_tht", "McTruthMatch & HitTag==1");
	double antiLambda0_uncut =  h_antiLambda0_tht_uncut->GetEntries();

	TH1D * h_antiLambda0_tht = new TH1D("h_antiLambda0_tht", "h_antiLambda0_tht", 100, 0,10);
	ntpAntiLambda0->Project("h_antiLambda0_tht", "antiLambda0_tht", "HitTag && "+cuts);
	int AntiLambda0 =  h_antiLambda0_tht->GetEntries();

	TH1D * h_antiLambda0_dp = new TH1D("h_antiLambda0_dp", "h_antiLambda0_dp", 250, -0.1,0.1);
	ntpAntiLambda0->Project("h_antiLambda0_dp", "(antiLambda0_p-McTruth_p)/McTruth_p", "HitTag && "+cuts);

	Double_t paramAL0[6] = jenny::GetFitParameterDoubleFit(h_antiLambda0_dp, false, 0.02,0.1, true);
	//jenny::CreateDrawAndSaveHistogramDoulbeFit(h_antiLambda0_dp, "","", false, false, false, 0.02,0.1, true);

	double ratio_antiLambda0_cut = AntiLambda0/antiLambda0_uncut;
	double ratio_antiLambda0_mc = AntiLambda0/nevents_mc;

	cout << "AntiLambda0|   " <<  antiLambda0_uncut << "|   " <<  AntiLambda0 << "(" << ratio_antiLambda0_cut*100 << ")|   " << ratio_antiLambda0_mc*100<< endl;// << "|   " << paramAL0[2]*100 << endl;


	//**** XiPlus1820
	TH1D * h_xiplus_tht_uncut = new TH1D("h_xiplus_tht_uncut", "h_xiplus_tht", 100, 0,10);
	ntpXiPlus1820->Project("h_xiplus_tht_uncut", "xiplus_tht", "McTruthMatch & HitTag==1");
	double xiplus_uncut =  h_xiplus_tht_uncut->GetEntries();

	TH1D * h_xiplus_tht = new TH1D("h_xiplus_tht", "h_xiplus_tht", 100, 0,10);
	ntpXiPlus1820->Project("h_xiplus_tht", "xiplus_tht", VtxCut);
	int XiPlus =  h_xiplus_tht->GetEntries();

	TH1D * h_xiplus_dp = new TH1D("h_xiplus_dp", "h_xiplus_dp", 250, -0.1,0.1);
	ntpXiPlus1820->Project("h_xiplus_dp", "(xiplus_p-MCTruth_p)/MCTruth_p", VtxCut);

	Double_t paramxip[6] = jenny::GetFitParameterDoubleFit(h_xiplus_dp, false, 0.02,0.1, true);
	//jenny::CreateDrawAndSaveHistogramDoulbeFit(h_xiplus_dp, "","", false, false, false, 0.02,0.1, true);

	double ratio_xiplus_cut = XiPlus/xiplus_uncut;
	double ratio_XiPlus_McTruth = XiPlus/nevents_mc;

	cout << "XiPlus1820|   " <<  xiplus_uncut << "|   " <<  XiPlus << "(" << ratio_xiplus_cut*100 << ")|   " << ratio_XiPlus_McTruth*100<< endl;// << "|   " << paramxip[2]*100 << endl;


	//**** XiMinus
	TH1D * h_XiMinus_tht_uncut = new TH1D("h_XiMinus_tht_uncut", "h_XiMinus_tht", 100, 0,10);
	ntpXiMinus->Project("h_XiMinus_tht_uncut", "VtxFit_tht", "McTruthMatch & HitTag==1");
	double XiMinus_uncut =  h_XiMinus_tht_uncut->GetEntries();

	TH1D * h_XiMinus_tht = new TH1D("h_XiMinus_tht", "h_XiMinus_tht", 100, 0,10);
	ntpXiMinus->Project("h_XiMinus_tht", "VtxFit_tht", cuts+"& HitTag==1");
	int XiMinus1820 =  h_XiMinus_tht->GetEntries();

	TH1D * h_xiMinus_dp = new TH1D("h_xiMinus_dp", "h_xiMinus_dp", 250, -0.1,0.1);
	ntpXiMinus->Project("h_xiMinus_dp", "VtxFit_p-MCTruth_p");

	//jenny::CreateDrawAndSaveHistogramDoulbeFit(h_xiMinus_dp, " ", " ", false, false, false, 0.02 , 0.1, true);
	Double_t paramxim[6] = jenny::GetFitParameterDoubleFit(h_xiMinus_dp, false, 0.02,0.1, true);

	double ratio_XiMinus_cut = XiMinus1820/XiMinus_uncut;
	double ratio_XiMinus_mc = XiMinus1820/nevents_mc;

	cout << "XiMinus|   " <<  XiMinus_uncut << "|   " <<  XiMinus1820 << "(" << ratio_XiMinus_cut*100 << ")|   " << ratio_XiMinus_mc*100 << endl;// "|   " << paramxim[2]*100 << endl;


	//**** XiSys
	TH1D * h_XiSys_tht_uncut = new TH1D("h_XiSys_tht_uncut", "h_XiSys_tht", 100, 0,10);
	ntpXiSys->Project("h_XiSys_tht_uncut", "XiSys_tht", "McTruthMatch");
	double XiSys_uncut =  h_XiSys_tht_uncut->GetEntries();

	TH1D * h_XiSys_tht = new TH1D("h_XiSys_tht", "h_XiSys_tht", 100, 0,10);
	ntpXiSys->Project("h_XiSys_tht", "XiSys_tht", cut4c);
	int XiSys =  h_XiSys_tht->GetEntries();

	TH1D * h_XiSys_dp = new TH1D("h_XiSys_dp", "h_XiSys_dp", 250, -0.1,0.1);
	ntpXiSys->Project("h_XiSys_dp", "(XiSys_p-McTruth_p)/McTruth_p", cut4c);

	Double_t paramxisys[6] = jenny::GetFitParameterDoubleFit(h_XiSys_dp, false, 0.02,0.1, true);


	double ratio_XiSys_cut = XiSys/XiSys_uncut;
	double ratio_XiSys_mc = XiSys/nevents_mc;

	cout << "XiSys|   " <<  XiSys_uncut << "|   " <<  XiSys << "(" << ratio_XiSys_cut*100 << ")|   " << ratio_XiSys_mc*100<< endl;// << "|   " << paramxisys[2]*100 << endl;


}
