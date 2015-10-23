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

void number_of_events_table(TString inFile=""){

	//*** In file

	TFile * input = new TFile(inFile, "READ");

	TTree * ntpMC = (TTree*) input->Get("ntpMC");
	TTree * ntpPiMinus = (TTree*) input->Get("ntpPiMinus");
	TTree * ntpPiPlus = (TTree*) input->Get("ntpPiPlus");
	TTree * ntpKaonMinus = (TTree*) input->Get("ntpKaonMinus");
	TTree * ntpProton = (TTree*) input->Get("ntpProton");
	TTree * ntpAntiProton = (TTree*) input->Get("ntpAntiProton");
	TTree * ntpLambda0 = (TTree*) input->Get("ntpLambda0");
	TTree * ntpAntiLambda0 = (TTree*) input->Get("ntpAntiLambda0");
	TTree * ntpXiPlus = (TTree*) input->Get("ntpXiPlus");
	TTree * ntpXiMinus1820 = (TTree*) input->Get("ntpXiMinus1820");
	TTree * ntpXiSys = (TTree*) input->Get("ntpXiSys");


	double nevents_mc = ntpMC->GetEntriesFast();
	TString cuts = " McTruthMatch && VtxFit_HowGood==1 && MassFit_HowGood>0";
	TString cut4c = "McTruthMatch && 4CFit_prob>0.01";


	cout << "particle|   #evts (uncut)|    #evts (ratio in %)|   MC ratio in %" << endl;


	//**** PiMinus
	TH1D * h_piminus_tht_uncut = new TH1D("h_piminus_tht_uncut", "h_piminus_tht", 100, 0,10);
	ntpPiMinus->Project("h_piminus_tht_uncut", "piminus_tht", "McTruthMatch && Mother==3122");
	double piminus_uncut =  h_piminus_tht_uncut->GetEntries();

	TH1D * h_piminus_tht = new TH1D("h_piminus_tht", "h_piminus_tht", 100, 0,10);
	ntpPiMinus->Project("h_piminus_tht", "piminus_tht", "McTruthMatch && piminus_HitTag && Mother==3122");
	int piminus =  h_piminus_tht->GetEntries();

	double ratio_piminus_cut = piminus/piminus_uncut;
	double ratio_piminus_mc = piminus/nevents_mc;

	cout << "PiMinus|   " <<  piminus_uncut << "|   " <<  piminus << "(" << ratio_piminus_cut*100 << ")|   " << ratio_piminus_mc*100 << endl;


	//**** PiPlus (AntiLambda0)
	TH1D * h_piplus_tht_uncut = new TH1D("h_piplus_tht_uncut", "h_piplus_tht", 100, 0,10);
	ntpPiPlus->Project("h_piplus_tht_uncut", "piplus_tht", "McTruthMatch && Mother==-3122");
	double piplus_uncut =  h_piplus_tht_uncut->GetEntries();

	TH1D * h_piplus_tht = new TH1D("h_piplus_tht", "h_piplus_tht", 100, 0,10);
	ntpPiPlus->Project("h_piplus_tht", "piplus_tht", "McTruthMatch && piplus_HitTag && Mother==-3122");
	int piplus =  h_piplus_tht->GetEntries();

	double ratio_piplus_cut = piplus/piplus_uncut;
	double ratio_piplus_mc = piplus/nevents_mc;

	cout << "PiPlus(AL0)|   " <<  piplus_uncut << "|   " <<  piplus << "(" << ratio_piplus_cut*100 << ")|   " << ratio_piplus_mc*100 << endl;



	//**** piplus (Xi+)
	TH1D * h_piplus2_tht_uncut = new TH1D("h_piplus2_tht_uncut", "h_piplus2_tht", 100, 0,10);
	ntpPiPlus->Project("h_piplus2_tht_uncut", "piplus_tht", "McTruthMatch && Mother==-3312");
	double piplus_uncut =  h_piplus2_tht_uncut->GetEntries();

	TH1D * h_piplus2_tht = new TH1D("h_piplus2_tht", "h_piplus2_tht", 100, 0,10);
	ntpPiPlus->Project("h_piplus2_tht", "piplus_tht", "McTruthMatch && piplus_HitTag && Mother==-3312");
	int piplus =  h_piplus2_tht->GetEntries();

	double ratio_piplus_cut = piplus/piplus_uncut;
	double ratio_piplus_mc = piplus/nevents_mc;

	cout << "PiPlus(Xi+)|" <<  piplus_uncut << "|" <<  piplus << "(" << ratio_piplus_cut*100 << ")|" << ratio_piplus_mc*100 << endl;


	//**** kaonMinus
	TH1D * h_kaonMinus_tht_uncut = new TH1D("h_kaonMinus_tht_uncut", "h_kaonMinus_tht", 100, 0,10);
	ntpKaonMinus->Project("h_kaonMinus_tht_uncut", "kaonminus_tht", "McTruthMatch");
	double kaonMinus_uncut =  h_kaonMinus_tht_uncut->GetEntries();

	TH1D * h_kaonMinus_tht = new TH1D("h_kaonMinus_tht", "h_kaonMinus_tht", 100, 0,10);
	ntpKaonMinus->Project("h_kaonMinus_tht", "kaonminus_tht", "McTruthMatch && kaonminus_HitTag");
	int kaonMinus =  h_kaonMinus_tht->GetEntries();

	double ratio_kaonMinus_cut = kaonMinus/kaonMinus_uncut;
	double ratio_kaonMinus_mc = kaonMinus/nevents_mc;

	cout << "KaonMinus|   " <<  kaonMinus_uncut << "|   " <<  kaonMinus << "(" << ratio_kaonMinus_cut*100 << ")|   " << ratio_kaonMinus_mc*100 << endl;



	//**** Proton
	TH1D * h_proton_tht_uncut = new TH1D("h_proton_tht_uncut", "h_proton_tht", 100, 0,10);
	ntpProton->Project("h_proton_tht_uncut", "proton_tht", "McTruthMatch && Mother==3122");
	double proton_uncut =  h_proton_tht_uncut->GetEntries();

	TH1D * h_proton_tht = new TH1D("h_proton_tht", "h_proton_tht", 100, 0,10);
	ntpProton->Project("h_proton_tht", "proton_tht", "McTruthMatch && proton_HitTag && Mother==3122");
	int proton =  h_proton_tht->GetEntries();

	double ratio_proton_cut = proton/proton_uncut;
	double ratio_proton_mc = proton/nevents_mc;

	cout << "proton|   " <<  proton_uncut << "|   " <<  proton << "(" << ratio_proton_cut*100 << ")|   " << ratio_proton_mc*100 << endl;




	//**** AntiProton
	TH1D * h_AntiProton_tht_uncut = new TH1D("h_AntiProton_tht_uncut", "h_AntiProton_tht", 100, 0,10);
	ntpAntiProton->Project("h_AntiProton_tht_uncut", "AntiProton_tht", "McTruthMatch");
	double AntiProton_uncut =  h_AntiProton_tht_uncut->GetEntries();

	TH1D * h_AntiProton_tht = new TH1D("h_AntiProton_tht", "h_AntiProton_tht", 100, 0,10);
	ntpAntiProton->Project("h_AntiProton_tht", "AntiProton_tht", "McTruthMatch && AntiProton_HitTag");
	int AntiProton =  h_AntiProton_tht->GetEntries();

	double ratio_AntiProton_cut = AntiProton/AntiProton_uncut;
	double ratio_AntiProton_mc = AntiProton/nevents_mc;

	cout << "AntiProton|   " <<  AntiProton_uncut << "|   " <<  AntiProton << "(" << ratio_AntiProton_cut*100 << ")|   " << ratio_AntiProton_mc*100 << endl;



	//**** lambda0
	TH1D * h_Lambda0_tht_uncut = new TH1D("h_Lambda0_tht_uncut", "h_Lambda0_tht", 100, 0,10);
	ntpLambda0->Project("h_Lambda0_tht_uncut", "Lambda0_tht", "McTruthMatch");
	double Lambda0_uncut =  h_Lambda0_tht_uncut->GetEntries();

	TH1D * h_Lambda0_tht = new TH1D("h_Lambda0_tht", "h_Lambda0_tht", 100, 0,10);
	ntpLambda0->Project("h_Lambda0_tht", "Lambda0_tht", "Lambda0_HitTag && "+cuts);
	int lambda0 =  h_Lambda0_tht->GetEntries();

	double ratio_Lambda0_cut = lambda0/Lambda0_uncut;
	double ratio_Lambda0_mc = lambda0/nevents_mc;

	cout << "lambda0|   " <<  Lambda0_uncut << "|   " <<  lambda0 << "(" << ratio_Lambda0_cut*100 << ")|   " << ratio_Lambda0_mc*100 << endl;



	//**** AntiLambda0
	TH1D * h_antiLambda0_tht_uncut = new TH1D("h_antiLambda0_tht_uncut", "h_antiLambda0_tht", 100, 0,10);
	ntpAntiLambda0->Project("h_antiLambda0_tht_uncut", "antiLambda0_tht", "McTruthMatch");
	double antiLambda0_uncut =  h_antiLambda0_tht_uncut->GetEntries();

	TH1D * h_antiLambda0_tht = new TH1D("h_antiLambda0_tht", "h_antiLambda0_tht", 100, 0,10);
	ntpAntiLambda0->Project("h_antiLambda0_tht", "antiLambda0_tht", "antiLambda0_HitTag && "+cuts);
	int AntiLambda0 =  h_antiLambda0_tht->GetEntries();

	double ratio_antiLambda0_cut = AntiLambda0/antiLambda0_uncut;
	double ratio_antiLambda0_mc = AntiLambda0/nevents_mc;

	cout << "AntiLambda0|   " <<  antiLambda0_uncut << "|   " <<  AntiLambda0 << "(" << ratio_antiLambda0_cut*100 << ")|   " << ratio_antiLambda0_mc*100 << endl;


	//**** XiPlus
	TH1D * h_xiplus_tht_uncut = new TH1D("h_xiplus_tht_uncut", "h_xiplus_tht", 100, 0,10);
	ntpXiPlus->Project("h_xiplus_tht_uncut", "xiplus_tht", "McTruthMatch");
	double xiplus_uncut =  h_xiplus_tht_uncut->GetEntries();

	TH1D * h_xiplus_tht = new TH1D("h_xiplus_tht", "h_xiplus_tht", 100, 0,10);
	ntpXiPlus->Project("h_xiplus_tht", "xiplus_tht", cuts);
	int XiPlus =  h_xiplus_tht->GetEntries();

	double ratio_xiplus_cut = XiPlus/xiplus_uncut;
	double ratio_xiplus_mc = XiPlus/nevents_mc;

	cout << "XiPlus|   " <<  xiplus_uncut << "|   " <<  XiPlus << "(" << ratio_xiplus_cut*100 << ")|   " << ratio_xiplus_mc*100 << endl;


	//**** XiMinus1820
	TH1D * h_XiMinus_tht_uncut = new TH1D("h_XiMinus_tht_uncut", "h_XiMinus_tht", 100, 0,10);
	ntpXiMinus1820->Project("h_XiMinus_tht_uncut", "XiMinus_tht", "McTruthMatch");
	double XiMinus_uncut =  h_XiMinus_tht_uncut->GetEntries();

	TH1D * h_XiMinus_tht = new TH1D("h_XiMinus_tht", "h_XiMinus_tht", 100, 0,10);
	ntpXiMinus1820->Project("h_XiMinus_tht", "XiMinus_tht", cuts);
	int XiMinus1820 =  h_XiMinus_tht->GetEntries();

	double ratio_XiMinus_cut = XiMinus1820/XiMinus_uncut;
	double ratio_XiMinus_mc = XiMinus1820/nevents_mc;

	cout << "XiMinus1820|   " <<  XiMinus_uncut << "|   " <<  XiMinus1820 << "(" << ratio_XiMinus_cut*100 << ")|   " << ratio_XiMinus_mc*100 << endl;


	//**** XiSys
	TH1D * h_XiSys_tht_uncut = new TH1D("h_XiSys_tht_uncut", "h_XiSys_tht", 100, 0,10);
	ntpXiSys->Project("h_XiSys_tht_uncut", "XiSys_tht", "McTruthMatch");
	double XiSys_uncut =  h_XiSys_tht_uncut->GetEntries();

	TH1D * h_XiSys_tht = new TH1D("h_XiSys_tht", "h_XiSys_tht", 100, 0,10);
	ntpXiSys->Project("h_XiSys_tht", "XiSys_tht", cut4c);
	int XiSys =  h_XiSys_tht->GetEntries();

	double ratio_XiSys_cut = XiSys/XiSys_uncut;
	double ratio_XiSys_mc = XiSys/nevents_mc;

	cout << "XiSys|   " <<  XiSys_uncut << "|   " <<  XiSys << "(" << ratio_XiSys_cut*100 << ")|   " << ratio_XiSys_mc*100 << endl;


}
