
#include "TFile.h"


void GemHit(RhoCandidate *c){

	int tag = 0;

	PndPidCandidate * pid = (PndPidCandidate*) c->GetRecoCandidate();

	if(pid){

		int gemhits = pid->GetGemHits();
		tag = (gemhits>0) ? 1: 0;
	}

	return tag;
}

void trackBranch(RhoCandidate *c){

	int branch = 0;

	PndPidCandidate * pid = (PndPidCandidate*) c->GetRecoCandidate();

	if(pid){
		branch = pid->GetTrackBranch();
	}

	return branch;
}

void hitTag(RhoCandidate *c){

	int tag = 0;

	PndPidCandidate * pid = (PndPidCandidate*) c->GetRecoCandidate();

	if(pid){

		int gemhits = pid->GetGemHits();
		int mvdhits = pid->GetMvdHits();
		int stthits = pid->GetSttHits();

		if (gemhits>3 || mvdhits>3 ||stthits>3) tag=1;
	}

	return tag;
}

void number_of_particles_leaving_GEM_hits(TString pre=""){

	TString inPIDFile = pre + "_pid_complete.root";
	TString inParFile = pre + "_simparams.root";

	TString OutputFile = pre + "how_many_in_gem.root";

	TString PIDParFile = TString( gSystem->Getenv("VMCWORKDIR")) + "/macro/params/all.par";

	//Initialization
	FairLogger::GetLogger()->SetLogToFile(kFALSE);
	FairRunAna* RunAna = new FairRunAna();
	FairRuntimeDb* rtdb = RunAna->GetRuntimeDb();
	RunAna->SetInputFile(inPIDFile);

	//setup parameter database
	FairParRootFileIo* parIo = new FairParRootFileIo();
	parIo->open(inParFile);
	FairParAsciiFileIo* parIoPID = new FairParAsciiFileIo();
	parIoPID->open(PIDParFile.Data(),"in");

	rtdb->setFirstInput(parIo);
	rtdb->setSecondInput(parIoPID);
	rtdb->setOutput(parIo);

	RunAna->SetOutputFile(OutputFile);
	RunAna->Init();

	//Create output file
	TFile *out = TFile::Open(pre + "_how_many_in_gem_output.root","RECREATE");

	// data reader Object
	PndAnalysis* theAnalysis = new PndAnalysis();
	int nevts = theAnalysis->GetEntries();


	//RhoCandLists for analysis
	RhoCandList piplus, piminus, proton, antiproton, kaonminus;


	TCanvas * c = new TCanvas("c", "c", 0,0,1200,700);
	TH1D * hist = new TH1D("hist", "particles per event going through GEM planes (all particles); particles per event; counts", 7,0,7);

	TH1D * hpim = new TH1D("hpim","particles per event going through GEM planes (#pi^{-}); ; counts", 2,0,2);
	TH1D * hpip1 = new TH1D("hpip1","particles per event going through GEM planes (#pi^{+}(#Lambda)); ; counts", 2,0,2);
	TH1D * hpip2 = new TH1D("hpip2","particles per event going through GEM planes (#pi^{+}(#bar{#Xi})); ; counts", 2,0,2);
	TH1D * hprot = new TH1D("hprot","particles per event going through GEM planes (proton); ; counts", 2,0,2);
	TH1D * haprot = new TH1D("haprot","particles per event going through GEM planes (antiproton); ; counts", 2,0,2);
	TH1D * hkaon = new TH1D("hkaon","particles per event going through GEM planes (K^{-}); ; counts", 2,0,2);


	int evt=-1;
	while (theAnalysis->GetEvent() && ++evt<nevts){

		if ((evt%100)==0) cout << "evt "<< evt <<endl;
//		cout << "Event number: " << evt << endl;

		theAnalysis->FillList(piminus, "PionBestMinus", "PidAlgoIdealCharged");
		theAnalysis->FillList(piplus, "PionBestPlus", "PidAlgoIdealCharged");
		theAnalysis->FillList(proton, "ProtonBestPlus", "PidAlgoIdealCharged");
		theAnalysis->FillList(antiproton, "ProtonBestMinus", "PidAlgoIdealCharged");
		theAnalysis->FillList(kaonminus, "KaonBestMinus", "PidAlgoIdealCharged");

		int sum=0;

		for(int pim=0; pim<piminus.GetLength(); pim++){

			int pimcount = 0;

			int truthmatch = theAnalysis->McTruthMatch(piminus[pim]);
			RhoCandidate * truth = piminus[pim]->GetMcTruth();
			int mother = (truth==0x0)? 88888: truth->TheMother()->PdgCode();
			int branch = trackBranch(piminus[pim]);
			int hittag = hitTag(piminus[pim]);

			int count = (truthmatch&& mother==3122&& branch==44 && hittag==1) ? GemHit(piminus[pim]): 0;

			pimcount += count;

			sum += count;
		}

		TString fillpim = (pimcount==0)? "no GEM hit": "GEM hits";
		hpim->Fill(fillpim, 1);

		for(int pip=0; pip<piplus.GetLength(); pip++){

			int pipcount1 = 0;
			int pipcount2 = 0;

			int truthmatch = theAnalysis->McTruthMatch(piplus[pip]);
			RhoCandidate * truth = piplus[pip]->GetMcTruth();
			int mother = (truth==0x0)? 88888: truth->TheMother()->PdgCode();
			int branch = trackBranch(piplus[pip]);
			int hittag = hitTag(piplus[pip]);


			int count1 = (truthmatch&& mother==-3122&& branch==44 && hittag==1)? GemHit(piplus[pip]): 0;
			int count2 = (truthmatch&& mother==-3312&& branch==44 && hittag==1)? GemHit(piplus[pip]): 0;

			int count = count1 + count2;

			pipcount1 += count1;
			pipcount2 += count2;


			sum += count;
		}

		TString fillpip1 = (pipcount1==0)? "no GEM hit": "GEM hits";
		hpip1->Fill(fillpip1, 1);


		TString fillpip2 = (pipcount2==0)? "no GEM hit": "GEM hits";
		hpip2->Fill(fillpip2, 1);

		for(int prot=0; prot<proton.GetLength(); prot++){

			int protcount=0;

			int truthmatch = theAnalysis->McTruthMatch(proton[prot]);
			RhoCandidate * truth = proton[prot]->GetMcTruth();
			int mother = (truth==0x0)? 88888: truth->TheMother()->PdgCode();
			int branch = trackBranch(proton[prot]);
			int hittag = hitTag(proton[prot]);

			int count = (truthmatch&& mother==3122&& branch==44 && hittag==1)? GemHit(proton[prot]): 0;

			protcount += count;

			sum += count;
		}
		TString fillprot = (protcount==0)? "no GEM hit": "GEM hits";
		hprot->Fill(fillprot, 1);

		for(int aprot=0; aprot<antiproton.GetLength(); aprot++){

			int aprotcount=0;

			int truthmatch = theAnalysis->McTruthMatch(antiproton[aprot]);
			RhoCandidate * truth = antiproton[aprot]->GetMcTruth();
			int mother = (truth==0x0)? 88888: truth->TheMother()->PdgCode();
			int branch = trackBranch(antiproton[aprot]);
			int hittag = hitTag(antiproton[aprot]);

			int count = (truthmatch&& mother==-3122&& branch==44 && hittag==1)? GemHit(antiproton[aprot]): 0;

			aprotcount+=count;


			sum += count;
		}
		TString fillaprot = (aprotcount==0)? "no GEM hit": "GEM hits";
		haprot->Fill(fillaprot, 1);

		for(int km=0; km<kaonminus.GetLength(); km++){

			int kaoncount = 0;

			int truthmatch = theAnalysis->McTruthMatch(kaonminus[km]);
			RhoCandidate * truth = kaonminus[km]->GetMcTruth();
			int mother = (truth==0x0)? 88888: truth->TheMother()->PdgCode();
			int branch = trackBranch(kaonminus[km]);
			int hittag = hitTag(kaonminus[km]);

			int count = (truthmatch&& mother==23314&& branch==44 && hittag==1)? GemHit(kaonminus[km]): 0;

			kaoncount+=count;

			sum += count;
		}
		TString fillkaon = (kaoncount==0)? "no GEM hit": "GEM hits";
		hkaon->Fill(fillkaon,1);


		hist->Fill(sum);


	}

	out->cd();

	hist->Draw();
	hist->Write();

	gStyle->SetOptStat(0);

	TCanvas * cdiv = new TCanvas("cdiv", "cdiv", 0,0,1500,1000);

	cdiv->Divide(2,3);
	cdiv->cd(1);
	hpim->Write("hpim");
	hpim->Draw();
	cdiv->cd(2);
	hpip1->Draw();
	hpip1->Write("hpip1");
	cdiv->cd(3);
	hpip2->Draw();
	hpip2->Write("hpip2");
	cdiv->cd(4);
	hprot->Draw();
	hprot->Write("hprot");
	cdiv->cd(5);
	haprot->Draw();
	haprot->Write("haprot");
	cdiv->cd(6);
	hkaon->Draw();
	hkaon->Write("hkaon");

	out->Save();

}
