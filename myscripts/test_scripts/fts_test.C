
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

void fts_test(){

	TString inPath = "~/panda/mysimulations/test/XiMinus1820_";

	TFile * simFile = new TFile(inPath+"sim_complete.root", "READ");
	TFile * recoFile = new TFile(inPath+"reco_complete.root", "READ");
	TFile * pidFile =  new TFile(inPath+"pid_complete.root", "READ");
	TFile * anaFile = new TFile(inPath+"output_ana.root", "READ");

	TTree * sim = (TTree*) simFile->Get("cbmsim");
	TTree * reco = (TTree*) recoFile->Get("cbmsim");
	TTree * pid = (TTree*) pidFile->Get("cbmsim");
	TTree * proton = (TTree*) anaFile->Get("ntpProton");
	TTree * AntiProton = (TTree*) anaFile->Get("ntpAntiProton");

	TBranch * mc = sim->GetBranch("MCTrack");
	TBranch * pid_branch = pid->GetBranch("PidChargedCand");
	TBranch * reco_fts = reco->GetBranch("FtsIdealGenTrack");



	TClonesArray * mcArray;
	TClonesArray * pidArray;
	TClonesArray * fts;

	mc->SetAddress(&mcArray);
	pid_branch->SetAddress(&pidArray);
	reco_fts->SetAddress(&fts);

	mc->GetEntry(45);
	reco->GetEntry(45);
	pid->GetEntry(45);
	proton->GetEntry(47);
	AntiProton->GetEntry(37);

	int mcIndex_proton=0;
	int mcIndex_antiproton=0;

	cout << "************************** MC **************************" << endl;
	for (int mc_entry=0; mc_entry<10; mc_entry++){
		PndMCTrack * MCTrack = (PndMCTrack*) mcArray->At(mc_entry);
		if(0x0==MCTrack) continue;
		int pdgCode = MCTrack->GetPdgCode();

		if (pdgCode==2212 || pdgCode==-2212){
			if (pdgCode==2212) mcIndex_proton = mc_entry;
			else mcIndex_antiproton=mc_entry;

			TLorentzVector * pMC = (TLorentzVector*) MCTrack->Get4Momentum();
			TVector3 * pos = (TVector3*) MCTrack->GetStartVertex();
			cout << "Particle PdgCode: " << pdgCode << ", Momentum: (" << pMC->Px() << ", " << pMC->Py() << ", " << pMC->Pz()  << "), Energy: " << pMC->E() << endl;
			cout << "Particle PdgCode: " << pdgCode << ", Position: (" << pos->X() << ", " << pos->Y() << ", " << pos->Z()  << ")" << endl;
			cout << "**************************" << endl;
		}
	}

	int npid = pidArray->GetEntries();
	int nFts = fts->GetEntries();
	int pdg = 88888;


	cout << "************************** Reco(Fts) **************************" << endl;
	for (int evtFts = 0 ; evtFts<nFts; evtFts++){
		PndTrack * ftsTrack = (PndTrack*) fts->At(evtFts);
		if(0x0==ftsTrack) continue;

		FairTrackParP firstparam = ftsTrack->GetParamFirst();
		FairTrackParP lastparam = ftsTrack->GetParamLast();

		TVector3 * pReco_first = (TVector3*) firstparam.GetMomentum();
		TVector3 * pReco_last = (TVector3*) lastparam.GetMomentum();

		TVector3 * pReco_Hitfirst = (TVector3*) firstparam.GetPosition();
		TVector3 * pReco_Hitlast = (TVector3*) lastparam.GetPosition();

		PndTrackCand * track_cand_fts = (PndTrackCand*) ftsTrack->GetTrackCandPtr();
		int trackid = track_cand_fts->getMcTrackId();
		if (trackid==mcIndex_proton){
			pdg = 2212;
		}
		else if (trackid == mcIndex_antiproton){
			pdg = -2212;
		}

		else{
			continue;
//			PndMCTrack * MCTrack = (PndMCTrack*) mcArray->At(track_cand_fts->getMcTrackId());
//			if(0x0==MCTrack) continue;
//			int pdg = MCTrack->GetPdgCode();
		}

		cout << "First Hit: Particle PdgCode: " << pdg << ", Momentum: (" << pReco_first->Px() << ", " << pReco_first->Py() << ", " << pReco_first->Pz()  << ")" << endl;
		cout << "First Hit: Particle PdgCode: " << pdg << ", Position: (" << pReco_Hitfirst->x() << ", " << pReco_Hitfirst->y() << ", " << pReco_Hitfirst->z()  << ")" << endl;
		cout << "*******" << endl;
		cout << "Last Hit: Particle PdgCode: " << pdg << ", Momentum: (" << pReco_last->Px() << ", " << pReco_last->Py() << ", " << pReco_last->Pz()  << ")" << endl;
		cout << "Last Hit: Particle PdgCode: " << pdg << ", Position: (" << pReco_Hitlast->x() << ", " << pReco_Hitlast->y() << ", " << pReco_Hitlast->z()  << ")" << endl;
		cout << "**************************" << endl;
	}

	cout << "************************** PidChargedCand **************************" << endl;

	for(int evtPid=0; evtPid<npid; evtPid++){
		PndPidCandidate * pidcand = (PndPidCandidate*) pidArray->At(evtPid);
		int mcindex = pidcand->GetMcIndex();

		if(mcindex==mcIndex_proton || mcindex==mcIndex_antiproton){
			int pdg = (mcindex==mcIndex_proton)? 2212: -2212;
			int charge = pidcand->GetCharge();
			TVector3 * mom = (TVector3*) pidcand->GetMomentum();
			TVector3 * firstHit = (TVector3*) pidcand->GetFirstHit();
			TVector3 * lastHit = (TVector3*) pidcand->GetLastHit();

			cout << "First Hit: Particle PdgCode: " << pdg << ", Momentum: (" << mom->Px() << ", " << mom->Py() << ", " << mom->Pz()  << ")" << endl;
			cout << "First Hit: Particle PdgCode: " << pdg << ", Position: (" << firstHit->x() << ", " << firstHit->y() << ", " << firstHit->z()  << ")" << endl;
			cout << "*******" << endl;
			cout << "Last Hit: Particle PdgCode: " << pdg << ", Position: (" << lastHit->x() << ", " << lastHit->y() << ", " << lastHit->z()  << ")" << endl;
			cout << "**************************" << endl;


		}
	}

	cout << "************************** ana **************************" << endl;

	TLeaf * pdg_proton = proton->GetLeaf("proton_pdg");

	TLeaf * Px_proton = proton->GetLeaf("proton_px");
	TLeaf * Py_proton = proton->GetLeaf("proton_py");
	TLeaf * Pz_proton = proton->GetLeaf("proton_pz");


	cout << "Particle PdgCode: " << pdg_proton->GetValue() << ", Momentum: (" << Px_proton->GetValue() << ", " << Py_proton->GetValue() << ", " << Pz_proton->GetValue()  << ")" << endl;

	TLeaf * x_proton = proton->GetLeaf("proton_x");
	TLeaf * y_proton = proton->GetLeaf("proton_y");
	TLeaf * z_proton = proton->GetLeaf("proton_z");

	cout << "Particle pdgCode: " << pdg_proton->GetValue() << ", Position: (" << x_proton->GetValue() << ", " << y_proton->GetValue() << ", " << z_proton->GetValue()  << ")" << endl;
	cout << "**************************" << endl;



	TLeaf * pdg_AntiProton = AntiProton->GetLeaf("AntiProton_pdg");

	TLeaf * Px_AntiProton = AntiProton->GetLeaf("AntiProton_px");
	TLeaf * Py_AntiProton = AntiProton->GetLeaf("AntiProton_py");
	TLeaf * Pz_AntiProton = AntiProton->GetLeaf("AntiProton_pz");


	cout << "Particle PdgCode: " << pdg_AntiProton->GetValue() << ", Momentum: (" << Px_AntiProton->GetValue() << ", " << Py_AntiProton->GetValue() << ", " << Pz_AntiProton->GetValue()  << ")" << endl;

	TLeaf * x_AntiProton = AntiProton->GetLeaf("AntiProton_x");
	TLeaf * y_AntiProton = AntiProton->GetLeaf("AntiProton_y");
	TLeaf * z_AntiProton = AntiProton->GetLeaf("AntiProton_z");

	cout << "Particle pdgCode: " << pdg_AntiProton->GetValue() << ", Position: (" << x_AntiProton->GetValue() << ", " << y_AntiProton->GetValue() << ", " << z_AntiProton->GetValue()  << ")" << endl;
	cout << "**************************" << endl;








}
