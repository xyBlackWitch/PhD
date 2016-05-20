

void select_crazy_tracks(){

	TFile * reco = new TFile("/home/ikp1/puetz/panda/fairsoft_mar15/pandaroot/mysimulations/test/tracking/test_2pi_100_events_third_run/reco_complete_with_kalman.root", "READ");
	TFile * mc = new TFile("/home/ikp1/puetz/panda/fairsoft_mar15/pandaroot/mysimulations/test/tracking/test_2pi_100_events_third_run/sim_complete.root", "READ");

	TTree * recoTree = (TTree*) reco->Get("cbmsim");
	TTree * mcTree = (TTree*) mc->Get("cbmsim");

	TBranch * IdealGenTrack = recoTree -> GetBranch("IdealGenTrack");
	TBranch * McTrack = mcTree ->GetBranch("MCTrack");


	TClonesArray * recoArray;
	IdealGenTrack -> SetAddress(&recoArray);

	TClonesArray * mcArray;
	McTrack -> SetAddress(&mcArray);

	int nevents_reco = recoTree ->GetEntriesFast();
	int nevents_mc = mcTree ->GetEntriesFast();

	if (nevents_reco != nevents_mc) cout << "Number of Events are not equal!" << endl;
	int counter = 0;

	for (int reco_evt = 0; reco_evt < nevents_reco; reco_evt++){
		double PxLast[2], PxFirst[2], PxMc[2], PFirst[2], PMc[2], QMc[2], QFirst[2];
		recoTree -> GetEntry(reco_evt);
		mcTree -> GetEntry(reco_evt);

		for (int reco_entry=0; reco_entry < 2; reco_entry++){

			PndTrack * recoTrack = (PndTrack*) recoArray->At(reco_entry);
			if (0x0 == recoTrack) continue;

			FairTrackParP FirstHitParamReco = recoTrack->GetParamFirst();
			FairTrackParP LastHitParamReco = recoTrack->GetParamLast();
			PxFirst[reco_entry] = FirstHitParamReco.GetPx();
			PFirst[reco_entry] = sqrt( FirstHitParamReco.GetPx()**2 + FirstHitParamReco.GetPy()**2+ FirstHitParamReco.GetPz()**2 );
			QFirst[reco_entry] = FirstHitParamReco.GetQ();
			//PxLast[reco_entry] = LastHitParamReco.GetPx();



		}

		for (int mc_entry=0; mc_entry < mcArray->GetEntries(); mc_entry++){

			PndMCTrack * mcTrack = (PndMCTrack*) mcArray -> At(mc_entry);
			if (0x0 == mcTrack) continue;

			int pdg = mcTrack->GetPdgCode();
			int mother = mcTrack->GetMotherID();

			if ((pdg==211 || pdg==-211) && mother==-1){
				PxMc[mc_entry-1] = mcTrack->GetMomentum().Px();
				PMc[mc_entry-1] = sqrt( mcTrack->GetMomentum().Px()**2 + mcTrack->GetMomentum().Py()**2+ mcTrack->GetMomentum().Pz()**2 );
				QMc[mc_entry-1] = (pdg == 211)? 1 : -1;


			}
		}

		int length_array = sizeof(PxMc)/sizeof(PxMc[0]);

		for (int array_entry=0; array_entry < length_array; array_entry++){


			float Px_diff = ((int)(PxMc[array_entry]*100)/100.) - ((int)(PxFirst[array_entry]*100)/100.);
			Px_diff = (int)(Px_diff*10)/10. ;

			float P_diff = ((int) (PMc[array_entry]*100)/100.) - ((int)(PFirst[array_entry]*100)/100.);
			P_diff = (int)(P_diff*10)/10. ;



			if (Px_diff == 0 ) continue;
			else{
				if (P_diff != 0 || fabs(Px_diff)<=0.11) continue;
				cout << "Event: " << reco_evt << " with charge : " << QFirst[array_entry] << endl;
				counter++;
			}

		}


	}
	cout << "Number of crazy track parameter: " << counter << endl;


}
