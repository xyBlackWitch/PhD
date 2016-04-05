class PndAnalysis;
class PndAnaPidSelector;
class RhoCandList;
class RhoTuple;
class PndTrackCand;


void Comparison_trackParameter(){

	//TStopWatch timer;

	//InputFile
	TFile * recoDataKalman = new TFile("/home/ikp1/puetz/panda/fairsoft_mar15/pandaroot/mysimulations/test/tracking/test_2pi_100_events_third_run/reco_complete_with_kalman.root", "READ");
	TFile * recoData = new TFile("/home/ikp1/puetz/panda/fairsoft_mar15/pandaroot/mysimulations/test/tracking/test_2pi_100_events_third_run/reco_complete_without_kalman.root");
	TFile * mcData = new TFile("/home/ikp1/puetz/panda/fairsoft_mar15/pandaroot/mysimulations/test/tracking/test_2pi_100_events_third_run/sim_complete.root");


	//****Get TrackParameter for IdealTrack with Kalman
	TTree * IdealTrackKalman = (TTree*)recoDataKalman->Get("cbmsim");
	TBranch * IdealGenTrack = IdealTrackKalman->GetBranch("IdealGenTrack");

	TClonesArray * arrayKalman;
	IdealGenTrack->SetAddress(&arrayKalman);



	//****GetTrackParamenter for IdealTrack without Kalman
	TTree * IdealTrackWoKalman = (TTree*)recoData->Get("cbmsim");
	TBranch * idealTrack = IdealTrackWoKalman->GetBranch("IdealTrack");


	TClonesArray * ideal;
	idealTrack->SetAddress(&ideal);




	//****Get Track parameter for MC
	TTree * MC = (TTree*) mcData->Get("cbmsim");
	TBranch * mcTrack = MC->GetBranch("MCTrack");


	TClonesArray * mcArray;
	mcTrack->SetAddress(&mcArray);


	int neventsKalman = IdealTrackKalman->GetEntriesFast();
	int neventsWoKalman = IdealTrackWoKalman->GetEntriesFast();
	int neventsMc = MC->GetEntriesFast();

	if (neventsKalman != neventsWoKalman || neventsMc != neventsKalman){
		cout << "Number of Events are different!" << endl;
		break;
	}
	else{
		int nevents = neventsMc;
	}



	cout << "Events for |" << " P |" << " Px |" << " Py |" <<" Pz |" << " Charge | " << "NHits | " << "X | " << "Y | " << "Z | " << endl;
	int selectedEvents[]= {15,40,67,93};
	int length = sizeof(selectedEvents)/sizeof(selectedEvents[0]);


	for (int evt=0; evt < length; evt ++){
		int eventnumber = selectedEvents[evt];
		IdealTrackKalman -> GetEntry(eventnumber);

		for (int entry=0; entry < 2; entry++){//arrayKalman->GetEntries(); entry++){

			PndTrack * trackKalman = (PndTrack*) arrayKalman->At(entry);
			if ( 0x0 == trackKalman ) continue;


			PndTrackCand * track_cand_kalman = (PndTrackCand*) trackKalman->GetTrackCandPtr();
			int nhits_kalman = track_cand_kalman->GetNHits();

			FairTrackParP fairTrackFirstParam_Kalman = trackKalman->GetParamFirst();
			FairTrackParP fairTrackLastParam_Kalman = trackKalman->GetParamLast();

			TVector3 * First_momentum_kalman = (TVector3*) fairTrackFirstParam_Kalman.GetMomentum();
			TVector3 * First_position_kalman = (TVector3*) fairTrackFirstParam_Kalman.GetPosition();
			TVector3 * Last_momentum_kalman = (TVector3*) fairTrackLastParam_Kalman.GetMomentum();
			TVector3 * Last_position_kalman = (TVector3*) fairTrackLastParam_Kalman.GetPosition();

			double P_first = sqrt( First_momentum_kalman->Px()**2 + First_momentum_kalman->Py()**2 + First_momentum_kalman->Pz()**2 );
			cout << "Event " << eventnumber << " Track " << entry << " with Kalman FirstParams: " << P_first << " |" << First_momentum_kalman->Px() << " |" << First_momentum_kalman->Py()
					<< " |" << First_momentum_kalman->Pz() << " |" << fairTrackFirstParam_Kalman.GetQ() << " |" << nhits_kalman << " |" << First_position_kalman->X() << " |"
					<< First_position_kalman->Y() << " |" << First_position_kalman->Z() <<  endl;


			double P_Last = sqrt( Last_momentum_kalman->Px()**2 + Last_momentum_kalman->Py()**2 + Last_momentum_kalman->Pz()**2 );
			cout << "Event " << eventnumber << " Track " << entry << " with Kalman LastParams: " << P_Last << " |" << Last_momentum_kalman->Px() << " |" << Last_momentum_kalman->Py()
					<< " |" << Last_momentum_kalman->Pz()<< " |" << fairTrackLastParam_Kalman.GetQ()<< " |" << nhits_kalman << " |" << Last_position_kalman->X() << " |"
					<< Last_position_kalman->Y() << " |" << Last_position_kalman->Z() <<  endl;
			cout << "------------------" << endl;
		}
		cout<< "------------------------------------------------------------------------------------------------" << endl;


		idealTrack -> GetEntry(eventnumber);

		for (int entry_wo_kalman = 0; entry_wo_kalman < 2 ; entry_wo_kalman++){//ideal->GetEntries(); entry_wo_kalman++){

			PndTrack * track = (PndTrack*) ideal->At(entry_wo_kalman);
			if (0x0 == track) continue;


			PndTrackCand * track_cand = (PndTrackCand*) track->GetTrackCandPtr();
			int nhits = track_cand->GetNHits();

			FairTrackParP fairTrackFirstParam = track->GetParamFirst();
			FairTrackParP fairTrackLastParam = track->GetParamLast();


			TVector3 * First_momentum = (TVector3*) fairTrackFirstParam.GetMomentum();
			TVector3 * First_position = (TVector3*) fairTrackFirstParam.GetPosition();
			TVector3 * Last_momentum = (TVector3*) fairTrackLastParam.GetMomentum();
			TVector3 * Last_position = (TVector3*) fairTrackLastParam.GetPosition();

			double P_first = sqrt( First_momentum->Px()**2 + First_momentum->Py()**2 + First_momentum->Pz()**2 );
			cout << "Event " << eventnumber << " Track " << entry_wo_kalman << " without Kalman (FirstParams): " << P_first << " |" << First_momentum->Px() << " |" << First_momentum->Py()
					<< " |" << First_momentum->Pz()<< " |" << fairTrackFirstParam.GetQ() << " |" <<  nhits << " |" << First_position->X() << " |"
					<< First_position->Y() << " |" << First_position->Z() <<  endl;


			double P_Last = sqrt( Last_momentum->Px()**2 + Last_momentum->Py()**2 + Last_momentum->Pz()**2 );
			cout << "Event " << eventnumber << " Track " << entry_wo_kalman << " without Kalman (LastParams): " << P_Last << " |" << Last_momentum->Px() << " |" << Last_momentum->Py()
					<< " |" << Last_momentum->Pz()<< " |" << fairTrackLastParam.GetQ() << " |" <<  nhits << " |" << Last_position->X() << " |"
					<< Last_position->Y() << " |" << Last_position->Z() <<  endl;

			cout << "------------------" << endl;
		}

		cout << "-----------------------------------------------------------------------------------------------------" << endl;

		mcTrack->GetEntry(eventnumber);

		for (int entry_mc = 0; entry_mc < 10; entry_mc++){

			PndMCTrack * track_mc = (PndMCTrack*) mcArray->At(entry_mc);
			if (0x0 == track_mc) continue;
			int pdg = track_mc->GetPdgCode();
			int mother = track_mc->GetMotherID();
			if ((pdg == 211 || pdg == -211) && mother==-1){
				int charge = (pdg==211) ? 1 : -1;
				TVector3 * momentum_mc = (TVector3*) track_mc->GetMomentum();
				TVector3 * position_mc = (TVector3*) track_mc->GetStartVertex();
				double P = sqrt( momentum_mc->Px()**2 + momentum_mc->Py()**2 + momentum_mc->Pz()**2 );
				cout << "Event " << eventnumber << " Track " << entry_mc << " MC: " << P << " |" << momentum_mc->Px() << " |" << momentum_mc->Py()
						<< " |" << momentum_mc->Pz()<< " |" << charge << " |" <<  "N/A" << " |" << position_mc->X() << " |"
						<< position_mc->Y() << " |" << position_mc->Z() <<  endl;

				cout << "------------------" << endl;
			}
			else{
				continue;
			}
		}

	}
}
