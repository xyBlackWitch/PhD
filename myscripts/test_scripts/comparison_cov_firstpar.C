class PndAnalysis;
class PndAnaPidSelector;
class RhoCandList;
class RhoTuple;

#include "../common_jenny.cpp"


void comparison_cov_firstpar(int nevts=0, TString Path=""){

  TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);

  TStopwatch timer;

  //Output File
  TString outPath = Path;
  TString OutputFile = "cov_output.root";

  //Input simulation Files
  TString DigiFile  = Path + "digi_complete.root";
  TString RecoFile  = Path + "reco_complete.root";
  TString inPIDFile = Path + "pid_complete.root";
  TString inParFile = Path + "simparams.root";
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

  RunAna->AddFriend(DigiFile);
  RunAna->AddFriend(RecoFile);
  RunAna->SetOutputFile(OutputFile);
  RunAna->Init();


  // data reader Object
  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts = theAnalysis->GetEntries();

  RhoCandList piminus, proton;


  TMatrixD CovPion(5,5);
  TMatrixD CovProton(5,5);

	int evt=-1;
	while (theAnalysis->GetEvent() && ++evt<nevts){
		cout << "evt: " << evt << endl;

		TString PidSelection = "PidAlgoIdealCharged";//"PidAlgoMvd;PidAlgoStt;PidAlgoDrc";


		//***Selection with no PID info
		theAnalysis->FillList(piminus, "PionAllMinus", PidSelection);
		theAnalysis->FillList(proton, "ProtonAllPlus", PidSelection);


		//Get piminus

		for (int j=0; j<piminus.GetLength(); ++j){

			bool truth = theAnalysis->McTruthMatch(piminus[j]);

			if(truth){

				TMatrixD piCov(7,7);
				piCov = piminus[j]->Cov7();


				FairTrackParP paramFirst = theAnalysis->GetFirstPar(piminus[j]);

				double cov[15];
				paramFirst.GetCov(cov);

				for (int i=0; i<15; i++){
					if(TMath::Abs(cov[i]<1e-6)) cov[i]=0;
				}

				CovPion[0][0] = cov[0];
				CovPion[0][1] = cov[1];
				CovPion[0][2] = cov[2];
				CovPion[0][3] = cov[3];
				CovPion[0][4] = cov[4];

				CovPion[1][1] = cov[5];
				CovPion[1][2] = cov[6];
				CovPion[1][3] = cov[7];
				CovPion[1][4] = cov[8];

				CovPion[2][2] = cov[9];
				CovPion[2][3] = cov[10];
				CovPion[2][4] = cov[11];

				CovPion[3][3] = cov[12];
				CovPion[3][4] = cov[13];

				CovPion[4][4] = cov[14];

				for (int i=0; i<5; i++){
					for(int j=0; j<5; j++){
						CovPion[j][i]=CovPion[i][j];
					}
				}

				cout << "Covariance matrix for pion: " << endl;

				for(int i=0; i<7; i++){
					for(int j=0; j<7; j++){
						if(TMath::Abs(piCov[i][j])<1e-6) piCov[i][j]=0;
					}
				}
				piCov.Print();

				cout << "First Par Covariance matrix for pion: " << endl;
				CovPion.Print();

			}

			//Get proton

			for (int j=0; j<proton.GetLength(); ++j){

				bool truth = theAnalysis->McTruthMatch(proton[j]);

				if(truth){
					TMatrixD protCov = proton[j]->Cov7();

					FairTrackParP paramFirst = theAnalysis->GetFirstPar(proton[j]);

					double cov[15];
					paramFirst.GetCov(cov);

					for (int i=0; i<15; i++){
						if(TMath::Abs(cov[i]<1e-6)) cov[i]=0;
					}

					CovProton[0][0] = cov[0];
					CovProton[0][1] = cov[1];
					CovProton[0][2] = cov[2];
					CovProton[0][3] = cov[3];
					CovProton[0][4] = cov[4];

					CovProton[1][1] = cov[5];
					CovProton[1][2] = cov[6];
					CovProton[1][3] = cov[7];
					CovProton[1][4] = cov[8];

					CovProton[2][2] = cov[9];
					CovProton[2][3] = cov[10];
					CovProton[2][4] = cov[11];

					CovProton[3][3] = cov[12];
					CovProton[3][4] = cov[13];

					CovProton[4][4] = cov[14];

					for (int i=0; i<5; i++){
						for(int j=0; j<5; j++){
							CovProton[j][i]=CovProton[i][j];
						}
					}

					cout << "Covariance matrix for proton: " << endl;

					for(int i=0; i<7; i++){
						for(int j=0; j<7; j++){
							if(TMath::Abs(protCov[i][j])<1e-6) protCov[i][j]=0;
						}
					}
					protCov.Print();

					cout << "First Par Covariance matrix for proton: " << endl;
					CovProton.Print();
				}
			}

		}
	}
}

