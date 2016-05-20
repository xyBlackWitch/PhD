class PndAnalysis;
class PndAnaPidSelector;
class RhoCandList;
class RhoTuple;

#include "../common_jenny.cpp"


void covariance_matrix_test( int nevts=0, TString pre=""){
  
  TDatabasePDG::Instance()-> AddParticle("pbarpSystem","pbarpSystem", 1.9, kFALSE, 0.1, 0,"", 88888);
  
  TStopwatch timer; 

  //Output File
  TString Path ="/home/ikp1/puetz/panda/mysimulations/test/boxgenerator/lambda0/10000_events/";
  TString outPath = Path;
  TString OutputFile = pre + "analysis_output.root";
  
  //Input simulation Files
  TString inPIDFile = Path + "pid_complete.root";
  TString inParFile = Path + "simparams.root";
  TString PIDParFile = TString( gSystem->Getenv("VMCWORKDIR")) + "/macro/params/all.par";
  

//  TFile * recofile = new TFile(Path+"reco_complete.root", "READ");
//  TTree * reco = (TTree*) recofile->Get("cbmsim");
//
//  TBranch * sttmvdgemgentrack = reco->GetBranch("SttMvdGemGenTrack");
//  TBranch * ftstrack = reco->GetBranch("FtsIdealGenTrack");
//
//  TClonesArray * sttmvdgem;
//  sttmvdgemgentrack->SetAddress(&sttmvdgem);
//
//  TClonesArray * fts;
//  ftstrack->SetAddress(&fts);


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


  // data reader Object
  PndAnalysis* theAnalysis = new PndAnalysis();
  if (nevts==0) nevts = theAnalysis->GetEntries();
  
  //RhoCandLists for analysis
  RhoCandList piplus, piminus, lambda0, proton, piplus_diag, piminus_diag, lambda0_diag, proton_diag;
  RhoCandList mclist, all, piplus_higher, piminus_higher, lambda0_higher, proton_higher;

  RhoCandidate * dummyCand = new RhoCandidate();

  //***Mass selector
  double m0_lambda0= TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
  cout<<"Mass of Lambda0: "<<m0_lambda0<<endl;
  RhoMassParticleSelector * lambdaMassSelector = new RhoMassParticleSelector("lambda0", m0_lambda0, 0.3);

  double pbarmom = 3;

  PndRhoTupleQA qa(theAnalysis, pbarmom);
  TH1D * h = new TH1D("h", "prob dist; prob; counts", 100,0,1);
  TH1D * h_diag = new TH1D("h_diag", "prob dist; prob; counts", 100,0,1);
  TH1D * h_higher = new TH1D("h_higher", "prob dist; prob; counts", 100,0,1);

  int evt=-1;
  while (theAnalysis->GetEvent() && ++evt<nevts){

//	reco->GetEntry(evt);

//    cout << "Event number: " << evt << endl;

    //***Setup event shape object

    TString PidSelection = "PidAlgoIdealCharged";//"PidAlgoMvd;PidAlgoStt;PidAlgoDrc";

    theAnalysis->FillList(all, "All", PidSelection);
//    PndEventShape evsh(all, ini, 0.05, 0.1);
    
    //***Selection with no PID info
    theAnalysis->FillList(piminus, "PionAllMinus", PidSelection);
    theAnalysis->FillList(piplus, "PionAllPlus", PidSelection);
    theAnalysis->FillList(proton, "ProtonAllPlus", PidSelection);

    theAnalysis->FillList(piminus_diag, "PionAllMinus", PidSelection);
    theAnalysis->FillList(piplus_diag, "PionAllPlus", PidSelection);
    theAnalysis->FillList(proton_diag, "ProtonAllPlus", PidSelection);

    theAnalysis->FillList(piminus_higher, "PionAllMinus", PidSelection);
    theAnalysis->FillList(piplus_higher, "PionAllPlus", PidSelection);
    theAnalysis->FillList(proton_higher, "ProtonAllPlus", PidSelection);




	//Get piminus

    for (int j=0; j<piminus.GetLength(); ++j){


    	TMatrixD cov_pim = piminus[j]->Cov7();
    	TMatrixD cov_diag(7,7);
    	TMatrixD cov_higher(cov_pim);

//    	FairRecoCandidate* recoCand_pim = piminus[j]->GetRecoCandidate();
//    	int trackid_pim = recoCand_pim->GetTrackIndex();
//    	int trackbranch_pim = recoCand_pim->GetTrackBranch();
//
//    	if(trackbranch_pim==48) PndTrack * recoTrack_pim = (PndTrack*) sttmvdgem->At(trackid_pim);
//    	if(trackbranch_pim==53) PndTrack * recoTrack_pim = (PndTrack*) fts->At(trackid_pim);
//
//    	int flag_pim = recoTrack_pim->GetFlag();
//
////    	if(flag_pim<=0)
//    	cout << "Flag " << flag_pim << endl;


//    	cout << "Cov piminus: " << endl;
//    	cov_pim.Print();

		for(int i=0; i<7; i++){
			cov_diag[i][i]=cov_pim[i][i];
		}
		piminus_diag[j]->SetCov7(cov_diag);

		for(int i=0; i<7; i++){
			for(int k=0; k<7; k++){
				if(i==k) continue;
				else{
					if(cov_higher[i][k]<1e-5) cov_higher[i][k]=0.5;
				}
			}
		}

		piminus_higher[j]->SetCov7(cov_higher);




    }

		//Get Proton
    	for (int j=0; j<proton.GetLength(); j++){
			TMatrixD cov_prot = proton[j]->Cov7();
			TMatrixD cov_protnew(7,7);
			TMatrixD cov_higher(cov_prot);


//	    	FairRecoCandidate* recoCand_prot = proton[j]->GetRecoCandidate();
//	    	int trackid_prot = recoCand_prot->GetTrackIndex();
//	    	int trackbranch_prot = recoCand_prot->GetTrackBranch();
//
//	    	if(trackbranch_prot==48) PndTrack * recoTrack_prot = (PndTrack*) sttmvdgem->At(trackid_prot);
//	    	if(trackbranch_prot==53) PndTrack * recoTrack_prot = (PndTrack*) fts->At(trackid_prot);
//
//	    	int flag_prot = recoTrack_prot->GetFlag();
//
//	    	if(flag_prot<=0) cout << "Flag Proton " << flag_prot << endl;

//			cout << "Cov proton: " << endl;
//        	cov_prot.Print();

			for(int i=0; i<7; i++){
			cov_protnew[i][i]=cov_prot[i][i];
				}

			proton_diag[j]->SetCov7(cov_protnew);

			for(int i=0; i<7; i++){
				for(int k=0; k<7; k++){
					if(i==k) continue;
					else{
						if(cov_higher[i][k]<1e-5) cov_higher[i][k]=0.5;
					}
				}
			}
			proton_higher[j]->SetCov7(cov_higher);



	}


   //***Lambda0 -> PiMinus + Proton
    lambda0.Combine(piminus,proton);
	lambda0.Select(lambdaMassSelector);
    lambda0.SetType(3122);



    for (int j=0; j<lambda0.GetLength(); ++j){

    	PndKinVtxFitter vtxfitter(lambda0[j]);
    	vtxfitter.Fit();
    	RhoCandidate * lambda0Fit = lambda0[j]->GetFit();

    	double prob = vtxfitter.GetProb();
    	RhoCandidate * d1 = lambda0Fit->Daughter(0);
    	RhoCandidate * d2 = lambda0Fit->Daughter(1);

    	TMatrixD covd1 = d1->Cov7();
    	TMatrixD covd2 = d2->Cov7();

//    	cout << "Probability: " << prob << endl;

    	h->Fill(prob);

//    	if (prob>0.9){
//
//
//    	}


    }

    lambda0_diag.Combine(piminus_diag,proton_diag);
	lambda0_diag.Select(lambdaMassSelector);
    lambda0_diag.SetType(3122);

    for (int j=0; j<lambda0_diag.GetLength(); ++j){

    	PndKinVtxFitter vtxfitter(lambda0_diag[j]);
    	vtxfitter.Fit();
    	RhoCandidate * lambda0_diagFit = lambda0_diag[j]->GetFit();

    	double prob = vtxfitter.GetProb();
    	RhoCandidate * d1 = lambda0_diagFit->Daughter(0);
    	RhoCandidate * d2 = lambda0_diagFit->Daughter(1);

    	TMatrixD covd1 = d1->Cov7();
    	TMatrixD covd2 = d2->Cov7();

//    	if (prob>0.9){
//    		if(covd1[0][1]<1e-7 && covd2[0][1]<1e-7){
    			h_diag->Fill(prob);
//				cout << "Cov of " << d1->PdgCode() << endl;
//				covd1.Print();
//				cout << "Cov of " << d2->PdgCode() << endl;
//				covd2.Print();
//
//				cout << "Prob(chi2): " << prob << endl;
//    		}
//    	}


    }


    lambda0_higher.Combine(piminus_higher,proton_higher);
	lambda0_higher.Select(lambdaMassSelector);
    lambda0_higher.SetType(3122);

    for (int j=0; j<lambda0_higher.GetLength(); ++j){

    	PndKinVtxFitter vtxfitter(lambda0_higher[j]);
    	vtxfitter.Fit();
    	RhoCandidate * lambda0_higherFit = lambda0_higher[j]->GetFit();

    	double prob = vtxfitter.GetProb();
    	RhoCandidate * d1 = lambda0_higherFit->Daughter(0);
    	RhoCandidate * d2 = lambda0_higherFit->Daughter(1);

    	TMatrixD covd1 = d1->Cov7();
    	TMatrixD covd2 = d2->Cov7();

//    	if (prob>0.9){
//    		if(covd1[0][1]<1e-7 && covd2[0][1]<1e-7){
    			h_higher->Fill(prob);
//				cout << "Cov of " << d1->PdgCode() << endl;
//				covd1.Print();
//				cout << "Cov of " << d2->PdgCode() << endl;
//				covd2.Print();
//
//				cout << "Prob(chi2): " << prob << endl;
//    		}
//    	}


    }

}//loop over all events

  gStyle->SetOptStat(0);

  TCanvas * c = new TCanvas("c", "c", 0,0,800,500);
  h->SetLineColor(kRed);
  h_diag->SetLineColor(kBlue);
  h_higher->SetLineColor(kBlack);

  TLegend * legend = new TLegend(0.7,0.62,0.86,0.795, "");
  legend->AddEntry(h, "original", "l");
  legend->AddEntry(h_diag, "only diagonal elements", "l");
  legend->AddEntry(h_higher, "higher errors", "l");

  h->Draw();
  h_diag->Draw("Same");
  h_higher->Draw("Same");
  legend->Draw();

//  timer.Stop();
//  Double_t rtime = timer.RealTime();
//  Double_t ctime = timer.CpuTime();
//
//  cout<<"Macro finisched successfully."<<endl;
//  cout<<"Realtime: "<<rtime<<" s, CPU time: "<<ctime<<" s"<<endl;
//  cout<<endl;
//
//
//  exit(0);
 
}
