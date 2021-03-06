class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;

// *** routine to only keep PID matched candidates in list
int SelectTruePid(PndAnalysis *ana, RhoCandList &l)
{
	int removed = 0;
	
	for (int ii=l.GetLength()-1;ii>=0;--ii)
	{
		if ( !(ana->McTruthMatch(l[ii])) )
		{
			l.Remove(l[ii]);
			removed++;
		}
	}
	
	return removed;
}

void printCand(RhoCandidate *c)
{
	TLorentzVector lv=c->P4();
	
	cout <<c->PdgCode()<<" ("<<lv.X()<<"/"<<lv.Y()<<"/"<<lv.Z()<<"/"<<lv.E()<<")"<<endl;
}

void countDoubles(RhoCandList &l, int &n1, int &n2, int &n3)
{
	int n_smc  = 0;
	int n_strk = 0;
	int n_both = 0;
	double d = 0.00001;
	
	for (int i=0;i<l.GetLength()-1;++i)
	{
		for (int j=i+1;j<l.GetLength();++j)
		{
			TLorentzVector dl = l[i]->P4() - l[j]->P4();
		
			bool chkmc = (l[i]->GetMcTruth()==l[j]->GetMcTruth());
			bool chktrk = (fabs(dl.X())<d) && (fabs(dl.Y())<d) && (fabs(dl.Z())<d) && (fabs(dl.E())<d);
			if (chkmc) n_smc++;
			if (chktrk) n_strk++;
			if (chktrk && chkmc) n_both++;
		}	
	}
	n1 = n_strk;
	n2 = n_smc;
	n3 = n_both;
}

void ana_complete(int nevts=0)
{
        TDatabasePDG::Instance()->AddParticle("pbarpSystem","pbarpSystem",1.9,kFALSE,0.1,0,"",88888);
        TStopwatch fTimer;
	// *** some variables
	int i=0,j=0, k=0, l=0;
	gStyle->SetOptFit(1011);
	
	// *** the output file for FairRunAna
	TString OutFile="output.root";  
					
	// *** the files coming from the simulation
	TString inPidFile  = "psi2s_jpsi2pi_jpsi_mumu_pid.root";    // this file contains the PndPidCandidates and McTruth
	TString inParFile  = "psi2s_jpsi2pi_jpsi_mumu_par.root";
	
	// *** PID table with selection thresholds; can be modified by the user
	TString pidParFile = TString(gSystem->Getenv("VMCWORKDIR"))+"/macro/params/all_day1.par";	
	
	// *** initialization
	FairLogger::GetLogger()->SetLogToFile(kFALSE);
	FairRunAna* fRun = new FairRunAna();
	FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
	fRun->SetInputFile(inPidFile);
	
	// *** setup parameter database 	
	FairParRootFileIo* parIO = new FairParRootFileIo();
	parIO->open(inParFile);
	FairParAsciiFileIo* parIOPid = new FairParAsciiFileIo();
	parIOPid->open(pidParFile.Data(),"in");
	
	rtdb->setFirstInput(parIO);
	rtdb->setSecondInput(parIOPid);
	rtdb->setOutput(parIO);  
	
	fRun->SetOutputFile(OutFile);
	fRun->Init(); 
	
        // *** create an output file for all histograms
	TFile *out = TFile::Open("output_ana.root","RECREATE");
	
	// *** create some histograms
	TH1F *hmomtrk    = new TH1F("hmomtrk","track momentum (all)",200,0,5);
	TH1F *hthttrk    = new TH1F("hthttrk","track theta (all)",200,0,3.1415);
	
	TH1F *hjpsim_all = new TH1F("hjpsim_all","J/#psi mass (all)",200,0,4.5);
	TH1F *hpsim_all  = new TH1F("hpsim_all","#psi(2S) mass (all)",200,0,5);
	
	TH1F *hjpsim_lpid = new TH1F("hjpsim_lpid","J/#psi mass (loose pid)",200,0,4.5);
	TH1F *hpsim_lpid  = new TH1F("hpsim_lpid","#psi(2S) mass (loose pid)",200,0,5);
	
	TH1F *hjpsim_tpid = new TH1F("hjpsim_tpid","J/#psi mass (tight pid)",200,0,4.5);
	TH1F *hpsim_tpid  = new TH1F("hpsim_tpid","#psi(2S) mass (tight pid)",200,0,5);
	
	TH1F *hjpsim_trpid = new TH1F("hjpsim_trpid","J/#psi mass (true pid)",200,0,4.5);
	TH1F *hpsim_trpid  = new TH1F("hpsim_trpid","#psi(2S) mass (true pid)",200,0,5);
	
	
	TH1F *hjpsim_ftm = new TH1F("hjpsim_ftm","J/#psi mass (full truth match)",200,0,4.5);
	TH1F *hpsim_ftm  = new TH1F("hpsim_ftm","#psi(2S) mass (full truth match)",200,0,5);
	
	TH1F *hjpsim_nm = new TH1F("hjpsim_nm","J/#psi mass (no truth match)",200,0,4.5);
	TH1F *hpsim_nm  = new TH1F("hpsim_nm","#psi(2S) mass (no truth match)",200,0,5);
	
	TH1F *hjpsim_diff = new TH1F("hjpsim_diff","J/#psi mass diff to truth",100,-2,2);
	TH1F *hpsim_diff  = new TH1F("hpsim_diff","#psi(2S) mass diff to truth",100,-2,2);
	
	
	TH1F *hjpsim_vf   = new TH1F("hjpsim_vf","J/#psi mass (vertex fit)",200,0,4.5);
	TH1F *hjpsim_4cf  = new TH1F("hjpsim_4cf","J/#psi mass (4C fit)",200,0,4.5);
	TH1F *hjpsim_mcf  = new TH1F("hjpsim_mcf","J/#psi mass (mass constraint fit)",200,0,4.5);
	
	TH1F *hjpsi_chi2_vf  = new TH1F("hjpsi_chi2_vf", "J/#psi: #chi^{2} vertex fit",100,0,10);
	TH1F *hpsi_chi2_4c   = new TH1F("hpsi_chi2_4c",  "#psi(2S): #chi^{2} 4C fit",100,0,250);
	TH1F *hjpsi_chi2_mf  = new TH1F("hjpsi_chi2_mf", "J/#psi: #chi^{2} mass fit",100,0,10);
	
	TH1F *hjpsi_prob_vf  = new TH1F("hjpsi_prob_vf", "J/#psi: Prob vertex fit",100,0,1);
	TH1F *hpsi_prob_4c   = new TH1F("hpsi_prob_4c",  "#psi(2S): Prob 4C fit",100,0,1);
	TH1F *hjpsi_prob_mf  = new TH1F("hjpsi_prob_mf", "J/#psi: Prob mass fit",100,0,1);
	
	TH2F *hvpos = new TH2F("hvpos","(x,y) projection of fitted decay vertex",100,-2,2,100,-2,2);
	
	//
	// Now the analysis stuff comes...
	//
	
	
	// *** the data reader object
	PndAnalysis* theAnalysis = new PndAnalysis();
	if (nevts==0) nevts= theAnalysis->GetEntries();
	
	// *** RhoCandLists for the analysis
	RhoCandList chrg, muplus, muminus, piplus, piminus, jpsi, psi2s;
	
	// *** Mass selector for the jpsi cands
	double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass();   // Get nominal PDG mass of the J/psi
	RhoMassParticleSelector *jpsiMassSel=new RhoMassParticleSelector("jpsi",m0_jpsi,1.0);
	
	// *** the lorentz vector of the initial psi(2S)
	TLorentzVector ini(0, 0, 6.231552, 7.240065);
	
	// ***
	// the event loop
	// ***
	
	int cntdbltrk=0, cntdblmc=0, cntdblboth=0, cnttrk=0, cnt_dbl_jpsi=0, cnt_dbl_psip=0;
	
	while (theAnalysis->GetEvent() && i++<nevts)
	{
		if ((i%100)==0) cout<<"evt " << i << endl;
				
		// *** Select with no PID info ('All'); type and mass are set 		
		theAnalysis->FillList(chrg,    "Charged");
		theAnalysis->FillList(muplus,  "MuonAllPlus");
		theAnalysis->FillList(muminus, "MuonAllMinus");
		theAnalysis->FillList(piplus,  "PionAllPlus");
		theAnalysis->FillList(piminus, "PionAllMinus");

		// *** momentum and theta histograms
		for (j=0;j<muplus.GetLength();++j) 
		{
			hmomtrk->Fill(muplus[j]->P());
			hthttrk->Fill(muplus[j]->P4().Theta());
		}
		for (j=0;j<muminus.GetLength();++j) 
		{
			hmomtrk->Fill(muminus[j]->P());
			hthttrk->Fill(muminus[j]->P4().Theta());
		}
		
		cnttrk += chrg.GetLength();
		
		int n1, n2, n3;
		
		countDoubles(chrg,n1,n2,n3);
		cntdbltrk  += n1;
		cntdblmc   += n2;
		cntdblboth += n3;		
		
		// *** combinatorics for J/psi -> mu+ mu-
		jpsi.Combine(muplus, muminus);
		
		
		// ***
		// *** do the TRUTH MATCH for jpsi
		// ***
		jpsi.SetType(443);
				
		int nm = 0;
		for (j=0;j<jpsi.GetLength();++j) 
		{
			hjpsim_all->Fill( jpsi[j]->M() );
			
			if (theAnalysis->McTruthMatch(jpsi[j]))
			{ 
				nm++;
				hjpsim_ftm->Fill( jpsi[j]->M() );
			 	hjpsim_diff->Fill( jpsi[j]->GetMcTruth()->M() - jpsi[j]->M() );
			}
			else 
				hjpsim_nm->Fill( jpsi[j]->M() );
		}
		
		if (nm>1) cnt_dbl_jpsi++;
		// ***
		// *** do VERTEX FIT (J/psi)
		// ***
		for (j=0;j<jpsi.GetLength();++j) 
		{
			PndKinVtxFitter vtxfitter(jpsi[j]);	// instantiate a vertex fitter
			vtxfitter.Fit();
			
			double chi2_vtx = vtxfitter.GetChi2();	// access chi2 of fit
			double prob_vtx = vtxfitter.GetProb();	// access probability of fit
			hjpsi_chi2_vf->Fill(chi2_vtx);
			hjpsi_prob_vf->Fill(prob_vtx);			
			
			if ( prob_vtx > 0.01 )				// when good enough, fill some histos
			{
				RhoCandidate *jfit = jpsi[j]->GetFit();	// access the fitted cand
				TVector3 jVtx=jfit->Pos();		// and the decay vertex position
				
				hjpsim_vf->Fill(jfit->M());            
				hvpos->Fill(jVtx.X(),jVtx.Y());
			}
		}
		
		// *** some rough mass selection
		jpsi.Select(jpsiMassSel);
		
		// *** combinatorics for psi(2S) -> J/psi pi+ pi-
		psi2s.Combine(jpsi, piplus, piminus);
		
		
		// ***
		// *** do the TRUTH MATCH for psi(2S)
		// ***
		psi2s.SetType(88888);

		nm = 0;
		for (j=0;j<psi2s.GetLength();++j) 
		{
			hpsim_all->Fill( psi2s[j]->M() );
			
			if (theAnalysis->McTruthMatch(psi2s[j])) 
			{
				nm++;
			 	hpsim_ftm->Fill( psi2s[j]->M() );
			 	hpsim_diff->Fill( psi2s[j]->GetMcTruth()->M() - psi2s[j]->M() );
			}
			else 
				hpsim_nm->Fill( psi2s[j]->M() );
		}			
		if (nm>1) cnt_dbl_psip++;

		
		// ***
		// *** do 4C FIT (initial psi(2S) system)
		// ***
		for (j=0;j<psi2s.GetLength();++j) 
		{
			PndKinFitter fitter(psi2s[j]);	// instantiate the kin fitter in psi(2S)
			fitter.Add4MomConstraint(ini);	// set 4 constraint
			fitter.Fit();		            // do fit
			
			double chi2_4c = fitter.GetChi2();	// get chi2 of fit
			double prob_4c = fitter.GetProb();	// access probability of fit
			hpsi_chi2_4c->Fill(chi2_4c);
			hpsi_prob_4c->Fill(prob_4c);			
			
			if ( prob_4c > 0.01 )			// when good enough, fill some histo
			{
				RhoCandidate *jfit = psi2s[j]->Daughter(0)->GetFit();	// get fitted J/psi
				
				hjpsim_4cf->Fill(jfit->M());
			}
		}		
		
		
		// ***
		// *** do MASS CONSTRAINT FIT (J/psi)
		// ***
		for (j=0;j<jpsi.GetLength();++j) 
		{
			PndKinFitter mfitter(jpsi[j]);		// instantiate the PndKinFitter in psi(2S)
			mfitter.AddMassConstraint(m0_jpsi);	// add the mass constraint
			mfitter.Fit();						// do fit
			
			double chi2_m = mfitter.GetChi2();	// get chi2 of fit
			double prob_m = mfitter.GetProb();	// access probability of fit
			hjpsi_chi2_mf->Fill(chi2_m);
			hjpsi_prob_mf->Fill(prob_m);			
			
			if ( prob_m > 0.01 )				// when good enough, fill some histo
			{
				RhoCandidate *jfit = jpsi[j]->GetFit();	// access the fitted cand
				hjpsim_mcf->Fill(jfit->M());
			}
		}		
		
		
		// ***
		// *** TRUE PID combinatorics
		// ***
		
		// *** do MC truth match for PID type
		SelectTruePid(theAnalysis, muplus);
		SelectTruePid(theAnalysis, muminus);
		SelectTruePid(theAnalysis, piplus);
		SelectTruePid(theAnalysis, piminus);
				
		// *** all combinatorics again with true PID
		jpsi.Combine(muplus, muminus);
		for (j=0;j<jpsi.GetLength();++j) hjpsim_trpid->Fill( jpsi[j]->M() );
		jpsi.Select(jpsiMassSel);
		
		psi2s.Combine(jpsi, piplus, piminus);
		for (j=0;j<psi2s.GetLength();++j) hpsim_trpid->Fill( psi2s[j]->M() );
		
		
		// ***
		// *** LOOSE PID combinatorics
		// ***
		
		// *** and again with PidAlgoMvd;PidAlgoStt;PidAlgoDrc and loose selection
		theAnalysis->FillList(muplus,  "MuonLoosePlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts");
		theAnalysis->FillList(muminus, "MuonLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts");
		theAnalysis->FillList(piplus,  "PionLoosePlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
		theAnalysis->FillList(piminus, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
		
		jpsi.Combine(muplus, muminus);
		for (j=0;j<jpsi.GetLength();++j) hjpsim_lpid->Fill( jpsi[j]->M() );
		jpsi.Select(jpsiMassSel);
		
		psi2s.Combine(jpsi, piplus, piminus);
		for (j=0;j<psi2s.GetLength();++j) hpsim_lpid->Fill( psi2s[j]->M() );
		
		
		// ***
		// *** TIGHT PID combinatorics
		// ***
		
		// *** and again with PidAlgoMvd;PidAlgoStt and tight selection
		theAnalysis->FillList(muplus,  "MuonTightPlus",  "PidAlgoMdtHardCuts");
		theAnalysis->FillList(muminus, "MuonTightMinus", "PidAlgoMdtHardCuts");
		theAnalysis->FillList(piplus,  "PionLoosePlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
		theAnalysis->FillList(piminus, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc");
		
		jpsi.Combine(muplus, muminus);
		for (j=0;j<jpsi.GetLength();++j) hjpsim_tpid->Fill( jpsi[j]->M() );
		jpsi.Select(jpsiMassSel);
		
		psi2s.Combine(jpsi, piplus, piminus);
		for (j=0;j<psi2s.GetLength();++j) hpsim_tpid->Fill( psi2s[j]->M() );
		
	}
		
	// *** write out all the histos
	out->cd();
	
	hmomtrk->Write();
	hthttrk->Write();
	
	hjpsim_all->Write();
	hpsim_all->Write();
	hjpsim_lpid->Write();
	hpsim_lpid->Write();
	hjpsim_tpid->Write();
	hpsim_tpid->Write();
	hjpsim_trpid->Write();
	hpsim_trpid->Write();
	
	hjpsim_ftm->Write();
	hpsim_ftm->Write();
	hjpsim_nm->Write();
	hpsim_nm->Write();
	
	hpsim_diff->Write();
	hjpsim_diff->Write();
	
	hjpsim_vf->Write();
	hjpsim_4cf->Write();
	hjpsim_mcf->Write();
	
	hjpsi_chi2_vf->Write();
	hpsi_chi2_4c->Write();
	hjpsi_chi2_mf->Write();
			
	hjpsi_prob_vf->Write();
	hpsi_prob_4c->Write();
	hjpsi_prob_mf->Write();
			
	hvpos->Write();
		
	out->Save();
        
	// Extract the maximal used memory an add is as Dart measurement
	// This line is filtered by CTest and the value send to CDash
	FairSystemInfo sysInfo;
	Float_t maxMemory=sysInfo.GetMaxMemory();
	cout << "<DartMeasurement name=\"MaxMemory\" type=\"numeric/double\">";
	cout << maxMemory;
	cout << "</DartMeasurement>" << endl;
  
	fTimer.Stop();
	Double_t rtime = fTimer.RealTime();
	Double_t ctime = fTimer.CpuTime();
  
	Float_t cpuUsage=ctime/rtime;
	cout << "<DartMeasurement name=\"CpuLoad\" type=\"numeric/double\">";
	cout << cpuUsage;
	cout << "</DartMeasurement>" << endl;
  
	cout << endl;
	cout << "Real time " << rtime << " s, CPU time " << ctime
		  << "s" << endl;
	cout << "CPU usage " << cpuUsage*100. << "%" << endl;
	cout << "Max Memory " << maxMemory << " MB" << endl;
   
	cout << "Macro finished successfully." << endl;

        exit(0);
	
}
