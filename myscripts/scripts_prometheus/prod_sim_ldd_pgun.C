// Macro for running Panda simulation  with Geant3  or Geant4 (M. Al-Turany)
// This macro is supposed to run the full simulation of the panda detector
// to run the macro:
// root  sim_complete.C  or in root session root>.x  sim_complete_stt.C
// to run with different options:(e.g more events, different momentum, Geant4)
// root  sim_complete.C"(100, "TGeant4",2)"

prod_sim_ldd_pgun(TString outpre="", Int_t nEvents = 100, Int_t PdgType=22 , Float_t mom = 10., double thtmin=0., double thtmax=90., int mult = 1  )
{
	

  //-----User Settings:-----------------------------------------------
  TString  SimEngine      ="TGeant3";
  TString  Workdir        =gSystem->Getenv("VMCWORKDIR");
  //TString  Decfile        =Workdir+"/tutorials/apr13/psi2s_jpsi2pi.dec";
  //TString  Resonance      ="psi(2S)";
  
  TString evtPdlFile = "/hera/panda/jpuetz/myscripts/evt.pdl";

  TString  OutputFile     = outpre+"_sim.root";
  TString  ParOutputfile  = outpre+"_par.root";
  Double_t BeamMomentum   = 15.0; // beam momentum ONLY for the scaling of the dipole field. For the generator use "mom"
  TString  MediaFile      = "media_pnd.geo";
  gDebug                  = 0;
  TString digiFile        = "all.par"; //The emc run the hit producer directly 
  // choose your event generator 
  Bool_t UseEvtGen	      = kFALSE;
  Bool_t UseEvtGenDirect  = kFALSE;     
  Bool_t UseDpm 	      = kFALSE;
  Bool_t UseBoxGenerator  = kTRUE;
  
  
  // DPM or EvtGen?
 // if (Decfile=="DPM") UseDpm=kTRUE;
 // else UseEvtGenDirect=kTRUE;

  //if (!UseBoxGenerator) BeamMomentum = mom;

  //------------------------------------------------------------------
  TLorentzVector fIni(0, 0, mom, sqrt(mom*mom+9.3827203e-01*9.3827203e-01)+9.3827203e-01);  
  TDatabasePDG::Instance()->AddParticle("pbarpSystem","pbarpSystem",fIni.M(),kFALSE,0.1,0, "",88888); 
  //------------------------------------------------------------------
  TStopwatch timer;
  timer.Start();
  gRandom->SetSeed(2704); 

  // Create the Simulation run manager--------------------------------
  FairRunSim *fRun = new FairRunSim();
  fRun->SetName(SimEngine.Data() );
  fRun->SetOutputFile(OutputFile.Data());
  fRun->SetGenerateRunInfo(kFALSE);
  fRun->SetBeamMom(BeamMomentum);
  fRun->SetMaterials(MediaFile.Data());
  fRun->SetUseFairLinks(kTRUE);			 //Added
  FairRuntimeDb *rtdb=fRun->GetRuntimeDb();
  
  // Set the parameters 
  //-------------------------------
  TString allDigiFile = gSystem->Getenv("VMCWORKDIR");
  allDigiFile += "/macro/params/";
  allDigiFile += digiFile;
 
 
  //-------Set the parameter output --------------------
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
  parIo1->open(allDigiFile.Data(),"in");
  rtdb->setFirstInput(parIo1);        

 //---------------------Set Parameter output      ---------- 
  Bool_t kParameterMerged=kTRUE;
  FairParRootFileIo* output=new FairParRootFileIo(kParameterMerged);
  output->open(ParOutputfile.Data());
  rtdb->setOutput(output);

   // Create and add detectors

 //-------------------------  CAVE      -----------------

  FairModule *Cave= new PndCave("CAVE");
  Cave->SetGeometryFileName("pndcave.geo");
  fRun->AddModule(Cave); 
  //-------------------------  Magnet   ----------------- 
  //FairModule *Magnet= new PndMagnet("MAGNET");
  //Magnet->SetGeometryFileName("FullSolenoid_V842.root");
  //Magnet->SetGeometryFileName("FullSuperconductingSolenoid_v831.root");
  //fRun->AddModule(Magnet);
  FairModule *Dipole= new PndMagnet("MAGNET");
  Dipole->SetGeometryFileName("dipole.geo");
  fRun->AddModule(Dipole);
  //-------------------------  Pipe     -----------------
  FairModule *Pipe= new PndPipe("PIPE");
  Pipe->SetGeometryFileName("beampipe_201309.root");
  fRun->AddModule(Pipe);
  //-------------------------  STT       -----------------
  FairDetector *Stt= new PndStt("STT", kTRUE);
  Stt->SetGeometryFileName("straws_skewed_blocks_35cm_pipe.geo");
  fRun->AddModule(Stt);
  //-------------------------  MVD       -----------------
  FairDetector *Mvd = new PndMvdDetector("MVD", kTRUE);
  Mvd->SetGeometryFileName("Mvd-2.1_FullVersion.root");
  fRun->AddModule(Mvd);
  //-------------------------Lambda Disks-----------------
  FairDetector *LamDisks = new PndLddDetector("LamDisks", kTRUE);
  LamDisks->SetGeometryFileName("LambdaDisksSeparatedSupport.root");
  fRun->AddModule(LamDisks);
  //-------------------------  GEM       -----------------
  FairDetector *Gem = new PndGemDetector("GEM", kTRUE);
  Gem->SetGeometryFileName("gem_3Stations_Tube.root");
  fRun->AddModule(Gem);
  //-------------------------  EMC       -----------------
  PndEmc *Emc = new PndEmc("EMC",kTRUE);
  Emc->SetGeometryVersion(1);
  Emc->SetStorageOfData(kFALSE);
  fRun->AddModule(Emc);
  //-------------------------  SCITIL    -----------------
  FairDetector *SciT = new PndSciT("SCIT",kTRUE);
  SciT->SetGeometryFileName("SciTil_201504.root");//SciT->SetGeometryFileName("barrel-SciTil_07022013.root");
  fRun->AddModule(SciT);
  //-------------------------  DRC       -----------------
  PndDrc *Drc = new PndDrc("DIRC", kTRUE);
  Drc->SetGeometryFileName("dirc_l0_p0_updated.root"); 
  Drc->SetRunCherenkov(kFALSE);
  fRun->AddModule(Drc); 
  //-------------------------  DISC      -----------------
  PndDsk* Dsk = new PndDsk("DSK", kTRUE);
  Dsk->SetStoreCerenkovs(kFALSE);
  Dsk->SetStoreTrackPoints(kFALSE);
  fRun->AddModule(Dsk);
  //-------------------------  MDT       -----------------
  PndMdt *Muo = new PndMdt("MDT",kTRUE);
  Muo->SetBarrel("fast");
  Muo->SetEndcap("fast");
  Muo->SetMuonFilter("fast");
  Muo->SetForward("fast");
  Muo->SetMdtMagnet(kTRUE);
  Muo->SetMdtCoil(kTRUE);
  Muo->SetMdtMFIron(kTRUE);
  fRun->AddModule(Muo);
  //-------------------------  FTS       -----------------
  FairDetector *Fts= new PndFts("FTS", kTRUE);
  Fts->SetGeometryFileName("fts.geo");
  fRun->AddModule(Fts); 
  //-------------------------  FTOF      -----------------
  FairDetector *FTof = new PndFtof("FTOF",kTRUE);
  FTof->SetGeometryFileName("ftofwall.root");
  fRun->AddModule(FTof);
  //-------------------------  RICH       ----------------
  FairDetector *Rich= new PndRich("RICH",kFALSE);
  Rich->SetGeometryFileName("rich_v2_shift.geo");
  fRun->AddModule(Rich);

  // Create and Set Event Generator
  //-------------------------------
  FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
  fRun->SetGenerator(primGen);
  
  if(UseBoxGenerator){	// Box Generator
     FairBoxGenerator* boxGen = new FairBoxGenerator(PdgType, mult); // 13 = muon; 1 = multipl.
     boxGen->SetPRange(mom,mom); // GeV/c
     boxGen->SetPhiRange(0., 360.); // Azimuth angle range [degree]
     boxGen->SetThetaRange(thtmin, thtmax); // Polar angle in lab system range [degree]
     boxGen->SetXYZ(0., 0., 0.); // cm
     primGen->AddGenerator(boxGen);
  }
  if(UseDpm){
  	  PndDpmDirect *Dpm= new PndDpmDirect(mom,1);
	  primGen->AddGenerator(Dpm);
  }
  if(UseEvtGenDirect){
      PndEvtGenDirect *EvtGen = new PndEvtGenDirect(Resonance, Decfile.Data(), mom, -1, "", evtPdlFile.Data());
	  EvtGen->SetStoreTree(kTRUE);
	  primGen->AddGenerator(EvtGen);
  }	
	 
 //---------------------Create and Set the Field(s)---------- 
  PndMultiField *fField= new PndMultiField("AUTO");
  fRun->SetField(fField);

 // EMC Hit producer
  //-------------------------------
  PndEmcHitProducer* emcHitProd = new PndEmcHitProducer();
  fRun->AddTask(emcHitProd);

  fRun->SetStoreTraj(kTRUE);

  
 //-------------------------  Initialize the RUN  -----------------  
  fRun->Init();
 //-------------------------  Run the Simulation  -----------------   
  fRun->Run(nEvents);
 //-------------------------  Save the parameters ----------------- 
  rtdb->saveOutput();
 //------------------------Print some info and exit----------------     
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
  
  cout << " Test passed" << endl;
  cout << " All ok " << endl;
  
  //exit(0);

}  
 
