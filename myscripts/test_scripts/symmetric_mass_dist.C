


void symmetric_mass_dist(TString prefix="", TString ntuple="ntpXiMinus1820", double min=1.75, double max=1.9, TString key="VtxFit"){

	TFile * data = new TFile(prefix, "READ");
	TTree * ntp = (TTree*) data->Get(ntuple);


	TString cut = "VtxFit_HowGood==1";

	int bins=500;

	TH1D * h_m = new TH1D("h_m", "Mass distribution;m/GeV/c^{2}; counts", bins, min, max);
	ntp->Project("h_m", key+"_m", cut);

	TCanvas * c = new TCanvas("c", "c", 0,0,800,500);
	h_m ->Draw();

	double entries = h_m->GetEntries();



	double mass[500];
	double delta[500];


	for( int bin=0; bin<bins; bin++){
		double e = h_m->GetBinWidth(bin);

		if(bin==0){
			int j=1;
		}
		else{
			int j=0;
		}
		double diff = 0;
		double center = h_m->GetBinCenter(bin);
//		double t1 = h_m->GetBinContent(bin);

		while (j<bins){
			double t2 = h_m->GetBinCenter(j);

			diff +=(t2-center)**2/e**2;
			j++;

		}
		mass[bin]=center;
		delta[bin]=diff/(entries);

	}
	TCanvas *c2 = new TCanvas("c2", "c2", 0,0,800,500);
	TGraph * graph = new TGraph(bins, mass, delta);
//	graph->GetYaxis()->SetRangeUser(0,6);

	graph->Draw();
}
