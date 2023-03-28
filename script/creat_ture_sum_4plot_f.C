void creat_ture_sum_4plot_f(TString mcfile, TString outfileroot, TString flux_file){

    double pi = TMath::Pi();
    double solidangle_factor =  50.*50./pi/25.;
    double GC_factor = 4.*pi; 	
    double total_remained;
    const int nbins = 1e2;
    double start = 5e-5;
    double end = 1;
    double step = pow(end/start,1./nbins);
    Double_t xEdges[nbins + 1]={0};
    for(int i=0;i<=nbins;i++){
        xEdges[i] = start*pow(step,i);
    } 

    TGraph *graw=new TGraph("../dat/dm_flux/flux_mass_5e-03.txt");

    TH1D *hraw=new TH1D("hraw","hraw",nbins,xEdges);
    for(int i=1;i<=nbins;i++){
        hraw->SetBinContent(i,graw->Eval(hraw->GetBinCenter(i)));
    }

    TFile *file;

    TChain *chain=new TChain("t3");
    chain->Add(mcfile+"*.root");
    double x[5000], y[5000], z[5000];
    int Hit, nstep;
    double T_before, T_reached;
    double x_before, y_before, z_before;

    chain->SetBranchAddress("x",x);	
    chain->SetBranchAddress("y",y);	
    chain->SetBranchAddress("z",z);	
    chain->SetBranchAddress("Hit",&Hit);	
    chain->SetBranchAddress("nstep",&nstep);	
    chain->SetBranchAddress("T_before",&T_before);	
    chain->SetBranchAddress("T_reached",&T_reached);	
    chain->SetBranchAddress("x_before",&x_before);	
    chain->SetBranchAddress("y_before",&y_before);	
    chain->SetBranchAddress("z_before",&z_before);	

    double Evt_reached = chain->GetEntries("Hit==1&&T_reached/1000.>5e-5&&((nstep==1&&((z[0]-1580)*(z_before-1580)<0))||(nstep>1&&((z[nstep-1]-1580)*(z[nstep-2]-1580)<0)))");
    double Evt_total = chain->GetEntries();
    cout << "total=" << Evt_total << ";reached=" << Evt_reached << endl;
    total_remained += Evt_reached;

    // reached spectrum
    TH1D *hreach=new TH1D("hreach","hreach",nbins,xEdges);

    chain->Draw("T_reached/1000.>>hreach","Hit==1&&((nstep==1&&((z[0]-1580)*(z_before-1580)<0))||(nstep>1&&((z[nstep-1]-1580)*(z[nstep-2]-1580)<0)))","");

    // scaled reached spectrum
    TH1D *hscale=new TH1D("hscale","hscale",nbins,xEdges);
    hscale->SetLineColor(2);
    hscale->SetTitle("Reached flux; Txl [GeV]; d#Phi/dTx [cm^{-2}s^{-1}GeV^{-1}]");

    double scale_factor = Evt_reached/Evt_total*hraw->Integral(0,nbins,"width")*GC_factor*solidangle_factor/hreach->Integral(0,nbins,"width");
    for(int i=1;i<=nbins;i++){
        if(hreach->GetBinContent(i)>0)
            hscale->SetBinContent(i,hreach->GetBinContent(i)*scale_factor);
        else hscale->SetBinContent(i,0);
    }
    hscale->Smooth(3);
    ofstream outfile;
    outfile.open(flux_file);
    for(int i=1;i<=nbins;i++){
        outfile <<hscale->GetBinCenter(i) << " " <<hscale->GetBinContent(i) << endl;
    }
    outfile.close();
    // Save
    file=new TFile(outfileroot,"recreate");
    hscale->Write();

    file->Close();

    chain->Delete();
    hscale->Delete();
    hreach->Delete();


}
