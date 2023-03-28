void Trans_Tr_form(double sigma,double m_x,double ms, TString file_f, TString file_b, TString file_outtxt){

    double m_N = 131.*931.49;		// xenon mass
    double A = 131.;			// xenon atom number
    double m_p = 938.27;			// proton mass
    double mu_xe = 141.5*141.5/2./m_N;
    double q02 = 3e-4;
    double ND = 1000.*1000./131.*6e23*3600.*24.*365.;		// xenon atom number per ton per year
    double pi = TMath::Pi();
    double sigma_xn = sigma*A*A*pow((m_N*(m_x+m_p))/(m_p*(m_x+m_N)),2);

    const int line = 1e5;
    //const int line = 1e4;
    const int line2 = 2000;

    TGraph *g;
    TH1D *Tx_in;
    TH1D *Tx_raw;
    TH1D *Tx_raw_back;
    TH1D *Tr_file;
    TGraph *gin;

    // DM form factor and Helm form factor
    TString F_sum = Form("((4.*pow(%.1f,2)+2.*%.1f*x)*pow(pow(%.1f,2)+%.2e,2)*(4.*pow(%.1f,2)+2.*%.1f*x)/(16.*pow(%1.f*%.1f*(pow(%.1f,2)+2.*%.1f*x),2)))*pow(3.*exp(-1.*pow(6.92e-3*sqrt(%.3f*x),2)/2.)*(sin(6.92e-3*sqrt(%.3f*x)*sqrt(pow(1.23*pow(%.3f,1/3.)-0.6,2)+7./3*pow(%.6f,2)*0.52*0.52-5.*0.9*0.9))-(6.92e-3*sqrt(%.3f*x)*sqrt(pow(1.23*pow(%.3f,1/3.)-0.6,2)+7./3*pow(%.6f,2)*0.52*0.52-5.*0.9*0.9))*cos(6.92e-3*sqrt(%.3f*x)*sqrt(pow(1.23*pow(%.3f,1/3.)-0.6,2)+7./3*pow(%.6f,2)*0.52*0.52-5.*0.9*0.9)))/pow(6.92e-3*sqrt(%.3f*x)*sqrt(pow(1.23*pow(%.3f,1/3.)-0.6,2)+7./3*pow(%.6f,2)*0.52*0.52-5.*0.9*0.9),3),2)",m_N*1e3,m_N*1e3,ms*1e3,q02*1e3,m_x*1e3,m_N*1e3,m_N*1e3,m_x*1e3,ms*1e3,m_N*1e3,A,A,A,pi,A,A,pi,A,A,pi,A,A,pi);	
    TF1 *fff=new TF1("fff",F_sum,0.1,1e3);

    double Tx_min[line2];
    double Tx[line2];
    double Tr[line2];
    double Tr_flux[line2] = {0};
    double Tr_max[line] = {0};

    /////////////////// Front contribution
    TFile *file=new TFile(file_f,"read");
    Tx_raw= (TH1D*)file->Get("hscale");
    /////////////////// Back contribution
    TFile *file_back=new TFile(file_b,"read");
    Tx_raw_back= (TH1D*)file_back->Get("hscale");

    Tx_in =new TH1D("Tx_in","Tx_in",line,5e-5,1);
    for(int k=1; k<=line; k++){
        Tr_max[k] = Tx_in->GetBinCenter(k)*(Tx_in->GetBinCenter(k)+2*m_x/1000.)/(m_N/1000./2.+Tx_in->GetBinCenter(k));
        Tx_in->SetBinContent(k,(Tx_raw->GetBinContent(k)+Tx_raw_back->GetBinContent(k))/Tr_max[k]);
    }

    Tr_file =new TH1D("Tr_file","Tr_file",line2,0,2e-4);
    for(int i=1; i<line2; i++){
        Tr[i] = Tr_file->GetBinCenter(i);
        Tx_min[i] = sqrt(pow(m_x/1000.-0.5*Tr[i],2)+m_N/1000./2.*Tr[i])-(m_x/1000.-0.5*Tr[i]);
        // find the Tr_min start bin number in the Tx spetrum
        double bin = Tx_in->FindBin(Tx_min[i]);
        Tr_flux[i] = ND*sigma_xn*fff->Eval(Tr[i]*1e6)*Tx_in->Integral(bin,line,"width");
        Tr_file->SetBinContent(i,Tr_flux[i]);
    }
    ofstream outfile;
    outfile.open(file_outtxt);
    for(int i=1; i<line2; i++){
        outfile << Tr_file->GetBinCenter(i) << " " << Tr_file->GetBinContent(i) << endl;
    }
    outfile.close();


    Tx_raw->Delete();
    Tx_raw_back->Delete();
    Tr_file->Delete();
}
