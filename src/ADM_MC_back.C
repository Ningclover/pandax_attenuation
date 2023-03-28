#include "../include/General_function.h"

int main(int argc, char** argv)
{
    if(argc!=6) {
        cout<<"syntax :"<<argv[0]<<" <xsec> <m_dm> <m_mediator> <input_energy> <outfile> "<<endl;
        exit(1);
    }
    double xsec = atof(argv[1]);
    double xsec1, xsec2;
    double m_dm = atof(argv[2]);
    double m_med = atof(argv[3]);
    double pi = TMath::Pi();
    double initial_sec = cal_sigma(m_dm/1000,m_med/1000);
    string mass = argv[2];

    gRandom->SetSeed(0);

    // The energy PDF
    TF1 *func_spec = get_eng_spec(TString(argv[4]));

    double mp = 938.;
    const double density = 2.8;
    // Basic Earth information
    double radius = 6371*1e3;
    double JP_deep = 2400.;
    double Core_radius = 3480*1e3;	

    const Int_t kMaxStep = 4096;
    Int_t nstep;
    Int_t stepno[kMaxStep];
    Double_t x[kMaxStep];
    Double_t y[kMaxStep];
    Double_t z[kMaxStep];
    Double_t theta[kMaxStep];
    Double_t T[kMaxStep];
    int Hit;
    Double_t T_reached;
    Double_t T_before;
    Double_t x_before;
    Double_t y_before;
    Double_t z_before;
    Double_t AA;
    Double_t aa;
    double X0, Y0, Z0;
    double step_m[kMaxStep];
    double theta_each[kMaxStep];
    int type[kMaxStep];

    TFile *outfile=new TFile(TString(argv[5]),"recreate");
    TTree *t3 = new TTree("t3","Reconst ntuple");
    t3->Branch("nstep",&nstep,"nstep/I");
    t3->Branch("T_reached",&T_reached,"T_reached/D");
    t3->Branch("T_before",&T_before,"T_before/D");
    t3->Branch("Hit",&Hit,"Hit/I");
    t3->Branch("x",x,"px[nstep]/D");
    t3->Branch("y",y,"py[nstep]/D");
    t3->Branch("z",z,"pz[nstep]/D");
    t3->Branch("theta",theta,"theta[nstep]/D");
    t3->Branch("T",T,"T[nstep]/D");
    t3->Branch("X0",&X0,"X0/D");
    t3->Branch("Y0",&Y0,"Y0/D");
    t3->Branch("Z0",&Z0,"Z0/D");
    t3->Branch("AA",&AA,"AA/D");
    t3->Branch("aa",&aa,"aa/D");
    t3->Branch("stepno",stepno,"stepno[nstep]/I");
    t3->Branch("x_before",&x_before,"x_before/D");
    t3->Branch("y_before",&y_before,"y_before/D");
    t3->Branch("z_before",&z_before,"z_before/D");
    t3->Branch("theta_each",theta_each,"theta_each[nstep]/D");
    t3->Branch("step_m",step_m,"step_m[nstep]/D");
    t3->Branch("type",type,"type[nstep]/I");

    // Get the cross section ratio file

    TString ms_string;
    if(m_med>99){
        ms_string = Form("%.1f",m_med/1000);
    }else if(m_med>9){
        ms_string = Form("%.2f",m_med/1000);
    }else{
        ms_string = Form("%.3f",m_med/1000);
    }

    TGraph *g_cs_ratio_Oxygen = new TGraph(Form("../dat/xsec_data/sigma_tot/Oxygen/sigma_ratio/logmchi-%.3f_mS-",log10(m_dm/1000))+ms_string+".txt");
    TGraph *g_cs_ratio_Iron = new TGraph(Form("../dat/xsec_data/sigma_tot/Iron/sigma_ratio/logmchi-%.3f_mS-",log10(m_dm/1000))+ms_string+".txt");
    // Get the elastic effective cross section
    TGraph *g_cs_elas_Oxygen = new TGraph(Form("../dat/xsec_data/sigma_tot/Oxygen/sigma_ela/logmchi-%.3f_mS-",log10(m_dm/1000))+ms_string+".txt");
    TGraph *g_cs_elas_Iron = new TGraph(Form("../dat/xsec_data/sigma_tot/Iron/sigma_ela/logmchi-%.3f_mS-",log10(m_dm/1000))+ms_string+".txt");
    // Get the QE effective cross section
    TGraph *g_cs_QE_Oxygen = new TGraph(Form("../dat/xsec_data/sigma_tot/Oxygen/sigma_QE/logmchi-%.3f_mS-",log10(m_dm/1000))+ms_string+".txt");
    TGraph *g_cs_QE_Iron = new TGraph(Form("../dat/xsec_data/sigma_tot/Iron/sigma_QE/logmchi-%.3f_mS-",log10(m_dm/1000))+ms_string+".txt");



    // Start the simulation
    int evtno(0);
    int N_max = 1e4;
    while(evtno<N_max){

        int cross = 0;	

        if(int(evtno/N_max*100)==evtno/N_max*100) cout << evtno/N_max*100. << "%;\n";
        // input energy
        double T_dm_in = 1000.*func_spec->GetRandom();



        Double_t T_dm_before = T_dm_in;
        T_before = T_dm_before;
        Double_t T_dm_after, theta_cm, theta_lab_relative;

        ////// input position
        double input_theta = gRandom->Uniform(0,pi);
        double input_phi = gRandom->Uniform(0,2.*pi);

        x_before = radius*sin(input_theta)*cos(input_phi);
        y_before = radius*sin(input_theta)*sin(input_phi);
        z_before = radius*cos(input_theta);

        TVector3 v_pos_before(x_before,y_before,z_before), v_pos_after(0,0,0); //position

        ///////////////////
        // input angle
        AA = gRandom->Uniform(-pi,pi);
        aa = asin(gRandom->Uniform(-1,1));

        // define the input angle
        Z0 = 1.*cos(pi/2.-aa);
        if(AA>0){
            X0 = 1.*sin(pi/2.-aa)*cos(AA);
            Y0 = 1.*sin(pi/2.-aa)*sin(AA);}
        else{
            X0 = 1.*sin(pi/2.-aa)*cos(2.*pi+AA);
            Y0 = 1.*sin(pi/2.-aa)*sin(2.*pi+AA);}

        TVector3 vec_before(X0,Y0,Z0), vec_after(0,0,1);//direction

        // start the steps
        if(T_dm_before<=0) cout<<evtno<<endl;  

        int counter(0);
        while(T_dm_before>0.05 && counter<kMaxStep){    

            // choose the scattering nucleus
            // And obtained the total cross section
            int A = gRandom->Uniform(0,1);
            double lambda;
            string element_one;
            // as no Mg24, so A^2 is used for the approximation that O16 is ~91% and Fe56 is !9%
            if(A>0.91){
                A = 56;
                lambda = 180.;
                element_one = "Oxygen";
                xsec1 = g_cs_QE_Iron->Eval(T_dm_before/1000.);
                xsec2 = g_cs_elas_Iron->Eval(T_dm_before/1000.);
            }
            else{
                A = 16;
                lambda = 250.6;
                element_one = "Iron";
                xsec1 = g_cs_QE_Oxygen->Eval(T_dm_before/1000.);
                xsec2 = g_cs_elas_Oxygen->Eval(T_dm_before/1000.);
            }
            double m_tgt = A*mp;
            double sigma = (xsec1+xsec2)*xsec/initial_sec;	

            // choose the scattering type, QE or elastic
            double noim = gRandom->Uniform(0,1);
            double near = 0;
            double elas_total_ratio;
            if(A==16) elas_total_ratio = g_cs_ratio_Oxygen->Eval(T_dm_before/1000.);
            else elas_total_ratio = g_cs_ratio_Iron->Eval(T_dm_before/1000.);
            if(noim>elas_total_ratio){
                type[counter] = 1;
                // find the input 2D PDF file (Tx_out:theta)
                double fi[30];
                for(int i=0;i<30;i++){
                    double noim1;
                    double noim2;
                    if((i+2)/2==(i+2)/2.) noim1 = 0.-(0.034+0.035)*i/2;
                    else noim1 = 0.-(0.034+0.035)*(i-1.)/2-0.034;

                    noim2 = pow(10,noim1);

                    if(T_dm_before/1000.-noim2>0) fi[i]=T_dm_before/1000.-noim2;
                    else fi[i]=-1.*T_dm_before/1000.+noim2;

                    if(i>0){
                        if(fi[i]<fi[i-1]) near = noim1;
                        else continue;
                    }
                    else continue;
                }
            }
            else{
                type[counter] = 0;
            }
            // the mean free path
            double step_in_m_correct;
            step_in_m_correct = A*mp/density/5.6e26/100./sigma;
            TString path = Form("exp(-1.*x/%.3e)/%.3e",step_in_m_correct,step_in_m_correct);
            TF1 *func_path = new TF1("func_path",path,0,step_in_m_correct*20.);

            // the elastic scattering angle PDF
            TF1 *f1 = new TF1("theta_pdf",ThetaCM_PDF,0,TMath::Pi(),7);
            f1->SetParameter(0,m_dm);
            f1->SetParameter(1,m_tgt);
            f1->SetParameter(2,m_dm);
            f1->SetParameter(3,m_tgt);
            f1->SetParameter(5,lambda);
            f1->SetParameter(6,m_med);
            // the QE 2D T_dm_out vs. theta PDF
            TString QE_name;
            if(m_med>99){

                QE_name = Form("../dat/xsec_data/dsigma_dTchipdtheta/%s/logmchi-%.3f_logTchi-%.3f_mS-%.1f.txt",element_one.c_str(),log10(m_dm/1000),near,m_med/1000);
            }else if(m_med>9){
                QE_name = Form("../dat/xsec_data/dsigma_dTchipdtheta/%s/logmchi-%.3f_logTchi-%.3f_mS-%.2f.txt",element_one.c_str(),log10(m_dm/1000),near,m_med/1000);
            }else{
                QE_name = Form("../dat/xsec_data/dsigma_dTchipdtheta/%s/logmchi-%.3f_logTchi-%.3f_mS-%.3f.txt",element_one.c_str(),log10(m_dm/1000),near,m_med/1000);
            }

            TGraph2D *QE_input = new TGraph2D(QE_name); 
            TH2D *QE_PDF = new TH2D("QE_PDF","",100,0,pow(10,near),100,0,180);
            for(int i=1;i<=100;i++){
                for(int j=1;j<=100;j++){
                    if(QE_input->Interpolate(i/100.*pow(10,near),180*j/100.)>0){
                        QE_PDF->SetBinContent(i,j,QE_input->Interpolate(i/100.*pow(10,near),180*j/100.));
                    }
                    else QE_PDF->SetBinContent(i,j,0);
                }
            }

            // The first scattering after one MFP
            if(counter==0){
                step_m[counter] = func_path->GetRandom();

                TVector3 disp = vec_before;
                disp.SetMag(step_m[counter]);
                v_pos_after = v_pos_before - disp;
                theta[counter] = vec_before.Theta();
                x[counter] = v_pos_after.X();
                y[counter] = v_pos_after.Y();
                z[counter] = v_pos_after.Z();
                T[counter] = T_dm_before;
                stepno[counter] = counter;
                vec_after = vec_before;
                T_dm_after = T_dm_before;
            }

            else{
                step_m[counter] = func_path->GetRandom();  

                // now determine the T_dm_out and the theta
                if(type[counter]==1){
                    QE_PDF->GetRandom2(T_dm_after,theta_cm);
                    T_dm_after = T_dm_after*1000.;
                }
                else{
                    f1->SetParameter(4,T_dm_before);
                    theta_cm = f1->GetRandom();
                    GetKinematics(m_dm, m_tgt, m_dm, m_tgt, T_dm_before, theta_cm, theta_lab_relative, T_dm_after);	
                }

                theta_each[counter] = theta_cm;
                theta_lab_relative = theta_cm;
                //now start from theta_lab_relative, calculate the real angle of this step in the lab
                vec_after = vec_before;
                TVector3 vec_z(0,0,1); //choose z as a convention
                TVector3 vec_axis = vec_z.Cross(vec_before); //rotation axis for theta
                if(!vec_axis.Mag()) vec_axis = TVector3(0,1,0);
                //now make a theta rotation again in the same xz plane by tb2
                vec_after.Rotate(theta_lab_relative, vec_axis); 
                //now make the azimuthal rotation by v2 by phi2
                vec_after.Rotate(gRandom->Uniform(0,TMath::Pi()*2), vec_before);

                //now assign values
                TVector3 disp = vec_after;
                disp.SetMag(step_m[counter]);
                v_pos_after = v_pos_before - disp;

                theta[counter] = vec_after.Theta();
                x[counter] = v_pos_after.X();
                y[counter] = v_pos_after.Y();
                z[counter] = v_pos_after.Z();
                T[counter] = T_dm_after;
                stepno[counter] = counter;

            }

            v_pos_before = v_pos_after;
            vec_before = vec_after;
            T_dm_before = T_dm_after;
            Hit = 0;
            T_reached = T[counter];     

            counter++;

            // break while the event over the earth surface
            if((pow(z[counter-1]*z[counter-1]+x[counter-1]*x[counter-1]+y[counter-1]*y[counter-1],0.5)>radius&&cross==0)||(pow(z[counter-1]*z[counter-1]+x[counter-1]*x[counter-1]+y[counter-1]*y[counter-1],0.5)<Core_radius)){
                f1->Delete();
                func_path->Delete();
                QE_input->Delete();
                QE_PDF->Delete();
                break;
            }
            // record Hit=1 while the event pass the inner circle first
            else if(pow(z[counter-1]*z[counter-1]+x[counter-1]*x[counter-1]+y[counter-1]*y[counter-1],0.5)<radius-JP_deep){
                f1->Delete();
                func_path->Delete();
                QE_input->Delete();
                QE_PDF->Delete();
                cross = 1;
                continue;
            }
            // recoid Hit=2 while the event back to the interlayer again
            else if(cross==1&&pow(z[counter-1]*z[counter-1]+x[counter-1]*x[counter-1]+y[counter-1]*y[counter-1],0.5)>radius-JP_deep){
                f1->Delete();
                func_path->Delete();
                QE_input->Delete();
                QE_PDF->Delete();
                cross = 2;
                break;
            }
            // else continue the circle
            else{
                f1->Delete();
                func_path->Delete();
                QE_input->Delete();
                QE_PDF->Delete();
                continue;
            }
        }
        Hit = cross;

        nstep = counter;
        outfile->cd();
        t3->Fill();

        evtno++;
    }
    outfile->cd();  
    t3->Write();
    outfile->Close();

    return 0;
}
