#pragma once
#include "Rtypes.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <TString.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <TApplication.h>
#include <TMath.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TF2.h>
#include <TMultiGraph.h>
#include <TFeldmanCousins.h>
#include <TRandom.h>
#include <Math/PdfFuncMathCore.h>
#include <TKey.h>

using namespace std;

// The mountain eval
TGraph2D *gr_mountain;
TGraph *gr_eng;

void GetNeutronKinematics(double T1, const double m1, const double m2,
        const double m3, const double m4, 
        const double theta_cm, const double phi_cm,
        const TVector3 vec_alpha,
        double& energy, TVector3& vec_neutron)
{
    double beta_cm = sqrt(T1*T1+2*m1*T1)/(T1+m1+m2);
    double w2 = pow(T1+m1+m2,2.)-(T1*T1+2*m1*T1);
    double w = sqrt(w2);
    double E3_cm = (w2+m3*m3-m4*m4)/2./w;
    double P3_cm = sqrt(E3_cm*E3_cm-m3*m3);
    double T3(0);
    TVector3 vec_neutron_in_alpha_coord(0,0,1);

    if(E3_cm<m3) {
        cout<<"below reaction threshold"<<endl;
        T3 = 0;
    } else {
        TVector3 vbeta(0,0,beta_cm);
        TVector3 p3vec(0,0,1);
        // magnitude = sqrt(x2+y2+z2)
        p3vec.SetMag(P3_cm);
        // set polar angle
        p3vec.SetTheta(theta_cm);
        // set azimuth angle
        p3vec.SetPhi(phi_cm);

        // Four dimensional vector, here for position and energy 
        TLorentzVector lvec(p3vec,E3_cm);
        // lorentz transfer with beta = v/c
        lvec.Boost(vbeta);
        // return the fourth parameter like energy in (px,py,pz,E) and Mag()
        T3 = lvec.E()-lvec.M();
        // lvec.Vect() return the three dimentional vector
        vec_neutron_in_alpha_coord = lvec.Vect(); //relative to alpha vector
    }
    energy = T3;

    //start with neutron vector the same as alpha vector
    vec_neutron = vec_alpha;
    TVector3 vec_z(0,0,1);
    TVector3 vec_axis = vec_z.Cross(vec_alpha);

    if(vec_axis.Mag()>0){
        //now make a theta rotation again in the same xz plane by tb2
        vec_neutron.Rotate(vec_neutron_in_alpha_coord.Theta(), vec_axis); 
        //now make the azimuthal rotation by v2 by phi2
        vec_neutron.Rotate(vec_neutron_in_alpha_coord.Phi(), vec_alpha);
    } else {
        TVector3 v(0,1,0); //around y axis
        vec_neutron.Rotate(vec_neutron_in_alpha_coord.Theta(), v);
        vec_neutron.Rotate(vec_neutron_in_alpha_coord.Phi(), vec_alpha);
    }
}

//use MeV as the "natural unit"
//Theta_cm, weighted by the form factor
Double_t ThetaCM_PDF(Double_t *x, Double_t *par)
{
    Double_t m1 = par[0];
    Double_t m2 = par[1];
    Double_t m3 = par[2];
    Double_t m4 = par[3];
    Double_t T1 = par[4];
    Double_t lambda = par[5];
    Double_t ms = par[6]; 
    Double_t theta_cm = x[0]; //argument in 0,2pi

    double beta_cm = sqrt(T1*T1+2*m1*T1)/(T1+m1+m2);
    double w2 = pow(T1+m1+m2,2.)-(T1*T1+2*m1*T1);
    double w = sqrt(w2);
    double E3_cm = (w2+m3*m3-m4*m4)/2./w;
    double P3_cm = sqrt(E3_cm*E3_cm-m3*m3);

    double T3(0);

    TVector3 v1(0,0,1);
    TVector3 v3(0,0,1);

    if(E3_cm<m3) {
        cout<<"below reaction threshold"<<endl;
        T3 = 0;
    } else {
        TVector3 vbeta(0,0,beta_cm);
        TVector3 p3vec(0,0,1);
        p3vec.SetMag(P3_cm);
        p3vec.SetTheta(theta_cm);

        TLorentzVector lvec(p3vec,E3_cm);
        lvec.Boost(vbeta);
        T3 = lvec.E()-lvec.M();
        v3 = lvec.Vect(); //final momentum vector
    }

    v1.SetMag(sqrt(T1*T1+2*T1*m1));
    TLorentzVector lv1(v1,T1+m1);

    v3.SetMag(sqrt(T3*T3+2*T3*m3));
    TLorentzVector lv3(v3,T3+m3);

    TLorentzVector four_mon_trans = lv1-lv3;

    Double_t Q = four_mon_trans.Mag();
    Double_t ff = 1./pow(1+Q*Q/lambda/lambda,2);
    // Add the DM form factor for up-philic scalar mediator
    Double_t f_dm= (4.*pow(m1,2)+Q*Q)*pow(pow(ms,2),2)*(4.*pow(m1,2)+Q*Q)/(16.*pow(m1*m2*(pow(ms,2)+Q*Q),2));
    ff = ff*f_dm;
    return ff;
}

void GetKinematics(const double m1, const double m2, const double m3, const double m4, double T1, const double theta_cm, double& theta_lab, double_t& t3)
{
    //repeat basic calculation
    double beta_cm = sqrt(T1*T1+2*m1*T1)/(T1+m1+m2);
    double w2 = pow(T1+m1+m2,2.)-(T1*T1+2*m1*T1);
    double w = sqrt(w2);
    double E3_cm = (w2+m3*m3-m4*m4)/2./w;
    double P3_cm = sqrt(E3_cm*E3_cm-m3*m3);


    TVector3 v1(0,0,1);
    TVector3 v3(0,0,1);

    if(E3_cm<m3) {
        cout<<"below reaction threshold"<<endl;
        t3 = 0;
    } else {
        TVector3 vbeta(0,0,beta_cm);
        TVector3 p3vec(0,0,1);
        p3vec.SetMag(P3_cm);
        p3vec.SetTheta(theta_cm);

        TLorentzVector lvec(p3vec,E3_cm);
        lvec.Boost(vbeta);
        t3 = lvec.E()-lvec.M();
        v3 = lvec.Vect(); //final momentum vector
    }
    theta_lab = v3.Angle(v1);
    cout<<endl;
    return;
}

double myfunc_eng(double *xx, double *)
{
    return gr_eng->Eval(xx[0]);
}

TF1* get_eng_spec(const TString filename)
{ 
    gr_eng = new TGraph(filename);
    TF1* func_eng = new TF1("f1",myfunc_eng, 5e-5, 1);
    func_eng->SetNpx(1000);
    return func_eng;
}

double myfunc_mountain(double *xx, double *){
    return gr_mountain->Interpolate(xx[0],xx[1]);
}

double cal_sigma(double mx, double ms, double gn=0.01, double gchi=0.01){
    // This is the initial cross section related to xsec_data
    // All in GeV
    double mp = 0.938;
    double GeV2cm =1/2.568*1e-27;
    double mun = mp*mx/(mp+mx);
    double sigma_n = mun*mun*pow(gn*gchi,2)/3.14/pow(ms,4)*GeV2cm;
    return sigma_n;
}
