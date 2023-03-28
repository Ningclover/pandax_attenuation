#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <regex>
#include <sstream>

void spec_changeBr(TString file, TString nname,double masstmp,double xsectmp,double mediator,double br,double gu)
{   
    TFile *f_p4eff = new TFile("../dat/eff_RDQ_graph.root","read");
    TGraph *geff = (TGraph*)f_p4eff->Get("eff_nr");
                TGraph *g1 = new TGraph(file);
                g1->SetName(Form("g_mass_%.2f_xsec_%.5e",masstmp,xsectmp));
                TFile *fout = new TFile(nname,"recreate");         
                g1->Write();
                TH1F *h_sd = new TH1F("h_Adm",Form("h_mass_%.2f_xsec_%.5e",masstmp,xsectmp),2000,0.5,200);
                TH1F *h_sd_eff = new TH1F("h_Adm_p4",Form("h_eff_mass_%.2f_xsec_%.5e",masstmp,xsectmp),2000,0.5,200);
                for(int i=1;i<2001;i++)
                {   
                    double val = g1->Eval(h_sd->GetBinCenter(i)*1e-6)*1e-6/1000/365*br/1e-5; // per (kg day keV) change br
                    h_sd->SetBinContent(i,val);
                    double eff = geff->Eval(h_sd_eff->GetBinCenter(i));
                    h_sd_eff->SetBinContent(i,val*eff);

                }
                double count = h_sd_eff->Integral("width")*2669*86; //2669 kg 86 day
                cout<<"br = "<<br<<"  mass = "<<masstmp<<"MeV  meditaor = "<<mediator<<"MeV  xsec = "<<xsectmp<<" counts = "<<count<<endl;

                h_sd->Write();
                h_sd_eff->Write();
                fout->Close();

}
