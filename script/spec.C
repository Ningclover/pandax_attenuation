#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <regex>
#include <sstream>

void spec(string file, TString nname)
{       
    std::regex pattern2("cs_([0-9]+.[0-9]+e-?[0-9]+)_([0-9]+.?[0-9]*)mev_ms([0-9]+.?[0-9]*)_([0-9]+)GeV");
    std::smatch match2;
    TFile *f_p4eff = new TFile("../dat/eff_RDQ_graph.root","read");
    TGraph *geff = (TGraph*)f_p4eff->Get("eff_nr");
    std::map<double,std::map<double,int>> masslist;
    if(std::regex_search(file, match2, pattern2)){
        if(match2.size() == 5){
            double masstmp = std::stof(match2[2].str());
            cout<<match2[1].str()<<endl;
            double xsectmp = std::stod(match2[1].str());
            double mediator = std::stof(match2[3].str());
            string th = match2[4];
            cout<<masstmp<<" "<<xsectmp<<endl;
            TGraph *g1 = new TGraph(file.c_str());
            g1->SetName(Form("g_mass_%.2f_xsec_%.5e",masstmp,xsectmp));
            TFile *fout = new TFile(nname,"recreate");         
            g1->Write();
            TH1F *h_sd = new TH1F("h_Adm",Form("h_mass_%.2f_xsec_%.5e",masstmp,xsectmp),2000,0.5,200);
            TH1F *h_sd_eff = new TH1F("h_Adm_p4",Form("h_eff_mass_%.2f_xsec_%.5e",masstmp,xsectmp),2000,0.5,200);
            for(int i=1;i<2001;i++)
            {   
                double val = g1->Eval(h_sd->GetBinCenter(i)*1e-6)*1e-6/1000/365; // per (kg day keV)
                h_sd->SetBinContent(i,val);
                double eff = geff->Eval(h_sd_eff->GetBinCenter(i));
                h_sd_eff->SetBinContent(i,val*eff);

            }
            double count = h_sd_eff->Integral("width")*2669*86; //2669 kg 86 day
            cout<<"mass = "<<masstmp<<"MeV xsec = "<<xsectmp<<" th = "<<th<<" counts = "<<count<<endl;
            if(count>0){
                masslist[masstmp][xsectmp] = 1; 
            }
            h_sd->Write();
            h_sd_eff->Write();
            fout->Close();
        }
    }

}
