{
    TH1F* h_tpc[9];
    TH1F* h_tpctof[9];
    TFile* file[9];
  
    //TF1* f1 = new TF1("erf", "[0]*exp(-pow([1]/x, [2]))", 0.2, 2.0);
    TF1* f1 = new TF1("erf", "[0]*TMath::Erf([1]*(x-[2]))+[3]", 0.15, 2.0);
    ofstream outfile("fit_para_proton_tof_eff.txt");
    Double_t par[4];
    for(Int_t i = 0; i < 9; i++){
        TString file_nm;
        file_nm.Form("../data/eff_proton_cen%d.root", i+1);
        file[i] = new TFile(file_nm.Data(), "read");
        TString h_tpc_nm; 
        TString h_tpctof_nm;
        h_tpc_nm.Form("h_tpc_cen%d", i+1);
        h_tpctof_nm.Form("h_tpctof_cen%d", i+1);
        h_tpc[i] = (TH1F*) file[i] -> Get(h_tpc_nm.Data());
        h_tpc[i] -> Sumw2();
        h_tpctof[i] = (TH1F*) file[i] -> Get(h_tpctof_nm.Data());
        h_tpctof[i] -> Sumw2();
        
        h_tpctof[i] -> Divide(h_tpc[i]);

        f1 -> SetParameter(0, 0.4);
        f1 -> SetParameter(1, 18);
        f1 -> SetParameter(3, 0);
        h_tpctof[i] -> Fit(f1, "ERWW0");
        f1 -> GetParameters(par);
        outfile << i << " " << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << std::endl;

        TString can_nm;
        can_nm.Form("can_%d", i);
        TCanvas* c_temp = new TCanvas("c_temp");
        c_temp -> SetLogy();
        h_tpctof[i] -> Draw("E");
        f1 -> SetParameters(par);
        f1 -> Draw("same");
        TString can_nm;
        can_nm.Form("tof_eff_cen%d.eps", i);
        c_temp -> SaveAs(can_nm.Data());
    
    }
}
