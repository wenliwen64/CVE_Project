void QA(){
    gROOT -> LoadMacro("v0dst.C+");
    TFile f1 = new TFile("auau200GeV_run11_ks_la.total.picodst.root", "read");
    TTree* tree = f1 -> Get("McV0PicoDst");

    v0dst dst(tree);

    TFile f2 = new TFile("QA_auau200GeV_run11_ks_la.total.histo.root", "recreate");
    f2 -> cd();

    TH1F* h_mcv0pt = new TH1F("h_mcv0pt", "h_mcv0pt", 100, 0, 10.0);
    TH1F* h_v0pt = new TH1F("h_v0pt", "h_v0pt", 100, 0, 10.0);

    TH1F* h_mcv0x = new TH1F("h_mcv0x", "h_mcv0x", 100, 0, 10);
    TH1F* h_v0x = new TH1F("h_v0x", "h_v0x", 100, 0, 10);

    TH1F* h_mcv0y = new TH1F("h_mcv0y", "h_mcv0y", 100, 0, 10);
    TH1F* h_v0y = new TH1F("h_v0y", "h_v0y", 100, 0, 10);

    TH1F* h_mcv0z = new TH1F("h_mcv0z", "h_mcv0z", 100, 0, 10);
    TH1F* h_v0z = new TH1F("h_v0z", "h_v0z", 100, 0, 10);

    TH1F* h_mcv0dca = new TH1F("h_mcv0dca", "h_mcv0dca", 100, 0, 5);
    TH1F* h_v0dca = new TH1F("h_v0dca", "h_v0dca", 100, 0, 5);

    int nEntries = tree -> nEntries();
    for(int i = 0; i < nEntries; i++){
        dst.LoadEntry(i);
        h_mcv0pt -> Fill(dst.mcv0pt);
        h_v0pt -> Fill(dst.v0pt);
   
        h_mcv0x -> Fill(dst.v0x);
        h_mcv0y -> Fill(dst.v0y);
        h_mcv0z -> Fill(dst.v0z);

        h_mcv0dca -> Fill();
    }

}
