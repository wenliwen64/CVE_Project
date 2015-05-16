void FitEff_cen(int cen=5){

gStyle->SetOptFit(1);
gStyle->SetOptDate(0);
gStyle->SetOptStat(0);
gStyle->SetFillColor(0);
char filename[200];
sprintf(filename,"cen%d.eff_200_p.root",cen);
TFile *ff = new TFile(filename);


TH1D* mc = (TH1*)ff->Get("Hist_mc_pt_cen");
TH1D* rc = (TH1*)ff->Get("Hist_rc_pt_cen");
rc->Divide(mc);
rc->SetTitle("");
rc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
rc->GetYaxis()->SetTitle("Tracking efficiency");
rc->GetXaxis()->SetRangeUser(0,8);

//TF1 *func = new TF1("func","[0]*TMath::Erf((x-[1])*[2])+[0]",0,5);
//TF1 *func = new TF1("func","[0]*tanh([1]*(pow(x,[2])-[3]))",0,6);
//func->SetParameters(0.8,250,0.002,1);
//func->SetParameters(0.81,370,0.002,1);

TF1 *func = new TF1("func","[0]*exp(-pow([1]/x,[2]))",0,6);
func->SetParameters(0.80,0.12,6);
//func->SetParameters(1.,0.10,2.0);

func->SetLineStyle(1);
func->SetLineColor(2);
rc->Fit("func","E","",0,6);
gPad->SetRightMargin(0.03);
gPad->SetTopMargin(0.03);

TF1 *func_39 = new TF1("func_39","[0]*exp(-pow([1]/x,[2]))",0,6);
func_39->SetParameters(0.862022,0.170073,2.23007);
func_39->SetLineColor(6);
func_39->SetLineStyle(4);
//func_39->Draw("same");

TF1 *func_run4 = new TF1("func_run4","[0]*exp(-pow([1]/x,[2]))",0,6);
func_run4->SetParameters(0.8232,0.09349,1.517);
func_run4->SetLineColor(4);
func_run4->SetLineStyle(2);
//func_run4->Draw("same");

TF1 *func_run7 = new TF1("func_run7","[0]*exp(-pow([1]/x,[2]))",0,6);
func_run7->SetParameters(0.723741,0.2886,1.63058);
func_run7->SetLineColor(8);
//func_run7->Draw("same");

   TLatex *   tex = new TLatex(0.4,1.0,"30% - 40%");
   tex->SetTextSize(0.06);
   tex->SetTextColor(1);
   tex->Draw();

   TLegend* legend = new TLegend(0.40, 0.17, 0.70, 0.50);
   legend->SetFillColor(0);
   legend->SetTextSize(0.055);
   legend->SetLineColor(0);
   legend->SetBorderSize(0.000001);
   legend->AddEntry(func_run4, "run4  200 GeV Au+Au","l");
   legend->AddEntry(func_run7, "run7  200 GeV Au+Au","l");
   legend->AddEntry(func_39, "run10 39 GeV Au+Au","l");
   legend->AddEntry(func, "run11 200 GeV Au+Au","l");
   // legend->Draw();
}
