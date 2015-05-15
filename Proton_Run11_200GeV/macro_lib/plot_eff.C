void plot_eff(string par_nm){
    //Get the data and kick out the zero data TODO: add some case switch code
    ifstream infile_fp("../data/weight_lambda_exp.txt");
    ifstream infile_fp_scale("../data/weight_lambda_exp_scale.txt");
    ifstream infile_exp("../data/weight_lambda_fp.txt");
    ifstream infile_exp_scale("../data/weight_lambda_fp_scale.txt");

    int centbin;
    int ptbin;
    float eff;
    float err;
    float ptbd[8] = { 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.15};
    float efficiency_fp[9][8];
    float efficiency_exp[9][8];
    float efficiency_fp_scale[9][8];
    float efficiency_exp_scale[9][8];
    float efferror_fp[9][8];
    float efferror_exp[9][8];
    float efferror_fp_scale[9][8];
    float efferror_exp_scale[9][8];

    float err_x[8] = {0, 0, 0, 0, 0, 0, 0, 0};

    while(infile_fp >> centbin){
	infile_fp >> ptbin;
	infile_fp >> eff;
	infile_fp >> err;
	efficiency_fp[centbin][ptbin-2] = eff;
	efferror_fp[centbin][ptbin-2] = err;
    }
    while(infile_exp >> centbin){
	infile_exp >> ptbin;
	infile_exp >> eff;
	infile_exp >> err;
	efficiency_exp[centbin][ptbin-2] = eff;
	efferror_exp[centbin][ptbin-2] = err;
    }
    while(infile_fp_scale >> centbin){
	infile_fp_scale >> ptbin;
	infile_fp_scale >> eff;
	infile_fp_scale >> err;
	efficiency_fp_scale[centbin][ptbin-2] = eff;
	efferror_fp_scale[centbin][ptbin-2] = err;
    }
    while(infile_exp_scale >> centbin){
	infile_exp_scale >> ptbin;
	infile_exp_scale >> eff;
	infile_exp_scale >> err;
	efficiency_exp_scale[centbin][ptbin-2] = eff;
	efferror_exp_scale[centbin][ptbin-2] = err;
    }

//======================================Plot 1. Exp. vs. Flat

    TCanvas* can = new TCanvas("c1", "c1", 800, 600);
    can -> SetLogy();
    TMultiGraph* mg = new TMultiGraph();

    for(int i = 0; i < 9; i++){
	TGraphErrors* gr_exp = new TGraphErrors(8, ptbd, efficiency_fp_scale[i], NULL, efferror_fp_scale[i]);
	gr_exp -> SetMarkerSize(1.0);
	gr_exp -> SetMarkerColor(2);
	gr_exp -> SetMarkerStyle(29);
	mg -> Add(gr_exp);
    }

    for(int i = 0; i < 9; i++){
	TGraphErrors* gr_fp = new TGraphErrors(8, ptbd, efficiency_exp_scale[i], NULL, efferror_exp_scale[i]);
	gr_fp -> SetMarkerSize(1.0);
	gr_fp -> SetMarkerColor(1);
	gr_fp -> SetMarkerStyle(29);
	mg -> Add(gr_fp);
    }

    can -> cd();
    char plot_title[100];
    if(par_nm == "Lambda"){
	mg -> SetTitle("#Lambda Efficiency in AuAu 200GeV");
    }
    else{
	mg -> SetTitle("#bar{#Lambda} Efficiency in AuAu 200GeV");
    }

    mg -> Draw("ap");
    mg->GetXaxis()->SetTitle("Pt(GeV/c)");
    mg->GetXaxis()->SetRangeUser(0.6, 2.5);
    mg->GetYaxis()->SetTitle("Efficiency");
    mg->GetYaxis()->SetRangeUser(pow(10, -11), 1000);
    gPad->Update();

    TGraphErrors* gr_fp = (TGraphErrors*)mg -> GetListOfGraphs() -> At(1);    //Get one of the exponentially distributed efficiency plot from multigraph list.
    TGraphErrors* gr_exp = (TGraphErrors*)mg -> GetListOfGraphs() -> At(12);    //Get one of the flat distributed efficiency plot from multigraph.
    TLegend* leg = new TLegend(0.45, 0.75, 0.7, 0.87);
    leg -> AddEntry(gr_exp, "exponential", "p");
    leg -> AddEntry(gr_fp, "flat", "p");
    leg -> SetBorderSize(0);
    leg -> SetFillColor(0);
    leg -> SetTextSize(0.045);
    leg -> Draw("same");

    char fig_nm_eps[100];
    char fig_nm_jpg[100];
    sprintf(fig_nm_eps, "plot_exp_vs_fp_%s.eps", par_nm.c_str());
    sprintf(fig_nm_jpg, "plot_exp_vs_fp_%s.jpg", par_nm.c_str());
    can -> SaveAs(fig_nm_eps);
    can -> SaveAs(fig_nm_jpg);

    can -> Update();

//===========================Plot 2. Fitting Exponential Result with Error Function

    TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);  
    c2 -> cd();
    c2 -> SetLogy();
    TMultiGraph* mg1 = new TMultiGraph(); 
    char mg1_title[100];
    if(par_nm == "Lambda"){
	mg1 -> SetTitle("#Lambda Efficiency in AuAu 200GeV(Exp.)");
    }
    else{
	mg1 -> SetTitle("#bar{#Lambda} Efficiency in AuAu 200GeV");
    }

    TF1* erf = new TF1("erf", "[0]*TMath::Erf(x-[1])+[2]", 0.6, 3.0);   
    double par[9][3];
    for(int i = 0; i < 9; i++){
	TGraphErrors* gr_exp = new TGraphErrors(8, ptbd, efficiency_exp[i], NULL, efferror_exp[i]);
	TGraphErrors* gr_exp_scale = new TGraphErrors(8, ptbd, efficiency_exp_scale[i], NULL, efferror_exp_scale[i]);
	gr_exp_scale -> SetMarkerSize(1.0);
	gr_exp_scale -> SetMarkerColor(2);
	gr_exp_scale -> SetMarkerStyle(29);
	mg1 -> Add(gr_exp_scale);
	gr_exp -> Fit(erf, "R0");
	erf -> GetParameters(par[i]); 
	std::cout<<"[0]="<<par[i][0]<<" [1]="<<par[i][1]<<" [2]="<<par[i][2]<<std::endl;
    }
    mg1 -> Draw("ap");
    mg1 -> GetXaxis() -> SetTitle("Pt(GeV/c)");
    mg1 -> GetXaxis() -> SetRangeUser(0.6, 2.5);
    mg1 -> GetYaxis() -> SetTitle("Efficiency");
    mg1 -> GetYaxis() -> SetRangeUser(pow(10, -11), 1000);
    gPad -> Update();

    char eff_fit_output[100];
    sprintf(eff_fit_output, "../output/eff_fit_par_%s.dat", par_nm.c_str());
    ofstream eff_fit(eff_fit_output);
    for(int i = 0; i < 9; i++){
	TF1* erf_temp = new TF1("erf_temp", "[0]*TMath::Erf(x-[1])+[2]", 0.6, 3.0);
	erf_temp -> SetParameter(0, pow(10, -i)*par[i][0]); 
	erf_temp -> SetParameter(1, par[i][1]); 
	erf_temp -> SetParameter(2, pow(10, -i)*par[i][2]); 
	erf_temp -> SetLineColor(1);
	erf_temp -> SetLineStyle(2);
	erf_temp -> Draw("sames");
	eff_fit << par[i][0] << " " << par[i][1] << " " << par[i][2] <<std::endl;
    }
    eff_fit.close();
    gPad -> Update();

    TLegend* leg1 = new TLegend(0.45, 0.75, 0.7, 0.87);
    TGraphErrors* graph_temp_1 = (TGraphErrors*)mg1 -> GetListOfGraphs() -> At(2);
    TF1* erf_temp = new TF1("erf_temp", "[0]*TMath::Erf(x-[1])+[2]", 0.5, 3.0);
    erf_temp -> SetLineColor(1);
    erf_temp -> SetLineStyle(2);

    leg1 -> AddEntry(graph_temp_1, "Efficiency", "p");
    leg1 -> AddEntry(erf_temp, "Error Function", "l");
    leg1 -> SetBorderSize(0);
    leg1 -> Draw("same");

    char fit_fig_nm_eps[100];
    char fit_fig_nm_jpg[100];
    sprintf(fit_fig_nm_eps, "../figures/%s_fit_plot.eps", par_nm.c_str()); 
    sprintf(fit_fig_nm_jpg, "../figures/%s_fit_plot.jpg", par_nm.c_str()); 
    c2 -> SaveAs(fit_fig_nm_eps);
    c2 -> SaveAs(fit_fig_nm_jpg);

}
