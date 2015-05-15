#include "TFile.h"
#include <fstream>
#include <iostream>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TNtuple.h"

void Calc_Eff_200(int cen) {
const float PI = TMath::Pi();
const float mean = 0.01166;
const float rms = 0.05312;
TFile *f1 = new TFile("cen6.ntuple_result_v2_11_subEP_dca1_pt015_eta1_noP_pionDca1_higherEP.root");
TH1* h_pt = (TH1*)f1->Get("Hist_Pt");

float real[9][7] = {{1.68441e+07,1.94857e+07,1.91433e+07,1.76339e+07,1.57031e+07,1.3709e+07,1.18246e+07},
			{3.87595e+07,4.43888e+07,4.34577e+07,4.00144e+07,3.57267e+07,3.13252e+07,2.71556e+07},
			{8.04378e+07,9.2009e+07,9.01896e+07,8.33911e+07,7.47927e+07,6.59201e+07,5.74707e+07},
			{1.28651e+08,1.47865e+08,1.45606e+08,1.3528e+08,1.2199e+08,1.08084e+08,9.47436e+07},
			{1.9998e+08,2.31542e+08,2.29316e+08,2.14124e+08,1.94076e+08,1.72786e+08,1.52079e+08},
			{2.91541e+08,3.40818e+08,3.39637e+08,3.1877e+08,2.90097e+08,2.59203e+08,2.29072e+08},
			{3.94517e+08,4.66796e+08,4.68155e+08,4.41467e+08,4.03254e+08,3.61425e+08,3.20203e+08},
			{2.3895e+08,2.85928e+08,2.88314e+08,2.72767e+08,2.49768e+08,2.24261e+08,1.98897e+08},
			{2.70925e+08,3.27361e+08,3.31383e+08,3.14285e+08,2.88208e+08,2.5899e+08,2.29898e+08}};
float emrc[9][7] = {{127.287,180.108,180.46,203.105,201.325,185.72,183.217},
			{111.057,137.61,163.772,151.6,161.841,157.728,178.696},
			{157.959,196.018,206.222,184.609,183.632,238.975,225.804},
			{294.14,355.613,355.514,396.844,397.958,385.601,425},
			{499,551,666,616,706,750,677},
			{715,916,988,992,1031,1047,1070},
			{1045,1247,1393,1414,1523,1589,1527},
			{653,853,898,928,1016,1073,996},
			{863,975,1199,1263,1262,1258,1314}};
const int cenDef[9] = {10, 22, 43, 76, 125, 193, 281, 396, 466};

float Eweight = 1;

        TChain* chain = new TChain("McV0PicoDst");
        int nfile = 0;
        nfile += chain->Add("output_p/*.root");
	nfile += chain->Add("output_pbar/*.root");
//	nfile += chain->Add("output_test/*.root");
//	nfile += chain->Add("output_km/*.root");
//        nfile += chain->Add("output_kp/*.root");
                cout <<"Added "<<nfile<<" files"<<endl;
                cout<<"# entries in chain: "<<chain->GetEntries()<<endl;

	char fname_out[200];
	sprintf(fname_out,"cen%d.eff_200_p.root",cen);
        TFile fout(fname_out,"RECREATE");

	TH1D* Hist_RefMult = new TH1D("Hist_RefMult","Hist_RefMult",500,-0.5,499.5);
	TH2D* Hist_mc_pt = new TH2D("Hist_mc_pt","Hist_mc_pt",28,0,560,200,0,10);
	TH1D* Hist_mc_pt_cen = new TH1D("Hist_mc_pt_cen","Hist_mc_pt_cen",200,0,10);
        TH1D* Hist_mc_pt_L = new TH1D("Hist_mc_pt_L","Hist_mc_pt_Lokesh",200,0,10);
	TH1D* Hist_mc_eta = new TH1D("Hist_mc_eta","Hist_mc_eta",200,-1,1);
	TH1D* Hist_rc_flag = new TH1D("Hist_rc_flag","Hist_rc_flag",3000,-1000+0.5,2000.5);
	TH2D* Hist_rc_pt = new TH2D("Hist_rc_pt","Hist_rc_pt",28,0,560,200,0,10);
	TH1D* Hist_rc_pt_cen = new TH1D("Hist_rc_pt_cen","Hist_rc_pt_cen",200,0,10);
	TH1D* Hist_rc_pt_scale = new TH1D("Hist_rc_pt_scale","Hist_rc_pt_scale",200,0,10);
	TH1D* Hist_rc_pt_L = new TH1D("Hist_rc_pt_L","Hist_rc_pt_Lokesh",200,0,10);
	TH1D* Hist_rc_eta = new TH1D("Hist_rc_eta","Hist_rc_eta",200,-1,1);
	TH1D* Hist_rc_dca = new TH1D("Hist_rc_dca","Hist_rc_dca",200,0,20);
	TH1D* Hist_rc_dca_scale = new TH1D("Hist_rc_dca_scale","Hist_rc_dca_scale",200,0,20);
	TH1D* Hist_rc_nfit = new TH1D("Hist_rc_nfit","Hist_rc_nfit",100,0,100);
	TH1D* Hist_rc_rat = new TH1D("Hist_rc_rat","Hist_rc_rat",100,0,2);
	TH1D* Hist_dphi = new TH1D("Hist_dphi","Hist_dphi",128,-3.2,3.2);
        TH1D* Hist_dphi1 = new TH1D("Hist_dphi1","Hist_dphi1",128,-3.2,3.2);
        TH1D* Hist_dphi2 = new TH1D("Hist_dphi2","Hist_dphi2",128,-3.2,3.2);
        TH1D* Hist_dphi3 = new TH1D("Hist_dphi3","Hist_dphi3",128,-3.2,3.2);
        TH1D* Hist_dphi4 = new TH1D("Hist_dphi4","Hist_dphi4",128,-3.2,3.2);
        TH1D* Hist_dphi_pion = new TH1D("Hist_dphi_pion","Hist_dphi_pion",128,-3.2,3.2);
        TH1D* Hist_ddphi = new TH1D("Hist_ddphi","Hist_ddphi",128,-3.2,3.2);
        TProfile* p_dphi = new TProfile("p_dphi","p_dphi",4,0.5,4.5,-100,100);
        TProfile* p_dphi_pt = new TProfile("p_dphi_pt","p_dphi_pt",200,0,10,-100,100);
        TProfile* p_dphi_L = new TProfile("p_dphi_L","p_dphi_L",4,0.5,4.5,-100,100);
        TProfile* p_dphi_LM= new TProfile("p_dphi_LM","p_dphi_LM",4,0.5,4.5,-100,100);
        TProfile* p_dphi_M = new TProfile("p_dphi_M","p_dphi_M",4,0.5,4.5,-100,100);
        TProfile* p_dphi_RM= new TProfile("p_dphi_RM","p_dphi_RM",4,0.5,4.5,-100,100);
        TProfile* p_dphi_R = new TProfile("p_dphi_R","p_dphi_R",4,0.5,4.5,-100,100);
	TProfile* p_ddphi = new TProfile("p_ddphi","p_ddphi",4,0.5,4.5,-100,100);
        TProfile* p_ddphi_L = new TProfile("p_ddphi_L","p_ddphi_L",4,0.5,4.5,-100,100);
        TProfile* p_ddphi_LM= new TProfile("p_ddphi_LM","p_ddphi_LM",4,0.5,4.5,-100,100);
        TProfile* p_ddphi_M = new TProfile("p_ddphi_M","p_ddphi_M",4,0.5,4.5,-100,100);
        TProfile* p_ddphi_RM= new TProfile("p_ddphi_RM","p_ddphi_RM",4,0.5,4.5,-100,100);
        TProfile* p_ddphi_R = new TProfile("p_ddphi_R","p_ddphi_R",4,0.5,4.5,-100,100);
        TH1D* Hist_asym = new TH1D("Hist_asym","Hist_asym",600,-1.5+0.0025,1.5+0.0025);
        TH1D* Hist_asym_L = new TH1D("Hist_asym_L","Hist_asym_L",600,-1.5+0.0025,1.5+0.0025);
        TH1D* Hist_asym_LM= new TH1D("Hist_asym_LM","Hist_asym_LM",600,-1.5+0.0025,1.5+0.0025);
        TH1D* Hist_asym_M = new TH1D("Hist_asym_M","Hist_asym_M",600,-1.5+0.0025,1.5+0.0025);
        TH1D* Hist_asym_RM= new TH1D("Hist_asym_RM","Hist_asym_RM",600,-1.5+0.0025,1.5+0.0025);
        TH1D* Hist_asym_R = new TH1D("Hist_asym_R","Hist_asym_R",600,-1.5+0.0025,1.5+0.0025);

	Int_t nentries = chain->GetEntries();
	for(int i = 0; i < nentries; i++){

                if((i+1)%1000==0) cout<<"Processing entry == "<< i+1 <<" == out of "<<nentries<<".\n";
                chain->GetEntry(i);

                TLeaf* leaf_RunId   = chain->GetLeaf("runnumber");
		TLeaf* leaf_RefMult = chain->GetLeaf("nrefmult");
		TLeaf* leaf_Px	    = chain->GetLeaf("primvertexX");
		TLeaf* leaf_Py      = chain->GetLeaf("primvertexY");
		TLeaf* leaf_Pz      = chain->GetLeaf("primvertexZ");
		TLeaf* leaf_ZDCrate = chain->GetLeaf("zdcrate");
		TLeaf* leaf_nmcv0   = chain->GetLeaf("nmcv0");
		TLeaf* leaf_nv0	    = chain->GetLeaf("nv0");
                TLeaf* leaf_Np      = chain->GetLeaf("numberP");
                TLeaf* leaf_Nn      = chain->GetLeaf("numberN");

		int Run 	    = leaf_RunId->GetValue(0);
		int RefMult	    = leaf_RefMult->GetValue(0);
		float PVtxz	    = leaf_Pz->GetValue(0);
		int NmcTracks	    = leaf_nmcv0->GetValue(0);
		int NrcTracks	    = leaf_nv0->GetValue(0);
                int Np              = leaf_Np->GetValue(0);
                int Nn              = leaf_Nn->GetValue(0);

		Hist_RefMult->Fill(RefMult);
                int Centrality  = 0;
                for(int j=0;j<9;j++) if(RefMult>cenDef[j]) Centrality = j+1;
                if(RefMult) {
                        float gM = 0.9995 + 21.89/(4.191*RefMult-18.17) - 2.723e-5*(4.191*RefMult-18.17);
                        Eweight = gM + 0.0009326*(gM-1)*PVtxz;
                }
                if(cen && Centrality != cen) continue;

                TLeaf* leaf_mc_Pt        = chain->GetLeaf("mcv0pt");
		TLeaf* leaf_mc_Pz        = chain->GetLeaf("mcv0pz");
		TLeaf* leaf_mc_Px        = chain->GetLeaf("mcv0px");
		TLeaf* leaf_mc_Py        = chain->GetLeaf("mcv0py");
		TLeaf* leaf_mc_id        = chain->GetLeaf("mcv0id");

		for(int trk = 0; trk < NmcTracks; trk++) {
			float mc_pt = leaf_mc_Pt->GetValue(trk);
			float mc_pz = leaf_mc_Pz->GetValue(trk);
			float mc_theta = atan2(mc_pt,mc_pz);
			float mc_eta = -log(tan(mc_theta/2));
			Hist_mc_eta->Fill(mc_eta, Eweight);
			if(mc_eta>-0.5 && mc_eta<0.5) {Hist_mc_pt->Fill(RefMult,mc_pt);}
			if(mc_eta>-1 && mc_eta<1) {Hist_mc_pt_cen->Fill(mc_pt);}
			if(mc_eta>-0.1 && mc_eta<0.1) Hist_mc_pt_L->Fill(mc_pt, Eweight);
		}

                TLeaf* leaf_rc_Pt        = chain->GetLeaf("v0pt");
		TLeaf* leaf_rc_Phi	 = chain->GetLeaf("v0phi");
		TLeaf* leaf_rc_Eta	 = chain->GetLeaf("v0eta");
                TLeaf* leaf_rc_Charge    = chain->GetLeaf("v0charge");
		TLeaf* leaf_rc_flag	 = chain->GetLeaf("v0flag");
		TLeaf* leaf_rc_dca	 = chain->GetLeaf("v0dca");
		TLeaf* leaf_rc_Nfithits  = chain->GetLeaf("v0nfithits");
		TLeaf* leaf_rc_Nmaxhits  = chain->GetLeaf("v0nmaxhits"); 
		TLeaf* leaf_rc_Ndedxhits = chain->GetLeaf("v0ndedxhits");
                TLeaf* leaf_rc_id        = chain->GetLeaf("v0mcid");

               for(int trk = 0; trk < NrcTracks; trk++) {
			int flag = leaf_rc_flag->GetValue(trk);
			Hist_rc_flag->Fill(flag,Eweight);
			if(flag<0 || flag>1000) continue;
			float dca = leaf_rc_dca->GetValue(trk);
			float rc_pt = leaf_rc_Pt->GetValue(trk);
			float eta = leaf_rc_Eta->GetValue(trk);
			int nFitHits = leaf_rc_Nfithits->GetValue(trk);
			int nMaxHits = leaf_rc_Nmaxhits->GetValue(trk);
			float ratio = float(nFitHits)/nMaxHits;
			int ndEdxHits = leaf_rc_Ndedxhits->GetValue(trk);
			Hist_rc_nfit->Fill(nFitHits, Eweight);
			Hist_rc_rat->Fill(ratio, Eweight);
			if(nFitHits<20 || nFitHits>50) continue;
			if(ratio<0.52 || ratio>1.05) continue;
				if(dca>2) continue;

				int ptbin = 0;
				if(rc_pt>0.2) ptbin = 1;
				if(rc_pt>0.25) ptbin = 2;
				if(rc_pt>0.3) ptbin = 3;
				if(rc_pt>0.35) ptbin = 4;
				if(rc_pt>0.4) ptbin = 5;
				if(rc_pt>0.45) ptbin = 6;
				float ptW = 0.001*real[cen-1][ptbin]/emrc[cen-1][ptbin];
//			if((i+1)%1000==0 && rc_pt<0.3 && rc_pt>0.25)  cout<<ptW<<endl;
			Hist_rc_eta->Fill(eta, Eweight);
			
			if(rc_pt>0.15 && rc_pt<0.5) {Hist_rc_dca->Fill(dca, Eweight);Hist_rc_dca_scale->Fill(dca, Eweight*ptW);}
                        if(eta>-0.5 && eta<0.5) {Hist_rc_pt->Fill(RefMult,rc_pt);}

                        if(dca<2 && ndEdxHits>14) {
                        if(eta>-1 && eta<1) Hist_rc_pt_cen->Fill(rc_pt);
			Hist_rc_pt_scale->Fill(rc_pt, Eweight*ptW);
}
			if(dca<3 && nFitHits>=25 &&ndEdxHits>15 && eta>-0.1 && eta<0.1) Hist_rc_pt_L->Fill(rc_pt, Eweight);
                }
///////ONLY ABOVE THIS LINE IS NEEDED FOR EFFICIENCY/////////////////////////////////////////////////// 	
		float Mpt[30], Mdca[30], MndEdxHits[30], Mdphi[30], Mweight[30];
		float integral = (float)h_pt->Integral(1,40);
                float bin_mean = integral/40.;
		int Mnmatch = 0;
		for(int trk = 0; trk < NrcTracks; trk++) {
                        int flag = leaf_rc_flag->GetValue(trk);
			if(flag<0 || flag>1000) continue;
                        float dca = leaf_rc_dca->GetValue(trk);
                        float rc_pt = leaf_rc_Pt->GetValue(trk);
			float rc_phi = leaf_rc_Phi->GetValue(trk);
			float rc_charge = leaf_rc_Charge->GetValue(trk);
                        float eta = leaf_rc_Eta->GetValue(trk);
			if(eta>1 || eta<-1) continue; 
                        int nFitHits = leaf_rc_Nfithits->GetValue(trk);
                        int nMaxHits = leaf_rc_Nmaxhits->GetValue(trk);
                        float ratio = float(nFitHits)/nMaxHits;
                        int ndEdxHits = leaf_rc_Ndedxhits->GetValue(trk);
                        if(nFitHits<20 || nFitHits>50) continue;
                        if(ratio<0.52 || ratio>1.05) continue;
			int rc_id = leaf_rc_id->GetValue(trk);
			int bin = h_pt->FindBin(rc_pt);
			float weight = h_pt->GetBinContent(bin)/bin_mean;

			for(int tr = 0; tr < NmcTracks; tr++) {
                        	float mc_px = leaf_mc_Px->GetValue(tr);
                        	float mc_py = leaf_mc_Py->GetValue(tr);
				float mc_phi = atan2(mc_py,mc_px);
				int mc_id = leaf_mc_id->GetValue(tr);
				if(mc_id==rc_id) {
					Mpt[Mnmatch] = rc_pt;
					Mdca[Mnmatch] = dca;
					MndEdxHits[Mnmatch] = ndEdxHits;
					float dphi = rc_phi - mc_phi;
					if(dphi > PI) dphi -= 2*PI;
					if(dphi <-PI) dphi += 2*PI;
					if(rc_pt>0.15 && rc_pt<2) Hist_dphi->Fill(dphi,weight);
					if(rc_pt>0.15 && rc_pt<0.5) Hist_dphi1->Fill(dphi,weight);
					else if(rc_pt<1) Hist_dphi2->Fill(dphi,weight);
                                        else if(rc_pt<1.5) Hist_dphi3->Fill(dphi,weight);
                                        else if(rc_pt<2) Hist_dphi4->Fill(dphi,weight);
					Mdphi[Mnmatch] = dphi;
					Mweight[Mnmatch] = weight;
					Mnmatch++;
                                        if(rc_charge>0) Np--;
                                        else Nn--;
				}
			}
		}
              float asym = (float(Np)-float(Nn))/(Np+Nn);
              Hist_asym->Fill(asym);
                if(asym<mean-rms) Hist_asym_L->Fill(asym);
                else if(asym<mean-0.3*rms) Hist_asym_LM->Fill(asym);
                else if(asym<mean+0.3*rms) Hist_asym_M->Fill(asym);
                else if(asym<mean+rms) Hist_asym_RM->Fill(asym);
                else Hist_asym_R->Fill(asym);
//cout<<"Nmatch = "<<Mnmatch<<endl;
		for(int ii=0;ii<Mnmatch;ii++) {
                                float dphi = Mdphi[ii];
                        if(Mdca[ii]>1) continue;
                        if(MndEdxHits[ii]<=10) continue;
                        p_dphi_pt->Fill(Mpt[ii],100*cos(2*dphi));
                        if(Mpt[ii]<0.15 || Mpt[ii]>0.5) continue;
                        Hist_dphi_pion->Fill(Mdphi[ii]);
                                p_dphi->Fill(1,100*cos(2*dphi));
                                p_dphi->Fill(2,100*sin(2*dphi));
                                if(asym<mean-rms) {p_dphi_L->Fill(3,100*cos(2*dphi));
                                                        p_dphi_L->Fill(1,100*cos(2*dphi), Mweight[ii]);
                                                        p_dphi_L->Fill(4,100*sin(2*dphi));
                                                        p_dphi_L->Fill(2,100*sin(2*dphi), Mweight[ii]);}
                                else if(asym<mean-0.3*rms) {p_dphi_LM->Fill(3,100*cos(2*dphi));
                                                        p_dphi_LM->Fill(1,100*cos(2*dphi), Mweight[ii]);
                                                        p_dphi_LM->Fill(4,100*sin(2*dphi));
                                                        p_dphi_LM->Fill(2,100*sin(2*dphi), Mweight[ii]);}
                                else if(asym<mean+0.3*rms) {p_dphi_M->Fill(3,100*cos(2*dphi));
                                                        p_dphi_M->Fill(1,100*cos(2*dphi), Mweight[ii]);
                                                        p_dphi_M->Fill(4,100*sin(2*dphi));
                                                        p_dphi_M->Fill(2,100*sin(2*dphi), Mweight[ii]);}
                                else if(asym<mean+rms) {p_dphi_RM->Fill(3,100*cos(2*dphi));
                                                        p_dphi_RM->Fill(1,100*cos(2*dphi), Mweight[ii]);
                                                        p_dphi_RM->Fill(4,100*sin(2*dphi));
                                                        p_dphi_RM->Fill(2,100*sin(2*dphi), Mweight[ii]);}
                                else {p_dphi_R->Fill(3,100*cos(2*dphi));
                                                        p_dphi_R->Fill(1,100*cos(2*dphi), Mweight[ii]);
                                                        p_dphi_R->Fill(4,100*sin(2*dphi));
                                                        p_dphi_R->Fill(2,100*sin(2*dphi), Mweight[ii]);}

			for(int jj=0;jj<Mnmatch;jj++) {
				if(ii==jj) continue;
				if(Mpt[jj]<0.15 || Mpt[jj]>2.0) continue;
				float ddphi = Mdphi[ii] - Mdphi[jj];
				Hist_ddphi->Fill(ddphi);
				p_ddphi->Fill(1,100*cos(2*ddphi));
				p_ddphi->Fill(2,100*sin(2*ddphi));
                                if(asym<mean-rms) {p_ddphi_L->Fill(3,100*cos(2*ddphi));
							p_ddphi_L->Fill(1,100*cos(2*ddphi), Mweight[ii]*Mweight[jj]);
							p_ddphi_L->Fill(4,100*sin(2*ddphi));
							p_ddphi_L->Fill(2,100*sin(2*ddphi), Mweight[ii]*Mweight[jj]);}
                                else if(asym<mean-0.3*rms) {p_ddphi_LM->Fill(3,100*cos(2*ddphi));
							p_ddphi_LM->Fill(1,100*cos(2*ddphi), Mweight[ii]*Mweight[jj]);
                                                        p_ddphi_LM->Fill(4,100*sin(2*ddphi));			
							p_ddphi_LM->Fill(2,100*sin(2*ddphi), Mweight[ii]*Mweight[jj]);}
                                else if(asym<mean+0.3*rms) {p_ddphi_M->Fill(3,100*cos(2*ddphi));
							p_ddphi_M->Fill(1,100*cos(2*ddphi), Mweight[ii]*Mweight[jj]);
                                                        p_ddphi_M->Fill(4,100*sin(2*ddphi));
							p_ddphi_M->Fill(2,100*sin(2*ddphi), Mweight[ii]*Mweight[jj]);}
                                else if(asym<mean+rms) {p_ddphi_RM->Fill(3,100*cos(2*ddphi));
							p_ddphi_RM->Fill(1,100*cos(2*ddphi), Mweight[ii]*Mweight[jj]);
                                                        p_ddphi_RM->Fill(4,100*sin(2*ddphi));
							p_ddphi_RM->Fill(2,100*sin(2*ddphi), Mweight[ii]*Mweight[jj]);}
                                else {p_ddphi_R->Fill(3,100*cos(2*ddphi));
							p_ddphi_R->Fill(1,100*cos(2*ddphi), Mweight[ii]*Mweight[jj]);
                                                        p_ddphi_R->Fill(4,100*sin(2*ddphi));
							p_ddphi_R->Fill(2,100*sin(2*ddphi), Mweight[ii]*Mweight[jj]);}
			}
		}
	}

        fout.Write();
	return;

}
