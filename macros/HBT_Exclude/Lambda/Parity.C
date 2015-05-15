/*
 * EId.C
 * macro to identify electrons from *.Tracks.root ntuple files
 * kurnadi, 2005/11/14
 */
using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TSystem.h"
#include "MyRef.h"
//#include "StThreeVectorD.hh"
//#include "StPhysicalHelixD.hh"
//#include "StLorentzVectorD.hh"

//const float PI = TMath::Pi();
const float PI = 3.14159265;
//cout<<"PI = "<<PI<<endl;
const float MM = 2/PI;
const int   opt_baryon = 0;
const float pt_trig_up = 2;
const float pt_trig_lo = .15;
const float pt_asso_up = 2.0;//1;
const float pt_asso_lo = 0.15;//0.15;
const float EtaCut = 1;
const float DcaCut = 2;
const Float_t Vz_offset = 0;
const Float_t Vz_cut = 40;    //40 for 39 Gev, 75 for 11 GeV and 70 for 7GeV
const int cenDef[9] = {7, 15, 28, 50, 81, 125, 185, 265, 316};	// 39 GeV
const int bad_Ref_day2[12] = {13703,13800,13802,14005,14107,14301,14307,14600,15901,16304,16407,16705};
const int Nrun_MB1 = 134;
const int Nrun_MB2 = 60;
const int Nrun_MB5 = 107;
const int Nrun_MB6 = 42;
const int bad_Ref_day3_MB1[Nrun_MB1] = {133010,133011,133022,133027,133028,134065,127019,128021,133020,137004,127039,132024,132025,132026,132032,127020,127048,127049,128028,128029,128030,128031,128032,132052,133041,136039,127002,132061,134017,134063,135004,135039,137020,127030,128007,132062,133005,133039,133040,133053,133054,134006,134007,134018,134026,134028,134038,134041,135019,135020,135045,135057,136007,136032,136064,137010,127009,127046,128038,132019,132022,132044,132045,133002,133019,133038,133052,134005,134040,134055,135002,135033,135048,135049,136069,137003,127021,127022,127023,127024,128024,128025,132020,132021,132023,132033,132048,132051,132057,132063,132065,133021,134008,135012,135021,135024,135030,135054,136044,136081,136086,127003,127010,127011,127017,127018,127032,132009,132034,132043,132066,132069,133018,134023,134057,136005,136006,136014,136017,136022,136023,136024,136025,136027,136028,136029,136030,136031,136034,136053,136054,136070,136071,138017};//MB1
const int bad_Ref_day3_MB2[Nrun_MB2] = {139032,139043,139044,139045,142002,139042,140021,140029,142063,142064,142065,144004,138081,138082,138087,138088,138089,138090,138091,139002,139003,139006,139007,139008,139009,139010,139015,139016,139017,139018,139021,142016,142033,142061,142062,144051,138092,139019,139020,140016,140020,141003,141004,141026,141062,141065,142001,142013,142023,142034,142046,142068,142076,143009,143024,143058,144016,144028,144033,145003};  //MB2
const int bad_Ref_day3_MB5[Nrun_MB5] = {155050,155056,158010,165028,154043,154044,154045,155058,158069,158070,158072,158073,164067,154046,154047,155008,155009,156015,156062,156063,158074,162015,154067,155002,155012,155047,156008,156009,157023,157030,157052,158006,159023,160021,161006,161015,161060,162004,162028,162034,163024,163058,164009,164056,164066,164089,165013,154048,154066,155011,155021,155038,155051,155060,155062,155064,156004,156035,156056,157012,157014,157031,157038,157051,158015,158021,158026,158040,158041,158051,158054,158056,158057,158058,158061,159005,159021,159022,159024,160016,160025,161007,161014,161017,161020,161022,161053,162017,162030,162035,162055,162056,162057,162058,163006,163008,163015,164011,164037,164043,164086,165001,165003,165005,165007,165026,165031}; //MB5
const int bad_Ref_day3_MB6[Nrun_MB6] = {166051,170016,167014,170034,170050,170051,171009,171015,167049,168010,169028,169032,170007,165042,166052,166059,167002,167040,167048,169031,169059,170009,170012,170018,170020,170031,171004,171014,166002,166003,167015,167024,168009,168022,168060,168077,169033,169034,170044,170045,170054,170056};//MB6
const Float_t MeanNetProton[9] = {0.24, 0.48, 0.90, 1.56, 2.49, 3.76, 5.38, 6.89, 8.13};

void Parity(int cen, int correction){	//main_function
 ifstream infile_la("../eff_fit_par_Lambda.dat");
    ifstream infile_tpc_p("../proton_tpc_eff_fit.dat");
    ifstream infile_tof_p("../proton_tof_eff_fit.dat");

    float par0_la, par1_la, par2_la;
    float par0_p, par1_p, par2_p;
    float par0_p_tof, par1_p_tof, par2_p_tof, par3_p_tof;
    for(int i = 0; i < 9; i++){
        if(i == (cen-1)){
            infile_la >> par0_la >> par1_la >> par2_la;
            infile_tpc_p >> par0_p >> par1_p >> par2_p;
            infile_tof_p >> par0_p_tof >> par1_p_tof >> par2_p_tof >> par3_p_tof;
        }
        else continue;
    }
    TF1* eff_la = new TF1("eff_la",  "[0]*TMath::Erf(x-[1])+[2]", 0, 3.0);
    TF1* eff_tpc_p = new TF1("eff_tpc_p", "[0]*exp(-pow([1]/x, [2]))", 0, 3.0);
    TF1* eff_tof_p = new TF1("eff_tof_p", "[0]*TMath::Erf([1]*(x-[2]))+[3]", 0, 3.0);

    eff_la -> SetParameter(0, par0_la);
    eff_la -> SetParameter(1, par1_la);
    eff_la -> SetParameter(2, par2_la);

    eff_tpc_p -> SetParameter(0, par0_p);
    eff_tpc_p -> SetParameter(1, par1_p);
    eff_tpc_p -> SetParameter(2, par2_p);

    eff_tof_p -> SetParameter(0, par0_p_tof);
    eff_tof_p -> SetParameter(1, par1_p_tof);
    eff_tof_p -> SetParameter(2, par2_p_tof);
    eff_tof_p -> SetParameter(3, par3_p_tof);

delete gRandom;
gRandom = new TRandom3(0);

struct StRefMultCorr refmultCorrUtil  = StRefMultCorr("refmult") ;

float Eweight = 1;
TFile *fMult = new TFile("MultWeight.root","READ");
TH2D *Mult_Weight;
if(fMult->IsOpen()) Mult_Weight = (TH2D*)fMult->Get("gRefMult_Weight");
else cout<<"no mult weight files!"<<endl;
const int Phibin = 80;
char fname[200];
char fname_pr[200];
float PsiShiftE1=0,PsiShiftE2=0,PsiShiftE3=0,PsiShiftE4=0,PsiShiftE5=0,PsiShiftE6=0,PsiShiftE7=0,PsiShiftE8=0;
float PsiShiftW1=0,PsiShiftW2=0,PsiShiftW3=0,PsiShiftW4=0,PsiShiftW5=0,PsiShiftW6=0,PsiShiftW7=0,PsiShiftW8=0;
float PsiShiftF1=0,PsiShiftF2=0,PsiShiftF3=0,PsiShiftF4=0,PsiShiftF5=0,PsiShiftF6=0,PsiShiftF7=0,PsiShiftF8=0;
const int order = 4;
float PhiMean[4*order]={0,};
float PhiMean_Lambda[4*order]={0,};
float PhiMean_Proton[4*order]={0,};
float PhiWgtFF[Phibin][4],PhiWgtRF[Phibin][4];
float PhiWgtFF_Lambda[Phibin][4],PhiWgtRF_Lambda[Phibin][4];
float PhiWgtFF_Proton[Phibin][4],PhiWgtRF_Proton[Phibin][4];
TProfile2D *TPCmean_FF, *TPCmean_RF, *TPCmean_FF_Lambda, *TPCmean_RF_Lambda, *TPCmean_FF_Proton, *TPCmean_RF_Proton;
TH2D *TPCPhi_FF, *TPCPhi_RF;
TH2D *TPCPhi_FF_Lambda, *TPCPhi_RF_Lambda;
TH2D *TPCPhi_FF_Proton, *TPCPhi_RF_Proton;
TProfile2D *Read_TPC_EP_full, *Read_TPC_EP_east, *Read_TPC_EP_west;
sprintf(fname,"cen%d.weight.root",cen);
TFile *fWgt=new TFile(fname,"READ");
if(!fWgt->IsOpen()) cout<<"no phi weight files!"<<endl;
if(fWgt->IsOpen()) {
        TPCmean_FF = (TProfile2D*)fWgt->Get("TPCmeanPhi_FF");
        TPCmean_RF = (TProfile2D*)fWgt->Get("TPCmeanPhi_RF");
	TPCmean_FF_Lambda = (TProfile2D*)fWgt->Get("TPCmeanPhi_FF_Lambda");
	TPCmean_RF_Lambda = (TProfile2D*)fWgt->Get("TPCmeanPhi_RF_Lambda");
//        TPCmean_FF_Proton = (TProfile2D*)fWgt->Get("TPCmeanPhi_FF_Proton");
//        TPCmean_RF_Proton = (TProfile2D*)fWgt->Get("TPCmeanPhi_RF_Proton");
	TPCPhi_FF = (TH2D*)fWgt->Get("Hist_Phi_FF");
	TPCPhi_RF = (TH2D*)fWgt->Get("Hist_Phi_RF");
	TPCPhi_FF_Lambda = (TH2D*)fWgt->Get("Hist_Phi_FF_Lambda");
	TPCPhi_RF_Lambda = (TH2D*)fWgt->Get("Hist_Phi_RF_Lambda");
//        TPCPhi_FF_Proton = (TH2D*)fWgt->Get("Hist_Phi_FF_Proton");
//        TPCPhi_RF_Proton = (TH2D*)fWgt->Get("Hist_Phi_RF_Proton");
	Read_TPC_EP_full = (TProfile2D*)fWgt->Get("pTPC_EP_full");
	Read_TPC_EP_east = (TProfile2D*)fWgt->Get("pTPC_EP_east");
	Read_TPC_EP_west = (TProfile2D*)fWgt->Get("pTPC_EP_west");
	for(int j=0;j<4;j++) {
                float PhiMeanFF = TPCPhi_FF->ProjectionX("",j+1,j+1)->GetSum()/(float)Phibin;
                float PhiMeanRF = TPCPhi_RF->ProjectionX("",j+1,j+1)->GetSum()/(float)Phibin;
		float PhiMeanFF_Lambda = TPCPhi_FF_Lambda->ProjectionX("",j+1,j+1)->GetSum()/(float)Phibin;
		float PhiMeanRF_Lambda = TPCPhi_RF_Lambda->ProjectionX("",j+1,j+1)->GetSum()/(float)Phibin;
//                float PhiMeanFF_Proton = TPCPhi_FF_Proton->ProjectionX("",j+1,j+1)->GetSum()/(float)Phibin;
//                float PhiMeanRF_Proton = TPCPhi_RF_Proton->ProjectionX("",j+1,j+1)->GetSum()/(float)Phibin;
		for(int i=0;i<Phibin;i++) {
                	PhiWgtFF[i][j] = (TPCPhi_FF->GetBinContent(i+1,j+1)>0)? PhiMeanFF/(TPCPhi_FF->GetBinContent(i+1,j+1)):1;
                        PhiWgtRF[i][j] = (TPCPhi_RF->GetBinContent(i+1,j+1)>0)? PhiMeanRF/(TPCPhi_RF->GetBinContent(i+1,j+1)):1;
			PhiWgtFF_Lambda[i][j] = (TPCPhi_FF_Lambda->GetBinContent(i+1,j+1)>0)? PhiMeanFF_Lambda/(TPCPhi_FF_Lambda->GetBinContent(i+1,j+1)):1;
			PhiWgtRF_Lambda[i][j] = (TPCPhi_RF_Lambda->GetBinContent(i+1,j+1)>0)? PhiMeanRF_Lambda/(TPCPhi_RF_Lambda->GetBinContent(i+1,j+1)):1;
//                        PhiWgtFF_Proton[i][j] = (TPCPhi_FF_Proton->GetBinContent(i+1,j+1)>0)? PhiMeanFF_Proton/(TPCPhi_FF_Proton->GetBinContent(i+1,j+1)):1;
//                        PhiWgtRF_Proton[i][j] = (TPCPhi_RF_Proton->GetBinContent(i+1,j+1)>0)? PhiMeanRF_Proton/(TPCPhi_RF_Proton->GetBinContent(i+1,j+1)):1;
//                	cout<<" PhiWgt= "<<PhiWgtFF[i][j];
		}
	}
}  
else {for(int j=0;j<4;j++) for(int i=0;i<Phibin;i++) { PhiWgtFF[i][j] = 1; PhiWgtRF[i][j] = 1; PhiWgtFF_Lambda[i][j] = 1; PhiWgtRF_Lambda[i][j] = 1;}
}

sprintf(fname_pr,"cen%d.weight_proton.root",cen);
TFile *fWgt_pr=new TFile(fname_pr,"READ");
if(!fWgt_pr->IsOpen()) cout<<"no proton phi weight files!"<<endl;
if(fWgt_pr->IsOpen()) {
        TPCmean_FF_Proton = (TProfile2D*)fWgt_pr->Get("TPCmeanPhi_FF_Proton");
        TPCmean_RF_Proton = (TProfile2D*)fWgt_pr->Get("TPCmeanPhi_RF_Proton");
        TPCPhi_FF_Proton = (TH2D*)fWgt_pr->Get("Hist_Phi_FF_Proton");
        TPCPhi_RF_Proton = (TH2D*)fWgt_pr->Get("Hist_Phi_RF_Proton");
        TPCPhi_FF_Proton = (TH2D*)fWgt_pr->Get("Hist_Phi_FF_Proton");
        TPCPhi_RF_Proton = (TH2D*)fWgt_pr->Get("Hist_Phi_RF_Proton");
        for(int j=0;j<4;j++) {
                float PhiMeanFF_Proton = TPCPhi_FF_Proton->ProjectionX("",j+1,j+1)->GetSum()/(float)Phibin;
                float PhiMeanRF_Proton = TPCPhi_RF_Proton->ProjectionX("",j+1,j+1)->GetSum()/(float)Phibin;
                for(int i=0;i<Phibin;i++) {
                        PhiWgtFF_Proton[i][j] = (TPCPhi_FF_Proton->GetBinContent(i+1,j+1)>0)? PhiMeanFF_Proton/(TPCPhi_FF_Proton->GetBinContent(i+1,j+1)):1;
                        PhiWgtRF_Proton[i][j] = (TPCPhi_RF_Proton->GetBinContent(i+1,j+1)>0)? PhiMeanRF_Proton/(TPCPhi_RF_Proton->GetBinContent(i+1,j+1)):1;
		}
        }
}
else {for(int j=0;j<4;j++) for(int i=0;i<Phibin;i++) { PhiWgtRF_Proton[i][j] = 1;PhiWgtFF_Proton[i][j] = 1;}
}


	TChain* chain = new TChain("StrangenessDst");
	int nfile = 0;

        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data154/*.lambda.picodst.root");
/*
        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data155/*.lambda.picodst.root");
        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data156/*.lambda.picodst.root");
        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data157/*.lambda.picodst.root");
        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data158/*.lambda.picodst.root");
        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data159/*.lambda.picodst.root");

        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data160/*.lambda.picodst.root");
        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data161/*.lambda.picodst.root");
        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data162/*.lambda.picodst.root");
        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data163/*.lambda.picodst.root");
        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data164/*.lambda.picodst.root");
        nfile += chain->Add("/media/Disk_Yan/liwen/Run11_200GeV/Data165/*.lambda.picodst.root");
*/
		cout <<"Added "<<nfile<<" files"<<endl;
		cout<<"# entries in chain: "<<chain->GetEntries()<<endl;

char fname_out[200];
sprintf(fname_out,"cen%d.ntuple_result_MSC_39_subEP.154.HBT.root",cen);
	TFile fout(fname_out,"RECREATE");
	//defining variables
	Float_t PVtxz, Bz, psi_E, psi_W, mod_E, mod_W;						//run, event info
	Int_t   Run, Day, Day2, Day3, Event, Trigger, RefMult, Centrality, NPTracks, zdcCoincidenceRate;	//
	Float_t Charge, ndEdx, nSigma_p, nSigma_pi, DCAGlobal, Eta, Theta, Phi, Pt;		//track info	

	//defining histograms
	TH1D* hEventTally = new TH1D("EventTally","Event Tally",10,0,1);
	hEventTally->SetBit(TH1::kCanRebin);
	hEventTally->SetStats(0);

	TH1D *hBz      = new TH1D("hBz","magnetic field",10, -10, 10);
	TH1D *hTrigger = new TH1D("hTrigger","hTrigger",200, 0.5, 200.5);
	TH1D *hCentrality = new TH1D("hCentrality","hCentrality",10,0,10);
	TH1D *hVertexZ = new TH1D("hVertexZ","hVertexZ",100,-100,100);
	TH2D *hMult_Vz = new TH2D("hMult_Vz","hMult_Vz",1000,-0.5,999.5,100,-100,100);
        TH2D *hMult_Vz_new = new TH2D("hMult_Vz_new","hMult_Vz_new",1000,-0.5,999.5,100,-100,100);

	TH2D* SigmaProtonPlus = new TH2D("SigmaProtonPlus","SigmaProtonPlus",500,0,2.5,1000,-10,10);
	TH2D* SigmaProtonMinus = new TH2D("SigmaProtonMinus","SigmaProtonMinus",500,0,2.5,1000,-10,10);

        TProfile2D *pTPCmeanPhi_FF = new TProfile2D("TPCmeanPhi_FF","TPCmeanPhi_FF",
                                        8*order,0.5,8*order+0.5,10000,8000,18000,-1,1,"");
        TProfile2D *pTPCmeanPhi_FF_Lambda = new TProfile2D("TPCmeanPhi_FF_Lambda","TPCmeanPhi_FF_Lambda",
                                        8*order,0.5,8*order+0.5,10000,8000,18000,-1,1,"");
        TProfile2D *pTPCmeanPhi_FF_Proton = new TProfile2D("TPCmeanPhi_FF_Proton","TPCmeanPhi_FF_Proton",
                                        8*order,0.5,8*order+0.5,10000,8000,18000,-1,1,"");
        TProfile2D *pTPCmeanPhi_RF = new TProfile2D("TPCmeanPhi_RF","TPCmeanPhi_RF",
                                        8*order,0.5,8*order+0.5,10000,8000,18000,-1,1,"");
        TProfile2D *pTPCmeanPhi_RF_Lambda = new TProfile2D("TPCmeanPhi_RF_Lambda","TPCmeanPhi_RF_Lambda",
                                        8*order,0.5,8*order+0.5,10000,8000,18000,-1,1,"");
        TProfile2D *pTPCmeanPhi_RF_Proton = new TProfile2D("TPCmeanPhi_RF_Proton","TPCmeanPhi_RF_Proton",
                                        8*order,0.5,8*order+0.5,10000,8000,18000,-1,1,"");

        TH1D* Hist_proton = new TH1D("Hist_proton","Hist_proton",50,-0.5,49.5);
        TH1D* Hist_pbar = new TH1D("Hist_pbar","Hist_pbar",50,-0.5,49.5);
	TH1D* Hist_netP = new TH1D("Hist_netP","Hist_netP",99,-49.5,49.5);
        TH2D *hEtaPtDist = new TH2D("EtaPtDist","EtaPtDist",26, -1.3, 1.3,300,0,15);
	TH2D *hEtaPt_Proton_Dist = new TH2D("EtaPt_Proton_Dist","EtaPt_Proton_Dist",26, -1.3, 1.3,300,0,15);
	TH2D* Hist_Phi = new TH2D("Hist_Phi","Hist_Phi",Phibin,-PI,PI,4,0.5,4.5);
	TH2D* Hist_TPC_EP_east = new TH2D("Hist_TPC_EP_east","Hist_TPC_EP_east",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_west = new TH2D("Hist_TPC_EP_west","Hist_TPC_EP_west",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_full = new TH2D("Hist_TPC_EP_full","Hist_TPC_EP_full",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_east_flat = new TH2D("Hist_TPC_EP_east_flat","Hist_TPC_EP_east_flat",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_west_flat = new TH2D("Hist_TPC_EP_west_flat","Hist_TPC_EP_west_flat",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_full_flat = new TH2D("Hist_TPC_EP_full_flat","Hist_TPC_EP_full_flat",36,0,PI,100,80,180);
	TProfile2D* pTPC_EP_east = new TProfile2D("pTPC_EP_east","pTPC_EP_east",8,0.5,8.5,100,80,180,-1,1,"");
        TProfile2D* pTPC_EP_west = new TProfile2D("pTPC_EP_west","pTPC_EP_west",8,0.5,8.5,100,80,180,-1,1,"");
        TProfile2D* pTPC_EP_full = new TProfile2D("pTPC_EP_full","pTPC_EP_full",8,0.5,8.5,100,80,180,-1,1,"");
	TH1F* Hist_dif_count = new TH1F("Hist_dif_count","Hist_dif_count",500,-250,250);
        TH1F* Hist_ful_count = new TH1F("Hist_ful_count","Hist_ful_count",1000,0,1000);


	TH1D* Hist_Pt = new TH1D("Hist_Pt","Hist_Pt",300,0,15);
        TH2D* Hist_Phi_FF = new TH2D("Hist_Phi_FF","Hist_Phi_FF",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_FF_Lambda = new TH2D("Hist_Phi_FF_Lambda","Hist_Phi_FF_Lambda",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_FF_Proton = new TH2D("Hist_Phi_FF_Proton","Hist_Phi_FF_Proton",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_RF = new TH2D("Hist_Phi_RF","Hist_Phi_RF",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_RF_Lambda = new TH2D("Hist_Phi_RF_Lambda","Hist_Phi_RF_Lambda",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_RF_Proton = new TH2D("Hist_Phi_RF_Proton","Hist_Phi_RF_Proton",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_FF_new = new TH2D("Hist_Phi_FF_new","Hist_Phi_FF_new",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_RF_new = new TH2D("Hist_Phi_RF_new","Hist_Phi_RF_new",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_FF_Lambda_new = new TH2D("Hist_Phi_FF_Lambda_new","Hist_Phi_FF_Lambda_new",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_RF_Lambda_new = new TH2D("Hist_Phi_RF_Lambda_new","Hist_Phi_RF_Lambda_new",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_FF_Proton_new = new TH2D("Hist_Phi_FF_Proton_new","Hist_Phi_FF_Proton_new",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_RF_Proton_new = new TH2D("Hist_Phi_RF_Proton_new","Hist_Phi_RF_Proton_new",Phibin,-PI,PI,4,0.5,4.5);
	TProfile *Hist_cos = new TProfile("Hist_cos","Hist_cos",3,0.5,3.5,-1,1,"");
       	TH1D* hDpt   = new TH1D("hDpt","hDpt",200,0,2);
	TH1D* hQinv  = new TH1D("hQinv","hQinv",1000,0,10);
	TH1D* hQinv2 = new TH1D("hQinv2","hQinv2",1100,-1,10);

   Float_t pdgmass = 1.115683;
   Float_t masswidth = 0.07+0.01;
   TH1F * hmInvMass = new TH1F("hmInvMass","Invariant mass", 200,pdgmass-masswidth,pdgmass+masswidth);

	TProfile *Hist_v2_pt_obs = new TProfile("Hist_v2_pt_obs","Hist_v2_pt_obs",300,0,15,-100,100,"");
        TProfile *Hist_v2_pt_pos_obs = new TProfile("Hist_v2_pt_pos_obs","Hist_v2_pt_pos_obs",300,0,15,-100,100,"");
        TProfile *Hist_v2_pt_neg_obs = new TProfile("Hist_v2_pt_neg_obs","Hist_v2_pt_neg_obs",300,0,15,-100,100,"");

	TH1D* Hist_dQ_in1  = new TH1D("Hist_dQ_in1","Hist_dQ_in1",201,-100.5,100.5);
	TH1D* Hist_dQ_out1 = new TH1D("Hist_dQ_out1","Hist_dQ_out1",201,-100.5,100.5);
        TH1D* Hist_dQ_in2  = new TH1D("Hist_dQ_in2","Hist_dQ_in2",201,-100.5,100.5);
        TH1D* Hist_dQ_out2 = new TH1D("Hist_dQ_out2","Hist_dQ_out2",201,-100.5,100.5);

	TProfile *pParity_int_MSC_same_in1 = new TProfile("pParity_int_MSC_same_in1","pParity_int_MSC_same_in1",201,-100.5,100.5,-100,100,"");
        TProfile *pParity_int_MSC_oppo_in1 = new TProfile("pParity_int_MSC_oppo_in1","pParity_int_MSC_oppo_in1",201,-100.5,100.5,-100,100,"");
        TProfile *pParity_int_MSC_same_out1 = new TProfile("pParity_int_MSC_same_out1","pParity_int_MSC_same_out1",201,-100.5,100.5,-100,100,"");
        TProfile *pParity_int_MSC_oppo_out1 = new TProfile("pParity_int_MSC_oppo_out1","pParity_int_MSC_oppo_out1",201,-100.5,100.5,-100,100,"");
        TProfile *pParity_int_MSC_same_in2 = new TProfile("pParity_int_MSC_same_in2","pParity_int_MSC_same_in2",201,-100.5,100.5,-100,100,"");
        TProfile *pParity_int_MSC_oppo_in2 = new TProfile("pParity_int_MSC_oppo_in2","pParity_int_MSC_oppo_in2",201,-100.5,100.5,-100,100,"");
        TProfile *pParity_int_MSC_same_out2 = new TProfile("pParity_int_MSC_same_out2","pParity_int_MSC_same_out2",201,-100.5,100.5,-100,100,"");
        TProfile *pParity_int_MSC_oppo_out2 = new TProfile("pParity_int_MSC_oppo_out2","pParity_int_MSC_oppo_out2",201,-100.5,100.5,-100,100,"");

        TProfile *pParity_int_MSC1_obs1 = new TProfile("pParity_int_MSC1_obs1","pParity_int_MSC1_obs1",4,0.5,4.5,-200,200,"");
        TProfile *pParity_int_MSC1_obs2 = new TProfile("pParity_int_MSC1_obs2","pParity_int_MSC1_obs2",4,0.5,4.5,-200,200,"");
        TProfile *pParity_int_MSC2_obs1 = new TProfile("pParity_int_MSC2_obs1","pParity_int_MSC2_obs1",4,0.5,4.5,-200,200,"");
        TProfile *pParity_int_MSC2_obs2 = new TProfile("pParity_int_MSC2_obs2","pParity_int_MSC2_obs2",4,0.5,4.5,-200,200,"");

//        TProfile *pParity_int_MSC_dQ_in_obs1 = new TProfile("pParity_int_MSC_dQ_in_obs1","pParity_int_MSC_dQ_in_obs1",4,0.5,4.5,-200,200,"");
//        TProfile *pParity_int_MSC_dQ_in_obs2 = new TProfile("pParity_int_MSC_dQ_in_obs2","pParity_int_MSC_dQ_in_obs2",4,0.5,4.5,-200,200,"");
//        TProfile *pParity_int_MSC_dQ_out_obs1 = new TProfile("pParity_int_MSC_dQ_out_obs1","pParity_int_MSC_dQ_out_obs1",4,0.5,4.5,-200,200,"");
//        TProfile *pParity_int_MSC_dQ_out_obs2 = new TProfile("pParity_int_MSC_dQ_out_obs2","pParity_int_MSC_dQ_out_obs2",4,0.5,4.5,-200,200,"");
	TProfile *pParity_int_obs1 = new TProfile("Parity_int_obs1","Parity_int_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_obs2 = new TProfile("Parity_int_obs2","Parity_int_obs2",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_r_obs1 = new TProfile("Parity_int_r_obs1","Parity_int_r_obs1",4,0.5,4.5,-1600,1600,"");
        TProfile *pParity_int_r_obs2 = new TProfile("Parity_int_r_obs2","Parity_int_r_obs2",4,0.5,4.5,-1600,1600,"");
        TProfile *pParity_int_rr_obs1 = new TProfile("Parity_int_rr_obs1","Parity_int_rr_obs1",4,0.5,4.5,-1600,1600,"");
        TProfile *pParity_int_rr_obs2 = new TProfile("Parity_int_rr_obs2","Parity_int_rr_obs2",4,0.5,4.5,-1600,1600,"");
        TProfile *pParity_int_s_obs1 = new TProfile("Parity_int_s_obs1","Parity_int_s_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_s_obs2 = new TProfile("Parity_int_s_obs2","Parity_int_s_obs2",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ss_obs1 = new TProfile("Parity_int_ss_obs1","Parity_int_ss_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ss_obs2 = new TProfile("Parity_int_ss_obs2","Parity_int_ss_obs2",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ss_obs2_eff = new TProfile("Parity_int_ss_obs2_eff","Parity_int_ss_obs2_eff",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_w_obs1 = new TProfile("Parity_int_w_obs1","Parity_int_w_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_w_obs2 = new TProfile("Parity_int_w_obs2","Parity_int_w_obs2",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ww_obs1 = new TProfile("Parity_int_ww_obs1","Parity_int_ww_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ww_obs2 = new TProfile("Parity_int_ww_obs2","Parity_int_ww_obs2",4,0.5,4.5,-100,100,"");

        TProfile *pParity_int_same_run = new TProfile("Parity_int_same_run","Parity_int_same_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_rr_same_run = new TProfile("Parity_int_rr_same_run","Parity_int_rr_same_run",10000,8000,18000,-100,100,"");
	TProfile *pParity_int_ss_same_run = new TProfile("Parity_int_ss_same_run","Parity_int_ss_same_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_ww_same_run = new TProfile("Parity_int_ww_same_run","Parity_int_ww_same_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_oppo_run = new TProfile("Parity_int_oppo_run","Parity_int_oppo_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_rr_oppo_run = new TProfile("Parity_int_rr_oppo_run","Parity_int_rr_oppo_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_ss_oppo_run = new TProfile("Parity_int_ss_oppo_run","Parity_int_ss_oppo_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_ww_oppo_run = new TProfile("Parity_int_ww_oppo_run","Parity_int_ww_oppo_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_ss_ran_obs1 = new TProfile("Parity_int_ss_ran_obs1","Parity_int_ss_ran_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ss_ran_obs2 = new TProfile("Parity_int_ss_ran_obs2","Parity_int_ss_ran_obs2",4,0.5,4.5,-100,100,"");
        TProfile2D *pParity_eta_ss_obs1 = new TProfile2D("Parity_eta_ss_obs1","Parity_eta_ss_obs1",12,0.5,12.5,20,-1,1,-100,100,"");
        TProfile2D *pParity_eta_ss_obs2 = new TProfile2D("Parity_eta_ss_obs2","Parity_eta_ss_obs2",12,0.5,12.5,20,-1,1,-100,100,"");
        TProfile2D *pParity_Deta_ss_obs1 = new TProfile2D("Parity_Deta_ss_obs1","Parity_Deta_ss_obs1",12,0.5,12.5,20,0,2,-100,100,"");
        TProfile2D *pParity_Deta_ss_obs2 = new TProfile2D("Parity_Deta_ss_obs2","Parity_Deta_ss_obs2",12,0.5,12.5,20,0,2,-100,100,"");
        TProfile2D *pParity_pt_ss_obs1  = new TProfile2D("Parity_pt_ss_obs1","Parity_pt_ss_obs1",12,0.5,12.5,20,0,2.0,-100,100,"");
        TProfile2D *pParity_pt_ss_obs2  = new TProfile2D("Parity_pt_ss_obs2","Parity_pt_ss_obs2",12,0.5,12.5,20,0,2.0,-100,100,"");
        TProfile2D *pParity_Dpt_ss_obs1 = new TProfile2D("Parity_Dpt_ss_obs1","Parity_Dpt_ss_obs1",12,0.5,12.5,200,0,2.0,-100,100,"");
        TProfile2D *pParity_Dpt_ss_obs2 = new TProfile2D("Parity_Dpt_ss_obs2","Parity_Dpt_ss_obs2",12,0.5,12.5,200,0,2.0,-100,100,"");
        TProfile2D *pParity_Q_ss_obs1   = new TProfile2D("Parity_Q_ss_obs1","Parity_Q_ss_obs1",4,0.5,4.5,500,0,5.0,-100,100,"");
        TProfile2D *pParity_Q_ss_obs2   = new TProfile2D("Parity_Q_ss_obs2","Parity_Q_ss_obs2",4,0.5,4.5,500,0,5.0,-100,100,"");

	Int_t nentries = chain->GetEntries();
//	nentries = 50000;
	//loop through events
	for(int i = 0; i < nentries; i++){
		if((i+1)%1000==0) cout<<"Processing entry == "<< i+1 <<" == out of "<<nentries<<".\n";
		chain->GetEntry(i);

		TLeaf* leaf_RunId   = chain->GetLeaf("mRunId");		
		TLeaf* leaf_EventId = chain->GetLeaf("mEventId");
		TLeaf* leaf_Trigger = chain->GetLeaf("mTrigger");
                TLeaf* leaf_Bz	    = chain->GetLeaf("mBz");
                TLeaf* leaf_PrimaryVertexZ = chain->GetLeaf("mPrimaryVertexZ");
                TLeaf* leaf_RefMult = chain->GetLeaf("mRefMult");
//		TLeaf* leaf_Centrality = chain->GetLeaf("mCentrality");
//		TLeaf* leaf_ZDC_EP  = chain->GetLeaf("mZDC_EP");
                TLeaf* leaf_NoTracks = chain->GetLeaf("mNoTracks");
		Run	= (int)leaf_RunId->GetValue(0);
		Event	= (int)leaf_EventId->GetValue(0);
		Trigger = (int)leaf_Trigger->GetValue(0);
		Bz	= leaf_Bz->GetValue(0);
		PVtxz	= leaf_PrimaryVertexZ->GetValue(0);
		RefMult = (int)leaf_RefMult->GetValue(0);
		NPTracks= (int)leaf_NoTracks->GetValue(0);
		zdcCoincidenceRate = 27000;
//		psi_E   = leaf_ZDC_EP->GetValue(0);
//		psi_W	= leaf_ZDC_EP->GetValue(2);
//		mod_E	= leaf_ZDC_EP->GetValue(1);
//                mod_W   = leaf_ZDC_EP->GetValue(3);
//		Centrality = leaf_Centrality->GetValue(0);
		Day 	= (int)((Run-12000000)/1000); 
		Day2    = (int)((Run-12000000)/10);
                Day3    = (int)((Run-12000000)/1);

		if ( refmultCorrUtil.isBadRun(Run) )
			continue;

		refmultCorrUtil.init(Run);
		refmultCorrUtil.initEvent(RefMult,PVtxz,zdcCoincidenceRate);
		Eweight = refmultCorrUtil.getWeight();
		Centrality = 1 + refmultCorrUtil.getCentralityBin9();
           
                // Check the day2 first 
                int bad_flag = 0;
                for(int jj = 0; jj < 12; jj++){
		    if(Day2 == bad_Ref_day2[jj]){
			bad_flag = 1; 
			break;
		    }
		}
                // Check the day3 
                if(Day3 <= 138024){
                    for(int jj = 0; jj < Nrun_MB1; jj++){
                        if(Day3 == bad_Ref_day3_MB1[jj]){
                            bad_flag = 1;
			    break; 
			}
		    }
                }
                else if(Day3 <= 145020){
                    for(int jj = 0; jj < Nrun_MB2; jj++){
                        if(Day3 == bad_Ref_day3_MB2[jj]){
                            bad_flag = 1;
			    break;
			}
		    }
		}
                else if(Day3 <= 154021){}
                else if(Day3 <= 165031){
                    for(int jj = 0; jj < Nrun_MB5; jj++){
                        if(Day3 == bad_Ref_day3_MB5[jj]){
			    bad_flag = 1;
			    break;
			}
		    }
		}
                else{
                    for(int jj = 0; jj < Nrun_MB6; jj++){
                        if(Day3 == bad_Ref_day3_MB6[jj]){
			    bad_flag = 1;
			    break;
			}
                    }
		}
		if(bad_flag) continue; 
//		if((Trigger%100)>=50) Eweight = Mult_Weight->GetBinContent(RefMult+1,int((PVtxz+30)/2)+1);
//		Eweight = 1.0/(1.0 - exp(-0.92*pow(RefMult,0.43)));
//		if(isnan(Eweight)) continue;
		hBz->Fill(Bz);
                hTrigger->Fill(Trigger);
		if((Trigger%100)>=0) {hMult_Vz->Fill(RefMult,PVtxz); hMult_Vz_new->Fill(RefMult,PVtxz,Eweight);}
		hCentrality->Fill(Centrality);
		hVertexZ->Fill(PVtxz);
		if(cen && Centrality != cen) continue;
		hEventTally->Fill("Total Event",1);
                if(TMath::Abs(PVtxz) > Vz_cut) continue;              //Z-vertex cut; track quality cut done in PkNtupleMaker
/*
//Corrections
                if(fWgt->IsOpen() && (TPCmean_FF->GetEntries() || TPCmean_RF->GetEntries())) {
                        for(int k=1;k<5;k++) {
                                if(Bz>0) PhiMean[k-1] = TPCmean_FF->GetBinContent(k,Day2-7999);
                                if(Bz<0) PhiMean[k-1] = TPCmean_RF->GetBinContent(k,Day2-7999);                        
			}
if((i+1)%1000==0) cout<<Day2<<"shift = "<<PhiMean[0]<<endl;
                }
*/

                TLeaf* leaf_PtV0        = chain->GetLeaf("fV0s.mPtprimaryV0");
                TLeaf* leaf_EtaV0       = chain->GetLeaf("fV0s.mEtaV0");
                TLeaf* leaf_PhiV0       = chain->GetLeaf("fV0s.mPhiV0");
                TLeaf* leaf_ChargeV0    = chain->GetLeaf("fV0s.mChargeV0");
                TLeaf* leaf_DCAglobalV0 = chain->GetLeaf("fV0s.mDCAglobalV0");
                TLeaf* leaf_ndEdxV0     = chain->GetLeaf("fV0s.mdEdxV0");
                TLeaf* leaf_nSigmaPV0   = chain->GetLeaf("fV0s.mnSigmaPV0");
                TLeaf* leaf_nSigmaPiV0  = chain->GetLeaf("fV0s.mnSigmaPiV0");
                TLeaf* leaf_nSigmaKV0   = chain->GetLeaf("fV0s.mnSigmaKV0");
		TLeaf* leaf_TofflagV0	= chain->GetLeaf("fV0s.mTofflagV0");
		TLeaf* leaf_TofV0   	= chain->GetLeaf("fV0s.mTofV0");
		TLeaf* leaf_PathlenV0	= chain->GetLeaf("fV0s.mPathlenV0");
		TLeaf* leaf_PV0		= chain->GetLeaf("fV0s.mPV0");
		TLeaf* leaf_trackIdV0   = chain->GetLeaf("fV0s.mTrackIdV0");
//Lambda Leaf
                TLeaf* leaf_PtW0        = chain->GetLeaf("fW0s.mPtglobalW0");
                TLeaf* leaf_PxW0        = chain->GetLeaf("fW0s.mPxglobalW0");
                TLeaf* leaf_PyW0        = chain->GetLeaf("fW0s.mPyglobalW0");
                TLeaf* leaf_PzW0        = chain->GetLeaf("fW0s.mPzglobalW0");
                TLeaf* leaf_EtaW0       = chain->GetLeaf("fW0s.mEtaW0");
                TLeaf* leaf_DCAglobalW0 = chain->GetLeaf("fW0s.mDcaW0");
                TLeaf* leaf_DecayLength = chain->GetLeaf("fW0s.mDecaylenW0");
                TLeaf* leaf_Dau1Dca     = chain->GetLeaf("fW0s.mDau1dcaW0");
                TLeaf* leaf_Dau2Dca     = chain->GetLeaf("fW0s.mDau2dcaW0");
                TLeaf* leaf_Dca1to2     = chain->GetLeaf("fW0s.mDca1to2W0");
		TLeaf* leaf_Dau1PrMatch = chain->GetLeaf("fW0s.mDau1PrMatchW0");
		TLeaf* leaf_Dau2PrMatch = chain->GetLeaf("fW0s.mDau2PrMatchW0");
		TLeaf* leaf_Dau1PrPx	= chain->GetLeaf("fW0s.mDau1px_prW0");
		TLeaf* leaf_Dau1PrPy    = chain->GetLeaf("fW0s.mDau1py_prW0");
		TLeaf* leaf_Dau1PrPz    = chain->GetLeaf("fW0s.mDau1pz_prW0");
		TLeaf* leaf_Dau1PrPt    = chain->GetLeaf("fW0s.mDau1pt_prW0");
		TLeaf* leaf_Dau2PrPx    = chain->GetLeaf("fW0s.mDau2px_prW0");
		TLeaf* leaf_Dau2PrPy    = chain->GetLeaf("fW0s.mDau2py_prW0");
		TLeaf* leaf_Dau2PrPz    = chain->GetLeaf("fW0s.mDau2pz_prW0");
		TLeaf* leaf_Dau2PrPt    = chain->GetLeaf("fW0s.mDau2pt_prW0");
                TLeaf* leaf_Dau1nSigma = chain->GetLeaf("fW0s.mDau1nSigmaW0");
		TLeaf* leaf_Dau2nSigma = chain->GetLeaf("fW0s.mDau2nSigmaW0");
		TLeaf* leaf_NoLambda 	= chain->GetLeaf("mNoLambdas");
		TLeaf* leaf_LambdaMass	= chain->GetLeaf("fW0s.mMassW0");
		TLeaf* leaf_Dau1Id      = chain->GetLeaf("fW0s.mDau1idW0");
		TLeaf* leaf_Dau2Id      = chain->GetLeaf("fW0s.mDau2idW0");
		int Np = 0, Npbar = 0;
		//TPC EP reconstruction
		TVector2 mQ, mQ1, mQ2, mLambdaPhi;
		Double_t mQx=0., mQy=0., mQx1=0., mQy1=0., mQx2=0., mQy2=0.;
		int Fcount = 0, Ecount = 0, Wcount =0;
		int NLambda = (int)leaf_NoLambda ->GetValue(0);
		for(int trk = 0; trk < NPTracks; trk++){
		    float EtaAsso   = leaf_EtaV0->GetValue(trk);
		    float PtAsso    = leaf_PtV0->GetValue(trk);
		    float PhiAsso   = leaf_PhiV0->GetValue(trk);
		    float DCAglAsso = leaf_DCAglobalV0->GetValue(trk);
		    float ChargeAsso= leaf_ChargeV0->GetValue(trk);
		    float nSigma_p  = leaf_nSigmaPV0->GetValue(trk);
		    float nSigma_K  = leaf_nSigmaKV0->GetValue(trk);
		    int trackID = leaf_trackIdV0 -> GetValue(trk);
		    int TOF_flag = leaf_TofflagV0 -> GetValue(trk);
		    //if(nSigma_p>-3)continue;
		    //		if(ChargeAsso>0)SigmaProtonPlus->Fill(PtAsso,nSigma_p);
		    //		if(ChargeAsso<0)SigmaProtonMinus->Fill(PtAsso,nSigma_p);
		    // also nSigna_e>1
		    //		if(PtAsso>0.4 && PtAsso<1 && DCAglAsso<1 && EtaAsso>-1 && EtaAsso<1 && nSigma_p>-2 && nSigma_p<2 && nSigma_K>2 && ChargeAsso >0) Np++;
		    //              if(PtAsso>0.4 && PtAsso<1 && DCAglAsso<1 && EtaAsso>-1 && EtaAsso<1 && nSigma_p>-2 && nSigma_p<2 && nSigma_K>2 && ChargeAsso <0) Npbar++;

		    if(PtAsso > pt_asso_up || PtAsso < pt_asso_lo) continue;
		    if(DCAglAsso > DcaCut) continue;
		    if(EtaAsso>EtaCut || EtaAsso<-EtaCut) continue;
		    if(nSigma_p > -2 && TOF_flag > 0) continue;
		    int kUse = 1;
		    for(int trkj = 0; trkj < NLambda; trkj++){
			int Dau1Id = leaf_Dau1Id -> GetValue(trkj);
			int Dau2Id = leaf_Dau2Id -> GetValue(trkj);
			if(trackID == Dau1Id || trackID == Dau2Id) kUse = 0;
		    }
		    if(kUse == 0) continue;
		    Fcount++;
		    int n = (int)((PhiAsso+PI)/2./PI*Phibin);
		    int fl = (PVtxz > 0)? ((ChargeAsso > 0)? 1:2):((ChargeAsso > 0)? 3:4);
		    Hist_Phi->Fill(PhiAsso,fl,PtAsso);
		    float W_phi = (EtaAsso>0)? PhiWgtFF[n][fl-1]:PhiWgtRF[n][fl-1];
		    mQx += W_phi*PtAsso * cos(PhiAsso * 2.);
		    mQy += W_phi*PtAsso * sin(PhiAsso * 2.);
		}

		int iTrack[Fcount], Scount = Fcount/2 -1;
		for(int q=0;q<Fcount;q++) iTrack[q] = q;
		random_shuffle(iTrack,iTrack+Fcount);
		Fcount = 0;
		for(int trk = 0; trk < NPTracks; trk++){
		    float EtaAsso   = leaf_EtaV0->GetValue(trk);
		    float PtAsso    = leaf_PtV0->GetValue(trk);
		    float PhiAsso   = leaf_PhiV0->GetValue(trk);
		    float DCAglAsso = leaf_DCAglobalV0->GetValue(trk);
		    float ChargeAsso= leaf_ChargeV0->GetValue(trk);
		    float nSigma_p  = leaf_nSigmaPV0->GetValue(trk);
		    int TOF_flag = leaf_TofflagV0 -> GetValue(trk);
                    int trackID = leaf_trackIdV0 -> GetValue(trk);
		    if(PtAsso > pt_asso_up || PtAsso < pt_asso_lo) continue;
		    if(DCAglAsso > DcaCut) continue;
		    if(EtaAsso>EtaCut || EtaAsso<-EtaCut) continue;
		    if(nSigma_p > -2 && TOF_flag > 0) continue;

		    int kUse = 1;
		    for(int trkj = 0; trkj < NLambda; trkj++){
			int Dau1Id = leaf_Dau1Id -> GetValue(trkj);
			int Dau2Id = leaf_Dau2Id -> GetValue(trkj);
			if(trackID == Dau1Id || trackID == Dau2Id) kUse = 0;
		    }
		    if(kUse == 0) continue;

		    //if(nSigma_p>-3)continue;
		    int n = (int)((PhiAsso+PI)/2./PI*Phibin);
		    int fl = (PVtxz > 0)? ((ChargeAsso > 0)? 1:2):((ChargeAsso > 0)? 3:4);
		    float W_phi = (EtaAsso>0)? PhiWgtFF[n][fl-1]:PhiWgtRF[n][fl-1];
		    if(iTrack[Fcount] > Scount) {mQx1 +=W_phi*PtAsso*cos(PhiAsso*2.); mQy1 +=W_phi*PtAsso*sin(PhiAsso * 2.); Ecount++;}
		    else {mQx2 += W_phi*PtAsso * cos(PhiAsso * 2.); mQy2 += W_phi*PtAsso * sin(PhiAsso * 2.); Wcount++;}
		    Fcount++;
		}
		mQ.Set(mQx, mQy); mQ1.Set(mQx1, mQy1); mQ2.Set(mQx2, mQy2);
		float TPC_EP_full = 0.5*mQ.Phi();
		float TPC_EP_east = 0.5*mQ1.Phi();
		float TPC_EP_west = 0.5*mQ2.Phi();
		Hist_TPC_EP_full->Fill(TPC_EP_full,Day);
		Hist_TPC_EP_east->Fill(TPC_EP_east,Day);
		Hist_TPC_EP_west->Fill(TPC_EP_west,Day);
		pTPC_EP_east->Fill(1,Day,cos(2*TPC_EP_east)); pTPC_EP_east->Fill(2,Day,sin(2*TPC_EP_east)); 
		pTPC_EP_east->Fill(3,Day,cos(4*TPC_EP_east)); pTPC_EP_east->Fill(4,Day,sin(4*TPC_EP_east));
		pTPC_EP_east->Fill(5,Day,cos(6*TPC_EP_east)); pTPC_EP_east->Fill(6,Day,sin(6*TPC_EP_east));
		pTPC_EP_east->Fill(7,Day,cos(8*TPC_EP_east)); pTPC_EP_east->Fill(8,Day,sin(8*TPC_EP_east));

		pTPC_EP_west->Fill(1,Day,cos(2*TPC_EP_west)); pTPC_EP_west->Fill(2,Day,sin(2*TPC_EP_west));
		pTPC_EP_west->Fill(3,Day,cos(4*TPC_EP_west)); pTPC_EP_west->Fill(4,Day,sin(4*TPC_EP_west));
		pTPC_EP_west->Fill(5,Day,cos(6*TPC_EP_west)); pTPC_EP_west->Fill(6,Day,sin(6*TPC_EP_west));
		pTPC_EP_west->Fill(7,Day,cos(8*TPC_EP_west)); pTPC_EP_west->Fill(8,Day,sin(8*TPC_EP_west));

		if(fWgt->IsOpen() && Read_TPC_EP_east->GetEntries()) {
		    PsiShiftE1 = Read_TPC_EP_east->GetBinContent(1,Day-79); PsiShiftE2 = Read_TPC_EP_east->GetBinContent(2,Day-79);
		    PsiShiftE3 = Read_TPC_EP_east->GetBinContent(3,Day-79); PsiShiftE4 = Read_TPC_EP_east->GetBinContent(4,Day-79);
		    PsiShiftE5 = Read_TPC_EP_east->GetBinContent(5,Day-79); PsiShiftE6 = Read_TPC_EP_east->GetBinContent(6,Day-79);
		    PsiShiftE7 = Read_TPC_EP_east->GetBinContent(7,Day-79); PsiShiftE8 = Read_TPC_EP_east->GetBinContent(8,Day-79);
		}
		if(fWgt->IsOpen() && Read_TPC_EP_west->GetEntries()) {
		    PsiShiftW1 = Read_TPC_EP_west->GetBinContent(1,Day-79); PsiShiftW2 = Read_TPC_EP_west->GetBinContent(2,Day-79);
		    PsiShiftW3 = Read_TPC_EP_west->GetBinContent(3,Day-79); PsiShiftW4 = Read_TPC_EP_west->GetBinContent(4,Day-79);
		    PsiShiftW5 = Read_TPC_EP_west->GetBinContent(5,Day-79); PsiShiftW6 = Read_TPC_EP_west->GetBinContent(6,Day-79);
		    PsiShiftW7 = Read_TPC_EP_west->GetBinContent(7,Day-79); PsiShiftW8 = Read_TPC_EP_west->GetBinContent(8,Day-79);
		}
		float TPC_EP_east_new = TPC_EP_east, TPC_EP_west_new = TPC_EP_west;
		TPC_EP_east_new += 2*(-PsiShiftE2*cos(2*TPC_EP_east)+PsiShiftE1*sin(2*TPC_EP_east))/(float)2 + 2*(-PsiShiftE4*cos(4*TPC_EP_east)+PsiShiftE3*sin(4*TPC_EP_east))/(float)4 + 2*(-PsiShiftE6*cos(6*TPC_EP_east)+PsiShiftE5*sin(6*TPC_EP_east))/(float)6 + 2*(-PsiShiftE8*cos(8*TPC_EP_east)+PsiShiftE7*sin(8*TPC_EP_east))/(float)8;
		TPC_EP_west_new += 2*(-PsiShiftW2*cos(2*TPC_EP_west)+PsiShiftW1*sin(2*TPC_EP_west))/(float)2 + 2*(-PsiShiftW4*cos(4*TPC_EP_west)+PsiShiftW3*sin(4*TPC_EP_west))/(float)4 + 2*(-PsiShiftW6*cos(6*TPC_EP_west)+PsiShiftW5*sin(6*TPC_EP_west))/(float)6 + 2*(-PsiShiftW8*cos(8*TPC_EP_west)+PsiShiftW7*sin(8*TPC_EP_west))/(float)8;
		if(TPC_EP_east_new>PI) TPC_EP_east_new -= PI;
		if(TPC_EP_east_new< 0) TPC_EP_east_new += PI;
		if(TPC_EP_west_new>PI) TPC_EP_west_new -= PI;
		if(TPC_EP_west_new< 0) TPC_EP_west_new += PI;
		Hist_TPC_EP_east_flat->Fill(TPC_EP_east_new,Day);
		Hist_TPC_EP_west_flat->Fill(TPC_EP_west_new,Day);
		mQx = mQ1.Mod()*cos(2*TPC_EP_east_new) + mQ2.Mod()*cos(2*TPC_EP_west_new);
		mQy = mQ1.Mod()*sin(2*TPC_EP_east_new) + mQ2.Mod()*sin(2*TPC_EP_west_new);
		mQ.Set(mQx, mQy);
		TPC_EP_full = 0.5*mQ.Phi();
		pTPC_EP_full->Fill(1,Day,cos(2*TPC_EP_full)); pTPC_EP_full->Fill(2,Day,sin(2*TPC_EP_full));
		pTPC_EP_full->Fill(3,Day,cos(4*TPC_EP_full)); pTPC_EP_full->Fill(4,Day,sin(4*TPC_EP_full));
		pTPC_EP_full->Fill(5,Day,cos(6*TPC_EP_full)); pTPC_EP_full->Fill(6,Day,sin(6*TPC_EP_full));
		pTPC_EP_full->Fill(7,Day,cos(8*TPC_EP_full)); pTPC_EP_full->Fill(8,Day,sin(8*TPC_EP_full));

		if(fWgt->IsOpen() && Read_TPC_EP_full->GetEntries()) {
		    PsiShiftF1 = Read_TPC_EP_full->GetBinContent(1,Day-79); PsiShiftF2 = Read_TPC_EP_full->GetBinContent(2,Day-79);
		    PsiShiftF3 = Read_TPC_EP_full->GetBinContent(3,Day-79); PsiShiftF4 = Read_TPC_EP_full->GetBinContent(4,Day-79);
		    PsiShiftF5 = Read_TPC_EP_full->GetBinContent(5,Day-79); PsiShiftF6 = Read_TPC_EP_full->GetBinContent(6,Day-79);
		    PsiShiftF7 = Read_TPC_EP_full->GetBinContent(7,Day-79); PsiShiftF8 = Read_TPC_EP_full->GetBinContent(8,Day-79);
		}
		float TPC_EP_full_new = TPC_EP_full;
		TPC_EP_full_new += 2*(-PsiShiftF2*cos(2*TPC_EP_full)+PsiShiftF1*sin(2*TPC_EP_full))/(float)2 + 2*(-PsiShiftF4*cos(4*TPC_EP_full)+PsiShiftF3*sin(4*TPC_EP_full))/(float)4 + 2*(-PsiShiftF6*cos(6*TPC_EP_full)+PsiShiftF5*sin(6*TPC_EP_full))/(float)6 + 2*(-PsiShiftF8*cos(8*TPC_EP_full)+PsiShiftF7*sin(8*TPC_EP_full))/(float)8;
		if(TPC_EP_full_new>PI) TPC_EP_full_new -= PI;
		if(TPC_EP_full_new< 0) TPC_EP_full_new += PI;
		Hist_TPC_EP_full_flat->Fill(TPC_EP_full_new,Day);
		mQx = mQ.Mod()*cos(2*TPC_EP_full_new);
		mQy = mQ.Mod()*sin(2*TPC_EP_full_new);
		Hist_cos->Fill(2,cos(2.*TPC_EP_east_new-2.*TPC_EP_west_new), Eweight);
		Hist_dif_count->Fill(Ecount - Wcount);
		Hist_ful_count->Fill(Ecount + Wcount);
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int N_L_P1 = 0, N_L_N1 = 0, N_R_P1 = 0, N_R_N1 = 0, N_T_P1 = 0, N_T_N1 = 0, N_B_P1 = 0, N_B_N1 = 0;
		int N_L_P2 = 0, N_L_N2 = 0, N_R_P2 = 0, N_R_N2 = 0, N_T_P2 = 0, N_T_N2 = 0, N_B_P2 = 0, N_B_N2 = 0;
		Fcount = 0;
		//**************************************************************************

		//***************************************************************************

		//loop through matched primary tracks
		for(int trki = 0; trki < NPTracks; trki++){
		    Pt	  = leaf_PtV0->GetValue(trki);
		    Eta	  = leaf_EtaV0->GetValue(trki);
		    Theta     = 2.*atan(exp(-Eta));
		    Charge	  = leaf_ChargeV0->GetValue(trki);
		    Phi	  = leaf_PhiV0->GetValue(trki);
		    ndEdx	  = leaf_ndEdxV0->GetValue(trki);
		    DCAGlobal = leaf_DCAglobalV0->GetValue(trki);
		    nSigma_p  = leaf_nSigmaPV0->GetValue(trki);
		    nSigma_pi = leaf_nSigmaPiV0->GetValue(trki);
		    int trackID     = leaf_trackIdV0->GetValue(trki);

		    if(nSigma_p<-3)
		    {
			float En  = sqrt(0.1396*0.1396+pow(Pt*cosh(Eta),2));

			float mQx_i = mQx, mQy_i = mQy;
			int n = (int)((Phi+PI)/2./PI*Phibin);
			int fl = (PVtxz > 0)? ((Charge > 0)? 1:2):((Charge > 0)? 3:4);
			float W_phi = (Eta>0)? PhiWgtFF[n][fl-1]:PhiWgtRF[n][fl-1];
			if(Pt > pt_asso_lo && Pt < pt_asso_up && Eta < EtaCut && Eta > -EtaCut && DCAGlobal < DcaCut) {
			    mQx_i -= W_phi*Pt * cos(Phi * 2.);
			    mQy_i -= W_phi*Pt * sin(Phi * 2.);
			}
			TVector2 mQ_i(mQx_i,mQy_i);
			float psi_F = 0.5*mQ_i.Phi();
			float v2 = cos(2*Phi - 2*psi_F)*100;
			Hist_v2_pt_obs->Fill(Pt,v2);
			if(Charge>0 && Eta > -1 && Eta <1 && ndEdx>10 && nSigma_pi>-2 && nSigma_pi<2) Hist_v2_pt_pos_obs->Fill(Pt,v2);
			if(Charge<0 && Eta > -1 && Eta <1 && ndEdx>10 && nSigma_pi>-2 && nSigma_pi<2) Hist_v2_pt_neg_obs->Fill(Pt,v2);

			Hist_Pt->Fill(Pt,Eweight);
			if(DCAGlobal > DcaCut) continue;
			hEtaPtDist->Fill(Eta,Pt,Eweight);

			if(Pt < pt_trig_lo || Pt > pt_trig_up) continue;
			if(Eta > EtaCut || Eta < -EtaCut) continue;

			//Corrections
			if(fWgt->IsOpen() && (TPCmean_FF->GetEntries() || TPCmean_RF->GetEntries())) {
			    for(int k=1;k<5;k++) {
				for(int kk=0;kk<order;kk++) {
				    if(Eta>0 && PVtxz>0) PhiMean[k-1+4*kk] = TPCmean_FF->GetBinContent(k+8*kk,Day2-7999);
				    if(Eta<0 && PVtxz>0) PhiMean[k-1+4*kk] = TPCmean_RF->GetBinContent(k+8*kk,Day2-7999);
				    if(Eta>0 && PVtxz<0) PhiMean[k-1+4*kk] = TPCmean_FF->GetBinContent(k+4+8*kk,Day2-7999);
				    if(Eta<0 && PVtxz<0) PhiMean[k-1+4*kk] = TPCmean_RF->GetBinContent(k+4+8*kk,Day2-7999);
				}
			    }
			}

			float cos1 =0, cos2=0, sin1=0, sin2=0, Phi_new = Phi;
			// Recentering parameters
			if(Eta>0 && PVtxz>0) {
			    if(Charge > 0) {
				for(int kk=0;kk<order;kk++) {
				    pTPCmeanPhi_FF->Fill(1+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
				    pTPCmeanPhi_FF->Fill(2+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
				}
				Hist_Phi_FF->Fill(Phi,1,Eweight);}
			    if(Charge < 0) {
				for(int kk=0;kk<order;kk++) {
				    pTPCmeanPhi_FF->Fill(3+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
				    pTPCmeanPhi_FF->Fill(4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
				}
				Hist_Phi_FF->Fill(Phi,2,Eweight);}
			}
			if(Eta<0 && PVtxz>0) {
			    if(Charge > 0) {
				for(int kk=0;kk<order;kk++) {
				    pTPCmeanPhi_RF->Fill(1+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
				    pTPCmeanPhi_RF->Fill(2+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
				}
				Hist_Phi_RF->Fill(Phi,1,Eweight);}
			    if(Charge < 0) {
				for(int kk=0;kk<order;kk++) {
				    pTPCmeanPhi_RF->Fill(3+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
				    pTPCmeanPhi_RF->Fill(4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
				}
				Hist_Phi_RF->Fill(Phi,2,Eweight);}
			}
			if(Eta>0 && PVtxz<0) {
			    if(Charge > 0) {
				for(int kk=0;kk<order;kk++) {
				    pTPCmeanPhi_FF->Fill(1+4+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
				    pTPCmeanPhi_FF->Fill(2+4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
				}
				Hist_Phi_FF->Fill(Phi,1+2,Eweight);}
			    if(Charge < 0) {
				for(int kk=0;kk<order;kk++) {
				    pTPCmeanPhi_FF->Fill(3+4+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
				    pTPCmeanPhi_FF->Fill(4+4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
				}
				Hist_Phi_FF->Fill(Phi,2+2,Eweight);}
			}
			if(Eta<0 && PVtxz<0) {
			    if(Charge > 0) {
				for(int kk=0;kk<order;kk++) {
				    pTPCmeanPhi_RF->Fill(1+4+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
				    pTPCmeanPhi_RF->Fill(2+4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
			  }
			  Hist_Phi_RF->Fill(Phi,1+2,Eweight);}
		      if(Charge < 0) {
			  for(int kk=0;kk<order;kk++) {
			      pTPCmeanPhi_RF->Fill(3+4+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
			      pTPCmeanPhi_RF->Fill(4+4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
			  }
			  Hist_Phi_RF->Fill(Phi,2+2,Eweight);}
		  }

		  if(Charge > 0) {cos1 = cos(Phi) - PhiMean[0]; sin1 = sin(Phi) - PhiMean[1];
		      double a2np = 1+PhiMean[0+4*(order-1)], a2nn = 1-PhiMean[0+4*(order-1)];
		      double lambda_2nsp = PhiMean[1+4*(order-1)]/a2np;
		      double lambda_2nsn = PhiMean[1+4*(order-1)]/a2nn;
		      double cos11 = (cos1 - lambda_2nsn*sin1)/(1-lambda_2nsn*lambda_2nsp);
		      double sin11 = (sin1 - lambda_2nsp*cos1)/(1-lambda_2nsn*lambda_2nsp);
		      cos1 = cos11/a2np;
		      sin1 = sin11/a2nn;
		      for(int jj=0;jj<order;jj++) Phi_new += -2*PhiMean[1+4*jj]*cos(jj*Phi+Phi)/int(jj+1) 
			  +2*PhiMean[0+4*jj]*sin(jj*Phi+Phi)/int(jj+1); 
		      if(Phi_new> PI) Phi_new -= 2*PI;
		      if(Phi_new<-PI) Phi_new += 2*PI;
		      if(Eta>0 && PVtxz>0) Hist_Phi_FF_new->Fill(Phi_new,1,Eweight);
		      if(Eta<0 && PVtxz>0) Hist_Phi_RF_new->Fill(Phi_new,1,Eweight);
		      if(Eta>0 && PVtxz<0) Hist_Phi_FF_new->Fill(Phi_new,1+2,Eweight);
		      if(Eta<0 && PVtxz<0) Hist_Phi_RF_new->Fill(Phi_new,1+2,Eweight);}
		  if(Charge < 0) {cos1 = cos(Phi) - PhiMean[2]; sin1 = sin(Phi) - PhiMean[3];
		      double a2np = 1+PhiMean[2+4*(order-1)], a2nn = 1-PhiMean[2+4*(order-1)];
		      double lambda_2nsp = PhiMean[3+4*(order-1)]/a2np;
		      double lambda_2nsn = PhiMean[3+4*(order-1)]/a2nn;
		      double cos11 = (cos1 - lambda_2nsn*sin1)/(1-lambda_2nsn*lambda_2nsp);
		      double sin11 = (sin1 - lambda_2nsp*cos1)/(1-lambda_2nsn*lambda_2nsp);
		      cos1 = cos11/a2np;
		      sin1 = sin11/a2nn;
		      for(int jj=0;jj<order;jj++) Phi_new += -2*PhiMean[3+4*jj]*cos(jj*Phi+Phi)/int(jj+1) 
			  +2*PhiMean[2+4*jj]*sin(jj*Phi+Phi)/int(jj+1);
		      if(Phi_new> PI) Phi_new -= 2*PI;
		      if(Phi_new<-PI) Phi_new += 2*PI;
		      if(Eta>0 && PVtxz>0) Hist_Phi_FF_new->Fill(Phi_new,2,Eweight);
		      if(Eta<0 && PVtxz>0) Hist_Phi_RF_new->Fill(Phi_new,2,Eweight);
		      if(Eta>0 && PVtxz<0) Hist_Phi_FF_new->Fill(Phi_new,2+2,Eweight);
		      if(Eta<0 && PVtxz<0) Hist_Phi_RF_new->Fill(Phi_new,2+2,Eweight);}

		  if(iTrack[Fcount] > Scount && Charge>0) {
		      if(sin(Phi_new - TPC_EP_west_new)>0) N_T_P1++;
		      if(sin(Phi_new - TPC_EP_west_new)<0) N_B_P1++;
		      if(cos(Phi_new - TPC_EP_west_new)>0) N_R_P1++;
		      if(cos(Phi_new - TPC_EP_west_new)<0) N_L_P1++;
		  }
		  if(iTrack[Fcount] > Scount && Charge<0) {
		      if(sin(Phi_new - TPC_EP_west_new)>0) N_T_N1++;
		      if(sin(Phi_new - TPC_EP_west_new)<0) N_B_N1++;
		      if(cos(Phi_new - TPC_EP_west_new)>0) N_R_N1++;
		      if(cos(Phi_new - TPC_EP_west_new)<0) N_L_N1++;
		  }
		  if(iTrack[Fcount] <= Scount && Charge>0) {
		      if(sin(Phi_new - TPC_EP_east_new)>0) N_T_P2++;
		      if(sin(Phi_new - TPC_EP_east_new)<0) N_B_P2++;
		      if(cos(Phi_new - TPC_EP_east_new)>0) N_R_P2++;
		      if(cos(Phi_new - TPC_EP_east_new)<0) N_L_P2++;
		  }
		  if(iTrack[Fcount] <= Scount && Charge<0) {
		      if(sin(Phi_new - TPC_EP_east_new)>0) N_T_N2++;
		      if(sin(Phi_new - TPC_EP_east_new)<0) N_B_N2++;
		      if(cos(Phi_new - TPC_EP_east_new)>0) N_R_N2++;
		      if(cos(Phi_new - TPC_EP_east_new)<0) N_L_N2++;
		  }
		  Fcount++;
	      }

	      //Second Loop
	      int TOF_flag = leaf_TofflagV0->GetValue(trki);
	      double TOF_proton = leaf_TofV0->GetValue(trki);
	      double TOF_path = leaf_PathlenV0->GetValue(trki);
	      double proton_p = leaf_PV0->GetValue(trki);
	      double mass2proton = 0;
	      if(TOF_flag>0)
	      {  mass2proton = proton_p*proton_p*(900.0*TOF_proton*TOF_proton/TOF_path/TOF_path - 1.0);  }

	      if(nSigma_p>-2&&nSigma_p<2&&mass2proton>0.8&&mass2proton<1&&Pt>0.4)
	      {
		  float En_proton  = sqrt(0.938272013*0.938272013+pow(Pt*cosh(Eta),2));
		  float Theta_proton     = 2.*atan(exp(-Eta));
		  float Phi_proton       = Phi;
		  if(Phi_proton> PI) Phi_proton -= 2*PI;
		  if(Phi_proton<-PI) Phi_proton += 2*PI;

		  if(DCAGlobal > DcaCut) continue;
		  hEtaPt_Proton_Dist->Fill(Eta,Pt,Eweight);
		  if(Pt > pt_trig_up) continue;
		  if(Eta > EtaCut || Eta < -EtaCut) continue;

		  if(fWgt_pr->IsOpen() && (TPCmean_FF_Proton->GetEntries() || TPCmean_RF_Proton->GetEntries())) {
		      for(int k=1;k<5;k++) {
			  for(int kk=0;kk<order;kk++) {
			      if(Eta>0 && PVtxz>0) PhiMean_Proton[k-1+4*kk] = TPCmean_FF_Proton->GetBinContent(k+8*kk,Day2-7999);
			      if(Eta<0 && PVtxz>0) PhiMean_Proton[k-1+4*kk] = TPCmean_RF_Proton->GetBinContent(k+8*kk,Day2-7999);
			      if(Eta>0 && PVtxz<0) PhiMean_Proton[k-1+4*kk] = TPCmean_FF_Proton->GetBinContent(k+4+8*kk,Day2-7999);
			      if(Eta<0 && PVtxz<0) PhiMean_Proton[k-1+4*kk] = TPCmean_RF_Proton->GetBinContent(k+4+8*kk,Day2-7999);
			  }
		      }
		  }
		  int n = (int)((Phi_proton+PI)/2./PI*Phibin);
		  int fl = (PVtxz > 0)? ((Charge > 0)? 1:2):((Charge > 0)? 3:4);
		  float W_phi1 = (Eta>0)? PhiWgtFF_Proton[n][fl-1]:PhiWgtRF_Proton[n][fl-1];

		  float cos1 =0, cos2=0, sin1=0, sin2=0, Phi_new = Phi_proton;

		  if(Eta>0 && PVtxz>0) {
		      if(Charge > 0) {
			  for(int kk=0;kk<order;kk++) {
			      pTPCmeanPhi_FF_Proton->Fill(1+8*kk,Day2,cos(kk*Phi_proton+Phi_proton),Eweight);
			      pTPCmeanPhi_FF_Proton->Fill(2+8*kk,Day2,sin(kk*Phi_proton+Phi_proton),Eweight);
			  }
			  Hist_Phi_FF_Proton->Fill(Phi_proton,1,Eweight);
		      }
		      if(Charge < 0) {
			  for(int kk=0;kk<order;kk++) {
			      pTPCmeanPhi_FF_Proton->Fill(3+8*kk,Day2,cos(kk*Phi_proton+Phi_proton),Eweight);
			      pTPCmeanPhi_FF_Proton->Fill(4+8*kk,Day2,sin(kk*Phi_proton+Phi_proton),Eweight);
			  }
			  Hist_Phi_FF_Proton->Fill(Phi_proton,2,Eweight);
		      }
		  }
		  if(Eta<0 && PVtxz>0) {
		      if(Charge > 0) {
			  for(int kk=0;kk<order;kk++) {
			      pTPCmeanPhi_RF_Proton->Fill(1+8*kk,Day2,cos(kk*Phi_proton+Phi_proton),Eweight);
			      pTPCmeanPhi_RF_Proton->Fill(2+8*kk,Day2,sin(kk*Phi_proton+Phi_proton),Eweight);
			  }
			  Hist_Phi_RF_Proton->Fill(Phi_proton,1,Eweight);
		      }
		      if(Charge < 0) {
			  for(int kk=0;kk<order;kk++) {
			      pTPCmeanPhi_RF_Proton->Fill(3+8*kk,Day2,cos(kk*Phi_proton+Phi_proton),Eweight);
			      pTPCmeanPhi_RF_Proton->Fill(4+8*kk,Day2,sin(kk*Phi_proton+Phi_proton),Eweight);
			  }
			  Hist_Phi_RF_Proton->Fill(Phi_proton,2,Eweight);
		      }
		  }
		  if(Eta>0 && PVtxz<0) {
		      if(Charge > 0) {
			  for(int kk=0;kk<order;kk++) {
			      pTPCmeanPhi_FF_Proton->Fill(1+4+8*kk,Day2,cos(kk*Phi_proton+Phi_proton),Eweight);
			      pTPCmeanPhi_FF_Proton->Fill(2+4+8*kk,Day2,sin(kk*Phi_proton+Phi_proton),Eweight);
			  }
			  Hist_Phi_FF_Proton->Fill(Phi_proton,1+2,Eweight);
		      }
		      if(Charge < 0) {
			  for(int kk=0;kk<order;kk++) {
			      pTPCmeanPhi_FF_Proton->Fill(3+4+8*kk,Day2,cos(kk*Phi_proton+Phi_proton),Eweight);
			      pTPCmeanPhi_FF_Proton->Fill(4+4+8*kk,Day2,sin(kk*Phi_proton+Phi_proton),Eweight);
			  }
			  Hist_Phi_FF_Proton->Fill(Phi_proton,2+2,Eweight);
		      }
		  }
		  if(Eta<0 && PVtxz<0) {
		      if(Charge > 0) {
			  for(int kk=0;kk<order;kk++) {
			      pTPCmeanPhi_RF_Proton->Fill(1+4+8*kk,Day2,cos(kk*Phi_proton+Phi_proton),Eweight);
			      pTPCmeanPhi_RF_Proton->Fill(2+4+8*kk,Day2,sin(kk*Phi_proton+Phi_proton),Eweight);
			  }
			  Hist_Phi_RF_Proton->Fill(Phi_proton,1+2,Eweight);
		      }
		      if(Charge < 0) {
			  for(int kk=0;kk<order;kk++) {
			      pTPCmeanPhi_RF_Proton->Fill(3+4+8*kk,Day2,cos(kk*Phi_proton+Phi_proton),Eweight);
			      pTPCmeanPhi_RF_Proton->Fill(4+4+8*kk,Day2,sin(kk*Phi_proton+Phi_proton),Eweight);
			  }
			  Hist_Phi_RF_Proton->Fill(Phi_proton,2+2,Eweight);
		      }
		  }
		  if(Charge > 0)
		  {
		      cos1 = cos(Phi_proton) - PhiMean_Proton[0]; sin1 = sin(Phi_proton) - PhiMean_Proton[1];
		      double a2np = 1+PhiMean_Proton[0+4*(order-1)], a2nn = 1-PhiMean_Proton[0+4*(order-1)];
		      double proton_2nsp = PhiMean_Proton[1+4*(order-1)]/a2np;
		      double proton_2nsn = PhiMean_Proton[1+4*(order-1)]/a2nn;
		      double cos11 = (cos1 - proton_2nsn*sin1)/(1-proton_2nsn*proton_2nsp);
		      double sin11 = (sin1 - proton_2nsp*cos1)/(1-proton_2nsn*proton_2nsp);
		      cos1 = cos11/a2np;
		      sin1 = sin11/a2nn;
		      for(int jj=0;jj<order;jj++) Phi_new += -2*PhiMean_Proton[1+4*jj]*cos(jj*Phi_proton+Phi_proton)/int(jj+1)
			  +2*PhiMean_Proton[0+4*jj]*sin(jj*Phi_proton+Phi_proton)/int(jj+1);
		      if(Phi_new> PI) Phi_new -= 2*PI;
		      if(Phi_new<-PI) Phi_new += 2*PI;
		      if(Eta>0 && PVtxz>0) Hist_Phi_FF_Proton_new->Fill(Phi_new,1,Eweight);
		      if(Eta<0 && PVtxz>0) Hist_Phi_RF_Proton_new->Fill(Phi_new,1,Eweight);
		      if(Eta>0 && PVtxz<0) Hist_Phi_FF_Proton_new->Fill(Phi_new,1+2,Eweight);
		      if(Eta<0 && PVtxz<0) Hist_Phi_RF_Proton_new->Fill(Phi_new,1+2,Eweight);
		  }
		  if(Charge < 0)
		  {
		      cos1 = cos(Phi_proton) - PhiMean_Proton[2]; sin1 = sin(Phi_proton) - PhiMean_Proton[3];
		      double a2np = 1+PhiMean_Proton[2+4*(order-1)], a2nn = 1-PhiMean_Proton[2+4*(order-1)];
		      double proton_2nsp = PhiMean_Proton[3+4*(order-1)]/a2np;
		      double proton_2nsn = PhiMean_Proton[3+4*(order-1)]/a2nn;
		      double cos11 = (cos1 - proton_2nsn*sin1)/(1-proton_2nsn*proton_2nsp);
		      double sin11 = (sin1 - proton_2nsp*cos1)/(1-proton_2nsn*proton_2nsp);
		      cos1 = cos11/a2np;
		      sin1 = sin11/a2nn;
		      for(int jj=0;jj<order;jj++) Phi_new += -2*PhiMean_Proton[3+4*jj]*cos(jj*Phi_proton+Phi_proton)/int(jj+1)
			  +2*PhiMean_Proton[2+4*jj]*sin(jj*Phi_proton+Phi_proton)/int(jj+1);
		      if(Phi_new> PI) Phi_new -= 2*PI;
		      if(Phi_new<-PI) Phi_new += 2*PI;
		      if(Eta>0 && PVtxz>0) Hist_Phi_FF_Proton_new->Fill(Phi_new,2,Eweight);
		      if(Eta<0 && PVtxz>0) Hist_Phi_RF_Proton_new->Fill(Phi_new,2,Eweight);
		      if(Eta>0 && PVtxz<0) Hist_Phi_FF_Proton_new->Fill(Phi_new,2+2,Eweight);
		      if(Eta<0 && PVtxz<0) Hist_Phi_RF_Proton_new->Fill(Phi_new,2+2,Eweight);
		  }

		  //Track 2: Lambda Candidates
				if(!correction){
			int NLambda= (int)leaf_NoLambda->GetValue(0);
			for(int trkj = 0; trkj < NLambda; trkj++) {
                                int Dau1Id      = leaf_Dau1Id->GetValue(trkj);
                                int Dau2Id      = leaf_Dau2Id->GetValue(trkj);
                                if(Dau1Id==trackID || Dau2Id==trackID) continue;
                        	float Pt2        = leaf_PtW0->GetValue(trkj);
                        	float Eta2       = leaf_EtaW0->GetValue(trkj);
				double Px_Lambda = leaf_PxW0->GetValue(trkj);
				double Py_Lambda = leaf_PyW0->GetValue(trkj);
				double Mass_Lambda = leaf_LambdaMass->GetValue(trkj);
				double Track1_nSigma = leaf_Dau1nSigma->GetValue(trkj);
				double Track2_nSigma = leaf_Dau2nSigma->GetValue(trkj);
				if(fabs(Mass_Lambda - 1.115683)>0.004)continue;
				if(fabs(Pt2-Pt)<0.15)continue;
				if(fabs(Eta2-Eta)<0.15)continue;

				TVector2 mLambdaPhi;
				mLambdaPhi.Set(Px_Lambda,Py_Lambda);
                        	float Phi2       = mLambdaPhi.Phi();
				if(Phi2> PI) Phi2 -= 2*PI;
				if(Phi2<-PI) Phi2 += 2*PI;

				float DCAGlobal2 = leaf_DCAglobalW0->GetValue(trkj);
				float En2	 = sqrt(1.115683*1.115683+pow(Pt2*cosh(Eta2),2));
				if(Pt2 < 0.5 || Pt2 > pt_asso_up) continue;
				if(Eta2 > EtaCut || Eta2 < -EtaCut) continue;
				if(DCAGlobal > 1.0) continue;

			if(fabs(Track1_nSigma) > 3.0) continue;
			if(fabs(Track2_nSigma) > 3.0) continue;
                        if((leaf_DecayLength->GetValue(trkj))<6)continue;
                        if((DCAGlobal2) > 0.6)continue;
                        if((leaf_Dau1Dca->GetValue(trkj))<0.6)continue;
                        if((leaf_Dau2Dca->GetValue(trkj))<1.8)continue;
                        if((leaf_Dca1to2->GetValue(trkj))>0.7)continue;

                        if(Pt2<0.6&&(leaf_Dau1Dca->GetValue(trkj))<0.7)continue;
                        if(Pt2<0.6&&(leaf_Dau2Dca->GetValue(trkj))<2.5)continue;

				float mQx_j = mQx, mQy_j = mQy;			
	                        int n2 = (int)((Phi2+PI)/2./PI*Phibin);
        	                int fl2 = (PVtxz > 0)? 1:3;
                	        float W_phi2 = (Eta2>0)? PhiWgtFF_Lambda[n2][fl2-1]:PhiWgtRF_Lambda[n2][fl2-1];
				int flag_dau1 = leaf_Dau1PrMatch->GetValue(trkj);
				float px_dau1 = leaf_Dau1PrPx->GetValue(trkj);
				float py_dau1 = leaf_Dau1PrPy->GetValue(trkj);
				float pz_dau1 = leaf_Dau1PrPz->GetValue(trkj);
				float pt_dau1 = leaf_Dau1PrPt->GetValue(trkj);
                                int flag_dau2 = leaf_Dau2PrMatch->GetValue(trkj);
                                float px_dau2 = leaf_Dau2PrPx->GetValue(trkj);
                                float py_dau2 = leaf_Dau2PrPy->GetValue(trkj);
                                float pz_dau2 = leaf_Dau2PrPz->GetValue(trkj);
				float pt_dau2 = leaf_Dau2PrPt->GetValue(trkj);
				TVector2 mDau1Phi;
				mDau1Phi.Set(px_dau1,py_dau1);
                                TVector2 mDau2Phi;
                                mDau2Phi.Set(px_dau2,py_dau2);
				float Phi_dau1 = mDau1Phi.Phi();
				float Phi_dau2 = mDau2Phi.Phi();
				float eta_dau1 = 0, eta_dau2 = 0;
				if(flag_dau1 == 1) eta_dau1 = 0.5*log( (sqrt(px_dau1*px_dau1+py_dau1*py_dau1+pz_dau1*pz_dau1)+pz_dau1)/(sqrt(px_dau1*px_dau1+py_dau1*py_dau1+pz_dau1*pz_dau1)-pz_dau1) );
				if(flag_dau2 == 1) eta_dau2 = 0.5*log( (sqrt(px_dau2*px_dau2+py_dau2*py_dau2+pz_dau2*pz_dau2)+pz_dau2)/(sqrt(px_dau2*px_dau2+py_dau2*py_dau2+pz_dau2*pz_dau2)-pz_dau2) );
				int fl_dau1 = (PVtxz > 0)? 1:3;
				int fl_dau2 = (PVtxz > 0)? 2:4;
				float W_phi_dau1 = (eta_dau1>0)? PhiWgtFF[n2][fl_dau1-1]:PhiWgtRF[n2][fl_dau1-1];
				float W_phi_dau2 = (eta_dau2>0)? PhiWgtFF[n2][fl_dau2-1]:PhiWgtRF[n2][fl_dau2-1];
                        	mQx_j -= (W_phi_dau1*pt_dau1*cos(Phi_dau1*2.0)*0.0 + W_phi_dau2*pt_dau2*cos(Phi_dau2*2.0)*flag_dau2 );
                        	mQy_j -= (W_phi_dau1*pt_dau1*sin(Phi_dau1*2.0)*0.0 + W_phi_dau2*pt_dau2*sin(Phi_dau2*2.0)*flag_dau2 );
				TVector2 mQ_j(mQx_j, mQy_j);
				float psi_F_new = 0.5*mQ_j.Phi();

	                        float Delt_phi1 = Phi_new - psi_F_new;
        	                if(Delt_phi1>PI) Delt_phi1 -= 2*PI;
                	        if(Delt_phi1<-PI) Delt_phi1 += 2*PI;
                        	Delt_phi1 = (Delt_phi1>0)? psi_F_new + (gRandom->Rndm())*PI:psi_F_new - (gRandom->Rndm())*PI;

				hDpt->Fill(fabs(Pt-Pt2),Eweight);
				float q_inv = pow(Pt*sin(Phi)-Pt2*sin(Phi2),2)+pow(Pt*cos(Phi)-Pt2*cos(Phi2),2)+pow(Pt*sinh(Eta)-Pt2*sinh(Eta2),2)-pow(En_proton-En2,2);
				hQinv2->Fill(q_inv,Eweight);
//				if(q_inv<0) continue;
				q_inv = (q_inv>0)? sqrt(q_inv):0;
				hQinv->Fill(q_inv,Eweight);
//Corrections
	                        if(fWgt->IsOpen() && (TPCmean_FF_Lambda->GetEntries() || TPCmean_RF_Lambda->GetEntries())) {
        	                        for(int k=1;k<5;k++) {
                                          for(int kk=0;kk<order;kk++) {
                                          	if(Eta2>0 && PVtxz>0) PhiMean_Lambda[k-1+4*kk] = TPCmean_FF_Lambda->GetBinContent(k+8*kk,Day2-7999);
                                          	if(Eta2<0 && PVtxz>0) PhiMean_Lambda[k-1+4*kk] = TPCmean_RF_Lambda->GetBinContent(k+8*kk,Day2-7999);
                                          	if(Eta2>0 && PVtxz<0) PhiMean_Lambda[k-1+4*kk] = TPCmean_FF_Lambda->GetBinContent(k+4+8*kk,Day2-7999);
                                          	if(Eta2<0 && PVtxz<0) PhiMean_Lambda[k-1+4*kk] = TPCmean_RF_Lambda->GetBinContent(k+4+8*kk,Day2-7999);
                                          }
                                	}
                        	}

				// Recentering parameters
				float Phi2_new = Phi2;
				 {cos2 = cos(Phi2) - PhiMean_Lambda[0]; sin2 = sin(Phi2) - PhiMean_Lambda[1];
                                       		 double a2np = 1+PhiMean_Lambda[0+4*(order-1)], a2nn = 1-PhiMean_Lambda[0+4*(order-1)];
                                        	 double lambda_2nsp = PhiMean_Lambda[1+4*(order-1)]/a2np;
                                        	 double lambda_2nsn = PhiMean_Lambda[1+4*(order-1)]/a2nn;
                                        	 double cos22 = (cos2 - lambda_2nsn*sin2)/(1-lambda_2nsn*lambda_2nsp);
                                        	 double sin22 = (sin2 - lambda_2nsp*cos2)/(1-lambda_2nsn*lambda_2nsp);
                                        	 cos2 = cos22/a2np;
                                        	 sin2 = sin22/a2nn;
                                        	for(int jj=0;jj<order;jj++) 
							Phi2_new += -2*PhiMean_Lambda[1+4*jj]*cos(jj*Phi2+Phi2)/int(jj+1)
                                                                    +2*PhiMean_Lambda[0+4*jj]*sin(jj*Phi2+Phi2)/int(jj+1);
						}
//                                if(Charge2 < 0) {cos2 = cos(Phi2) - PhiMean_Lambda[2]; sin2 = sin(Phi2) - PhiMean_Lambda[3];
//                                        	double a2np = 1+PhiMean_Lambda[2+4], a2nn = 1-PhiMean_Lambda[2+4];
//                                        	double lambda_2nsp = PhiMean_Lambda[3+4]/a2np;
//                                        	double lambda_2nsn = PhiMean_Lambda[3+4]/a2nn;
//                                        	double cos22 = (cos2 - lambda_2nsn*sin2)/(1-lambda_2nsn*lambda_2nsp);
//                                        	double sin22 = (sin2 - lambda_2nsp*cos2)/(1-lambda_2nsn*lambda_2nsp);
//                                        	cos2 = cos22/a2np;
//                                        	sin2 = sin22/a2nn;
//                                        	for(int jj=0;jj<order;jj++) 
//							Phi2_new += -2*PhiMean_Lambda[3+4*jj]*cos(jj*Phi2+Phi2)/int(jj+1)
//                                                                   +2*PhiMean_Lambda[2+4*jj]*sin(jj*Phi2+Phi2)/int(jj+1);
//						}
                        	float Delt_phi2 = Phi2_new - psi_F_new;
                        	if(Delt_phi2>PI) Delt_phi2 -= 2*PI;
                        	if(Delt_phi2<-PI) Delt_phi2 += 2*PI;
                        	Delt_phi2 = (Delt_phi2>0)? psi_F_new + (gRandom->Rndm())*PI:psi_F_new - (gRandom->Rndm())*PI;
				float correlator0 = cos(Phi + Phi2 - 2*psi_F_new);
				float correlator1 = (cos1*cos(Phi2)-sin1*sin(Phi2))*cos(2*psi_F_new)
                                                        + (sin1*cos(Phi2)+cos1*sin(Phi2))*sin(2*psi_F_new);
				float correlator2 = (cos1*cos2-sin1*sin2)*cos(2*psi_F_new)
							+ (sin1*cos2+cos1*sin2)*sin(2*psi_F_new);
				float correlator3 = cos(Phi_new + Phi2 - 2*psi_F_new);
				float correlator4 = cos(Phi_new + Phi2_new - 2*psi_F_new);
				// Evaluate the efficiency
				float eff_la_val = eff_la -> Eval(Pt2);
                                float eff_tpc_p_val = eff_tpc_p -> Eval(Pt);
                                float eff_tof_p_val = eff_tof_p -> Eval(Pt);
                                float eff_comb = eff_la_val * eff_tpc_p_val * eff_tof_p_val;

				float correlator5 = cos(Delt_phi1 + Delt_phi2 +PI - 2*psi_F_new);
				float correlator6 = cos(Phi_new-psi_F_new)*cos(Phi2_new-psi_F_new);
                                float correlator7 = sin(Phi_new-psi_F_new)*sin(Phi2_new-psi_F_new);
//			cout<<"Phi = "<<Phi<<" ; Phi2 = "<<Phi2<<" ; PhiMean = "<<PhiMean[0]<<" "<<PhiMean[1]<<" "<<PhiMean[2]<<" "<<PhiMean[3]<<" ; PhiMean_Lambda = "<<PhiMean_Lambda[0]<<" "<<PhiMean_Lambda[1]<<" "<<PhiMean_Lambda[2]<<" "<<PhiMean_Lambda[3]<<" "<<endl;
//			cout<<"cos1 = "<<cos1<<" ; sin1 = "<<sin1<<" ; cos2 = "<<cos2<<" ; sin2 = "<<sin2<<endl;
//			cout<<"Phi_new = "<<Phi_new<<" ; Phi2_new = "<<Phi2_new<<" ; psi_F_new ="<<psi_F_new<<endl;
//			cout<<correlator0<<" ; "<<correlator1<<" ; "<<correlator2<<" ; "<<correlator3<<" ; "<<correlator4<<" ; "<<correlator5<<" ; "<<correlator6<<" ; "<<correlator7<<endl;
				if(Charge>0) {
		if(!isnan(correlator0)){	pParity_int_obs1->Fill(1,100*correlator0);
					pParity_int_obs2->Fill(1,100*correlator0,Eweight);		
                                        pParity_int_w_obs1->Fill(1,100*correlator0,W_phi1);
                                        pParity_int_w_obs2->Fill(1,100*correlator0,Eweight*W_phi1);
                                        pParity_int_ww_obs1->Fill(1,100*correlator0,W_phi1*W_phi2);
                                        pParity_int_ww_obs2->Fill(1,100*correlator0,Eweight*W_phi1*W_phi2);}
                if(!isnan(correlator1)){        pParity_int_r_obs1->Fill(1,100*correlator1);
                                        pParity_int_r_obs2->Fill(1,100*correlator1,Eweight);}
                if(!isnan(correlator2)){        pParity_int_rr_obs1->Fill(1,100*correlator2);
                                        pParity_int_rr_obs2->Fill(1,100*correlator2,Eweight);}
                if(!isnan(correlator3)){        pParity_int_s_obs1->Fill(1,100*correlator3);
                                        pParity_int_s_obs2->Fill(1,100*correlator3,Eweight);}
                if(!isnan(correlator4)){        pParity_int_ss_obs1->Fill(1,100*correlator4);
                                        pParity_int_ss_obs2 -> Fill(1,100*correlator4,Eweight);
                                        pParity_int_ss_obs2_eff -> Fill(1,100*correlator4,Eweight/eff_comb);}
                if(!isnan(correlator5)){        pParity_int_ss_ran_obs1->Fill(1,100*correlator5);
                                        pParity_int_ss_ran_obs2->Fill(1,100*correlator5,Eweight);}
		if(!isnan(correlator4)){	pParity_eta_ss_obs1->Fill(1,0.5*(Eta+Eta2),100*correlator4);
                                        pParity_eta_ss_obs2->Fill(1,0.5*(Eta+Eta2),100*correlator4,Eweight);
                                        pParity_Deta_ss_obs1->Fill(1,fabs(Eta-Eta2),100*correlator4);
                                        pParity_Deta_ss_obs2->Fill(1,fabs(Eta-Eta2),100*correlator4,Eweight);
                                        pParity_pt_ss_obs1->Fill(1,0.5*(Pt+Pt2),100*correlator4);
                                        pParity_pt_ss_obs2->Fill(1,0.5*(Pt+Pt2),100*correlator4,Eweight);
                                        pParity_Dpt_ss_obs1->Fill(1,fabs(Pt-Pt2),100*correlator4);
                                        pParity_Dpt_ss_obs2->Fill(1,fabs(Pt-Pt2),100*correlator4,Eweight);}
                if(!isnan(correlator6)){        pParity_eta_ss_obs1->Fill(1+4,0.5*(Eta+Eta2),100*correlator6);
                                        pParity_eta_ss_obs2->Fill(1+4,0.5*(Eta+Eta2),100*correlator6,Eweight);
                                        pParity_Deta_ss_obs1->Fill(1+4,fabs(Eta-Eta2),100*correlator6);
                                        pParity_Deta_ss_obs2->Fill(1+4,fabs(Eta-Eta2),100*correlator6,Eweight);
                                        pParity_pt_ss_obs1->Fill(1+4,0.5*(Pt+Pt2),100*correlator6);
                                        pParity_pt_ss_obs2->Fill(1+4,0.5*(Pt+Pt2),100*correlator6,Eweight);
                                        pParity_Dpt_ss_obs1->Fill(1+4,fabs(Pt-Pt2),100*correlator6);
                                        pParity_Dpt_ss_obs2->Fill(1+4,fabs(Pt-Pt2),100*correlator6,Eweight);}
                if(!isnan(correlator7)){        pParity_eta_ss_obs1->Fill(1+8,0.5*(Eta+Eta2),100*correlator7);
                                        pParity_eta_ss_obs2->Fill(1+8,0.5*(Eta+Eta2),100*correlator7,Eweight);
                                        pParity_Deta_ss_obs1->Fill(1+8,fabs(Eta-Eta2),100*correlator7);
                                        pParity_Deta_ss_obs2->Fill(1+8,fabs(Eta-Eta2),100*correlator7,Eweight);
                                        pParity_pt_ss_obs1->Fill(1+8,0.5*(Pt+Pt2),100*correlator7);
                                        pParity_pt_ss_obs2->Fill(1+8,0.5*(Pt+Pt2),100*correlator7,Eweight);
                                        pParity_Dpt_ss_obs1->Fill(1+8,fabs(Pt-Pt2),100*correlator7);
                                        pParity_Dpt_ss_obs2->Fill(1+8,fabs(Pt-Pt2),100*correlator7,Eweight);}
		if(!isnan(correlator4)){	pParity_Q_ss_obs1->Fill(1,q_inv,100*correlator4);
                                        pParity_Q_ss_obs2->Fill(1,q_inv,100*correlator4,Eweight);}
				}
                                if(Charge<0) {
                if(!isnan(correlator0)){        pParity_int_obs1->Fill(2,100*correlator0);
                                        pParity_int_obs2->Fill(2,100*correlator0,Eweight);
                                        pParity_int_w_obs1->Fill(2,100*correlator0,W_phi1);
                                        pParity_int_w_obs2->Fill(2,100*correlator0,Eweight*W_phi1);
                                        pParity_int_ww_obs1->Fill(2,100*correlator0,W_phi1*W_phi2);
                                        pParity_int_ww_obs2->Fill(2,100*correlator0,Eweight*W_phi1*W_phi2);}
                if(!isnan(correlator1)){        pParity_int_r_obs1->Fill(2,100*correlator1);
                                        pParity_int_r_obs2->Fill(2,100*correlator1,Eweight);}
                if(!isnan(correlator2)){        pParity_int_rr_obs1->Fill(2,100*correlator2);
                                        pParity_int_rr_obs2->Fill(2,100*correlator2,Eweight);}
                if(!isnan(correlator3)){        pParity_int_s_obs1->Fill(2,100*correlator3);
                                        pParity_int_s_obs2->Fill(2,100*correlator3,Eweight);}
                if(!isnan(correlator4)){        pParity_int_ss_obs1->Fill(2,100*correlator4);
                                        pParity_int_ss_obs2->Fill(2,100*correlator4,Eweight);
                                        pParity_int_ss_obs2_eff -> Fill(2,100*correlator4,Eweight/eff_comb);}
                if(!isnan(correlator5)){        pParity_int_ss_ran_obs1->Fill(2,100*correlator5);
                                        pParity_int_ss_ran_obs2->Fill(2,100*correlator5,Eweight);}
                if(!isnan(correlator4)){        pParity_eta_ss_obs1->Fill(2,0.5*(Eta+Eta2),100*correlator4);
                                        pParity_eta_ss_obs2->Fill(2,0.5*(Eta+Eta2),100*correlator4,Eweight);
                                        pParity_Deta_ss_obs1->Fill(2,fabs(Eta-Eta2),100*correlator4);
                                        pParity_Deta_ss_obs2->Fill(2,fabs(Eta-Eta2),100*correlator4,Eweight);
                                        pParity_pt_ss_obs1->Fill(2,0.5*(Pt+Pt2),100*correlator4);
                                        pParity_pt_ss_obs2->Fill(2,0.5*(Pt+Pt2),100*correlator4,Eweight);
                                        pParity_Dpt_ss_obs1->Fill(2,fabs(Pt-Pt2),100*correlator4);
                                        pParity_Dpt_ss_obs2->Fill(2,fabs(Pt-Pt2),100*correlator4,Eweight);}
                if(!isnan(correlator6)){        pParity_eta_ss_obs1->Fill(2+4,0.5*(Eta+Eta2),100*correlator6);
                                        pParity_eta_ss_obs2->Fill(2+4,0.5*(Eta+Eta2),100*correlator6,Eweight);
                                        pParity_Deta_ss_obs1->Fill(2+4,fabs(Eta-Eta2),100*correlator6);
                                        pParity_Deta_ss_obs2->Fill(2+4,fabs(Eta-Eta2),100*correlator6,Eweight);
                                        pParity_pt_ss_obs1->Fill(2+4,0.5*(Pt+Pt2),100*correlator6);
                                        pParity_pt_ss_obs2->Fill(2+4,0.5*(Pt+Pt2),100*correlator6,Eweight);
                                        pParity_Dpt_ss_obs1->Fill(2+4,fabs(Pt-Pt2),100*correlator6);
                                        pParity_Dpt_ss_obs2->Fill(2+4,fabs(Pt-Pt2),100*correlator6,Eweight);}
                if(!isnan(correlator7)){        pParity_eta_ss_obs1->Fill(2+8,0.5*(Eta+Eta2),100*correlator7);
                                        pParity_eta_ss_obs2->Fill(2+8,0.5*(Eta+Eta2),100*correlator7,Eweight);
                                        pParity_Deta_ss_obs1->Fill(2+8,fabs(Eta-Eta2),100*correlator7);
                                        pParity_Deta_ss_obs2->Fill(2+8,fabs(Eta-Eta2),100*correlator7,Eweight);
                                        pParity_pt_ss_obs1->Fill(2+8,0.5*(Pt+Pt2),100*correlator7);
                                        pParity_pt_ss_obs2->Fill(2+8,0.5*(Pt+Pt2),100*correlator7,Eweight);
                                        pParity_Dpt_ss_obs1->Fill(2+8,fabs(Pt-Pt2),100*correlator7);
                                        pParity_Dpt_ss_obs2->Fill(2+8,fabs(Pt-Pt2),100*correlator7,Eweight);}
		if(!isnan(correlator4)){	pParity_Q_ss_obs1->Fill(2,q_inv,100*correlator4);
                                        pParity_Q_ss_obs2->Fill(2,q_inv,100*correlator4,Eweight);}
                                }
			} // 2nd track Lambda Candidates
}

//***************************Mark Here

		}//Second Loop	

		}  //Track

                float PP_in1  = 100.*(N_L_P1*(N_L_P1-1) + N_R_P1*(N_R_P1-1) -2*N_L_P1*N_R_P1)/(N_L_P1+N_R_P1)/(N_L_P1+N_R_P1-1);
                float PP_out1 = 100.*(N_T_P1*(N_T_P1-1) + N_B_P1*(N_B_P1-1) -2*N_T_P1*N_B_P1)/(N_T_P1+N_B_P1)/(N_T_P1+N_B_P1-1);
                float NN_in1  = 100.*(N_L_N1*(N_L_N1-1) + N_R_N1*(N_R_N1-1) -2*N_L_N1*N_R_N1)/(N_L_N1+N_R_N1)/(N_L_N1+N_R_N1-1);
                float NN_out1 = 100.*(N_T_N1*(N_T_N1-1) + N_B_N1*(N_B_N1-1) -2*N_T_N1*N_B_N1)/(N_T_N1+N_B_N1)/(N_T_N1+N_B_N1-1);
                float PN_in1  = 100.*(N_L_P1*N_L_N1 + N_R_P1*N_R_N1 - N_L_P1*N_R_N1 - N_L_N1*N_R_P1)/(N_L_P1+N_R_P1)/(N_L_N1+N_R_N1);
                float PN_out1 = 100.*(N_T_P1*N_T_N1 + N_B_P1*N_B_N1 - N_T_P1*N_B_N1 - N_T_N1*N_B_P1)/(N_T_P1+N_B_P1)/(N_T_N1+N_B_N1);
                float NP_in1  = 100.*(N_L_N1*N_L_P1 + N_R_N1*N_R_P1 - N_L_N1*N_R_P1 - N_L_P1*N_R_N1)/(N_L_N1+N_R_N1)/(N_L_P1+N_R_P1);
                float NP_out1 = 100.*(N_T_N1*N_T_P1 + N_B_N1*N_B_P1 - N_T_N1*N_B_P1 - N_T_P1*N_B_N1)/(N_T_N1+N_B_N1)/(N_T_P1+N_B_P1);
                float PP_in2  = 100.*(N_L_P2*(N_L_P2-1) + N_R_P2*(N_R_P2-1) -2*N_L_P2*N_R_P2)/(N_L_P2+N_R_P2)/(N_L_P2+N_R_P2-1);
                float PP_out2 = 100.*(N_T_P2*(N_T_P2-1) + N_B_P2*(N_B_P2-1) -2*N_T_P2*N_B_P2)/(N_T_P2+N_B_P2)/(N_T_P2+N_B_P2-1);
                float NN_in2  = 100.*(N_L_N2*(N_L_N2-1) + N_R_N2*(N_R_N2-1) -2*N_L_N2*N_R_N2)/(N_L_N2+N_R_N2)/(N_L_N2+N_R_N2-1);
                float NN_out2 = 100.*(N_T_N2*(N_T_N2-1) + N_B_N2*(N_B_N2-1) -2*N_T_N2*N_B_N2)/(N_T_N2+N_B_N2)/(N_T_N2+N_B_N2-1);
                float PN_in2  = 100.*(N_L_P2*N_L_N2 + N_R_P2*N_R_N2 - N_L_P2*N_R_N2 - N_L_N2*N_R_P2)/(N_L_P2+N_R_P2)/(N_L_N2+N_R_N2);
                float PN_out2 = 100.*(N_T_P2*N_T_N2 + N_B_P2*N_B_N2 - N_T_P2*N_B_N2 - N_T_N2*N_B_P2)/(N_T_P2+N_B_P2)/(N_T_N2+N_B_N2);
                float NP_in2  = 100.*(N_L_N2*N_L_P2 + N_R_N2*N_R_P2 - N_L_N2*N_R_P2 - N_L_P2*N_R_N2)/(N_L_N2+N_R_N2)/(N_L_P2+N_R_P2);
                float NP_out2 = 100.*(N_T_N2*N_T_P2 + N_B_N2*N_B_P2 - N_T_N2*N_B_P2 - N_T_P2*N_B_N2)/(N_T_N2+N_B_N2)/(N_T_P2+N_B_P2);

                int Nin1  = (N_L_P1-N_L_N1)-(N_R_P1-N_R_N1);
                int Nout1 = (N_T_P1-N_T_N1)-(N_B_P1-N_B_N1);
                int Nin2  = (N_L_P2-N_L_N2)-(N_R_P2-N_R_N2);
                int Nout2 = (N_T_P2-N_T_N2)-(N_B_P2-N_B_N2);
                Hist_dQ_in1->Fill(Nin1,Eweight);
                Hist_dQ_out1->Fill(Nout1,Eweight);
                Hist_dQ_in2->Fill(Nin2,Eweight);
                Hist_dQ_out2->Fill(Nout2,Eweight);

                                if(!isnan(PP_in1)) pParity_int_MSC_same_in1->Fill(Nin1,MM*PP_in1,Eweight);
                                if(!isnan(NN_in1)) pParity_int_MSC_same_in1->Fill(Nin1,MM*NN_in1,Eweight);
                                if(!isnan(PN_in1)) pParity_int_MSC_oppo_in1->Fill(Nin1,MM*PN_in1,Eweight);
                                if(!isnan(NP_in1)) pParity_int_MSC_oppo_in1->Fill(Nin1,MM*NP_in1,Eweight);
                                if(!isnan(PP_out1)) pParity_int_MSC_same_out1->Fill(Nout1,MM*PP_out1,Eweight);
                                if(!isnan(NN_out1)) pParity_int_MSC_same_out1->Fill(Nout1,MM*NN_out1,Eweight);
                                if(!isnan(PN_out1)) pParity_int_MSC_oppo_out1->Fill(Nout1,MM*PN_out1,Eweight);
                                if(!isnan(NP_out1)) pParity_int_MSC_oppo_out1->Fill(Nout1,MM*NP_out1,Eweight);
                                if(!isnan(PP_in2)) pParity_int_MSC_same_in2->Fill(Nin2,MM*PP_in2,Eweight);
                                if(!isnan(NN_in2)) pParity_int_MSC_same_in2->Fill(Nin2,MM*NN_in2,Eweight);
                                if(!isnan(PN_in2)) pParity_int_MSC_oppo_in2->Fill(Nin2,MM*PN_in2,Eweight);
                                if(!isnan(NP_in2)) pParity_int_MSC_oppo_in2->Fill(Nin2,MM*NP_in2,Eweight);
                                if(!isnan(PP_out2)) pParity_int_MSC_same_out2->Fill(Nout2,MM*PP_out2,Eweight);
                                if(!isnan(NN_out2)) pParity_int_MSC_same_out2->Fill(Nout2,MM*NN_out2,Eweight);
                                if(!isnan(PN_out2)) pParity_int_MSC_oppo_out2->Fill(Nout2,MM*PN_out2,Eweight);
                                if(!isnan(NP_out2)) pParity_int_MSC_oppo_out2->Fill(Nout2,MM*NP_out2,Eweight);

                if(!isnan(PP_in1) && !isnan(PP_out1)) pParity_int_MSC1_obs1->Fill(1, MM*(PP_in1 - PP_out1));
                if(!isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs1->Fill(2, MM*(NN_in1 - NN_out1));
                if(!isnan(PP_in1) && !isnan(PP_out1) && !isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs1->Fill(3, MM*(PP_in1 - PP_out1));
                if(!isnan(PP_in1) && !isnan(PP_out1) && !isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs1->Fill(3, MM*(NN_in1 - NN_out1));
                if(!isnan(PN_in1) && !isnan(PN_out1) && !isnan(NP_in1) && !isnan(NP_out1)) pParity_int_MSC1_obs1->Fill(4, MM*(PN_in1 - PN_out1));
                if(!isnan(PN_in1) && !isnan(PN_out1) && !isnan(NP_in1) && !isnan(NP_out1)) pParity_int_MSC1_obs1->Fill(4, MM*(NP_in1 - NP_out1));
                if(!isnan(PP_in1) && !isnan(PP_out1)) pParity_int_MSC1_obs2->Fill(1, MM*(PP_in1 - PP_out1),Eweight);
                if(!isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs2->Fill(2, MM*(NN_in1 - NN_out1),Eweight);
                if(!isnan(PP_in1) && !isnan(PP_out1) && !isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs2->Fill(3, MM*(PP_in1 - PP_out1),Eweight);
                if(!isnan(PP_in1) && !isnan(PP_out1) && !isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs2->Fill(3, MM*(NN_in1 - NN_out1),Eweight);
                if(!isnan(PN_in1) && !isnan(PN_out1) && !isnan(NP_in1) && !isnan(NP_out1)) pParity_int_MSC1_obs2->Fill(4, MM*(PN_in1 - PN_out1),Eweight);
                if(!isnan(PN_in1) && !isnan(PN_out1) && !isnan(NP_in1) && !isnan(NP_out1)) pParity_int_MSC1_obs2->Fill(4, MM*(NP_in1 - NP_out1),Eweight);
                if(!isnan(PP_in2) && !isnan(PP_out2)) pParity_int_MSC2_obs1->Fill(1, MM*(PP_in2 - PP_out2));
                if(!isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs1->Fill(2, MM*(NN_in2 - NN_out2));
                if(!isnan(PP_in2) && !isnan(PP_out2) && !isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs1->Fill(3, MM*(PP_in2 - PP_out2));
                if(!isnan(PP_in2) && !isnan(PP_out2) && !isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs1->Fill(3, MM*(NN_in2 - NN_out2));
                if(!isnan(PN_in2) && !isnan(PN_out2) && !isnan(NP_in2) && !isnan(NP_out2)) pParity_int_MSC2_obs1->Fill(4, MM*(PN_in2 - PN_out2));
                if(!isnan(PN_in2) && !isnan(PN_out2) && !isnan(NP_in2) && !isnan(NP_out2)) pParity_int_MSC2_obs1->Fill(4, MM*(NP_in2 - NP_out2));
                if(!isnan(PP_in2) && !isnan(PP_out2)) pParity_int_MSC2_obs2->Fill(1, MM*(PP_in2 - PP_out2),Eweight);
                if(!isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs2->Fill(2, MM*(NN_in2 - NN_out2),Eweight);
                if(!isnan(PP_in2) && !isnan(PP_out2) && !isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs2->Fill(3, MM*(PP_in2 - PP_out2),Eweight);
                if(!isnan(PP_in2) && !isnan(PP_out2) && !isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs2->Fill(3, MM*(NN_in2 - NN_out2),Eweight);
                if(!isnan(PN_in2) && !isnan(PN_out2) && !isnan(NP_in2) && !isnan(NP_out2)) pParity_int_MSC2_obs2->Fill(4, MM*(PN_in2 - PN_out2),Eweight);
                if(!isnan(PN_in2) && !isnan(PN_out2) && !isnan(NP_in2) && !isnan(NP_out2)) pParity_int_MSC2_obs2->Fill(4, MM*(NP_in2 - NP_out2),Eweight);
/*
		if((N_L_N - N_L_P) == (N_R_N - N_R_P)) {
                if(!isnan(PP_in)) pParity_int_MSC_dQ_in_obs1->Fill(1, MM*PP_in);
                if(!isnan(NN_in)) pParity_int_MSC_dQ_in_obs1->Fill(2, MM*NN_in);
                if(!isnan(PP_in) && !isnan(NN_in)) pParity_int_MSC_dQ_in_obs1->Fill(3, MM*PP_in);
                if(!isnan(PP_in) && !isnan(NN_in)) pParity_int_MSC_dQ_in_obs1->Fill(3, MM*NN_in);
                if(!isnan(PN_in) && !isnan(NP_in)) pParity_int_MSC_dQ_in_obs1->Fill(4, MM*PN_in);
                if(!isnan(PN_in) && !isnan(NP_in)) pParity_int_MSC_dQ_in_obs1->Fill(4, MM*NP_in);
                if(!isnan(PP_in)) pParity_int_MSC_dQ_in_obs2->Fill(1, MM*PP_in,Eweight);
                if(!isnan(NN_in)) pParity_int_MSC_dQ_in_obs2->Fill(2, MM*NN_in,Eweight);
                if(!isnan(PP_in) && !isnan(NN_in)) pParity_int_MSC_dQ_in_obs2->Fill(3, MM*PP_in,Eweight);
                if(!isnan(PP_in) && !isnan(NN_in)) pParity_int_MSC_dQ_in_obs2->Fill(3, MM*NN_in,Eweight);
                if(!isnan(PN_in) && !isnan(NP_in)) pParity_int_MSC_dQ_in_obs2->Fill(4, MM*PN_in,Eweight);
                if(!isnan(PN_in) && !isnan(NP_in)) pParity_int_MSC_dQ_in_obs2->Fill(4, MM*NP_in,Eweight);
		}
                if((N_T_N - N_T_P) == (N_B_N - N_B_P)) {
                if(!isnan(PP_out)) pParity_int_MSC_dQ_out_obs1->Fill(1, MM*PP_out);
                if(!isnan(NN_out)) pParity_int_MSC_dQ_out_obs1->Fill(2, MM*NN_out);
                if(!isnan(PP_out) && !isnan(NN_out)) pParity_int_MSC_dQ_out_obs1->Fill(3, MM*PP_out);
                if(!isnan(PP_out) && !isnan(NN_out)) pParity_int_MSC_dQ_out_obs1->Fill(3, MM*NN_out);
                if(!isnan(PN_out) && !isnan(NP_out)) pParity_int_MSC_dQ_out_obs1->Fill(4, MM*PN_out);
                if(!isnan(PN_out) && !isnan(NP_out)) pParity_int_MSC_dQ_out_obs1->Fill(4, MM*NP_out);
                if(!isnan(PP_out)) pParity_int_MSC_dQ_out_obs2->Fill(1, MM*PP_out,Eweight);
                if(!isnan(NN_out)) pParity_int_MSC_dQ_out_obs2->Fill(2, MM*NN_out,Eweight);
                if(!isnan(PP_out) && !isnan(NN_out)) pParity_int_MSC_dQ_out_obs2->Fill(3, MM*PP_out,Eweight);
                if(!isnan(PP_out) && !isnan(NN_out)) pParity_int_MSC_dQ_out_obs2->Fill(3, MM*NN_out,Eweight);
                if(!isnan(PN_out) && !isnan(NP_out)) pParity_int_MSC_dQ_out_obs2->Fill(4, MM*PN_out,Eweight);
                if(!isnan(PN_out) && !isnan(NP_out)) pParity_int_MSC_dQ_out_obs2->Fill(4, MM*NP_out,Eweight);
                }
*/

//**********************************************************************************************************************
//Lambda
if(correction){
                int NLambda= (int)leaf_NoLambda->GetValue(0);
                for(int trki = 0; trki < NLambda; trki++){
                        double Pt_Lambda = leaf_PtW0->GetValue(trki);
			double Px_Lambda = leaf_PxW0->GetValue(trki);
			double Py_Lambda = leaf_PyW0->GetValue(trki);
			double Mass_Lambda = leaf_LambdaMass->GetValue(trki);
			double Track1_nSigma = leaf_Dau1nSigma->GetValue(trki);
			double Track2_nSigma = leaf_Dau2nSigma->GetValue(trki);
//			if(fabs(Mass_Lambda - 1.115683)>0.004)continue;

                        Eta       = leaf_EtaW0->GetValue(trki);
                        Theta     = 2.*atan(exp(-Eta));
			TVector2 mLambdaPhi;
			mLambdaPhi.Set(Px_Lambda,Py_Lambda);
                        Phi       = mLambdaPhi.Phi();
                        if(Phi> PI) Phi -= 2*PI;
                        if(Phi<-PI) Phi += 2*PI;
                        DCAGlobal = leaf_DCAglobalW0->GetValue(trki);
			if(fabs(Track1_nSigma)>3) continue;
			if(fabs(Track2_nSigma)>3) continue;
                	if((leaf_DecayLength->GetValue(trki))<6)continue;
                	if((DCAGlobal)>0.6)continue;
                	if((leaf_Dau1Dca->GetValue(trki))<0.6)continue;
                	if((leaf_Dau2Dca->GetValue(trki))<1.8)continue;
                	if((leaf_Dca1to2->GetValue(trki))>0.7)continue;

                	if(Pt_Lambda<0.6&&(leaf_Dau1Dca->GetValue(trki))<0.7)continue;
                	if(Pt_Lambda<0.6&&(leaf_Dau2Dca->GetValue(trki))<2.5)continue;

                        if(DCAGlobal > DcaCut) continue;
                        hEtaPtDist->Fill(Eta,Pt,Eweight);
                        if(Pt_Lambda < 0.5 || Pt_Lambda > pt_trig_up) continue;
                        if(Eta > EtaCut || Eta < -EtaCut) continue;
			hmInvMass->Fill(Mass_Lambda);
                        if(fabs(Mass_Lambda - 1.115683)>0.004)continue;

                        if(fWgt->IsOpen() && (TPCmean_FF_Lambda->GetEntries() || TPCmean_RF_Lambda->GetEntries())) {
                                for(int k=1;k<5;k++) {
                                        for(int kk=0;kk<order;kk++) {
                                          if(Eta>0 && PVtxz>0) PhiMean_Lambda[k-1+4*kk] = TPCmean_FF_Lambda->GetBinContent(k+8*kk,Day2-7999);
                                          if(Eta<0 && PVtxz>0) PhiMean_Lambda[k-1+4*kk] = TPCmean_RF_Lambda->GetBinContent(k+8*kk,Day2-7999);
                                          if(Eta>0 && PVtxz<0) PhiMean_Lambda[k-1+4*kk] = TPCmean_FF_Lambda->GetBinContent(k+4+8*kk,Day2-7999);
                                          if(Eta<0 && PVtxz<0) PhiMean_Lambda[k-1+4*kk] = TPCmean_RF_Lambda->GetBinContent(k+4+8*kk,Day2-7999);
                                        }
                                }
                        }
                        float cos1 =0, cos2=0, sin1=0, sin2=0, Phi_new = Phi;
                        // Recentering parameters
                        if(Eta>0 && PVtxz>0) {
                        	{
                                        for(int kk=0;kk<order;kk++) {
                                                pTPCmeanPhi_FF_Lambda->Fill(1+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
                                                pTPCmeanPhi_FF_Lambda->Fill(2+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
                                        }
                                        Hist_Phi_FF_Lambda->Fill(Phi,1,Eweight);}
                        }
                        if(Eta<0 && PVtxz>0) {
                        	{
                                        for(int kk=0;kk<order;kk++) {
                                                pTPCmeanPhi_RF_Lambda->Fill(1+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
                                                pTPCmeanPhi_RF_Lambda->Fill(2+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
                                        }
                                        Hist_Phi_RF_Lambda->Fill(Phi,1,Eweight);}
                        }
                        if(Eta>0 && PVtxz<0) {
                        	 {
                                        for(int kk=0;kk<order;kk++) {
                                                pTPCmeanPhi_FF_Lambda->Fill(1+4+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
                                                pTPCmeanPhi_FF_Lambda->Fill(2+4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
                                        }
                                        Hist_Phi_FF_Lambda->Fill(Phi,1+2,Eweight);}
                        }
                        if(Eta<0 && PVtxz<0) {
                        
                                        for(int kk=0;kk<order;kk++) {
                                                pTPCmeanPhi_RF_Lambda->Fill(1+4+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
                                                pTPCmeanPhi_RF_Lambda->Fill(2+4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
                                        }
                                        Hist_Phi_RF_Lambda->Fill(Phi,1+2,Eweight);
                        }
				{
                        		cos1 = cos(Phi) - PhiMean_Lambda[0]; sin1 = sin(Phi) - PhiMean_Lambda[1];
                                        double a2np = 1+PhiMean_Lambda[0+4*(order-1)], a2nn = 1-PhiMean_Lambda[0+4*(order-1)];
                                        double lambda_2nsp = PhiMean_Lambda[1+4*(order-1)]/a2np;
                                        double lambda_2nsn = PhiMean_Lambda[1+4*(order-1)]/a2nn;
                                        double cos11 = (cos1 - lambda_2nsn*sin1)/(1-lambda_2nsn*lambda_2nsp);
                                        double sin11 = (sin1 - lambda_2nsp*cos1)/(1-lambda_2nsn*lambda_2nsp);
                                        cos1 = cos11/a2np;
                                        sin1 = sin11/a2nn;
                                        for(int jj=0;jj<order;jj++) Phi_new += -2*PhiMean_Lambda[1+4*jj]*cos(jj*Phi+Phi)/int(jj+1)
                                                                               +2*PhiMean_Lambda[0+4*jj]*sin(jj*Phi+Phi)/int(jj+1);
                                        if(Phi_new> PI) Phi_new -= 2*PI;
                                        if(Phi_new<-PI) Phi_new += 2*PI;
                                        if(Eta>0 && PVtxz>0) Hist_Phi_FF_Lambda_new->Fill(Phi_new,1,Eweight);
                                        if(Eta<0 && PVtxz>0) Hist_Phi_RF_Lambda_new->Fill(Phi_new,1,Eweight);
                                        if(Eta>0 && PVtxz<0) Hist_Phi_FF_Lambda_new->Fill(Phi_new,1+2,Eweight);
                                        if(Eta<0 && PVtxz<0) Hist_Phi_RF_Lambda_new->Fill(Phi_new,1+2,Eweight);}
		}//Lambda
//************************************************************************************************************************************
}
        } // Event

	fout.Write();
char fname_new[200];
sprintf(fname_new,"cen%d.weight_new.root",cen);
TFile *fWgtNew = new TFile(fname_new,"RECREATE");
pTPCmeanPhi_FF->Write();
pTPCmeanPhi_RF->Write();
Hist_Phi_FF->Write();
Hist_Phi_RF->Write();
pTPC_EP_east->Write();
pTPC_EP_west->Write();
pTPC_EP_full->Write();

pTPCmeanPhi_FF_Lambda->Write();
pTPCmeanPhi_RF_Lambda->Write();
pTPCmeanPhi_FF_Proton->Write();
pTPCmeanPhi_RF_Proton->Write();
Hist_Phi_FF_Lambda->Write();
Hist_Phi_RF_Lambda->Write();
Hist_Phi_FF_Proton->Write();
Hist_Phi_RF_Proton->Write();
fWgtNew->Close();
	return;
}

