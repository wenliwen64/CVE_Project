#include "TLeaf.h"
#include "TH1F.h"
#include "TChain.h"
#include <iostream>
#include <stdio.h> 
#include "TFile.h"
void cal_eff_protonII(Int_t cen){
    //define the centreality bin based on the multiplicity
    const int cenDef[9] = {10, 22, 43, 76, 125, 193, 281, 396, 466};

    TChain* chain = new TChain("StrangenessDst");
    chain -> Add("/media/PingYuan/Run11_200GeV/Data154/*.lambda.picodst.root");
    Int_t nEntries = chain -> GetEntries();
    cout<<nEntries<<" events added"<<endl;

    char h_tpc_nm[100];
    char h_tpctof_nm[100];
    char outfile_nm[100];
    sprintf(h_tpc_nm, "h_tpc_cen%d", cen);
    sprintf(h_tpctof_nm, "h_tpctof_cen%d", cen);
    sprintf(outfile_nm, "eff_proton_cen%d", cen);
   
    TFile* outfile = new TFile(outfile_nm, "recreate"); 
    TH1F* h_tpc = new TH1F(h_tpc_nm, h_tpctof_nm, 50, 0, 2.5);
    TH1F* h_tpctof = new TH1F(h_tpctof_nm, h_tpctof_nm, 50, 0, 2.5);

    for(int i = 0; i < nEntries; i++){
        if(i%10000 == 0) cout<<"Processing entry == "<<i<< "== out of "<<nEntries<<endl;
        chain -> GetEntry(i);

        //Eventwise
        TLeaf* leaf_RunId = chain -> GetLeaf("mRunId");
        TLeaf* leaf_EventId = chain -> GetLeaf("mEventId");
        TLeaf* leaf_Trigger = chain -> GetLeaf("mTrigger");
        TLeaf* leaf_Bz = chain -> GetLeaf("mBz");
        TLeaf* leaf_PrimaryVertexZ = chain -> GetLeaf("mPrimaryVertexZ");
        TLeaf* leaf_RefMult = chain -> GetLeaf("mRefMult");
        TLeaf* leaf_NoTracks = chain -> GetLeaf("mNoTracks");

        Int_t refmult = leaf_RefMult -> GetValue(0);
        Int_t cenbin = 0;

        for(int ii = 0; ii < 9; ii++){ 
            if(refmult > cenDef[ii]) cenbin++;
        }

        if(cenbin != cen) continue;
        //Trackwise
        TLeaf* leaf_PtV0 = chain -> GetLeaf("fV0s.mPtprimaryV0");
        TLeaf* leaf_EtaV0 = chain -> GetLeaf("fV0s.mEtaV0");
        TLeaf* leaf_PhiV = chain -> GetLeaf("fV0s.mEtaV");
        TLeaf* leaf_ChargeV0 = chain -> GetLeaf("fV0s.mChargeV0");
        TLeaf* leaf_DCAglobalV0 =  chain -> GetLeaf("fV0s.mDCAglobalV0");
        TLeaf* leaf_ndEdxV0 = chain -> GetLeaf("fV0s.mdEdxV0");
        TLeaf* leaf_nSigmaPV0 = chain -> GetLeaf("fV0s.mnSigmaPV0"); 
        TLeaf* leaf_nSigmaPiV0 = chain -> GetLeaf("fV0s.mnSigmaPiV0");
        TLeaf* leaf_nsigmaKV0 = chain -> GetLeaf("fV0s.mnSigmaKV0"); 
        TLeaf* leaf_TofflagV0 = chain -> GetLeaf("fV0s.mTofflagV0");
        TLeaf* leaf_TofV0 = chain -> GetLeaf("fV0s.mTofV0");
        TLeaf* leaf_PathlenV0 = chain -> GetLeaf("fV0s.mPathlenV0");
        TLeaf* leaf_PV0 = chain -> GetLeaf("fV0s.mPV0");
        TLeaf* leaf_trackIdV0 = chain -> GetLeaf("fV0s.mTrackIdV0");

        Int_t no_tracks = leaf_NoTracks -> GetValue(0);

        for(int j = 0; j < no_tracks; j++){
            Float_t DCAGlobal = leaf_DCAglobalV0 -> GetValue(j);
            Float_t Eta = leaf_EtaV0 -> GetValue(j);
            Float_t Pt = leaf_PtV0 -> GetValue(j);
            Float_t nSigmaP = leaf_nSigmaPV0 -> GetValue(j);

            Float_t proton_p = leaf_PV0 -> GetValue(j); 
            Int_t TofFlag = leaf_TofflagV0 -> GetValue(j);
            Float_t TOF_proton = leaf_TofV0 -> GetValue(j);
            Float_t TOF_path = leaf_PathlenV0 -> GetValue(j);
          
	    if(abs(nSigmaP) > 2.0) continue;
	    if(abs(Eta) > 1.0) continue;
            if(DCAGlobal > 2.0) continue;
            if(Pt < 0.15 || Pt > 2.0) continue;
//TODO: add some dedx check plots to see if it is correct
            h_tpc -> Fill(Pt);
            
            Float_t mass2proton = 0;
            if(TofFlag > 0){
		mass2proton = proton_p*proton_p*(900.0*TOF_proton*TOF_proton/TOF_path/TOF_path - 1.0); 
            }
            if(mass2proton > 0.8 && mass2proton < 1.0){
                h_tpctof -> Fill(Pt);
            }
	}
    }
    outfile -> Write();
    outfile -> Close();
}
