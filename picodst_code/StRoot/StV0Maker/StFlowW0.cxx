#include "StFlowW0.h"

#include "TObject.h"

ClassImp(StFlowW0)

StFlowW0::StFlowW0(StFlowW0* flowW0) : TObject ()
{
        mMassW0 = flowW0->massW0();
        mRapidityW0 = flowW0->rapidityW0();
        mEtaW0 = flowW0->etaW0();
        mPtglobalW0 = flowW0->ptglobalW0();
        mPxglobalW0 = flowW0->pxglobalW0();
        mPyglobalW0 = flowW0->pyglobalW0();
        mPzglobalW0 = flowW0->pzglobalW0();
        mOxglobalW0 = flowW0->oxglobalW0();
        mOyglobalW0 = flowW0->oyglobalW0();
        mOzglobalW0 = flowW0->ozglobalW0();
        mDecaylenW0 = flowW0->decaylenW0();
        mDcaW0 = flowW0->dcaW0();
        mDca2DW0 = flowW0->dca2DW0();
        mPathlenW0 = flowW0->pathlenW0();

        mDau1idW0 = flowW0->dau1idW0();
        mDau1dcaW0 = flowW0->dau1dcaW0();
        mDau1dca2DW0 = flowW0->dau1dca2DW0();
        mDau1nhitsW0 = flowW0->dau1nhitsW0();
        mDau1dedxW0 = flowW0->dau1dedxW0();
        mDau1nSigmaW0 = flowW0->dau1nSigmaW0();
        mDau1etaW0 = flowW0->dau1etaW0();
        mDau1ptW0 = flowW0->dau1ptW0();
        mDau1pxW0 = flowW0->dau1pxW0();
        mDau1pyW0 = flowW0->dau1pyW0();
        mDau1pzW0 = flowW0->dau1pzW0();
        mDau1tofflagW0 = flowW0->dau1tofflagW0();
        mDau1tofW0 = flowW0->dau1tofW0();
        mDau1pathlenW0 = flowW0->dau1pathlenW0();
        mDau1pt_prW0 = flowW0->dau1pt_prW0();
        mDau1px_prW0 = flowW0->dau1px_prW0();
        mDau1py_prW0 = flowW0->dau1py_prW0();
        mDau1pz_prW0 = flowW0->dau1pz_prW0();

        mDca1to2W0 = flowW0->dca1to2W0();

        mDau2idW0 = flowW0->dau2idW0();
        mDau2dcaW0 = flowW0->dau2dcaW0();
        mDau2dca2DW0 = flowW0->dau2dca2DW0();
        mDau2nhitsW0 = flowW0->dau2nhitsW0();
        mDau2dedxW0 = flowW0->dau2dedxW0();
        mDau2nSigmaW0 = flowW0->dau2nSigmaW0();
        mDau2etaW0 = flowW0->dau2etaW0();
        mDau2ptW0 = flowW0->dau2ptW0();
        mDau2pxW0 = flowW0->dau2pxW0();
        mDau2pyW0 = flowW0->dau2pyW0();
        mDau2pzW0 = flowW0->dau2pzW0();
        mDau2tofflagW0 = flowW0->dau2tofflagW0();
        mDau2tofW0 = flowW0->dau2tofW0();
        mDau2pathlenW0 = flowW0->dau2pathlenW0();
        mDau2pt_prW0 = flowW0->dau2pt_prW0();
        mDau2px_prW0 = flowW0->dau2px_prW0();
        mDau2py_prW0 = flowW0->dau2py_prW0();
        mDau2pz_prW0 = flowW0->dau2pz_prW0();
}
