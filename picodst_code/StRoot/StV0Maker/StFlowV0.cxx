#include "StFlowV0.h"

#include "TObject.h"

ClassImp(StFlowV0)

StFlowV0::StFlowV0(StFlowV0* flowV0) : TObject ()
{
	mTrackTypeV0 = flowV0->tracktypeV0();
        mTrackFlagV0 = flowV0->trackflagV0();
        mTrackIdV0 = flowV0->trackIdV0();
	mChi2zV0= flowV0->chi2zV0();
	mChi2V0 = flowV0->chi2V0();
	mPV0	= flowV0->pV0();        
	mEtaV0	= flowV0->etaV0();    
	mPhiV0  = flowV0->phiV0();
	mChargeV0=flowV0->chargeV0();
	mDCAglobalV0=flowV0->dcaglobalV0();
	mdEdxV0 = flowV0->dedxV0(); 
        mnSigmaPiV0 = flowV0->nSigmaPiV0();
        mnSigmaKV0 = flowV0->nSigmaKV0();
        mnSigmaPV0 = flowV0->nSigmaPV0();
	mnSigmaEV0 = flowV0->nSigmaEV0();
        mPtprimaryV0= flowV0->ptprimaryV0();
	mPxprimaryV0= flowV0->pxprimaryV0();
        mPyprimaryV0= flowV0->pyprimaryV0();
        mPzprimaryV0= flowV0->pzprimaryV0();
	mTofflagV0= flowV0->tofflagV0();
	mTofV0= flowV0->tofV0();
	mPathlenV0= flowV0->pathlenV0();
}
