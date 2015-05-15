#ifndef StFlowV0_hh
#define StFlowV0_hh
#include "TObject.h"

class StFlowV0 : public TObject {
public:
	StFlowV0() { }
	StFlowV0(StFlowV0* flowV0);
	
	virtual       ~StFlowV0() { }
	
	Int_t	tracktypeV0()	const;	// track type
	Int_t	trackflagV0()	const;  // track flag
        Int_t   trackIdV0()   const;  // track ID
        Float_t chi2zV0()  	const;  // chi2z
	Float_t chi2V0()	const;  // chi2
	Float_t	pV0()		const;	// Total momentum
	Float_t	etaV0()		const;	// Pseudorapidity
	Float_t phiV0()		const;  // azimuth
	Int_t   chargeV0()	const;  // charge
	Float_t dcaglobalV0()	const;  // dca global
	Float_t dedxV0()	const;  // dedx
        Float_t nSigmaPiV0()     const;  // nSigma of electron
        Float_t nSigmaKV0()     const;  // nSigma of electron
        Float_t nSigmaPV0()     const;  // nSigma of electron
	Float_t nSigmaEV0()	const;  // nSigma of electron
        Float_t ptprimaryV0()    const;  // pt of primary track
	Float_t pxprimaryV0()	const;  // px of primary track
        Float_t pyprimaryV0()    const;  // py of primary track
        Float_t pzprimaryV0()    const;  // pz of primary track
	Int_t	tofflagV0()	const;	// tof flag of primary track
	Float_t tofV0()		const;	//Time of flight of primary track
	Float_t pathlenV0()	const;	//pathlength of primary track

	void 	SetTrackTypeV0	( Int_t	  ttS)			{ mTrackTypeV0	= ttS;}
	void    SetTrackFlagV0  ( Int_t   ttS)                  { mTrackFlagV0  = ttS;}
        void    SetTrackIdV0  ( Int_t   pId)                  { mTrackIdV0  = pId;}
        void    SetChi2zV0	( Float_t chi2zS)		{ mChi2zV0	= chi2zS;}
	void 	SetChi2V0	( Float_t chi2S)		{ mChi2V0	= chi2S;}
	void	SetPV0		( Float_t ptS)			{ mPV0		= ptS; }
	void	SetEtaV0	( Float_t pseudorapS)		{ mEtaV0	= pseudorapS; }
        void    SetPhiV0	( Float_t phiS)			{ mPhiV0	= phiS; }
        void    SetChargeV0	( Int_t chargeS)		{ mChargeV0	= chargeS; }
	void 	SetDCAglobalV0  ( Float_t dcaG)			{ mDCAglobalV0	= dcaG;}
	void 	SetDedxV0	( Float_t dedxS)		{ mdEdxV0	= dedxS;}
	void	SetNSigmaPiV0	( Float_t nSigPi)		{ mnSigmaPiV0	= nSigPi;}
        void    SetNSigmaKV0    ( Float_t nSigK)                { mnSigmaKV0    = nSigK;}
        void    SetNSigmaPV0    ( Float_t nSigP)                { mnSigmaPV0    = nSigP;}
        void    SetNSigmaEV0    ( Float_t nSigE)                { mnSigmaEV0    = nSigE;}
        void    SetPtprimaryV0  ( Float_t ptG)                  { mPtprimaryV0  = ptG;}
	void	SetPxprimaryV0	( Float_t pG)			{ mPxprimaryV0	= pG;}
        void    SetPyprimaryV0   ( Float_t pG)                   { mPyprimaryV0   = pG;}
        void    SetPzprimaryV0   ( Float_t pG)                   { mPzprimaryV0   = pG;}
	void	SetTofflagV0	( Int_t	ttF)			{ mTofflagV0  = ttF;}
	void	SetTofV0	( Float_t tof)			{ mTofV0 = tof;}
	void	SetPathlenV0	( Float_t pathlength)		{ mPathlenV0 = pathlength;}
 
 private:
        Int_t   mTrackTypeV0;
        Int_t   mTrackFlagV0;
        Int_t   mTrackIdV0;
        Float_t mChi2zV0;
	Float_t mChi2V0;
	Float_t	mPV0;
	Float_t	mEtaV0;
        Float_t mPhiV0;
	Int_t   mChargeV0;
	Float_t mDCAglobalV0;
	Float_t mdEdxV0;
        Float_t mnSigmaPiV0;
        Float_t mnSigmaKV0;
        Float_t mnSigmaPV0;
	Float_t mnSigmaEV0;
        Float_t mPtprimaryV0;
	Float_t mPxprimaryV0;
        Float_t mPyprimaryV0;
        Float_t mPzprimaryV0;
	Int_t mTofflagV0;
	Float_t mTofV0;
	Float_t mPathlenV0;

ClassDef(StFlowV0,1)
};

inline          Int_t	       StFlowV0::tracktypeV0()	   const   { return mTrackTypeV0;}
inline          Int_t          StFlowV0::trackflagV0()     const   { return mTrackFlagV0;}
inline          Int_t          StFlowV0::trackIdV0()     const   { return mTrackIdV0;}
inline          Float_t        StFlowV0::chi2zV0()	   const   { return mChi2zV0;}
inline          Float_t        StFlowV0::chi2V0()	   const   { return mChi2V0;}
inline		Float_t	       StFlowV0::pV0()	           const   { return mPV0; }
inline		Float_t	       StFlowV0::etaV0()	   const   { return mEtaV0; }
inline          Float_t        StFlowV0::phiV0()           const   { return mPhiV0; }
inline          Int_t	       StFlowV0::chargeV0()	   const   { return mChargeV0; }
inline          Float_t        StFlowV0::dcaglobalV0()	   const   { return mDCAglobalV0;}
inline          Float_t        StFlowV0::dedxV0() 	   const   { return mdEdxV0;}
inline          Float_t        StFlowV0::nSigmaPiV0()       const   { return mnSigmaPiV0;}
inline          Float_t        StFlowV0::nSigmaKV0()       const   { return mnSigmaKV0;}
inline          Float_t        StFlowV0::nSigmaPV0()       const   { return mnSigmaPV0;}
inline          Float_t        StFlowV0::nSigmaEV0()	   const   { return mnSigmaEV0;}
inline          Float_t        StFlowV0::ptprimaryV0()     const   { return mPtprimaryV0;}
inline          Float_t        StFlowV0::pxprimaryV0()	   const   { return mPxprimaryV0;}
inline          Float_t        StFlowV0::pyprimaryV0()      const   { return mPyprimaryV0;}
inline          Float_t        StFlowV0::pzprimaryV0()      const   { return mPzprimaryV0;}
inline          Int_t          StFlowV0::tofflagV0()	   const    { return mTofflagV0;}
inline          Float_t        StFlowV0::tofV0()	   const    { return mTofV0;}
inline          Float_t        StFlowV0::pathlenV0()	   const    { return mPathlenV0;}
#endif
