#ifndef StFlowW0_hh
#define StFlowW0_hh
#include "TObject.h"

class StFlowW0 : public TObject {
public:
	StFlowW0() { }
	StFlowW0(StFlowW0* flowW0);
	
	virtual       ~StFlowW0() { }

	Float_t massW0()	const;	// InvMass
	Float_t	rapidityW0()	const;	// Rapidity
	Float_t etaW0()         const;  // Pseudorapidity
	Float_t ptglobalW0()    const;  // pt of global track
        Float_t pxglobalW0()    const;  // px of global track
        Float_t pyglobalW0()    const;  // py of global track
        Float_t pzglobalW0()    const;  // pz of global track
        Float_t oxglobalW0()    const;  // origin x of global track
        Float_t oyglobalW0()    const;  // origin y of global track
        Float_t ozglobalW0()    const;  // origin z of global track
	Float_t	decaylenW0()	const;	// decaylength
        Float_t dcaW0()         const;  // dca
        Float_t dca2DW0()	const;  // 2D dca
        Float_t pathlenW0()	const;  // pathlength

	Int_t	dau1idW0()	const;	// dau1id
	Float_t	dau1dcaW0()	const;	// dau1dca
	Float_t	dau1dca2DW0()	const;	// dau1dca2D
	Int_t	dau1nhitsW0()	const;	// dau1 nhits
	Float_t	dau1dedxW0()	const;	// dau1 dedx
	Float_t	dau1nSigmaW0()	const;	// dau1 nSigma
	Float_t	dau1etaW0()	const;	// dau1 eta
	Float_t dau1ptW0()	const;	// dau1 pt
	Float_t	dau1pxW0()	const;	// dau1 px
	Float_t	dau1pyW0()	const;	// dau1	py
	Float_t	dau1pzW0()	const;	// dau1	pz
	Int_t	dau1tofflagW0()	const;	// dau1 tof flag
	Float_t	dau1tofW0()	const;	// dau1 tof
	Float_t	dau1pathlenW0()	const;	// dau1 pathlength
	Int_t	dau1PrMatchW0()	const;	// dau1 primary Match
	Float_t	dau1pt_prW0()	const;	// dau1	primary pt
	Float_t	dau1px_prW0()	const;	// dau1 primary px
	Float_t	dau1py_prW0()	const;	// dau1 primary py
	Float_t	dau1pz_prW0()	const;	// dau1 primary pz

	Float_t dca1to2W0()	const;	// dca proton to pion

        Int_t   dau2idW0()      const;  // dau2id
        Float_t dau2dcaW0()     const;  // dau2dca
        Float_t dau2dca2DW0()   const;  // dau2dca2D
        Int_t   dau2nhitsW0()   const;  // dau2 nhits
        Float_t dau2dedxW0()    const;  // dau2 dedx
        Float_t dau2nSigmaW0()  const;  // dau2 nSigma
        Float_t dau2etaW0()     const;  // dau2 eta
        Float_t dau2ptW0()      const;  // dau2 pt
        Float_t dau2pxW0()      const;  // dau2 px
        Float_t dau2pyW0()      const;  // dau2 py
        Float_t dau2pzW0()      const;  // dau2 pz
        Int_t   dau2tofflagW0() const;  // dau2 tof flag
        Float_t dau2tofW0()     const;  // dau2 tof
        Float_t dau2pathlenW0() const;  // dau2 pathlength
        Int_t   dau2PrMatchW0() const;  // dau2 primary Match
        Float_t dau2pt_prW0()   const;  // dau2 primary pt
        Float_t dau2px_prW0()   const;  // dau2 primary px
        Float_t dau2py_prW0()   const;  // dau2 primary py
        Float_t dau2pz_prW0()   const;  // dau2 primary pz


        void	SetMassW0	( Float_t  massG)		{ mMassW0	= massG;}
        void	SetRapidityW0	( Float_t  rapidityG)		{ mRapidityW0	= rapidityG;}
        void 	SetEtaW0	( Float_t  etaG)		{ mEtaW0	= etaG;}
        void	SetPtglobalW0	( Float_t  ptG)			{ mPtglobalW0	= ptG;}
        void    SetPxglobalW0   ( Float_t pG)                   { mPxglobalW0   = pG;}
        void    SetPyglobalW0   ( Float_t pG)                   { mPyglobalW0   = pG;}
        void    SetPzglobalW0   ( Float_t pG)                   { mPzglobalW0   = pG;}
        void    SetOxglobalW0   ( Float_t oG)                   { mOxglobalW0   = oG;}
        void    SetOyglobalW0   ( Float_t oG)                   { mOyglobalW0   = oG;}
        void    SetOzglobalW0   ( Float_t oG)                   { mOzglobalW0   = oG;}
        void	SetDecaylenW0	( Float_t  declenG)		{ mDecaylenW0	= declenG;}
        void	SetDcaW0	( Float_t  dcaG)		{ mDcaW0	= dcaG;}
        void	SetDca2DW0	( Float_t  dca2DG)		{ mDca2DW0	= dca2DG;}
        void	SetPathlenW0	( Float_t  plenG)		{ mPathlenW0	= plenG;}

        void	SetDau1idW0	( Int_t	  ttS)			{ mDau1idW0	= ttS;}
        void	SetDau1dcaW0	( Float_t dcaG)			{ mDau1dcaW0	= dcaG;}
        void	SetDau1dca2DW0	( Float_t dcaG)			{ mDau1dca2DW0	= dcaG;}
        void	SetDau1nhitsW0	( Float_t nhitsG)		{ mDau1nhitsW0	= nhitsG;}
        void	SetDau1dedxW0	( Float_t dedxG)		{ mDau1dedxW0	= dedxG;}
        void	SetDau1nSigmaW0	( Float_t nSigmaG)		{ mDau1nSigmaW0	= nSigmaG;}
        void	SetDau1etaW0	( Float_t etaG)			{ mDau1etaW0	= etaG;}
        void	SetDau1ptW0	( Float_t ptG)			{ mDau1ptW0	= ptG;}
        void	SetDau1pxW0	( Float_t pG)			{ mDau1pxW0	= pG;}
        void	SetDau1pyW0	( Float_t pG)			{ mDau1pyW0	= pG;}
        void	SetDau1pzW0	( Float_t pG)			{ mDau1pzW0	= pG;}
        void	SetDau1tofflagW0( Int_t	  ttF)			{ mDau1tofflagW0= ttF;}
        void	SetDau1tofW0	( Float_t tof)			{ mDau1tofW0	= tof;}
        void	SetDau1pathlenW0( Float_t pathlength)		{ mDau1pathlenW0= pathlength;}
	void	SetDau1PrMatchW0( Int_t	  ttS)			{ mDau1PrMatchW0= ttS;}
        void	SetDau1pt_prW0	( Float_t ptS)			{ mDau1pt_prW0	= ptS;}
        void	SetDau1px_prW0	( Float_t pS)			{ mDau1px_prW0	= pS;}
        void	SetDau1py_prW0	( Float_t pS)			{ mDau1py_prW0	= pS;}
        void	SetDau1pz_prW0	( Float_t pS)			{ mDau1pz_prW0	= pS;}

	void	SetDca1to2W0	( Float_t dcaG)			{ mDca1to2W0	= dcaG;}

        void    SetDau2idW0     ( Int_t   ttS)                  { mDau2idW0     = ttS;}
        void    SetDau2dcaW0    ( Float_t dcaG)                 { mDau2dcaW0    = dcaG;}
        void    SetDau2dca2DW0  ( Float_t dcaG)                 { mDau2dca2DW0  = dcaG;}
        void    SetDau2nhitsW0  ( Float_t nhitsG)               { mDau2nhitsW0  = nhitsG;}
        void    SetDau2dedxW0   ( Float_t dedxG)                { mDau2dedxW0   = dedxG;}
        void    SetDau2nSigmaW0 ( Float_t nSigmaG)              { mDau2nSigmaW0 = nSigmaG;}
        void    SetDau2etaW0    ( Float_t etaG)                 { mDau2etaW0    = etaG;}
        void    SetDau2ptW0     ( Float_t ptG)                  { mDau2ptW0     = ptG;} 
        void    SetDau2pxW0     ( Float_t pG)                   { mDau2pxW0     = pG;}
        void    SetDau2pyW0     ( Float_t pG)                   { mDau2pyW0     = pG;}
        void    SetDau2pzW0     ( Float_t pG)                   { mDau2pzW0     = pG;}
        void    SetDau2tofflagW0( Int_t   ttF)                  { mDau2tofflagW0= ttF;}
        void    SetDau2tofW0    ( Float_t tof)                  { mDau2tofW0    = tof;}
        void    SetDau2pathlenW0( Float_t pathlength)           { mDau2pathlenW0= pathlength;}
        void    SetDau2PrMatchW0( Int_t   ttS)                  { mDau2PrMatchW0= ttS;}
        void    SetDau2pt_prW0  ( Float_t ptS)                  { mDau2pt_prW0  = ptS;}
        void    SetDau2px_prW0  ( Float_t pS)                   { mDau2px_prW0  = pS;}
        void    SetDau2py_prW0  ( Float_t pS)                   { mDau2py_prW0  = pS;}
        void    SetDau2pz_prW0  ( Float_t pS)                   { mDau2pz_prW0  = pS;}

 private:
        Float_t mMassW0;
        Float_t mRapidityW0;
        Float_t mEtaW0;
        Float_t mPtglobalW0;
        Float_t mPxglobalW0;
        Float_t mPyglobalW0;
        Float_t mPzglobalW0;
        Float_t mOxglobalW0;
        Float_t mOyglobalW0;
        Float_t mOzglobalW0;
        Float_t mDecaylenW0;
        Float_t mDcaW0;
        Float_t mDca2DW0;
        Float_t mPathlenW0;

        Int_t   mDau1idW0;
        Float_t mDau1dcaW0;
        Float_t mDau1dca2DW0;
        Float_t mDau1nhitsW0;
        Float_t mDau1dedxW0;
        Float_t mDau1nSigmaW0;
        Float_t mDau1etaW0;
        Float_t mDau1ptW0;
        Float_t mDau1pxW0;
        Float_t mDau1pyW0;
        Float_t mDau1pzW0;
        Int_t   mDau1tofflagW0;
        Float_t mDau1tofW0;
        Float_t mDau1pathlenW0;
	Int_t	mDau1PrMatchW0;
        Float_t mDau1pt_prW0;
        Float_t mDau1px_prW0;
        Float_t mDau1py_prW0;
        Float_t mDau1pz_prW0;

        Float_t mDca1to2W0;

        Int_t   mDau2idW0;
        Float_t mDau2dcaW0;
        Float_t mDau2dca2DW0;
        Float_t mDau2nhitsW0;
        Float_t mDau2dedxW0;
        Float_t mDau2nSigmaW0;
        Float_t mDau2etaW0;
        Float_t mDau2ptW0;
        Float_t mDau2pxW0;
        Float_t mDau2pyW0;
        Float_t mDau2pzW0;
        Int_t   mDau2tofflagW0;
        Float_t mDau2tofW0;
        Float_t mDau2pathlenW0;
	Int_t	mDau2PrMatchW0;
        Float_t mDau2pt_prW0;
        Float_t mDau2px_prW0;
        Float_t mDau2py_prW0;
        Float_t mDau2pz_prW0;

ClassDef(StFlowW0,1)
};

inline		Float_t		StFlowW0::massW0()	const	{ return mMassW0;}
inline		Float_t 	StFlowW0::rapidityW0()	const	{ return mRapidityW0;}
inline		Float_t		StFlowW0::etaW0()	const	{ return mEtaW0;}
inline		Float_t		StFlowW0::ptglobalW0()	const	{ return mPtglobalW0;}
inline		Float_t		StFlowW0::pxglobalW0()	const	{ return mPxglobalW0;}
inline		Float_t		StFlowW0::pyglobalW0()	const	{ return mPyglobalW0;}
inline		Float_t		StFlowW0::pzglobalW0()	const	{ return mPzglobalW0;}
inline		Float_t		StFlowW0::oxglobalW0()	const	{ return mOxglobalW0;}
inline		Float_t		StFlowW0::oyglobalW0()	const	{ return mOyglobalW0;}
inline		Float_t 	StFlowW0::ozglobalW0()	const	{ return mOzglobalW0;}
inline		Float_t 	StFlowW0::decaylenW0()	const	{ return mDecaylenW0;}
inline		Float_t 	StFlowW0::dcaW0()	const	{ return mDcaW0;}
inline		Float_t 	StFlowW0::dca2DW0()	const	{ return mDca2DW0;}
inline		Float_t 	StFlowW0::pathlenW0()	const	{ return mPathlenW0;}

inline		Int_t		StFlowW0::dau1idW0()	const	{ return mDau1idW0;}
inline		Float_t 	StFlowW0::dau1dcaW0()	const	{ return mDau1dcaW0;}
inline		Float_t 	StFlowW0::dau1dca2DW0()	const	{ return mDau1dca2DW0;}
inline		Int_t   	StFlowW0::dau1nhitsW0()	const	{ return mDau1nhitsW0;}
inline		Float_t 	StFlowW0::dau1dedxW0()	const	{ return mDau1dedxW0;}
inline		Float_t 	StFlowW0::dau1nSigmaW0()const	{ return mDau1nSigmaW0;}
inline		Float_t 	StFlowW0::dau1etaW0()	const	{ return mDau1etaW0;}
inline		Float_t 	StFlowW0::dau1ptW0()	const	{ return mDau1ptW0;}
inline		Float_t 	StFlowW0::dau1pxW0()	const	{ return mDau1pxW0;}
inline		Float_t		StFlowW0::dau1pyW0()	const	{ return mDau1pyW0;}
inline		Float_t 	StFlowW0::dau1pzW0()	const	{ return mDau1pzW0;}
inline		Int_t   	StFlowW0::dau1tofflagW0()const	{ return mDau1tofflagW0;}
inline		Float_t 	StFlowW0::dau1tofW0()	const	{ return mDau1tofW0;}
inline		Float_t 	StFlowW0::dau1pathlenW0()const	{ return mDau1pathlenW0;}
inline		Int_t		StFlowW0::dau1PrMatchW0()const	{ return mDau1PrMatchW0;}
inline		Float_t 	StFlowW0::dau1pt_prW0()	const	{ return mDau1pt_prW0;}
inline		Float_t 	StFlowW0::dau1px_prW0() const	{ return mDau1px_prW0;}
inline		Float_t 	StFlowW0::dau1py_prW0()	const	{ return mDau1py_prW0;}
inline		Float_t 	StFlowW0::dau1pz_prW0() const	{ return mDau1pz_prW0;}

inline		Float_t 	StFlowW0::dca1to2W0()	const	{ return mDca1to2W0;}

inline          Int_t           StFlowW0::dau2idW0()    const   { return mDau2idW0;}
inline          Float_t         StFlowW0::dau2dcaW0()   const   { return mDau2dcaW0;}
inline          Float_t         StFlowW0::dau2dca2DW0() const   { return mDau2dca2DW0;}
inline          Int_t           StFlowW0::dau2nhitsW0() const   { return mDau2nhitsW0;}
inline          Float_t         StFlowW0::dau2dedxW0()  const   { return mDau2dedxW0;}
inline          Float_t         StFlowW0::dau2nSigmaW0()const   { return mDau2nSigmaW0;}
inline          Float_t         StFlowW0::dau2etaW0()   const   { return mDau2etaW0;}
inline          Float_t         StFlowW0::dau2ptW0()    const   { return mDau2ptW0;}
inline          Float_t         StFlowW0::dau2pxW0()    const   { return mDau2pxW0;}
inline          Float_t         StFlowW0::dau2pyW0()    const   { return mDau2pyW0;}
inline          Float_t         StFlowW0::dau2pzW0()    const   { return mDau2pzW0;}
inline          Int_t           StFlowW0::dau2tofflagW0()const  { return mDau2tofflagW0;}
inline          Float_t         StFlowW0::dau2tofW0()   const   { return mDau2tofW0;}
inline          Float_t         StFlowW0::dau2pathlenW0()const  { return mDau2pathlenW0;}
inline          Int_t           StFlowW0::dau2PrMatchW0()const  { return mDau2PrMatchW0;}
inline          Float_t         StFlowW0::dau2pt_prW0() const   { return mDau2pt_prW0;}
inline          Float_t         StFlowW0::dau2px_prW0() const   { return mDau2px_prW0;}
inline          Float_t         StFlowW0::dau2py_prW0() const   { return mDau2py_prW0;}
inline          Float_t         StFlowW0::dau2pz_prW0() const   { return mDau2pz_prW0;}

#endif
