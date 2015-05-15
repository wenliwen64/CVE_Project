#ifndef StFlowStrangenessEvent_hh
#define StFlowStrangenessEvent_hh

//#include <iostream.h>
#include "TObject.h"
#include "StThreeVectorF.hh"

#include "StFlowV0.h"
#include "StFlowW0.h"

class zdcTriggerDetector;
class TClonesArray;

class StFlowStrangenessEvent : public TObject {
public:
	StFlowStrangenessEvent();
        virtual ~StFlowStrangenessEvent();
	

        TClonesArray* v0(){return fV0s;}
        TClonesArray* w0(){return fW0s;}

	Int_t 	eventId()			const;
	Int_t 	runId()				const;
	Int_t	numberOfV0s()			const;
        Int_t   numberOfW0s()                   const;
	UShort_t	refMult()		const;
        UShort_t        refMultTOF()               const;
        Float_t primaryVertexX()        const;
        Float_t primaryVertexY()        const;
	Float_t	primaryVertexZ()	const;
        Float_t   vpdEWdiff()     const;
        Double_t bbcEast() const;
        Double_t bbcWest() const;
        Double_t bbcCoin() const;
	Double_t zdcCoin() const;	
	Int_t trgId()	const;
	Float_t magn()	const;


	void	Clear		(Option_t* option ="");
	void	AddV0	(StFlowV0* v0);
        void    AddW0   (StFlowW0* w0);
	void	SetEventId	(const Int_t eventidS)		{ mEventId	= eventidS; }
	void	SetRunId(const Int_t runidS)			{ mRunId	= runidS; }
        void    SetTrigger	(Int_t trigger)			{ mTrigger	= trigger; }
	void	SetPs_1		(Float_t Ps)			{ mPreScale_1	= Ps;}
        void    SetPs_2         (Float_t Ps)                    { mPreScale_2   = Ps;}
	void	SetBz		(Float_t B)			{ mBz		= B;}
        void    SetPrimaryVertexX  (Float_t x)                  { mPrimaryVertexX = x; }
        void    SetPrimaryVertexY  (Float_t y)                  { mPrimaryVertexY = y; }
	void	SetPrimaryVertexZ  (Float_t z)		        { mPrimaryVertexZ = z; }
        void    SetVPDEWdiff(Float_t d) {mVPDEWdiff = d;}
	void	SetRefMult	(const UShort_t refMultS)	{ mRefMult	= refMultS; }
        void    SetRefMultTOF      (const UShort_t refMultTOFS)       { mRefMultTOF      = refMultTOFS; }
        void    SetBBCE (Double_t e)    {mBBCEast=e;}
        void    SetBBCW (Double_t w)    {mBBCWest=w;}
        void    SetBBCX (Double_t x)    {mBBCCoin=x;}
        void    SetZDCX (Double_t zx)    {mZDCCoin=zx;}
	void	SetTrgId (const Int_t triggerIdS)		{ mTriggerId	= triggerIdS;}
	void	SetMagn	(Float_t magnfieldS)			{ mMagn		= magnfieldS;}
	void	SetNumV0 (Int_t	NumV0)	{mNoTracks = NumV0;}
        void    SetNumW0 (Int_t NumW0)  {mNoLambdas = NumW0;}

//        static Int_t        numCtorCalled;
//        static Int_t        numDtorCalled;



private:
        Int_t   mRunId;
	Int_t	mEventId;
	Int_t	mTrigger;
        Float_t mPreScale_1;
	Float_t mPreScale_2;
	Float_t mBz;
        Float_t mPrimaryVertexX;
        Float_t mPrimaryVertexY;
        Float_t mPrimaryVertexZ;
        Float_t mVPDEWdiff;
        UShort_t mRefMult;
        UShort_t mRefMultTOF;
        Double_t mBBCEast;
        Double_t mBBCWest;
        Double_t mBBCCoin;
	Double_t mZDCCoin;
	Int_t	mTriggerId;
	Float_t	mMagn;

	Int_t	mNoTracks;
	TClonesArray*		fV0s;
	static TClonesArray*	fgV0s;
        Int_t   mNoLambdas;
        TClonesArray*           fW0s;
        static TClonesArray*    fgW0s;

ClassDef(StFlowStrangenessEvent,1)
};
inline	Int_t 	StFlowStrangenessEvent::eventId()				const	{ return mEventId; }
inline	Int_t 	StFlowStrangenessEvent::runId()					const	{ return mRunId; }
inline	Int_t	StFlowStrangenessEvent::numberOfV0s()			const	{ return mNoTracks; }
inline  Int_t   StFlowStrangenessEvent::numberOfW0s()                   const   { return mNoLambdas; }
inline  Float_t StFlowStrangenessEvent::primaryVertexX()                const   { return mPrimaryVertexX;}
inline  Float_t StFlowStrangenessEvent::primaryVertexY()                const   { return mPrimaryVertexY;}
inline	Float_t	StFlowStrangenessEvent::primaryVertexZ()		const	{ return mPrimaryVertexZ;}
inline  Float_t   StFlowStrangenessEvent::vpdEWdiff()                   const   { return mVPDEWdiff;}
inline	UShort_t		StFlowStrangenessEvent::refMult()				const	{ return mRefMult; }
inline  UShort_t                StFlowStrangenessEvent::refMultTOF()                               const   { return mRefMultTOF; }
inline Double_t StFlowStrangenessEvent::bbcEast()                       const   {return mBBCEast;}
inline Double_t StFlowStrangenessEvent::bbcWest()                       const   {return mBBCWest;}
inline Double_t StFlowStrangenessEvent::bbcCoin()                       const   {return mBBCCoin;}
inline Double_t StFlowStrangenessEvent::zdcCoin()                       const   {return mZDCCoin;}
inline	Int_t	StFlowStrangenessEvent::trgId()				const	{ return mTriggerId;}
inline	Float_t	StFlowStrangenessEvent::magn()				const	{ return mMagn;}
#endif
