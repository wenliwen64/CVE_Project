#include "TObject.h"
#include "TClonesArray.h"

#include "StFlowV0.h"
#include "StFlowW0.h"
#include "StFlowStrangenessEvent.h"


ClassImp(StFlowStrangenessEvent)

TClonesArray *StFlowStrangenessEvent::fgV0s = 0;
TClonesArray *StFlowStrangenessEvent::fgW0s = 0;
//Int_t StFlowStrangenessEvent::numCtorCalled=0;
//Int_t StFlowStrangenessEvent::numDtorCalled=0;


StFlowStrangenessEvent::StFlowStrangenessEvent() : TObject()
{
       if (!fgV0s)	fgV0s = new TClonesArray("StFlowV0", 3000);
       fV0s  = fgV0s;
       mNoTracks  = 0;

       if (!fgW0s)      fgW0s = new TClonesArray("StFlowW0", 2000);
       fW0s  = fgW0s;
       mNoLambdas  = 0;
}


StFlowStrangenessEvent::~StFlowStrangenessEvent()
{
  Clear();
}


void	StFlowStrangenessEvent::AddV0( StFlowV0* v0 ) 
{
  	TClonesArray &v0s = *fV0s;
  	new(v0s[mNoTracks++]) StFlowV0(v0);
}


void    StFlowStrangenessEvent::AddW0( StFlowW0* w0 )
{
        TClonesArray &w0s = *fW0s;
        new(w0s[mNoLambdas++]) StFlowW0(w0);
}

void	StFlowStrangenessEvent::Clear(Option_t *option)
{
  if(fV0s)   fV0s->Clear(option);     mNoTracks = 0;
  if(fW0s)   fW0s->Clear(option);     mNoLambdas = 0;
}
