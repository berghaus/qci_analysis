#ifndef CONTROL_FRAME_HPP
#define CONTROL_FRAME_HPP

#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>

class ControlFrame : public TGMainFrame {

private:
   TGCompositeFrame *fCframe;
   TGTextButton     *fStart, *fPause, *fExit;
   Bool_t            start, pause;

public:
   ControlFrame(const TGWindow *p, UInt_t w, UInt_t h);
   virtual ~ControlFrame();
   // slots
   void ChangeStartLabel();
   void ChangePauseLabel();

   ClassDef(ControlFrame, 0)
};


#endif // CONTROL_FRAME_HPP
