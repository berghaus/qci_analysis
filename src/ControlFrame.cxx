#include "ControlFrame.hpp"

using namespace std;

ClassImp(ControlFrame);

void ControlFrame::ChangeStartLabel()
{
  // Slot connected to the Clicked() signal. 
  // It will toggle labels "Start" and "Stop".
  
  fStart->SetState(kButtonDown);
  if (!start) {
     fStart->SetText("&Stop");
     start = kTRUE;
  } else {
     fStart->SetText("&Start");
     start = kFALSE;
  }
  fStart->SetState(kButtonUp);
}

void ControlFrame::ChangePauseLabel()
{
  // Slot connected to the Clicked() signal. 
  // It will toggle labels "Resume" and "Pause".
  
  fPause->SetState(kButtonDown);
  if (!pause) {
     fPause->SetText("&Resume");
     pause = kTRUE;
  } else {
     fPause->SetText("&Pause");
     pause = kFALSE;
  }
  fPause->SetState(kButtonUp);
}

ControlFrame::ControlFrame(const TGWindow *p, UInt_t w, UInt_t h) :
  TGMainFrame(p, w, h)
{
   // Create a horizontal frame containing buttons
   fCframe = new TGCompositeFrame(this, 170, 20, kHorizontalFrame|kFixedWidth);
   
   fStart = new TGTextButton(fCframe, "&Start");
   fStart->Connect("Clicked()", "ControlFrame", this, "ChangeStartLabel()");
   fCframe->AddFrame(fStart, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 
                                               3, 2, 2, 2));
   fStart->SetToolTipText("Click to toggle the button label (Start/Stop)");
   start = kFALSE;
   
   fPause = new TGTextButton(fCframe, "&Pause");
   fPause->Connect("Clicked()", "ControlFrame", this, "ChangePauseLabel()");
   fPause->SetToolTipText("Click to toggle the button label (Pause/Resume)");
   fCframe->AddFrame(fPause, new TGLayoutHints(kLHintsTop | kLHintsExpandX,
                                               3, 2, 2, 2));
   pause = kFALSE;
   
   AddFrame(fCframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 5, 1));

   fExit = new TGTextButton(this, "&Exit ","gApplication->Terminate(0)");
   AddFrame(fExit, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,2,2));
   
   SetWindowName("Change Labels");
   
   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();
}


ControlFrame::~ControlFrame()
{
   // Clean up all widgets, frames and layouthints that were used
   fCframe->Cleanup();
   Cleanup();
}


