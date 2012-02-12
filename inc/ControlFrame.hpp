#ifndef CONTROL_FRAME_HPP
#define CONTROL_FRAME_HPP

#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>

class ControlFrame: public TGMainFrame {

private:
  TGCompositeFrame *_cFrame;
  TGTextButton *_start, *_pause, *fExit;
  Bool_t _isRunning, _isPaused;

public:
  ControlFrame( const TGWindow *p, UInt_t w, UInt_t h );
  virtual ~ControlFrame();
  // slots
  void ChangeStartLabel();
  void ChangePauseLabel();

ClassDef(ControlFrame, 0)
};

#endif // CONTROL_FRAME_HPP
