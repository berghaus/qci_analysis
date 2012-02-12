#include "ControlFrame.hpp"

using namespace std;

ClassImp( ControlFrame );

void ControlFrame::ChangeStartLabel() {
  // Slot connected to the Clicked() signal. 
  // It will toggle labels "Start" and "Stop".

  if ( !_isRunning ) {
    _start->SetText( "&Stop" );
    _start->SetState( kButtonDown );
    _isRunning = kTRUE;
  } else {
    _start->SetText( "&Start" );
    _start->SetState( kButtonUp );
    _isRunning = kFALSE;
  }

}

void ControlFrame::ChangePauseLabel() {
  // Slot connected to the Clicked() signal. 
  // It will toggle labels "Resume" and "Pause".

  if ( !_isPaused ) {
    _pause->SetText( "&Resume" );
    _pause->SetState( kButtonDown );
    _isPaused = kTRUE;
  } else {
    _pause->SetText( "&Pause" );
    _pause->SetState( kButtonUp );
    _isPaused = kFALSE;
  }

}

ControlFrame::ControlFrame( const TGWindow *p, UInt_t w, UInt_t h ) :
    TGMainFrame( p, w, h ) {

  // Create a horizontal frame containing buttons
  _cFrame = new TGCompositeFrame( this, 170, 20, kHorizontalFrame | kFixedWidth );

  _start = new TGTextButton( _cFrame, "&Start" );
  _start->Connect( "Clicked()", "ControlFrame", this, "ChangeStartLabel()" );
  _cFrame->AddFrame( _start, new TGLayoutHints( kLHintsTop | kLHintsExpandX, 3, 2, 2, 2 ) );
  _start->SetToolTipText( "Click to toggle the button label (Start/Stop)" );
  _isRunning = kFALSE;

  _pause = new TGTextButton( _cFrame, "&Pause" );
  _pause->Connect( "Clicked()", "ControlFrame", this, "ChangePauseLabel()" );
  _pause->SetToolTipText( "Click to toggle the button label (Pause/Resume)" );
  _cFrame->AddFrame( _pause, new TGLayoutHints( kLHintsTop | kLHintsExpandX, 3, 2, 2, 2 ) );
  _isPaused = kFALSE;

  AddFrame( _cFrame, new TGLayoutHints( kLHintsCenterX, 2, 2, 5, 1 ) );

  fExit = new TGTextButton( this, "&Exit ", "gApplication->Terminate(0)" );
  AddFrame( fExit, new TGLayoutHints( kLHintsTop | kLHintsExpandX, 5, 5, 2, 2 ) );

  SetWindowName( "Change Labels" );

  MapSubwindows();
  Resize( GetDefaultSize() );
  MapWindow();
}

ControlFrame::~ControlFrame() {
  // Clean up all widgets, frames and layouthints that were used
  _cFrame->Cleanup();
  Cleanup();
}

