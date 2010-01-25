// -*- C++ -*-
//
// Package:     Core
// Class  :     FWGUISubviewArea
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Fri Feb 15 14:13:33 EST 2008
// $Id: FWGUISubviewArea.cc,v 1.33 2009/10/07 19:02:31 amraktad Exp $
//

// system include files
#include <assert.h>
#include <stdexcept>
#include <iostream>
#include <boost/bind.hpp>

#include "TSystem.h"
#include "TGButton.h"
#include "TGSplitFrame.h"
#include "TGFont.h"
#include "TGLabel.h"
#include "TEveWindow.h"

#include "Fireworks/Core/interface/FWGUISubviewArea.h"
#include "Fireworks/Core/interface/FWGUIManager.h"
#include "Fireworks/Core/src/FWCheckBoxIcon.h"

//
// constructors and destructor
//
FWGUISubviewArea::FWGUISubviewArea(TEveCompositeFrame* ef, TGCompositeFrame* parent, Int_t height) :
   TGHorizontalFrame(parent, 20, height),
   m_frame (ef),
   m_swapButton(0),
   m_undockButton(0), m_dockButton(0),
   m_closeButton(0),
   m_infoButton(0)
{
   UInt_t lh = kLHintsNormal | kLHintsExpandX | kLHintsExpandY;

   // info
   m_infoButton = new TGPictureButton(this,infoIcon());
   m_infoButton->ChangeOptions(kRaisedFrame);
   m_infoButton->SetDisabledPicture(infoDisabledIcon());
   AddFrame(m_infoButton, new TGLayoutHints(lh));
   m_infoButton->AllowStayDown(kTRUE);
   m_infoButton->Connect("Pressed()","FWGUISubviewArea",this,"selectButtonDown()");
   m_infoButton->Connect("Released()","FWGUISubviewArea",this,"selectButtonUp()");
   m_infoButton->SetToolTipText("Edit View");

   //swap
   m_swapButton = new TGPictureButton(this, swapIcon());
   m_swapButton->SetDisabledPicture(swapDisabledIcon());
   m_swapButton->SetToolTipText("Swap with the first view or select window to swap with click on toolbar.");
   m_swapButton->ChangeOptions(kRaisedFrame);
   m_swapButton->SetHeight(height);
   AddFrame(m_swapButton, new TGLayoutHints(lh));
   m_swapButton->Connect("Clicked()","FWGUISubviewArea",this,"swap()");

   // dock 
   if (dynamic_cast<const TGMainFrame*>(ef->GetParent()))
   {  
      m_dockButton = new TGPictureButton(this, dockIcon());
      m_dockButton->ChangeOptions(kRaisedFrame);
      m_dockButton->ChangeOptions(kRaisedFrame);
      m_dockButton->SetToolTipText("Dock view");
      m_dockButton->SetHeight(height);
      AddFrame(m_dockButton, new TGLayoutHints(lh));
      m_dockButton->Connect("Clicked()", "FWGUISubviewArea",this,"dock()");
   }
   else
   {
      // undock  
      m_undockButton = new TGPictureButton(this, undockIcon());
      m_undockButton->ChangeOptions(kRaisedFrame);
      m_undockButton->SetDisabledPicture(undockDisabledIcon());
      m_undockButton->SetToolTipText("Undock view to own window");
      m_undockButton->SetHeight(height);
      AddFrame(m_undockButton, new TGLayoutHints(lh));
      m_undockButton->Connect("Clicked()", "FWGUISubviewArea",this,"undock()");
   }
   // destroy
   m_closeButton = new TGPictureButton(this,closeIcon());
   m_closeButton->ChangeOptions(kRaisedFrame);
   m_closeButton->SetToolTipText("Close view");
   m_closeButton->SetHeight(height);
   AddFrame(m_closeButton, new TGLayoutHints(lh));
   m_closeButton->Connect("Clicked()", "FWGUISubviewArea",this,"destroy()");

   FWGUIManager::getGUIManager()->connectSubviewAreaSignals(this);
}

FWGUISubviewArea::~FWGUISubviewArea()
{
   //std::cout <<"IN dstr FWGUISubviewArea"<<std::endl;
   m_closeButton->Disconnect("Clicked()", this,"destroy()");
   m_infoButton->Disconnect("Pressed()",this,"selectButtonDown()");
   m_infoButton->Disconnect("Released()",this,"selectButtonUp()");
}

//______________________________________________________________________________
//
// actions
//

void
FWGUISubviewArea::selectButtonDown()
{
   selected_(this);
}

void
FWGUISubviewArea::selectButtonUp()
{
   unselected_(this);
}

void
FWGUISubviewArea::unselect()
{
   m_infoButton->SetState(kButtonUp);
}

void
FWGUISubviewArea::setSwapIcon(bool isOn)
{
  m_swapButton->SetEnabled(isOn);
}

void
FWGUISubviewArea::swap()
{
   swap_(this);
}

void
FWGUISubviewArea::destroy()
{
   goingToBeDestroyed_(this);
}

void
FWGUISubviewArea::undock()
{
    TTimer::SingleShot(50, m_frame->GetEveWindow()->ClassName(), m_frame->GetEveWindow(), "UndockWindow()");
}

void
FWGUISubviewArea::dock()
{
   TGWindow* w = (TGWindow*)(m_frame->GetParent());
   w->UnmapWindow();
   TTimer::SingleShot(0, m_frame->ClassName(), m_frame, "MainFrameClosed()");
}

//
// const member functions
//
bool
FWGUISubviewArea::isSelected() const
{
   return m_infoButton->IsDown();
}


//______________________________________________________________________________
TEveWindow*
FWGUISubviewArea::getEveWindow()
{
   return m_frame->GetEveWindow();
}

//______________________________________________________________________________
// static member functions
//
const TGPicture *
FWGUISubviewArea::swapIcon()
{
   static const TGPicture* s_icon = 0;
   if(0== s_icon) {
      const char* cmspath = gSystem->Getenv("CMSSW_BASE");
      if(0 == cmspath) {
         throw std::runtime_error("CMSSW_BASE environment variable not set");
      }
      s_icon = gClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"moveup.png");
   }
   return s_icon;
}

const TGPicture *
FWGUISubviewArea::swapDisabledIcon()
{
   static const TGPicture* s_icon = 0;
   if(0== s_icon) {
      const char* cmspath = gSystem->Getenv("CMSSW_BASE");
      if(0 == cmspath) {
         throw std::runtime_error("CMSSW_BASE environment variable not set");
      }
      s_icon = gClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"moveup-disabled.png");
   }
   return s_icon;
}

const TGPicture *
FWGUISubviewArea::closeIcon()
{
   static const TGPicture* s_icon = 0;
   if(0== s_icon) {
      const char* cmspath = gSystem->Getenv("CMSSW_BASE");
      if(0 == cmspath) {
         throw std::runtime_error("CMSSW_BASE environment variable not set");
      }
      s_icon = gClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"delete.png");
   }
   return s_icon;
}

const TGPicture *
FWGUISubviewArea::closeDisabledIcon()
{
   static const TGPicture* s_icon = 0;
   if(0== s_icon) {
      const char* cmspath = gSystem->Getenv("CMSSW_BASE");
      if(0 == cmspath) {
         throw std::runtime_error("CMSSW_BASE environment variable not set");
      }
      s_icon = gClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"delete-disabled.png");
   }
   return s_icon;
}


const TGPicture *
FWGUISubviewArea::undockIcon()
{
   static const TGPicture* s_icon = 0;
   if(0== s_icon) {
      const char* cmspath = gSystem->Getenv("CMSSW_BASE");
      if(0 == cmspath) {
         throw std::runtime_error("CMSSW_BASE environment variable not set");
      }
      s_icon = gClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"expand.png");
   }
   return s_icon;
}

const TGPicture *
FWGUISubviewArea::dockIcon()
{
   static const TGPicture* s_icon = 0;
   if(0== s_icon) {
      const char* cmspath = gSystem->Getenv("CMSSW_BASE");
      if(0 == cmspath) {
         throw std::runtime_error("CMSSW_BASE environment variable not set");
      }
      s_icon = gClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"dock.png");
   }
   return s_icon;
}

const TGPicture *
FWGUISubviewArea::undockDisabledIcon()
{
   static const TGPicture* s_icon = 0;
   if(0== s_icon) {
      const char* cmspath = gSystem->Getenv("CMSSW_BASE");
      if(0 == cmspath) {
         throw std::runtime_error("CMSSW_BASE environment variable not set");
      }
      s_icon = gClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"expand-disabled.png");
   }
   return s_icon;
}

const TGPicture *
FWGUISubviewArea::infoIcon()
{
   static const TGPicture* s_icon = 0;
   if(0== s_icon) {
      const char* cmspath = gSystem->Getenv("CMSSW_BASE");
      if(0 == cmspath) {
         throw std::runtime_error("CMSSW_BASE environment variable not set");
      }
      s_icon = gClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"info.png");
   }
   return s_icon;
}

const TGPicture *
FWGUISubviewArea::infoDisabledIcon()
{
   static const TGPicture* s_icon = 0;
   if(0== s_icon) {
      const char* cmspath = gSystem->Getenv("CMSSW_BASE");
      if(0 == cmspath) {
         throw std::runtime_error("CMSSW_BASE environment variable not set");
      }
      s_icon = gClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"info-disabled.png");
   }
   return s_icon;
}

void
FWGUISubviewArea::setInfoButton(bool downp)
{
   m_infoButton->SetState(downp ? kButtonEngaged : kButtonUp, false);
}

FWGUISubviewArea*
FWGUISubviewArea::getToolBarFromWindow(TEveWindow* w)
{
   // horizontal frame
   TGFrameElement *el = (TGFrameElement*) w->GetEveFrame()->GetList()->First();
   TGCompositeFrame* hf = (TGCompositeFrame*)el->fFrame;
   // subview last in the horizontal frame
   el = (TGFrameElement*)hf->GetList()->Last();
   FWGUISubviewArea* ar = dynamic_cast<FWGUISubviewArea*>(el->fFrame);
   return ar;
}
