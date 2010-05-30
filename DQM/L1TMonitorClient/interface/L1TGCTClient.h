#ifndef DQM_L1TMONITORCLIENT_L1TGCTCLIENT_H
#define DQM_L1TMONITORCLIENT_L1TGCTCLIENT_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

class L1TGCTClient: public edm::EDAnalyzer {

 public:

  /// Constructor
  L1TGCTClient(const edm::ParameterSet& ps);
  
  /// Destructor
  virtual ~L1TGCTClient();
 
 protected:

  /// BeginJob
  void beginJob(void);

  /// BeginRun
  void beginRun(const edm::Run& r, const edm::EventSetup& c);

  /// Analyze
  void analyze(const edm::Event& e, const edm::EventSetup& c) ;

  void beginLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
                            const edm::EventSetup& context) ;

  /// DQM Client Diagnostic
  void endLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
                          const edm::EventSetup& c);

  /// EndRun
  void endRun(const edm::Run& r, const edm::EventSetup& c);

  /// Endjob
  void endJob();

 private:

  void makeXProjection(TH2D* input, MonitorElement* output);
  void makeYProjection(TH2D* input, MonitorElement* output);

  DQMStore* dbe_;   
  std::string monitorDir_; 
  int counterLS_;      ///counter 
  int counterEvt_;     ///counter 
  int prescaleLS_;     ///units of lumi sections 
  int prescaleEvt_;    ///prescale on number of events

  MonitorElement *l1GctIsoEmOccEta_;
  MonitorElement *l1GctIsoEmOccPhi_;
  MonitorElement *l1GctNonIsoEmOccEta_;
  MonitorElement *l1GctNonIsoEmOccPhi_;
  MonitorElement *l1GctAllJetsOccEta_;
  MonitorElement *l1GctAllJetsOccPhi_;
  MonitorElement *l1GctCenJetsOccEta_;
  MonitorElement *l1GctCenJetsOccPhi_;
  MonitorElement *l1GctForJetsOccEta_;
  MonitorElement *l1GctForJetsOccPhi_;
  MonitorElement *l1GctTauJetsOccEta_;
  MonitorElement *l1GctTauJetsOccPhi_;

};

#endif
