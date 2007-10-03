#ifndef RecoBTag_JetTagComputer_JetTagProducer_h
#define RecoBTag_JetTagComputer_JetTagProducer_h

// system include files
#include <string>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"

class JetTagProducer : public edm::EDProducer {
    public:
	explicit JetTagProducer(const edm::ParameterSet&);
	~JetTagProducer();

	virtual void produce(edm::Event&, const edm::EventSetup&);

    private:
	void setup(edm::Event&);

	const JetTagComputer			*m_computer;
	std::string				m_jetTagComputer;
	std::map<std::string, edm::InputTag>	m_tagInfoLabels;
	std::vector<edm::InputTag>		m_tagInfos;
};

#endif // RecoBTag_JetTagComputer_JetTagProducer_h
