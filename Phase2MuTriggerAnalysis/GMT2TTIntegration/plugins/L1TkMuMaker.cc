/*
 * L1TkMuMaker.cc
 *
 *  Created on: Jul 31, 2014
 *      Author: Austin
 */

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"

#include "etaPhiPropagators.h"

using namespace std;
using namespace edm;
using namespace l1extra;

class L1TkMuMaker : public edm::EDProducer {
public:
  typedef math::XYZTLorentzVectorF LorentzVector;
  typedef TTTrack< Ref_PixelDigi_ >                     L1Track;
  typedef std::vector< L1Track >        L1TrackCollection;
  explicit L1TkMuMaker (const edm::ParameterSet&);
  ~L1TkMuMaker() {};

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  edm::InputTag l1tkmuInputTag;
  std::string aliasprefix_;
  unique_ptr<etaPhiProp> propTFS2;
};


L1TkMuMaker::L1TkMuMaker(const edm::ParameterSet& iConfig) {
  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix", "l1tkmus");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos)
    branchprefix.replace(branchprefix.find("_"),1,"");

  produces<std::vector<float> > (branchprefix+"muextpt" ).setBranchAlias(aliasprefix_+"_muext_pt" );
  produces<std::vector<float> > (branchprefix+"muexteta" ).setBranchAlias(aliasprefix_+"_muext_eta" );
  produces<std::vector<float> > (branchprefix+"muextphi" ).setBranchAlias(aliasprefix_+"_muext_phi" );
  produces<std::vector<int> > (branchprefix+"muextq" ).setBranchAlias(aliasprefix_+"_muext_q" );

  produces<std::vector<float> > (branchprefix+"csctfpt" ).setBranchAlias(aliasprefix_+"_csctf_pt" );
	produces<std::vector<float> > (branchprefix+"csctfeta" ).setBranchAlias(aliasprefix_+"_csctf_eta" );
	produces<std::vector<float> > (branchprefix+"csctfphi" ).setBranchAlias(aliasprefix_+"_csctf_phi" );
	produces<std::vector<int> > (branchprefix+"csctfq" ).setBranchAlias(aliasprefix_+"_csctf_q" );


  // input tags
  l1tkmuInputTag = iConfig.getParameter<edm::InputTag>("l1tkmuInputTag");
  propTFS2 = std::unique_ptr<etaPhiProp>( new etaPhiToStation2Run1TF());
}


void L1TkMuMaker::produce(edm::Event& ev, const edm::EventSetup& es){
  auto_ptr<vector<float> >	muext_pt				(new vector<float>  );
  auto_ptr<vector<float> >	muext_eta				(new vector<float>  );
  auto_ptr<vector<float> >	muext_phi				(new vector<float>  );
  auto_ptr<vector<int> >		muext_q					(new vector<int>		);

  auto_ptr<vector<float> >	csctf_pt				(new vector<float>  );
	auto_ptr<vector<float> >	csctf_eta				(new vector<float>  );
	auto_ptr<vector<float> >	csctf_phi				(new vector<float>  );
	auto_ptr<vector<int> >		csctf_q					(new vector<int>		);

  Handle<L1TkMuonParticleCollection> ttH;
  ev.getByLabel(l1tkmuInputTag, ttH);
  const L1TkMuonParticleCollection& tts(*ttH.product());

  for (const L1TkMuonParticle & l1tkmu : tts){
  	const L1MuonParticleExtendedRef & l1mu_ref = l1tkmu.getMuExtendedRef();
  	const L1MuRegionalCand & cscCand = l1mu_ref.get()->cscCand();

  	bool is_good = true;

  	if(l1mu_ref.get()->eta() >= 1.1) {
  		is_good = false;
  	}

  	if(is_good) {
			muext_pt->push_back(l1mu_ref.get()->pt());
			muext_eta->push_back(l1mu_ref.get()->eta());
			muext_phi->push_back(l1mu_ref.get()->phi());
			muext_q->push_back(l1mu_ref.get()->charge()>0? 1: -1);
  	}
  	else {
  		muext_pt->push_back(-10);
			muext_eta->push_back(-10);
			muext_phi->push_back(-10);
			muext_q->push_back(-10);
  	}

    if(cscCand.chargeValid() && is_good) {
    	csctf_pt->push_back(cscCand.ptValue());
    	csctf_eta->push_back(cscCand.etaValue());
    	csctf_phi->push_back(cscCand.phiValue());
    	csctf_q->push_back(cscCand.chargeValue());
    }
    else {
    	csctf_pt->push_back(-10);
			csctf_eta->push_back(-10);
			csctf_phi->push_back(-10);
			csctf_q->push_back(-10);
    }
  }

  ev.put(muext_pt, aliasprefix_+"muext_pt");
	ev.put(muext_eta, aliasprefix_+"muext_eta");
	ev.put(muext_phi, aliasprefix_+"muext_phi");
	ev.put(muext_q , aliasprefix_+"muext_q");

	ev.put(csctf_pt, aliasprefix_+"csctf_pt");
	ev.put(csctf_eta, aliasprefix_+"csctf_eta");
	ev.put(csctf_phi, aliasprefix_+"csctf_phi");
	ev.put(csctf_q , aliasprefix_+"csctf_q");

}

DEFINE_FWK_MODULE(L1TkMuMaker);

