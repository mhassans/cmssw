#ifndef RecoTauTag_RecoTau_MVADMDiscriminationProducerBase_H_
#define RecoTauTag_RecoTau_MVADMDiscriminationProducerBase_H_

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"


#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PATTauDiscriminator.h"


template<class TauType, class TauDiscriminator>
class MVADMDiscriminationProducerBase : public edm::stream::EDProducer<> {
  public:
    // setup framework types for this tautype
    typedef std::vector<TauType>        TauCollection;
    typedef edm::Ref<TauCollection>     TauRef;
    typedef edm::RefProd<TauCollection> TauRefProd;

    // standard constructor from PSet
    explicit MVADMDiscriminationProducerBase(const edm::ParameterSet& iConfig);

    // default constructor must not be called - it will throw an exception
    // derived!  classes must call the parameterset constructor.
    MVADMDiscriminationProducerBase();

    ~MVADMDiscriminationProducerBase() override {}

    void produce(edm::Event&, const edm::EventSetup&) override;

    // called at the beginning of every event - override if necessary.
    virtual void beginEvent(const edm::Event&, const edm::EventSetup&) {}

    // abstract functions implemented in derived classes.
    virtual std::vector<float> discriminate(const TauRef& tau) const = 0;

    // called at the end of event processing - override if necessary.
    virtual void endEvent(edm::Event&) {}

    struct TauDiscInfo {
      edm::InputTag label;
      edm::Handle<TauDiscriminator> handle;
      edm::EDGetTokenT<TauDiscriminator> disc_token;
      double cut;
      void fill(const edm::Event& evt) { 
	evt.getByToken(disc_token, handle); 
      };
    };

    //std::unique_ptr<pat::PATTauDiscriminator> output1;
    //std::unique_ptr<pat::PATTauDiscriminator> output2;
    //std::unique_ptr<pat::PATTauDiscriminator> output3;
    //std::unique_ptr<pat::PATTauDiscriminator> output4;
    //std::unique_ptr<pat::PATTauDiscriminator> output5;
    //std::unique_ptr<pat::PATTauDiscriminator> output6;
    //std::unique_ptr<pat::PATTauDiscriminator> output7;
    

  protected:

    edm::InputTag TauProducer_;

    std::string moduleLabel_;
    edm::EDGetTokenT<TauCollection> Tau_token;

    // current tau
    size_t tauIndex_;

  private:

};

// define our implementations
typedef MVADMDiscriminationProducerBase<pat::Tau, pat::PATTauDiscriminator>
  PATTauMVADMDiscriminationProducerBase;


#endif
