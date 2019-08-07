#include "RecoTauTag/RecoTau/interface/MVADMDiscriminationProducerBase.h"

#include <string>

using namespace reco;

// default constructor; must not be called
template<class TauType, class TauDiscriminator>
MVADMDiscriminationProducerBase<TauType, TauDiscriminator>::MVADMDiscriminationProducerBase()
{
   throw cms::Exception("MVADMDiscriminationProducerBase") << " -- default ctor called; derived classes must call " <<
      "MVADMDiscriminationProducerBase(const ParameterSet&)";
}

//--- standard constructor from PSet
template<class TauType, class TauDiscriminator>
MVADMDiscriminationProducerBase<TauType, TauDiscriminator>::MVADMDiscriminationProducerBase(const edm::ParameterSet& iConfig)
  : moduleLabel_(iConfig.getParameter<std::string>("@module_label"))
{
   // tau collection to discriminate
   TauProducer_        = iConfig.getParameter<edm::InputTag>("PATTauProducer");
   Tau_token= consumes<TauCollection>(TauProducer_);

   // register product
   produces<TauDiscriminator>("DMother");
   produces<TauDiscriminator>("DM0");
   produces<TauDiscriminator>("DM1");
   produces<TauDiscriminator>("DM2");
   produces<TauDiscriminator>("DM10");
   produces<TauDiscriminator>("DM11");
   produces<TauDiscriminator>();
   
}

template<class TauType, class TauDiscriminator>
void MVADMDiscriminationProducerBase<TauType, TauDiscriminator>::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
   tauIndex_=0;
   // setup function - does nothing in base, but can be overridden to retrieve PV or other stuff
   beginEvent(event, eventSetup);

   // retrieve the tau collection to discriminate
   edm::Handle<TauCollection> taus;
   event.getByToken(Tau_token, taus);

   // output product
   auto output1 = std::make_unique<TauDiscriminator>(TauRefProd(taus)); //(MVA) dm="other"
   auto output2 = std::make_unique<TauDiscriminator>(TauRefProd(taus)); //(MVA) dm=0
   auto output3 = std::make_unique<TauDiscriminator>(TauRefProd(taus)); //(MVA) dm=1
   auto output4 = std::make_unique<TauDiscriminator>(TauRefProd(taus)); //(MVA) dm=2
   auto output5 = std::make_unique<TauDiscriminator>(TauRefProd(taus)); //(MVA) dm=10
   auto output6 = std::make_unique<TauDiscriminator>(TauRefProd(taus)); //(MVA) dm=11
   // output7 defines MVA dm as integer value based on size of class scores
   auto output7 = std::make_unique<TauDiscriminator>(TauRefProd(taus));

   size_t nTaus = taus->size();

   // loop over taus
   for( size_t iTau = 0; iTau < nTaus; ++iTau )
   {
      // get reference to tau
      TauRef tauRef(taus, iTau);

      // this tau passes the prereqs, call our implemented discrimination function
      std::vector<float> result = discriminate(tauRef); ++tauIndex_;

      // store the result of this tau into our new discriminator
      if(result.size()<3) {
        output1->setValue(iTau, 0.);
        output2->setValue(iTau, 0.);
        output3->setValue(iTau, 0.);
        output4->setValue(iTau, 0.);
        output5->setValue(iTau, 0.);
        output6->setValue(iTau, 0.);
        output7->setValue(iTau, -1);
      }
      if(result.size()==4) {
        output1->setValue(iTau, result[0]);
        output2->setValue(iTau, result[2]);
        output3->setValue(iTau, result[1]);
        output4->setValue(iTau, result[3]); 
        output5->setValue(iTau, 0.0);
        output6->setValue(iTau, 0.0);

        if     (result[2]>result[0]&&result[2]>result[1]&&result[2]>result[3]) output7->setValue(iTau, 0);
        else if(result[1]>result[0]&&result[1]>result[2]&&result[1]>result[3]) output7->setValue(iTau, 1);
        else if(result[3]>result[0]&&result[3]>result[1]&&result[3]>result[2]) output7->setValue(iTau, 2);
        else output7->setValue(iTau, -1);
      }
      if(result.size()==3){
        output1->setValue(iTau, result[0]);
        output2->setValue(iTau, 0.0);
        output3->setValue(iTau, 0.0);
        output4->setValue(iTau, 0.0);
        output5->setValue(iTau, result[1]);
        output6->setValue(iTau, result[2]);
        if     (result[1]>result[0] && result[1]>result[2]) output7->setValue(iTau, 10);
        else if(result[2]>result[0] && result[2]>result[1]) output7->setValue(iTau, 11); 
        else output7->setValue(iTau, -1);
      }
   }
   event.put(std::move(output1), (std::string)"DMother");
   event.put(std::move(output2), (std::string)"DM0");
   event.put(std::move(output3), (std::string)"DM1");
   event.put(std::move(output4), (std::string)"DM2"); 
   event.put(std::move(output5), (std::string)"DM10");
   event.put(std::move(output6), (std::string)"DM11");
   event.put(std::move(output7));
   // function to put additional information into the event - does nothing in base, but can be overridden in derived classes
   endEvent(event);
}


// compile our desired types and make available to linker
template class MVADMDiscriminationProducerBase<pat::Tau, pat::PATTauDiscriminator>;
