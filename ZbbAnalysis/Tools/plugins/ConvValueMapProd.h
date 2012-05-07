#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


class ConvValueMapProd : public edm::EDProducer {
    public:
        explicit ConvValueMapProd(const edm::ParameterSet&);
        ~ConvValueMapProd();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag gsfLabel_;
        edm::InputTag tkLabel_;

};
