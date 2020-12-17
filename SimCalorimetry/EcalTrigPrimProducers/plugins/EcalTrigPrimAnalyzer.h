// -*- C++ -*-
//
// Class:      EcalTrigPrimAnalyzer
//
/**\class EcalTrigPrimAnalyzer

 Description: rereads the result of the EcalTrigPrimProducer

*/
//
// Original Author:  Ursula Berthon
//         Created:  Thu Jul 4 11:38:38 CEST 2005
//
//

// system include files
//#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "CalibCalorimetry/EcalTPGTools/interface/EcalTPGScale.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/EcalBarrelGeometryRecord.h"
#include "Geometry/Records/interface/EcalEndcapGeometryRecord.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TTree.h>
#include <string>
#include <vector>

//
// class declaration
//

class EcalTrigPrimAnalyzer : public edm::one::EDAnalyzer<> {
public:
  explicit EcalTrigPrimAnalyzer(const edm::ParameterSet &);
  ~EcalTrigPrimAnalyzer() override;

  void analyze(const edm::Event &, const edm::EventSetup &) override;
  void endJob() override;

private:
  // for histos of nr of hits
  std::vector<std::string> ecal_parts_;
  TH1I *ecal_et_[2];
  TH1I *ecal_tt_[2];
  TH1I *ecal_fgvb_[2];
  TH1I *histEndc, *histBar;
  TFile *histfile_;
  TH2F *hTPvsRechit_;
  TH1F *hTPoverRechit_;
  TTree *tree_;

  int iphi_, ieta_, tpgADC_, ttf_, fg_;
  float eRec_, tpgGeV_;

  edm::InputTag label_;

  edm::InputTag rechits_labelEB_;
  edm::InputTag rechits_labelEE_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;
  edm::ESGetToken<CaloSubdetectorGeometry, EcalEndcapGeometryRecord> endcapGeomToken_;
  edm::ESGetToken<CaloSubdetectorGeometry, EcalBarrelGeometryRecord> barrelGeomToken_;
  edm::ESGetToken<EcalTrigTowerConstituentsMap, IdealGeometryRecord> eTTmapToken_;
  EcalTPGScale::Tokens tokens_;

  bool recHits_;
};
