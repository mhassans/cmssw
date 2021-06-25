
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTowerGeometryHelper.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

HGCalTriggerTowerGeometryHelper::HGCalTriggerTowerGeometryHelper(const edm::ParameterSet& conf)
    : doNose_(conf.getParameter<bool>("doNose")),
      minEta_(conf.getParameter<double>("minEta")),
      maxEta_(conf.getParameter<double>("maxEta")),
      minPhi_(conf.getParameter<double>("minPhi")),
      maxPhi_(conf.getParameter<double>("maxPhi")),
      nBinsEta_(conf.getParameter<int>("nBinsEta")),
      nBinsPhi_(conf.getParameter<int>("nBinsPhi")),
      binsEta_(conf.getParameter<std::vector<double> >("binsEta")),
      binsPhi_(conf.getParameter<std::vector<double> >("binsPhi")),
      splitModuleSum_(conf.getParameter<bool>("splitModuleSum")),
      moduleTowerMapping_(conf.getParameter<edm::FileInPath>("moduleTowerMapping")),
      splitDivisorSilic_(conf.getParameter<int>("splitDivisorSilic")),
      splitDivisorScint_(conf.getParameter<int>("splitDivisorScint")){
  if (!binsEta_.empty() && ((unsigned int)(binsEta_.size()) != nBinsEta_ + 1)) {
    throw edm::Exception(edm::errors::Configuration, "Configuration")
        << "HGCalTriggerTowerGeometryHelper nBinsEta for the tower map not consistent with binsEta size" << std::endl;
  }

  if (!binsPhi_.empty() && ((unsigned int)(binsPhi_.size()) != nBinsPhi_ + 1)) {
    throw edm::Exception(edm::errors::Configuration, "Configuration")
        << "HGCalTriggerTowerGeometryHelper nBinsPhi for the tower map not consistent with binsPhi size" << std::endl;
  }

  // if the bin vecctor is empty we assume the bins to be regularly spaced
  if (binsEta_.empty()) {
    for (unsigned int bin1 = 0; bin1 != nBinsEta_ + 1; bin1++) {
      binsEta_.push_back(minEta_ + bin1 * ((maxEta_ - minEta_) / nBinsEta_));
    }
  }

  // if the bin vecctor is empty we assume the bins to be regularly spaced
  if (binsPhi_.empty()) {
    for (unsigned int bin2 = 0; bin2 != nBinsPhi_ + 1; bin2++) {
      binsPhi_.push_back(minPhi_ + bin2 * ((maxPhi_ - minPhi_) / nBinsPhi_));
    }
  }

  for (int zside = -1; zside <= 1; zside += 2) {
    for (unsigned int bin1 = 0; bin1 != nBinsEta_; bin1++) {
      for (unsigned int bin2 = 0; bin2 != nBinsPhi_; bin2++) {
        l1t::HGCalTowerID towerId(doNose_, zside, bin1, bin2);
        tower_coords_.emplace_back(towerId.rawId(),
                                   zside * ((binsEta_[bin1 + 1] + binsEta_[bin1]) / 2),
                                   (binsPhi_[bin2 + 1] + binsPhi_[bin2]) / 2);
      }
    }
  }

  if (conf.getParameter<bool>("readMappingFile")) {
    // We read the TC to TT mapping from file,
    // otherwise we derive the TC to TT mapping on the fly from eta-phi coord. of the TCs
    std::ifstream l1tTriggerTowerMappingStream(conf.getParameter<edm::FileInPath>("L1TTriggerTowerMapping").fullPath());
    if (!l1tTriggerTowerMappingStream.is_open()) {
      throw cms::Exception("MissingDataFile") << "Cannot open HGCalTriggerGeometry L1TTriggerTowerMapping file\n";
    }

    unsigned trigger_cell_id = 0;
    unsigned short iEta = 0;
    unsigned short iPhi = 0;

    for (; l1tTriggerTowerMappingStream >> trigger_cell_id >> iEta >> iPhi;) {
      if (iEta >= nBinsEta_ || iPhi >= nBinsPhi_) {
        throw edm::Exception(edm::errors::Configuration, "Configuration")
            << "HGCalTriggerTowerGeometryHelper warning inconsistent mapping TC : " << trigger_cell_id
            << " to TT iEta: " << iEta << " iPhi: " << iPhi << " when max #bins eta: " << nBinsEta_
            << " phi: " << nBinsPhi_ << std::endl;
      }
      l1t::HGCalTowerID towerId(doNose_, triggerTools_.zside(DetId(trigger_cell_id)), iEta, iPhi);
      cells_to_trigger_towers_[trigger_cell_id] = towerId.rawId();
    }
    l1tTriggerTowerMappingStream.close();
  }
}

const std::vector<l1t::HGCalTowerCoord>& HGCalTriggerTowerGeometryHelper::getTowerCoordinates() const {
  return tower_coords_;
}

unsigned short HGCalTriggerTowerGeometryHelper::getTriggerTowerFromEtaPhi(const float& eta, const float& phi) const {
  auto bin_eta_l = std::lower_bound(binsEta_.begin(), binsEta_.end(), fabs(eta));
  unsigned int bin_eta = 0;
  // we add a protection for TCs in Hadron part which are outside the boundaries and possible rounding effects
  if (bin_eta_l == binsEta_.end()) {
    if (fabs(eta) < minEta_) {
      bin_eta = 0;
    } else if (fabs(eta) >= maxEta_) {
      bin_eta = nBinsEta_;
    } else {
      edm::LogError("HGCalTriggerTowerGeometryHelper")
          << " did not manage to map eta " << eta << " to any Trigger Tower\n";
    }
  } else {
    bin_eta = bin_eta_l - binsEta_.begin() - 1;
  }

  auto bin_phi_l = std::lower_bound(binsPhi_.begin(), binsPhi_.end(), phi);
  unsigned int bin_phi = 0;
  if (bin_phi_l == binsPhi_.end()) {
    if (phi < minPhi_) {
      bin_phi = nBinsPhi_;
    } else if (phi >= maxPhi_) {
      bin_phi = 0;
    } else {
      edm::LogError("HGCalTriggerTowerGeometryHelper")
          << " did not manage to map phi " << phi << " to any Trigger Tower\n";
    }
  } else {
    bin_phi = bin_phi_l - binsPhi_.begin() - 1;
  }
  int zside = eta < 0 ? -1 : 1;

  return l1t::HGCalTowerID(doNose_, zside, bin_eta, bin_phi).rawId();
}

std::unordered_map<unsigned short, float> HGCalTriggerTowerGeometryHelper::getTriggerTower(const l1t::HGCalTriggerCell& thecell) const {
  std::unordered_map<unsigned short, float> binIDandShares = {}; 
  unsigned int trigger_cell_id = thecell.detId();
  // NOTE: if the TC is not found in the map than it is mapped via eta-phi coords.
  // this can be considered dangerous (silent failure of the map) but it actually allows to save
  // memory mapping explicitly only what is actually needed
  auto tower_id_itr = cells_to_trigger_towers_.find(trigger_cell_id);
  if (tower_id_itr != cells_to_trigger_towers_.end()){
    binIDandShares.insert({tower_id_itr->second, 1.0});
    return binIDandShares;
  }
  binIDandShares.insert({getTriggerTowerFromEtaPhi(thecell.position().eta(), thecell.position().phi()), 1.0});
  return binIDandShares;
}

std::unordered_map<unsigned short, float> HGCalTriggerTowerGeometryHelper::getTriggerTower(const l1t::HGCalTriggerSums& thesum) const {
  std::unordered_map<unsigned short, float> binIDandShares = {}; 
  if(!splitModuleSum_){
    binIDandShares.insert({getTriggerTowerFromEtaPhi(thesum.position().eta(), thesum.position().phi()), 1.0});
    return binIDandShares;
  } else{
    HGCalModuleDetId detid(thesum.detId());
    int moduleU = detid.moduleU();
    int moduleV = detid.moduleV();
    int layer = detid.layer();
    int sector = detid.sector();
    int zside = detid.zside();

    int subdet = -1;
    int splitDivisor = -1;
    if(detid.isHScintillator()){ //FIXME HFNose not handled
      subdet = 2;
      layer += 28;
      splitDivisor = splitDivisorScint_; 
    } else if(detid.isEE()){
      subdet = 0;
      splitDivisor = splitDivisorSilic_;
    } else if(detid.isHSilicon()) {
      subdet = 1;
      layer += 28;
      splitDivisor = splitDivisorSilic_;
    }
    
    std::ifstream moduleTowerMappingStream(moduleTowerMapping_.fullPath());
    if (!moduleTowerMappingStream.is_open()) {
      throw cms::Exception("MissingDataFile") << "Cannot open HGCalTowerMapProducer moduleTowerMapping file\n";
    }

    std::vector<std::string> result; // the line including towers and module sum shares for this module
    std::string line;
    getline(moduleTowerMappingStream, line); //To skip first row
    for( std::string line; getline(moduleTowerMappingStream, line ); ){
      std::stringstream ss(line);
      while(ss.good()){
        std::string substr;
        std::getline(ss, substr, ' ');
        result.push_back(substr);
      }
      if ((std::stoi(result.at(0))==subdet) && (std::stoi(result.at(1))==layer) && (std::stoi(result.at(2))==moduleU) && (std::stoi(result.at(3))==moduleV)){
        break;
      } else{
        result.clear();
      }
    }
    if(result.empty()){ // for modules not found in the mapping file (currently a few partial modules) use the traditional method.
      binIDandShares.insert({getTriggerTowerFromEtaPhi(thesum.position().eta(), thesum.position().phi()), 1.0});
      return binIDandShares;
    }

    int towerEta;
    int towerPhi;
    for (int i=1; i<=std::stoi(result.at(4)); i++){
      towerEta = 2 + std::stoi(result.at(3*i+2)); // shift by two to avoid negative eta
      towerPhi = (std::stoi(result.at(3*i+3)) + sector*int(nBinsPhi_)/3 + int(nBinsPhi_)/2) % int(nBinsPhi_);//move to the correct sector
      if(zside==1){
        towerPhi = (3*int(nBinsPhi_)/2 - towerPhi - 1) % int(nBinsPhi_); //correct x -> -x in z>0
      }
      binIDandShares.insert( {l1t::HGCalTowerID(doNose_, zside, towerEta, towerPhi).rawId(),  std::stod(result.at(3*i+4))/splitDivisor } );
    }
    return binIDandShares; 
  }
}
