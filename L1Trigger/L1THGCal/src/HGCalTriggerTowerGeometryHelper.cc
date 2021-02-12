
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
    : minEta_(conf.getParameter<double>("minEta")),
      maxEta_(conf.getParameter<double>("maxEta")),
      minPhi_(conf.getParameter<double>("minPhi")),
      maxPhi_(conf.getParameter<double>("maxPhi")),
      nBinsEta_(conf.getParameter<int>("nBinsEta")),
      nBinsPhi_(conf.getParameter<int>("nBinsPhi")),
      binsEta_(conf.getParameter<std::vector<double> >("binsEta")),
      binsPhi_(conf.getParameter<std::vector<double> >("binsPhi")) {
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
        l1t::HGCalTowerID towerId(zside, bin1, bin2);
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
      l1t::HGCalTowerID towerId(triggerTools_.zside(DetId(trigger_cell_id)), iEta, iPhi);
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

//TESSSSSSSSSSSSSSSSSSSTTTTTTTTTTTTTTTTTTTTTTTTT
    std::cout<<"*******"<<std::endl;
    std::cout<<"eta="<<eta<<",  phi="<<phi<<std::endl;
    std::cout<<"thierbin_eta="<<bin_eta<<",  theirbin_phi="<<bin_phi<<",  rawID="<<l1t::HGCalTowerID(zside, bin_eta, bin_phi).rawId()<<std::endl;
    std::cout<<"zside="<<zside<<std::endl;
    std::cout<<"*******"<<std::endl;
//TESSSSSSSSSSSSSSSSSSSTTTTTTTTTTTTTTTTTTTTTTTTT

  return l1t::HGCalTowerID(zside, bin_eta, bin_phi).rawId();
}

unsigned short HGCalTriggerTowerGeometryHelper::getTriggerTower(const l1t::HGCalTriggerCell& thecell) const {
  unsigned int trigger_cell_id = thecell.detId();
  // NOTE: if the TC is not found in the map than it is mapped via eta-phi coords.
  // this can be considered dangerous (silent failure of the map) but it actually allows to save
  // memory mapping explicitly only what is actually needed
  auto tower_id_itr = cells_to_trigger_towers_.find(trigger_cell_id);
  if (tower_id_itr != cells_to_trigger_towers_.end())
    return tower_id_itr->second;

  return getTriggerTowerFromEtaPhi(thecell.position().eta(), thecell.position().phi());
}

//Original function
//unsigned short HGCalTriggerTowerGeometryHelper::getTriggerTower(const l1t::HGCalTriggerSums& thesum) const {
//  return getTriggerTowerFromEtaPhi(thesum.position().eta(), thesum.position().phi());
//}

//unsigned HGCalTriggerTowerGeometryHelper::uvMapping(unsigned layer, std::pair<int,int> &uv) {
//  unsigned sector(0);
//  int offset;
//  
//  if(layer<=28) { // CE-E    
//    if(uv.first>0 && uv.second>=0) return sector;
//    
//    offset=0;
//    if(uv.first>=uv.second && uv.second<0) sector=2;
//    else sector=1;
//    
//  } else if((layer%2)==1) { // CE-H Odd
//    if(uv.first>=0 && uv.second>=0) return sector;
//    
//    offset=-1;    
//    if(uv.first>uv.second && uv.second<0) sector=2;
//    else sector=1;
//    
//  } else { // CE-H Even
//    if(uv.first>=1 && uv.second>=1) return sector;
//    
//    offset=1;
//    if(uv.first>=uv.second && uv.second<1) sector=2;
//    else sector=1;
//  }
//  
//  int up,vp;
//  
//  if(sector==1) {
//    up=uv.second-uv.first;
//    vp=-uv.first+offset;    
//    
//  } else {
//    up=-uv.second+offset;
//    vp=uv.first-uv.second+offset;
//  }
//  
//  uv.first=up;
//  uv.second=vp;
//  return sector;
//}


unsigned short HGCalTriggerTowerGeometryHelper::getTriggerTower(const l1t::HGCalTriggerSums& thesum) const {
  HGCSiliconDetId detid(thesum.detId());
  int waferu = detid.waferU();
  int waferv = detid.waferV();
  int layer = detid.layer();
  
  if(triggerTools_.isSilicon(detid) && waferu>=1 && waferv>=1){//if u or v non-negative, etc.  
 //   std::ifstream input("inputFile/tower_per_module.txt");
 //   std::vector<std::string> result;
 //   for( std::string line; getline( input, line ); ){
 //       std::stringstream ss(line);
 //       while(ss.good()){
 //           std::string substr;
 //           std::getline(ss, substr, ' ');
 //           result.push_back(substr);
 //      }

 //      if (std::stoi(result[0])==layer && std::stoi(result[1])==waferu && std::stoi(result[2])==waferv){
 //         break;
 //      }
 //      else{
 //         result.clear();
 //      }
 //   }
 //   
 //   for (int i=1; i<=std::stoi(result[3]); i++){
 //       int zside = thesum.position().eta() < 0 ? -1 : 1;
 //       std::cout<<"iter="<<i<<", myeta="<<std::stoi(result[3*i+1])<<", myphi="<<std::stoi(result[3*i+2])<<", frac="<<std::stoi(result[3*i+3])<<std::endl;
 //       std::cout<<"myTT:"<<l1t::HGCalTowerID( zside, std::stoi(result[3*i+1]), std::stoi(result[3*i+2]) ).rawId()<<std::endl;
 //       std::cout<<"theirTT:"<<getTriggerTowerFromEtaPhi(thesum.position().eta(), thesum.position().phi());
 //       std::cout<<"------------------"<<std::endl;
 //   }
 //   std::cout<<"++++++++++FINISH+++++++++++++"<<std::endl;
    return 0;//getTriggerTowerFromEtaPhi(thesum.position().eta(), thesum.position().phi());
  }
  else{
    return 0;
  }
  
  //std::cout<<", waferu="<<waferu<<",  waferv="<<waferv<<",  layer="<<layer<<std::endl;
  //std::cout<<"silicon="<<triggerTools_.isSilicon(detid)<<std::endl;
  
 // if(layer==5 && waferu==-4 && waferv==1){
 // std::cout<<"*******"<<std::endl;
 // std::cout<<""<<std::endl;
 // std::cout<<"*******"<<std::endl;
 // }
  
  //return getTriggerTowerFromEtaPhi(thesum.position().eta(), thesum.position().phi());
}

//std::unordered_map<unsigned short, unsigned short> HGCalTriggerTowerGeometryHelper::getTriggerTowerSplit(const l1t::HGCalTriggerSums& thesum){}
