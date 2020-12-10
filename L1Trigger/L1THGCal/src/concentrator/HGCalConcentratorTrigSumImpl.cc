#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorTrigSumImpl.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToModule.h"

HGCalConcentratorTrigSumImpl::HGCalConcentratorTrigSumImpl(const edm::ParameterSet& conf) {}

void HGCalConcentratorTrigSumImpl::doSum(uint32_t module_id,
                                         const std::vector<l1t::HGCalTriggerCell>& trigCellVecInput,
                                         std::vector<l1t::HGCalTriggerSums>& trigSumsVecOutput) const {
  double ptsum = 0;
  double mipptsum = 0;
  double hwptsum = 0;

  for (const auto& trigCell : trigCellVecInput) {
    // detId selection is already done in HGCalConcentratorProcessorSelection:
    // here we do not worry about it and assume all cells are from the same module
    ptsum += trigCell.pt();
    mipptsum += trigCell.mipPt();
    hwptsum += trigCell.hwPt();
  }
  if (!trigCellVecInput.empty()) {
    uint32_t id0 = trigCellVecInput[0].detId();
    uint32_t tsid = module_id;
    if(triggerTools_.isSilicon(id0)) {
      HGCalTriggerDetId id0si(id0);
      DetId::Detector det = (id0si.subdet() == HGCalTriggerSubdetector::HGCalEETrigger) ? DetId::HGCalEE : DetId::HGCalHSi;
      tsid = HGCSiliconDetId(det, id0si.zside(), id0si.type(), id0si.layer(), id0si.waferU(), id0si.waferV(), 0, 0);
    }
    GlobalPoint module_pos = triggerTools_.getTriggerGeometry()->getModulePosition(module_id);

    math::PtEtaPhiMLorentzVector p4(ptsum, module_pos.eta(), module_pos.phi(), 0);
    l1t::HGCalTriggerSums ts;
    ts.setP4(p4);
    ts.setDetId(tsid);
    ts.setPosition(module_pos);
    ts.setMipPt(mipptsum);
    ts.setHwPt(hwptsum);
    trigSumsVecOutput.push_back(ts);
  }
}
