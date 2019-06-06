#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include <Math/VectorUtil.h>

/* class PATTauDiscriminationMVADM
 *
 *  Seperate tau decay types using MVA
 *  Returns MVA score for each class
 *
 */

namespace {

template <class T, class U>
bool sortStrips (std::pair<T,U> i, std::pair<T,U> j) {
  return (i.first.pt() > j.first.pt());
}



class PATTauDiscriminationMVADM final : public PATTauDiscriminationProducerBase  {
  public:
    explicit PATTauDiscriminationMVADM(const edm::ParameterSet& iConfig)
        :PATTauDiscriminationProducerBase(iConfig){
          targetDM_        = iConfig.getParameter<int>("targetDM");
        }
    ~PATTauDiscriminationMVADM() override{}
    double discriminate(const TauRef& tau) const override;


  private:
    int targetDM_;
    mutable reco::CandidatePtrVector gammas_;

    typedef ROOT::Math::PtEtaPhiEVector Vector;

    Vector GetPi0 (reco::CandidatePtrVector gammas, bool leadEtaPhi) const {
      Vector pi0;
      if(gammas.size()>0) {
        double E = 0.;
        double phi = 0.;
        double eta = 0.;
        for(auto g: gammas) {
          E+=g->energy();
          phi+=g->energy()*g->phi();
          eta+=g->energy()*g->eta();
        }
        eta/=E;
        phi/=E;
  
        if(leadEtaPhi){
          // if true sets the eta and phi of the pi0 to that of the leading gamma rather than using the weighted average
          eta = gammas[0]->eta();
          phi = gammas[0]->phi();
        }
  
        double mass = 0.1349;
        double p = sqrt(E*E-mass*mass);
        double theta = atan(exp(-eta))*2;
        double pt = p*sin(theta);
        pi0 = Vector(pt,eta,phi,E);
      }
      return pi0;
    }


    std::vector<std::pair<Vector,reco::CandidatePtrVector>> HPSGammas (reco::CandidatePtrVector cands) const {
      std::vector<std::pair<Vector,reco::CandidatePtrVector>> strips;   
      while(!cands.empty()) {
  
        reco::CandidatePtrVector Associated = {};
        reco::CandidatePtrVector notAssociated = {};
  
        Vector stripVector(0,0,0,0);
        stripVector=cands[0]->p4();
        Associated.push_back(cands[0]);
 
        bool repeat = true;
        while (repeat) {
          repeat = false;
          for(unsigned int i=1;i<cands.size();++i) {
            double etaAssociationDistance = 0.20*pow(cands[i]->pt(),-0.66) + 0.20*pow(stripVector.Pt(),-0.66);
            double phiAssociationDistance = 0.35*pow(cands[i]->pt(),-0.71) + 0.35*pow(stripVector.Pt(),-0.71);
            etaAssociationDistance = std::min(etaAssociationDistance, 0.15);
            etaAssociationDistance = std::max(etaAssociationDistance, 0.05);
            phiAssociationDistance = std::min(phiAssociationDistance, 0.30);
            phiAssociationDistance = std::max(phiAssociationDistance, 0.05);
  
            if(fabs(cands[i]->eta()-stripVector.eta())<etaAssociationDistance &&
              fabs(ROOT::Math::VectorUtil::DeltaPhi(cands[i]->p4(),stripVector))<phiAssociationDistance) {
              stripVector+=cands[i]->p4();
              Associated.push_back(cands[i]);
              repeat = true;
            }
            else {
              notAssociated.push_back(cands[i]);
            }
          }
          cands.swap(notAssociated);
          notAssociated.clear(); 
        }

        Vector strip = GetPi0(Associated, false);
        strips.push_back(std::make_pair(strip, Associated));
 
      }
      std::sort(strips.begin(), strips.end(), sortStrips<Vector, reco::CandidatePtrVector>);
  
      return strips;
    }


    std::pair<Vector,Vector> GetRho (const TauRef& tau, double gammas_pt_cut) const {
      Vector pi;
      Vector pi0;
      gammas_.clear();

      reco::CandidatePtrVector gammas;
      for (auto g: tau->signalGammaCands()) if(g->pt()>gammas_pt_cut) gammas.push_back(g);
      reco::CandidatePtrVector hads = tau->signalChargedHadrCands();

      if(hads.size()>0) pi = hads[0]->p4();

      double cone_size = std::max(std::min(0.1, 3./tau->pt()),0.05);
      std::vector<std::pair<Vector, reco::CandidatePtrVector>> strip_pairs = HPSGammas(gammas);
      std::vector<std::pair<Vector, reco::CandidatePtrVector>> strips_incone;
      for(auto s : strip_pairs) {
        if(std::fabs(ROOT::Math::VectorUtil::DeltaR(s.first,tau->p4()))<cone_size) strips_incone.push_back(s);
      }
      if(tau->decayMode()==0) {
        if(strips_incone.size()>0) {
          gammas = strips_incone[0].second;
          for (auto s : strips_incone) {
            for (auto g : s.second) gammas_.push_back(g);
          }
        } else if(strip_pairs.size()>0) {
          double min_dR = 0.4;
          std::pair<Vector,reco::CandidatePtrVector> closest_strip;
          for (auto s : strip_pairs) {
            double dR = ROOT::Math::VectorUtil::DeltaR(s.first,tau->p4());
            if(dR<min_dR) {
              min_dR = dR;
              closest_strip = s;
            }
          }
          gammas = closest_strip.second;
          gammas_ = gammas;
        }
      } else {
        gammas_ = gammas;
      }
      pi0 = GetPi0(gammas, true);
    
      return std::make_pair(pi,pi0);
    }


//    float Egamma1_tau_;
//    float Egamma2_tau_;
//    float Epi_tau_;
//    float rho_dEta_tau_;
//    float rho_dphi_tau_;
//    float gammas_dEta_tau_;
//    float gammas_dR_tau_;
//    float DeltaR2WRTtau_tau_;
//    float eta_;
//    float pt_;
//    float Epi0_;
//    float Epi_;
//    float rho_dEta_;
//    float rho_dphi_;
//    float gammas_dEta_;
//    float tau_decay_mode_;
//    float Mrho_;
//    float Mpi0_;
//    float DeltaR2WRTtau_;
//    float Mpi0_TwoHighGammas_;
//    float Mrho_OneHighGammas_;
//    float Mrho_TwoHighGammas_;
//    float Mrho_subleadingGamma_;
//    float strip_pt_; 
//    float E_;
//    float E1_;
//    float E1_overEa1_;
//    float E1_overEtau_;
//    float E2_;
//    float E2_overEa1_;
//    float E2_overEtau_;
//    float E3_;
//    float E3_overEtau_;
//    float a1_pi0_dEta_;
//    float a1_pi0_dEta_timesEtau_;
//    float a1_pi0_dphi_;
//    float a1_pi0_dphi_timesEtau_;
//    float h1_h2_dEta_;
//    float h1_h2_dEta_timesE12_;
//    float h1_h2_dphi_;
//    float h1_h2_dphi_timesE12_;
//    float h1_h3_dEta_;
//    float h1_h3_dEta_timesE13_;
//    float h1_h3_dphi_;
//    float h1_h3_dphi_timesE13_;
//    float h2_h3_dEta_;
//    float h2_h3_dEta_timesE23_;
//    float h2_h3_dphi_;
//    float h2_h3_dphi_timesE23_;
//    float mass0_;
//    float mass1_;
//    float mass2_;
//    float strip_E_;
 
};


double PATTauDiscriminationMVADM::discriminate(const TauRef& tau) const {
  float gammas_pt_cut = 0.5; // change this for 94X samples
  gammas_.clear();
  // define all variables used by MVA
  float tau_decay_mode = tau->decayMode();

  if(tau_decay_mode>1&&tau_decay_mode<10) return -1;

  if(tau_decay_mode>=10) {
    ;
  } else {
    for (auto g: tau->signalGammaCands()) if(g->pt()>gammas_pt_cut) gammas_.push_back(g);
    std::pair<Vector,Vector> rho = GetRho (tau, gammas_pt_cut);
    if(rho.first.mass()>0) ;
  }


        std::vector<float> inputs = {};

        if(tau_decay_mode<2) {
          inputs.resize(24);

          //inputs[0] = Egamma1_tau;
          //inputs[1] = Egamma2_tau;
          //inputs[2] = Epi_tau;
          //inputs[3] = rho_dEta_tau;
          //inputs[4] = rho_dphi_tau;
          //inputs[5] = gammas_dEta_tau;
          //inputs[6] = gammas_dR_tau;
          //inputs[7] = DeltaR2WRTtau_tau;
          inputs[8] = tau_decay_mode;
          //inputs[9] = eta;
          //inputs[10] = pt;
          //inputs[11] = Epi0;
          //inputs[12] = Epi;
          //inputs[13] = rho_dEta;
          //inputs[14] = rho_dphi;
          //inputs[15] = gammas_dEta;
          //inputs[16] = Mrho;
          //inputs[17] = Mpi0;
          //inputs[18] = DeltaR2WRTtau;
          //inputs[19] = Mpi0_TwoHighGammas;
          //inputs[20] = Mrho_OneHighGammas;
          //inputs[21] = Mrho_TwoHighGammas;
          //inputs[22] = Mrho_subleadingGamma;
          //inputs[23] = strip_pt;
        }


  double mvaScore = -1;
  if((targetDM_<=2&&targetDM_>=0) || targetDM_==10 || targetDM_==11) mvaScore = targetDM_;
  return mvaScore;
}

}
DEFINE_FWK_MODULE(PATTauDiscriminationMVADM);
