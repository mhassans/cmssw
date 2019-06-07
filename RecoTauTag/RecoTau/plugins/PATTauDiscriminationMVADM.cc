#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include <Math/VectorUtil.h>
#include "TMVA/Reader.h"

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
          TString input_name_dm_10_applytoeven = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_dm_10_applytoeven.xml";
          TString input_name_dm_10_applytoodd = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_dm_10_applytoeven.xml";
          TString input_name_dm_0_1_applytoeven = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_dm_0_1_applytoeven.xml";
          TString input_name_dm_0_1_applytoodd = (std::string)getenv("CMSSW_BASE") + "/src/RecoTauTag/RecoTau/TrainingFiles/data/MVADM/mvadm_dm_0_1_applytoodd.xml"; 
          reader_even_ = new TMVA::Reader();
          reader_odd_ = new TMVA::Reader();
          reader_dm10_even_ = new TMVA::Reader();
          reader_dm10_odd_ = new TMVA::Reader();
   
          for(unsigned i=0; i<(unsigned)var_names_.size(); ++i){
            reader_even_->AddVariable( var_names_[i], &(vars_[i]) );
          }
 
        }
    ~PATTauDiscriminationMVADM() override{}
    double discriminate(const TauRef& tau) const override;


  private:

    TMVA::Reader *reader_even_;
    TMVA::Reader *reader_odd_;
    TMVA::Reader *reader_dm10_even_;
    TMVA::Reader *reader_dm10_odd_;

    std::vector<float> vars_; 
    std::vector<float> vars_dm10_;
    //vars_.resize(24);
    //vars_dm10_.resize(40);

    std::vector<TString> var_names_ = {};
    std::vector<TString> var_names_dm10_ = {};

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

  std::pair<std::vector<Vector>, Vector> GetA1 (const TauRef& tau, float gammas_pt_cut) const {
    std::vector<Vector> prongs;
    Vector pi0;
    reco::CandidatePtrVector hads = tau->signalChargedHadrCands();
    if(hads.size()==3) {
      // arrange hadrons so the oppositly charged hadron is contained in the first element
      if(hads[1]->charge()!=hads[0]->charge()&&hads[1]->charge()!=hads[2]->charge()){
        auto temp = hads[1];
        hads[1] = hads[0];
        hads[0] = temp;
      }
      else if(hads[2]->charge()!=hads[0]->charge()&&hads[2]->charge()!=hads[1]->charge()){
        auto temp = hads[2];
        hads[2] = hads[0];
        hads[0] = temp;
      } 
      // from the two same sign hadrons place the one that gives the mass most similar to the rho meson as the second element
      double rho_mass = 0.7755;
      double dM1 = std::fabs((hads[0]->p4()+hads[1]->p4()).M()-rho_mass);
      double dM2 = std::fabs((hads[0]->p4()+hads[2]->p4()).M()-rho_mass);
      if(dM2<dM1){
        auto temp = hads[2];
        hads[2] = hads[1];
        hads[1] = temp;
      }
    }

    reco::CandidatePtrVector gammas;
    for (auto g: tau->isolationGammaCands()) if(g->pt()>gammas_pt_cut) gammas.push_back(g); // change to signal gammas for 94X
    double cone_size = std::max(std::min(0.1, 3./tau->pt()),0.05);
    std::vector<std::pair<Vector, reco::CandidatePtrVector>> strip_pairs = HPSGammas(gammas);
    std::vector<std::pair<Vector, reco::CandidatePtrVector>> strips_incone; 
    for(auto s : strip_pairs) if(std::fabs(ROOT::Math::VectorUtil::DeltaR(s.first,tau->p4()))<cone_size) strips_incone.push_back(s);
     
    reco::CandidatePtrVector signal_gammas = {};
    if(strips_incone.size()>0) {
      pi0 = GetPi0(strips_incone[0].second, true);
      for (auto s : strips_incone) {
        for (auto g : s.second) signal_gammas.push_back(g);
      } 
    } else if(strip_pairs.size()>0) {
      double min_dR = 0.4;
      std::pair<Vector, reco::CandidatePtrVector> closest_strip;
      for (auto s : strip_pairs) {
        double dR = ROOT::Math::VectorUtil::DeltaR(s.first,tau->p4());
        if(dR<min_dR) {
          min_dR = dR;
          closest_strip = s;
        }
      }
      pi0 = GetPi0(closest_strip.second, true);
      for (auto g : closest_strip.second) signal_gammas.push_back(g);
    }
    gammas_ = signal_gammas;
    for (auto h : hads) prongs.push_back((Vector)h->p4());  
    return std::make_pair(prongs, pi0);
  }

 
};


double PATTauDiscriminationMVADM::discriminate(const TauRef& tau) const {
  float gammas_pt_cut = 0.5; // change this for 94X samples
  gammas_.clear();
  // define all variables used by MVA
  float tau_decay_mode = tau->decayMode();

  if (tau_decay_mode>11 || (tau_decay_mode>1&&tau_decay_mode<10)) return 0.;

  Vector pi0;
  Vector pi;
  std::pair<Vector,Vector> rho;
  std::vector<Vector> a1_daughters = {};

  if(tau_decay_mode>1&&tau_decay_mode<10) return -1;

  if(tau_decay_mode>=10) {
    std::pair<std::vector<Vector>, Vector>  a1 = GetA1(tau, gammas_pt_cut);
    a1_daughters  = a1.first;
    pi0 = a1.second;
  } else {
    for (auto g: tau->signalGammaCands()) if(g->pt()>gammas_pt_cut) gammas_.push_back(g);
    rho = GetRho (tau, gammas_pt_cut);
    pi0 = rho.second;
  }
 
  float strip_pt = pi0.pt();
  float E = tau->energy();

  float E1=-1;
  float E2=-1;
  float E3=-1;
  float a1_pi0_dEta=-1;
  float a1_pi0_dphi=-1;
  float a1_pi0_dEta_timesEtau=-1;
  float a1_pi0_dphi_timesEtau=-1;
  float h1_h2_dEta=-1;
  float h1_h2_dphi=-1;
  float h1_h3_dEta=-1;
  float h1_h3_dphi=-1;
  float h2_h3_dEta=-1;
  float h2_h3_dphi=-1;
  float h1_h2_dphi_timesE12=-1;
  float h1_h3_dphi_timesE13=-1;
  float h2_h3_dphi_timesE23=-1;
  float h1_h2_dEta_timesE12=-1;
  float h1_h3_dEta_timesE13=-1;
  float h2_h3_dEta_timesE23=-1;
  float mass0=-1;
  float mass1=-1;
  float mass2=-1;
  float strip_E=-1;
  float E1_overEa1=-1;
  float E2_overEa1=-1;
  float E1_overEtau=-1;
  float E2_overEtau=-1;
  float E3_overEtau=-1;

  if(tau_decay_mode>9) {
    strip_E = pi0.energy();
    mass0 = (a1_daughters[0] + a1_daughters[1] + a1_daughters[2]).M();
    mass1 = (a1_daughters[0] + a1_daughters[1]).M();
    mass2 = (a1_daughters[0] + a1_daughters[2]).M();
    E1 = a1_daughters[0].energy();
    E2 = a1_daughters[1].energy();
    E3 = a1_daughters[2].energy();

    float Ea1 = E1+E2+E3;
    E1_overEa1 = E1/Ea1;
    E2_overEa1 = E2/Ea1;
    float Etau = Ea1+strip_E;
    E1_overEtau = E1/Etau;
    E2_overEtau = E2/Etau;
    E3_overEtau = E3/Etau;

    if(strip_pt>0) {
      a1_pi0_dEta = std::fabs(pi0.eta()-tau->eta());
      a1_pi0_dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(pi0,tau->p4()));
    }

    a1_pi0_dEta_timesEtau=a1_pi0_dEta*Etau;
    a1_pi0_dphi_timesEtau=a1_pi0_dphi*Etau;

    h1_h2_dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(a1_daughters[0],a1_daughters[1]));
    h1_h3_dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(a1_daughters[0],a1_daughters[2]));
    h2_h3_dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(a1_daughters[1],a1_daughters[2]));
    h1_h2_dEta = std::fabs(a1_daughters[0].eta()-a1_daughters[1].eta());
    h1_h3_dEta = std::fabs(a1_daughters[0].eta()-a1_daughters[2].eta());
    h2_h3_dEta = std::fabs(a1_daughters[1].eta()-a1_daughters[2].eta());

    h1_h2_dphi_timesE12=h1_h2_dphi*(E1+E2);
    h1_h3_dphi_timesE13=h1_h3_dphi*(E1+E3);
    h2_h3_dphi_timesE23=h2_h3_dphi*(E2+E3);
    h1_h2_dEta_timesE12=h1_h2_dEta*(E1+E2);
    h1_h3_dEta_timesE13=h1_h3_dEta*(E1+E3);
    h2_h3_dEta_timesE23=h2_h3_dEta*(E2+E3);
  }

  if (tau_decay_mode<12) {
    pi = rho.first; 
  }

  float Egamma1=-1, Egamma2=-1;
  float Epi = pi.energy();
  float Epi0 = pi0.energy();

  if(gammas_.size()>=1) Egamma1 = gammas_[0]->energy();
  if(gammas_.size()>=2) Egamma2 = gammas_[1]->energy();

  float Egamma1_tau = Egamma1/E;
  float Egamma2_tau = Egamma2/E;

  float Epi_tau = Epi/E;

  float pt = tau->pt();
  float eta = tau->eta();

  float rho_dEta=-1, rho_dphi=-1, gammas_dEta = -1., gammas_dphi = -1.;

  if(Epi0>0) {
    rho_dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(pi,pi0));
    rho_dEta = std::fabs(pi.eta()-pi0.eta());
  }
  float rho_dEta_tau = rho_dEta*E;
  float rho_dphi_tau = rho_dphi*E;

  if(gammas_.size()>1) {
    gammas_dphi =  std::fabs(ROOT::Math::VectorUtil::DeltaPhi(gammas_[0]->p4(),gammas_[1]->p4()));
    gammas_dEta =  std::fabs(gammas_[0]->eta()-gammas_[1]->eta());
  }
  float gammas_dEta_tau = gammas_dEta* E;
  float gammas_dR_tau =  sqrt(gammas_dEta*gammas_dEta + gammas_dphi*gammas_dphi)*E;

  float Mpi0=-1, Mpi0_TwoHighGammas=-1;
  Vector gammas_vector;
  for (auto g : gammas_) gammas_vector+=g->p4();
  Mpi0 = gammas_vector.M();
  if(gammas_.size()>=2) Mpi0_TwoHighGammas = (gammas_[0]->p4() + gammas_[1]->p4()).M();

  float Mrho=-1, Mrho_OneHighGammas=-1, Mrho_TwoHighGammas=-1, Mrho_subleadingGamma=-1;
  Mrho = (pi + pi0).M();
  if(gammas_.size()>=1) Mrho_OneHighGammas=(pi + gammas_[0]->p4() ).M();
  if(gammas_.size()>=2) Mrho_TwoHighGammas=(pi + gammas_[0]->p4() + gammas_[1]->p4()).M();
  if(gammas_.size()>=2) Mrho_subleadingGamma=(pi + gammas_[1]->p4()).M();

  float DeltaR2WRTtau=-999;
  if(gammas_.size()>=1){
    DeltaR2WRTtau=0;
    double SumPt=0;
    DeltaR2WRTtau=std::pow(ROOT::Math::VectorUtil::DeltaR(pi,tau->p4()),2)*std::pow(pi.pt(),2);
    SumPt=std::pow(pi.pt(),2);
    for(auto g : gammas_){
      DeltaR2WRTtau+=std::pow(ROOT::Math::VectorUtil::DeltaR(g->p4(),tau->p4()),2)*std::pow(g->pt(),2);
      SumPt+=std::pow(g->pt(),2);
    }
    DeltaR2WRTtau/=SumPt;
  }
  float DeltaR2WRTtau_tau = DeltaR2WRTtau*E*E;

  // once the variables are computed they need to be stored in the order expected by TMVA reader
  std::vector<float> inputs = {};
  if(tau_decay_mode<2) {
    inputs.resize(24);

    inputs[0] = Egamma1_tau;
    inputs[1] = Egamma2_tau;
    inputs[2] = Epi_tau;
    inputs[3] = rho_dEta_tau;
    inputs[4] = rho_dphi_tau;
    inputs[5] = gammas_dEta_tau;
    inputs[6] = gammas_dR_tau;
    inputs[7] = DeltaR2WRTtau_tau;
    inputs[8] = tau_decay_mode;
    inputs[9] = eta;
    inputs[10] = pt;
    inputs[11] = Epi0;
    inputs[12] = Epi;
    inputs[13] = rho_dEta;
    inputs[14] = rho_dphi;
    inputs[15] = gammas_dEta;
    inputs[16] = Mrho;
    inputs[17] = Mpi0;
    inputs[18] = DeltaR2WRTtau;
    inputs[19] = Mpi0_TwoHighGammas;
    inputs[20] = Mrho_OneHighGammas;
    inputs[21] = Mrho_TwoHighGammas;
    inputs[22] = Mrho_subleadingGamma;
    inputs[23] = strip_pt;
  }
  if(tau_decay_mode>9) {
    inputs.resize(40);

    inputs[0] = E1_overEa1;
    inputs[1] = E2_overEa1;
    inputs[2] = E1_overEtau;
    inputs[3] = E2_overEtau;
    inputs[4] = E3_overEtau;
    inputs[5] = a1_pi0_dEta_timesEtau;
    inputs[6] = a1_pi0_dphi_timesEtau;
    inputs[7] = h1_h2_dphi_timesE12;
    inputs[8] = h1_h2_dEta_timesE12;
    inputs[9] = h1_h3_dphi_timesE13;
    inputs[10] = h1_h3_dEta_timesE13;
    inputs[11] = h2_h3_dphi_timesE23;
    inputs[12] = h2_h3_dEta_timesE23;
    inputs[13] = gammas_dEta_tau;
    inputs[14] = gammas_dR_tau;
    inputs[15] = tau_decay_mode;
    inputs[16] = mass0;
    inputs[17] = mass1;
    inputs[18] = mass2;
    inputs[19] = E1;
    inputs[20] = E2;
    inputs[21] = E3;
    inputs[22] = strip_E;
    inputs[23] = a1_pi0_dEta;
    inputs[24] = a1_pi0_dphi;
    inputs[25] = strip_pt;
    inputs[26] = pt;
    inputs[27] = eta;
    inputs[28] = E;
    inputs[29] = h1_h2_dphi;
    inputs[30] = h1_h3_dphi;
    inputs[31] = h2_h3_dphi;
    inputs[32] = h1_h2_dEta;
    inputs[33] = h1_h3_dEta;
    inputs[34] = h2_h3_dEta;
    inputs[35] = Egamma1;
    inputs[36] = Egamma2;
    inputs[37] = gammas_dEta;
    inputs[38] = Mpi0;
    inputs[39] = Mpi0_TwoHighGammas;
  }

  double mvaScore = -1;
  if((targetDM_<=2&&targetDM_>=0) || targetDM_==10 || targetDM_==11) mvaScore = targetDM_;
  return mvaScore;
}

}
DEFINE_FWK_MODULE(PATTauDiscriminationMVADM);
