#pragma once

// ROOT includes
#include "TVector3.h"
#include "TLorentzVector.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Constants.hh"

class KICalculator {

public:

  enum CalcType { kOpt1, kOpt2, kOpt3, kOpt4, nOpts };

  KICalculator( const TVector3& p3_lep, const TVector3& p3_had,
    const CalcType calc_opt = kOpt1, const double m_lep = MUON_MASS,
    const double m_had = PROTON_MASS );

  inline const TVector3& p3_lep() const { return p3_lep_; }
  inline const TVector3& p3_had() const { return p3_had_; }
  inline const double m_lep() const { return m_lep_; }
  inline const double m_had() const { return m_had_; }
  inline const CalcType calc_opt() const { return calc_opt_; }

  inline double k_miss() const { return k_miss_; }
  inline double E_miss() const { return E_miss_; }
  inline double p_miss_minus() const { return p_miss_minus_; }
  inline double p_miss() const { return p_miss_; }
  inline double pT() const { return pT_; }
  inline double pL() const { return pL_; }
  inline double pn() const { return pn_; }
  inline double delta_alphaT() const { return delta_alphaT_; }
  inline double delta_alpha3D_q() const { return delta_alpha3D_q_; }
  inline double delta_alpha3D_mu() const { return delta_alpha3D_mu_; }
  inline double delta_phiT() const { return delta_phiT_; }
  inline double delta_phi3D() const { return delta_phi3D_; }
  inline double Ecal() const { return Ecal_; }
  inline double Ecal_MB() const { return Ecal_MB_; }
  inline double E_QE() const { return E_QE_; }
  inline double Q2() const { return Q2_; }
  inline double A() const { return A_; }
  inline double pTx() const { return pTx_; }
  inline double pTy() const { return pTy_; }
  inline double pn_perp() const { return pn_perp_; }
  inline double pn_perpx() const { return pn_perpx_; }
  inline double pn_perpy() const { return pn_perpy_; }
  inline double pn_par() const { return pn_par_; }
  inline double theta_lep_had() const { return theta_lep_had_; }

protected:

  // Store copies of the inputs used by the constructor
  const TVector3 p3_lep_;
  const TVector3 p3_had_;
  const double m_lep_;
  const double m_had_;
  const CalcType calc_opt_;

  double k_miss_;
  double E_miss_;
  double p_miss_minus_;
  double p_miss_;
  double pT_;
  double pL_;
  double pn_;
  double delta_alphaT_;
  double delta_alpha3D_q_;
  double delta_alpha3D_mu_;
  double delta_phiT_;
  double delta_phi3D_;
  double Ecal_;
  double Ecal_MB_;
  double E_QE_;
  double Q2_;
  double A_;
  double pTx_;
  double pTy_;
  double pn_perp_;
  double pn_perpx_;
  double pn_perpy_;
  double pn_par_;
  double theta_lep_had_;
};
