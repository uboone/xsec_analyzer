// Afroditi Papadopoulou <afropapp13@outlook.com>
// Daniel Barrow <daniel.barrow@physics.ox.ac.uk>
// Steven Gardiner <gardiner@fnal.gov>
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/KICalculator.hh"

KICalculator::KICalculator( const TVector3& p3_lep, const TVector3& p3_had,
  const KICalculator::CalcType calc_opt, const double m_lep,
  const double m_had ) : p3_lep_( p3_lep ), p3_had_( p3_had ), m_lep_( m_lep ),
  m_had_( m_had ), calc_opt_( calc_opt )
{
  double E_lep = real_sqrt( p3_lep.Mag2() + m_lep*m_lep );
  double E_had = real_sqrt( p3_had.Mag2() + m_had*m_had );

  double binding_energy;

  switch ( calc_opt ) {
  case kOpt1:
  case kOpt4:
    // This value is the shell-occupancy-weighted mean of the $E_{\alpha}$
    // values listed for 40Ar in Table II of arXiv:1609.03530. MINERvA uses an
    // identical procedure for 12C to obtain the binding energy value of 27.13
    // MeV, which is adopted in their STV analysis described in
    // arXiv:1910.08658
    binding_energy = 0.02478; // GeV
    break;
  case kOpt2:
    // For the calculation of the excitation energies
    // https://doi.org/10.1140/epjc/s10052-019-6750-3
    binding_energy = 0.0309; // GeV
    break;
  case kOpt3:
    // https://tinyurl.com/binding-stv-tools
    binding_energy = 0.04; // GeV
  default:
    throw std::runtime_error( "Unrecognized option for KICalculator "
      + std::to_string( calc_opt_ ) );
  }

  double m2_diff_np = NEUTRON_MASS*NEUTRON_MASS - PROTON_MASS*PROTON_MASS;

  TVector3 p3_lep_T( p3_lep.X(), p3_lep.Y(), 0.);

  TLorentzVector p4_lep( p3_lep, E_lep );

  TVector3 p3_had_T( p3_had.X(), p3_had.Y(), 0. );

  TLorentzVector p4_had( p3_had, E_had );
  double had_KE = E_had - m_had;

  TVector3 p3_T = p3_lep_T + p3_had_T;

  pT_ = p3_T.Mag();
  TVector2 p2_T = ( p3_lep + p3_had ).XYvector();

  delta_alphaT_ = std::acos( ( -1. * p3_lep_T * p3_T )
    / ( p3_lep_T.Mag() * pT_ ) ) * 180. / M_PI;
  if ( delta_alphaT_ > 180. ) delta_alphaT_ -= 180.;
  if ( delta_alphaT_ < 0. ) delta_alphaT_ += 180.;

  delta_phiT_ = std::acos( (- p3_lep_T * p3_had_T)
    / ( p3_lep_T.Mag() * p3_had_T.Mag() ) ) * 180. / M_PI;
  if ( delta_phiT_ > 180. ) { delta_phiT_ -= 180.; }
  if ( delta_phiT_ < 0. ) { delta_phiT_ += 180.; }

  // Calorimetric Energy Reconstruction
  Ecal_ = E_lep + had_KE + binding_energy; // GeV

  // QE Energy Reconstruction

  double E_QE_numer = 2.*( NEUTRON_MASS - binding_energy )*E_lep
    - ( binding_energy*binding_energy - 2.*NEUTRON_MASS*binding_energy
    + MUON_MASS*MUON_MASS + m2_diff_np );
  double E_QE_denom = 2.*( NEUTRON_MASS - binding_energy - E_lep
    + p3_lep.Mag() * p3_lep.CosTheta() );
  E_QE_ = E_QE_numer / E_QE_denom;

  // Reconstructed Q2
  TLorentzVector p4_nu( 0., 0., Ecal_, Ecal_ );
  TLorentzVector q4 = p4_nu - p4_lep;
  Q2_ = -1. * q4.Mag2();

  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.101.092001

  // Just the magnitudes
  // pTx_ = pT_ * std::sin( delta_alphaT * M_PI / 180. );
  // pTy_ = pT_ * std::cos( delta_alphaT * M_PI / 180. );

  TVector3 z_unit_vec( 0., 0., 1. );

  TVector2 xT_unit_vec = z_unit_vec.Cross( p3_lep ).XYvector().Unit();
  pTx_ = xT_unit_vec.X()*p2_T.X() + xT_unit_vec.Y()*p2_T.Y();

  TVector2 yT_unit_vec = ( -p3_lep ).XYvector().Unit();
  pTy_ = yT_unit_vec.X()*p2_T.X() + yT_unit_vec.Y()*p2_T.Y();

  // JLab light cone variables

  TLorentzVector p4_miss = p4_lep + p4_had - p4_nu;

  E_miss_ = std::abs( p4_miss.E() );
  p_miss_ = p4_miss.Vect().Mag();

  // Avoid Ecal assumption (suggestion from Jackson)
  p_miss_minus_ = ( E_lep - p3_lep.Z() ) + ( E_had - p3_had.Z() );

  double k_miss_numer = pT_*pT_ + m_had*m_had;
  double k_miss_denom = p_miss_minus_ * ( 2.*m_had - p_miss_minus_ );

  // Expression from Jackson's GlueX note
  double k_miss2 = m_had*m_had * k_miss_numer / k_miss_denom - m_had*m_had;

  k_miss_ = real_sqrt( k_miss2 );

  A_ = p_miss_minus_ / m_had;

  // MINERvA longitudinal and total variables

  // See https://journals.aps.org/prc/pdf/10.1103/PhysRevC.95.065501
  double mA = 22.*NEUTRON_MASS + 18.*PROTON_MASS - 0.34381; // GeV

  // For the calculation of the excitation energies, see Table 7 from
  // https://doi.org/10.1140/epjc/s10052-019-6750-3
  double mA_prime = mA - NEUTRON_MASS + binding_energy; // GeV

  // For the calculation of p_n, see Eq. (8) from
  // https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.121.022504
  double R = mA + p3_lep.Z() + p3_had.Z() - E_lep - E_had;

  // Beyond the transverse variables
  // Based on Andy F's xsec meeting presentation
  // https://tinyurl.com/uboone-beyond-stv

  // GeV, after discussion with Andy F who got the numbers from Jan S
  Ecal_MB_ = E_lep + had_KE + binding_energy;

  TLorentzVector p4_nu_MB( 0., 0., Ecal_MB_, Ecal_MB_ );
  TLorentzVector q4_MB = p4_nu_MB - p4_lep;

  switch (calc_opt) {
  case kOpt1:
  case kOpt2:
  case kOpt3:
    pL_ = 0.5*R - ( mA_prime*mA_prime + pT_*pT_ ) / ( 2.*R );
    break;
  case kOpt4:
    pL_ = p3_lep.Z() + p3_had.Z() - Ecal_MB_;
    break;
  default:
    throw std::runtime_error( "Unrecognized option for KICalculator "
      + std::to_string( calc_opt_ ) );
  }

  pn_ = real_sqrt( pT_*pT_ + pL_*pL_ );

  TVector3 p3_pn( p3_T.X(), p3_T.Y(), pL_ );

  TVector3 q3_MB = q4_MB.Vect();
  TVector3 q3T_MB( q3_MB.X(), q3_MB.Y(), 0. );

  delta_alpha3D_q_ = std::acos( ( q3_MB * p3_pn )
    / ( q3_MB.Mag() * pn_ ) ) * 180. / M_PI;
  if ( delta_alpha3D_q_ > 180. ) delta_alpha3D_q_ -= 180.;
  if ( delta_alpha3D_q_ < 0. ) delta_alpha3D_q_ += 180.;

  delta_alpha3D_mu_ = std::acos( -( p3_lep * p3_pn )
    / ( p3_lep.Mag() * pn_ ) ) * 180. / M_PI;
  if ( delta_alpha3D_mu_ > 180. ) delta_alpha3D_mu_ -= 180.;
  if ( delta_alpha3D_mu_ < 0. ) delta_alpha3D_mu_ += 180.;

  delta_phi3D_ = std::acos( ( q3_MB * p3_had )
    / ( q3_MB.Mag() * p3_had.Mag() ) ) * 180. / M_PI;
  if ( delta_phi3D_ > 180. ) delta_phi3D_ -= 180.;
  if ( delta_phi3D_ < 0. ) delta_phi3D_ += 180.;

  // Magnitudes
  pn_perp_ = pn_ * std::sin( delta_alpha3D_q_ * M_PI / 180. );
  pn_par_ = pn_ * std::cos( delta_alpha3D_q_ * M_PI / 180. );

  pn_perpx_ = ( q3T_MB.Unit().Cross( z_unit_vec ) ).Dot( p3_pn );
  pn_perpy_ = ( q3_MB.Unit().Cross(
    q3T_MB.Unit().Cross( z_unit_vec ) ) ).Dot( p3_pn );

  // Opening angle between the lepton and the hadronic system
  theta_lep_had_ = std::acos( p3_lep.Dot( p3_had ) / p3_lep.Mag()
    / p3_had.Mag() ) * 180. / M_PI;
}
