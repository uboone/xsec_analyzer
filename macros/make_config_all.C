#include "ConfigMakerUtils.hh"
#include "HistUtils.hh"
#include "UniverseMaker.hh"
#include "SliceBinning.hh"

// Placeholder value for the block index for bins in which it is irrelevant
constexpr int DUMMY_BLOCK_INDEX = -1;

// Sideband selection cuts
/*const std::string sel_dirt =
  " !sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && sel_topo_cut_passed && sel_has_muon_candidate"
  " && sel_muon_contained && sel_muon_quality_ok"
  " && sel_muon_passed_mom_cuts && sel_no_reco_showers"
  " && sel_has_p_candidate && sel_protons_contained"
  " && sel_passed_proton_pid_cut && sel_lead_p_passed_mom_cuts";

const std::string sel_NC =
  " sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV && sel_topo_cut_passed"
  " && !sel_has_muon_candidate && muon_candidate_idx != -1"
  " && sel_no_reco_showers && sel_has_p_candidate"
  " && sel_protons_contained && sel_passed_proton_pid_cut"
  " && sel_lead_p_passed_mom_cuts";

const std::string sel_CCNpi =
  " sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && sel_topo_cut_passed && sel_has_muon_candidate"
  " && sel_muon_contained && sel_muon_quality_ok"
  " && sel_muon_passed_mom_cuts && sel_no_reco_showers"
  " && sel_has_p_candidate && sel_protons_contained"
  " && Sum$( pfp_generation_v == 2 && trk_llr_pid_score_v > 0.2 ) > 1"
  " && sel_lead_p_passed_mom_cuts";*/

const std::string sel_combined =
  "((" + sel_dirt + ") || (" + sel_NC + ") || (" + sel_CCNpi + "))";

struct EdgeDef {

  EdgeDef( std::map< double, std::vector<double> >* edges, bool use_overflow,
    const std::string& var_to_use, const std::string& act_var_name,
    const std::string& oth_var_name) : bin_edges_2d_( edges ),
    needs_overflow_( use_overflow ), branch_var_name_( var_to_use ),
    active_var_name_( act_var_name ), other_var_name_( oth_var_name ) {}

  std::map< double, std::vector<double> >* bin_edges_2d_;
  bool needs_overflow_;
  std::string branch_var_name_;
  std::string active_var_name_;
  std::string other_var_name_;
};

// cos_theta_neutron in 1D
std::vector< double > cos_theta_n_1D_edges = { -1., -0.5, 0.0, 0.27,
   0.45, 0.62, 0.76, 0.86, 0.94, 1.0 };

void make_config_all() {

  // Set up an initially empty container to hold the slice definitions. We'll
  // populate it in parallel with defining the bins themselves.
  SliceBinning sb;

  // Set the variables to use when defining phase-space slices
  sb.slice_vars_ = {
    { "cos#theta_{n}", "\\cos\\theta_{\\n}", "", "" },
    { "bin number", "\\text{ bin number}", "", "" }
  };

  std::string selection = "sel_CC1muNnXp0pi";
  std::string signal_def = "mc_is_signal";

  // By construction, MC event categories 5-11 contain all beam-correlated
  // backgrounds. This list is therefore comprehensive apart from cosmic
  // overlay stuff which is directly measured in a dedicated sample.
  std::vector< std::string > background_defs = {
    "category == 5", "category == 6", "category == 7", "category == 8",
    "category == 9", "category == 10", "category == 11"
  };

  std::vector< TrueBin > true_bins;
  std::vector< RecoBin > reco_bins;

  int block_index = -1;

//// 1D cos_theta_n block

  // Increment the block index for the current set of bin edges
  ++block_index;

  // Get the index of the first analysis bin in this block
  first_block_bin_idx = reco_bins.size();

  // Get the index for the cos_theta_mu variable definition. We will use
  // it to make a new slice while also defining the 1D bins.
  int cos_var_idx = find_slice_var_index( "cos#theta_{n}",
    sb.slice_vars_ );

  size_t num_cos_1D_edges = cos_theta_n_1D_edges.size();
  size_t num_cos_1D_bins = 0u;
  if ( num_cos_1D_edges >= 2u ) num_cos_1D_bins
    = num_cos_1D_edges - 1u;

  // Before defining each bin, make a new Slice object and set up the
  // corresponding ROOT histogram within it
  auto& cos_slice = add_slice( sb, cos_theta_mn_1D_edges, cos_var_idx );

  for ( size_t b = 0u; b < num_cos_1D_bins; ++b ) {

    double cos_low = cos_theta_n_1D_edges.at( b );
    double cos_high = cos_theta_n_1D_edges.at( b + 1 );

    std::stringstream true_ss;
    true_ss << signal_def << " && mc_p3_mu.CosTheta() >= " << cos_low
     << " && mc_p3_mu.CosTheta() < " << cos_high;

    std::string true_bin_def = true_ss.str();

    true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_index );

    std::stringstream reco_ss;
    reco_ss << selection << " && p3_mu.CosTheta() >= " << cos_low
     << " && p3_mu.CosTheta() < " << cos_high;

    std::string reco_bin_def = reco_ss.str();

    // Here we use a trick: the current analysis bin index is equal
    // to the size of the reco_bins vector before we add the new element.
    size_t ana_bin_idx = reco_bins.size();
    // Here's another trick: the call to operator[]() below will create
    // a new map entry if needed. We then insert the current analysis
    // bin index into the map entry.
    cos_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

    // Define the new bin and add it to the vector of reco bins
    reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_index );

  } // loop over 1D cos_theta_mu bins

//// Final definitions

  int num_bins = reco_bins.size();

  //// Create a slice showing just the sideband bins as a function of bin number
  //auto& sideband_slice = add_slice( sb, num_ord_bins, num_ord_bins, num_bins,
  //  bin_number_var_idx );
  //for ( int ab = num_ord_bins; ab < num_bins; ++ab ) {
  //  // The ROOT histogram bins are one-based, so we correct for this here
  //  sideband_slice.bin_map_[ ab + 1 ].insert( ab );
  //}

  // Create a slice showing all results together as a function of bin number
  auto& bin_num_slice = add_slice( sb, num_bins, 0, num_bins,
    bin_number_var_idx );
  for ( int ab = 0; ab < num_bins; ++ab ) {
    // The ROOT histogram bins are one-based, so we correct for this here
    bin_num_slice.bin_map_[ ab + 1 ].insert( ab );
  }

  // Add true bins for the background categories of interest. We'll use a
  // dummy block index now since it's only important for signal bins.
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin, DUMMY_BLOCK_INDEX );
  }

  // Dump this information to the output files
  std::ofstream out_file( "myconfig_all.txt" );
  out_file << "ALL\n";
  out_file << "stv_tree\n";
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  // Also write a SliceBinning configuration file
  std::ofstream sb_file( "mybins_all.txt" );
  sb_file << sb;
}
