// Standard library includes
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <set>
#include <string>
#include <sstream>
#include <iostream>

// ROOT includes
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TPad.h"

// STV analysis includes
#include "XSecAnalyzer/Selections/EventCategoriesNp1pi.hh"
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/FilePropertiesManager.hh"
#include "XSecAnalyzer/HistUtils.hh"
#include "XSecAnalyzer/PlotUtils.hh"
#include "XSecAnalyzer/EventCategoryInterpreter.hh" 

// Abbreviation to make using the enum class easier
using NFT = NtupleFileType;
const int FontStyle = 132;
std::vector<std::string> eventCategoryNames = {"EXT ", "kUknown", "kSignal CCQE ", "kSignal MEC ", "kSignal RES ", "kSignal DIS ", "kSignal COH ", "kSignal Other ", "kNumuCCpi0 ", "kNuMuCC0piXp ", "kNuMuCCOther ", "kNuECC ", "kNC ", "kOOFV ", "kOther "};

TH2D* make_tally_histogram( const std::string& hist_name_prefix,
  const std::string& sel, const std::set<int>& runs,
  const std::string& mc_event_weight = DEFAULT_MC_EVENT_WEIGHT )
{
  // Keep the expression used as the selection inside a pair of
  // parentheses. This will ensure that it "plays well" with the
  // other pieces of the TTree::Draw() commands below.
  std::string selection = '(' + sel + ')';

  int num_category_bins = static_cast<int>( EventCategoryNp1pi::kOther ) + 3;

  // The y-axis indicates whether an event was selected, while the x-axis gives
  // the true EventCategory. An x value in the bin with low edge -1 is reserved
  // for EXT events, while -2 is reserved for BNB events.
  std::string tally_hist_name = hist_name_prefix + "_tally_hist";
  TH2D* event_tally_hist = new TH2D( tally_hist_name.c_str(),
    "; event category; selected; events", num_category_bins, -2.,
    num_category_bins - 2, 2, 0., 2 );

  // Get access to the singleton utility classes that we'll need
  const EventCategoryInterpreter& eci = EventCategoryInterpreter::Instance();
  const FilePropertiesManager& fpm = FilePropertiesManager::Instance();

  // Consider samples for data taken with the beam on, data taken with the beam
  // off, and CV MC samples for numus, intrinsic nues, and dirt events
  constexpr std::array< NFT, 5 > file_types = { NFT::kOnBNB, NFT::kExtBNB,
    NFT::kNumuMC, NFT::kIntrinsicNueMC, NFT::kDirtMC };

  // Similar array that includes only the CV MC samples
  constexpr std::array< NFT, 3 > mc_file_types = { NFT::kNumuMC,
    NFT::kIntrinsicNueMC, NFT::kDirtMC };

  // Prepare TChains needed to loop over the event ntuples to be analyzed. Also
  // prepare maps to keep track of the corresponding POT normalizations and
  // total number of triggers (the latter of these is actually used only for
  // data samples).
  std::map< NFT, std::unique_ptr<TChain> > tchain_map;
  std::map< NFT, double > pot_map;
  std::map< NFT, long > trigger_map;
  for ( const auto& type : file_types ) {
    tchain_map.emplace( std::make_pair(type, new TChain("stv_tree")) );
    pot_map[ type ] = 0.;
    trigger_map[ type ] = 0;
  }

  // Add files for each of the selected runs to the appropriate TChain. Also
  // update the corresponding POT normalizations. Use the FilePropertiesManager
  // to find the right ntuple files for each run.
  const auto& ntuple_map = fpm.ntuple_file_map();
  const auto& data_norm_map = fpm.data_norm_map();

  for ( const int& run : runs ) {
    // Get the map storing the ntuple file names for the current run
    const auto& run_map = ntuple_map.at( run );

    for ( const auto& type : file_types ) {

      // Get the set of ntuple files for the current run and sample type
      const auto& ntuple_files = run_map.at( type );

      // Get access to the corresponding TChain, total POT value, and total
      // number of triggers that we want to use to handle these files
      auto* tchain = tchain_map.at( type ).get();
      double& total_pot = pot_map.at( type );
      long& total_triggers = trigger_map.at( type );

      for ( const auto& file_name : ntuple_files ) {
        // Add the current file to the appropriate TChain
        tchain->Add( file_name.c_str() );

        // For data samples, get normalization information from the
        // FilePropertiesManager and add it to the total (it's not stored in
        // the files themselves)
        if ( type == NFT::kOnBNB || type == NFT::kExtBNB ) {
          const auto& norm_info = data_norm_map.at( file_name );
          total_triggers += norm_info.trigger_count_;
          // This will just be zero for beam-off data. We will calculate an
          // effective value using the trigger counts below.
          total_pot += norm_info.pot_;
        }
        // For MC samples, extract the POT normalization from the TParameter
        // stored in the file
        else if ( type == NFT::kNumuMC || type == NFT::kIntrinsicNueMC
          || type == NFT::kDirtMC )
        {
          TFile temp_file( file_name.c_str(), "read" );
          TParameter<float>* temp_pot = nullptr;
          temp_file.GetObject( "summed_pot", temp_pot );
          double pot = temp_pot->GetVal();
          total_pot += pot;
        }
      } // file names
    } // ntuple types
  } // runs

  // We need to scale the beam-off data based on the ratio of the total trigger
  // counts for beam-off and beam-on data. Calculate the ratio here.
  double trigs_on = trigger_map.at( NFT::kOnBNB );
  double trigs_off = trigger_map.at( NFT::kExtBNB );

  double ext_scale_factor = trigs_on / trigs_off;
  std::string ext_scale_factor_str = std::to_string( ext_scale_factor );

  // Fill the beam-off data portion of the histogram using the matching TChain.
  // Use the EXT scale factor as an event weight so that we match the POT
  // normalization of the BNB data.
  TChain* off_chain = tchain_map.at( NFT::kExtBNB ).get();
  off_chain->Draw( (selection + " : -1 >>+" + tally_hist_name).c_str(),
    ext_scale_factor_str.c_str(), "goff" );

  // Fill the beam-on data portion of the histogram using the matching TChain.
  // The BNB events have unit weight.
  TChain* on_chain = tchain_map.at( NFT::kOnBNB ).get();
  on_chain->Draw( (selection + " : -2 >>+" + tally_hist_name).c_str(),
    "1", "goff" );

  // Loop over the different MC samples and collect their contributions. We
  // have to handle them separately in order to get the POT normalization
  // correct.
  for ( const auto& type : mc_file_types ) {

    TChain* mc_ch = tchain_map.at( type ).get();
    double on_pot = pot_map.at( NFT::kOnBNB );
    double mc_pot = pot_map.at( type );

    // Add this sample's contribution to the tally histogram by MC event
    // category
    for ( const auto& pair : eci.label_map() ) {

      EventCategoryNp1pi ec = pair.first;
      std::string ec_str = std::to_string( ec );
      double mc_scale_factor = on_pot / mc_pot;
      std::string mc_scale_factor_str = std::to_string( mc_scale_factor );

      mc_ch->Draw( (selection + " : " + ec_str
        + " >>+" + tally_hist_name).c_str(),
        (mc_event_weight + " * " + mc_scale_factor_str + " * (EventCategory == "
        + ec_str + ')').c_str(), "goff" );

    } // event categories

  } // MC samples

  //// Plot the result and return the finished histogram
  //auto* c1 = new TCanvas;
  //event_tally_hist->Draw( "colz" );
  return event_tally_hist;
}


void calc_effpur( const TH2D& tally_hist, double& efficiency,
  double& purity, double& mc_purity )
{
  // Offset that shifts an EventCategoryNp1pi value to the corresponding x-axis bin
  // index in the tally histogram. This accounts for the extra bins allocated
  // for BNB and EXT data
  const int category_bin_offset = 3;

  const int first_signal_bin = EventCategoryNp1pi::kSignalCCQE + category_bin_offset;
  const int last_signal_bin = EventCategoryNp1pi::kSignalOther + category_bin_offset;

  // Note that we want to skip the bin with the BNB data since it shouldn't
  // factor into our calculation of the total number of selected MC+EXT events
  const int first_nonBNB_category_bin = 2;

  const int first_mc_category_bin = category_bin_offset;
  int last_category_bin = tally_hist.GetNbinsX();
  
  // These reflect how we've defined the y-axis of the tally histogram
  const int missed_bin = 1;
  const int selected_bin = 2;

  double all_signal = tally_hist.Integral( first_signal_bin, last_signal_bin,
    missed_bin, selected_bin );
  
  double selected_signal = tally_hist.Integral( first_signal_bin,
    last_signal_bin, selected_bin, selected_bin );

  double all_selected = tally_hist.Integral( first_nonBNB_category_bin,
    last_category_bin, selected_bin, selected_bin );

  double all_mc_selected = tally_hist.Integral( first_mc_category_bin,
    last_category_bin, selected_bin, selected_bin );

  for (int i = 0; i < 15; ++i) {
	//std::cout << eventCategoryNames[i] << tally_hist.Integral(first_nonBNB_category_bin + i, first_nonBNB_category_bin + (i+1), selected_bin, selected_bin ) << std::endl;
       std::cout << eventCategoryNames[i] << " " << tally_hist.GetBinContent(i+2, 2) << std::endl; 
  }
 

  std::cout << "Total selected (signal): " << selected_signal << std::endl;
  std::cout << "Total selected: " << all_selected << std::endl; 
  efficiency = selected_signal / all_signal;
  purity = selected_signal / all_selected;
  mc_purity = selected_signal / all_mc_selected;
}

std::string runsToString(const std::set<int>& runs) {
    std::ostringstream oss;
    for (auto it = runs.begin(); it != runs.end(); ++it) {
        if (it != runs.begin()) {
            oss << ", ";
        }
        oss << *it;
    }
    return oss.str();
}

void effpur_cutflow() {

  const std::vector< std::string > selection_defs = { "1",
  "reco_vertex_in_FV",
  "reco_vertex_in_FV && pfp_starts_in_PCV", 
  "reco_vertex_in_FV && pfp_starts_in_PCV && topo_cut_passed",
  "nu_mu_cc",
  "nu_mu_cc && no_reco_showers",
  "nu_mu_cc && no_reco_showers && min_3_tracks",
  "nu_mu_cc && no_reco_showers && min_3_tracks && 2_non_proton",
  "nu_mu_cc && no_reco_showers && min_3_tracks && 2_non_proton && all_pfp_contained",
  "nu_mu_cc && no_reco_showers && min_3_tracks && 2_non_proton && all_pfp_contained && all_pfp_in_vtx_proximity",
  "nu_mu_cc && no_reco_showers && min_3_tracks && 2_non_proton && all_pfp_contained && all_pfp_in_vtx_proximity && muon_passed_mom_cuts",
  "nu_mu_cc && no_reco_showers && min_3_tracks && 2_non_proton && all_pfp_contained && all_pfp_in_vtx_proximity && muon_passed_mom_cuts && muon_quality_ok",
  "nu_mu_cc && no_reco_showers && min_3_tracks && 2_non_proton && all_pfp_contained && all_pfp_in_vtx_proximity && muon_passed_mom_cuts && muon_quality_ok && lead_p_passed_mom_cuts",
  "nu_mu_cc && no_reco_showers && min_3_tracks && 2_non_proton && all_pfp_contained && all_pfp_in_vtx_proximity && muon_passed_mom_cuts && muon_quality_ok && lead_p_passed_mom_cuts && pion_passed_mom_cuts"
  };

  // New ntuple flags introduced at each stage of the selection. Used
  // solely for pgfplots output (hence the weird \char 95 to get literal
  // underscores)
  const std::vector< std::string > selection_flags = { "",
    "95reco\\char 95vertex\\char 95in\\char 95FV",
    "95pfp\\char 95starts\\char 95in\\char 95PCV", 
    "95topo\\char 95cut\\char 95passed",
    "95has\\char 95muon\\char 95candidate",
    "95no\\char 95reco\\char 95showers",
    "95min\\char 95\\char 51\\char 95tracks",
    "95\\char 50\\char 95non\\char 95proton",
    "95all\\char 95pfp\\char 95contained",
    "95all\\char 95pfp\\char 95in\\char 95vtx\\char 95proximity",
    "95muon\\char 95passed\\char 95mom\\char 95cuts",
    "95muon\\char 95quality\\char 95ok",
    "95lead\\char 95p\\char 95passed\\char 95mom\\char 95cuts",
    "95pion\\char 95passed\\char 95mom\\char 95cuts",

  };

  const std::vector< std::string > bin_labels = { "no cuts",
  "reco_vtx_in_FV", "pfp_start_in_PCV", "topo_cut_passed",  "nu_mu_cc", "no_reco_showers", "min_3_tracks", "2_non_proton", "all_pfp_contained", "all_pfp_in_vtx_proximity",
  "muon_passed_mom_cuts", "muon_quality_ok",
  "lead_p_passed_mom_cuts", "pion_passed_mom_cuts" };

  const std::vector< std::string > latex_bin_labels = { "no cuts",
  "in FV", "pfps in PCV", "topo score",  "CC Inclusive", "no showers", "min 3 tracks", "2 non proton", "pfps fully contained", "vtx proximity",
  "$\\mu$ momentum limits", "$\\mu$ mom quality ok",
  "$p$ momentum limits", "$\\pi$ momentum limits" };

  const std::set<int>& runs = {1,2,3};

  size_t num_points = selection_defs.size();
  TGraph* eff_graph = new TGraph( num_points );
  TGraph* pur_graph = new TGraph( num_points );
  TGraph* mc_pur_graph = new TGraph( num_points );
  gStyle->SetTitleFontSize(0.);
  eff_graph->SetTitle(" ");
  double eff, pur, mc_pur;

  for ( size_t k = 0u; k < num_points; ++k  ) {

    const auto& selection = selection_defs.at( k );

    std::string k_str = std::to_string( k );
    auto* tally_hist = make_tally_histogram( k_str, selection, runs );
 
    std::cout << "selection = " << selection << '\n';
    calc_effpur( *tally_hist, eff, pur, mc_pur );

    eff_graph->SetPoint( k, k + 1, eff );
    pur_graph->SetPoint( k, k + 1, pur );
    mc_pur_graph->SetPoint( k, k + 1, mc_pur );

    std::cout << "eff = " << eff << '\n';
    std::cout << "pur = " << pur << '\n';
    std::cout << "mc_pur = " << mc_pur << '\n';
    std::cout << "\n\n";

    //TCanvas* c = new TCanvas;
    //tally_hist->Draw( "colz" );
    delete tally_hist;

  } // selection definitions


  TCanvas* c1 = new TCanvas;
  c1->SetBottomMargin(0.3);
  eff_graph->SetLineColor(kBlue);
  eff_graph->SetMarkerColor(kBlue);
  eff_graph->SetLineWidth(3);
  eff_graph->SetMarkerStyle(20);
  eff_graph->GetYaxis()->SetRangeUser( 0., 1. );
  eff_graph->GetYaxis()->SetTitleFont( FontStyle);
  eff_graph->GetXaxis()->SetTitleFont( FontStyle);
  eff_graph->GetYaxis()->SetLabelFont( FontStyle );
  eff_graph->GetXaxis()->SetLabelFont( FontStyle );
  eff_graph->Draw( "alp" );

  for ( int b = 1; b <= bin_labels.size(); ++b ) {
    eff_graph->GetHistogram()->GetXaxis()
      ->SetBinLabel( eff_graph->GetHistogram()->FindBin(b),
      bin_labels.at(b - 1).c_str() );
  }
  eff_graph->Draw( "same" );

  pur_graph->SetLineColor(kRed);
  pur_graph->SetMarkerColor(kRed);
  pur_graph->SetLineWidth(3);
  pur_graph->SetMarkerStyle(20);
  pur_graph->Draw("same lp");

  mc_pur_graph->SetLineColor(kBlack);
  mc_pur_graph->SetMarkerColor(kBlack);
  mc_pur_graph->SetLineWidth(3);
  mc_pur_graph->SetMarkerStyle(20);
  mc_pur_graph->Draw("same lp");

  TLegend* lg = new TLegend(0.1,0.92,0.9,0.99);
  lg->SetNColumns(3);
  lg->SetBorderSize( 0 ); 
  lg->AddEntry( eff_graph, "Efficiency", "lp" );
  lg->AddEntry( pur_graph, "Purity", "lp" );
  lg->AddEntry( mc_pur_graph, "MC purity", "lp" );
  lg->SetTextSize(0.06);
  lg->SetTextFont(FontStyle);
  lg->Draw("same");

  std::string runs_str = runsToString(runs);
  c1->SaveAs( ("/exp/uboone/app/users/kwresilo/xsec_analyzer/plots/eff_pur/cuts_evolution_runs" + runs_str + ".pdf").c_str() );
  // We've finished making the ROOT plot. Dump the results to an input file for
  // plotting offline with pgfplots. The keys of the map are column names,
  // while the values are the corresponding entries in each row.
  std::map< std::string, std::vector<std::string> > pgfplots_table;
  dump_tgraph( "eff", *eff_graph, pgfplots_table );
  dump_tgraph( "pur", *pur_graph, pgfplots_table );
  dump_tgraph( "mcpur", *mc_pur_graph, pgfplots_table );

  // Add the x labels in their own column
  pgfplots_table[ "x_label" ] = std::vector< std::string >();
  for ( const auto& label : latex_bin_labels ) {
    std::string entry = '{' + label + '}';
    pgfplots_table.at( "x_label" ).push_back( entry );
  }

  // Add the selection flags in their own column
  pgfplots_table[ "sel_flag" ] = std::vector< std::string >();
  for ( const auto& flag : selection_flags ) {
    std::string entry = '{' + flag + '}';
    pgfplots_table.at( "sel_flag" ).push_back( entry );
  }

  write_pgfplots_file( ("/exp/uboone/app/users/kwresilo/xsec_analyzer/plots/eff_pur/cutflow_runs" + runs_str + ".txt"), pgfplots_table );
}
