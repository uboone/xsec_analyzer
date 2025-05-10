// Standard library includes
#include <map>
#include <string>

// ROOT includes
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TPad.h"

// STV analysis includes
#include "XSecAnalyzer/FilePropertiesManager.hh"
#include "XSecAnalyzer/HistUtils.hh"
#include "XSecAnalyzer/PlotUtils.hh"
#include "XSecAnalyzer/QuickPlotter.hh"
#include "XSecAnalyzer/Selections/SelectionFactory.hh"

// Abbreviation to make using the enum class easier
using NFT = NtupleFileType;

QuickPlotter::QuickPlotter( const std::set< int >& runs,
  const std::string& sel_name ) : runs_( runs )
{
  SelectionFactory sf;
  sel_for_categories_.reset(
    sf.CreateSelection( sel_name )
  );
}

void QuickPlotter::plot( const std::string& branchexpr,
  const std::string& cutexpr, const std::vector< double >& bin_low_edges,
  const std::vector< std::string >& extra_tree_names,
  const std::string& x_axis_label, const std::string& y_axis_label,
  const std::string& title ) const
{
  // Get the number of bins to use in histograms
  int Nbins = bin_low_edges.size() - 1;

  // Make a counter that iterates each time this function is called. We'll use
  // it to avoid duplicate histogram names (which can confuse ROOT).
  static long plot_counter = -1;
  ++plot_counter;

  // Get access to the singleton FilePropertiesManager
  const FilePropertiesManager& fpm = FilePropertiesManager::Instance();

  // Consider samples for data taken with the beam on, data taken with the beam
  // off, and CV MC samples for numus, intrinsic nues, and dirt events
  constexpr std::array< NFT, 5 > file_types = { NFT::kOnBNB, NFT::kExtBNB,
    NFT::kNumuMC, NFT::kIntrinsicNueMC, NFT::kDirtMC };

  // Similar array that includes only the CV MC samples
  constexpr std::array< NFT, 3 > mc_file_types = { NFT::kNumuMC,
    NFT::kIntrinsicNueMC, NFT::kDirtMC };

  std::string tree_name = "XSecAnalyzer/" + sel_for_categories_->name();

  // Prepare TChains needed to loop over the event ntuples to be analyzed. Also
  // prepare maps to keep track of the corresponding POT normalizations and
  // total number of triggers (the latter of these is actually used only for
  // data samples).
  std::map< NFT, std::vector< std::unique_ptr< TChain > > > tchain_map;
  std::map< NFT, double > pot_map;
  std::map< NFT, long > trigger_map;
  for ( const auto& type : file_types ) {
    std::vector< std::unique_ptr< TChain > > chain_vec;
    chain_vec.push_back( std::make_unique< TChain >( tree_name.c_str() ) );
    for ( const auto& tr_name : extra_tree_names ) {
      chain_vec.push_back( std::make_unique< TChain >( tr_name.c_str() ) );
    }
    tchain_map.emplace( std::make_pair(type, std::move( chain_vec ) ) );
    pot_map[ type ] = 0.;
    trigger_map[ type ] = 0;
  }

  // Add files for each of the selected runs to the appropriate TChain. Also
  // update the corresponding POT normalizations. Use the FilePropertiesManager
  // to find the right ntuple files for each run.
  const auto& ntuple_map = fpm.ntuple_file_map();
  const auto& data_norm_map = fpm.data_norm_map();

  for ( const int& run : runs_ ) {
    // Get the map storing the ntuple file names for the current run
    const auto& run_map = ntuple_map.at( run );

    for ( const auto& type : file_types ) {

      // Get the set of ntuple files for the current run and sample type
      const auto& ntuple_files = run_map.at( type );

      // Get access to the corresponding TChain, total POT value, and total
      // number of triggers that we want to use to handle these files
      auto& tchain_vec = tchain_map.at( type );
      double& total_pot = pot_map.at( type );
      long& total_triggers = trigger_map.at( type );

      for ( const auto& file_name : ntuple_files ) {

        // Get the first TChain to add the later ones as friends
        auto& first_ch = tchain_vec.front();

        // Add the current file to the appropriate TChain
        for ( auto& tchain : tchain_vec ) {
          tchain->Add( file_name.c_str() );
          if ( tchain != first_ch ) first_ch->AddFriend( tchain.get() );
        }

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
          temp_file.GetObject( "XSecAnalyzer/summed_pot", temp_pot );
          double pot = temp_pot->GetVal();
          total_pot += pot;
        }
      } // file names
    } // ntuple types
  } // runs

  // Prepare strings used by multiple histograms below
  std::string hist_name_prefix = "hist_plot" + std::to_string( plot_counter );
  std::string plot_title = title + "; " + x_axis_label + "; " + y_axis_label;

  // Fill the beam-off data histogram using the matching TChain
  std::string off_data_hist_name = hist_name_prefix + "_ext";
  TH1D* off_data_hist = new TH1D( off_data_hist_name.c_str(),
    plot_title.c_str(), Nbins, bin_low_edges.data() );

  auto& off_chain = tchain_map.at( NFT::kExtBNB ).front();
  off_chain->Draw( (branchexpr + " >> " + off_data_hist_name).c_str(),
    cutexpr.c_str(), "goff" );
  //off_data_hist->SetDirectory( nullptr );

  // We need to scale the beam-off data to an effective POT based on the ratio
  // of the total trigger counts for beam-off and beam-on data. Do that here.
  double pot_on = pot_map.at( NFT::kOnBNB );
  double trigs_on = trigger_map.at( NFT::kOnBNB );
  double trigs_off = trigger_map.at( NFT::kExtBNB );

  // Compute the effective POT and store it in the map
  double ext_effective_pot = trigs_off * pot_on / trigs_on;
  pot_map[ NFT::kExtBNB ] = ext_effective_pot;

  // Scale the beam-off data based on the effective POT
  off_data_hist->Scale( pot_on / ext_effective_pot );
  off_data_hist->SetFillColor( 28 );
  off_data_hist->SetLineColor( 28 );
  off_data_hist->SetLineWidth( 2 );
  off_data_hist->SetFillStyle( 3005 );
  off_data_hist->SetStats( false );

  // Fill the beam-on data histogram using the matching TChain
  std::string on_data_hist_name = hist_name_prefix + "_on";
  TH1D* on_data_hist = new TH1D( on_data_hist_name.c_str(),
    plot_title.c_str(), Nbins, bin_low_edges.data() );

  auto& on_chain = tchain_map.at( NFT::kOnBNB ).front();
  on_chain->Draw( (branchexpr + " >> " + on_data_hist_name).c_str(),
    cutexpr.c_str(), "goff" );

  //on_data_hist->SetDirectory( nullptr );
  on_data_hist->Scale( 1. );

  on_data_hist->SetLineColor( kBlack );
  on_data_hist->SetLineWidth( 3 );
  on_data_hist->SetMarkerStyle( kFullCircle );
  on_data_hist->SetMarkerSize( 0.8 );
  on_data_hist->SetStats( false );

  on_data_hist->GetXaxis()->SetTitleOffset( 0.0 );
  on_data_hist->GetXaxis()->SetTitleSize( 0.0 );
  on_data_hist->GetYaxis()->SetTitleSize( 0.05 );
  on_data_hist->GetYaxis()->CenterTitle( true );
  on_data_hist->GetXaxis()->SetLabelSize( 0.0 );

  // This prevents the first y-axis label label (0) to be clipped by the
  // ratio plot
  on_data_hist->SetMinimum( 1e-3 );

  // Initialize empty stacked histograms organized by MC event category (an int
  // code)
  std::map< int, TH1D* > mc_hists;
  // Loop over all MC event categories
  for ( const auto& pair : sel_for_categories_->category_map() ) {
    int cat = pair.first;
    std::string cat_label = pair.second.name_;
    int cat_color_code = pair.second.color_;

    std::string temp_mc_hist_name = hist_name_prefix + "_mc"
      + std::to_string( cat );

    TH1D* temp_mc_hist = new TH1D( temp_mc_hist_name.c_str(),
      plot_title.c_str(), Nbins, bin_low_edges.data() );

    mc_hists[ cat ] = temp_mc_hist;

    //temp_mc_hist->SetDirectory( nullptr );
    temp_mc_hist->SetFillColor( cat_color_code );
    temp_mc_hist->SetLineColor( cat_color_code );
    temp_mc_hist->SetStats( false );
  }

  // Loop over the different MC samples and collect their contributions. We
  // have to handle them separately in order to get the POT normalization
  // correct.

  // Counter that avoids duplicate temporary MC histogram names. This is used
  // to avoid annoying ROOT warnings.
  static int dummy_counter = 0;
  for ( const auto& type : mc_file_types ) {

    auto& mc_ch = tchain_map.at( type ).front();
    double on_pot = pot_map.at( NFT::kOnBNB );
    double mc_pot = pot_map.at( type );

    // Add this sample's contribution to the stacked histograms by MC event
    // category
    std::string category_branch_name = sel_for_categories_->name()
      + "_Category";

    for ( const auto& pair : sel_for_categories_->category_map() ) {

      int cat = pair.first;

      std::string temp_mc_hist_name = hist_name_prefix + "_temp_mc"
        + std::to_string( cat ) + "_number" + std::to_string( dummy_counter );

      ++dummy_counter;

      TH1D* temp_mc_hist = new TH1D( temp_mc_hist_name.c_str(),
        plot_title.c_str(), Nbins, bin_low_edges.data() );

      mc_ch->Draw( (branchexpr + " >> " + temp_mc_hist_name).c_str(),
        (mc_event_weight_ + "*(" + cutexpr + " && " + category_branch_name
        + " == " + std::to_string(cat) + ')').c_str(), "goff" );

      // Scale to the same exposure as the beam-on data
      temp_mc_hist->Scale( on_pot / mc_pot );

      // Add this histogram's contribution (now properly scaled) to the total
      mc_hists.at( cat )->Add( temp_mc_hist );

      // We don't need the temporary histogram anymore, so just get rid of it
      delete temp_mc_hist;

    } // event categories

  } // MC samples

  // All the input histograms are now ready. Prepare the plot.
  auto* c1 = new TCanvas;
  //c1->SetLeftMargin( 0.12 );
  //c1->SetBottomMargin( 1.49 );

  TPad* pad1 = new TPad( "pad1", "", 0.0, 0.23, 1.0, 1.0 );
  pad1->SetBottomMargin( 0 );
  pad1->SetRightMargin( 0.06 );
  pad1->SetLeftMargin( 0.13 );
  pad1->SetGridx();
  pad1->Draw();
  pad1->cd();

  on_data_hist->Draw( "E1" );

  // Stack of categorized MC predictions plus extBNB contribution
  THStack* stacked_hist = new THStack( "mc", "" );

  // Sum all contributions into this TH1D so that we can get the overall
  // statistical uncertainty easily
  TH1D* stat_err_hist = new TH1D(
    ("stat_err_hist_" + hist_name_prefix).c_str(), "",
    Nbins, bin_low_edges.data()
  );

  stacked_hist->Add( off_data_hist );
  stat_err_hist->Add( off_data_hist );

  for ( auto citer = mc_hists.crbegin(); citer != mc_hists.crend();
    ++citer )
  {
    TH1D* hist = citer->second;
    stacked_hist->Add( hist );
    stat_err_hist->Add( hist );
  }

  stacked_hist->Draw( "hist same" );

  on_data_hist->Draw( "E1 same" );

  stat_err_hist->SetFillColor( kBlack );
  stat_err_hist->SetLineColor( kBlack );
  stat_err_hist->SetLineWidth( 2 );
  stat_err_hist->SetFillStyle( 3004 );

  stat_err_hist->Draw( "E2 same" );

  // Adjust y-axis range for stacked plot. Check both the data and the
  // stacked histograms (via the combined stat_err_hist)
  double ymax = stat_err_hist->GetBinContent( stat_err_hist->GetMaximumBin() );
  double ymax2 = on_data_hist->GetBinContent( on_data_hist->GetMaximumBin() );
  if ( ymax < ymax2 ) ymax = ymax2;

  // Redraw the histograms with the updated y-axis range
  on_data_hist->GetYaxis()->SetRangeUser( 0., 1.05*ymax );
  on_data_hist->Draw( "E1 same" );

  // Prepare the plot legend
  TLegend* lg = new TLegend( 0.64, 0.32, 0.94, 0.85 );

  std::string legend_title = get_legend_title( pot_on );
  lg->SetHeader( legend_title.c_str(), "C" );

  lg->AddEntry( on_data_hist, "Data (beam on)", "lp" );
  lg->AddEntry( stat_err_hist, "Statistical uncertainty", "f" );

  double total_events = stat_err_hist->Integral();
  for ( const auto& pair : sel_for_categories_->category_map() ) {
    int cat = pair.first;
    std::string label = pair.second.name_;

    TH1* category_hist = mc_hists.at( cat );

    // Use TH1::Integral() to account for CV reweighting correctly
    double events_in_category = category_hist->Integral();
    double category_percentage = events_in_category * 100. / total_events;

    std::string cat_pct_label = Form( "%.2f%#%", category_percentage );

    lg->AddEntry( category_hist, (label + ", " + cat_pct_label).c_str(), "f" );
  }

  double beam_off_events = off_data_hist->Integral();
  double beam_off_percentage = beam_off_events * 100. / total_events;

  std::string off_pct_label = Form( "%.2f%#%", beam_off_percentage );

  lg->AddEntry( off_data_hist, ("Data (beam off), "
    + off_pct_label).c_str(), "f" );

  lg->SetBorderSize( 0 );

  // Increase the font size for the legend header
  // (see https://root-forum.cern.ch/t/tlegend-headers-font-size/14434)
  TLegendEntry* lg_header = dynamic_cast< TLegendEntry* >(
    lg->GetListOfPrimitives()->First() );
  lg_header->SetTextSize( 0.03 );

  lg->Draw( "same" );

  // Ratio plot
  c1->cd(); // Go back from pad1 to main canvas c1

  TPad* pad2 = new TPad( "pad2", "", 0, 0.01, 1.0, 0.23 );
  pad2->SetTopMargin( 0 );
  pad2->SetFrameFillStyle( 4000 );
  pad2->SetBottomMargin( 0.38 );
  pad2->SetRightMargin( 0.06 );
  pad2->SetLeftMargin( 0.13 );

  pad2->SetGridx();
  pad2->Draw();
  pad2->cd(); // change current pad to pad2

  // Ratio plot
  TH1D* h_ratio = dynamic_cast<TH1D*>( on_data_hist->Clone("h_ratio") );
  h_ratio->SetStats( false );
  h_ratio->Divide( stat_err_hist );
  h_ratio->SetLineWidth( 2 );
  h_ratio->SetLineColor( kBlack );
  h_ratio->SetMarkerStyle( kFullCircle );
  h_ratio->SetMarkerSize( 0.8 );
  h_ratio->SetTitle( "" );

  // x-axis
  h_ratio->GetXaxis()->SetTitle( on_data_hist->GetXaxis()->GetTitle() );
  h_ratio->GetXaxis()->CenterTitle( true );
  h_ratio->GetXaxis()->SetLabelSize( 0.12 );
  h_ratio->GetXaxis()->SetTitleSize( 0.18 );
  h_ratio->GetXaxis()->SetTickLength( 0.05 );
  h_ratio->GetXaxis()->SetTitleOffset( 0.9 );

  // y-axis
  h_ratio->GetYaxis()->SetTitle( "ratio" ); //"#frac{Beam ON}{Beam OFF + MC}" );
  h_ratio->GetYaxis()->CenterTitle( true );
  h_ratio->GetYaxis()->SetLabelSize( 0.08);
  h_ratio->GetYaxis()->SetTitleSize( 0.15 );
  h_ratio->GetYaxis()->SetTitleOffset( 0.35 );

  h_ratio->Draw( "E1" );

  gStyle->SetGridColor( 17 );

  // Adjust y-axis
  double ratio_max = h_ratio->GetBinContent( h_ratio->GetMaximumBin() );
  double ratio_min = h_ratio->GetBinContent( h_ratio->GetMinimumBin() );

  h_ratio->SetMaximum( ratio_max + ratio_max*0.15 );
  h_ratio->SetMinimum( ratio_min - ratio_min*0.2 );

  gPad->Update();

  // Draw a horizontal dashed line at ratio == 1
  TLine* line = new TLine( h_ratio->GetXaxis()->GetXmin(), 1.0,
    h_ratio->GetXaxis()->GetXmax(), 1.0 );
  line->SetLineColor( kBlack );
  line->SetLineStyle( 9 ); // dashed
  line->Draw();

  c1->Update();

  // Ensure the canvas is drawn even when we're working in a ROOT-flavored
  // Jupyter notebook
  c1->Draw();
}

// Overloaded version with constant-width binning
void QuickPlotter::plot( const std::string& branchexpr,
  const std::string& cutexpr, double xmin, double xmax, int Nbins,
  const std::vector< std::string >& extra_tree_names,
  const std::string& x_axis_label, const std::string& y_axis_label,
  const std::string& title ) const
{
  // Generates a vector of bin low edges equivalent to the approach used by the
  // TH1 constructor that takes xmin and xmax in addition to the number of bins
  auto bin_low_edges = get_bin_low_edges( xmin, xmax, Nbins );

  this->plot( branchexpr, cutexpr, bin_low_edges, extra_tree_names,
    x_axis_label, y_axis_label, title );
}
