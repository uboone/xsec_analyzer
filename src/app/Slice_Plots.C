// Standard library includes
#include <algorithm>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/FilePropertiesManager.hh"
#include "XSecAnalyzer/MCC9SystematicsCalculator.hh"
#include "XSecAnalyzer/PlotUtils.hh"
#include "XSecAnalyzer/SliceBinning.hh"
#include "XSecAnalyzer/SliceHistogram.hh"

using NFT = NtupleFileType;

//#define USE_FAKE_DATA ""

namespace {

  void set_mc_histogram_style( int event_category, TH1* mc_hist, int color ) {
    mc_hist->SetFillColor( color );
    mc_hist->SetLineColor( color );
    mc_hist->SetStats( false );
  }

  void set_ext_histogram_style( TH1* ext_hist ) {
    ext_hist->SetFillColor( 28 );
    ext_hist->SetLineColor( 28 );
    ext_hist->SetLineWidth( 2 );
    ext_hist->SetFillStyle( 3005 );
    ext_hist->SetStats( false );
  }

  void set_bnb_data_histogram_style( TH1* bnb_hist ) {

    bnb_hist->SetLineColor( kBlack );
    bnb_hist->SetLineWidth( 3 );
    bnb_hist->SetMarkerStyle( kFullCircle );
    bnb_hist->SetMarkerSize( 0.8 );
    bnb_hist->SetStats( false );

    bnb_hist->GetXaxis()->SetTitleOffset( 0.0 );
    bnb_hist->GetXaxis()->SetTitleSize( 0.0 );
    bnb_hist->GetYaxis()->SetTitleSize( 0.05 );
    bnb_hist->GetYaxis()->CenterTitle( true );
    bnb_hist->GetXaxis()->SetLabelSize( 0.0 );

    // This prevents the first y-axis label label (0) to be clipped by the
    // ratio plot
    bnb_hist->SetMinimum( 1e-3 );
  }

  void set_stat_err_histogram_style( TH1* stat_err_hist ) {
    stat_err_hist->SetFillColor( kBlack );
    stat_err_hist->SetLineColor( kBlack );
    stat_err_hist->SetLineWidth( 2 );
    stat_err_hist->SetFillStyle( 3004 );
  }

} // anonymous namespace

void tutorial_slice_plots(std::string FPM_Config, std::string SYST_Config, std::string SLICE_Config, std::string Univ_Output, std::string Plot_OutputDir) {

  // Counter to ensure plots aren't overwritten
  uint FileNameCounter = 0;
  std::string Plot_Prefix = "SlicePlots";
  std::string Plot_Suffix = ".pdf";
  std::string PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + Form("_%i",FileNameCounter) + Plot_Suffix;

  std::cout << "\nRunning Slice_Plots with options:" << std::endl;
  std::cout << "\tFPM_Config: " << FPM_Config << std::endl;
  std::cout << "\tSYST_Config: " << SYST_Config << std::endl;
  std::cout << "\tSLICE_Config: " <<  SLICE_Config << std::endl;
  std::cout << "\tUniv_Output: " << Univ_Output << std::endl;
  std::cout << "\tPlot_OutputDir: " << Plot_OutputDir << std::endl;
  std::cout << "\t\tWith filename: " << PlotFileName << std::endl;
  std::cout << "\n" << std::endl;

#ifdef USE_FAKE_DATA
  // Initialize the FilePropertiesManager and tell it to treat the NuWro
  // MC ntuples as if they were data
  auto& fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties( FPM_Config );
#endif

  // Check that we can read the universe output file
  TFile* temp_file = new TFile(Univ_Output.c_str(), "read");
  if (!temp_file || temp_file->IsZombie()) {
    std::cerr << "Could not read file: " << Univ_Output << std::endl;
    throw;
  }
  delete temp_file;

  auto* syst_ptr = new MCC9SystematicsCalculator(Univ_Output, SYST_Config);
  auto& syst = *syst_ptr;

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the
  // slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

  #ifdef USE_FAKE_DATA
    // Add the EXT to the "data" when working with fake data
    reco_bnb_hist->Add( reco_ext_hist );
  #endif

  TH2D* category_hist = syst.cv_universe().hist_categ_.get();

  // Total MC+EXT prediction in reco bin space. Start by getting EXT.
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  auto* sb_ptr = new SliceBinning( SLICE_Config );
  auto& sb = *sb_ptr;

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {

    const auto& slice = sb.slices_.at( sl_idx );

    // We now have all of the reco bin space histograms that we need as input.
    // Use them to make new histograms in slice space.
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
      *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

    SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
      *reco_ext_hist, slice, &matrix_map.at("EXTstats") );

    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

    //auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
    //std::cout << "Slice " << sl_idx << ": \u03C7\u00b2 = "
    //  << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bins,"
    //  << " p-value = " << chi2_result.p_value_ << '\n';

    // Build a stack of categorized central-value MC predictions plus the
    // extBNB contribution in slice space
    set_ext_histogram_style( slice_ext->hist_.get() );

    THStack* slice_pred_stack = new THStack( "mc+ext", "" );
    slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

    const auto& sel_for_cat = syst.get_selection_for_categories();
    const auto& cat_map = sel_for_cat.category_map();

    // Go in reverse so that, if the signal is defined first in the map, it
    // ends up on top. Note that this index is one-based to match the ROOT
    // histograms
    int cat_bin_index = cat_map.size();
    for ( auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter )
    {
      int cat = iter->first;
      int color = iter->second.second;
      TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
        cat_bin_index, cat_bin_index );
      temp_mc_hist->SetDirectory( nullptr );

      SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
        *temp_mc_hist, slice  );

      set_mc_histogram_style( cat, temp_slice_mc->hist_.get(), color );

      slice_pred_stack->Add( temp_slice_mc->hist_.get() );

      std::string cat_col_prefix = "MC" + std::to_string( cat );

      --cat_bin_index;
    }

    TCanvas* c1 = new TCanvas;
    slice_bnb->hist_->SetLineColor( kBlack );
    slice_bnb->hist_->SetLineWidth( 3 );
    slice_bnb->hist_->SetMarkerStyle( kFullCircle );
    slice_bnb->hist_->SetMarkerSize( 1.2 );
    slice_bnb->hist_->SetStats( false );
    double ymax = std::max( slice_bnb->hist_->GetMaximum(),
      slice_mc_plus_ext->hist_->GetMaximum() ) * 1.07;
    slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., ymax );

    slice_bnb->hist_->Draw( "e" );

    slice_pred_stack->Draw( "hist same" );

    slice_mc_plus_ext->hist_->SetLineWidth( 3 );
    slice_mc_plus_ext->hist_->SetLineColor(kRed);
    slice_mc_plus_ext->hist_->Draw( "same hist e" );

    slice_bnb->hist_->Draw( "same e" );

    PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + Form("_%i",FileNameCounter) + Plot_Suffix;
    c1->SaveAs(PlotFileName.c_str());
    FileNameCounter += 1;

    // Get the binning and axis labels for the current slice by cloning the
    // (empty) histogram owned by the Slice object
    TH1* slice_hist = dynamic_cast< TH1* >(
      slice.hist_->Clone("slice_hist") );

    slice_hist->SetDirectory( nullptr );

    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map< std::string, TH1* >();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    // Show fractional uncertainties computed using these covariance matrices
    // in the ROOT plot. All configured fractional uncertainties will be
    // included in the output pgfplots file regardless of whether they appear
    // in this vector.
    const std::vector< std::string > cov_mat_keys = { "total",
      "detVar_total", "flux", "reint", "xsec_total", "POT", "numTargets",
      "MCstats", "EXTstats", "BNBstats"
    };

    // Loop over the various systematic uncertainties
    int color = 0;
    for ( const auto& pair : matrix_map ) {

      const auto& key = pair.first;
      const auto& cov_matrix = pair.second;

      SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
        *reco_mc_plus_ext_hist, slice, &cov_matrix );

      // The SliceHistogram object already set the bin errors appropriately
      // based on the slice covariance matrix. Just change the bin contents
      // for the current histogram to be fractional uncertainties. Also set
      // the "uncertainties on the uncertainties" to zero.
      // TODO: revisit this last bit, possibly assign bin errors here
      for ( const auto& bin_pair : slice.bin_map_ ) {
        int global_bin_idx = bin_pair.first;
        double y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
        double err = slice_for_syst->hist_->GetBinError( global_bin_idx );
        double frac = 0.;
        if ( y > 0. ) frac = err / y;
        slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
        slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
      }

      // Check whether the current covariance matrix name is present in
      // the vector defined above this loop. If it isn't, don't bother to
      // plot it, and just move on to the next one.
      auto cbegin = cov_mat_keys.cbegin();
      auto cend = cov_mat_keys.cend();
      auto iter = std::find( cbegin, cend, key );
      if ( iter == cend ) continue;

      frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

      if ( color <= 9 ) ++color;
      if ( color == 5 ) ++color;
      if ( color >= 10 ) color += 10;

      slice_for_syst->hist_->SetLineColor( color );
      slice_for_syst->hist_->SetLineWidth( 3 );
    }

    TCanvas* c2 = new TCanvas;
    TLegend* lg2 = new TLegend( 0.7, 0.7, 0.9, 0.9 );

    auto* total_frac_err_hist = frac_uncertainty_hists.at( "total" );
    total_frac_err_hist->SetStats( false );
    total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
      total_frac_err_hist->GetMaximum() * 1.05 );
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineWidth( 3 );
    total_frac_err_hist->Draw( "hist" );

    lg2->AddEntry( total_frac_err_hist, "total", "l" );

    for ( auto& pair : frac_uncertainty_hists ) {
      const auto& name = pair.first;
      TH1* hist = pair.second;
      // We already plotted the "total" one above
      if ( name == "total" ) continue;

      lg2->AddEntry( hist, name.c_str(), "l" );
      hist->Draw( "same hist" );

      std::cout << name << " frac err in bin #1 = "
        << hist->GetBinContent( 1 )*100. << "%\n";
    }

    lg2->Draw( "same" );

    std::cout << "Total frac error in bin #1 = "
      << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";

    PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + Form("_%i",FileNameCounter) + Plot_Suffix;
    c2->SaveAs(PlotFileName.c_str());
    FileNameCounter += 1;

  } // slices

}

int main(int argc, char* argv[]) {
  if ( argc != 6 ) {
    std::cout << "Usage: Slice_Plots FPM_CONFIG"
	      << " SYST_Config SLICE_Config Univ_Output Plot_OutputDir\n";
    return 1;
  }

  std::string list_file_name( argv[1] );
  std::string univmake_config_file_name( argv[2] );
  std::string output_file_name( argv[3] );

  std::string FPM_Config( argv[1] );
  std::string SYST_Config( argv[2] );
  std::string SLICE_Config( argv[3] );
  std::string Univ_Output( argv[4] );
  std::string Plot_OutputDir( argv[5] );

  tutorial_slice_plots(FPM_Config, SYST_Config, SLICE_Config, Univ_Output, Plot_OutputDir);
  return 0;
}
