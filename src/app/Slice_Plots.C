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


// Define custom color map for known covariance matrix labels
const std::map<std::string, int> color_map = {
  {"total", kBlack},
  {"flux", kRed+1},
  {"reint", kGreen+2},
  //{"neutron_reint", kViolet+1},
  {"xsec_total", kBlue+1},
  {"detVar_total", kOrange+2},
  {"POT", kGray+2},
  {"numTargets", kPink+6},
  {"MCstats", kTeal+3},
  {"EXTstats", kMagenta+2},
  {"BNBstats", kAzure+2}
};

std::map<std::string, int> category_colors = {
  { "Signal", kGreen },
  { "Background", kRed },
  { "OOFV", kOrange+3 },
  { "EXT", kBlack }
};

const std::map<std::string, std::vector<int>> grouped_categories = {
  { "Signal", {5,6,8,9} },
  { "Background", {1,2,3,4,7,10,11,12,13,15} },
  { "OOFV", {14} },
  { "EXT", {0} }
};

namespace {

  // Burke Function
  void make_fractional_uncertainty_barchart(const std::map<std::string, TH1*>& frac_uncertainty_hists, const std::map<std::string, TH1*>& detvar_subcomponents) {
    std::vector<std::string> keys;
    for (const auto& p : frac_uncertainty_hists) {
        if (p.first != "total") keys.push_back(p.first);
    }
    keys.push_back("total");

    const int n = keys.size();
    TH1D* h_bar = new TH1D("h_bar", ";;Fractional Uncertainty", n, 0, n);

    int bin = 1;
    for (const auto& key : keys) {
        double frac = frac_uncertainty_hists.at(key)->GetBinContent(1);
        h_bar->SetBinContent(bin, frac);
        h_bar->GetXaxis()->SetBinLabel(bin, key.c_str());
        ++bin;
    }

    h_bar->LabelsOption("v", "X");
    h_bar->SetStats(false);
    h_bar->SetFillColor(kAzure + 1);
    h_bar->SetBarWidth(0.7);
    h_bar->SetBarOffset(0.15);

    // Build detVar breakdown bar chart
    std::vector<std::string> det_keys;
    for (const auto& p : detvar_subcomponents) {
        det_keys.push_back(p.first);
    }
    const int ndet = det_keys.size();
    TH1D* h_det = new TH1D("h_det", ";;DetVar Subcomponent Uncertainty", ndet, 0, ndet);

    bin = 1;
    for (const auto& key : det_keys) {
        double frac = detvar_subcomponents.at(key)->GetBinContent(1);
        h_det->SetBinContent(bin, frac);
        h_det->GetXaxis()->SetBinLabel(bin, key.c_str());
        ++bin;
    }

    h_det->LabelsOption("v", "X");
    h_det->SetStats(false);
    h_det->SetFillColor(kTeal + 1);
    h_det->SetBarWidth(0.7);
    h_det->SetBarOffset(0.15);

    TCanvas* c = new TCanvas("c_frac_unc_bar", "Fractional Uncertainty Bar Chart", 1000, 1200);
    c->Divide(1,2);
    c->cd(1);
    gPad->SetBottomMargin(0.25);
    h_bar->Draw("bar");
    c->cd(2);
    gPad->SetBottomMargin(0.25);
    h_det->Draw("bar");

    c->SaveAs("fractional_uncertainty_barchart.pdf");
  }


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
    bnb_hist->SetLineStyle(2);
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
    //stat_err_hist->SetFillColor( kBlack );
    //stat_err_hist->SetLineColor( kBlack );
    //stat_err_hist->SetLineWidth( 2 );
    //stat_err_hist->SetFillStyle( 3004 );
    stat_err_hist->SetFillColor( kBlack );
    stat_err_hist->SetLineColor( kBlack );
    stat_err_hist->SetLineWidth( 1 );
    stat_err_hist->SetFillStyle( 3345 );
  }

  // Grouped mode
  TLegend* add_legend(
    THStack* stack,
    const std::map<std::string, TH1D*>& grouped_hists,
    const std::map<std::string, int>& category_colors,
    double x1 = 0.47, double y1 = 0.6,
    double x2 = 0.75, double y2 = 0.88)
  {
    TLegend* legend = new TLegend(x1, y1, x2, y2);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);

    double total = 0.0;
    for (const auto& pair : grouped_hists) {
      total += pair.second->Integral();
    }

    for (const auto& pair : grouped_hists) {
      const std::string& group_name = pair.first;
      TH1D* hist = pair.second;
      double frac = (total > 0.) ? hist->Integral() / total : 0.0;
      std::string label = group_name + Form(" (%.1f%%)", frac * 100.0);

      legend->AddEntry(hist, label.c_str(), "f");
    }

    // Optionally add BNB Data marker-style entry (you can conditionally pass it in if needed)
    for (const auto& hobj : *stack->GetHists()) {
      TH1D* h = dynamic_cast<TH1D*>(hobj);
      if (!h) continue;
      if (std::string(h->GetName()).find("bnb") != std::string::npos) {
        legend->AddEntry(h, "BNB Data", "lep");
        break;
      }
    }

    return legend;
  }

  // Stacked Mode
    TLegend* add_legend(
    THStack* stack,
    const std::map<int, std::pair<std::string, int>>& cat_map,
    double x1 = 0.5, double y1 = 0.4,
    double x2 = 0.88, double y2 = 0.88)
  {
    TLegend* legend = new TLegend(x1, y1, x2, y2);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);

    double total = 0.0;
    for (const auto& hobj : *stack->GetHists()) {
      TH1D* h = dynamic_cast<TH1D*>(hobj);
      if (!h) continue;
      total += h->Integral();
    }

    for (const auto& hobj : *stack->GetHists()) {
      TH1D* hist = dynamic_cast<TH1D*>(hobj);
      if (!hist) continue;

      std::string name = hist->GetName();
      if (name.find("bnb") != std::string::npos || name.find("data") != std::string::npos) {
        continue;
      }

      int cat_id = -1;
      if (name.rfind("MC", 0) == 0) {
        try {
          cat_id = std::stoi(name.substr(2));
        } catch (...) {
          continue;
        }
      } else {
        continue;
      }

      std::string label = Form("Category %d", cat_id);
      auto it = cat_map.find(cat_id);
      if (it != cat_map.end()) {
        label = it->second.first;
      }

      double frac = total > 0. ? hist->Integral() / total : 0.0;
      label += Form(" (%.1f%%)", frac * 100.0);
      legend->AddEntry(hist, label.c_str(), "f");
    }

    for (const auto& hobj : *stack->GetHists()) {
      TH1D* h = dynamic_cast<TH1D*>(hobj);
      if (!h) continue;
      std::string name = h->GetName();
      if (name.find("bnb") != std::string::npos || name.find("data") != std::string::npos) {
        legend->AddEntry(h, "BNB Data", "lep");
        break;
      }
    }

    return legend;
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

    // Burke Debugging for unfolder.C
    std::cout << "[BNB Data: Bin Content]" << std::endl;
    for (int bin = 1; bin <= slice_bnb->hist_->GetNbinsX(); ++bin) {
      double content = slice_bnb->hist_->GetBinContent(bin);
      std::cout << "  Bin " << bin << ": Value = " << content << std::endl;
    }
    // End Burke Debugging

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

    //THStack* slice_pred_stack = new THStack( "mc+ext", "" );
    //slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

    const auto& sel_for_cat = syst.get_selection_for_categories();
    const auto& cat_map = sel_for_cat.CategoryMap();

    // Begin Burke Edits
    bool use_grouped_mode = false;
    // Set to `true` to enable minimal plotting mode:
    // - Only draws the MC+EXT prediction and the BNB data histogram
    // - Disables stacked category histograms
    // - Useful for simple data vs. prediction comparisons (e.g., paper drafts)
    const bool minimal_mode = false; 

    if (use_grouped_mode && minimal_mode) {
      throw std::runtime_error("Cannot use both grouped plotting mode and minimal plotting mode at the same time.");
    }

    std::map<std::string, TH1D*> grouped_hists;

    THStack* slice_pred_stack = nullptr;

    if (!minimal_mode) {
      slice_pred_stack = new THStack("mc+ext", "");
      slice_pred_stack->Add(slice_ext->hist_.get()); // extBNB
    }

    if (use_grouped_mode) {
      grouped_hists["EXT"] = dynamic_cast<TH1D*>(slice_ext->hist_.get());
    }

    if (use_grouped_mode) {
      std::map<std::string, TH1D*> temp_grouped;

      int cat_bin_index = cat_map.size();
      for (auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter) {
        int cat = iter->first;
        int color = iter->second.second;

        TH1D* temp_mc_hist = category_hist->ProjectionY("temp_mc_hist", cat_bin_index, cat_bin_index);
        temp_mc_hist->SetDirectory(nullptr);
  
        SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(*temp_mc_hist, slice);
        temp_slice_mc->hist_->SetName(("MC" + std::to_string(cat)).c_str());
        temp_slice_mc->hist_->SetDirectory(nullptr);

        for (const auto& group : grouped_categories) {
          const std::string& group_name = group.first;
          const auto& cat_ids = group.second;
          if (std::find(cat_ids.begin(), cat_ids.end(), cat) != cat_ids.end()) {
            if (temp_grouped.count(group_name)) {
              temp_grouped[group_name]->Add(temp_slice_mc->hist_.get());
            } else {
              TH1D* h_clone = dynamic_cast<TH1D*>(temp_slice_mc->hist_->Clone(("group_" + group_name).c_str()));
              h_clone->SetDirectory(nullptr);
              temp_grouped[group_name] = h_clone;
            }
            break;
          }
        }

        --cat_bin_index;
      }
  
      for (const auto& pair : temp_grouped) {
        const std::string& group_name = pair.first;
        TH1D* hist = pair.second;
	if (group_name == "EXT") {
          // Use the original slice_ext histogram (with hatched styling)
          continue;  // already added to both slice_pred_stack and grouped_hists
	}

        int color = category_colors.at(group_name);
        set_mc_histogram_style(0, hist, color);
        slice_pred_stack->Add(hist);
        grouped_hists[group_name] = hist;
      }
 
    } else if (!minimal_mode) {
      int cat_bin_index = cat_map.size();
      for (auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter) {
        int cat = iter->first;
        int color = iter->second.second;

        TH1D* temp_mc_hist = category_hist->ProjectionY("temp_mc_hist", cat_bin_index, cat_bin_index);
        temp_mc_hist->SetDirectory(nullptr);

        SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(*temp_mc_hist, slice);
        temp_slice_mc->hist_->SetName(("MC" + std::to_string(cat)).c_str());
  
        set_mc_histogram_style(cat, temp_slice_mc->hist_.get(), color);
        slice_pred_stack->Add(temp_slice_mc->hist_.get());

        --cat_bin_index;
      }
    }

    // End Burke Edits


    // Go in reverse so that, if the signal is defined first in the map, it
    // ends up on top. Note that this index is one-based to match the ROOT
    // histograms
    /*int cat_bin_index = cat_map.size();
    for ( auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter )
    {
      int cat = iter->first;
      int color = iter->second.second;
      TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
        cat_bin_index, cat_bin_index );
      temp_mc_hist->SetDirectory( nullptr );

      SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
        *temp_mc_hist, slice  );

      temp_slice_mc->hist_->SetName(("MC" + std::to_string(cat)).c_str());

      set_mc_histogram_style( cat, temp_slice_mc->hist_.get(), color );

      slice_pred_stack->Add( temp_slice_mc->hist_.get() );

      std::string cat_col_prefix = "MC" + std::to_string( cat );

      --cat_bin_index;
    }*/

    TCanvas* c1 = new TCanvas;
    slice_bnb->hist_->SetLineColor( kBlack );
    slice_bnb->hist_->SetLineWidth( 3 );
    slice_bnb->hist_->SetMarkerStyle( kFullCircle );
    slice_bnb->hist_->SetMarkerSize( 1.2 );
    slice_bnb->hist_->SetStats( false );
    slice_bnb->hist_->SetTitle("");
    double ymax = std::max( slice_bnb->hist_->GetMaximum(),
      slice_mc_plus_ext->hist_->GetMaximum() ) * 1.25;
    slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., ymax );

    slice_bnb->hist_->Draw( "e" );

    if (minimal_mode) {
      slice_mc_plus_ext->hist_->SetLineWidth( 3 );
      slice_mc_plus_ext->hist_->SetLineColor(kRed);
      slice_mc_plus_ext->hist_->SetTitle("");
      //set_stat_err_histogram_style(slice_mc_plus_ext->hist_.get());
      //slice_mc_plus_ext->hist_->Draw( "same hist e" );
      slice_mc_plus_ext->hist_->Draw( "same hist" );
      // Clone and draw the MC+EXT stat error band
      TH1D* mc_err_band = dynamic_cast<TH1D*>(slice_mc_plus_ext->hist_->Clone("mc_err_band"));
      mc_err_band->SetDirectory(nullptr);
      mc_err_band->SetTitle("");
      set_stat_err_histogram_style(mc_err_band);
      mc_err_band->Draw("E2 same");  // Hatched fill

      TLegend* lg_min = new TLegend(0.4, 0.75, 0.6, 0.88);
      lg_min->SetBorderSize(0);
      lg_min->SetFillStyle(0);
      lg_min->SetTextFont(42);
      lg_min->SetTextSize(0.03);

      lg_min->AddEntry(slice_mc_plus_ext->hist_.get(), "MC + EXT (CV)", "l");
      lg_min->AddEntry(slice_bnb->hist_.get(), "Nuwro Fake Data", "lep");
      lg_min->Draw("same");
    } else {
      slice_pred_stack->Draw("hist same");

      TH1D* mc_err_band = dynamic_cast<TH1D*>(slice_mc_plus_ext->hist_->Clone("mc_err_band"));
      mc_err_band->SetDirectory(nullptr);
      set_stat_err_histogram_style(mc_err_band);
      mc_err_band->Draw("E2 same");
    }

    //slice_pred_stack->Draw( "hist same" );

    // Begin Burke Edit
    //TH1D* mc_err_band = dynamic_cast<TH1D*>(slice_mc_plus_ext->hist_->Clone("mc_err_band"));
    //set_stat_err_histogram_style(mc_err_band);
    //mc_err_band->Draw("E2 same");  // E2 draws the error band with fill

    if (!minimal_mode) {
      if (use_grouped_mode) {
        TLegend* lg = add_legend(slice_pred_stack, grouped_hists, category_colors);
        lg->AddEntry(slice_bnb->hist_.get(), "Nuwro F.D. (1.15e21 POT)", "lep");
        lg->Draw();
      }
      else {
        TLegend* lg = add_legend(slice_pred_stack, cat_map);
        lg->AddEntry(slice_bnb->hist_.get(), "Nuwro F.D. (1.15e21 POT)", "lep");
        lg->Draw();
      }
    }
    // End Burke Edit



    // Get number of bins from the first histogram in the stack
    /*int num_bins = slice_pred_stack->GetHistogram()->GetNbinsX(); 
    
    for (int bin_idx = 1; bin_idx <= num_bins; ++bin_idx) {  // ROOT bins start at 1
      double total_bin_value = 0.0;
      
      // Sum over all histograms in the stack
      TIter next(slice_pred_stack->GetHists());
      TH1* hist;
      while ((hist = (TH1*)next())) {
        total_bin_value += hist->GetBinContent(bin_idx);
      }
      std::cout << "  Bin " << bin_idx << ": " << total_bin_value << std::endl;
    }
    */
    ///////////////////////////////////

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
    
    //Burke Edit is adding "neutron_reint" to this list
    const std::vector< std::string > cov_mat_keys = { "total",
      "detVar_total", "flux", "reint", "neutron_reint", "xsec_total", "POT", "numTargets",
      "MCstats", "EXTstats", "BNBstats"
    };
    /*const std::vector< std::string > cov_mat_keys = { "total",
      "detVar_total", "flux", "reint", "xsec_total", "POT", "numTargets",
      "MCstats", "EXTstats", "BNBstats"
    };*/

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

      // Begin Burke Edit
      //if ( color <= 9 ) ++color;
      //if ( color == 5 ) ++color;
      //if ( color >= 10 ) color += 10;
      auto col_it = color_map.find(key);
      int line_color = (col_it != color_map.end()) ? col_it->second : kGray+1;
      slice_for_syst->hist_->SetLineColor(line_color);

      // End Burke Edit
      //slice_for_syst->hist_->SetLineColor( color );
      slice_for_syst->hist_->SetLineWidth( 3 );
    }

    TCanvas* c2 = new TCanvas;
    TLegend* lg2 = new TLegend( 0.7, 0.7, 0.9, 0.9 );
    
    // Begin Burke Edit
    auto* total_frac_err_hist = frac_uncertainty_hists.at( "total" );
    auto* bnbstats_hist = frac_uncertainty_hists.count("BNBstats") ? frac_uncertainty_hists.at("BNBstats") : nullptr;
    double total_max = total_frac_err_hist->GetMaximum();
    double bnb_max = bnbstats_hist ? bnbstats_hist->GetMaximum() : 0.0;

    double ymax_burke_edit = std::max(total_max, bnb_max) * 1.05;

    total_frac_err_hist->SetStats( false );
    //total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
    //  total_frac_err_hist->GetMaximum() * 1.05 );
    total_frac_err_hist->GetYaxis()->SetRangeUser(0., ymax_burke_edit);
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineWidth( 3 );
    total_frac_err_hist->Draw( "hist" );
    // End Burke Edit

    lg2->AddEntry( total_frac_err_hist, "total", "l" );

    for ( auto& pair : frac_uncertainty_hists ) {
      const auto& name = pair.first;
      TH1* hist = pair.second;

      // We already plotted the "total" one above
      if ( name == "total" ) continue;

      if ( name == "BNBstats" ) {hist->SetLineStyle(2);}

      lg2->AddEntry( hist, name.c_str(), "l" );
      hist->Draw( "same hist" );

      std::cout << name << " frac err in bin #1 = "
        << hist->GetBinContent( 1 )*100. << "%\n";
    }

    lg2->Draw( "same" );

    // Burke Edits
    
    if (matrix_map.count("neutron_reint")) {
        const auto& neutron_cov = matrix_map.at("neutron_reint");

        SliceHistogram* neutron_slice = SliceHistogram::make_slice_histogram(
            *reco_mc_plus_ext_hist, slice, &neutron_cov );

	// remake the nominal histogram
        SliceHistogram* nominal_slice = SliceHistogram::make_slice_histogram(
            *reco_mc_plus_ext_hist, slice );

	// Slice the reinteraction unisim histogram
        TH1D& reint_hist = *(dynamic_cast<TH1D*>(
            syst.rw_universes_.at("weight_neutron_reint").at(0)->hist_reco_.get()));

	TH1D* ext_clone = dynamic_cast<TH1D*>(reco_ext_hist->Clone("ext_clone_for_reint"));
	ext_clone->SetDirectory(nullptr);

        reint_hist.Add(ext_clone);

        SliceHistogram* reint_slice = SliceHistogram::make_slice_histogram(
            reint_hist, slice );

        std::cout << "[Neutron Reinteraction fractional errors and error per bin:]" << std::endl;

        for (const auto& bin_pair : slice.bin_map_) {
            int bin_idx = bin_pair.first;
	    double variance = std::pow(neutron_slice->hist_->GetBinError(bin_idx), 2);
	    double nominal_val = nominal_slice->hist_->GetBinContent(bin_idx);
	    std::cout << "  Bin " << bin_idx
              << " | Variance: " << variance
	      << " | Frac Unc.: " << std::sqrt(variance) / nominal_val
	      << " | nominal value: " << nominal_val
              << std::endl;
        }
	delete nominal_slice;
        delete neutron_slice;
    }

    std::cout << "Total frac error in bin #1 = "
      << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";

    PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + Form("_%i",FileNameCounter) + Plot_Suffix;
    c2->SaveAs(PlotFileName.c_str());
    FileNameCounter += 1;

    // Begin Burke Addition (detvar ratio plots per slice)
    std::map<std::string, TH1*> detvar_subcomponents;
    std::vector<std::string> detvar_components = {
      "detVarLYatten", "detVarLYdown", "detVarLYrayl", "detVarRecomb2",
      "detVarSCE", "detVarWMAngleXZ", "detVarWMAngleYZ",
      "detVarWMX", "detVarWMYZ"
    };

    TCanvas* c3 = new TCanvas("c3", "DetVar Ratio Breakdown", 800, 600);
    TLegend* lg3 = new TLegend(0.6, 0.55, 0.9, 0.88);
    lg3->SetBorderSize(0);
    lg3->SetFillStyle(0);
    lg3->SetTextFont(42);
    lg3->SetTextSize(0.03);

    std::vector<int> detvar_colors = {
      kRed+1, kBlue+1, kGreen+2, kOrange+2, kViolet+1, kAzure+2, kPink+6, kTeal+3
    };

    int detvar_color_index = 0;
    bool first_detvar_draw = true;

    for (const auto& key : detvar_components) {
      auto it = matrix_map.find(key);
      if (it == matrix_map.end()) continue;

      const auto& cov_matrix = it->second;
      SliceHistogram* detvar_slice = SliceHistogram::make_slice_histogram(
        *reco_mc_plus_ext_hist, slice, &cov_matrix );

      for (const auto& bin_pair : slice.bin_map_) {
        int bin_idx = bin_pair.first;
        double y = detvar_slice->hist_->GetBinContent(bin_idx);
        double err = detvar_slice->hist_->GetBinError(bin_idx);
        double frac = (y > 0.) ? err / y : 0.;
        detvar_slice->hist_->SetBinContent(bin_idx, frac);
        detvar_slice->hist_->SetBinError(bin_idx, 0.);
      }

      int color = detvar_colors[detvar_color_index % detvar_colors.size()];
      ++detvar_color_index;

      detvar_slice->hist_->SetLineColor(color);
      detvar_slice->hist_->SetLineWidth(2);
      detvar_slice->hist_->SetStats(false);

      if (first_detvar_draw) {
        detvar_slice->hist_->SetMinimum(0.);
	detvar_slice->hist_->SetMaximum(.6);
        detvar_slice->hist_->Draw("hist");
        first_detvar_draw = false;
      } else {
        detvar_slice->hist_->Draw("hist same");
      }

      lg3->AddEntry(detvar_slice->hist_.get(), key.c_str(), "l");

      detvar_subcomponents[key] = detvar_slice->hist_.get();
    }

    // Now add the total detVar line (black)
    if (matrix_map.count("detVar_total")) {
      const auto& total_cov = matrix_map.at("detVar_total");
      SliceHistogram* total_detvar_slice = SliceHistogram::make_slice_histogram(
        *reco_mc_plus_ext_hist, slice, &total_cov );

      for (const auto& bin_pair : slice.bin_map_) {
        int bin_idx = bin_pair.first;
        double y = total_detvar_slice->hist_->GetBinContent(bin_idx);
        double err = total_detvar_slice->hist_->GetBinError(bin_idx);
        double frac = (y > 0.) ? err / y : 0.;
        total_detvar_slice->hist_->SetBinContent(bin_idx, frac);
        total_detvar_slice->hist_->SetBinError(bin_idx, 0.);
      }

      total_detvar_slice->hist_->SetLineColor(kBlack);
      total_detvar_slice->hist_->SetLineWidth(3);
      total_detvar_slice->hist_->SetLineStyle(2);
      total_detvar_slice->hist_->Draw("hist same");

      lg3->AddEntry(total_detvar_slice->hist_.get(), "detVar_total", "l");
    }

    lg3->Draw();

    std::string ratio_outname = Plot_OutputDir + "/DetVar_Ratio_Slice_" + std::to_string(sl_idx) + ".png";
    c3->SaveAs(ratio_outname.c_str());
    // End Burke Addition (detvar ratio plots per slice)

    // Burke Bar Chart Edit
    make_fractional_uncertainty_barchart(frac_uncertainty_hists,detvar_subcomponents);
    // End Burke Bar Chart Edit

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
