// Standard library includes
#include <iomanip>
#include <iostream>
#include <sstream>

#include <algorithm>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/CrossSectionExtractor.hh"
#include "XSecAnalyzer/PGFPlotsDumpUtils.hh"
#include "XSecAnalyzer/SliceBinning.hh"
#include "XSecAnalyzer/SliceHistogram.hh"


struct GeneratorInfo {
  std::string path;
  int lineColor;
  int lineStyle;
  int lineWidth;
  std::string name;
  float scaling;
};

// Helper function that dumps a lot of the results to simple text files.
// The events_to_xsec_factor is a constant that converts expected true event
// counts to a total cross section (10^{-39} cm^2 / Ar) via multiplication.
void dump_overall_results( const UnfoldedMeasurement& result,
  const std::map< std::string, std::unique_ptr<TMatrixD> >& unf_cov_matrix_map,
  double events_to_xsec_factor,
  const std::map< std::string, std::unique_ptr< PredictedTrueEvents > >& pred_map )
{
  // Dump the unfolded flux-averaged total cross sections (by converting
  // the units on the unfolded signal event counts)
  TMatrixD unf_signal = *result.unfolded_signal_;
  unf_signal *= events_to_xsec_factor;
  dump_text_column_vector( "../unfold_output/vec_table_unfolded_signal.txt", unf_signal );

  // Dump similar tables for each of the theoretical predictions (and the fake
  // data truth if applicable). Note that this function expects that the
  // additional smearing matrix A_C has not been applied to these predictions.
  for ( const auto& gen_pair : pred_map ) {
    std::string gen_short_name = gen_pair.second->name();
    TMatrixD temp_gen = gen_pair.second->get_prediction();
    temp_gen *= events_to_xsec_factor;
    dump_text_column_vector( "../unfold_output/vec_table_" + gen_short_name + ".txt",
      temp_gen );
  }

  // No unit conversions are necessary for the unfolding, error propagation,
  // and additional smearing matrices since they are dimensionless
  dump_text_matrix( "../unfold_output/mat_table_unfolding.txt", *result.unfolding_matrix_ );
  dump_text_matrix( "../unfold_output/mat_table_err_prop.txt", *result.err_prop_matrix_ );
  dump_text_matrix( "../unfold_output/mat_table_add_smear.txt", *result.add_smear_matrix_ );

  // Convert units on the covariance matrices one-by-one and dump them
  for ( const auto& cov_pair : unf_cov_matrix_map ) {
    const auto& name = cov_pair.first;
    TMatrixD temp_cov_matrix = *cov_pair.second;
    // Note that we need to square the unit conversion factor for the
    // covariance matrix elements
    temp_cov_matrix *= std::pow( events_to_xsec_factor, 2 );

    std::cout << name << ": " << temp_cov_matrix[0][0] << std::endl;

    dump_text_matrix( "../unfold_output/mat_table_cov_" + name + ".txt", temp_cov_matrix );
  }
}

std::string toLatexScientific(double value) {
  std::stringstream stream;
  stream << std::scientific << std::setprecision(2) << value;
  std::string str = stream.str();
  size_t pos = str.find('e');
  if (pos != std::string::npos) {
      str.replace(pos, 2, " #times 10^{");
      str += "}";
  }
  return str;
}

TH2D TMatrixDToTH2D(const TMatrixD & mat, const char* name, const char* title, double xlow, double xup, double ylow, double yup) {
  int nX = mat.GetNrows();
  int nY = mat.GetNcols();
  TH2D hist(name, title, nX, xlow, xup, nY, ylow, yup);
  for (int i = 0; i < nX; ++i) {
      for (int j = 0; j < nY; ++j) {
          hist.SetBinContent(i+1, j+1, mat[i][j]);
      }
  }
  return hist;
}

TH1D* get_generator_hist(const TString& filePath, const unsigned int sl_idx, const float scaling = 1.f )
{
    // Open the file
    TFile* file = new TFile(filePath, "readonly");

    // Check if the file was successfully opened
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return nullptr;
    }

    // Define the plot names
    std::vector<TString> plotNames = {
        "TrueElectronEnergyPlot",
        "TrueElectronCosBetaPlot",
        "TruePionCosBetaPlot",
        "TrueElectronPionOpeningAnglePlot",
        "TrueTotalPlot"
    };

    // special case for combined histogram, all bins
    if (sl_idx == 5) {

      // retrive all histograms
      std::vector<TH1D*> hists_set;

      for (int idx = 0; idx < plotNames.size(); idx++) {
        TH1D* hist = (TH1D*)file->Get(plotNames[idx]);
        hists_set.push_back(hist);
      }

      // construct new histogram from these
      TH1D* hist = new TH1D("", "", 21, 0, 21);

      for (int h_idx = 0; h_idx < hists_set.size(); h_idx++) {

          int n_bins = hists_set[h_idx]->GetNbinsX();

          for (int bin_idx = 1; bin_idx < n_bins + 1; bin_idx++) {  // root counts from 1

            double bin_content = hists_set[h_idx]->GetBinContent(bin_idx);
            double bin_error = hists_set[h_idx]->GetBinError(bin_idx);

            // not general, just testing
            int combined_bin_idx = h_idx*5 + bin_idx;

            // set bin
            hist->SetBinContent(combined_bin_idx, bin_content);
            hist->SetBinError(combined_bin_idx, bin_error);

          }
      }

      return hist;
    }
    else {

      // Check if the index is valid
      if (sl_idx >= plotNames.size()) {
          std::cerr << "Invalid slice index: " << sl_idx << std::endl;
          return nullptr;
      }

      // Get the histogram from the file
      TH1D* hist = (TH1D*)file->Get(plotNames[sl_idx]);

      // Check if the histogram was successfully retrieved
      if (!hist) {
          std::cerr << "Failed to retrieve histogram: " << plotNames[sl_idx] << std::endl;
          return nullptr;
      }

      // Scale the histogram
      hist->Scale(scaling);

      return hist;
    }

    // fail safe
    std::cerr << "Failed to retrieve generator prediction, sl_idx = " << sl_idx << std::endl;
    return nullptr;
}

void multiply_1d_hist_by_matrix(TMatrixD *mat, TH1 *hist)
{
    // Copy the histogram contents into a column vector
    int num_bins = mat->GetNcols();
    TMatrixD hist_mat(num_bins, 1);
    for (int r = 0; r < num_bins; ++r)
    {
        hist_mat(r, 0) = hist->GetBinContent(r + 1);
    }

    // Multiply the column vector by the input matrix
    // TODO: add error handling here related to matrix dimensions
    TMatrixD hist_mat_transformed(*mat, TMatrixD::EMatrixCreatorsOp2::kMult,
                                  hist_mat);

    // Update the input histogram contents with the new values
    for (int r = 0; r < num_bins; ++r)
    {
        double val = hist_mat_transformed(r, 0);
        hist->SetBinContent(r + 1, val);
    }
}

void UnfolderNuMI(std::string XSEC_Config, std::string SLICE_Config, std::string OutputDirectory, std::string OutputFileName) {

  // set to using fake data
  bool using_fake_data = false;
  bool total_only = false;

  std::cout << "\nRunning Unfolder.C with options:" << std::endl;
  std::cout << "\tXSEC_Config: " << XSEC_Config << std::endl;
  std::cout << "\tSLICE_Config: " << SLICE_Config << std::endl;
  std::cout << "\tOutputDirectory: " << OutputDirectory << std::endl;
  std::cout << "\tOutputFileName: " << OutputFileName << std::endl;
  std::cout << "\n" << std::endl;

  // Use a CrossSectionExtractor object to handle the systematics and unfolding
  auto extr = std::make_unique< CrossSectionExtractor >( XSEC_Config );

  // get unfolded results
  auto* sb_ptr = new SliceBinning( SLICE_Config );
  auto& sb = *sb_ptr;

  auto xsec = extr->get_unfolded_events();
  double conv_factor = extr->conversion_factor();
  const auto& pred_map = extr->get_prediction_map();
  double total_pot = extr->get_data_pot();

  double A_C_total = 1;
  if (total_only) {
    const TMatrixD &A_C_temp =  *xsec.result_.add_smear_matrix_;
    A_C_total = A_C_temp(0,0);
    std::cout << "Total Only Mode: A_C matrix element = " << A_C_total << std::endl;
  }

  // Truth Generator Plots
  const std::string genPath = "/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/generatorFiles/";

  // Format: path, lineColor, lineStyle, lineWidth, name, scaling
  std::vector<GeneratorInfo> generators = {
    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //{genPath + "FlatTreeAnalyzerOutput_NuWro_Combined.root", kTeal+3, 1, 3, "NuWro 21.09.2", 1.f},
    //{genPath + "FlatTreeAnalyzerOutput_NEUT_Combined.root", kOrange+8, 1, 3, "NEUT 5.4.0.1", 1.f},
    //{genPath + "FlatTreeAnalyzerOutput_GiBUU_Combined.root", kMagenta-3, 1, 3, "GiBUU 2023", 1.f},
    //{genPath + "FlatTreeAnalyzerOutput_GiBUU2025_Combined.root", kMagenta-3, 1, 3, "GiBUU 2025", 1.f},
    //{genPath + "FlatTreeAnalyzerOutput_Genie_Combined.root", kRed+1, 1, 3, "GENIE 3.4.2 AR23", 1.f},
  };

  for (int sl_idx = 0; sl_idx < 2; sl_idx++) {

    const auto& slice = sb.slices_.at( sl_idx ); // only considering single slice

    // Make a histogram showing the unfolded true event counts in the current slice
    SliceHistogram* slice_unf = SliceHistogram::make_slice_histogram(
      (*xsec.result_.unfolded_signal_), slice, xsec.result_.cov_matrix_.get() );

    // Temporary copies of the unfolded true event count slices with
    // different covariance matrices
    std::map< std::string, std::unique_ptr<SliceHistogram> > sh_cov_map;
    for ( const auto& uc_pair : xsec.unfolded_cov_matrix_map_ ) {
      const auto& uc_name = uc_pair.first;
      const auto& uc_matrix = uc_pair.second;

      auto& uc_ptr = sh_cov_map[ uc_name ];
      uc_ptr.reset(
        SliceHistogram::make_slice_histogram( (*xsec.result_.unfolded_signal_),
          slice, uc_matrix.get() )
      );
    }

    // Also use the GENIE CV model to do the same
    auto genie_cv_it = pred_map.find("MicroBooNETune");
    TMatrixD genie_cv_truth = genie_cv_it->second->get_prediction();

    for ( const auto& gen_pair : pred_map ) {
      std::cout << "Key: " << gen_pair.first << std::endl;
    }

    SliceHistogram* slice_cv = SliceHistogram::make_slice_histogram(
      genie_cv_truth, slice, nullptr );

    // If present, also use the truth information from the fake data to do the same
    SliceHistogram* slice_truth = nullptr;
    if ( using_fake_data ) {
      auto fake_data_truth_it = pred_map.find("Fakedata");
      TMatrixD fake_data_truth = fake_data_truth_it->second->get_prediction();
      TMatrixD fake_data_truth_cov(fake_data_truth.GetNrows(), fake_data_truth.GetNrows());      

      slice_truth = SliceHistogram::make_slice_histogram( fake_data_truth,
        slice, &fake_data_truth_cov );
    }

    // Keys are legend labels, values are SliceHistogram objects containing
    // true-space predictions from the corresponding generator models
    auto* slice_gen_map_ptr = new std::map< std::string, SliceHistogram* >();
    auto& slice_gen_map = *slice_gen_map_ptr;

    slice_gen_map[ "unfolded data" ] = slice_unf;
    if ( using_fake_data ) {
      slice_gen_map[ "truth" ] = slice_truth;
    }
    slice_gen_map[ "MicroBooNE Tune" ] = slice_cv;

    int var_count = 0;
    std::string diff_xsec_denom;
    std::string name_latex;
    std::string diff_xsec_units_denom;
    std::string diff_xsec_denom_latex;
    std::string diff_xsec_units_denom_latex;
    double other_var_width = 1.;
    for ( const auto& ov_spec : slice.other_vars_ ) {
      double high = ov_spec.high_bin_edge_;
      double low = ov_spec.low_bin_edge_;
      const auto& var_spec = sb.slice_vars_.at( ov_spec.var_index_ );
      if ( high != low && std::abs(high - low) < BIG_DOUBLE ) {
        ++var_count;
        other_var_width *= ( high - low );
        diff_xsec_denom += 'd' + var_spec.name_;
        diff_xsec_denom_latex += " d" + var_spec.latex_name_;
        const std::string& temp_units = var_spec.units_;
        if ( !temp_units.empty() ) {
          diff_xsec_units_denom += " / " + temp_units;
          diff_xsec_units_denom_latex += " / " + var_spec.latex_units_;
        }
      }
    }

    for ( size_t av_idx : slice.active_var_indices_ ) {
      const auto& var_spec = sb.slice_vars_.at( av_idx );
      const std::string& temp_name = var_spec.name_;
      if ( temp_name != "true bin number" ) {
        var_count += slice.active_var_indices_.size();
        name_latex += var_spec.name_;
        diff_xsec_denom += 'd' + var_spec.name_;
        diff_xsec_denom_latex += " d" + var_spec.latex_name_;

        if ( !var_spec.units_.empty() ) {
          diff_xsec_units_denom += " / " + var_spec.units_;
          diff_xsec_units_denom_latex += " / " + var_spec.latex_units_;
        }
      }
    }

    // NOTE: This currently assumes that each slice is a 1D histogram
    int num_slice_bins = slice_unf->hist_->GetNbinsX();
    TMatrixD trans_mat( num_slice_bins, num_slice_bins );
    for ( int b = 0; b < num_slice_bins; ++b ) {
      double width = slice_unf->hist_->GetBinWidth( b + 1 );
      width *= other_var_width;
      trans_mat( b, b ) = ( 1 / (conv_factor*width));
    }

    if (total_only) {
      const TMatrixD &A_C_temp =  *xsec.result_.add_smear_matrix_;
      double A_C_total = A_C_temp(0,0);
      trans_mat(0,0) = trans_mat(0,0) * 1/A_C_total;
      std::cout << "Ac total element: " << A_C_total << std::endl;
    }

    std::string slice_y_title;
    std::string slice_y_latex_title;
    if ( var_count > 0 && sl_idx != 4 && !total_only) {
      slice_y_title += "d";
      slice_y_latex_title += "{$d";
      if ( var_count > 1 ) {
        slice_y_title += "^{" + std::to_string( var_count ) + "}";
        slice_y_latex_title += "^{" + std::to_string( var_count ) + "}";
      }
      slice_y_title += "#sigma/" + diff_xsec_denom;
      slice_y_latex_title += "\\sigma / " + diff_xsec_denom_latex;
    }
    else {
      slice_y_title += "#sigma";
      slice_y_latex_title += "\\sigma";
    }
    slice_y_title += " [10^{-39} cm^{2}" + diff_xsec_units_denom + " / nucleon]";
    slice_y_latex_title += "\\text{ }(10^{-39}\\text{ cm}^{2}"
      + diff_xsec_units_denom_latex + " / \\mathrm{nucleon})$}";

    // Convert all slice histograms from true event counts to differential
    // cross-section units
    for ( auto& pair : slice_gen_map ) {
      auto* slice_h = pair.second;
      slice_h->transform( trans_mat );
      slice_h->hist_->GetYaxis()->SetTitle( slice_y_title.c_str() );
      slice_h->hist_->GetYaxis()->SetTitleSize( 0.055 );
      slice_h->hist_->GetYaxis()->SetTitleOffset( 0.75 );

      slice_h->hist_->GetXaxis()->SetTitleSize( 0.055 );
      slice_h->hist_->GetXaxis()->SetTitleOffset( 1 );

      slice_h->hist_->GetXaxis()->SetLabelSize( 0.05 );
      slice_h->hist_->GetYaxis()->SetLabelSize( 0.05 );
    }

    // Also transform all of the unfolded data slice histograms which have
    // specific covariance matrices
    for ( auto& sh_cov_pair : sh_cov_map ) {
      auto& slice_h = sh_cov_pair.second;
      slice_h->transform( trans_mat );
    }

    // Keys are generator legend labels, values are the results of a chi^2
    // test compared to the unfolded data (or, in the case of the unfolded
    // data, to the fake data truth)
    std::map< std::string, SliceHistogram::Chi2Result > chi2_map;
    std::cout << '\n';
    for ( const auto& pair : slice_gen_map ) {
      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      // Decide what other slice histogram should be compared to this one,
      // then calculate chi^2
      SliceHistogram* other = nullptr;
      // We don't need to compare the unfolded data to itself, so just skip to
      // the next SliceHistogram and leave a dummy Chi2Result object in the map
      if ( name == "unfolded data" ) {
        chi2_map[ name ] = SliceHistogram::Chi2Result();
        continue;
      }
      // Compare all other distributions to the unfolded data
      else {
        other = slice_gen_map.at( "unfolded data" );
      }

      // Store the chi^2 results in the map
      const auto& chi2_result = chi2_map[ name ] = slice_h->get_chi2( *other );

      std::cout << name << ": \u03C7\u00b2 = "
        << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bin";
      if ( chi2_result.num_bins_ > 1 ) std::cout << 's';
      if ( chi2_result.num_bins_ > 1 ) std::cout << ", p-value = " << chi2_result.p_value_ << '\n';
    }

    TCanvas *c1 = new TCanvas(std::to_string(sl_idx).c_str(), std::to_string(sl_idx).c_str(), 2100, 900);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetCanvasPreferGL(0);

    bool noRatioPlot = true;
    const auto rightMargin = 0.33;

    TPad* pad1 = new TPad(("pad1 slice "+std::to_string(sl_idx)).c_str(), "", 0.0, noRatioPlot ? 0.05 : 0.3, 1.0, 1.0);
    pad1->SetBottomMargin(0.125);
    pad1->SetTopMargin(0.1);
    pad1->SetRightMargin(rightMargin);
    pad1->Draw();
    pad1->cd();

    slice_unf->hist_->SetLineColor( kBlack );
    slice_unf->hist_->SetLineWidth( 3 );
    slice_unf->hist_->SetMarkerStyle( kFullCircle );
    slice_unf->hist_->SetMarkerSize( 0.7 );
    slice_unf->hist_->SetStats( false );

    slice_unf->hist_->SetTitle("");

    if (total_only) {
      slice_unf->hist_->GetXaxis()->SetLabelSize(0);
      slice_unf->hist_->GetXaxis()->SetTickLength(0);
    }

    double ymax = -DBL_MAX;
    slice_unf->hist_->GetYaxis()->SetRangeUser( 0., slice_unf->hist_->GetMaximum()*1.5 );
    slice_unf->hist_->Draw( "e" );

    slice_cv->hist_->SetStats( false );
    slice_cv->hist_->SetLineColor( kAzure - 7 );
    slice_cv->hist_->SetLineWidth( 3 );
    slice_cv->hist_->SetLineStyle( 5 );
    slice_cv->hist_->Draw( "hist same" );

    if ( using_fake_data ) {
      slice_truth->hist_->SetStats( false );
      slice_truth->hist_->SetLineColor( kOrange );
      slice_truth->hist_->SetLineWidth( 5 );
      slice_truth->hist_->Draw( "hist same" );
    }

    // Print Values
    for (int i = 1; i <= slice_unf->hist_->GetNbinsX(); ++i) {
        double bin_content = slice_unf->hist_->GetBinContent(i);
        double bin_error = slice_unf->hist_->GetBinError(i);
        std::cout << "Bin: " << i << ", Val = " << bin_content << ", Unc = " << bin_error << std::endl;
    }

    // Draw generator predictions
    // Ac matrix
    const TMatrixD &A_C =  *xsec.result_.add_smear_matrix_;
    // Find bins for this slice
    size_t start = std::numeric_limits<size_t>::max();
    size_t stop = 0;
    for (const auto& entry : slice.bin_map_) {
        const auto& set = entry.second;
        if (set.size() != 1) {
            throw std::runtime_error("Error: set in bin_map_ has more or less than 1 entry");
        }
        start = std::min(start, *set.begin());
        stop = std::max(stop, *set.rbegin());
    }

    TMatrixD ac_hist_slice(stop - start + 1, stop - start + 1);
    for (int i = start; i <= stop; i++) {
        for (int j = start; j <= stop; j++) {
            ac_hist_slice(i - start, j - start) = A_C.operator()(i, j);
        }
    }

    TLegend *lg = new TLegend(1 - rightMargin + 0.02, 0.09, 1 - 0.02, 0.93);

    // Add in generator predictions
    // loop through generators
    for(const auto& generator : generators) {
      int h_idx = sl_idx;
      if (total_only) h_idx = 4;
      const auto gen_hist = get_generator_hist(generator.path, h_idx, generator.scaling);
      if (gen_hist) {

        // Multiple by AC matrix
        if (!total_only) multiply_1d_hist_by_matrix(&ac_hist_slice, gen_hist);

        // Normalise by bin width
        for (int i = 1; i <= gen_hist->GetNbinsX(); ++i) {
            double bin_content = gen_hist->GetBinContent(i);
            double bin_error = gen_hist->GetBinError(i);
            double bin_width = slice_unf->hist_->GetXaxis()->GetBinWidth(i); // Get the bin width from the slice_unf histogram
            gen_hist->SetBinContent(i, bin_content / bin_width);
            gen_hist->SetBinError(i, bin_error / bin_width);
        }

        gen_hist->SetLineColor(generator.lineColor); // Set the line color
        gen_hist->SetLineWidth(generator.lineWidth); // Set the line width
        gen_hist->SetLineStyle(generator.lineStyle); // Set the line style
        gen_hist->Draw( "hist same" );

        // Calculate a Chi2 and p-value
        SliceHistogram *gen_slice_h = SliceHistogram::slice_histogram_from_histogram(*gen_hist);
        const auto& chi2_result = gen_slice_h->get_chi2( *slice_unf );
        //std::cout << chi2_result.chi2_ << ", p-value = " << chi2_result.p_value_ << std::endl;

        std::ostringstream oss;
        oss << "#splitline{" << generator.name << "}{"
            << "#chi^{2} = " << (chi2_result.chi2_>= 0.01 && chi2_result.chi2_ < 100 ? std::fixed : std::scientific) << std::setprecision(2) << chi2_result.chi2_ << " / " << chi2_result.num_bins_ << " bin" << (chi2_result.num_bins_ > 1 ? "s" : "") << "}";
        std::string label = oss.str();

        lg->AddEntry(gen_hist, label.c_str(), "l");
      }
    }

    for ( const auto& pair : slice_gen_map ) {
      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      const auto& chi2_result = chi2_map.at( name );

      std::string name_clean = name;

      if (name_clean == "truth") name_clean = "NuWro Truth";
      //if (name_clean == "truth") name_clean = "Truth";
      if (name_clean == "MicroBooNE Tune") name_clean = "GENIE 3.0.6 G18 #muB"; // _10a_02_11a
      //if (label == "unfolded data") label = "Unfolded Fake Data";
      if (name_clean == "unfolded data") name_clean = "Unfolded Data";

      std::ostringstream oss;

      if (name != "unfolded data") {
        oss << "#splitline{" << name_clean << "}{"
            << "#chi^{2} = " << (chi2_result.chi2_>= 0.01 && chi2_result.chi2_ < 100 ? std::fixed : std::scientific) << std::setprecision(2) << chi2_result.chi2_ << " / " << chi2_result.num_bins_ << " bin" << (chi2_result.num_bins_ > 1 ? "s" : "") << "}";
      }
      else {
        oss << name_clean;
      }
      std::string label = oss.str();

      if (name == "unfolded data") lg->AddEntry( slice_h->hist_.get(), label.c_str(), "lep" );
      else lg->AddEntry( slice_h->hist_.get(), label.c_str(), "l" );
    }

    // redraw data points to put on top
    slice_unf->hist_->Draw( "e same" );

    // Draw Legend
    lg->Draw( "same" );

     // Create the label text with the POT value
    TLatex label;
    label.SetTextAlign(12); // Set text alignment (left-aligned)
    label.SetNDC(); // Set position in normalized coordinates
    std::string labelText1( "MicroBooNE NuMI Data" );
    std::string labelText2( "2.2#times10^{20} POT" );
    label.SetTextSize(0.045);

    if (sl_idx == 0 && !total_only) {
      label.DrawLatex(0.4, 0.85, labelText1.c_str() );
      label.DrawLatex(0.4, 0.80, labelText2.c_str() );
    }
    else {
      label.DrawLatex(0.135, 0.85, labelText1.c_str() );
      label.DrawLatex(0.135, 0.80, labelText2.c_str() );
    }

    // write to file
    std::string plot_name = "../unfold_output/plot_slice_" + std::to_string(sl_idx) + ".pdf";
    c1->SaveAs(plot_name.c_str());

  }

  // create plot of A_C matrix
  const Int_t n = 2;
  Double_t bins[n+1] = {0, 7, 8};
  const Char_t *labels[n] = {"E_{e}", "Total"};

  // Convert TMatrixD to TH2D
  TMatrixD temp_ac = *xsec.result_.add_smear_matrix_;
  TH2D h_A_C = TMatrixDToTH2D(temp_ac, "h_A_C", "Regularization Matrix", 0, temp_ac.GetNcols(), 0, temp_ac.GetNrows());

  gStyle->SetPalette(kBird);

  TCanvas *c_ac = new TCanvas("c_ac","A_C Matrix",200,10,1920,1080);
  c_ac->SetRightMargin(0.15);
  c_ac->SetTopMargin(0.125);
  h_A_C.SetStats(0); // Disable the statistics box
  //h_A_C.GetZaxis()->SetRangeUser(-0.5, 1.5); // Set the z range
  h_A_C.Draw("colz");
  h_A_C.GetXaxis()->SetTitle("Bin Number");
  h_A_C.GetYaxis()->SetTitle("Bin Number");
  h_A_C.GetZaxis()->SetTitle("Regularization");
  c_ac->Update();

  // Draw vertical and horizontal lines at the bin edges
  for (Int_t i = 1; i < n; i++) {
      TLine *vline = new TLine(bins[i], 0, bins[i], h_A_C.GetNbinsY());
      vline->SetLineColor(kBlack);
      vline->Draw();

      TLine *hline = new TLine(0, bins[i], h_A_C.GetNbinsX(), bins[i]);
      hline->SetLineColor(kBlack);
      hline->Draw();
  }

  for (Int_t i = 1; i <= n; i++) {
      // Draw white dotted lines from bins[i-1] to bins[i]
      TLine *vline_dotted1 = new TLine(bins[i], bins[i-1], bins[i], bins[i]);
      vline_dotted1->SetLineColor(kWhite);
      vline_dotted1->SetLineStyle(2); // Set line style to dotted
      vline_dotted1->Draw();

      TLine *hline_dotted1 = new TLine(bins[i-1], bins[i], bins[i], bins[i]);
      hline_dotted1->SetLineColor(kWhite);
      hline_dotted1->SetLineStyle(2); // Set line style to dotted
      hline_dotted1->Draw();

      if(i<n)
      {
          TLine *vline_dotted2 = new TLine(bins[i], bins[i+1], bins[i], bins[i]);
          vline_dotted2->SetLineColor(kWhite);
          vline_dotted2->SetLineStyle(2); // Set line style to dotted
          vline_dotted2->Draw();

          TLine *hline_dotted2 = new TLine(bins[i+1], bins[i], bins[i], bins[i]);
          hline_dotted2->SetLineColor(kWhite);
          hline_dotted2->SetLineStyle(2); // Set line style to dotted
          hline_dotted2->Draw();
      }
  }

  // Add labels in the middle of the intervals
  for (Int_t i = 0; i < n; i++) {
      Double_t midPoint = (bins[i] + bins[i+1]) / 2.0;
      TLatex *text = new TLatex(midPoint, 1.03*h_A_C.GetNbinsY(), labels[i]);
      text->SetTextSize(0.03); // Set text size to something smaller
      text->SetTextAlign(22); // Center alignment
      text->Draw();
      // delete text;
  }

  // add labels in each bin
  for (int i = 0; i < h_A_C.GetNbinsX(); i++) {
      for (int j = 0; j < h_A_C.GetNbinsY(); j++) {

          double bin_content = h_A_C.GetBinContent(i+1, j+1);
          if (bin_content == 0) continue;

          TLatex* latex = new TLatex(h_A_C.GetXaxis()->GetBinCenter(i+1), h_A_C.GetYaxis()->GetBinCenter(j+1), Form("%.3f",bin_content));
          latex->SetTextFont(42);
          latex->SetTextSize(0.02);
          latex->SetTextAlign(22);
          latex->Draw();
      }
  }

  h_A_C.Write();

  c_ac->SaveAs("../unfold_output/plot_regularization_matrix.pdf");

}

int main( int argc, char* argv[] ) {

  if ( argc != 4 ) {
    std::cout << "Usage: Unfolder.C XSEC_Config"
      << " SLICE_Config OUTPUT_FILE\n";
    return 1;
  }

  std::string XSEC_Config( argv[1] );
  std::string SLICE_Config( argv[2] );
  std::string OutputFile( argv[3] );

  //Take the output directory from the file handed as the expected output
  //Only used for dumping to text or plot, if that option is requested in the hardcoded options at start of file
  std::string OutputDirectory = OutputFile.substr(0, OutputFile.find_last_of("/") + 1);
  std::string OutputFileName = OutputFile.substr(OutputFile.find_last_of("/") + 1);

  UnfolderNuMI(XSEC_Config, SLICE_Config, OutputDirectory, OutputFileName);
  return 0;
}
