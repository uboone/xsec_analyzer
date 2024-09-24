// Standard library includes
#include <iomanip>
#include <iostream>
#include <sstream>

#include <algorithm>

// ROOT includes
#include "TCanvas.h"
#include "TLegend.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/CrossSectionExtractor.hh"
#include "XSecAnalyzer/PGFPlotsDumpUtils.hh"
#include "XSecAnalyzer/SliceBinning.hh"
#include "XSecAnalyzer/SliceHistogram.hh"

//Useful DEBUG options which can be turned on/off
std::string PlotExtension = ".pdf";
std::string TextExtension = ".txt";
bool DumpToText = false;
bool DumpToPlot = false;

void Unfolder(std::string XSEC_Config, std::string SLICE_Config, std::string OutputDirectory, std::string OutputFileName) {

  std::cout << "\nRunning Unfolder.C with options:" << std::endl;
  std::cout << "\tXSEC_Config: " << XSEC_Config << std::endl;
  std::cout << "\tSLICE_Config: " << SLICE_Config << std::endl;
  std::cout << "\tOutputDirectory: " << OutputDirectory << std::endl;
  std::cout << "\tOutputFileName: " << OutputFileName << std::endl;
  std::cout << "\n" << std::endl;

  // Use a CrossSectionExtractor object to handle the systematics and unfolding
  auto extr = std::make_unique< CrossSectionExtractor >( XSEC_Config );

  // Plot slices of the unfolded result
  auto* sb_ptr = new SliceBinning( SLICE_Config );
  auto& sb = *sb_ptr;

  auto xsec = extr->get_unfolded_events();
  double conv_factor = extr->conversion_factor();
  const auto& pred_map = extr->get_prediction_map();

  std::cout << "\n\nSaving results -----------------" << std::endl;
  OutputFileName = OutputDirectory+"/"+OutputFileName;

  std::cout << "Output file - " << OutputFileName << std::endl;
  if (DumpToText) std::cout << "\tDumping plots to " << TextExtension << " files" << std::endl;
  if (DumpToPlot) std::cout << "\tDumping plots to " << PlotExtension << " files" << std::endl;
  std::cout << "\n" << std::endl;

  // Check that the output file can be written to
  TFile* File = new TFile(OutputFileName.c_str(), "RECREATE");
  if (!File || File->IsZombie()) {
    std::cerr << "Could not write to output file:" << OutputFileName << std::endl;
    throw;
  }

  // Make results in both Event Count units and then the XSec units
  std::vector<std::string> ResultTypes(2);
  ResultTypes[0] = "EventCountUnits";
  ResultTypes[1] = "XsecUnits";
  size_t nResultTypes = ResultTypes.size();

  //Loop over the two ResultTypes (Event Counts and Xsec Units)
  for (size_t iRT=0;iRT<nResultTypes;iRT++) {

    std::string RT = ResultTypes[iRT];
    File->cd();
    File->mkdir(RT.c_str());
    File->cd(RT.c_str());

    //======================================================================================
    //Loop over all covariance matrices stored in the unfolded result and save them

    File->cd(RT.c_str());
    File->mkdir((RT+"/Covariances").c_str());
    File->cd((RT+"/Covariances").c_str());

    // Convert units on the covariance matrices one-by-one and dump them
    for ( const auto& cov_pair : xsec.unfolded_cov_matrix_map_ ) {
      const auto& name = cov_pair.first;
      TMatrixD temp_cov_matrix = *cov_pair.second;
      // Note that we need to square the unit conversion factor for the
      // covariance matrix elements

      if (RT == "XsecUnits") {
	temp_cov_matrix *= std::pow( 1.0 / conv_factor, 2 );
      }

      temp_cov_matrix.Write(name.c_str());
      if (DumpToText) dump_text_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_" + name + TextExtension, temp_cov_matrix );
      if (DumpToPlot) draw_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_" + name + PlotExtension, temp_cov_matrix, (name+" matrix").c_str(), "Bin Number", "Bin Number", "COLZ");
    }

    // No unit conversions are necessary for the unfolding, error propagation,
    // and additional smearing matrices since they are dimensionless
    TMatrixD temp_unfolding_matrix = *xsec.result_.unfolding_matrix_;
    TMatrixD temp_err_prop_matrix = *xsec.result_.err_prop_matrix_;
    TMatrixD temp_add_smear_matrix = *xsec.result_.add_smear_matrix_;
    temp_unfolding_matrix.Write("UnfoldingMatrix");
    temp_err_prop_matrix.Write("ErrorPropagationMatrix");
    temp_add_smear_matrix.Write("AdditionalSmearingMatrix");

    if (DumpToText) {
      dump_text_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_unfolding"+TextExtension, temp_unfolding_matrix );
      dump_text_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_err_prop"+TextExtension, temp_err_prop_matrix );
      dump_text_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_add_smear"+TextExtension, temp_add_smear_matrix );
    }
    if (DumpToPlot) {
      draw_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_unfolding"+PlotExtension, temp_unfolding_matrix, "Unfolding matrix", "Bin Number", "Bin Number", "COLZ");
      draw_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_err_prop"+PlotExtension, temp_err_prop_matrix, "Error Propagation matrix", "Bin Number", "Bin Number", "COLZ");
      draw_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_add_smear"+PlotExtension, temp_add_smear_matrix, "Addition Smearing matrix", "Bin Number", "Bin Number", "COLZ");
    }

    //======================================================================================
    //Loop over all the slices taken from the config and save the unfolded distribution/generator prediction

    for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
      auto& Slice = sb.slices_.at( sl_idx );
      //The following line will fall over if several active variables are used per Slice
      auto& SliceVar = sb.slice_vars_.at( Slice.active_var_indices_.front() );

      std::string SliceVariableName = SliceVar.name_;
      SliceVariableName.erase(std::remove(SliceVariableName.begin(), SliceVariableName.end(), ' '), SliceVariableName.end());

      File->cd(RT.c_str());
      File->mkdir((RT+"/"+SliceVariableName).c_str());
      File->cd((RT+"/"+SliceVariableName).c_str());

      //======================================================================================
      //Save unfolded distribution for each covariance matrix

      // Make a histogram showing the unfolded counts in the current slice
      // for a particular covariance matrix being used to define the uncertainties
      for ( const auto& uc_pair : xsec.unfolded_cov_matrix_map_ ) {
	const auto& uc_name = uc_pair.first;
	const auto& uc_matrix = uc_pair.second;

	SliceHistogram* Slice_unf = SliceHistogram::make_slice_histogram( *xsec.result_.unfolded_signal_, Slice, uc_matrix.get() );
	TH1* SliceHist = Slice_unf->hist_.get();
	if (RT == "XsecUnits") {
	  SliceHist->Scale(1.0 / conv_factor);
	}
	SliceHist->Write((SliceVariableName+"_"+uc_name).c_str());

	if (DumpToText) dump_text_column_vector( OutputDirectory+"/"+RT+"_vec_table_unfolded_signal_"+uc_name+TextExtension, *xsec.result_.unfolded_signal_ );
	if (DumpToPlot) draw_column_vector( OutputDirectory+"/"+RT+"_vec_table_unfolded_signal_"+uc_name+PlotExtension, *xsec.result_.unfolded_signal_, "Unfolded Signal", "Bin Number", "Cross Section [#times 10^{-38} cm^{2}]");
      }

      //======================================================================================
      //Loop over all generator predictions and save them to the same output

      /*
      //DB Still need to check this loop as I don't currently have generator prediction files for tutorial binning scheme
      for ( const auto& gen_pair : extr->get_prediction_map()) {
	std::string gen_short_name = gen_pair.second->name();
	TMatrixD temp_gen = gen_pair.second->get_prediction();
	if (RT == "XsecUnits") {
	  temp_gen *= (1.0 / conv_factor);
	}

	TH1D* temp_gen_hist = Matrix_To_TH1(temp_gen,gen_short_name,SliceVariableName,"Events");
	temp_gen_hist->Write(("GenPred_"+SliceVariableName+"_"+gen_short_name).c_str());

	if (DumpToText) dump_text_column_vector( OutputDirectory+"/"+RT+"_vec_table_" + gen_short_name + TextExtension, temp_gen );
	if (DumpToPlot) draw_column_vector( OutputDirectory+"/"+RT+"_vec_table_" + gen_short_name + PlotExtension, temp_gen, (gen_short_name + " Prediction").c_str(), "Bin Number", "Cross Section [#times 10^{-38} cm^{2}]");
      }
      */

    }
  }

  File->Close();
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

  Unfolder(XSEC_Config, SLICE_Config, OutputDirectory, OutputFileName);
  return 0;
}
