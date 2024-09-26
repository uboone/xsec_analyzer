#pragma once

// Standard library includes
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>

// ROOT includes
#include "TDecompQRH.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"

constexpr double DEFAULT_MATRIX_INVERSION_TOLERANCE = 1e-4;

std::unique_ptr< TMatrixD > invert_matrix( const TMatrixD& mat,
  const double inversion_tolerance = DEFAULT_MATRIX_INVERSION_TOLERANCE );

void dump_text_matrix( const std::string& output_file_name,
  const TMatrixD& matrix );

void dump_text_column_vector( const std::string& output_file_name,
  const TMatrixD& matrix );

void draw_column_vector( const std::string& output_file_name, const TMatrixD& matrix, const std::string& matrix_title, const std::string& xaxis_title="", const std::string& yaxis_title="", const std::string& draw_options="" );

void draw_matrix( const std::string& output_file_name, const TMatrixD& matrix, const std::string& matrix_title, const std::string& xaxis_title="", const std::string& yaxis_title="", const std::string& draw_options="" );

TH1D* Matrix_To_TH1(TMatrixD matrix, std::string matrix_title, std::string xaxis_title, std::string yaxis_title);

// Load a TMatrixD object saved in a text file by a previous call to
// dump_text_matrix() or dump_text_column_vector()
TMatrixD load_matrix( const std::string& input_file_name );

// Compute the direct sum of a vector of input TMatrixD objects
TMatrixD direct_sum( const std::vector< const TMatrixD* >& matrices );

// Overloaded version for a pair of input matrices
TMatrixD direct_sum( const TMatrixD& m1, const TMatrixD& m2 );
