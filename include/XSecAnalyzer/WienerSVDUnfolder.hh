#pragma once

// ROOT includes
#include "TDecompChol.h"
#include "TDecompSVD.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Unfolder.hh"

// Implementation of the Wiener-SVD unfolding method
// W. Tang et al., J. Instrum. 12, P10002 (2017)
// https://arxiv.org/abs/1705.03568
class WienerSVDUnfolder : public Unfolder {

  public:

    enum RegularizationMatrixType { kIdentity, kFirstDeriv, kSecondDeriv };

    inline WienerSVDUnfolder( bool use_wiener_filter = true,
      RegularizationMatrixType type = kIdentity ) : Unfolder(),
      use_filter_( use_wiener_filter ), reg_type_( type ) {}

    // Trick taken from https://stackoverflow.com/a/18100999
    using Unfolder::unfold;

    virtual UnfoldedMeasurement unfold( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal ) const override;

    inline bool use_filter() const { return use_filter_; }
    inline void set_use_filter( bool use_filter ) { use_filter_ = use_filter; }

    inline RegularizationMatrixType get_regularization_type() const
      { return reg_type_; }

    inline void set_regularization_type( const RegularizationMatrixType& type )
      { reg_type_ = type; }

  protected:

    // Helper function that sets the contents of the regularization matrix
    // based on the current value of reg_type_
    void set_reg_matrix( TMatrixD& C ) const;

    // Flag indicating whether the Wiener filter should be used. If it is
    // false, then the usual expression will be replaced with an identity
    // matrix
    bool use_filter_ = true;

    // Enum that determines the form to use for the regularization matrix C
    RegularizationMatrixType reg_type_ = kIdentity;
};
