#pragma once

// Standard library includes
#include <cmath>
#include <iostream>
#include <stdexcept>

// ROOT includes
#include "TVectorD.h"

// STV analysis includes
#include "Unfolder.hh"

// Implementation of the iterative D'Agostini unfolding method
// G. D'Agostini, Nucl. Instrum. Methods Phys. Res. A 362, 487-498 (1995)
// https://hep.physics.utoronto.ca/~orr/wwwroot/Unfolding/d-agostini.pdf.
//
// Uncertainties on the measurement are propagated through the unfolding
// procedure following a corrected expression given in Eq. (20) from
// https://arxiv.org/abs/1806.03350. The original paper by D'Agostini
// propagates the measurement uncertainties as if the unfolding matrix is
// independent of the measured data points, but this is only true for the first
// iteration.
class DAgostiniUnfolder : public Unfolder {

  public:

    // The maximum number of iterations to perform (to avoid infinite loops)
    static constexpr unsigned int DAGOSTINI_MAX_ITERATIONS = 99;

    // Enumerated type describing the rule to use to
    // stop the iterations
    enum ConvergenceCriterion {
      FixedIterations = 0, // Stop after a predefined iteration count
      FigureOfMerit = 1, // Stop once a figure of merit falls below threshold
    };

    inline DAgostiniUnfolder( unsigned int default_iterations = 1u )
      : Unfolder(), conv_criter_( ConvergenceCriterion::FixedIterations ),
      num_iterations_( default_iterations ) {}

    inline DAgostiniUnfolder( ConvergenceCriterion cc, double fig_merit_target )
      : Unfolder(), conv_criter_( cc ),
      num_iterations_( DAGOSTINI_MAX_ITERATIONS ),
      fig_merit_target_( fig_merit_target ) {}

    // Trick taken from https://stackoverflow.com/a/18100999
    using Unfolder::unfold;

    virtual UnfoldedMeasurement unfold( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal ) const override;

    inline unsigned int get_iterations() const { return num_iterations_; }
    inline void set_iterations( unsigned int iters )
      { num_iterations_ = iters; }

    // Sets the value of the flag which determines whether MC statistical
    // uncertainties on the response (smearceptance) matrix elements should
    // be propagated through the unfolding procedure. This is more correct,
    // but it is often a small effect and relatively slow to compute.
    inline void set_include_respmat_covariance( bool do_it )
      { include_respmat_covariance_ = do_it; }

  protected:

    // Calculates the "figure of merit" used to determine convergence of the
    // iterations when using the ConvergenceCriterion::FigureOfMerit option
    double calc_figure_of_merit( const TMatrixD& old_true_signal,
      const TMatrixD& new_true_signal ) const;

    ConvergenceCriterion conv_criter_;
    unsigned int num_iterations_;
    double fig_merit_target_;
    bool include_respmat_covariance_ = false;
};
