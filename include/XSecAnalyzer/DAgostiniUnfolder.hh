#pragma once

// XSecAnalyzer includes
#include "Unfolder.hh"

// Implementation of the iterative D'Agostini unfolding method
// G. D'Agostini, Nucl. Instrum. Methods Phys. Res. A 362, 487-498 (1995)
// https://hep.physics.utoronto.ca/~orr/wwwroot/Unfolding/d-agostini.pdf.
//
// Two approaches to propagating the uncertainties are implemented, with
// the first being the default. One can switch between these by setting
// the member boolean variable use_A_C_ to true (option #1) or false
// (option #2). The default option #1 is recommended.
//
//   (1) The unfolding matrix is used to transform both the central-value
//       data points and the input reconstructed-space covariance matrix.
//       To avoid the need for computing additional unfolding-related
//       uncertainties when assessing goodness-of-fit, a "regularization
//       matrix" or "additional smearing matrix" A_C is also computed
//       and included in the output. This matrix, originally defined
//       for the Wiener-SVD unfolding technique, should be applied to
//       a vector of theoretical predictions before a quantitative comparison
//       (e.g., evaluation of the chi-squared statistic) is made with the
//       unfolded data. Since the A_C matrix accounts for unfolding-related
//       bias, the dependence of the unfolding matrix on the
//       measured data for iterations beyond the first is ignored.
//
//   (2) While the same unfolding matrix from option #1 is used to transform
//       the central-value data points, a distinct "error propagation matrix"
//       is used to transform the reconstructed-space covariance matrix. Under
//       this approach, a regularization matrix A_C is not computed and should
//       not be used. In this case, the error propagation matrix differs from
//       the unfolding matrix by some additional terms that are nonzero for
//       iterations beyond the first. These terms account for the dependence of
//       the unfolding matrix on the measured data points themselves. When A_C
//       is not used, the unfolding matrix represents an estimate of the
//       transformation from reconstructed to true space. Uncertainties on the
//       original measurement thus lead to an uncertainty on the unfolding
//       matrix itself. This dependence is accounted for in the code using
//       the expression given in Eq. (20) from https://arxiv.org/abs/1806.03350.
//       Note also that the original paper by D'Agostini propagates the
//       measurement uncertainties as if the unfolding matrix is independent of
//       the measured data points (similarly to option #1). However, when A_C
//       is not used (as is the case for D'Agostini's original treatment),
//       this is only true for the first iteration.
//
//   Further discussion of related issues for D'Agostini unfolding can be found
//   in https://arxiv.org/abs/2401.04065.
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

    // Sets the value of the flag that determines whether A_C should be
    // computed and stored in the output. A value of true corresponds to
    // the "option #1" approach to uncertainty propagation described above.
    // A value of false corresponds to option #2.
    inline void set_use_AC( bool use_it ) { use_A_C_ = use_it; }

    // Sets the value of the flag that determines whether MC statistical
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
    bool use_A_C_ = true;
};
