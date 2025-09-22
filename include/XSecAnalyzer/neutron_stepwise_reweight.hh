// neutron_stepwise_reweight.hh
#pragma once

#include <vector>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include "TH1D.h"
#include "TVector3.h"
#include "TFile.h"
#include <random>

// Global constant for inverse nucleon number density (cm^3 / nucleon)
const double inverse_n = 39.95 / (1.394 * 6.023e23);

// Container holding pointers to nominal and MiniCAPTAIN cross section data
struct CrossSectionData {
  const TH1D* scaled_hist;
  const TH1D* nominal_hist;
  const TGraphAsymmErrors* miniCaptainGraph;
};

// Global pointer to the neutron cross section histogram
inline CrossSectionData LoadXSecHistogramOnce(const std::string& filename) {
  static std::unique_ptr<TH1D> scaled_hist_ptr;
  static std::unique_ptr<TH1D> nominal_hist_ptr;
  static std::unique_ptr<TGraphAsymmErrors> graph_ptr;

  if (!scaled_hist_ptr || !nominal_hist_ptr || !graph_ptr) {
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
      throw std::runtime_error("Failed to open xsec file: " + filename);
    }

    // Scaled CV histogram (green curve)
    TH1D* scaled = dynamic_cast<TH1D*>(file->Get("scaled_xsec_histogram"));
    if (!scaled) {
      throw std::runtime_error("Failed to find scaled_xsec_histogram in file: " + filename);
    }
    scaled->SetDirectory(nullptr);
    scaled_hist_ptr.reset(scaled);

    // Nominal GEANT4 histogram (red curve)
    TH1D* nominal = dynamic_cast<TH1D*>(file->Get("xsec_histogram"));
    if (!nominal) {
      throw std::runtime_error("Failed to find xsec_histogram in file: " + filename);
    }
    nominal->SetDirectory(nullptr);
    nominal_hist_ptr.reset(nominal);

    // MiniCAPTAIN data graph
    TGraphAsymmErrors* graph = dynamic_cast<TGraphAsymmErrors*>(file->Get("miniCaptain_graph"));
    if (!graph) {
      throw std::runtime_error("Failed to find miniCaptain_graph in file: " + filename);
    }
    graph_ptr.reset(graph);

    file->Close();
    delete file;
  }

  return {scaled_hist_ptr.get(), nominal_hist_ptr.get(), graph_ptr.get()};
}

// Describes a single neutron segment's passage through detector
struct Segment {
  double KE;             // kinetic energy in MeV
  double length;         // segment length in cm
  std::string end_process; // final process (e.g. neutronInelastic, protonInelastic, NCapture)
};

// Holds the results of recursive tracking
struct NeutronTrackResult {
  std::vector<Segment> segments; // neutron path broken into segments
  std::string final_status;  // how the neutron ends
};

// Get cross section from histogram
double get_geant4_bin_xsec(const TH1D* hist, double KE) {
  if (!hist) throw std::runtime_error("Null xsec histogram");
  int bin = hist->FindFixBin(KE);

  // Ensure KE is within valid histogram range
  if (bin <= 0 || bin > hist->GetNbinsX()) {
    std::cout << "[DEBUG] KE=" << KE << " is outside valid xsec histogram range." << std::endl;
    return 0.0;
  }
  double xsec = hist->GetBinContent(bin);
  return std::max(0.0, xsec);
}

// Interpolate MiniCAPTAIN cross section at given KE
double get_minicaptain_xsec(const TGraphAsymmErrors* graph, double KE) {
  if (!graph) throw std::runtime_error("Null minicaptain graph");

  int n = graph->GetN();
  const double* x = graph->GetX();

  // Find closest KE point
  int closest = 0;
  double min_dist = std::abs(KE - x[0]);
  for (int i = 1; i < n; ++i) {
    double dist = std::abs(KE - x[i]);
    if (dist < min_dist) {
      min_dist = dist;
      closest = i;
    }
  }

  // Clamp to last bin if KE is beyond the final point
  if (KE > x[n - 1]) closest = n - 1;

  double y = graph->GetY()[closest];
  return std::max(0.0, y);
}

// Perform stepwise reweighting for a single neutron track (used in the multisim for neutron reinteraction uncertainty calculation)
// inputs are: 	- vector of neutron segments (each time a neutron scatters in the detector volume, it becomes a new segment)
// 		- the scaled central value cross section histogram (sometimes called upon multiple times because the neutron can change energy significantly as it rescatters)
// 		- Minicaptain data error bars (used to set the gaussian width)
// 		- random number
// 		- step length
// output is: event weight
inline double stepwise_weight(const Segment& seg,
                              const TH1D* nominal_xsec,                 // scaled CV curve
                              const TGraphAsymmErrors* minicaptain_data,
                              double z_fixed,
                              double dL = 0.1)
{
  // Trivial rejects
  if (seg.KE <= 0.0 || seg.KE < 0.1 || seg.length <= 0.0 || dL <= 0.0) return 1.0;

  // Lookup scaled CV sigma at this KE
  double nom_xsec = get_geant4_bin_xsec(nominal_xsec, seg.KE);
  if (!std::isfinite(nom_xsec) || nom_xsec < 0.0) nom_xsec = 0.0;

  // Build the gaussian throw width around the scaled CV using Mini-CAPTAIN band
  double sigma_throw = 0.0;
  double mini_xsec = 0.0;
  if (minicaptain_data && minicaptain_data->GetN() > 0) {
    const int n = minicaptain_data->GetN();
    const double* x       = minicaptain_data->GetX();
    const double* ey_high = minicaptain_data->GetEYhigh();
    const double* ey_low  = minicaptain_data->GetEYlow();

    // nearest MCaptain point
    int closest = 0;
    double min_dist = std::abs(seg.KE - x[0]);
    for (int i = 1; i < n; ++i) {
      const double d = std::abs(seg.KE - x[i]);
      if (d < min_dist) { min_dist = d; closest = i; }
    }

    mini_xsec = get_minicaptain_xsec(minicaptain_data, seg.KE);
    const double mini_err = 0.5*(ey_high[closest] + ey_low[closest]);  // symmetric approx
    sigma_throw = 2.0 * std::abs(mini_xsec - nom_xsec) + mini_err;
    if (!std::isfinite(sigma_throw) || sigma_throw < 0.0) sigma_throw = 0.0;
  }

  // Create varied cross section around the nominal
  double var_xsec = nom_xsec + z_fixed * sigma_throw;

  const double delta_sigma = var_xsec - nom_xsec;

  // --- Integrate along the path without a vanishing micro-step ---
  constexpr double tiny = 1e-12;   // cm
  constexpr double epsP = 1e-300;  // guard against 0 denominator

  const double L = seg.length;
  const int    N_full = static_cast<int>(std::floor(L / dL + tiny));
  double used  = static_cast<double>(N_full) * dL;
  double rem   = L - used;
  if (rem < tiny) rem = 0.0;

  // Do survival over (N_full - 1) full steps; apply interaction (or survival) on the final step
  const int    n_surv_steps = std::max(0, N_full - 1);
  const double last_step    = (N_full > 0 ? dL : 0.0) + rem;

  double total_w = 1.0;

  // Survival over the first chunk (vectorized)
  if (n_surv_steps > 0) {
    const double f_step = std::exp(-delta_sigma * dL / inverse_n);
    if (std::isfinite(f_step)) total_w *= std::pow(f_step, n_surv_steps);
  }

  // Final step
  if (last_step > tiny) {
    const bool interacted = (seg.end_process == "neutronInelastic");

    if (interacted) {
      // Numerically-stable interaction probabilities
      const double x_nom =  nom_xsec   * last_step / inverse_n;
      const double x_var =  var_xsec   * last_step / inverse_n;
      const double P_nom = -std::expm1(-x_nom);   // 1 - exp(-x_nom)
      const double P_var = -std::expm1(-x_var);

      if (!(P_nom == 0.0 && P_var == 0.0)) {
        const double denom = (P_nom == 0.0 ? epsP : P_nom);
        const double w_int = P_var / denom;
        if (std::isfinite(w_int)) total_w *= w_int;
        // else: neutralize this step instead of poisoning the weight
      }
      // If both P's are exactly zero, treat as neutral (no ratio info).
    } else {
      const double f_last = std::exp(-delta_sigma * last_step / inverse_n);
      if (std::isfinite(f_last)) total_w *= f_last;
    }
  }

  return total_w;
}

// Perform stepwise reweighting for a single neutron track (used to scale the nominal GEANT4 neutron-argon total inelastic cross section up to the scaled CV curve shown in docdb 44856-v2 slide 7)
// inputs are:  - vector of neutron segments (each time a neutron scatters in the detector volume, it becomes a new segment)
//              - the central value cross section histogram (red curve in docdb 44856-v2 slide 7)
//              - the scaled central value cross section histogram (green curve in docdb 44856-v2 slide 7)
//              - step length
// output is: event weight
inline double neutron_cv_scale_weight(const std::vector<Segment>& segments,
                                      const TH1D* geant4_total_xsec,   // red
                                      const TH1D* scaled_CV_xsec,      // green
                                      double dL = 0.1) {
  if (!geant4_total_xsec || !scaled_CV_xsec || segments.empty() || dL <= 0.0) return 1.0;

  const Segment& seg = segments.front();
  if (seg.length <= 0.0) return 1.0;

  // Cross-section lookups at this KE
  const double sigma_red   = get_geant4_bin_xsec(geant4_total_xsec, seg.KE);
  const double sigma_green = get_geant4_bin_xsec(scaled_CV_xsec,    seg.KE);

  // Numerical knobs
  constexpr double tiny = 1e-12;      // cm; treat anything below as zero-length
  constexpr double epsP = 1e-300;     // guard for 0 denominator in ratios

  // Partition L into full steps + remainder, but fold a tiny remainder
  // into the last full step so there is no micro-step at the end.
  const double L = seg.length;
  const int    N_full = static_cast<int>(std::floor(L / dL + tiny)); // bias up a hair
  double used  = static_cast<double>(N_full) * dL;
  double rem   = L - used;
  if (rem < tiny) rem = 0.0;

  // Steps:
  // - First do survival over (N_full - 1) full steps of size dL
  // - Final step length = (N_full > 0 ? dL : 0) + rem. This is to guard for the edge
  // case where the neutron travels for less than a full dL
  const int    n_surv_steps = std::max(0, N_full - 1);
  const double last_step    = (N_full > 0 ? dL : 0.0) + rem;

  double total_w = 1.0;

  // Survival over the first chunk (vectorized via pow for stability)
  if (n_surv_steps > 0) {
    const double delta  = sigma_green - sigma_red;
    const double f_step = std::exp(-delta * dL / inverse_n);
    total_w *= std::pow(f_step, n_surv_steps);
  }

  // Final step: either interaction ratio (if inelastic) or one more survival
  if (last_step > tiny) {
    const bool interacted = (seg.end_process == "neutronInelastic");

    if (interacted) {
      // Use expm1 for numerical stability at tiny arguments
      const double x_red =  sigma_red   * last_step / inverse_n;
      const double x_grn =  sigma_green * last_step / inverse_n;
      const double P_red = -std::expm1(-x_red);  // == 1 - exp(-x_red)
      const double P_grn = -std::expm1(-x_grn);

      // If both are ~0 (no interaction per unit length in both models),
      // treat as neutral instead of 0/0 -> NaN.
      if (!(P_red == 0.0 && P_grn == 0.0)) {
        const double denom = (P_red == 0.0 ? epsP : P_red);
        const double w     = P_grn / denom;
        if (std::isfinite(w)) total_w *= w;
        // else: neutralize this step (skip) to avoid poisoning the event
      }
    } else {
      const double delta_last = sigma_green - sigma_red;
      total_w *= std::exp(-delta_last * last_step / inverse_n);
    }
  }

  //if (seg.KE > 0.3 && seg.KE < 0.35) {std::cout << "Neutron between 300 and 350 MeV. Track Weight: " << total_w << ", path length: " << seg.length << std::endl; }

  return total_w;
}

// This function is strictly for validation purposes
// This function calculates the neutron weight for each individual segment of a neutrons path. Each time a neutron interacts, a new weight is calculated for the next section of the neutrons path.
// This is actually what is happening in both of the reweighting functions directly above, however for each neutron, the product of all the individual segments is taken at the end of the path
// and the resultant weight is assigned to the neutron. Neutrons scatter multiple times on their way out of the detector, so each section of its travel gets a different weight depending on the energy
// and path length of the neutron. To validate that the weight are being calcualted properly, you can plot the track length for the nominal neutron xsec and the scaled neutron xsec and confirm that
// the mean free path changes characteristically with the xsec scaling
//
// Inputs: - vector of neutron segments (each time a neutron scatters in the detector volume, it becomes a new segment)
// 	   - the central value cross section histogram (red curve in docdb 44856-v2 slide 7)
// 	   - the scaled central value cross section histogram (green curve in docdb 44856-v2 slide 7)
// 	   - step length
// Output: - 2 vectors. The first contains the original segment lengths, and the second contains the scaled segment length
inline std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>
segment_length_scaling(const std::vector<Segment>& segments,
                       const TH1D* nominal_xsec_hist,     // red
                       const TH1D* scaled_xsec_hist,      // green
                       double dL = 0.1) {
  std::vector<double> nominal_lengths;
  std::vector<double> scaled_lengths;
  std::vector<double> kinetic_energies;
  std::vector<double> nom_xsecs;
  std::vector<double> scaled_xsecs;
  std::vector<double> weights;

  for (const auto& seg : segments) {
    if (seg.KE <= 0) continue;
    if (seg.KE < 0.1) break; // below 100 MeV, stop

    double L_remaining = seg.length;
    double sigma_nom   = get_geant4_bin_xsec(nominal_xsec_hist, seg.KE);
    double sigma_scaled= get_geant4_bin_xsec(scaled_xsec_hist, seg.KE);

    if (sigma_nom <= 0) continue;

    // Loop through steps within this segment
    double total_seg_weight = 1.0;
    while (L_remaining > 0) {
      double step = std::min(dL, L_remaining);
      L_remaining -= step;

      bool is_final    = (L_remaining <= 0);
      bool interacted  = is_final && (seg.end_process == "neutronInelastic");

      if (interacted) {
        double P_int_nom    = 1.0 - std::exp(-sigma_nom    * step / inverse_n);
        double P_int_scaled = 1.0 - std::exp(-sigma_scaled * step / inverse_n);
        double w = P_int_scaled / P_int_nom;
        total_seg_weight *= w;
      } else {
        double delta = sigma_scaled - sigma_nom;
        total_seg_weight *= std::exp(-delta * step / inverse_n);
      }
    }

    // Store lengths
    nominal_lengths.push_back(seg.length);
    scaled_lengths.push_back(seg.length * total_seg_weight);
    kinetic_energies.push_back(seg.KE);
    nom_xsecs.push_back(sigma_nom);
    scaled_xsecs.push_back(sigma_scaled);
    weights.push_back(total_seg_weight);
  }

  return {nominal_lengths, scaled_lengths, kinetic_energies, nom_xsecs, scaled_xsecs, weights};
} 

// Determines whether a point is inside the detector volume
inline bool isContained(const TVector3& pos) {
  return (pos.X() > 0 && pos.X() < 256.35 &&
          pos.Y() > -116.5 && pos.Y() < 116.5 &&
          pos.Z() > 0 && pos.Z() < 1036.8);
}

// Calculates the length from the neuton vertex to the detector boundary that the neutron exits through
inline double lengthToEdge(const TVector3& start, const TVector3& end) {
  TVector3 dir = (end - start).Unit();
  double step = 0.1;  // cm
  TVector3 probe = start;
  double traveled = 0.0;

  while (isContained(probe)) {
    probe += dir * step;
    traveled += step;
    if (traveled > 2000.0) break; // failsafe
  }

  return traveled;
}

// Recursively tracks a neutron and its inelastic descendants
NeutronTrackResult TrackNeutron(int start_index,
                                const std::vector<int>& pdg,
                                const std::vector<TVector3>& start,
                                const std::vector<TVector3>& end,
                                const std::vector<double>& energy,
                                const std::vector<std::string>& end_process) {

  NeutronTrackResult result;

  if (pdg[start_index] != 2112) {
    std::cout << "[WARNING] particle submitted to TrackNeutron not a neutron" << std::endl;
    return result;
  }

  constexpr double kM_n = 0.939565; // GeV
  const bool start_in = isContained(start[start_index]);
  const bool end_in   = isContained(end[start_index]);

  const double KE = energy[start_index] - kM_n; // GeV
  double len_in = 0.0;
  std::string label;

  if (start_in && end_in) {
    // Fully inside
    len_in = (end[start_index] - start[start_index]).Mag();
    label  = end_process[start_index];
  } else if (start_in && !end_in) {
    // Inside -> outside: distance to boundary along this segment
    len_in = lengthToEdge(start[start_index], end[start_index]);
    label  = "escaped detector";
  } else if (!start_in && end_in) {
    // Outside -> inside: chord from boundary to the inside endpoint
    len_in = lengthToEdge(end[start_index], start[start_index]);
    label  = end_process[start_index];
  } else {
    // Outside -> outside: contributes no in-detector path
    label  = "never entered detector";
  }

  if (len_in > 0.0) {
    result.segments.push_back({KE, len_in, label});
  }
  result.final_status = label;

  return result;
}
