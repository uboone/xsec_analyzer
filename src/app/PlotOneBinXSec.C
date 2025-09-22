#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMatrixD.h"

#include "XSecAnalyzer/CrossSectionExtractor.hh"

void PlotOneBinXSec(const std::string& xsec_config, const std::string& output_name = "one_bin_xsec_plot.pdf") {
  gStyle->SetOptStat(0);

  // Load unfolded result
  auto extractor = std::make_unique<CrossSectionExtractor>(xsec_config);
  auto unfolded = extractor->get_unfolded_events();

  // Conversion factor from signal events to cross section [cm^2]
  double conv_factor = extractor->conversion_factor();

  // Assume only one bin in the unfolded signal
  const TMatrixD& unfolded_sig = *unfolded.result_.unfolded_signal_;
  double val_evt = unfolded_sig(0, 0);
  double val_xsec = val_evt * conv_factor;

  // Uncertainty from covariance matrix
  const TMatrixD& cov = *unfolded.result_.cov_matrix_;
  double var_evt = cov(0, 0);
  double err_xsec = std::sqrt(var_evt) * conv_factor;

  // Validation stanza
  TH1D* reco_data_hist = extractor->get_data_histogram();
  std::cout << "Raw selected data events (reco bins):" << std::endl;

  const auto& pred_map = extractor->get_prediction_map();
  TMatrixD pred = pred_map.at("MicroBooNETune")->get_prediction();

  std::cout << "GENIE truth prediction (true bins):" << std::endl;
  for (int b = 0; b < pred.GetNrows(); ++b) {
    std::cout << "  True bin " << b << ": " << pred(b, 0) << std::endl;
  }

  std::cout << "Cross section = " << val_xsec << " ± " << err_xsec << " cm^2\n";

  // Plot
  TCanvas* c = new TCanvas("c", "One-bin Cross Section", 800, 600);

  TH1D* h = new TH1D("h", "", 1, 0, 1);
  h->SetBinContent(1, val_xsec);
  h->SetBinError(1, err_xsec);

  h->SetMinimum(0);
  h->SetMaximum(1.5 * (val_xsec + err_xsec));
  h->SetTitle("");
  h->GetXaxis()->SetLabelSize(0);
  h->GetXaxis()->SetTickLength(0);
  h->GetXaxis()->SetTitle("Total Cross Section");
  h->GetYaxis()->SetTitle("#sigma [cm^{2}]");

  h->SetLineColor(kBlue+1);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.2);

  h->Draw("E1");

  // Legend
  TLegend* leg = new TLegend(0.6, 0.75, 0.88, 0.88);
  leg->AddEntry(h, "MicroBooNE Unfolded #sigma", "lep");
  leg->Draw();

  // Label
  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.04);
  label.DrawLatex(0.2, 0.85, "MicroBooNE Simulation");

  c->SaveAs(output_name.c_str());
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: PlotOneBinXSec XSEC_CONFIG [output_name.pdf]\n";
    return 1;
  }

  std::string xsec_config(argv[1]);
  std::string output_name = (argc > 2) ? argv[2] : "one_bin_xsec_plot.pdf";

  PlotOneBinXSec(xsec_config, output_name);
  return 0;
}
