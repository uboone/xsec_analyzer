import ROOT

# Define categories, titles, and colors as structured data
categories = [
    ("kUnknown",         "Unknown",               ROOT.kGray),
    ("kNuMuCC0p0pi0n",   "CC#nu_{#mu}0p0#pi0n",    ROOT.kMagenta),
    ("kNuMuCC0p0pi1n",   "CC#nu_{#mu}0p0#pi1n",    ROOT.kGreen),
    ("kNuMuCC1p0pi0n",   "CC#nu_{#mu}1p0#pi0n",    ROOT.kMagenta + 2),
    ("kNuMuCC1p0pi1n",   "CC#nu_{#mu}1p0#pi1n",    ROOT.kGreen + 2),
    ("kNuMuCC0p0piNn",   "CC#nu_{#mu}0p0#piNn",    ROOT.kBlue - 4),
    ("kNuMuCC1p0piNn",   "CC#nu_{#mu}1p0#piNn",    ROOT.kBlue),
    ("kNuMuCCNp0pi0n",   "CC#nu_{#mu}Np0#pi0n",    ROOT.kMagenta + 4),
    ("kNuMuCCNp0pi1n",   "CC#nu_{#mu}Np0#pi1n",    ROOT.kGreen + 4),
    ("kNuMuCCNp0piNn",   "CC#nu_{#mu}Np0#piNn",    ROOT.kBlue + 2),
    ("kNuMuCCNpi",       "CC#nu_{#mu}N#pi",        ROOT.kOrange - 4),
    ("kNuMuCCOther",     "CC#nu_{#mu}Other",       ROOT.kOrange - 3),
    ("kNuECC",           "CC#nu_{e}",              ROOT.kRed - 7),
    ("kNC",              "NC",                     ROOT.kOrange),
    ("kOOFV",            "OOFV",                   ROOT.kRed + 3),
    ("kOther",           "Other",                  ROOT.kRed + 1),
]

# Create histograms from definitions
histos = []
for name, title, color in categories:
    h = ROOT.TH1F(name, "", 1, 0, 1)
    h.SetLineColor(color)
    h.SetFillColor(color)
    histos.append((h, title))

# Setup canvas and legend
c = ROOT.TCanvas()
leg = ROOT.TLegend(0.05, 0.05, 0.95, 0.95)
leg.SetLineColor(ROOT.kBlack)

for h, title in histos:
    leg.AddEntry(h, title, "f")

leg.Draw()
c.SaveAs("XSEC_Legend.png")
c.SaveAs("XSEC_Legend.pdf")
input("Press Enter to exit...")
c.Clear()

