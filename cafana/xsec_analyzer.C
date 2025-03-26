/**
 * @file example.C
 * @brief Convert flat CAFs into xsec_analyzer format.
 * @details This macro slices CAFAna files for use in xsec_analyzer. The macro configures
 * some basic variables and cuts, then runs the analysis over a single sample.
 * No cuts are applied as this is done upstream automatically.
 * @author brindenc@fnal.gov
*/

#include "include/pandora/slice_variables.h"
#include "include/pandora/slice_cuts.h"
#include "include/pandora/definitions.h"
#include "include/pandora/syst.h"
#include "include/analysis.h"
#include "include/preprocessor.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"


std::string out_fname = "/exp/sbnd/data/users/brindenc/xsec_analysis/test/sbruce_med/test1/xsec_analyzer_med";

void xsec_analyzer()
{
    /**
     * @brief Create an instance of the Analysis class.
     * @details This creates an instance of the Analysis class, which is used
     * to run the analysis on the specified samples. The name of the analysis,
     * and therefor the name of the output file, is specified as an argument to
     * the constructor.
     */
    //ana::Analysis analysis("/exp/sbnd/data/users/brindenc/xsec_analysis/test/sbruce/xsec_analyzer");
    ana::Analysis analysis(out_fname);

    /**
     * @brief Add a sample to the analysis.
     * @details This adds a sample to the analysis by creating a SpectrumLoader
     * object and adding it to the Analysis class. The SpectrumLoader object
     * represents the sample in the analysis, and is used to load the data from
     * the ROOT file and apply the cuts and variables. The name passed to the
     * AddLoader function is used to create a directory in the output ROOT file
     * to store the results of the analysis.
     */
    //ana::SpectrumLoader var00("/exp/sbnd/data/users/brindenc/xsec_analysis/test/bnb_small/small_syst.root");
    ana::SpectrumLoader var00("/exp/sbnd/data/users/brindenc/ML/bnb_cosmics/v09_89_01/syst/caf/caf_bnb_cosmics_dirt_0.flat.caf.root");
    analysis.AddLoader("nominal", &var00, true);
    
    /**
     * @brief Add a set of variables for slices to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    std::map<std::string, ana::SpillMultiVar> slice;

    // Double variables
    slice.insert({"leading_muon_momentum", ana::SpillMultiVar(PANDORAVAR(slicevars::kLeadingMuonMom,slicecuts::kNoCut))});
    slice.insert({"leading_muon_costheta", ana::SpillMultiVar(PANDORAVAR(slicevars::kLeadingMuonCostheta,slicecuts::kNoCut))});
    slice.insert({"true_leading_muon_momentum", ana::SpillMultiVar(PANDORAVAR(slicevars::kTrueLeadingMuonMom,slicecuts::kNoCut))});
    slice.insert({"true_leading_muon_costheta", ana::SpillMultiVar(PANDORAVAR(slicevars::kTrueLeadingMuonCostheta,slicecuts::kNoCut))});

    // Integer variables
    slice.insert({"event_type", ana::SpillMultiVar(PANDORAVAR(defs::kEventType,slicecuts::kNoCut))});
    slice.insert({"is_signal", ana::SpillMultiVar(PANDORAVAR(slicecuts::kIsSignal,slicecuts::kNoCut))});
    slice.insert({"Slice", ana::SpillMultiVar(PANDORAVAR(slicevars::kSliceIndex,slicecuts::kNoCut))});

    analysis.AddTree("slice", slice, true);

    // Add the header tree
    const SpillVar kIsMC = SIMPLESPILLVAR(hdr.ismc);
    std::vector<SpillVar> header_vars = {kIsMC};
    Tree hdrtree("header", {"is_mc"}, var00, header_vars, kNoSpillCut, true);

    // Add systematics

    // Get systematic weights and names
    std::vector<std::pair<int, int>> min_max;
    for(size_t i = 0; i < syst::genie_names.size(); ++i) min_max.emplace_back(-3, 3);

    // Initlaize storage for the trees
    std::vector<std::string> multisim_names;
    std::vector<std::vector<Var>> univsKnobs;
    std::vector<unsigned int> nuniverses;
    //std::vector<std::string> detsyst_names;
    //std::vector<std::pair<int,int>> min_max_det;
    // Load the knobs
    for(const auto& name: syst::flux_names) {
        // Hard code nuniv at 1000 for now since asking UniverseOracle breaks for 
        // piplus, piminus, kplus, kminus, and kzero.
        if (name != syst::flux_names[0]) continue;
        multisim_names.push_back(name);
        size_t nuniv = 1000;
        nuniverses.push_back(nuniv);
        univsKnobs.emplace_back();
        for(size_t i = 0; i < nuniv; ++i) 
        univsKnobs.back().push_back(GetUniverseWeight(name, i));
    }
    for(const auto& name: syst::xsec_multisim_names) {
        if (name != syst::flux_names[0]) continue;
        multisim_names.push_back(name);
        size_t nuniv = 100;
        nuniverses.push_back(nuniv);
        univsKnobs.emplace_back();
        for(size_t i = 0; i < nuniv; ++i) 
        univsKnobs.back().push_back(GetUniverseWeight(name, i));
    }
    //TODO: Add detsysts
    // for(size_t i = 0; i < detsysts.size(); ++i)
    // {
    //     min_max_det.emplace_back(-3,3);
    //     const ISyst* detsyst = detsysts[i];
    //     detsyst_names.push_back(detsyst->ShortName()+"_multisigma");
    // }

    // Add the trees
    //NSigmasTree detsysttree("detsystTree", detsyst_names, var00, detsysts, min_max_det, kNoSpillCut, kNoCut, kNoShift, true, true);
    //NSigmasTree nsigtree("multisigmaTree", syst::genie_names, var00, syst::genie_systs, min_max, kNoSpillCut, kNoCut, kNoShift, true, true);
    //NUniversesTree nunivtree("multisimTree", multisim_names, var00, univsKnobs, nuniverses, kNoSpillCut, kNoCut, kNoShift, true, true);
    NUniversesTree nunivtree("multisimTree", std::vector<std::string>{syst::flux_names[0]}, var00, univsKnobs, nuniverses, kNoSpillCut, kNoCut, kNoShift, true, true);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();

    //Write the trees
    std::string out_fname_root = out_fname + ".root"; //why tf is c++ so dumb?
    TFile fout(out_fname_root.c_str(), "UPDATE");
    TDirectory* dir = fout.mkdir("systs");
    //nsigtree.SaveTo(dir);
    nunivtree.SaveTo(dir);
    //detsysttree.SaveTo(dir);

    //Write the header tree to event/nominal folder
    TDirectory* dir2 = fout.GetDirectory("events/nominal");
    if (!dir2) {
        // Handle the error, directory does not exist
        std::cerr << "Directory 'events/nominal' does not exist." << std::endl;
        return;
    }
    hdrtree.SaveTo(dir2);
}