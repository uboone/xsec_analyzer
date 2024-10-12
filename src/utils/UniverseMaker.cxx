// XSecAnalyzer includes
#include "XSecAnalyzer/TreeUtils.hh"
#include "XSecAnalyzer/UniverseMaker.hh"

// Define this static member of the Universe class
size_t Universe::num_categories_;

UniverseMaker::UniverseMaker( const std::string& config_file_name ) {
  std::ifstream in_file( config_file_name );
  this->init( in_file );
}

UniverseMaker::UniverseMaker( std::istream& config_stream ) {
  this->init( config_stream );
}

void UniverseMaker::init( std::istream& in_file ) {

  // Load the root TDirectoryFile name to use when writing the universes to an
  // output ROOT file
  in_file >> output_directory_name_;

  // Load the TTree name to use when processing ntuple input files
  std::string ttree_name;
  in_file >> ttree_name;

  // Initialize the owned input TChain with the configured TTree name
  input_chain_.SetName( ttree_name.c_str() );

  // Load the name of the selection whose event category definitions
  // will be used for populating the category histogram
  std::string sel_categ_name;
  in_file >> sel_categ_name;

  // Instantiate the requested selection and store it in this object for later
  // use
  SelectionFactory sel_fact;
  //SelectionBase* temp_sel = sel_fact.CreateSelection( sel_categ_name );
  //sel_for_categories_.reset( temp_sel );
  //FIXME: using normal pointer to avoid invalid pointer error
  sel_for_categories_ = sel_fact.CreateSelection( sel_categ_name);

  // Load the true bin definitions
  size_t num_true_bins;
  in_file >> num_true_bins;

  for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
    TrueBin temp_bin;
    in_file >> temp_bin;

    /*
    // DEBUG
    std::cout << "tb = " << tb << '\n';
    std::cout << temp_bin << '\n';
    */

    true_bins_.push_back( temp_bin );
  }

  // Load the reco bin definitions
  size_t num_reco_bins;
  in_file >> num_reco_bins;
  for ( size_t rb = 0u; rb < num_reco_bins; ++rb ) {
    RecoBin temp_bin;
    in_file >> temp_bin;

    /*
    // DEBUG
    std::cout << "rb = " << rb << '\n';
    std::cout << temp_bin << '\n';
    */

    reco_bins_.push_back( temp_bin );
  }

}

void UniverseMaker::add_input_file( const std::string& input_file_name )
{
  // Check to make sure that the input file contains the expected ntuple
  TFile temp_file( input_file_name.c_str(), "read" );

  // Temporary storage
  TTree* temp_tree;

  std::string tree_name = input_chain_.GetName();
  temp_file.GetObject( tree_name.c_str(), temp_tree );
  if ( !temp_tree ) throw std::runtime_error( "Missing ntuple TTree "
    + tree_name + " in the input ntuple file " + input_file_name );

  // If we've made it here, then the input file has passed all of the checks.
  // Add it to the input TChain.
  input_chain_.AddFile( input_file_name.c_str() );
}

void UniverseMaker::prepare_formulas() {

  // Remove any pre-existing TTreeFormula objects from the owned vectors
  true_bin_formulas_.clear();
  reco_bin_formulas_.clear();
  category_formulas_.clear();

  // Create one TTreeFormula for each true bin definition
  for ( size_t tb = 0u; tb < true_bins_.size(); ++tb ) {
    const auto& bin_def = true_bins_.at( tb );
    std::string formula_name = "true_formula_" + std::to_string( tb );

    auto tbf = std::make_unique< TTreeFormula >( formula_name.c_str(),
      bin_def.signal_cuts_.c_str(), &input_chain_ );

    tbf->SetQuickLoad( true );

    true_bin_formulas_.emplace_back( std::move(tbf) );
  }

  // Create one TTreeFormula for each reco bin definition
  for ( size_t rb = 0u; rb < reco_bins_.size(); ++rb ) {
    const auto& bin_def = reco_bins_.at( rb );
    std::string formula_name = "reco_formula_" + std::to_string( rb );

    auto rbf = std::make_unique< TTreeFormula >( formula_name.c_str(),
      bin_def.selection_cuts_.c_str(), &input_chain_ );

    rbf->SetQuickLoad( true );

    reco_bin_formulas_.emplace_back( std::move(rbf) );
  }

  // Create one TTreeFormula for each true event category
  const auto& category_map = sel_for_categories_->category_map();
  Universe::set_num_categories( category_map.size() );
  for ( const auto& category_pair : category_map ) {

    int cur_category = static_cast< int >( category_pair.first );
    std::string str_category = std::to_string( cur_category );

    std::string category_formula_name = "category_formula_" + str_category;

    std::string category_cuts = sel_for_categories_->name()
      + "_EventCategory == " + str_category;

    auto cbf = std::make_unique< TTreeFormula >(
      category_formula_name.c_str(), category_cuts.c_str(), &input_chain_ );

    cbf->SetQuickLoad( true );

    category_formulas_.emplace_back( std::move(cbf) );
  }

}

void UniverseMaker::build_universes(
  const std::vector<std::string>& universe_branch_names )
{
  return this->build_universes( &universe_branch_names );
}

void UniverseMaker::build_universes(
  const std::vector<std::string>* universe_branch_names )
{
  int num_input_files = input_chain_.GetListOfFiles()->GetEntries();
  if ( num_input_files < 1 ) {
    std::cout << "ERROR: The UniverseMaker object has not been"
      " initialized with any input files yet.\n";
    return;
  }

  WeightHandler wh;
  wh.set_branch_addresses( input_chain_, universe_branch_names );

  // Make sure that we always have branches set up for the CV correction
  // weights, i.e., the spline and tune weights. Don't throw an exception if
  // these are missing in the input TTree (we could be working with real data)
  wh.add_branch( input_chain_, SPLINE_WEIGHT_NAME, false );
  wh.add_branch( input_chain_, TUNE_WEIGHT_NAME, false );

  this->prepare_formulas();

  // Set up storage for the "is_mc" boolean flag branch. If we're not working
  // with MC events, then we shouldn't do anything with the true bin counts.
  bool is_mc;
  input_chain_.SetBranchAddress( "is_mc", &is_mc );

  // Get the first TChain entry so that we can know the number of universes
  // used in each vector of weights
  input_chain_.GetEntry( 0 );

  // Now prepare the vectors of Universe objects with the correct sizes
  this->prepare_universes( wh );

  int treenumber = 0;
  for ( long long entry = 0; entry < input_chain_.GetEntries(); ++entry ) {
    // Load the TTree for the current TChain entry
    input_chain_.LoadTree( entry );

    // If the current entry is in a new TTree, then have all of the
    // TTreeFormula objects make the necessary updates
    if ( treenumber != input_chain_.GetTreeNumber() ) {
      treenumber = input_chain_.GetTreeNumber();
      for ( auto& tbf : true_bin_formulas_ ) tbf->Notify();
      for ( auto& rbf : reco_bin_formulas_ ) rbf->Notify();
      for ( auto& cbf : category_formulas_ ) cbf->Notify();
    }

    // Find the reco bin(s) that should be filled for the current event
    std::vector< FormulaMatch > matched_reco_bins;
    for ( size_t rb = 0u; rb < reco_bin_formulas_.size(); ++rb ) {
      auto& rbf = reco_bin_formulas_.at( rb );
      int num_formula_elements = rbf->GetNdata();
      for ( int el = 0; el < num_formula_elements; ++el ) {
        double formula_wgt = rbf->EvalInstance( el );
        if ( formula_wgt ) matched_reco_bins.emplace_back( rb, formula_wgt );
      }
    }

    // Find the EventCategory label(s) that apply to the current event
    std::vector< FormulaMatch > matched_category_indices;
    for ( size_t c = 0u; c < category_formulas_.size(); ++c ) {
      auto& cbf = category_formulas_.at( c );
      int num_formula_elements = cbf->GetNdata();
      for ( int el = 0; el < num_formula_elements; ++el ) {
        double formula_wgt = cbf->EvalInstance( el );
        if ( formula_wgt ) {
          matched_category_indices.emplace_back( c, formula_wgt );
        }
      }
    }

    input_chain_.GetEntry( entry );

    std::vector< FormulaMatch > matched_true_bins;
    double spline_weight = 0.;
    double tune_weight = 0.;

    // If we're working with an MC sample, then find the true bin(s)
    // that should be filled for the current event
    if ( is_mc ) {
      for ( size_t tb = 0u; tb < true_bin_formulas_.size(); ++tb ) {
        auto& tbf = true_bin_formulas_.at( tb );
        int num_formula_elements = tbf->GetNdata();
        for ( int el = 0; el < num_formula_elements; ++el ) {
          double formula_wgt = tbf->EvalInstance( el );
          if ( formula_wgt ) matched_true_bins.emplace_back( tb, formula_wgt );
        }
      } // true bins

      // If we have event weights in the map at all, then get the current
      // event's CV correction weights here for potentially frequent re-use
      // below
      auto& wm = wh.weight_map();
      if ( wm.size() > 0u ) {
        spline_weight = wm.at( SPLINE_WEIGHT_NAME )->front();
        tune_weight = wm.at( TUNE_WEIGHT_NAME )->front();
      }
    } // MC event

    for ( const auto& pair : wh.weight_map() ) {
      const std::string& wgt_name = pair.first;
      const auto& wgt_vec = pair.second;

      auto& u_vec = universes_.at( wgt_name );

      for ( size_t u = 0u; u < wgt_vec->size(); ++u ) {

        // No need to use the slightly slower "at" here since we're directly
        // looping over the weight vector
        double w = wgt_vec->operator[]( u );

        // Multiply by any needed CV correction weights
        apply_cv_correction_weights( wgt_name, w, spline_weight, tune_weight );

        // Deal with NaNs, etc. to make a "safe weight" in all cases
        double safe_wgt = safe_weight( w );

        // Get the universe object that should be filled with the processed
        // event weight
        auto& universe = u_vec.at( u );

        for ( const auto& tb : matched_true_bins ) {
          // TODO: consider including the TTreeFormula weight(s) in the check
          // applied via safe_weight() above
          universe.hist_true_->Fill( tb.bin_index_, tb.weight_ * safe_wgt );
          for ( const auto& rb : matched_reco_bins ) {
            universe.hist_2d_->Fill( tb.bin_index_, rb.bin_index_,
              tb.weight_ * rb.weight_ * safe_wgt );
          } // reco bins

          for ( const auto& other_tb : matched_true_bins ) {
            universe.hist_true2d_->Fill( tb.bin_index_, other_tb.bin_index_,
              tb.weight_ * other_tb.weight_ * safe_wgt );
          } // true bins

        } // true bins

        for ( const auto& rb : matched_reco_bins ) {
          universe.hist_reco_->Fill( rb.bin_index_, rb.weight_ * safe_wgt );

          for ( const auto& c : matched_category_indices ) {
            universe.hist_categ_->Fill( c.bin_index_, rb.bin_index_,
              c.weight_ * rb.weight_ * safe_wgt );
          }

          for ( const auto& other_rb : matched_reco_bins ) {
            universe.hist_reco2d_->Fill( rb.bin_index_, other_rb.bin_index_,
              rb.weight_ * other_rb.weight_ * safe_wgt );
          }
        } // reco bins
      } // universes
    } // weight names

    // Fill the unweighted histograms now that we're done with the
    // weighted ones. Note that "unweighted" in this context applies to
    // the universe event weights, but that any implicit weights from
    // the TTreeFormula evaluations will still be applied.
    auto& univ = universes_.at( UNWEIGHTED_NAME ).front();
    for ( const auto& tb : matched_true_bins ) {
      univ.hist_true_->Fill( tb.bin_index_, tb.weight_ );
      for ( const auto& rb : matched_reco_bins ) {
        univ.hist_2d_->Fill( tb.bin_index_, rb.bin_index_,
          tb.weight_ * rb.weight_ );
      } // reco bins

      for ( const auto& other_tb : matched_true_bins ) {
        univ.hist_true2d_->Fill( tb.bin_index_, other_tb.bin_index_,
          tb.weight_ * other_tb.weight_ );
      } // true bins

    } // true bins

    for ( const auto& rb : matched_reco_bins ) {

      univ.hist_reco_->Fill( rb.bin_index_, rb.weight_ );

      for ( const auto& c : matched_category_indices ) {
        univ.hist_categ_->Fill( c.bin_index_, rb.bin_index_,
          c.weight_ * rb.weight_ );
      }

      for ( const auto& other_rb : matched_reco_bins ) {
        univ.hist_reco2d_->Fill( rb.bin_index_, other_rb.bin_index_,
          rb.weight_ * other_rb.weight_ );
      }

    } // reco bins

  } // TChain entries

  input_chain_.ResetBranchAddresses();
}

void UniverseMaker::prepare_universes( const WeightHandler& wh ) {

  size_t num_true_bins = true_bins_.size();
  size_t num_reco_bins = reco_bins_.size();

  for ( const auto& pair : wh.weight_map() ) {
    const std::string& weight_name = pair.first;
    size_t num_universes = pair.second->size();

    std::vector< Universe > u_vec;

    for ( size_t u = 0u; u < num_universes; ++u ) {
      u_vec.emplace_back( weight_name, u, num_true_bins, num_reco_bins );
    }

    universes_[ weight_name ] = std::move( u_vec );
  }

  // Add the special "unweighted" universe unconditionally
  std::vector< Universe > temp_uvec;
  temp_uvec.emplace_back( UNWEIGHTED_NAME, 0, num_true_bins, num_reco_bins );
  universes_[ UNWEIGHTED_NAME ] = std::move( temp_uvec );

}

void UniverseMaker::save_histograms(
  const std::string& output_file_name,
  const std::string& subdirectory_name,
  bool update_file )
{
  // Decide whether to overwrite the output file or simply update the contents.
  // This difference is only important if the output file already exists before
  // this function is called.
  std::string tfile_option( "recreate" );
  if ( update_file ) {
    tfile_option = "update";
  }
  TFile out_file( output_file_name.c_str(), tfile_option.c_str() );

  // Navigate to the subdirectory within the output ROOT file where the
  // universe histograms will be saved. Create new TDirectoryFile objects as
  // needed.
  TDirectoryFile* root_tdir = nullptr;
  TDirectoryFile* sub_tdir = nullptr;

  out_file.GetObject( output_directory_name_.c_str(), root_tdir );
  if ( !root_tdir ) {
    // TODO: add error handling for a forward slash in the root TDirectoryFile
    // name
    root_tdir = new TDirectoryFile( output_directory_name_.c_str(),
      "universes", "", &out_file );
  }

  // Save the configuration settings for this class to the main
  // TDirectoryFile before moving on to the appropriate subdirectory. If
  // these settings have already been saved, then double-check that they
  // match the current configuration. In the event of a mismatch, throw an
  // exception to avoid data corruption.
  std::string tree_name, true_bin_spec, reco_bin_spec;
  tree_name = this->input_chain().GetName();

  std::ostringstream oss_true, oss_reco;

  for ( const auto& tbin : true_bins_ ) {
    oss_true << tbin << '\n';
  }

  for ( const auto& rbin : reco_bins_ ) {
    oss_reco << rbin << '\n';
  }

  true_bin_spec = oss_true.str();
  reco_bin_spec = oss_reco.str();

  std::string* saved_tree_name = nullptr;
  std::string* saved_tb_spec = nullptr;
  std::string* saved_rb_spec = nullptr;
  std::string* saved_sel_for_categ_name = nullptr;
  root_tdir->GetObject( "ntuple_name", saved_tree_name );
  root_tdir->GetObject( TRUE_BIN_SPEC_NAME.c_str(), saved_tb_spec );
  root_tdir->GetObject( RECO_BIN_SPEC_NAME.c_str(), saved_rb_spec );
  root_tdir->GetObject( "sel_for_categ", saved_sel_for_categ_name );

  if ( saved_tree_name ) {
    if ( tree_name != *saved_tree_name ) {
      throw std::runtime_error( "Tree name mismatch: " + tree_name
        + " vs. " + *saved_tree_name );
    }
  }
  else {
    root_tdir->WriteObject( &tree_name, "ntuple_name" );
  }

  if ( saved_tb_spec ) {
    if ( true_bin_spec != *saved_tb_spec ) {
      throw std::runtime_error( "Inconsistent true bin specification!" );
    }
  }
  else {
    root_tdir->WriteObject( &true_bin_spec, TRUE_BIN_SPEC_NAME.c_str() );
  }

  if ( saved_rb_spec ) {
    if ( reco_bin_spec != *saved_rb_spec ) {
      throw std::runtime_error( "Inconsistent reco bin specification!" );
    }
  }
  else {
    root_tdir->WriteObject( &reco_bin_spec, RECO_BIN_SPEC_NAME.c_str() );
  }

  const std::string& sel_for_categ_name = sel_for_categories_->name();
  if ( saved_sel_for_categ_name ) {
    if ( sel_for_categ_name != *saved_sel_for_categ_name ) {
      throw std::runtime_error( "Inconsistent selections configured for event"
        " categorization" );
    }
  }
  else {
    root_tdir->WriteObject( &sel_for_categ_name, "sel_for_categ" );
  }

  std::string subdir_name = ntuple_subfolder_from_file_name(
    subdirectory_name );

  root_tdir->GetObject( subdir_name.c_str(), sub_tdir );
  if ( !sub_tdir ) {
    sub_tdir = new TDirectoryFile( subdir_name.c_str(), "universes",
      "", root_tdir );
  }

  // Now we've found (or created) the TDirectoryFile where the output
  // will be saved. Ensure that it is the active file here before writing
  // out the histograms.
  sub_tdir->cd();

  for ( auto& pair : universes_ ) {
    auto& u_vec = pair.second;
    for ( auto& univ : u_vec ) {
      // Always save the reco histograms
      univ.hist_reco_->Write();
      univ.hist_reco2d_->Write();

      // Save the others if the true histogram was filled at least once
      // (used to infer that we have MC truth information)
      if ( univ.hist_true_->GetEntries() > 0. ) {
        univ.hist_true_->Write();
        univ.hist_2d_->Write();
        univ.hist_categ_->Write();
        univ.hist_true2d_->Write();
      }
    } // universes
  } // weight names
}
