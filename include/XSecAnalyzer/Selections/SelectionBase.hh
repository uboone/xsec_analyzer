#pragma once

// Standard library includes
#include <string>
#include <type_traits>

// ROOT includes
#include "TTree.h"
#include "TVector3.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/AnalysisEvent.hh"
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/STVTools.hh"

class SelectionBase {

public:

  SelectionBase( const std::string& sel_name );

  void setup( TTree* out_tree, bool create_branches = true );
  void apply_selection( AnalysisEvent* event );
  void summary();

  virtual void final_tasks() {};

  inline bool is_event_mc_signal() { return mc_signal_; }
  inline bool is_event_selected() { return selected_; }

  inline const std::string& name() const { return selection_name_; }

  inline const std::map< int, std::pair< std::string, int > >&
    category_map() const { return categ_map_; }

protected:

  // Sets the branch address for output TTree variables managed by this
  // SelectionBase object
  template < typename T, typename S > void set_branch( T*& var, S var_name )
  {
    std::string var_name_str( var_name );
    std::string full_name = selection_name_ + '_' + var_name_str;
    std::string leaf_list = full_name;

    // The use of if constexpr here is a C++17 feature
    if constexpr ( std::is_same_v< T, bool > ) leaf_list += "/O";
    else if constexpr ( std::is_same_v< T, double > ) leaf_list += "/D";
    else if constexpr ( std::is_same_v< T, float > ) leaf_list += "/F";
    else if constexpr ( std::is_same_v< T, int > ) leaf_list += "/I";
    else if constexpr ( std::is_same_v< T, unsigned int > ) leaf_list += "/i";
    else leaf_list = "";

    // Branches for objects do not use a leaf list and use a
    // pointer-to-a-pointer to set the address
    if ( leaf_list.empty() ) {
      if ( need_to_create_branches_ ) {
        out_tree_->Branch( var_name_str.c_str(), &var );
      }
      else {
        out_tree_->SetBranchAddress( var_name_str.c_str(), &var );
      }
    }
    // Branches for simple types need to specify a leaf list upon creation
    // and use a regular pointer to set the address
    else {
      if ( need_to_create_branches_ ) {
        out_tree_->Branch( var_name_str.c_str(), var, leaf_list.c_str() );
      }
      else {
        out_tree_->SetBranchAddress( var_name_str.c_str(), var );
      }
    }
  }

  // Overloaded version that takes a MyPointer argument instead
  template < typename T, typename S >
    void set_branch( MyPointer<T>& u_ptr, S var_name )
  {
    T*& address = u_ptr.get_bare_ptr();
    this->set_branch( address, var_name );
  }

  // Overloaded version that takes an rvalue reference to a pointer argument
  // instead (see, e.g., https://stackoverflow.com/a/5465371/4081973)
  template < typename T, typename S >
    void set_branch( T*&& rval_ptr, S var_name )
  {
    T* address = rval_ptr;
    this->set_branch( address, var_name );
  }

  void setup_tree();
  void reset_base();

  inline int get_event_number() { return event_number_; }

  inline void define_true_FV( double XMin, double XMax, double YMin,
    double YMax, double ZMin, double ZMax )
  {
    fv_true_ = { XMin, XMax, YMin, YMax, ZMin, ZMax };
    set_fv_true_ = true;
  }

  inline FiducialVolume true_FV() {
    if ( !set_fv_true_ ) {
      std::cerr << "True Fiducial volume has not been defined"
        << " for selection:" << selection_name_ << '\n';
      throw;
    }
    return fv_true_;
  }

  inline void define_reco_FV( double XMin, double XMax, double YMin,
    double YMax, double ZMin, double ZMax )
  {
    fv_reco_ = { XMin, XMax, YMin, YMax, ZMin, ZMax };
    set_fv_reco_ = true;
  }

  inline FiducialVolume reco_FV() {
    if ( !set_fv_reco_ ) {
      std::cerr << "Reco Fiducial volume has not been defined"
        << " for selection:" << selection_name_ << '\n';
      throw;
    }
    return fv_reco_;
  }

  virtual bool selection( AnalysisEvent* event ) = 0;
  virtual int categorize_event( AnalysisEvent* event ) = 0;
  virtual void compute_reco_observables( AnalysisEvent* event ) = 0;
  virtual void compute_true_observables( AnalysisEvent* event ) = 0;
  virtual void define_output_branches() = 0;
  virtual bool define_signal( AnalysisEvent* event ) = 0;
  virtual void define_constants() = 0;
  virtual void define_category_map() = 0;
  virtual void reset() = 0;
  void define_additional_input_branches() {};

  TTree* out_tree_;
  bool need_to_create_branches_;

  STVTools stv_tools_;

protected:

  std::map< int, std::pair< std::string, int > > categ_map_;

private:

  std::string selection_name_;
  int num_passed_events_;

  int event_category_;

  bool selected_;
  bool mc_signal_;

  FiducialVolume fv_true_;
  FiducialVolume fv_reco_;
  bool set_fv_true_ = false;
  bool set_fv_reco_ = false;

  int event_number_;

};
