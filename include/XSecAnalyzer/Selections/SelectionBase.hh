#pragma once

// Standard library includes
#include <map>
#include <memory>
#include <set>
#include <string>

// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"

// Forward-declare required classes
class AnalysisEvent;

class SelectionBase {

public:

  using CategoryMap = std::map< int, std::pair< std::string, int > >;
  using TreeNameSet = std::set< std::string >;
  using SelectionConfig = std::map< std::string,
    std::pair< std::shared_ptr< TreeNameSet >,
    std::shared_ptr< CategoryMap > > >;

  SelectionBase( const std::string& sel_name );

  inline virtual ~SelectionBase() {};

  void apply_selection( bool is_mc, AnalysisEvent& event );

  inline const std::string& name() const { return selection_name_; }
  inline int passed_events() const { return num_passed_events_; }

  const CategoryMap& category_map() const;
  const TreeNameSet& input_tree_names() const;

  const FiducialVolume& get_FV() const;

  // Resets the internal state of the selection object to be ready for
  // processing a new event
  virtual void reset() {}

protected:

  void define_FV( double x_min, double x_max, double y_min,
    double y_max, double z_min, double z_max );

  // This group of virtual functions provides an interface
  // that needs to be defined in all concrete derived classes
  virtual bool is_selected( AnalysisEvent& ev ) = 0;
  virtual bool is_signal( AnalysisEvent& ev ) = 0;
  virtual std::string categorize_event( AnalysisEvent& ev ) = 0;
  virtual void compute_reco_observables( AnalysisEvent& ev ) = 0;
  virtual void compute_true_observables( AnalysisEvent& ev ) = 0;

private:

  // Load selection information from the global configuration file
  void load_selection_config();

  // Looks up the integer category code corresponding to a string
  // returned by categorize_event()
  int get_category_code( const std::string& categ_name ) const;

  // Helper function for category_map() and input_tree_names()
  const std::pair< std::shared_ptr< TreeNameSet >,
    std::shared_ptr< CategoryMap > >& load_info() const;

  std::string selection_name_;
  int num_passed_events_;

  int event_category_;

  bool selected_;
  bool mc_signal_;

  std::unique_ptr< FiducialVolume > fv_;

  static std::unique_ptr< SelectionConfig > selection_config_;

  static constexpr int DUMMY_CATEGORY_CODE = -9999;
  static constexpr int DEFAULT_CATEGORY_CODE = 0;
  static constexpr int DEFAULT_CATEGORY_COLOR = 920; // ROOT's kGray
  static constexpr char DEFAULT_CATEGORY_NAME[] = "Unknown";
};
