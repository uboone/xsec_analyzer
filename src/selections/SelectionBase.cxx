// Standard library includes
#include <iostream>

// XSecAnalyzer includes
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/Selections/SelectionBase.hh"

SelectionBase::SelectionBase( const std::string& sel_name ) {

  selection_name_ = sel_name;
  num_passed_events_ = 0;

  event_number_ = 0;

  fv_true_ = { BOGUS, BOGUS, BOGUS, BOGUS, BOGUS, BOGUS };
  fv_reco_ = { BOGUS, BOGUS, BOGUS, BOGUS, BOGUS, BOGUS };

}

void SelectionBase::setup( TTree* out_tree, bool create_branches ) {

  out_tree_ = out_tree;
  need_to_create_branches_ = create_branches;
  this->setup_tree();
  this->define_category_map();
  this->define_constants();

}

void SelectionBase::apply_selection( AnalysisEvent* event ) {
  this->reset_base();
  this->reset();

  mc_signal_ = this->define_signal( event );
  selected_ = this->selection( event );
  event_category_ = this->categorize_event( event );

  this->compute_reco_observables( event );

  // Note that event->is_mc_ is set in CategorizeEvent() above
  if ( event->is_mc_ ) {
    this->compute_true_observables( event );
  }

  if ( selected_ ) {
    ++num_passed_events_;
  }

  ++event_number_;
}

void SelectionBase::summary() {
  std::cout << selection_name_ << " has " << num_passed_events_
    << " events which passed\n";
}

void SelectionBase::setup_tree() {

  this->set_branch( &selected_, "Selected" );
  this->set_branch( &mc_signal_, "MC_Signal" );
  this->set_branch( &event_category_, "EventCategory" );

  this->define_additional_input_branches();
  this->define_output_branches();

}

void SelectionBase::reset_base() {
  selected_ = false;
  mc_signal_ = false;
  event_category_ = BOGUS_INDEX;
}
