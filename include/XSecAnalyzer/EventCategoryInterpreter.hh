#pragma once

#include "TH1.h"

#include "XSecAnalyzer/Selections/EventCategoriesNp1pi.hh"
class EventCategoryInterpreter {

  public:

    // This is a singleton class, so we'll delete the copy constructor
    // the move constructor, and the assignment operators
    EventCategoryInterpreter( const EventCategoryInterpreter& ) = delete;
    EventCategoryInterpreter( EventCategoryInterpreter&& ) = delete;
    EventCategoryInterpreter& operator=( const EventCategoryInterpreter& )
      = delete;
    EventCategoryInterpreter& operator=( EventCategoryInterpreter&& )
      = delete;

    // Get a const reference to the singleton instance of the
    // EventCategoryInterpreter
    inline static const EventCategoryInterpreter& Instance() {

      // Create the EventCategoryInterpreter object using a static variable.
      // This ensures that the singleton instance is only created once.
      static std::unique_ptr<EventCategoryInterpreter>
        the_instance( new EventCategoryInterpreter() );

      // Return a reference to the singleton instance
      return *the_instance;
    }

    // Function to get the label map
    inline const std::map<EventCategoryNp1pi, std::string>& label_map() const {
        return event_category_to_label_map_;
    }

    inline std::string label( EventCategoryNp1pi ec ) const
      { return event_category_to_label_map_.at( ec ); }

    inline int color_code( EventCategoryNp1pi ec ) const
      { return event_category_to_color_map_.at( ec ); }

    inline void set_mc_histogram_style( EventCategoryNp1pi ec, TH1* mc_hist ) const
    {
      int color = color_code( ec );
      mc_hist->SetFillColor( color );
      mc_hist->SetLineColor( color );
      mc_hist->SetStats( false );
    }

    inline void set_ext_histogram_style( TH1* ext_hist ) const {
      ext_hist->SetFillColor( 28 );
      ext_hist->SetLineColor( 28 );
      ext_hist->SetLineWidth( 2 );
      ext_hist->SetFillStyle( 3005 );
      ext_hist->SetStats( false );
    }

    inline void set_bnb_data_histogram_style( TH1* bnb_hist ) const {

      bnb_hist->SetLineColor( kBlack );
      bnb_hist->SetLineWidth( 3 );
      bnb_hist->SetMarkerStyle( kFullCircle );
      bnb_hist->SetMarkerSize( 0.8 );
      bnb_hist->SetStats( false );

      bnb_hist->GetXaxis()->SetTitleOffset( 0.0 );
      bnb_hist->GetXaxis()->SetTitleSize( 0.0 );
      bnb_hist->GetYaxis()->SetTitleSize( 0.05 );
      bnb_hist->GetYaxis()->CenterTitle( true );
      bnb_hist->GetXaxis()->SetLabelSize( 0.0 );

      // This prevents the first y-axis label label (0) to be clipped by the
      // ratio plot
      bnb_hist->SetMinimum( 1e-3 );
    }

    inline void set_stat_err_histogram_style( TH1* stat_err_hist ) const {
      stat_err_hist->SetFillColor( kBlack );
      stat_err_hist->SetLineColor( kBlack );
      stat_err_hist->SetLineWidth( 2 );
      stat_err_hist->SetFillStyle( 3004 );
    }


  private:

    // Private constructor to prevent instantiation
    EventCategoryInterpreter() {
        for (const auto& entry : CC1muNp1pi_MAP) {
            event_category_to_label_map_[static_cast<EventCategoryNp1pi>(entry.first)] = entry.second.first;
            event_category_to_color_map_[static_cast<EventCategoryNp1pi>(entry.first)] = entry.second.second;

        }
    }
    std::map<EventCategoryNp1pi, std::string> event_category_to_label_map_;
    std::map<EventCategoryNp1pi, int> event_category_to_color_map_;

};
