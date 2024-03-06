#pragma once

#include <map>
#include "TH1.h"

// Enum used to label event categories of interest for analysis plots
enum EventCategory {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = 0,

  // Signal events broken down by underlying reaction mode
  kNuMuCC0p0pi_CCQE = 1,
  kNuMuCC0p0pi_CCMEC = 2,
  kNuMuCC0p0pi_CCRES = 3,
  kNuMuCC0p0pi_Other = 4,
  
  kNuMuCC1p0pi_CCQE = 5,
  kNuMuCC1p0pi_CCMEC = 6,
  kNuMuCC1p0pi_CCRES = 7,
  kNuMuCC1p0pi_Other = 8,

  kNuMuCC2p0pi_CCQE = 9,
  kNuMuCC2p0pi_CCMEC = 10,
  kNuMuCC2p0pi_CCRES = 11,
  kNuMuCC2p0pi_Other = 12,

  // M = >2
  kNuMuCCMp0pi_CCQE = 13,
  kNuMuCCMp0pi_CCMEC = 14,
  kNuMuCCMp0pi_CCRES = 15,
  kNuMuCCMp0pi_Other = 16,

  // True numu CC event with at least one final-state pion above threshold
  kNuMuCCNpi = 17,

  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  kNuMuCCOther = 18,

  // True nue CC event
  kNuECC = 19,

  // True neutral current event for any neutrino flavor
  kNC = 20,

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 21,

  // All events that do not fall within any of the other categories (e.g.,
  // numubar CC)
  kOther = 22,

  nCategories
};

// Singleton class that helps manipulate EventCategory enum values
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

    inline const std::map< EventCategory, std::string >& label_map() const
      { return event_category_to_label_map_; }

    inline std::string label( EventCategory ec ) const
      { return event_category_to_label_map_.at( ec ); }

    inline int color_code( EventCategory ec ) const
      { return event_category_to_color_map_.at( ec ); }

    inline void set_mc_histogram_style( EventCategory ec, TH1* mc_hist ) const
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

  inline int get_number_categories() const {return nCategories;}

  private:

    EventCategoryInterpreter() {}

    std::map< EventCategory, std::string > event_category_to_label_map_ = {
      { kUnknown, "Unknown" },
      { kNuMuCC0p0pi_CCQE, "CCmu0p0pi (CCQE)" },
      { kNuMuCC0p0pi_CCMEC, "CCmu0p0pi (CCMEC)" },
      { kNuMuCC0p0pi_CCRES, "CCmu0p0pi (CCRES)" },
      { kNuMuCC0p0pi_Other, "CCmu0p0pi (Other)" },
      { kNuMuCC1p0pi_CCQE, "CCmu1p0pi (CCQE)" },
      { kNuMuCC1p0pi_CCMEC, "CCmu1p0pi (CCMEC)" },
      { kNuMuCC1p0pi_CCRES, "CCmu1p0pi (CCRES)" },
      { kNuMuCC1p0pi_Other, "CCmu1p0pi (Other)" },
      { kNuMuCC2p0pi_CCQE, "CCmu2p0pi (CCQE)" },
      { kNuMuCC2p0pi_CCMEC, "CCmu2p0pi (CCMEC)" },
      { kNuMuCC2p0pi_CCRES, "CCmu2p0pi (CCRES)" },
      { kNuMuCC2p0pi_Other, "CCmu2p0pi (Other)" },
      { kNuMuCCMp0pi_CCQE, "CCmuMp0pi (CCQE)" },
      { kNuMuCCMp0pi_CCMEC, "CCmuMp0pi (CCMEC)" },
      { kNuMuCCMp0pi_CCRES, "CCmuMp0pi (CCRES)" },
      { kNuMuCCMp0pi_Other, "CCmuMp0pi (Other)" },
      { kNuMuCCNpi, "#nu_{#mu} CCN#pi" },
      { kNuMuCCOther, "Other #nu_{#mu} CC" },
      { kNuECC, "#nu_{e} CC" },
      { kNC, "NC" },
      { kOOFV, "Out FV" },
      { kOther, "Other" }
    };

    std::map< EventCategory, int > event_category_to_color_map_ = {
      { kUnknown, kGray },
      { kNuMuCC0p0pi_CCQE, kBlue -2 },
      { kNuMuCC0p0pi_CCMEC, kBlue - 6 },
      { kNuMuCC0p0pi_CCRES, kBlue - 9 },
      { kNuMuCC0p0pi_Other, kBlue - 10 },
      { kNuMuCC1p0pi_CCQE, kOrange + 4 },
      { kNuMuCC1p0pi_CCMEC, kOrange + 5 },
      { kNuMuCC1p0pi_CCRES, kOrange + 6 },
      { kNuMuCC1p0pi_Other, kOrange + 7 },
      { kNuMuCC2p0pi_CCQE, kCyan - 3 },
      { kNuMuCC2p0pi_CCMEC, kCyan - 6 },
      { kNuMuCC2p0pi_CCRES, kCyan - 4 },
      { kNuMuCC2p0pi_Other, kCyan - 9 },
      { kNuMuCCMp0pi_CCQE, kGreen },
      { kNuMuCCMp0pi_CCMEC, kGreen + 1 },
      { kNuMuCCMp0pi_CCRES, kGreen + 2 },
      { kNuMuCCMp0pi_Other, kGreen + 3 },
      { kNuMuCCNpi, kAzure - 2 },
      { kNuMuCCOther, kAzure },
      { kNuECC, kViolet },
      { kNC, kOrange },
      { kOOFV, kRed + 3 },
      { kOther, kRed + 1 }
    };
};
