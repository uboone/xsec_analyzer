#pragma once

// Standard library includes
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// ROOT includes
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TStyle.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/HistUtils.hh"
#include "XSecAnalyzer/FilePropertiesManager.hh"

class Block: public TNamed {

  public:

    Block( const std::string& name, const std::string& title,
      const int i ) : fName( name ), fTitle( title ), binType( i ) {}
    virtual ~Block() = default;

    virtual int GetNBinsX() = 0;
    virtual int GetNBinsY(const int idx) = 0;
    virtual std::string GetBinDef(const int idx)  = 0;
    virtual std::string GetBinDef(const int x, const int y)  = 0;
    virtual int GetBinType() = 0;
    virtual Double_t GetBinXLow(const int i) = 0;
    virtual Double_t GetBinXHigh(const int i) = 0;

    virtual std::string GetXName() = 0;
    virtual std::string GetXNameUnit() = 0;
    virtual std::string GetXTitle() = 0;
    virtual std::string GetXTitleUnit() = 0;
    virtual std::string GetXTexTitle() = 0;
    virtual std::string GetXTexTitleUnit() = 0;
    virtual std::string GetYName() = 0;
    virtual std::string GetYNameUnit()  = 0;
    virtual std::string GetYTitle() = 0;
    virtual std::string GetYTitleUnit() = 0;
    virtual std::string GetYTexTitle() = 0;
    virtual std::string GetYTexTitleUnit() = 0;
    virtual std::string GetSelection() = 0;

    virtual std::vector<double> GetVector() = 0;
    virtual std::vector<double> GetVector( const int i ) = 0;
    virtual bool Is1D() = 0;
    virtual bool Is2D() = 0;

  protected:

    std::string fName;
    std::string fTitle;
    int binType;
    std::vector<std::string> binDef;
    std::string xName;
    std::string xNameUnit;
    std::string xTitle;
    std::string xTitleUnit;
    std::string xTexTitle;
    std::string xTexTitleUnit;
    std::string yName;
    std::string yNameUnit;
    std::string yTitle;
    std::string yTitleUnit;
    std::string yTexTitle;
    std::string yTexTitleUnit;
};


class Block1D: public Block {
  public:
    Block1D( const std::string& name, const std::string& title,
      const std::vector<double>& block1d, const std::string& selection,
      const int i ) : Block( name, title, i ), fName( name ), fTitle( title ),
      fblock( block1d ), fselection( selection ), binType( i )
      { this->Init(); }

    Block1D( const std::string& name, const std::string& title,
      const std::string& textitle, const std::vector<double>& block1d,
      const std::string& selection, const int i ) : Block( name, title, i ),
      fName( name ), fTitle( title ), fTexTitle( textitle ), fblock( block1d ),
      fselection( selection ), binType( i )
      { this->Init(); }

    virtual ~Block1D() = default;
    int GetNBinsX(){ return binDef.size(); }
    int GetNBinsY( const int idx ) { return 0; }
    std::string GetBinDef( const int idx ){ return binDef[idx]; }
    std::string GetBinDef( const int x, const int y ){ return GetBinDef( x ); }
    int GetBinType() { return int(binType); }

    Double_t GetBinXLow( const int i ) { return -9999999; }
    Double_t GetBinXHigh( const int i ) { return -9999999; }

    std::string GetXName() { return xName; }
    std::string GetXNameUnit() { return xNameUnit; }
    std::string GetXTitle() { return xTitle; }
    std::string GetXTitleUnit() { return xTitleUnit; }
    std::string GetXTexTitle() { return xTexTitle; }
    std::string GetXTexTitleUnit() { return xTexTitleUnit; }
    std::string GetYName() { return yName; }
    std::string GetYNameUnit() { return yNameUnit; }
    std::string GetYTitle() { return yTitle; }
    std::string GetYTitleUnit() { return yTitleUnit; }
    std::string GetYTexTitle() { return yTexTitle; }
    std::string GetYTexTitleUnit() { return yTexTitleUnit; }

    std::string GetSelection() { return fselection; }

    std::vector<double> GetVector() { return fblock; }
    std::vector<double> GetVector( const int /*i*/ ) { return GetVector(); }
    bool Is1D() { return true; }
    bool Is2D() { return false; }

  protected:

    std::string fName;
    std::string fTitle;
    std::string fTexTitle;
    std::vector<double> fblock;
    std::string fselection;
    int binType;

    std::string xName;
    std::string xNameUnit;
    std::string xTitle;
    std::string xTitleUnit;
    std::string xTexTitle;
    std::string xTexTitleUnit;

    std::string yName;
    std::string yNameUnit;
    std::string yTitle;
    std::string yTitleUnit;
    std::string yTexTitle;
    std::string yTexTitleUnit;

    std::vector< std::string > binDef;

  private:
    void Init();
    void SetTitle( const std::string& title );
    void SetTexTitle( const std::string& title );
    void SetName( const std::string& name );
};

class Block2D : public Block {
  public:
    Block2D( const std::string& name, const std::string& title,
      const std::map< double, std::vector<double> >& block2d,
      const std::string& selection, int i ) : Block( name, title, i ),
      fName( name ), fTitle( title ), fblock( block2d ),
      fselection( selection ), binType( i )
      { this->Init(); }

    Block2D(const std::string& name, const std::string& title,
     const std::string& textitle,
     const std::map< double, std::vector<double> >& block2d,
     const std::string& selection, int i ) : Block( name, title, i ),
     fName( name ), fTitle( title ), fTexTitle( textitle ),
     fblock( block2d ), fselection( selection ), binType( i )
     { this->Init(); }

    virtual ~Block2D() = default;
    int GetNBinsX(){ return fblock_vv_.size(); }
    int GetNBinsY( const int idx ) { return (fblock_vv_[idx].size() - 1); }
    std::string GetBinDef( const int x, const int y ) {
      int index = 0;
      for( auto i = 0; i < x; i++ ) {
        index += (fblock_vv_[i].size() - 1);
      }
      index += y;
      return this->GetBinDef( index );
    }
    std::string GetBinDef( const int idx ){ return binDef[idx]; }
    int GetBinType(){ return int(binType); }

    Double_t GetBinXLow( const int i ) { return xbin[i]; }
    Double_t GetBinXHigh( const int i ) { return xbin[i+1]; }

    std::string GetXName() { return xName; }
    std::string GetXNameUnit() { return xNameUnit; }
    std::string GetXTitle() { return xTitle; }
    std::string GetXTitleUnit() { return xTitleUnit; }
    std::string GetXTexTitle() { return xTexTitle; }
    std::string GetXTexTitleUnit() { return xTexTitleUnit; }
    std::string GetYName() { return yName; }
    std::string GetYNameUnit() { return yNameUnit; }
    std::string GetYTitle() { return yTitle; }
    std::string GetYTitleUnit() { return yTitleUnit; }
    std::string GetYTexTitle() { return yTexTitle; }
    std::string GetYTexTitleUnit() { return yTexTitleUnit; }
    std::string GetSelection() { return fselection; }

    std::vector<double> GetVector() { return GetVector( 0u ); }
    std::vector<double> GetVector( const int i ) { return fblock_vv_[i]; }
    std::map<double, std::vector<double>> GetMap() { return fblock; }
    bool Is1D() { return false; }
    bool Is2D() { return true; }

  private:

    std::string fName;
    std::string fTitle;
    std::string fTexTitle;
    std::map< double, std::vector<double> >  fblock;
    std::vector< std::vector<double> > fblock_vv_;
    std::vector<double> xbin;
    std::string fselection;
    int binType;

    std::string xName;
    std::string xNameUnit;
    std::string xTitle;
    std::string xTitleUnit;
    std::string xTexTitle;
    std::string xTexTitleUnit;

    std::string yName;
    std::string yNameUnit;
    std::string yTitle;
    std::string yTitleUnit;
    std::string yTexTitle;
    std::string yTexTitleUnit;

    std::vector< std::string > binDef;

  private:
    void Init();
    void SetTitle(const std::string& title);
    void SetTexTitle(const std::string& title);
    void SetName( const std::string& name );
};


class BlockTrueReco {
  public:
    BlockTrueReco( Block* btrue, Block* breco) : block_true_( btrue ),\
      block_reco_( breco ) {}
    virtual ~BlockTrueReco() = default;
    Block* block_true_;
    Block* block_reco_;
};


class MakeConfig {

  public:

    MakeConfig();
    virtual ~MakeConfig();

    void BinScheme();
    void Print();
    void ResPlots();

  private:

    std::vector< BlockTrueReco > vect_block;
    std::string TREE = "stv_tree";
    std::string BIN_CONFIG;
    std::string SELECTION;
    std::string DIRECTORY;
    std::string CATEGORY;
    std::vector<int> background_index;

    // The anticipated POT to use when scaling the MC prediction in the
    // expected reco events plot. This will help ensure that all choices of
    // reco binning are informed by the expected statistical uncertainties when
    // the full dataset is analyzed.
    //
    // Full dataset for Runs 1-3
    static constexpr double EXPECTED_POT = 6.9314e21;//4.133e19;//1.28e21;

    // Number of true bins to use when plotting true distributions in a given
    // reco bin
    static constexpr int DEFAULT_TRUE_BINS = 100;

    // ROOT integer code for Arial font
    static constexpr int FONT_STYLE = 62; // Arial

    // When bins have zero content, set them to this very small value so that
    // the colz style will still paint them
    static constexpr double REALLY_SMALL = 1e-11;

    // Dummy counter used to ensure that each histogram generated by this
    // function has a unique name to use with TTree::Draw()
    static int hist_count;

    std::set< int > RUNS;

    void make_res_plots( const std::string& branchexpr,
        const std::string& variable_title, const std::string& selection,
        const std::set<int>& runs, std::vector<double> bin_low_edges,
        bool show_bin_plots = true,
        bool show_smear_numbers = false,
        int num_true_bins = DEFAULT_TRUE_BINS,
        const std::string& mc_branchexpr = "",
        const std::string& signal_cuts = "mc_is_signal",
        const std::string& mc_event_weight = DEFAULT_MC_EVENT_WEIGHT );


    // Overloaded version that uses a fixed number of equal-width bins
    void make_res_plots( const std::string& branchexpr,
        const std::string& variable_title, const std::string& selection,
        const std::set<int>& runs,
        double xmin, double xmax, int Nbins,
        bool show_bin_plots = true,
        bool show_smear_numbers = false,
        int num_true_bins = DEFAULT_TRUE_BINS,
        const std::string& mc_branchexpr = "",
        const std::string& signal_cuts = "mc_is_signal",
        const std::string& mc_event_weight = DEFAULT_MC_EVENT_WEIGHT );

    void make_res_plots( std::istream& in_file,
        const std::set<int>& runs,
        const std::string& universe_branch_name = "TunedCentralValue_UBGenie",
        size_t universe_index = 0u,
        bool show_smear_numbers = false );
};
