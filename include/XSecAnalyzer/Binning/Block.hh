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

class Block: public TNamed {

  public:

    Block( const std::string& name, const std::string& title,
      const int i ) : fName( name ), fTitle( title ), binType( i ) {}
    virtual ~Block() = default;

    virtual int GetNBinsX() = 0;
    virtual int GetNBinsY(const int idx) = 0;
    virtual int GetNBinsZ(const int idx, const int idy) = 0;
    virtual std::string GetBinDef(const int idx)  = 0;
    virtual std::string GetBinDef(const int x, const int y) = 0;
    virtual std::string GetBinDef(const int x, const int y, const int z) = 0;
    virtual int GetBinType() = 0;
    virtual Double_t GetBinXLow(const int i) = 0;
    virtual Double_t GetBinXHigh(const int i) = 0;
    virtual Double_t GetBinYLow(const int i, const int j) = 0;
    virtual Double_t GetBinYHigh(const int i, const int j) = 0;

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
    virtual std::string GetZName() = 0;
    virtual std::string GetZNameUnit()  = 0;
    virtual std::string GetZTitle() = 0;
    virtual std::string GetZTitleUnit() = 0;
    virtual std::string GetZTexTitle() = 0;
    virtual std::string GetZTexTitleUnit() = 0;
    virtual std::string GetSelection() = 0;

    virtual std::vector<double> GetVector() = 0;
    virtual std::vector<double> GetVector( const int x ) = 0;
    virtual std::vector<double> GetVector( const int x, const int y ) = 0;
    virtual bool Is1D() = 0;
    virtual bool Is2D() = 0;
    virtual bool Is3D() = 0;

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
    std::string zName;
    std::string zNameUnit;
    std::string zTitle;
    std::string zTitleUnit;
    std::string zTexTitle;
    std::string zTexTitleUnit;
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
    int GetNBinsZ(const int idx, const int idy) {return 0;}
    std::string GetBinDef( const int idx ){ return binDef[idx]; }
    std::string GetBinDef( const int x, const int y ){ return GetBinDef( x ); }
    std::string GetBinDef( const int x, const int y, const int z ){ return GetBinDef( x ); }
    int GetBinType() { return int(binType); }

    Double_t GetBinXLow( const int i ) { return -9999999; }
    Double_t GetBinXHigh( const int i ) { return -9999999; }
    Double_t GetBinYLow( const int i, const int j ) { return -9999999; }
    Double_t GetBinYHigh( const int i, const int j ) { return -9999999; }

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
    std::string GetZName() { return zName; }
    std::string GetZNameUnit() { return zNameUnit; }
    std::string GetZTitle() { return zTitle; }
    std::string GetZTitleUnit() { return zTitleUnit; }
    std::string GetZTexTitle() { return zTexTitle; }
    std::string GetZTexTitleUnit() { return zTexTitleUnit; }

    std::string GetSelection() { return fselection; }

    std::vector<double> GetVector() { return fblock; }
    std::vector<double> GetVector( const int /*x*/ ) { return GetVector(); }
    std::vector<double> GetVector( const int /*x*/, const int /*y*/ ) { return GetVector(); }
    bool Is1D() { return true; }
    bool Is2D() { return false; }
    bool Is3D() { return false; }

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
    std::string zName;
    std::string zNameUnit;
    std::string zTitle;
    std::string zTitleUnit;
    std::string zTexTitle;
    std::string zTexTitleUnit;

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
    int GetNBinsZ( const int idx, const int idy ) { return 0; }
    std::string GetBinDef( const int x, const int y, const int z ) {return this->GetBinDef(x, y);}
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
    Double_t GetBinYLow( const int i, const int j ) { return -9999999; }
    Double_t GetBinYHigh( const int i, const int j ) { return -9999999; }

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
    std::string GetZName() { return zName; }
    std::string GetZNameUnit() { return zNameUnit; }
    std::string GetZTitle() { return zTitle; }
    std::string GetZTitleUnit() { return zTitleUnit; }
    std::string GetZTexTitle() { return zTexTitle; }
    std::string GetZTexTitleUnit() { return zTexTitleUnit; }
    std::string GetSelection() { return fselection; }

    std::vector<double> GetVector() { return GetVector( 0u ); }
    std::vector<double> GetVector( const int x ) { return fblock_vv_[x]; }
    std::vector<double> GetVector( const int x, const int y ) { return fblock_vv_[x]; }
    std::map<double, std::vector<double>> GetMap() { return fblock; }
    bool Is1D() { return false; }
    bool Is2D() { return true; }
    bool Is3D() { return false; }

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
    std::string zName;
    std::string zNameUnit;
    std::string zTitle;
    std::string zTitleUnit;
    std::string zTexTitle;
    std::string zTexTitleUnit;

    std::vector< std::string > binDef;

  private:
    void Init();
    void SetTitle(const std::string& title);
    void SetTexTitle(const std::string& title);
    void SetName( const std::string& name );
};

class Block3D : public Block {
  public:
    Block3D( const std::string& name, const std::string& title,
      const std::map< double, std::map< double, std::vector<double> > >& block3d,
      const std::string& selection, int i ) : Block( name, title, i ),
      fName( name ), fTitle( title ), fblock( block3d ),
      fselection( selection ), binType( i )
      { this->Init(); }

    Block3D(const std::string& name, const std::string& title,
     const std::string& textitle,
     const std::map< double, std::map< double, std::vector<double> > >& block3d,
     const std::string& selection, int i ) : Block( name, title, i ),
     fName( name ), fTitle( title ), fTexTitle( textitle ),
     fblock( block3d ), fselection( selection ), binType( i )
     { this->Init(); }

    virtual ~Block3D() = default;
    int GetNBinsX(){ return fblock_map_.size(); }
    int GetNBinsY( const int idx ) { return fblock_map_[idx].size(); }
    int GetNBinsZ( const int idx, const int idy ) { return fblock_map_[idx][idy].size() - 1; }
    std::string GetBinDef( const int x, const int y, const int z ) {
      int index = 0;
      for( auto i = 0; i < x; i++ ) {
        for(auto j = 0; j < fblock_map_[i].size(); j++){
          index+= (fblock_map_[i][j].size() - 1);
        }
      }
        for( auto j = 0; j < y; j++ ) {
          index += (fblock_map_[x][j].size() - 1);
        }
      index += z;
      return this->GetBinDef( index );
    }
    std::string GetBinDef( const int x, const int y ) {
      int index = 0;
      for( auto i = 0; i < x; i++ ) {
        for( auto j = 0; j < y; j++ ) {
          index += (fblock_map_[i][j].size() - 1);
        }
      }
      return this->GetBinDef( index );
    }



    std::string GetBinDef( const int idx ){ 
      std::cout << "index : " << idx << "  " << binDef[idx] << std::endl;
      return binDef[idx]; }
    int GetBinType(){ return int(binType); }

    Double_t GetBinXLow( const int i ) { return xbin[i]; }
    Double_t GetBinXHigh( const int i ) { return xbin[i+1]; }
    Double_t GetBinYLow( const int i, const int j ) { return ybin_map_.at(i).at(j); }
    Double_t GetBinYHigh( const int i, const int j ) { return ybin_map_.at(i).at(j+1); }

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
    std::string GetZName() { return zName; }
    std::string GetZNameUnit() { return zNameUnit; }
    std::string GetZTitle() { return zTitle; }
    std::string GetZTitleUnit() { return zTitleUnit; }
    std::string GetZTexTitle() { return zTexTitle; }
    std::string GetZTexTitleUnit() { return zTexTitleUnit; }
    std::string GetSelection() { return fselection; }

    std::vector<double> GetVector() { return GetVector( 0u, 0u ); }
    std::vector<double> GetVector( const int x ) { return GetVector(x, 0u); }
    std::vector<double> GetVector( const int x, const int y ) { return fblock_map_.at(x).at(y); }
 //   std::map<double, std::vector<double>> GetMap() { return fblock; }
    bool Is1D() { return false; }
    bool Is2D() { return false; } // for Block more than 2D, using the same
    bool Is3D() { return true; } // for Block more than 2D, using the same

  private:

    std::string fName;
    std::string fTitle;
    std::string fTexTitle;
    std::map< double, std::map< double, std::vector<double> > > fblock;
    // using the index of x and y to store the bin edges
    std::map< int, std::map<int, std::vector<double> > >  fblock_map_;
    std::map< int, std::map<int, double > >  ybin_map_;
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

    std::string zName;
    std::string zNameUnit;
    std::string zTitle;
    std::string zTitleUnit;
    std::string zTexTitle;
    std::string zTexTitleUnit;

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

class BlockReco {
	public:

	BlockReco(Block* breco): block_reco_( breco ) {}
	virtual ~BlockReco() = default;
	Block* block_reco_;
};
