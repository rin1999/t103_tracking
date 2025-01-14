/*
  ConfMan.hh

  2012/5 K.Shirotori
*/

#ifndef Confman_h
#define Confman_h 1

#include <string>

class VEvent;
class TrGeomMan;
class FieldMan;

class ConfMan
{

private:
  //Conf
  std::string ConfFileName_;
  static ConfMan *confManager_;

  //Analysis mode
  bool anaMode_;

  //Geometry files  
  std::string TrGeomFileName_;
  TrGeomMan *TrGeomManager_;

  //Fiels map file
  std::string FieldMapFileName_;
  double SpecFieldNMR_;
  double SpecFieldCalc_;

  //Beam resolution
  double BeamResol_;

  //Tracker position resolution
  double bSSDResol_;
  double sSSD1Resol_;
  double sSSD2Resol_;
  double BFTResol_;
  double SFTResol_;
  double AFTResol_;
  double IT1Resol_;
  double IT2Resol_;
  double ST1Resol_;
  double ST2Resol_;

  //Counter time resolution
  double T0Resol_;
  double TofResol_;
  double ITofResol1_;
  double ITofResol2_;
  double PADResol_;

public:
  ConfMan( const std::string & filename );
  ~ConfMan();

  static ConfMan *GetConfManager( void ) { return confManager_; }
  bool Initialize( void );
  void SetFileName( const std::string & filename ) { ConfFileName_=filename; }
  VEvent* EventAllocator();

  bool AnaMode( void ) const { return anaMode_; }

  //Tr
  TrGeomMan *GetTrGeomManager( void ) { return TrGeomManager_; }

  double SpecFieldNMR( void ) const { return SpecFieldNMR_; }
  double SpecFieldCalc( void ) const { return SpecFieldCalc_; }

  //Resolution
  double GetBeamResol( void ) const { return BeamResol_; }

  double GetbSSDResol( void ) const { return bSSDResol_; }
  double GetsSSD1Resol( void ) const { return sSSD1Resol_; }
  double GetsSSD2Resol( void ) const { return sSSD2Resol_; }
  double GetBFTResol( void ) const { return BFTResol_; }
  double GetSFTResol( void ) const { return SFTResol_; }
  double GetAFTResol( void ) const { return AFTResol_; }
  double GetIT1Resol( void ) const { return IT1Resol_; }
  double GetIT2Resol( void ) const { return IT2Resol_; }
  double GetST1Resol( void ) const { return ST1Resol_; }
  double GetST2Resol( void ) const { return ST2Resol_; }

  double GetT0Resol( void ) const { return T0Resol_; }
  double GetTofResol( void ) const { return TofResol_; }
  double GetITofResol1( void ) const { return ITofResol1_; }
  double GetITofResol2( void ) const { return ITofResol2_; }
  double GetPADResol( void ) const { return PADResol_; }


private:
  ConfMan( const ConfMan & );
  ConfMan & operator = ( const ConfMan & );

public:
  bool InitializeHistograms();
  bool InitializeParameterFiles();

};

#endif
