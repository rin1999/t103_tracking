/*
  TrGeomMan.cc

  2012/5  K.Shirotori
*/

#include "TrGeomMan.hh"
#include "TrGeomRecord.hh"

#include <string>
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cstdlib>

const int MaxChar = 200;

TrGeomMan *TrGeomMan::geomMan_=0;

TrGeomMan::TrGeomMan()
{}

TrGeomMan::~TrGeomMan()
{}


TrGeomMan & TrGeomMan::GetInstance( void )
{
  if( !geomMan_ ){
    geomMan_ = new TrGeomMan();
  }
  return *geomMan_;
}

double TrGeomMan::GetLocalZ( int lnum ) const
{
  static const std::string funcname = "[TrGeomMan::GetLocalZ(int)]"; 
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->length_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double TrGeomMan::GetResolution( int lnum ) const
{
  static const std::string funcname = "[TrGeomMan::GetResolution(int)]"; 
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->resol_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double TrGeomMan::GetTiltAngle( int lnum ) const
{
  static const std::string funcname = "[TrGeomMan::GetTiltAngle(int)]";
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->tiltAngle_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double TrGeomMan::GetRotAngle1( int lnum ) const
{
  static const std::string funcname = "[TrGeomMan::GetRotAngle1(int)]";
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->rotAngle1_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double TrGeomMan::GetRotAngle2( int lnum ) const
{
  static const std::string funcname = "[TrGeomMan::GetRotAngle2(int)]";
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->rotAngle2_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

const ThreeVector & TrGeomMan::GetGlobalPosition( int lnum ) const
{
  static const std::string funcname = "[TrGeomMan::GetGlobalPosition(int)]";
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->pos_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

ThreeVector TrGeomMan::NormalVector( int lnum ) const
{
  static const std::string funcname = "[TrGeomMan::NormalVector(int)]";
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->NormalVector();
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

ThreeVector TrGeomMan::UnitVector( int lnum ) const
{
  static const std::string funcname = "[TrGeomMan::UnitVector(int)]";
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->UnitVector();
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

const TrGeomRecord * TrGeomMan::GetRecord( int lnum ) const
{
  static const std::string funcname = "[TrGeomMan::GetRecord(int)]";
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}
 
  
double TrGeomMan::calcWirePosition( int lnum, double wire ) const
{
  static const std::string funcname = "[TrGeomMan::calcWirePosition()]";
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ){
    return pGeo->WirePos(wire);
  }
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

int TrGeomMan::calcWireNumber( int lnum, double pos ) const
{
  static const std::string funcname = "[TrGeomMan::calcWireNumber()]";
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ){
    return pGeo->WireNumber(pos);
  }
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}
 

void TrGeomMan::clearElements( void )
{
  //  for_each( geomRecord_.begin(), geomRecord_.end(), DeleteObject() );
  std::map <int, TrGeomRecord *>::iterator itr;
  for( itr=geomRecord_.begin(); itr!=geomRecord_.end(); ++itr ){
    delete itr->second;
  }
  geomRecord_.clear();
}


bool TrGeomMan::Initialize( void )
{
  static const std::string funcname = "[TrGeomMan::Initialize]";
  char str[MaxChar];
  char cname[MaxChar];
  int id;
  double xs, ys, zs, ta, ra1, ra2, l, res, w0, dd, ofs;

  FILE *fp;

  if( ( fp = fopen( filename_.c_str(), "r" ) ) == 0 ){
    throw std::invalid_argument(funcname+": file open fail");
  }

  clearElements();

  while( fgets( str, MaxChar, fp ) != 0 ){
    if( str[0]!='#' ){
      if( sscanf( str, "%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &id, cname, &xs, &ys, &zs, &ta, &ra1, &ra2, &l, &res,
		  &w0, &dd, &ofs )
	  == 13 ){
	TrGeomRecord *pRec = 
	  new TrGeomRecord( id, cname, xs, ys, zs, 
			    ta, ra1, ra2, l, res,
			    w0, dd, ofs );

	TrGeomRecord *pOld = geomRecord_[id];
	geomRecord_[id] = pRec;

	if( pOld ){
	  std::cerr << funcname << ": duplicated id number. "
		    << " following record is deleted." << std::endl;
	  std::cerr << "Id=" << pOld->id_ << " " << pOld->pos_
		    << " ) ... " << std::endl;
	  delete pOld;
	}
      }
      else {
	std::string strtemp=str;
	std::cerr << funcname << ": Invalid format " << strtemp << std::endl;
      }
    }
  }

  fclose(fp);

  return true;
}


std::vector <int> TrGeomMan::GetDetectorIDList( void ) const
{
  std::vector<int> vlist;
  vlist.reserve(geomRecord_.size());
  std::map <int, TrGeomRecord *>::const_iterator 
    itr=geomRecord_.begin(), end=geomRecord_.end();

  for(; itr!=end; ++itr ){
    vlist.push_back( itr->first );
  }

  return vlist;
}

ThreeVector TrGeomMan::Local2GlobalPos( int lnum, 
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[TrGeomMan::Local2GlobalPos(ThreeVecor &)]";

  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ) 
    throw std::out_of_range(funcname+": No record" );

  double x = pGeo->dxds_*in.x() + pGeo->dxdt_*in.y()
    + pGeo->dxdu_*in.z() + pGeo->pos_.x();
  double y = pGeo->dyds_*in.x() + pGeo->dydt_*in.y()
    + pGeo->dydu_*in.z() + pGeo->pos_.y();
  double z = pGeo->dzds_*in.x() + pGeo->dzdt_*in.y()
    + pGeo->dzdu_*in.z() + pGeo->pos_.z();

  return ThreeVector( x, y, z );
}

ThreeVector TrGeomMan::Global2LocalPos( int lnum,
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[TrGeomMan::Global2LocalPos(ThreeVecor &)]";

  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ){
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }

  double x 
    = pGeo->dsdx_*(in.x()-pGeo->pos_.x())
    + pGeo->dsdy_*(in.y()-pGeo->pos_.y())
    + pGeo->dsdz_*(in.z()-pGeo->pos_.z());
  double y 
    = pGeo->dtdx_*(in.x()-pGeo->pos_.x())
    + pGeo->dtdy_*(in.y()-pGeo->pos_.y())
    + pGeo->dtdz_*(in.z()-pGeo->pos_.z());
  double z 
    = pGeo->dudx_*(in.x()-pGeo->pos_.x())
    + pGeo->dudy_*(in.y()-pGeo->pos_.y())
    + pGeo->dudz_*(in.z()-pGeo->pos_.z());

  return ThreeVector( x, y, z );
}

ThreeVector TrGeomMan::Local2GlobalDir( int lnum,
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[TrGeomMan::Local2GlobalDir(ThreeVecor &)]";

  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ){
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }

  double x = pGeo->dxds_*in.x() + pGeo->dxdt_*in.y()
    + pGeo->dxdu_*in.z();
  double y = pGeo->dyds_*in.x() + pGeo->dydt_*in.y()
    + pGeo->dydu_*in.z();
  double z = pGeo->dzds_*in.x() + pGeo->dzdt_*in.y()
    + pGeo->dzdu_*in.z();

  return ThreeVector( x, y, z );
}

ThreeVector TrGeomMan::Global2LocalDir( int lnum,
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[TrGeomMan::Global2LocalDir(ThreeVecor &)]";

  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ){
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }

  double x = pGeo->dsdx_*in.x() + pGeo->dsdy_*in.y()+ pGeo->dsdz_*in.z();
  double y = pGeo->dtdx_*in.x() + pGeo->dtdy_*in.y()+ pGeo->dtdz_*in.z();
  double z = pGeo->dudx_*in.x() + pGeo->dudy_*in.y()+ pGeo->dudz_*in.z();

  return ThreeVector( x, y, z );
}

void TrGeomMan::SetResolution( int lnum,  double res ) const
{
  static const std::string funcname = "[TrGeomMan::SetResolution(int, double)]"; 
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) {
    pGeo->resol_=res;
    return;
  }
  else{
    std::cerr << funcname << ": No record. Layer#=" 
              << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

void TrGeomMan::SetVertex( int lnum,  ThreeVector vertex ) const
{
  static const std::string funcname = "[TrGeomMan::SetVertex]"; 
  TrGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) {
    pGeo->pos_.setX(vertex.x());
    pGeo->pos_.setY(vertex.y());
    pGeo->pos_.setZ(vertex.z());

    return;
  }
  else{
    std::cerr << funcname << ": No record. Layer#=" 
              << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

int TrGeomMan::GetDetectorId( const std::string &detName ) const
{
  const std::string funcName = "TrGeomMan::GetDetectorId";

  std::map <int, TrGeomRecord *>::const_iterator 
    itr=geomRecord_.begin(), end=geomRecord_.end();

  for(; itr!=end; ++itr ){
    if (itr->second->name_ == detName)
      return itr->second->id_;
  }

  std::cerr << funcName << " : No such detector " << detName << std::endl;
  exit(-1);
}
