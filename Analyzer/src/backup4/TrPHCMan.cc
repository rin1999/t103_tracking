/*
  TrPHCMan.cc

  2019/2  K.Shirotori
*/

#include "TrPHCMan.hh"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <limits>
#include <fstream>
#include <cstdlib>

const double EPS = std::numeric_limits<double>::epsilon();

const int SegMask  = 0x03FF;
const int CidMask  = 0x00FF;
const int PlidMask = 0x00FF;
const int UdMask   = 0x0003;

const int SegShift  =  0;
const int CidShift  = 11;
const int PlidShift = 19;
const int UdShift   = 27;

inline int KEY( int cid, int pl, int seg, int ud )
{
  return (((cid&CidMask)<<CidShift) | ((pl&PlidMask)<<PlidShift)
          | ((seg&SegMask)<<SegShift) | ((ud&UdMask)<<UdShift) );
}

const int BufSize = 144;

TrPHCParam::TrPHCParam( int type, int np, double *parlist )
  : Type(type),NPar(np)
{
  if(np>0){
    ParList = new double [NPar];
    for( int i=0; i<NPar; ++i ) ParList[i]=parlist[i];
  }
  else
    ParList=0;
}

TrPHCParam::~TrPHCParam()
{
  if(ParList)
    delete [] ParList;
}

double TrPHCParam::DoPHC( double time, double de )
{
  static const std::string funcname = "[TrPHCParam::DoPHC]";
  double retval=time;

  switch(Type){
  case 0:
    retval=time; break;
  case 1:
    retval=type1Correction(time,de); break;
  case 2:
    retval=type2Correction(time,de); break;
  default:
    std::cerr << funcname << ": No Correction Method. Type=" 
              << Type << std::endl;
  }
  return retval;
}

double TrPHCParam::DoRPHC( double time, double de )
{
  static const std::string funcname = "[TrPHCParam::DoRPHC]";
  double retval=time;

  switch(Type){
  case 0:
    retval=time; break;
  case 1:
    retval=type1RCorrection(time,de); break;
  case 2:
    retval=type2RCorrection(time,de); break;
  default:
    std::cerr << funcname << ": No Correction Method. Type=" 
              << Type << std::endl;
  }
  return retval;
}

double TrPHCParam::type1Correction( double time, double de )
{
  if(fabs(de-ParList[1])<EPS) de=ParList[1]+EPS;
  return time-ParList[0]/sqrt(fabs(de-ParList[1]))+ParList[2];
}

double TrPHCParam::type2Correction( double time, double w )
{
  if(fabs(w-ParList[1])<EPS) w=ParList[1]+EPS;
  return time-(ParList[0]*w*w + ParList[1]*w + ParList[2]);
}

double TrPHCParam::type1RCorrection( double time, double de )
{
  if(fabs(de-ParList[1])<EPS) de=ParList[1]+EPS;
  return time+ParList[0]/sqrt(fabs(de-ParList[1]))-ParList[2];
}

double TrPHCParam::type2RCorrection( double time, double w )
{
  if(fabs(w-ParList[1])<EPS) w=ParList[1]+EPS;
  return time-(ParList[0]*w*w + ParList[1]*w + ParList[2]);
}

TrPHCMan::TrPHCMan( const std::string & filename )
  : PHCFileName(filename)
{
}

TrPHCMan::~TrPHCMan()
{
  clearMap();
}

void TrPHCMan::clearMap( void )
{
 for(PhcPIterator i=Container.begin(); i!=Container.end(); i++)
    delete i->second;
 Container.clear();
}

bool TrPHCMan::Initialize( void )
{
  static const std::string funcname = "[TrPHCMan::Initialize]";

  //  FILE *fp;
  std::ifstream f(PHCFileName.c_str());
  //  char buf[BufSize];

  if(f.fail()){
    std::cerr << funcname << ": file open fail :" <<
      PHCFileName << std::endl;
    exit(-1);
  }

  int cid, plid, seg, ud, type, np;
  double par[10];

  while(f.good() ){
    std::string buf;
    std::getline(f,buf);
    if(buf.empty())
      continue;
    if(buf[0]!='#'){
      if( sscanf(buf.c_str(),
                 "%d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &cid, &plid, &seg, &ud, &type, &np, &par[0], &par[1], &par[2],
                 &par[3], &par[4], &par[5], &par[6], &par[7], &par[8], 
                 &par[9] ) >= 6 ){
        if(np>10) np=10;
        int key=KEY(cid,plid,seg,ud);
        TrPHCParam *p=new TrPHCParam(type,np,par);
        if(p) Container[key]=p;
        else{
          std::cerr << funcname << ": new fail." << std::endl;
          std::cerr << " PlId=" << std::setw(2) << plid
                    << " Seg=" << std::setw(3) << seg
                    << " Ud=" << std::setw(1) << ud << std::endl;
        }
      }
      else{
        std::cerr << funcname << ": Invalid format" << std::endl;
        std::cerr << " ===> " << buf << std::endl;
      } /* if( sscanf... ) */ 
    } /* if(buf[0]...) */ 
  } /* while( fgets... ) */

  //  fclose(fp);
  f.close();
  std::cout << funcname << ": Initialization finished." << std::endl;
  return true;
}

bool TrPHCMan::doCorrection( int cid, int plid, int seg, int ud,
			     double time, double de, double & ctime )
{
  TrPHCParam *map=GetMap(cid,plid,seg,ud);
  if(!map){ ctime=time; return false; }
  ctime=map->DoPHC(time,de);
  return true;
}

bool TrPHCMan::doRCorrection( int cid, int plid, int seg, int ud,
			      double time, double de, double & ctime )
{
  TrPHCParam *map=GetMap(cid,plid,seg,ud);
  if(!map){ ctime=time; return false; }
  ctime=map->DoRPHC(time,de);
  return true;
}

TrPHCParam * TrPHCMan::GetMap( int cid, int plid, int seg, int ud )
{
  int key=KEY(cid,plid,seg,ud);
  TrPHCParam *map=0;
  PhcPIterator i=Container.find(key);
  if( i != Container.end() ) map=i->second;
  return map;
}
