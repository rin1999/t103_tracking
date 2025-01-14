/*
  TrTdcCalibMan.cc

  2019/2  K.Shirotori
*/

#include "TrTdcCalibMan.hh"

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstdlib>

const int BufSize = 200;

struct TrTdcCalMap {
  TrTdcCalMap( double q0, double q1 )
    : p0(q0), p1(q1)
  {}
  double p0, p1;
};

TrTdcCalibMan::TrTdcCalibMan( const std::string & filename )
  : MapFileName_(filename)
{}

TrTdcCalibMan::~TrTdcCalibMan()
{
  clearElements();
}

void TrTdcCalibMan::clearElements( void )
{
  std::map <unsigned int, TrTdcCalMap *>::iterator itr;
  for( itr=Cont_.begin(); itr!=Cont_.end(); ++itr ){
    delete itr->second;
  }
  Cont_.clear();
}

inline unsigned int SortKey( int plane, double fiber)
{
  return (plane<<10) | int(fiber);
}


bool TrTdcCalibMan::Initialize( void )
{
  static const std::string funcname = "[TrTdcCalibMan::Initialize]";

  int pid, wid;
  double p0, p1;
  FILE *fp;
  char buf[BufSize];

  if((fp=fopen(MapFileName_.c_str(),"r"))==0){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }

  while(fgets(buf,BufSize,fp)){
    if(buf[0]!='#'){
      if(sscanf(buf,"%d %d %lf %lf",&pid,&wid,&p1,&p0)==4){
	unsigned int key=SortKey(pid,wid);
	TrTdcCalMap *p=new TrTdcCalMap(p0,p1);
	if(p){
	  if(Cont_[key]) delete Cont_[key];
	  Cont_[key]=p;
	}
	else{
	  std::cerr << funcname << ": new fail. "
		    << " Plane=" << std::setw(3) << std::dec << pid 
		    << " Fiber=" << std::setw(3) << std::dec << wid 
		    << std::endl;
	}
      }
      else{
	std::cerr << funcname << ": Bad format => "
		  << std::string(buf) << std::endl;
      }
    }
  }
  fclose(fp);

  std::cout << funcname << ": Initialization finished" << std::endl;
  return true;
}

TrTdcCalMap * TrTdcCalibMan::getMap( int PlaneId, double FiberId ) const
{
  unsigned int key=SortKey(PlaneId,FiberId);
  return Cont_[key];
}  

bool TrTdcCalibMan::GetTime( int PlaneId, double FiberId,
			     int tdc, double & time ) const
{
  static const std::string funcname = "[TrTdcaCalibMan::GetTime]";
  TrTdcCalMap *p=getMap(PlaneId,FiberId);
  if(p){ 
    time=(p->p0)+(p->p1)*tdc;
    return true;
  }
  else{
    //    std::cerr << funcname << ": No record. "
    //	      << " PlaneId=" << std::setw(3) << std::dec << PlaneId
    //	      << " FiberId=" << std::setw(3) << std::dec << FiberId 
    //	      << std::endl;
    return false;
  }
}

bool TrTdcCalibMan::GetTdc( int PlaneId, double FiberId,
			    double time, int &tdc ) const
{
  static const std::string funcname = "[TrTdcaCalibMan::GetTdc]";
  TrTdcCalMap *p=getMap(PlaneId,FiberId);
  if(p){
    tdc=int((time-(p->p0))/(p->p1));
    return true;
  }
  else{
    //    std::cerr << funcname << ": No record. "
    //	      << " PlaneId=" << std::setw(3) << std::dec << PlaneId
    //	      << " FiberId=" << std::setw(3) << std::dec << FiberId 
    //	      << std::endl;
    return false;
  }
}
