/*
  DCDriftParamMan.cc

  2018/12  K.Shirotori
*/

#include "DCDriftParamMan.hh"
#include "DetectorID.hh"
#include "ConfMan.hh"

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cstdlib>

const int BufSize = 200;

const int MaxParam = 15;

struct DCDriftParamRecord
{
  int type, np;
  double p[MaxParam];
};

DCDriftParamMan::DCDriftParamMan( const std::string & filename )
  : ParameterFileName_(filename)
{}

DCDriftParamMan::~DCDriftParamMan()
{
  clearElements();
}

void DCDriftParamMan::clearElements( void )
{
  std::map <unsigned int, DCDriftParamRecord *>::iterator itr;
  for( itr=Cont_.begin(); itr!=Cont_.end(); ++itr ){
    delete itr->second;
  }
  Cont_.clear();
}

inline unsigned int SortKey( int plane, double wire )
{
  return (plane<<10) | int(wire);
}

bool DCDriftParamMan::Initialize( void )
{
  static const std::string funcname = "[DCDriftParamMan::Initialize]";

  int pid, wid, type, np, n;
  FILE *fp;
  char buf[BufSize];
  double q[MaxParam];

  if((fp=fopen(ParameterFileName_.c_str(),"r"))==0){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }

  while(fgets(buf,BufSize,fp)){
    if(buf[0]!='#'){
      if((n=sscanf(buf,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   &pid, &wid, &type, &np, &q[0], &q[1], &q[2], &q[3], &q[4],
		   &q[5], &q[6], &q[7], &q[8], &q[9], &q[10]))>=6 ){
	// wid is not used at present
	// reserved for future updates
	unsigned int key=SortKey(pid,0);
	DCDriftParamRecord *ptr = new DCDriftParamRecord();
	if(ptr){
	  if(np>MaxParam) np=MaxParam;
	  ptr->type=type; ptr->np=np;
	  for(int i=0; i<np; ++i ) ptr->p[i]=q[i];
	  if(Cont_[key]) delete Cont_[key];
	  Cont_[key]=ptr;
	}
	else{
	  std::cerr << funcname << ": new fail. "
                    << " Plane=" << std::setw(3) << pid 
		    << " Wire=" << std::setw(3) << wid << std::endl;
	}
      }
      else{
	std::cerr << funcname << ": Bad format =>\n"
		  << std::string(buf) << std::endl;
      } /* if((n=sscanf... */
    } /* if(buf[0]...) */
  } /* while(fgets... */
  fclose(fp);

  std::cout << funcname << ": Initialization finished" << std::endl;
  return true;
}

DCDriftParamRecord * 
DCDriftParamMan::getParameter( int PlaneId, double WireId ) const
{
  WireId=0;
  unsigned int key=SortKey(PlaneId,WireId);
  return Cont_[key];
}

double DCDriftParamMan::DriftLength1( double dt, double vel )
{
  return dt*vel;
}

double DCDriftParamMan::
DriftLength2( double dt, double p1, double p2, double p3,
	      double st, double p5, double p6 )
{
  double dtmax=12.+p2+1./p6;
  double dl;
  double alph=-0.5*p5*st+0.5*st*p3*p5*(p1-st)
    +0.5*p3*p5*(p1-st)*(p1-st);
  if( dt<-10. || dt>dtmax+12. )
    dl = -500.;
  else if( dt<st )
    dl = 0.5*(p5+p3*p5*(p1-st))/st*dt*dt;
  else if( dt<p1 )
    dl = alph+p5*dt-0.5*p3*p5*(p1-dt)*(p1-dt);
  else if( dt<p2 )
    dl = alph+p5*dt;
  else
    dl = alph+p5*dt-0.5*p6*p5*(dt-p2)*(dt-p2);
  return dl;
} 

double DCDriftParamMan::DriftLength3( double dt, double p1, double p2 ,int PlaneId)
{
  if (PlaneId>=1 && PlaneId <=4) {
    if (dt > 130.)
      return 999.9;
  } else if (PlaneId == 5) {
    if (dt > 500.)
      ;
  } else if (PlaneId>=6 && PlaneId <=11) {
    if (dt > 80.)
      return 999.9;
  } else if (PlaneId>=101 && PlaneId <=124) {
    if (dt > 80.)
      return 999.9;
  }

    return dt*p1+dt*dt*p2;
}

double DCDriftParamMan::DriftLength4( double dt, double p1, double p2 ,double p3)
{
  if (dt < p2)
    return dt*p1;
  else
    return p3*(dt-p2) + p1*p2;
}

double DCDriftParamMan::DriftLength5( double dt, double p1, double p2 ,double p3, double p4, double p5)
{
  if (dt < p2)
    return dt*p1;
  else if (dt >= p2 && dt < p4)
    return p3*(dt-p2) + p1*p2;
  else
    return p5*(dt-p4) + p3*(p4-p2) + p1*p2;
}

double DCDriftParamMan::DriftLength6(int PlaneId, double dt, double p1, double p2 ,double p3, double p4, double p5)
{
  ConfMan *confMan = ConfMan::GetConfManager();
  double dtcut1 = confMan->DTCut1();
  double dtcut2 = confMan->DTCut2();
  double dtcut3 = confMan->DTCut3();
  double dtcut4 = confMan->DTCut4();
  double dlcut1 = confMan->DLCut1();
  double dlcut2 = confMan->DLCut2();
  double dlcut3 = confMan->DLCut3();
  double dlcut4 = confMan->DLCut4();

  //start
/*
  if (PlaneId-PlOffsDC>=1 && PlaneId-PlOffsDC<=3) {

    if ( PlaneId==201 && (dt < -20.0) || (dt > 400.0) )
      return 999.9;
    if ( PlaneId==202 && (dt < -20.0) || (dt > 400.0) )
    return 999.9;
    if ( PlaneId==203 && (dt < -20.0) || (dt > 400.0) )
    return 999.9;
    //end
*/

    
  if (PlaneId-PlOffsDC>=1 && PlaneId-PlOffsDC<=4) {
    if ( dt<-20 || dt>400 )
      return 999.9;
    
    double dl = dt*p1+dt*dt*p2+p3*pow(dt, 3.0)+p4*pow(dt, 4.0)+p5*pow(dt, 5.0);

    //For Test DC
    
    if ( PlaneId==201 && (dl > dlcut1 || dt > dtcut1) )
      return dlcut1;
    if ( PlaneId==202 && (dl > dlcut2 || dt > dtcut2) )
      return dlcut2;
    if ( PlaneId==203 && (dl > dlcut3 || dt > dtcut3) )
      return dlcut3;
    if ( PlaneId==204 && (dl > dlcut4 || dt > dtcut4) )
      return dlcut4;
    if (dl < 0.0)
      return 0.0;
    else
      return dl;
  }
  else {
    std::cerr << "DCDriftParamMan::DriftLength6 Invalid planeId="
              << PlaneId << std::endl;
    exit(1);
  }
}

bool DCDriftParamMan::calcDrift( int PlaneId, double WireId, double ctime, 
				 double & dt, double & dl ) const
{
  static const std::string funcname = "DCDriftParamMan::calcDrift";
  DCDriftParamRecord *p=getParameter(PlaneId,WireId);
  if(p){
    dt=p->p[0]-ctime;
    switch(p->type){
    case 1:
      dl=DriftLength1( dt, p->p[1] );
      return true;
    case 2:
      dl=DriftLength2( dt, p->p[1], p->p[2], p->p[3], p->p[4], p->p[5], p->p[6] );
      return true;
    case 3:
      dl=DriftLength3( dt, p->p[1], p->p[2], PlaneId );
      return true;
    case 4:
      dl=DriftLength4( dt, p->p[1], p->p[2], p->p[3] );
      return true;
    case 5:
      dl=DriftLength5( dt, p->p[1], p->p[2], p->p[3], p->p[4], p->p[5] );
      return true;
    case 6:
      dl=DriftLength6( PlaneId, dt, p->p[1], p->p[2], p->p[3], p->p[4], p->p[5] );
      return true;
    default:
      std::cerr << funcname << ": No Type. Type=" 
		<< p->type << std::endl;
      return false;
    }
  }
  else{
#if 0
    std::cerr << funcname << ": No record. "
              << " PlaneId=" << std::setw(3) << PlaneId
              << " WireId=" << std::setw(3) << WireId << std::endl;
#endif
    return false;
  }
}


