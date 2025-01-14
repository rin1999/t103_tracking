/*
  UserSpecTracking.cc

  2016/2 K.Shirotori 
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>

#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "Particle.hh"
#include "RawData.hh"
#include "PrimInfo.hh"
#include "SpecLib.hh"

#define check 1

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventSpecTracking 
  : public VEvent
{

private:
  RawData *rawData;
  TrAnalyzer *TrAna;
  HodoAnalyzer *HodoAna;

public:
  EventSpecTracking();
  ~EventSpecTracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( std::ifstream & );
  bool InitializeHistograms(); 
  void InitializeEvent();
};

EventSpecTracking::EventSpecTracking()
  : VEvent(),
    rawData(0),
    TrAna(new TrAnalyzer()),
    HodoAna(new HodoAnalyzer())
{
}

EventSpecTracking::~EventSpecTracking()
{
  if (HodoAna)   delete HodoAna;
  if (TrAna)   delete TrAna;
  if (rawData) delete rawData;
}

#ifndef MaxHits 
#define MaxHits 30
#endif

struct Event{
  //Primary
  int    prinhits;
  double priposx, priposy, priposz;
  double pbeam, ubeam, vbeam, abmom;
  double prim1, prip1;
  double pritheta1, priphi1, prithetacm1, priphicm1;
  double prim2, prip2;
  double pritheta2, priphi2, prithetacm2, priphicm2;

  //Generated Beam
  int    gbnhits;
  std::vector<int>    gbid, gbtype;
  std::vector<double> gbp, gbpx, gbpy, gbpz;
  std::vector<int>    gbpid;
  std::vector<double> gbvx, gbvy;

  //Generated Scat
  int    gsnhits;
  std::vector<int>    gsid, gstype;
  std::vector<double> gsp, gspx, gspy, gspz;
  std::vector<int>    gspid;
  std::vector<double> gsvx, gsvy;

  //Beam Tracking
  int    ntb;
  std::vector<int>    idb, typeb;
  std::vector<int>    layerb;
  std::vector<double> chisqrb;
  std::vector<double> x0b, y0b;
  std::vector<double> u0b, v0b;
  std::vector<double> posb, resb;

  //Scat Tracking
  std::vector<int>    nts;
  std::vector<int>    ids, types;
  std::vector<int>    layers;
  std::vector<double> chisqrs;
  std::vector<double> p, theta, phi;
  std::vector<double> dp;
  std::vector<double> path, m2;
  std::vector<double> x0s, y0s;
  std::vector<double> u0s, v0s;
  std::vector<double> poss, ress;

  //PID
  int ntKp,  ntKm;
  int ntPip, ntPim;
  int ntP,   ntPb;
  int ntMup, ntMum;
  int ntEp,  ntEm;

  std::vector<double> pkp,  ukp,  vkp,  thetakp,  phikp;  
  std::vector<double> pkm,  ukm,  vkm,  thetakm,  phikm;  
  std::vector<double> ppip, upip, vpip, thetapip, phipip;;
  std::vector<double> ppim, upim, vpim, thetapim, phipim;;
  std::vector<double> pp,   up,   vp,   thetap,   phip;   
  std::vector<double> ppb,  upb,  vpb,  thetapb,  phipb;  
  std::vector<double> pmup, umup, vmup, thetamup, phimup;
  std::vector<double> pmum, umum, vmum, thetamum, phimum;
  std::vector<double> pep,  uep,  vep,  thetaep,  phiep;  
  std::vector<double> pem,  uem,  vem,  thetaem,  phiem;

  //PID flag
  std::vector<int> ntdet;
  int ntrich,  ntpid1, ntpid2, nttof, ntitof, ntpad, ntmf;
  std::vector<int> fkp, fkm, fpip, fpim, fp, fpb, fmup, fmum, fep, fem;

  //Output
  double pkpx[MaxHits],  pkpy[MaxHits],  pkpz[MaxHits];
  double pkmx[MaxHits],  pkmy[MaxHits],  pkmz[MaxHits];
  double ppipx[MaxHits], ppipy[MaxHits], ppipz[MaxHits];
  double ppimx[MaxHits], ppimy[MaxHits], ppimz[MaxHits];
  double ppx[MaxHits],   ppy[MaxHits],   ppz[MaxHits];
  double ppbx[MaxHits],  ppby[MaxHits],  ppbz[MaxHits];
  double pmupx[MaxHits], pmupy[MaxHits], pmupz[MaxHits];
  double pmumx[MaxHits], pmumy[MaxHits], pmumz[MaxHits];
  double pepx[MaxHits],  pepy[MaxHits],  pepz[MaxHits];
  double pemx[MaxHits],  pemy[MaxHits],  pemz[MaxHits];

  int ntB;
  double pB[MaxHits],  uB[MaxHits],  vB[MaxHits];
  double pBx[MaxHits], pBy[MaxHits], pBz[MaxHits];

  int flkp[MaxHits], flkm[MaxHits];
  int flpip[MaxHits], flpim[MaxHits];
  int flp[MaxHits], flpb[MaxHits];
  int flmup[MaxHits], flmum[MaxHits];
  int flep[MaxHits], flem[MaxHits];
};
static Event event;

bool EventSpecTracking::ProcessingBegin()
{
 return true;
}

bool EventSpecTracking::ProcessingNormal( std::ifstream &In )
{
  const std::string funcname = "ProcessingNormal";

  rawData = new RawData;
  if( !rawData->DecodeRawHits(In) ) return false;

#if check
  std::cout << "*************************************************" << std::endl; 
#endif

  const TrGeomMan & geomMan=TrGeomMan::GetInstance();

  //**************************************************************************
  //******************RawData

  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  double X0, Y0, Z0;
  double BeamMom;
  //PrimaryInfo
  {
    const PrimInfoContainer &cont=rawData->GetPrimHC();  
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      PrimInfo *hit=cont[i];
      event.priposx = hit->GetVertX();
      event.priposy = hit->GetVertY();
      event.priposz = hit->GetVertZ();
      X0=hit->GetVertX();
      Y0=hit->GetVertY();
      Z0=hit->GetVertZ();
      
      event.pbeam = hit->GetBeamMom();
      event.ubeam = hit->GetBeamU();
      event.vbeam = hit->GetBeamV();
      event.abmom = hit->GetAnaBeamMom();
      BeamMom=hit->GetAnaBeamMom();
      double ubeam = hit->GetBeamU();
      double vbeam = hit->GetBeamV();

      event.prim1 = hit->GetMass1();
      event.prip1 = hit->GetMom1();
      event.pritheta1 = hit->GetTheta1();
      event.priphi1 = hit->GetPhi1();
      event.prithetacm1 = hit->GetThetaCM1();
      event.priphicm1 = hit->GetPhiCM1();
      
      event.prim2 = hit->GetMass2();
      event.prip2 = hit->GetMom2();
      event.pritheta2 = hit->GetTheta2();
      event.priphi2 = hit->GetPhi2();
      event.prithetacm2 = hit->GetThetaCM2();
      event.priphicm2 = hit->GetPhiCM2();
    }
  }  
  //ThreeVector IniVert(0.0, 0.0, 0.0);
  ThreeVector IniVert(1400.0-Z0, 0.0, 0.0);

  //Beam
  std::vector <ThreeVector> BeamPCont, BeamXCont;
  ThreeVector primPos(X0, Y0, Z0);
  {
    const s_BeamRHitContainer &cont =rawData->GetsBeamRHC();
    int nhB = cont.size();
    event.gbnhits = nhB;
    for(int i=0; i<nhB; i++){
      s_BeamRawHit *hit=cont[i];
      event.gbid.push_back(hit->TrackId());
      event.gbtype.push_back(hit->TrackType());
      event.gbp.push_back(hit->GetMom());
      event.gbpx.push_back(hit->GetMomX());
      event.gbpy.push_back(hit->GetMomY());
      event.gbpz.push_back(hit->GetMomZ());
      event.gbpid.push_back(hit->GetPid());
      event.gbvx.push_back(hit->GetVertX());
      event.gbvy.push_back(hit->GetVertY());

      //////////////////////////Tracking
      event.idb.push_back(hit->TrackId());
      event.typeb.push_back(hit->TrackType());

      TrAna->DecodesBRawHits( hit, hit->TrackType() );
      TrAna->TrackSearchBFTT();
     
      int nt=TrAna->GetNtracksBFTT();
      for( int it=0; it<nt; ++it ){
	TrLocalTrack *tp=TrAna->GetTrackBFTT(it);
	double chisqr=tp->GetChiSquare();
	double x0=tp->GetX0(), y0=tp->GetY0();
	double u0=tp->GetU0(), v0=tp->GetV0();
	double xtgt=tp->GetX( Z0 ), ytgt=tp->GetY( Z0 );
	double utgt=u0, vtgt=v0;

	double pt=BeamMom/sqrt(1.+u0*u0+v0*v0);
	ThreeVector Mom( pt*u0, pt*v0, pt );
	ThreeVector Pos( X0, Y0, Z0 );
	BeamPCont.push_back(Mom); 
	BeamXCont.push_back(Pos);      
      }
      TrAna->TrackSearchBeam( IniVert );

      int ntBeam=TrAna->GetNTracksBeam();
      event.ntb = ntBeam;
      for( int it=0; it<ntBeam; ++it ){
	BeamTrack *tp=TrAna->GetBeamTrack(it);
	if(!tp) continue;
	int nh=tp->GetNHits();	

	for( int j=0; j<nh; ++j ){
	  TrackHit *hit=tp->GetHit(j);
	  if(!hit) continue;
	  int layerId=hit->GetLayer();
	  event.layerb.push_back(layerId);  
	  double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	  event.posb.push_back(pos);
	  event.resb.push_back(res);
	}

	double chisqr=tp->chisqr();  
	ThreeVector Ppos=tp->PrimaryPosition();
	ThreeVector Pmom=tp->PrimaryMomentum();
	double x=Ppos.x(), y=Ppos.y(), z=Ppos.z();
	double p=Pmom.mag();
	double pz=Pmom.z();
	double u=Pmom.x()/pz, v=Pmom.y()/pz;
	
	// std::cout<< "V1= " << ThreeVector(X0, Y0, 0.0) << std::endl;
	// std::cout<< "V2= " << Ppos << std::endl;
	
	event.chisqrb.push_back(chisqr);
	event.x0b.push_back(x);
	event.y0b.push_back(y);
	event.u0b.push_back(u);
	event.v0b.push_back(v);   

	// double pt=BeamMom/sqrt(1.+u*u+v*v);
	// ThreeVector Mom( pt*u, pt*v, pt );
	// ThreeVector Pos( X0, Y0, Z0 );
	// BeamPCont.push_back(Mom); 
	// BeamXCont.push_back(Pos);      
      }
    }
  }
  // if( BeamPCont.size()==0 ){
  //   return true;
  // }
#if check
  std::cout<< "Beam: " << BeamPCont.size() << ", pB= "<< BeamMom << std::endl;
#endif

  //Scat
  double IniP;

  //PID
  Particle *part = new Particle();

  int ntRich=0, ntPid1=0, ntPid2=0; 
  int ntTof=0, ntItof=0, ntPad=0, ntMF=0;
  {
    const s_ScatRHitContainer &cont =rawData->GetsScatRHC();
    int nhS = cont.size();
    event.gsnhits = nhS;
    for(int i=0; i<nhS; i++){
      s_ScatRawHit *hit=cont[i];

      event.gsid.push_back(hit->TrackId());
      event.gstype.push_back(hit->TrackType());
      event.gsp.push_back(hit->GetMom());
      event.gspx.push_back(hit->GetMomX());
      event.gspy.push_back(hit->GetMomY());
      event.gspz.push_back(hit->GetMomZ());
      event.gspid.push_back(hit->GetPid());
      event.gsvx.push_back(hit->GetVertX());
      event.gsvy.push_back(hit->GetVertY());

      //////////////////////////Tracking
      event.ids.push_back(hit->TrackId());
      event.types.push_back(hit->TrackType());
      IniP=hit->GetMom();

      TrAna->DecodesSRawHits( hit, hit->TrackType() );
      
      TrAna->TrackSearchSFTT();
      TrAna->TrackSearchIT1T();
      TrAna->TrackSearchIT2RT();
      TrAna->TrackSearchIT2LT();
      TrAna->TrackSearchST1T();
      TrAna->TrackSearchST2T();
      
      TrAna->TrackSearchPreIn( IniP );
      TrAna->TrackSearchPreOut( IniP );
      TrAna->TrackSearchPreOut2( IniP );

      int ntSFT     = TrAna->GetNtracksSFTT();
      int ntIT1     = TrAna->GetNtracksIT1T();
      int ntIT2R    = TrAna->GetNtracksIT2RT();
      int ntIT2L    = TrAna->GetNtracksIT2LT();
      int ntST1     = TrAna->GetNtracksST1T();
      int ntST2     = TrAna->GetNtracksST2T();
      int ntPreIn   = TrAna->GetNTracksPreIn();
      int ntPreOut  = TrAna->GetNTracksPreOut();
      int ntPreOut2 = TrAna->GetNTracksPreOut2();

#if check
      std::cout<< "** Type: " << hit->TrackType() <<std::endl;
      std::cout<< "ntSFT= " << ntSFT <<std::endl;
      std::cout<< "ntIT1= " << ntIT1 <<std::endl;
      std::cout<< "ntIT2R= " << ntIT2R <<std::endl;
      std::cout<< "ntIT2L= " << ntIT2L <<std::endl;
      std::cout<< "ntST1= " << ntST1 <<std::endl;
      std::cout<< "ntST2= " << ntST2 <<std::endl;
      std::cout<< "ntPreIn= " << ntPreIn <<std::endl;
      std::cout<< "ntPreOut= " << ntPreOut <<std::endl;
      std::cout<< "ntPreOut2= " << ntPreOut2 <<std::endl;
#endif

      int Type=0, ntScat=0;
      if( hit->TrackType() == TrackTypeScat1 && ntPreIn!=0 && ntPreOut!=0 ){
      	TrAna->TrackSearchScat1( IniP, IniVert );
      	ntScat=TrAna->GetNTracksScat1();
      	Type=TrackTypeScat1;	
      }
      if( hit->TrackType() == TrackTypeScat2 && ntPreIn!=0 ){
      	if( ntST1!=0 || ntIT2R!=0 || ntIT2L!=0 ){
      	  TrAna->TrackSearchScat2( IniP, IniVert, hit );
      	  ntScat=TrAna->GetNTracksScat2();
      	  Type=TrackTypeScat2;	
      	}
      	if( ntST1==0 && ntIT2R==0 && ntIT2L==0 ){
      	  TrAna->TrackSearchScat3( IniP, IniVert );
      	  ntScat=TrAna->GetNTracksScat3();
      	  Type=TrackTypeScat2A;	
      	}
      }
      if( hit->TrackType() == TrackTypeScat3 ){
      	TrAna->TrackSearchScat3( IniP, IniVert );
      	ntScat=TrAna->GetNTracksScat3();
      	Type=TrackTypeScat3;	
      }

      event.nts.push_back(ntScat);
      for( int it=0; it<ntScat; ++it ){
	int nh;
	double chisqr;  
	ThreeVector Ppos;
	ThreeVector Pmom;
	double pathL;

	if( Type==TrackTypeScat1 ){
	  Scat1Track *tp=TrAna->GetScat1Track(it);
	  if(!tp) continue;
	  nh=tp->GetNHits();
	  chisqr=tp->chisqr();  
	  Ppos=tp->PrimaryPosition();
	  Pmom=tp->PrimaryMomentum();
	  pathL=tp->PathLengthToTOF();
	  
	  for( int j=0; j<nh; ++j ){
	    TrackHit *hit=tp->GetHit(j);
	    if(!hit) continue;
	    int layerId=hit->GetLayer();
	    event.layers.push_back(layerId);  
	    double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	    event.poss.push_back(pos);
	    event.ress.push_back(res);
	  }
	}
	if( Type==TrackTypeScat2 ){
	  Scat2Track *tp=TrAna->GetScat2Track(it);
	  if(!tp) continue;
	  nh=tp->GetNHits();
	  chisqr=tp->chisqr();  
	  Ppos=tp->PrimaryPosition();
	  Pmom=tp->PrimaryMomentum();
	  pathL=tp->PathLengthToTOF();
	  
	  for( int j=0; j<nh; ++j ){
	    TrackHit *hit=tp->GetHit(j);
	    if(!hit) continue;
	    int layerId=hit->GetLayer();
	    event.layers.push_back(layerId);  
	    double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	    event.poss.push_back(pos);
	    event.ress.push_back(res);
	  }
	}
	if( Type==TrackTypeScat2A ){
	  Scat3Track *tp=TrAna->GetScat3Track(it);
	  if(!tp) continue;
	  nh=tp->GetNHits();
	  chisqr=tp->chisqr();  
	  Ppos=tp->PrimaryPosition();
	  Pmom=tp->PrimaryMomentum();
	  pathL=tp->PathLengthToTOF();
	  
	  for( int j=0; j<nh; ++j ){
	    TrackHit *hit=tp->GetHit(j);
	    if(!hit) continue;
	    int layerId=hit->GetLayer();
	    event.layers.push_back(layerId);  
	    double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	    event.poss.push_back(pos);
	    event.ress.push_back(res);
	  }
	}
	if( Type==TrackTypeScat3 ){
	  Scat3Track *tp=TrAna->GetScat3Track(it);
	  if(!tp) continue;
	  nh=tp->GetNHits();
	  chisqr=tp->chisqr();  
	  Ppos=tp->PrimaryPosition();
	  Pmom=tp->PrimaryMomentum();
	  pathL=tp->PathLengthToTOF();
	  
	  for( int j=0; j<nh; ++j ){
	    TrackHit *hit=tp->GetHit(j);
	    if(!hit) continue;
	    int layerId=hit->GetLayer();
	    event.layers.push_back(layerId);  
	    double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	    event.poss.push_back(pos);
	    event.ress.push_back(res);
	  }
	}
	
	double x=Ppos.x(), y=Ppos.y();
	double p=Pmom.mag();
	double pz=Pmom.z();
	double u, v;
	double theta, phi;	  
	ThreeVector Mom;
	
	if(pz>0){
	  u=Pmom.x()/pz, v=Pmom.y()/pz;
	  theta = Pmom.theta();	  
	  phi = Pmom.phi(); 
	  Mom.setX(Pmom.x());
	  Mom.setY(Pmom.y());
	  Mom.setZ(Pmom.z());
	  //std::cout<< Pmom.mag() << " " << Pmom << std::endl;
	}
	if(pz<0){
	  u=Pmom.x()/pz, v=Pmom.y()/pz;
	  theta = (-1.*Pmom).theta();	  
	  phi = (-1.*Pmom).phi();	  
	  Mom.setX(-1.*Pmom.x());
	  Mom.setY(-1.*Pmom.y());
	  Mom.setZ(-1.*Pmom.z());
	  //std::cout<< Pmom.mag() << " " << Pmom << std::endl;
	}
	
	// std::cout<< "V1= " << ThreeVector(X0, Y0, 0.0) << std::endl;
	// std::cout<< "V2= " << Ppos << std::endl;
	
	event.chisqrs.push_back(chisqr);
	event.p.push_back(p);
	event.theta.push_back(theta*Rad2Deg);
	event.phi.push_back(phi*Rad2Deg);
	event.path.push_back(pathL);
	event.x0s.push_back(x);
	event.y0s.push_back(y);
	event.u0s.push_back(u);
	event.v0s.push_back(v);   
	
	event.dp.push_back(p-IniP);
#if check
	std::cout<< "****** dP= " << p-IniP <<std::endl;
#endif
	
	//PID
	if( hit->TrackType() == TrackTypeScat1 ){
	  bool flag_rich =false;
	  bool flag_pid1 =false;
	  
	  //TOF
	  const s_HodoRHitContainer &contTOF =hit->GetsTOFRHC();
	  int nhTOF = contTOF.size();
	  for(int j=0; j<nhTOF; j++){
	    s_HodoRawHit *hitTOF=contTOF[j];
	    
	    int nt = hitTOF->GetSize();
	    for( int k=0; k<nt; k++ ){
	      double time = hitTOF->GetTime(k);
	      double path = hitTOF->GetPath(k);
	      double mom = hitTOF->GetMom(k);
	      
	      double C = 2.99792458E+8;
	      double T0, V0, B0, m2_0;
	      
	      V0 = (path/1000.)/(time*1.0E-9);
	      B0 = V0/C;
	      m2_0 = (mom/B0)*(mom/B0)*(1.-B0*B0);
	      
	      event.m2.push_back(m2_0); 
	    }
	  }
	  
	  const s_HodoRHitContainer &contRICH =hit->GetsRICHRHC();
	  int nhRICH = contRICH.size();
	  if( nhRICH!=0 ) flag_rich=true;
	  const s_HodoRHitContainer &contPID1 =hit->GetsPID1RHC();
	  int nhPID1 = contPID1.size();
	  if( nhPID1!=0 ) flag_pid1=true;
	  
	  if( flag_rich ){
#if check
	    std::cout<<"******1: RICH"<<std::endl;
	    std::cout<< "id= " << hit->TrackId()
		     << " PID= " << hit->GetPid()
		     << " Mom= " << p 
		     <<std::endl;
#endif
	    part->SetId(hit->TrackId());
	    part->SetPid(hit->GetPid());
	    part->SetMomX(Mom.x());
	    part->SetMomY(Mom.y());
	    part->SetMomZ(Mom.z());
	    part->SetDetFlag(1);
	    event.ntdet.push_back(1);
	    ntRich++;
	  }
	  if( flag_pid1 && !flag_rich ){
#if check
	    std::cout<<"******2: PID1"<<std::endl;
	    std::cout<< "id= " << hit->TrackId()
		     << " PID= " << hit->GetPid()
		     << " Mom= " << p 
		     <<std::endl;
#endif
	    part->SetId(hit->TrackId());
	    part->SetPid(hit->GetPid());
	    part->SetMomX(Mom.x());
	    part->SetMomY(Mom.y());
	    part->SetMomZ(Mom.z());
	    part->SetDetFlag(2);
	    event.ntdet.push_back(2);
	    ntPid1++;
	  }
	  if( !flag_pid1 && !flag_rich ){
#if check
	    std::cout<<"******3: TOF"<<std::endl;
	    std::cout<< "id= " << hit->TrackId()
		     << " PID= " << hit->GetPid()
		     << " Mom= " << p 
		     <<std::endl;
#endif
	    part->SetId(hit->TrackId());
	    part->SetPid(hit->GetPid());
	    part->SetMomX(Mom.x());
	    part->SetMomY(Mom.y());
	    part->SetMomZ(Mom.z());
	    part->SetDetFlag(3);
	    event.ntdet.push_back(3);
	    ntTof++;
	  }
	}

	if( hit->TrackType() == TrackTypeScat2 ){
	  bool flag_tof =false;
	  bool flag_pid1 =false;
	  bool flag_pid2 =false;
	  
	  //ITOF
	  const s_HodoRHitContainer &contITOF =hit->GetsITOFRHC();
	  int nhITOF = contITOF.size();
	  for(int j=0; j<nhITOF; j++){
	    s_HodoRawHit *hitITOF=contITOF[j];
	    
	    int nt = hitITOF->GetSize();
	    for( int k=0; k<nt; k++ ){
	      double time = hitITOF->GetTime(k);
	      double path = hitITOF->GetPath(k);
	      double mom = hitITOF->GetMom(k);
	      
	      double C = 2.99792458E+8;
	      double T0, V0, B0, m2_0;
	      
	      V0 = (path/1000.)/(time*1.0E-9);
	      B0 = V0/C;
	      m2_0 = (mom/B0)*(mom/B0)*(1.-B0*B0);
	      
	      event.m2.push_back(m2_0); 
	    }
	  }
	  
	  const s_HodoRHitContainer &contTOF =hit->GetsTOFRHC();
	  int nhTOF = contTOF.size();
	  if( nhTOF!=0 ) flag_tof=true;
	  const s_HodoRHitContainer &contPID1 =hit->GetsPID1RHC();
	  int nhPID1 = contPID1.size();
	  if( nhPID1!=0 ) flag_pid1=true;
	  const s_HodoRHitContainer &contPID2 =hit->GetsPID2RHC();
	  int nhPID2 = contPID2.size();
	  if( nhPID2!=0 ) flag_pid2=true;
	  
	  if( !flag_tof && flag_pid1 ){
#if check
	    std::cout<<"******2: PID1"<<std::endl;
	    std::cout<< "id= " << hit->TrackId()
		     << " PID= " << hit->GetPid()
		     << " Mom= " << p 
		     <<std::endl;
#endif
	    part->SetId(hit->TrackId());
	    part->SetPid(hit->GetPid());
	    part->SetMomX(Mom.x());
	    part->SetMomY(Mom.y());
	    part->SetMomZ(Mom.z());
	    part->SetDetFlag(2);
	    event.ntdet.push_back(2);
	    ntPid1++;
	  }
	  if( !flag_tof && !flag_pid1 && flag_pid2 ){
#if check
	    std::cout<<"******4: PID2"<<std::endl;
	    std::cout<< "id= " << hit->TrackId()
		     << " PID= " << hit->GetPid()
		     << " Mom= " << p 
		     <<std::endl;
#endif
	    part->SetId(hit->TrackId());
	    part->SetPid(hit->GetPid());
	    part->SetMomX(Mom.x());
	    part->SetMomY(Mom.y());
	    part->SetMomZ(Mom.z());
	    part->SetDetFlag(4);
	    event.ntdet.push_back(4);
	    ntPid2++;
	  }
	  if( !flag_tof && !flag_pid1 && !flag_pid2 ){
#if check
	    std::cout<<"******5: ITOF"<<std::endl;
	    std::cout<< "id= " << hit->TrackId()
		     << " PID= " << hit->GetPid()
		     << " Mom= " << p 
		     <<std::endl;
#endif
	    part->SetId(hit->TrackId());
	    part->SetPid(hit->GetPid());
	    part->SetMomX(Mom.x());
	    part->SetMomY(Mom.y());
	    part->SetMomZ(Mom.z());
	    part->SetDetFlag(5);
	    event.ntdet.push_back(5);
	    ntItof++;
	  }
	}
	
	if( hit->TrackType() == TrackTypeScat3 ){
	  bool flag_tof =false;
	  bool flag_itof =false;
	  bool flag_rich =false;
	  bool flag_pid1 =false;
	  bool flag_pid2 =false;
	  
	  //PAD
	  const s_HodoRHitContainer &contPAD =hit->GetsPADRHC();
	  int nhPAD = contPAD.size();
	  for(int j=0; j<nhPAD; j++){
	    s_HodoRawHit *hitPAD=contPAD[j];
	    
	    int nt = hitPAD->GetSize();
	    for( int k=0; k<nt; k++ ){
	      double time = hitPAD->GetTime(k);
	      double path = hitPAD->GetPath(k);
	      double mom = hitPAD->GetMom(k);
	      
	      double C = 2.99792458E+8;
	      double T0, V0, B0, m2_0;
	      
	      V0 = (path/1000.)/(time*1.0E-9);
	      B0 = V0/C;
	      m2_0 = (mom/B0)*(mom/B0)*(1.-B0*B0);
	      
	      event.m2.push_back(m2_0); 
	    }
	  }
	  
	  const s_HodoRHitContainer &contTOF =hit->GetsTOFRHC();
	  int nhTOF = contTOF.size();
	  if( nhTOF!=0 ) flag_tof=true;
	  const s_HodoRHitContainer &contITOF =hit->GetsITOFRHC();
	  int nhITOF = contITOF.size();
	  if( nhITOF!=0 ) flag_itof=true;
	  const s_HodoRHitContainer &contRICH =hit->GetsRICHRHC();
	  int nhRICH = contRICH.size();
	  if( nhRICH!=0 ) flag_rich=true;
	  const s_HodoRHitContainer &contPID1 =hit->GetsPID1RHC();
	  int nhPID1 = contPID1.size();
	  if( nhPID1!=0 ) flag_pid1=true;
	  const s_HodoRHitContainer &contPID2 =hit->GetsPID2RHC();
	  int nhPID2 = contPID2.size();
	  if( nhPID2!=0 ) flag_pid2=true;
	  
	  if( !flag_tof && !flag_itof && !flag_rich && !flag_pid1 && !flag_pid2 ){
#if check
	    std::cout<<"******6: PAD"<<std::endl;
	    std::cout<< "id= " << hit->TrackId()
		     << " PID= " << hit->GetPid()
		     << " Mom= " << p 
		     <<std::endl;
#endif
	    part->SetId(hit->TrackId());
	    part->SetPid(hit->GetPid());
	    part->SetMomX(Mom.x());
	    part->SetMomY(Mom.y());
	    part->SetMomZ(Mom.z());
	    part->SetDetFlag(6);
	    event.ntdet.push_back(6);
	    ntPad++;
	  }
	}//PID END	  
      }
    }
  }//Scat analysis

  event.ntrich=ntRich;
  event.ntpid1=ntPid1;
  event.ntpid2=ntPid2;
  event.nttof =ntTof;
  event.ntitof=ntItof;
  event.ntpad =ntPad;
  event.ntmf  =ntMF;

  std::vector <ThreeVector> KaonPCont, KaonPXCont; 
  std::vector <ThreeVector> KaonMCont, KaonMXCont; 
  std::vector <ThreeVector> PionPCont, PionPXCont;
  std::vector <ThreeVector> PionMCont, PionMXCont;
  std::vector <ThreeVector> ProtonCont, ProtonXCont;
  std::vector <ThreeVector> ProtonbCont, ProtonbXCont;
  std::vector <ThreeVector> MuonPCont, MuonPXCont; 
  std::vector <ThreeVector> MuonMCont, MuonMXCont; 
  std::vector <ThreeVector> PositronCont, PositronXCont; 
  std::vector <ThreeVector> ElectronCont, ElectronXCont; 

  int ikp=0, ikm=0;
  int ipip=0, ipim=0;
  int ip=0, ipb=0;
  int iep=0, iem=0;
  int imup=0, imum=0;
  
  int ntPart = part->GetSize();
  for( int i=0; i<ntPart; i++ ){
    int pid=part->GetPid(i);
    double px=part->GetMomX(i), py=part->GetMomY(i), pz=part->GetMomZ(i);
    ThreeVector mom(px, py, pz);
    double p=mom.mag(), u=mom.x()/p, v=mom.y()/p;
    double theta = mom.theta()*Rad2Deg, phi=mom.phi()*Rad2Deg;
    int fldet = part->GetDetFlag(i);
    
    if( pid==4 ){//K+
      event.pkp.push_back(p);
      event.ukp.push_back(u);
      event.vkp.push_back(v);
      event.thetakp.push_back(theta);	    
      event.phikp.push_back(phi);
      event.pkpx[ikp]=mom.x();
      event.pkpy[ikp]=mom.y();
      event.pkpz[ikp]=mom.z();
      
      event.fkp.push_back(fldet);
      event.flkp[ikp]=fldet;
      
      KaonPCont.push_back(mom);
      KaonPXCont.push_back(primPos);
      ikp++;
    }
    if( pid==3 ){//K-
      event.pkm.push_back(p);
      event.ukm.push_back(u);
      event.vkm.push_back(v);
      event.thetakm.push_back(theta);	    
      event.phikm.push_back(phi);
      event.pkmx[ikm]=mom.x();
      event.pkmy[ikm]=mom.y();
      event.pkmz[ikm]=mom.z();
      
      event.fkm.push_back(fldet);
      event.flkm[ikm]=fldet;
      
      KaonMCont.push_back(mom);
      KaonMXCont.push_back(primPos);
      ikm++;
    }
    if( pid==2 ){//Pi+
      event.ppip.push_back(p);
      event.upip.push_back(u);
      event.vpip.push_back(v);
      event.thetapip.push_back(theta);	    
      event.phipip.push_back(phi);
      event.ppipx[ipip]=mom.x();
      event.ppipy[ipip]=mom.y();
      event.ppipz[ipip]=mom.z();
      
      event.fpip.push_back(fldet);
      event.flpip[ipip]=fldet;
      
      PionPCont.push_back(mom);
      PionPXCont.push_back(primPos);
      ipip++;
    }
    if( pid==1 ){//Pi-
      event.ppim.push_back(p);
      event.upim.push_back(u);
      event.vpim.push_back(v);
      event.thetapim.push_back(theta);	    
      event.phipim.push_back(phi);
      event.ppimx[ipim]=mom.x();
      event.ppimy[ipim]=mom.y();
      event.ppimz[ipim]=mom.z();
      
      event.fpim.push_back(fldet);
      event.flpim[ipim]=fldet;
      
      PionMCont.push_back(mom);
      PionMXCont.push_back(primPos);
      ipim++;
    }
    if( pid==10 ){//Proton
      event.pp.push_back(p);
      event.up.push_back(u);
      event.vp.push_back(v);
      event.thetap.push_back(theta);	    
      event.phip.push_back(phi);
      event.ppx[ip]=mom.x();
      event.ppy[ip]=mom.y();
      event.ppz[ip]=mom.z();
      
      event.fp.push_back(fldet);
      event.flp[ip]=fldet;
      
      ProtonCont.push_back(mom);
      ProtonXCont.push_back(primPos);
      ip++;
    }
    if( pid==14 ){//Anti-Proton
      event.ppb.push_back(p);
      event.upb.push_back(u);
      event.vpb.push_back(v);
      event.thetapb.push_back(theta);	    
      event.phipb.push_back(phi);
      event.ppbx[ipb]=mom.x();
      event.ppby[ipb]=mom.y();
      event.ppbz[ipb]=mom.z();
      
      event.fpb.push_back(fldet);
      event.flpb[ipb]=fldet;
      
      ProtonbCont.push_back(mom);
      ProtonbXCont.push_back(primPos);
      ipb++;
    }
    if( pid==6 ){ //Mu+
      event.pmup.push_back(p);
      event.umup.push_back(u);
      event.vmup.push_back(v);
      event.thetamup.push_back(theta);	    
      event.phimup.push_back(phi);
      event.pmupx[imup]=mom.x();
      event.pmupy[imup]=mom.y();
      event.pmupz[imup]=mom.z();
      
      event.fmup.push_back(fldet);
      event.flmup[imup]=fldet;
      
      MuonPCont.push_back(mom);
      MuonPXCont.push_back(primPos);
      imup++;
    }
    if( pid==5 ){ //Mu-
      event.pmum.push_back(p);
      event.umum.push_back(u);
      event.vmum.push_back(v);
      event.thetamum.push_back(theta);	    
      event.phimum.push_back(phi);
      event.pmumx[imum]=mom.x();
      event.pmumy[imum]=mom.y();
      event.pmumz[imum]=mom.z();
      
      event.fmum.push_back(fldet);
      event.flmum[imum]=fldet;
      
      MuonMCont.push_back(mom);
      MuonMXCont.push_back(primPos);
      imum++;
    }
    if( pid==8 ){ //e+
      event.pep.push_back(p);
      event.uep.push_back(u);
      event.vep.push_back(v);
      event.thetaep.push_back(theta);	    
      event.phiep.push_back(phi);
      event.pepx[iep]=mom.x();
      event.pepy[iep]=mom.y();
      event.pepz[iep]=mom.z();
      
      event.fep.push_back(fldet);
      event.flep[iep]=fldet;
      
      PositronCont.push_back(mom);
      PositronXCont.push_back(primPos);
      iep++;
    }
    if( pid==7 ){ //e-
      event.pem.push_back(p);
      event.uem.push_back(u);
      event.vem.push_back(v);
      event.thetaem.push_back(theta);	    
      event.phiem.push_back(phi);
      event.pemx[iem]=mom.x();
      event.pemy[iem]=mom.y();
      event.pemz[iem]=mom.z();
      
      event.fem.push_back(fldet);
      event.flem[iem]=fldet;
      
      ElectronCont.push_back(mom);
      ElectronXCont.push_back(primPos);
      iem++;
    }
  }
  
  //Number of Particles
  int nKp=KaonPCont.size();
  int nKm=KaonMCont.size();
  int nPip=PionPCont.size();  
  int nPim=PionMCont.size();
  int nP=ProtonCont.size();
  int nPb=ProtonbCont.size();
  int nMup=MuonPCont.size();  
  int nMum=MuonMCont.size();
  int nEp=PositronCont.size();  
  int nEm=ElectronCont.size();
  
  event.ntKp = nKp; 
  event.ntKm = nKm; 
  event.ntPip = nPip; 
  event.ntPim = nPim; 
  event.ntP = nP; 
  event.ntPb = nPb; 
  event.ntMup = nMup; 
  event.ntMum = nMum; 
  event.ntEp = nEp; 
  event.ntEm = nEm; 
  
  //Beam
  int ibe=0;
  int nB=BeamPCont.size();
  event.ntB = nB;
  for( int ib=0; ib<nB; ++ib ){
    ThreeVector pb=BeamPCont[ib];
    double ub=pb.x()/pb.z(), vb=pb.y()/pb.z();
    event.pB[ibe]=pb.mag();
    event.uB[ibe]=ub;
    event.vB[ibe]=vb;
    event.pBx[ibe]=pb.x();
    event.pBy[ibe]=pb.y();
    event.pBz[ibe]=pb.z();
    ibe++;
  }

  tree->Fill();

  delete part;

  return true;
}

void EventSpecTracking::InitializeEvent( void )
{
  //PriInfo
  event.priposx = -9999.0;
  event.priposy = -9999.0;
  event.priposz = -9999.0;
  event.pbeam = -9999.0;
  event.ubeam = -9999.0;
  event.vbeam = -9999.0;
  event.abmom = -9999.0;
  event.prim1 = -9999.0;
  event.prip1 = -9999.0;
  event.pritheta1 = -9999.0;
  event.priphi1 = -9999.0;
  event.prithetacm1 = -9999.0;
  event.priphicm1 = -9999.0;
  event.prim2 = -9999.0;
  event.prip2 = -9999.0;
  event.pritheta2 = -9999.0;
  event.priphi2 = -9999.0;
  event.prithetacm2 = -9999.0;
  event.priphicm2 = -9999.0;

  //Generated Beam
  event.gbnhits = -1;
  event.gbid.clear();
  event.gbtype.clear();
  event.gbp.clear();
  event.gbpx.clear();
  event.gbpy.clear();
  event.gbpz.clear();
  event.gbpid.clear();
  event.gbvx.clear();
  event.gbvy.clear();

  //Generated Scat
  event.gsnhits = -1;
  event.gsid.clear();
  event.gstype.clear();
  event.gsp.clear();
  event.gspx.clear();
  event.gspy.clear();
  event.gspz.clear();
  event.gspid.clear();
  event.gsvx.clear();
  event.gsvy.clear();

  //Beam Tracking
  event.idb.clear();
  event.typeb.clear();
  event.ntb = -1;
  event.layerb.clear();
  event.chisqrb.clear();
  event.x0b.clear();
  event.y0b.clear();
  event.u0b.clear();
  event.v0b.clear();
  event.posb.clear();
  event.resb.clear();

  //Scat Tracking
  event.ids.clear();
  event.types.clear();
  event.nts.clear();
  event.layers.clear();
  event.chisqrs.clear();
  event.p.clear();
  event.theta.clear();
  event.phi.clear();
  event.dp.clear();
  event.path.clear();
  event.m2.clear();
  event.x0s.clear();
  event.y0s.clear();
  event.u0s.clear();
  event.v0s.clear();
  event.poss.clear();
  event.ress.clear();

  //PID
  event.ntKp = -1;
  event.ntKm = -1;
  event.ntPip = -1;
  event.ntPim = -1;
  event.ntP = -1;
  event.ntPb = -1;
  event.ntMup = -1;
  event.ntMum = -1;
  event.ntEp = -1;
  event.ntEm = -1;

  event.pkp.clear();
  event.ukp.clear();
  event.vkp.clear();
  event.thetakp.clear();
  event.phikp.clear();

  event.pkm.clear();
  event.ukm.clear();
  event.vkm.clear();
  event.thetakm.clear();
  event.phikm.clear();
  
  event.ppip.clear();
  event.upip.clear();
  event.vpip.clear();
  event.thetapip.clear();
  event.phipip.clear();

  event.ppim.clear();
  event.upim.clear();
  event.vpim.clear();
  event.thetapim.clear();
  event.phipim.clear();

  event.pp.clear();
  event.up.clear();
  event.vp.clear();
  event.thetap.clear();
  event.phip.clear();

  event.ppb.clear();
  event.upb.clear();
  event.vpb.clear();
  event.thetapb.clear();
  event.phipb.clear();

  event.pmup.clear();
  event.umup.clear();
  event.vmup.clear();
  event.thetamup.clear();
  event.phimup.clear();

  event.pmum.clear();
  event.umum.clear();
  event.vmum.clear();
  event.thetamum.clear();
  event.phimum.clear();
  
  event.pep.clear();
  event.uep.clear();
  event.vep.clear();
  event.thetaep.clear();
  event.phiep.clear();

  event.pem.clear();
  event.uem.clear();
  event.vem.clear();
  event.thetaem.clear();
  event.phiem.clear();

  //PID flag
  event.ntdet.clear();
  event.ntrich = -1;
  event.ntpid1 = -1;
  event.ntpid2 = -1;
  event.nttof = -1;
  event.ntitof = -1;
  event.ntpad = -1;
  event.ntmf = -1;

  event.fkp.clear();
  event.fkm.clear();
  event.fpip.clear();
  event.fpim.clear();
  event.fp.clear();
  event.fpb.clear();
  event.fmup.clear();
  event.fmum.clear();
  event.fep.clear();
  event.fem.clear();

  //Output
  for( int it=0; it<MaxHits; it++ ){
    event.pkpx[it]= -9999.0;
    event.pkpy[it]= -9999.0;
    event.pkpz[it]= -9999.0;

    event.pkmx[it]= -9999.0;
    event.pkmy[it]= -9999.0;
    event.pkmz[it]= -9999.0;

    event.ppipx[it]= -9999.0;
    event.ppipy[it]= -9999.0;
    event.ppipz[it]= -9999.0;

    event.ppimx[it]= -9999.0;
    event.ppimy[it]= -9999.0;
    event.ppimz[it]= -9999.0;

    event.ppx[it]= -9999.0;
    event.ppy[it]= -9999.0;
    event.ppz[it]= -9999.0;

    event.ppbx[it]= -9999.0;
    event.ppby[it]= -9999.0;
    event.ppbz[it]= -9999.0;

    event.pmupx[it]= -9999.0;
    event.pmupy[it]= -9999.0;
    event.pmupz[it]= -9999.0;

    event.pmumx[it]= -9999.0;
    event.pmumy[it]= -9999.0;
    event.pmumz[it]= -9999.0;

    event.pepx[it]= -9999.0;
    event.pepy[it]= -9999.0;
    event.pepz[it]= -9999.0;

    event.pemx[it]= -9999.0;
    event.pemy[it]= -9999.0;
    event.pemz[it]= -9999.0;
  }

  event.ntB = -1;
  for( int it=0; it<MaxHits; it++ ){
    event.pB[it]= -9999.0;
    event.uB[it]= -9999.0;
    event.vB[it]= -9999.0;
    event.pBx[it]= -9999.0;
    event.pBy[it]= -9999.0;
    event.pBz[it]= -9999.0;
  }

  for( int it=0; it<MaxHits; it++ ){
    event.flkp[it]= -1;
    event.flkm[it]= -1;
    event.flpip[it]= -1;
    event.flpim[it]= -1;
    event.flp[it]= -1;
    event.flpb[it]= -1;
    event.flmup[it]= -1;
    event.flmum[it]= -1;
    event.flep[it]= -1;
    event.flem[it]= -1;
  }
}


bool EventSpecTracking::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventSpecTracking;
}

bool ConfMan:: InitializeHistograms()
{  
  HBTree("tree","tree of Spec");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //PriInfo  
  tree->Branch("priposx", &event.priposx);
  tree->Branch("priposy", &event.priposy);
  tree->Branch("priposz", &event.priposz);
  tree->Branch("pbeam", &event.pbeam);
  tree->Branch("ubeam", &event.ubeam);
  tree->Branch("vbeam", &event.vbeam);
  tree->Branch("abmom", &event.abmom);
  tree->Branch("prim1", &event.prim1);
  tree->Branch("prip1", &event.prip1);
  tree->Branch("pritheta1", &event.pritheta1);
  tree->Branch("priphi1", &event.priphi1);
  tree->Branch("prithetacm1", &event.prithetacm1);
  tree->Branch("priphicm1", &event.priphicm1);
  tree->Branch("prim2", &event.prim2);
  tree->Branch("prip2", &event.prip2);
  tree->Branch("pritheta2", &event.pritheta2);
  tree->Branch("priphi2", &event.priphi2);
  tree->Branch("prithetacm2", &event.prithetacm2);
  tree->Branch("priphicm2", &event.priphicm2);

  //Generated Scat
  tree->Branch("gsnhits", &event.gsnhits);
  tree->Branch("gsid", &event.gsid);
  tree->Branch("gstype", &event.gstype);
  tree->Branch("gsp", &event.gsp);
  tree->Branch("gspx", &event.gspx);
  tree->Branch("gspy", &event.gspy);
  tree->Branch("gspz", &event.gspz);
  tree->Branch("gspid", &event.gspid);
  tree->Branch("gsvx", &event.gsvx);
  tree->Branch("gsvy", &event.gsvy);

  //Generated Beam
  tree->Branch("gbnhits", &event.gbnhits);
  tree->Branch("gbid", &event.gbid);
  tree->Branch("gbtype", &event.gbtype);
  tree->Branch("gbp", &event.gbp);
  tree->Branch("gbpx", &event.gbpx);
  tree->Branch("gbpy", &event.gbpy);
  tree->Branch("gbpz", &event.gbpz);
  tree->Branch("gbpid", &event.gbpid);
  tree->Branch("gbvx", &event.gbvx);
  tree->Branch("gbvy", &event.gbvy);

  //Beam Tracking
  tree->Branch("idb", &event.idb);
  tree->Branch("typeb", &event.typeb);
  tree->Branch("ntb", &event.ntb);
  tree->Branch("layerb", &event.layerb);
  tree->Branch("chisqrb", &event.chisqrb);
  tree->Branch("x0b", &event.x0b);
  tree->Branch("y0b", &event.y0b);
  tree->Branch("u0b", &event.u0b);
  tree->Branch("v0b", &event.v0b);
  tree->Branch("posb", &event.posb);
  tree->Branch("resb", &event.resb);

  //Scat Tracking
  tree->Branch("ids", &event.ids);
  tree->Branch("types", &event.types);
  tree->Branch("nts", &event.nts);
  tree->Branch("layers", &event.layers);
  tree->Branch("chisqrs", &event.chisqrs);
  tree->Branch("p", &event.p);
  tree->Branch("theta", &event.theta);
  tree->Branch("phi", &event.phi);
  tree->Branch("dp", &event.dp);
  tree->Branch("path", &event.path);
  tree->Branch("m2", &event.m2);
  tree->Branch("x0s", &event.x0s);
  tree->Branch("y0s", &event.y0s);
  tree->Branch("u0s", &event.u0s);
  tree->Branch("v0s", &event.v0s);
  tree->Branch("poss", &event.poss);
  tree->Branch("ress", &event.ress);

  //PID
  tree->Branch("ntKp", &event.ntKp, "ntKp/I");
  tree->Branch("ntKm", &event.ntKm, "ntKm/I");
  tree->Branch("ntPip", &event.ntPip, "ntPip/I");
  tree->Branch("ntPim", &event.ntPim, "ntPim/I");
  tree->Branch("ntP", &event.ntP, "ntP/I");
  tree->Branch("ntPb", &event.ntPb, "ntPb/I");
  tree->Branch("ntMup", &event.ntMup, "ntMup/I");
  tree->Branch("ntMum", &event.ntMum, "ntMum/I");
  tree->Branch("ntEp", &event.ntEp, "ntEp/I");
  tree->Branch("ntEm", &event.ntEm, "ntEm/I");

  tree->Branch("pkp", &event.pkp);
  tree->Branch("ukp", &event.ukp);
  tree->Branch("vkp", &event.vkp);
  tree->Branch("thetakp", &event.thetakp);
  tree->Branch("phikp", &event.phikp);

  tree->Branch("pkm", &event.pkm);
  tree->Branch("ukm", &event.ukm);
  tree->Branch("vkm", &event.vkm);
  tree->Branch("thetakm", &event.thetakm);
  tree->Branch("phikm", &event.phikm);

  tree->Branch("ppip", &event.ppip);
  tree->Branch("upip", &event.upip);
  tree->Branch("vpip", &event.vpip);
  tree->Branch("thetapip", &event.thetapip);
  tree->Branch("phipip", &event.phipip);

  tree->Branch("ppim", &event.ppim);
  tree->Branch("upim", &event.upim);
  tree->Branch("vpim", &event.vpim);
  tree->Branch("thetapim", &event.thetapim);
  tree->Branch("phipim", &event.phipim);

  tree->Branch("pp", &event.pp);
  tree->Branch("up", &event.up);
  tree->Branch("vp", &event.vp);
  tree->Branch("thetap", &event.thetap);
  tree->Branch("phip", &event.phip);

  tree->Branch("ppb", &event.ppb);
  tree->Branch("upb", &event.upb);
  tree->Branch("vpb", &event.vpb);
  tree->Branch("thetapb", &event.thetapb);
  tree->Branch("phipb", &event.phipb);

  tree->Branch("pmup", &event.pmup);
  tree->Branch("umup", &event.umup);
  tree->Branch("vmup", &event.vmup);
  tree->Branch("thetamup", &event.thetamup);
  tree->Branch("phimup", &event.phimup);

  tree->Branch("pmum", &event.pmum);
  tree->Branch("umum", &event.umum);
  tree->Branch("vmum", &event.vmum);
  tree->Branch("thetamum", &event.thetamum);
  tree->Branch("phimum", &event.phimum);

  tree->Branch("pep", &event.pep);
  tree->Branch("uep", &event.uep);
  tree->Branch("vep", &event.vep);
  tree->Branch("thetaep", &event.thetaep);
  tree->Branch("phiep", &event.phiep);

  tree->Branch("pem", &event.pem);
  tree->Branch("uem", &event.uem);
  tree->Branch("vem", &event.vem);
  tree->Branch("thetaem", &event.thetaem);
  tree->Branch("phiem", &event.phiem);

  //PID flag
  tree->Branch("ntdet", &event.ntdet);
  tree->Branch("ntrich", &event.ntrich, "ntrich/I");
  tree->Branch("ntpid1", &event.ntpid1, "ntpid1/I");
  tree->Branch("ntpid2", &event.ntpid2, "ntpid2/I");
  tree->Branch("nttof", &event.nttof, "nttof/I");
  tree->Branch("ntitof", &event.ntitof, "ntitof/I");
  tree->Branch("ntpad", &event.ntpad, "ntpad/I");
  tree->Branch("ntmf", &event.ntmf, "ntmf/I");

  tree->Branch("fkp", &event.fkp);
  tree->Branch("fkm", &event.fkm);
  tree->Branch("fpip", &event.fpip);
  tree->Branch("fpim", &event.fpim);
  tree->Branch("fp", &event.fp);
  tree->Branch("fpb", &event.fpb);
  tree->Branch("fmup", &event.fmup);
  tree->Branch("fmum", &event.fmum);
  tree->Branch("fep", &event.fep);
  tree->Branch("fem", &event.fem);

  //Output
  tree->Branch("pkpx", event.pkpx, "pkpx[ntKp]/D");
  tree->Branch("pkpy", event.pkpy, "pkpy[ntKp]/D");
  tree->Branch("pkpz", event.pkpz, "pkpz[ntKp]/D");

  tree->Branch("pkmx", event.pkmx, "pkmx[ntKm]/D");
  tree->Branch("pkmy", event.pkmy, "pkmy[ntKm]/D");
  tree->Branch("pkmz", event.pkmz, "pkmz[ntKm]/D");

  tree->Branch("ppipx", event.ppipx, "ppipx[ntPip]/D");
  tree->Branch("ppipy", event.ppipy, "ppipy[ntPip]/D");
  tree->Branch("ppipz", event.ppipz, "ppipz[ntPip]/D");

  tree->Branch("ppimx", event.ppimx, "ppimx[ntPim]/D");
  tree->Branch("ppimy", event.ppimy, "ppimy[ntPim]/D");
  tree->Branch("ppimz", event.ppimz, "ppimz[ntPim]/D");

  tree->Branch("ppx", event.ppx, "ppx[ntP]/D");
  tree->Branch("ppy", event.ppy, "ppy[ntP]/D");
  tree->Branch("ppz", event.ppz, "ppz[ntP]/D");

  tree->Branch("ppbx", event.ppbx, "ppbx[ntPb]/D");
  tree->Branch("ppby", event.ppby, "ppby[ntPb]/D");
  tree->Branch("ppbz", event.ppbz, "ppbz[ntPb]/D");

  tree->Branch("pmupx", event.pmupx, "pmupx[ntMup]/D");
  tree->Branch("pmupy", event.pmupy, "pmupy[ntMup]/D");
  tree->Branch("pmupz", event.pmupz, "pmupz[ntMup]/D");

  tree->Branch("pmumx", event.pmumx, "pmumx[ntMum]/D");
  tree->Branch("pmumy", event.pmumy, "pmumy[ntMum]/D");
  tree->Branch("pmumz", event.pmumz, "pmumz[ntMum]/D");

  tree->Branch("pepx", event.pepx, "pepx[ntEp]/D");
  tree->Branch("pepy", event.pepy, "pepy[ntEp]/D");
  tree->Branch("pepz", event.pepz, "pepz[ntEp]/D");

  tree->Branch("pemx", event.pemx, "pemx[ntEm]/D");
  tree->Branch("pemy", event.pemy, "pemy[ntEm]/D");
  tree->Branch("pemz", event.pemz, "pemz[ntEm]/D");

  tree->Branch("ntB", &event.ntB, "ntB/I");
  tree->Branch("pB", event.pB, "pB[ntB]/D");
  tree->Branch("uB", event.uB, "uB[ntB]/D");
  tree->Branch("vB", event.vB, "vB[ntB]/D");
  tree->Branch("pBx", event.pBx, "pBx[ntB]/D");
  tree->Branch("pBy", event.pBy, "pBy[ntB]/D");
  tree->Branch("pBz", event.pBz, "pBz[ntB]/D");

  tree->Branch("flkp", event.flkp, "flkp[ntKp]/I");
  tree->Branch("flkm", event.flkm, "flkm[ntKm]/I");
  tree->Branch("flpip", event.flpip, "flpip[ntPip]/I");
  tree->Branch("flpim", event.flpim, "flpim[ntPim]/I");
  tree->Branch("flp", event.flp, "flp[ntP]/I");
  tree->Branch("flpb", event.flpb, "flp[ntPb]/I");
  tree->Branch("flmup", event.flmup, "flmup[ntMup]/I");
  tree->Branch("flmum", event.flmum, "flmum[ntMum]/I");
  tree->Branch("flep", event.flep, "flep[ntEp]/I");
  tree->Branch("flem", event.flem, "flem[ntEm]/I");

  return true;
}
