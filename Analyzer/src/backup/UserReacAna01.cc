/*
  UserReacAna01.cc

  2024/05 K.Shirotori 
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>

#include "ConfMan.hh"
#include "RootHelper.hh"
#include "DetectorInfo.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "TrRawHit.hh"
#include "HodoRawHit.hh"

#include "Hodo2Hit.hh"
#include "Hodo1Hit.hh"
#include "HodoCluster.hh"
#include "HodoAnalyzer.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"

#include "DCRawHit.hh"
#include "DCAnalyzer.hh"
#include "DCHit.hh"
#include "GeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "DCDriftParamMan.hh"
#include "DCLocalTrack.hh"
#include "TrackHit.hh"
#include "DCLTrackHit.hh"

#include "ThreeVector.hh"
#include "Kinematics.hh"

#include "TFile.h"
#include "TTree.h"

#define check 0

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventReacAna01 : public VEvent
{
public:
   EventReacAna01();
   ~EventReacAna01();
   
   bool ProcessingBegin();
   bool ProcessingEnd();
   bool ProcessingNormal( TFile*, int evnum );
   void InitializeEvent();
   
private:
   RawData *rawData;
   HodoAnalyzer *hodoAna;
   DCAnalyzer *DCAna;
};

EventReacAna01::EventReacAna01()
   : VEvent(),
     rawData(0),
     hodoAna(new HodoAnalyzer()),
     DCAna(new DCAnalyzer())
{
}

EventReacAna01::~EventReacAna01()
{
   if (hodoAna) delete hodoAna; 
   if (DCAna)   delete DCAna;
   if (rawData) delete rawData;
}

#ifndef MaxHits 
#define MaxHits 32
#endif

struct Event{
   //////Timing counter
   double de0;

   int nhits_beam_tof;
   int nhits_scat_tof;

   double beam_tof[MaxHits];
   double scat_tof[MaxHits];

   int seg_ltof[MaxHits];
   double x_ltof[MaxHits];
   double y_ltof[MaxHits];
      
   //Local tracking
   int    nt;
   std::vector<int> layer;
   double chisqr[MaxHits];
   double x0[MaxHits], u0[MaxHits];
   double y0[MaxHits], v0[MaxHits];
   double x1[MaxHits], u1[MaxHits];
   double y1[MaxHits], v1[MaxHits];

   //Reaction events (2 tracks analysis)
   double theta1, phi1;
   double theta2, phi2;

   double pathl1_vt, pathl2_vt;
   double pathl1_0, pathl2_0;

   double dir1_x, dir1_y, dir1_z;
   double dir2_x, dir2_y, dir2_z;

   double vtx, vty, vtz;

   double cld, op_angle;

   double dy1_t1, dy2_t1;
   double dy1_t2, dy2_t2;

   double time1, time2;
   double de1, de2;

   double beta1, beta2;
};
static Event event;

bool EventReacAna01::ProcessingBegin()
{
 return true;
}

////////////////////////////
////////////////////////////
//Event selection parameters
const double beam_cut_low  =  -10.0;
const double beam_cut_high =   10.0;
const double scat_cut_low  = -20.0;
const double scat_cut_high = 100.0;

const double tracking_chisqr = 30.0;

const double y_ltof_cut = 105.0;

////////////////////////////
//Constants
const double z_ltof = 2195.;

const double LightVel = 299.792458;
const double time_ofs = z_ltof/LightVel;

bool EventReacAna01::ProcessingNormal( TFile* iFile, int evnum )
{
   const std::string funcname = "ProcessingNormal";
   
   ConfMan *confMan = ConfMan::GetConfManager();
   
   TTree *tree = (TTree*)iFile->Get("tree");
   
   rawData = new RawData;
   
   int n_entry = tree->GetEntries();
   for(int i_entry = 0; i_entry < n_entry; i_entry++){
      
      if( i_entry%10000 == 0 ){
         std::cout << "EventNum = " << i_entry << std::endl;
      }
      if( i_entry == evnum ){
         std::cout << "# of analyzed events: " << std::dec << i_entry << std::endl;
         break;
      }

      if( !rawData->DecodeRawHits(iFile, i_entry) ) return false;
      
      TTree *tree = static_cast<TTree *>(gFile->Get("tree"));
      
      InitializeEvent();
      
      HF1( 1, 0. );
      
      //**************************************************************************
      //******************Timing Counter
#if check
      std::cout << "Timing Counter***********************" << std::endl;
#endif
      hodoAna->DecodeUTOFHits( rawData );
      hodoAna->DecodeLTOFHits( rawData );
      hodoAna->DecodeT1Hits( rawData );
      
      bool flag_beam_tof = false;
      bool flag_scat_tof = false;
      double MT_UTOF = -9999.0;
      
      // T1-UTOF
      {
         int nh1=hodoAna->GetNHitsUTOF();
         int nh2=hodoAna->GetNHitsT1();
         int nhits_b=0;
         for( int i1=0; i1<nh1; ++i1 ){
            Hodo2Hit *utof=hodoAna->GetHitUTOF(i1);
            double mt_utof=utof->MeanTime();
            double de_utof=utof->DeltaE();
            MT_UTOF = mt_utof;
            event.de0 = de_utof;
            //std::cout << nh1 << " : " << mt_utof << std::endl;
            if(!utof) continue;
            for(int i2=0; i2<nh2; ++i2 ){
               Hodo2Hit *t1=hodoAna->GetHitT1(i2);
               if(!t1) continue;
               double mt_t1=t1->MeanTime();
               
               event.beam_tof[nhits_b] = mt_t1 - mt_utof;
               nhits_b++;
               
               if(beam_cut_low<(mt_t1-mt_utof)&&(mt_t1-mt_utof)<beam_cut_high)
                  flag_beam_tof = true;
            }
         }
         event.nhits_beam_tof = nhits_b;
      }
      // LTOF-UTOF
      {
         int nh1=hodoAna->GetNHitsUTOF();
         int nh2=hodoAna->GetNHitsLTOF();
         int nhits_s=0;
         for( int i1=0; i1<nh1; ++i1 ){
            Hodo2Hit *utof=hodoAna->GetHitUTOF(i1);
            double mt_utof=utof->MeanTime();
            if(!utof) continue;
            for(int i2=0; i2<nh2; ++i2 ){
               Hodo2Hit *ltof=hodoAna->GetHitLTOF(i2);
               if(!ltof) continue;
               double mt_ltof=ltof->MeanTime();
               
               event.scat_tof[nhits_s] = mt_ltof - mt_utof;
               nhits_s++;
               
               if(scat_cut_low<(mt_ltof-mt_utof)&&(mt_ltof-mt_utof)<scat_cut_high)
                  flag_scat_tof = true;
            }
         }
         event.nhits_scat_tof = nhits_s;
      }   
      
      {
         int nh=hodoAna->GetNHitsLTOF();
         for( int i=0; i<nh; ++i ){
            Hodo2Hit *hit=hodoAna->GetHitLTOF(i);
            int seg=hit->SegmentId();
            event.seg_ltof[i] = seg;
            
            double tl = hit->GetTLeft();
            double tr = hit->GetTRight();
            double tdif = tr-tl;
            double x=(tdif)*40.;
            double y=(3.5-seg)*200.;
            event.x_ltof[i] = x;
            event.y_ltof[i] = y;
         }
      }

      if( !(flag_beam_tof&&flag_scat_tof) ) continue;
      HF1( 1, 1. );
      
      //**************************************************************************
      //KLDC Tracking
#if check
      std::cout << "KLDC Tracking***********************" << std::endl;
#endif
      DCAna->DecodeKLDCRawHits( rawData );  
      DCAna->TrackSearchKLDC();
      bool flag_track_chisqr = false;
      {
         int nt=DCAna->GetNtracksKLDC();
         event.nt=nt;
         
         for( int it=0; it<nt; ++it ){
            DCLocalTrack *tp=DCAna->GetTrackKLDC(it);
            int nh=tp->GetNHit();
            double chisqr=tp->GetChiSquare();
            double x0=tp->GetX0();
            double u0=tp->GetU0();
            double y0=tp->GetY0();
            double v0=tp->GetV0();
            
            double xtgt=tp->GetX( 0.);
            double utgt=u0;
            double ytgt=tp->GetY( 0.);
            double vtgt=v0;
            
            double xltof=tp->GetX( z_ltof );
            double ultof=u0;
            double yltof=tp->GetY( z_ltof );
            double vltof=v0;
            
            event.chisqr[it] = chisqr;
            event.x0[it] = xtgt;
            event.u0[it] = utgt;
            event.y0[it] = ytgt;
            event.v0[it] = vtgt;
            
            event.x1[it] = xltof;
            event.u1[it] = ultof;
            event.y1[it] = yltof;
            event.v1[it] = vltof;
            
            if( chisqr < tracking_chisqr) flag_track_chisqr=true;
         }
      }
      
      if( !flag_track_chisqr ) continue;
      HF1( 1, 2. );
      
      //**************************************************************************
      //2 Track Analysis
#if check
      std::cout << "2 Track Analysis***********************" << std::endl;
#endif
      bool flag_tracks = false;
      bool flag_times = false;
      {
         int nt=DCAna->GetNtracksKLDC();
         if(nt==2){
            flag_tracks = true;

            ThreeVector pos0_t1, pos1_t1;
            ThreeVector pos0_t2, pos1_t2;
            
            ThreeVector dir0_t1, dir1_t1;
            ThreeVector dir0_t2, dir1_t2;
            
            DCLocalTrack *tp1=DCAna->GetTrackKLDC(0);
            DCLocalTrack *tp2=DCAna->GetTrackKLDC(1);
            
            double x0t1 = tp1->GetX0();
            double y0t1 = tp1->GetY0();
            double u0t1 = tp1->GetU0();
            double v0t1 = tp1->GetV0();
            
            double x0t2 = tp2->GetX0();
            double y0t2 = tp2->GetY0();
            double u0t2 = tp2->GetU0();
            double v0t2 = tp2->GetV0();
            
            double x1t1 = tp1->GetX(z_ltof);
            double y1t1 = tp1->GetY(z_ltof);
            double u1t1 = tp1->GetU0();
            double v1t1 = tp1->GetV0();
            
            double x1t2 = tp2->GetX(z_ltof);
            double y1t2 = tp2->GetY(z_ltof);
            double u1t2 = tp2->GetU0();
            double v1t2 = tp2->GetV0();
            
            //Track1
            pos0_t1.setX(x0t1);
            pos0_t1.setY(y0t1);
            pos0_t1.setZ(0.0);
            
            dir0_t1.setX(u0t1);
            dir0_t1.setY(v0t1);
            dir0_t1.setZ(1.0/pow(1.0+u0t1*u0t1+v0t1*v0t1,0.5));
            
            pos1_t1.setX(x1t1);
            pos1_t1.setY(y1t1);
            pos1_t1.setZ(z_ltof);
            
            dir1_t1.setX(u1t1);
            dir1_t1.setY(v1t1);
            dir1_t1.setZ(1.0/pow(1.0+u1t1*u1t1+v1t1*v1t1,0.5));
            
            //Track2
            pos0_t2.setX(x0t2);
            pos0_t2.setY(y0t2);
            pos0_t2.setZ(0.0);
            
            dir0_t2.setX(u0t2);
            dir0_t2.setY(v0t2);
            dir0_t2.setZ(1.0/pow(1.0+u0t2*u0t2+v0t2*v0t2,0.5));
            
            pos1_t2.setX(x1t2);
            pos1_t2.setY(y1t2);
            pos1_t2.setZ(z_ltof);
            
            dir1_t2.setX(u1t2);
            dir1_t2.setY(v1t2);
            dir1_t2.setZ(1.0/pow(1.0+u1t2*u1t2+v1t2*v1t2,0.5));
            
            event.theta1=dir0_t1.theta()*Rad2Deg;
            event.theta2=dir0_t2.theta()*Rad2Deg;
            event.phi1=dir0_t1.phi()*Rad2Deg;
            event.phi2=dir0_t2.phi()*Rad2Deg;
            
            event.dir1_x=dir0_t1.x();
            event.dir1_y=dir0_t1.y();
            event.dir1_z=dir0_t1.z();
            
            event.dir2_x=dir0_t2.x();
            event.dir2_y=dir0_t2.y();
            event.dir2_z=dir0_t2.z();
            
            // std::cout<< "Vectors = "
            //          << dir1_0.mag() << "|"
            //          << dir2_0.mag() <<std::endl; 
            
            ThreeVector vert=VertexPoint( pos0_t1, pos0_t2, dir0_t1, dir0_t2 );
            double closedist=closeDist( pos0_t1, pos0_t2, dir0_t1, dir0_t2 );
            
            double path1 = (pos1_t1-pos0_t1).mag();
            double path2 = (pos1_t2-pos0_t2).mag();
            
            event.pathl1_0 = path1;
            event.pathl2_0 = path2;
            
            event.pathl1_vt = (pos1_t1-vert).mag();
            event.pathl2_vt = (pos1_t2-vert).mag();
            
            event.vtx=vert.x();
            event.vty=vert.y();
            event.vtz=vert.z();
            event.cld=closedist;
            
            event.op_angle=acos((dir0_t1*dir0_t2)/(dir0_t1.mag()*dir0_t2.mag()))*Rad2Deg;
            
            int nh=hodoAna->GetNHitsLTOF();
            if(nh==2){   
               Hodo2Hit *hit1=hodoAna->GetHitLTOF(0);
               Hodo2Hit *hit2=hodoAna->GetHitLTOF(1);
               
               int seg1=hit1->SegmentId();
               int seg2=hit2->SegmentId();
               
               // std::cout<< "Segs = "
               //          << seg1 << "|"
               //          << seg2 <<std::endl; 
               
               double y1=(3.5-seg1)*200.;
               double y2=(3.5-seg2)*200.;
               
               double dy1_t1 = fabs(y1-pos1_t1.y());
               double dy2_t1 = fabs(y2-pos1_t1.y());
               
               double dy1_t2 = fabs(y1-pos1_t2.y());
               double dy2_t2 = fabs(y2-pos1_t2.y());
               
               event.dy1_t1 = dy1_t1;
               event.dy2_t1 = dy2_t1;
               
               event.dy1_t2 = dy1_t2;
               event.dy2_t2 = dy2_t2;
               
               bool f_dy1_t1 = false;
               bool f_dy2_t1 = false;
               
               bool f_dy1_t2 = false;
               bool f_dy2_t2 = false;
               
               int counts=0;
               if( dy1_t1<y_ltof_cut && y_ltof_cut<dy1_t2 ){
                  double time = hit1->MeanTime()-MT_UTOF;
                  double de = hit1->DeltaE();
                  
                  event.time1 = time;
                  event.de1 = de;
                  event.beta1 = (path1/(time+time_ofs))/LightVel;
                  
                  f_dy1_t1 = true;
                  counts++;
               }
               if( dy2_t1<y_ltof_cut && y_ltof_cut<dy2_t2 ){
                  double time = hit2->MeanTime()-MT_UTOF;
                  double de = hit2->DeltaE();
                  
                  event.time1 = time;
                  event.de1 = de;
                  event.beta1 = (path1/(time+time_ofs))/LightVel;
                  
                  f_dy2_t1 = true;
                  counts++;
               }
               if( dy1_t2<y_ltof_cut && y_ltof_cut<dy1_t1 ){
                  double time= hit1->MeanTime()-MT_UTOF;
                  double de = hit1->DeltaE();
                  
                  event.time2 = time;
                  event.de2 = de;
                  event.beta2 = (path2/(time+time_ofs))/LightVel;
                  
                  f_dy1_t2 = true;
                  counts++;
               }
               if( dy2_t2<y_ltof_cut && y_ltof_cut<dy2_t1 ){
                  double time = hit2->MeanTime()-MT_UTOF;
                  double de = hit2->DeltaE();
                  
                  event.time2 = time;
                  event.de2 = de;
                  event.beta2 = (path2/(time+time_ofs))/LightVel;
                  
                  f_dy2_t2 = true;
                  counts++;
               }
               if(counts>=2) flag_times = true;

               // std::cout<< "Flags = "
               //          << f_dy1_t1 << ":"
               //          << f_dy2_t1 << "|"
               //          << f_dy1_t2 << ":"
               //          << f_dy2_t2 << " : " << counts <<std::endl; 
            }
         }
      }       

      if( !(flag_tracks && flag_times) ) continue;
      HF1( 1, 3. );

      tree->Fill();
   }
   
   return true;
}

void EventReacAna01::InitializeEvent( void )
{
   //////Timing counter
   event.de0  = -9999.0;
   
   event.nhits_beam_tof = -1;
   event.nhits_scat_tof = -1;

   for( int i=0; i<MaxHits; i++){
      event.beam_tof[i]   = -9999.0;
      event.scat_tof[i]   = -9999.0;
   }

   for( int i=0; i<MaxHits; i++){
      event.seg_ltof[i] = -9999.0;
      event.x_ltof[i] = -9999.0;
      event.y_ltof[i] = -9999.0;
   }
   
   //Local Tracking
   event.nt = -1;
   event.layer.clear();

   for( int i=0; i<MaxHits; i++){
      event.chisqr[i] = -9999.0;
      event.x0[i] = -9999.0;
      event.u0[i] = -9999.0;
      event.y0[i] = -9999.0;
      event.v0[i] = -9999.0;
      event.x1[i] = -9999.0;
      event.u1[i] = -9999.0;
      event.y1[i] = -9999.0;
      event.v1[i] = -9999.0;
   }

   //Reaction events (2 tracks analysis)
   event.theta1 = -9999.0;
   event.theta2 = -9999.0;
   event.phi1 = -9999.0;
   event.phi2 = -9999.0;
   
   event.pathl1_vt = -9999.0;
   event.pathl2_vt = -9999.0;
   event.pathl1_0 = -9999.0;
   event.pathl2_0 = -9999.0;
   
   event.dir1_x = -9999.0;
   event.dir1_y = -9999.0;
   event.dir1_z = -9999.0;
   
   event.dir2_x = -9999.0;
   event.dir2_y = -9999.0;
   event.dir2_z = -9999.0;
   
   event.vtx = -9999.0;
   event.vty = -9999.0;
   event.vtz = -9999.0;
   
   event.cld = -9999.0;
   event.op_angle = -9999.0;
   
   event.dy1_t1 = -9999.0;
   event.dy2_t1 = -9999.0;
   event.dy1_t2 = -9999.0;
   event.dy2_t2 = -9999.0;
   
   event.time1 = -9999.0;
   event.time2 = -9999.0;

   event.de1 = -9999.0;
   event.de2 = -9999.0;

   event.beta1 = -9999.0;
   event.beta2 = -9999.0;
}

bool EventReacAna01::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();

  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventReacAna01;
}

bool ConfMan:: InitializeHistograms()
{
   HB1( 1, "Status", 60, 0., 60. );
   
   HBTree("tree","tree");
   TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
   
   //////Timing counter
   tree->Branch("de0", &event.de0, "de0/D");  

   tree->Branch("nhits_beam_tof", &event.nhits_beam_tof, "nhits_beam_tof/I");
   tree->Branch("beam_tof",        event.beam_tof,       "beam_tof[32]/D");  
   
   tree->Branch("nhits_scat_tof", &event.nhits_scat_tof, "nhits_scat_tof/I");
   tree->Branch("scat_tof",        event.scat_tof,       "scat_tof[32]/D");  

   tree->Branch("seg_ltof", event.seg_ltof, "seg_ltof[32]/I");  
   tree->Branch("x_ltof", event.x_ltof, "x_ltof[32]/D");  
   tree->Branch("y_ltof", event.y_ltof, "y_ltof[32]/D");  

   //Tracking
   tree->Branch("nt", &event.nt, "nt/I");
   tree->Branch("layer", &event.layer);

   tree->Branch("chisqr", event.chisqr, "chisqr[32]/D");
   tree->Branch("x0", &event.x0, "x0[32]/D");
   tree->Branch("u0", &event.u0, "u0[32]/D");
   tree->Branch("y0", &event.y0, "y0[32]/D");
   tree->Branch("v0", &event.v0, "v0[32]/D");
   tree->Branch("x1", &event.x1, "x1[32]/D");
   tree->Branch("u1", &event.u1, "u1[32]/D");
   tree->Branch("y1", &event.y1, "y1[32]/D");
   tree->Branch("v1", &event.v1, "v1[32]/D");

   //Reaction events (2 tracks analysis)
   tree->Branch("theta1", &event.theta1, "theta1/D");
   tree->Branch("theta2", &event.theta2, "theta2/D");
   tree->Branch("phi1", &event.phi1, "phi1/D");
   tree->Branch("phi2", &event.phi2, "phi2/D");

   tree->Branch("pathl1_vt", &event.pathl1_vt, "pathl1_vt/D");
   tree->Branch("pathl2_vt", &event.pathl2_vt, "pathl2_vt/D");
   tree->Branch("pathl1_0", &event.pathl1_0, "pathl1_0/D");
   tree->Branch("pathl2_0", &event.pathl2_0, "pathl2_0/D");

   tree->Branch("dir1_x", &event.dir1_x, "dir1_x/D");
   tree->Branch("dir1_y", &event.dir1_y, "dir1_y/D");
   tree->Branch("dir1_z", &event.dir1_z, "dir1_z/D");

   tree->Branch("dir2_x", &event.dir2_x, "dir2_x/D");
   tree->Branch("dir2_y", &event.dir2_y, "dir2_y/D");
   tree->Branch("dir2_z", &event.dir2_z, "dir2_z/D");

   tree->Branch("vtx", &event.vtx, "vtx/D");
   tree->Branch("vty", &event.vty, "vty/D");
   tree->Branch("vtz", &event.vtz, "vtz/D");

   tree->Branch("cld", &event.cld, "cld/D");
   tree->Branch("op_angle", &event.op_angle, "op_angle/D");

   tree->Branch("dy1_t1", &event.dy1_t1, "dy1_t1/D");
   tree->Branch("dy2_t1", &event.dy2_t1, "dy2_t1/D");
   tree->Branch("dy1_t2", &event.dy1_t2, "dy1_t2/D");
   tree->Branch("dy2_t2", &event.dy2_t2, "dy2_t2/D");

   tree->Branch("time1", &event.time1, "time1/D");
   tree->Branch("time2", &event.time2, "time2/D");

   tree->Branch("de1", &event.de1, "de1/D");
   tree->Branch("de2", &event.de2, "de2/D");

   tree->Branch("beta1", &event.beta1, "beta1/D");
   tree->Branch("beta2", &event.beta2, "beta2/D");

   return true;
}
