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

struct Event{
   //////Timing counter
   std::vector<double> de0;

   std::vector<double> beam_tof;
   std::vector<double> scat_tof;

   std::vector<int> seg_ltof, seg_ltof1;
   std::vector<double> x_ltof, x_ltof1;
   std::vector<double> y_ltof, y_ltof1;
      
   //Local tracking
   int    nt;
   std::vector<int> layer;
   std::vector<double> chisqr;
   std::vector<double> x0, u0;
   std::vector<double> y0, v0;
   std::vector<double> x1, u1;
   std::vector<double> y1, v1;

   //Reaction events
   std::vector<double> theta, phi;
   std::vector<double> pathl, pathl_vt;

   std::vector<double> dir_x, dir_y, dir_z;

   std::vector<double> vtx, vty, vtz;

   std::vector<double> cld, op_angle;

   std::vector<double> dy;

   std::vector<double>  time, de;
   std::vector<double> ctime;
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

const double y_ltof_cut = 120.0;

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
      double time0 = -9999.0;
      
      // T1-UTOF
      {
         int nh1=hodoAna->GetNHitsUTOF();
         int nh2=hodoAna->GetNHitsT1();
         for( int i1=0; i1<nh1; ++i1 ){
            Hodo2Hit *utof=hodoAna->GetHitUTOF(i1);
            double mt_utof=utof->MeanTime();
            double de_utof=utof->DeltaE();
            time0 = mt_utof;
            event.de0.push_back(de_utof);
            //std::cout << nh1 << " : " << mt_utof << std::endl;
            if(!utof) continue;
            for(int i2=0; i2<nh2; ++i2 ){
               Hodo2Hit *t1=hodoAna->GetHitT1(i2);
               if(!t1) continue;
               double mt_t1=t1->MeanTime();
               
               event.beam_tof.push_back(mt_t1 - mt_utof);
               
               if(beam_cut_low<(mt_t1-mt_utof)&&(mt_t1-mt_utof)<beam_cut_high)
                  flag_beam_tof = true;
            }
         }
      }
      // LTOF-UTOF
      {
         int nh1=hodoAna->GetNHitsUTOF();
         int nh2=hodoAna->GetNHitsLTOF();
         for( int i1=0; i1<nh1; ++i1 ){
            Hodo2Hit *utof=hodoAna->GetHitUTOF(i1);
            double mt_utof=utof->MeanTime();
            if(!utof) continue;
            for(int i2=0; i2<nh2; ++i2 ){
               Hodo2Hit *ltof=hodoAna->GetHitLTOF(i2);
               if(!ltof) continue;
               double mt_ltof=ltof->MeanTime();
               
               event.scat_tof.push_back(mt_ltof - mt_utof);
               
               if(scat_cut_low<(mt_ltof-mt_utof)&&(mt_ltof-mt_utof)<scat_cut_high)
                  flag_scat_tof = true;
            }
         }
      }   
      
      {
         int nh=hodoAna->GetNHitsLTOF();
         for( int i=0; i<nh; ++i ){
            Hodo2Hit *hit=hodoAna->GetHitLTOF(i);
            int seg=hit->SegmentId();
            event.seg_ltof.push_back(seg);
            
            double tl = hit->GetTLeft();
            double tr = hit->GetTRight();
            double tdif = tr-tl;
            double x=(tdif)*40.;
            double y=(3.5-seg)*200.;
            event.x_ltof.push_back(x);
            event.y_ltof.push_back(y);
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
            
            event.chisqr.push_back(chisqr);
            event.x0.push_back(xtgt);
            event.u0.push_back(utgt);
            event.y0.push_back(ytgt);
            event.v0.push_back(vtgt);
            
            event.x1.push_back(xltof);
            event.u1.push_back(ultof);
            event.y1.push_back(yltof);
            event.v1.push_back(vltof);
            
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
      bool flag_track = false;
      bool flag_time = false;
      {
         std::vector<ThreeVector> pos0, dir0;
         std::vector<ThreeVector> pos1, dir1;
         int nt=DCAna->GetNtracksKLDC();
         
         for( int it=0; it<nt; ++it ){
            DCLocalTrack *tp=DCAna->GetTrackKLDC(it);
            
            double x0 = tp->GetX0();
            double y0 = tp->GetY0();
            double u0 = tp->GetU0();
            double v0 = tp->GetV0();
            
            double x1 = tp->GetX(z_ltof);
            double y1 = tp->GetY(z_ltof);
            double u1 = tp->GetU0();
            double v1 = tp->GetV0();

            ThreeVector pos_a, pos_b, dir;
            pos_a.setX(x0);
            pos_a.setY(y0);
            pos_a.setZ(0.0);
            pos0.push_back(pos_a);

            pos_b.setX(x1);
            pos_b.setY(y1);
            pos_b.setZ(z_ltof);
            pos1.push_back(pos_b);

            double w = 1.0/std::sqrt(1.0+u0*u0+v0*v0);
            dir.setX(u0*w);
            dir.setY(v0*w);
            dir.setZ(w);
            dir0.push_back(dir);

            // std::cout<< "Direction : "
            //          << dir.mag() << " <- " << dir
            //          <<std::endl; 
         }
         
         int nhvt=0;
         for( int i1=0; i1<nt; ++i1 ){
            for( int i2=i1+1; i2<nt; ++i2 ){
               ThreeVector vert=VertexPoint( pos0.at(i1), pos0.at(i2), dir0.at(i1), dir0.at(i2) );
               double closedist=closeDist( pos0.at(i1), pos0.at(i2), dir0.at(i1), dir0.at(i2) );

               // std::cout<< nt << " : *** = " << i1 << " : " << i2 << std::endl;
               // std::cout<< pos0.at(i1).x() << " : " << pos0.at(i2).x() << std::endl;
               event.vtx.push_back(vert.x());
               event.vty.push_back(vert.y());
               event.vtz.push_back(vert.z());
               event.cld.push_back(closedist);
               
               event.pathl_vt.push_back((pos1.at(i1)-vert).mag());
               
               event.op_angle.push_back( acos((dir0.at(i1)*dir0.at(i2))/(dir0.at(i1).mag()*dir0.at(i2).mag()))*Rad2Deg );
            
               nhvt++;
            }
         }
         if( nhvt>0 ) flag_track = true;

         int nh=hodoAna->GetNHitsLTOF();

         if( !(nh>1) ) continue;

         //std::cout << "***********" << std::endl;
         int nhtime=0;
         for( int i=0; i<nh; ++i ){
            Hodo2Hit *hit=hodoAna->GetHitLTOF(i);
            int seg=hit->SegmentId();
            double tl = hit->GetTLeft();
            double tr = hit->GetTRight();
            double tdif = tr-tl;
            double x=(tdif)*40.;
            double y=(3.5-seg)*200.;
            
            int nhy = pos1.size();
            for( int j=0; j<nhy; ++j ){
               double dy = fabs(y-pos1.at(j).y());
               event.dy.push_back(dy);

               //std::cout << dy << std::endl;
               
               if( dy<y_ltof_cut ){
                  double time = hit->MeanTime()-time0;
                  double de = hit->DeltaE();
                  
                  event.time.push_back(time);
                  event.de.push_back(de);

                  event.seg_ltof1.push_back(seg);
                  event.x_ltof1.push_back(x);
                  event.y_ltof1.push_back(y);
                  
                  //std::cout << time << std::endl;
                  //std::cout << i << " : " << j << std::endl;
                  
                  double ctime = time+2.3*de-2.39;
                  event.ctime.push_back(ctime);
                              
                  double theta = dir0.at(j).theta()*Rad2Deg;
                  double phi = dir0.at(j).phi()*Rad2Deg;
                  event.theta.push_back(theta);
                  event.phi.push_back(phi);
                  
                  double path = (pos1.at(j)-pos0.at(j)).mag();
                  event.pathl.push_back(path);
                  
                  event.dir_x.push_back(dir0.at(j).x());
                  event.dir_y.push_back(dir0.at(j).y());
                  event.dir_z.push_back(dir0.at(j).z());
                  
                  nhtime++;
               }
            }
            if( nhtime>0 ) flag_time = true;
         }
      }
      
      if( !(flag_track && flag_time) ) continue;
      HF1( 1, 3. );

      tree->Fill();
   }
   
   return true;
}

void EventReacAna01::InitializeEvent( void )
{
   //////Timing counter
   event.de0.clear();

   event.beam_tof.clear();
   event.scat_tof.clear();
   
   event.seg_ltof.clear();
   event.x_ltof.clear();
   event.y_ltof.clear();

   event.seg_ltof1.clear();
   event.x_ltof1.clear();
   event.y_ltof1.clear();
   
   //Local Tracking
   event.nt = -1;
   event.layer.clear();
   
   event.chisqr.clear();
   event.x0.clear();
   event.u0.clear();
   event.y0.clear();
   event.v0.clear();
   event.x1.clear();
   event.u1.clear();
   event.y1.clear();
   event.v1.clear();

   //Reaction events
   event.theta.clear();
   event.phi.clear();
   event.pathl.clear();
   event.pathl_vt.clear();
   
   event.dir_x.clear();
   event.dir_y.clear();
   event.dir_z.clear();
   
   event.vtx.clear();
   event.vty.clear();
   event.vtz.clear();
   
   event.cld.clear();
   event.op_angle.clear();
   
   event.dy.clear();
   
   event.time.clear();
   event.de.clear();
   
   event.ctime.clear();
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
   tree->Branch("de0", &event.de0 );  
   tree->Branch("beam_tof", &event.beam_tof );  
   tree->Branch("scat_tof", &event.scat_tof );  

   tree->Branch("seg_ltof", &event.seg_ltof );  
   tree->Branch("x_ltof", &event.x_ltof );  
   tree->Branch("y_ltof", &event.y_ltof );  

   tree->Branch("seg_ltof1", &event.seg_ltof1 );  
   tree->Branch("x_ltof1", &event.x_ltof1 );  
   tree->Branch("y_ltof1", &event.y_ltof1 );  

   //Tracking
   tree->Branch("nt", &event.nt, "nt/I");
   tree->Branch("layer", &event.layer);

   tree->Branch("chisqr", &event.chisqr );
   tree->Branch("x0", &event.x0 );
   tree->Branch("u0", &event.u0 );
   tree->Branch("y0", &event.y0 );
   tree->Branch("v0", &event.v0 );
   tree->Branch("x1", &event.x1 );
   tree->Branch("u1", &event.u1 );
   tree->Branch("y1", &event.y1 );
   tree->Branch("v1", &event.v1 );

   //Reaction events
   tree->Branch("theta", &event.theta );
   tree->Branch("phi", &event.phi );

   tree->Branch("pathl", &event.pathl );
   tree->Branch("pathl_vt", &event.pathl_vt );

   tree->Branch("dir_x", &event.dir_x );
   tree->Branch("dir_y", &event.dir_y );
   tree->Branch("dir_z", &event.dir_z );

   tree->Branch("vtx", &event.vtx );
   tree->Branch("vty", &event.vty );
   tree->Branch("vtz", &event.vtz );

   tree->Branch("cld", &event.cld );
   tree->Branch("op_angle", &event.op_angle );

   tree->Branch("dy", &event.dy );

   tree->Branch("time", &event.time );
   tree->Branch("de", &event.de );

   tree->Branch("ctime", &event.ctime );

   return true;
}
