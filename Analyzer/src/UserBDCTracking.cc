/*
  UserBDCTracking.cc

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

class EventBDCTracking : public VEvent
{
public:
   EventBDCTracking();
   ~EventBDCTracking();
   
   bool ProcessingBegin();
   bool ProcessingEnd();
   bool ProcessingNormal( TFile*, int evnum );
   void InitializeEvent();
   
private:
   RawData *rawData;
   HodoAnalyzer *hodoAna;
   DCAnalyzer *DCAna;
};

EventBDCTracking::EventBDCTracking()
   : VEvent(),
     rawData(0),
     hodoAna(new HodoAnalyzer()),
     DCAna(new DCAnalyzer())
{
}

EventBDCTracking::~EventBDCTracking()
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
   int nhits_utof;
   int id_utof[MaxHits];
   double ltdc_utof_l[NumOfSegUTOF];
   double ltdc_utof_r[NumOfSegUTOF];
   double tot_utof_l[NumOfSegUTOF];
   double tot_utof_r[NumOfSegUTOF];
   double mt_utof[NumOfSegUTOF];
   double de_utof[NumOfSegUTOF];
   
   int nhits_ltof;
   int id_ltof[MaxHits];
   double ltdc_ltof_l[NumOfSegLTOF];
   double ltdc_ltof_r[NumOfSegLTOF];
   double tot_ltof_l[NumOfSegLTOF];
   double tot_ltof_r[NumOfSegLTOF];
   double mt_ltof[NumOfSegLTOF];
   double de_ltof[NumOfSegLTOF];

   int nhits_t1;
   int id_t1[MaxHits];
   double ltdc_t1_l[NumOfSegT1];
   double ltdc_t1_r[NumOfSegT1];
   double tot_t1_l[NumOfSegT1];
   double tot_t1_r[NumOfSegT1];
   double mt_t1[NumOfSegT1];
   double de_t1[NumOfSegT1];

   int nhits_bht;
   int id_bht[MaxHits];
   double ltdc_bht[NumOfSegBref];
   double tot_bht[NumOfSegBref];

   int nhits_beam_tof;
   double beam_tof[MaxHits];
   
   int nhits_scat_tof;
   double scat_tof[MaxHits];
   int nhits_scat_tof2;
   double scat_tof2[MaxHits];
   int nhits_scat_tof3;
   double scat_tof3[MaxHits];
   
   //BDC
   std::vector<int> id_bdc_l1, id_bdc_l2, id_bdc_l3, id_bdc_l4;
   std::vector<int> id_bdc_l5, id_bdc_l6, id_bdc_l7, id_bdc_l8;
   std::vector<std::vector<double>> tot_bdc_l1, tot_bdc_l2, tot_bdc_l3, tot_bdc_l4;  
   std::vector<std::vector<double>> tot_bdc_l5, tot_bdc_l6, tot_bdc_l7, tot_bdc_l8;  
   std::vector<std::vector<double>> ltdc_bdc_l1, ltdc_bdc_l2, ltdc_bdc_l3, ltdc_bdc_l4;  
   std::vector<std::vector<double>> ltdc_bdc_l5, ltdc_bdc_l6, ltdc_bdc_l7, ltdc_bdc_l8;  

   //Local tracking
   int    nt;
   std::vector<int>    layer;
   std::vector<double> chisqr;
   std::vector<double> x0, u0;
   std::vector<double> y0, v0;
   std::vector<double> x1, u1;
   std::vector<double> y1, v1;

   std::vector<std::vector<double>> dt_l1, dt_l2, dt_l3, dt_l4;  
   std::vector<std::vector<double>> dt_l5, dt_l6, dt_l7, dt_l8;  
   std::vector<std::vector<double>> dl_l1, dl_l2, dl_l3, dl_l4;  
   std::vector<std::vector<double>> dl_l5, dl_l6, dl_l7, dl_l8;  
   std::vector<std::vector<double>> pos_l1, pos_l2, pos_l3, pos_l4;  
   std::vector<std::vector<double>> pos_l5, pos_l6, pos_l7, pos_l8;  
   std::vector<std::vector<double>> res_l1, res_l2, res_l3, res_l4;  
   std::vector<std::vector<double>> res_l5, res_l6, res_l7, res_l8;  
};
static Event event;

bool EventBDCTracking::ProcessingBegin()
{
 return true;
}

bool EventBDCTracking::ProcessingNormal( TFile* iFile, int evnum )
{
   const std::string funcname = "ProcessingNormal";

   ConfMan *confMan = ConfMan::GetConfManager();
    
   TTree *tree = (TTree*)iFile->Get("tree");
   
   rawData = new RawData;
   
   int n_entry = tree->GetEntries();
   for(int i_entry = 0; i_entry < n_entry; i_entry++){
      
      if( i_entry%1000 == 0 ){
         std::cout << "EventNum = " << i_entry << std::endl;
      }
      if( i_entry == evnum ){
         std::cout << "# of analyzed events: " << std::dec << i_entry << std::endl;
         break;
      }
       if( !rawData->DecodeRawHits(iFile, i_entry) ) return false;
       
       TTree *tree = static_cast<TTree *>(gFile->Get("tree"));
       
       InitializeEvent();
       
       //**************************************************************************
       //******************NormalizedData
       bool f_dtof_veto = false;
         
       //UTOF
       hodoAna->DecodeUTOFHits( rawData );
       {
          int nh=hodoAna->GetNHitsUTOF();
          int nh2=0;
          for( int i=0; i<nh; ++i ){
             Hodo2Hit *hit=hodoAna->GetHitUTOF(i);
             
             int seg=hit->SegmentId();
             event.id_utof[i] = seg;
             event.ltdc_utof_l[seg-1] = hit->GetTLeft();
             event.ltdc_utof_r[seg-1] = hit->GetTRight();
             event.tot_utof_l[seg-1] = hit->GetALeft();
             event.tot_utof_r[seg-1] = hit->GetARight();
             event.mt_utof[seg-1] = hit->MeanTime();
             event.de_utof[seg-1] = hit->DeltaE();
             nh2++;
          }
          event.nhits_utof = nh2;
       }
       //LTOF
       hodoAna->DecodeLTOFHits( rawData );
       {
          int nh=hodoAna->GetNHitsLTOF();
          int nh2=0;
          for( int i=0; i<nh; ++i ){
             Hodo2Hit *hit=hodoAna->GetHitLTOF(i);
             
             int seg=hit->SegmentId();
             event.id_ltof[i] = seg;
             event.ltdc_ltof_l[seg-1] = hit->GetTLeft();
             event.ltdc_ltof_r[seg-1] = hit->GetTRight();
             event.tot_ltof_l[seg-1] = hit->GetALeft();
             event.tot_ltof_r[seg-1] = hit->GetARight();
             event.mt_ltof[seg-1] = hit->MeanTime();
             event.de_ltof[seg-1] = hit->DeltaE();
             nh2++;
          }
          event.nhits_ltof = nh2;
       }
       //T1
       hodoAna->DecodeT1Hits( rawData );
       {
          int nh=hodoAna->GetNHitsT1();
          int nh2=0;
          for( int i=0; i<nh; ++i ){
             Hodo2Hit *hit=hodoAna->GetHitT1(i);
             
             int seg=hit->SegmentId();
             event.id_t1[i] = seg;
             event.ltdc_t1_l[seg-1] = hit->GetTLeft();
             event.ltdc_t1_r[seg-1] = hit->GetTRight();
             event.tot_t1_l[seg-1] = hit->GetALeft();
             event.tot_t1_r[seg-1] = hit->GetARight();
             event.mt_t1[seg-1] = hit->MeanTime();
            event.de_t1[seg-1] = hit->DeltaE();
            nh2++;
          }
          event.nhits_t1 = nh2;
       }
       // T1-UTOF
       {
          int nh1=hodoAna->GetNHitsUTOF();
          int nh2=hodoAna->GetNHitsT1();
          int nhits_b=0;
          for( int i1=0; i1<nh1; ++i1 ){
             Hodo2Hit *utof=hodoAna->GetHitUTOF(i1);
             double mt_utof=utof->MeanTime();
             if(!utof) continue;
             for(int i2=0; i2<nh2; ++i2 ){
                Hodo2Hit *t1=hodoAna->GetHitT1(i2);
                if(!t1) continue;
                double mt_t1=t1->MeanTime();
                
                event.beam_tof[nhits_b] = mt_t1 - mt_utof;
                nhits_b++;
             }
          }
          event.nhits_beam_tof = nhits_b;
       }
       // LTOF-UTOF
      {
         int nh1=hodoAna->GetNHitsUTOF();
         int nh2=hodoAna->GetNHitsLTOF();
         int nhits_s=0;
         int nhits_s2=0;
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
               if(f_dtof_veto){
                  event.scat_tof2[nhits_s2] = mt_ltof - mt_utof;
                  nhits_s2++;
               }
            }
         }
         event.nhits_scat_tof = nhits_s;
         event.nhits_scat_tof2 = nhits_s2;
      }
      // LTOF-T1
      {
         int nh1=hodoAna->GetNHitsT1();
         int nh2=hodoAna->GetNHitsLTOF();
         int nhits_s=0;
         for( int i1=0; i1<nh1; ++i1 ){
            Hodo2Hit *t1=hodoAna->GetHitT1(i1);
            double mt_t1=t1->MeanTime();
            if(!t1) continue;
            for(int i2=0; i2<nh2; ++i2 ){
               Hodo2Hit *ltof=hodoAna->GetHitLTOF(i2);
               if(!ltof) continue;
               double mt_ltof=ltof->MeanTime();
               
               event.scat_tof3[nhits_s] = mt_ltof - mt_t1;
               nhits_s++;
            }
         }
         event.nhits_scat_tof3 = nhits_s;
      }   
      
     //**************************************************************************
     //BDC Rawdata
#if check
       std::cout << "BDC Rawdata***********************" << std::endl;
#endif
       event.tot_bdc_l1.resize(NumOfWireBDCX);
       event.tot_bdc_l2.resize(NumOfWireBDCX);
       event.tot_bdc_l3.resize(NumOfWireBDCUV);
       event.tot_bdc_l4.resize(NumOfWireBDCUV);
       event.tot_bdc_l5.resize(NumOfWireBDCX);
       event.tot_bdc_l6.resize(NumOfWireBDCX);
       event.tot_bdc_l7.resize(NumOfWireBDCUV);
       event.tot_bdc_l8.resize(NumOfWireBDCUV);
       
       event.ltdc_bdc_l1.resize(NumOfWireBDCX);
       event.ltdc_bdc_l2.resize(NumOfWireBDCX);
       event.ltdc_bdc_l3.resize(NumOfWireBDCUV);
       event.ltdc_bdc_l4.resize(NumOfWireBDCUV);
       event.ltdc_bdc_l5.resize(NumOfWireBDCX);
       event.ltdc_bdc_l6.resize(NumOfWireBDCX);
       event.ltdc_bdc_l7.resize(NumOfWireBDCUV);
       event.ltdc_bdc_l8.resize(NumOfWireBDCUV);
       
       {
          for( int layer=1; layer<=NumOfLayersBDC; ++layer ){
             const DCRHitContainer &cont =rawData->GetBDCRHC(layer);
             int nh=cont.size();
             for( int i=0; i<nh; ++i ){
                DCRawHit *hit=cont[i];
                int layerId = hit->LayerId();
                int wireId = hit->WireId();
                
                if(layerId==1){
                   event.id_bdc_l1.push_back(wireId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it )
                      event.ltdc_bdc_l1[wireId-1].push_back(hit->GetlTdc(it));
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bdc_l1[wireId-1].push_back(hit->GetTot(it));
                }
                if(layerId==2){
                   event.id_bdc_l2.push_back(wireId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it )
                      event.ltdc_bdc_l2[wireId-1].push_back(hit->GetlTdc(it));
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bdc_l2[wireId-1].push_back(hit->GetTot(it));
                }
                if(layerId==3){
                   event.id_bdc_l3.push_back(wireId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it )
                      event.ltdc_bdc_l3[wireId-1].push_back(hit->GetlTdc(it));
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bdc_l3[wireId-1].push_back(hit->GetTot(it));
                }
                if(layerId==4){
                   event.id_bdc_l4.push_back(wireId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it )
                      event.ltdc_bdc_l4[wireId-1].push_back(hit->GetlTdc(it));
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bdc_l4[wireId-1].push_back(hit->GetTot(it));
                }
                if(layerId==5){
                   event.id_bdc_l5.push_back(wireId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it )
                      event.ltdc_bdc_l5[wireId-1].push_back(hit->GetlTdc(it));
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bdc_l5[wireId-1].push_back(hit->GetTot(it));
                }
                if(layerId==6){
                   event.id_bdc_l6.push_back(wireId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it )
                      event.ltdc_bdc_l6[wireId-1].push_back(hit->GetlTdc(it));
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bdc_l6[wireId-1].push_back(hit->GetTot(it));
                }
                if(layerId==7){
                   event.id_bdc_l7.push_back(wireId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it )
                      event.ltdc_bdc_l7[wireId-1].push_back(hit->GetlTdc(it));
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bdc_l7[wireId-1].push_back(hit->GetTot(it));
                }
                if(layerId==8){
                   event.id_bdc_l8.push_back(wireId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it )
                      event.ltdc_bdc_l8[wireId-1].push_back(hit->GetlTdc(it));
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bdc_l8[wireId-1].push_back(hit->GetTot(it));
                }
             }
          }
       }

       //**************************************************************************
       //BDC Decode
#if check
       std::cout << "BDC Decode***********************" << std::endl;
#endif

       const int dctdclow = confMan->BDCTRangeLow();
       const int dctdchigh = confMan->BDCTRangeHigh();

       DCAna->DecodeBDCRawHits( rawData );  

       int multi_DC[NumOfLayersBDC];
       for( int i=0; i<NumOfLayersBDC; i++ ) multi_DC[i]=0;
       {
          for( int layer=1; layer<=NumOfLayersBDC; ++layer ){
             const DCHitContainer &cont = DCAna->GetBDCHC(layer);
             int nh = cont.size();
             
             for( int i=0; i<nh; ++i ){
                DCHit *hit = cont[i];
                double wire = hit->GetWire();
                int nhtdc = hit->GetTdcSize();
                bool tdcflag = false;
                int tdc1st = 99999;
                
                for( int k=0; k<nhtdc; k++ ){
                   int tdc = hit->GetTdcVal(k);
                   if( dctdclow< tdc && tdc < dctdchigh ) tdcflag = true;
                   if( (dctdclow< tdc && tdc < dctdchigh) && tdc < tdc1st ) tdc1st=tdc;
                }
                HF1( 100*layer+2, tdc1st );
                HF1( 10000*layer+int(wire), tdc1st );
                
                if( tdcflag ){
                   HF1( 100*layer+1, wire-0.5 );
                   multi_DC[layer-1]++;
                }
                
                int nhdt = hit->GetDriftTimeSize();
                for( int k=0; k<nhdt; k++ ){
                   double dt = hit->GetDriftTime(k);
                   HF1( 100*layer+3, dt );
                   HF1( 10000*layer+1000+int(wire), dt );
                }
                int nhdl = hit->GetDriftTimeSize();
                for( int k=0; k<nhdl; k++ ){
                   double dl = hit->GetDriftLength(k);
                   HF1( 100*layer+4, dl );
                   HF1( 10000*layer+2000+int(wire), dl );
                }
             }
             HF1( 100*layer, multi_DC[layer-1] );
          }
       }

       //**************************************************************************
       //BDC Tracking
#if check
       std::cout << "BDC Tracking***********************" << std::endl;
#endif
       event.dt_l1.resize(NumOfWireBDCX);
       event.dt_l2.resize(NumOfWireBDCX);
       event.dt_l3.resize(NumOfWireBDCUV);
       event.dt_l4.resize(NumOfWireBDCUV);
       event.dt_l5.resize(NumOfWireBDCX);
       event.dt_l6.resize(NumOfWireBDCX);
       event.dt_l7.resize(NumOfWireBDCUV);
       event.dt_l8.resize(NumOfWireBDCUV);

       event.dl_l1.resize(NumOfWireBDCX);
       event.dl_l2.resize(NumOfWireBDCX);
       event.dl_l3.resize(NumOfWireBDCUV);
       event.dl_l4.resize(NumOfWireBDCUV);
       event.dl_l5.resize(NumOfWireBDCX);
       event.dl_l6.resize(NumOfWireBDCX);
       event.dl_l7.resize(NumOfWireBDCUV);
       event.dl_l8.resize(NumOfWireBDCUV);

       event.pos_l1.resize(NumOfWireBDCX);
       event.pos_l2.resize(NumOfWireBDCX);
       event.pos_l3.resize(NumOfWireBDCUV);
       event.pos_l4.resize(NumOfWireBDCUV);
       event.pos_l5.resize(NumOfWireBDCX);
       event.pos_l6.resize(NumOfWireBDCX);
       event.pos_l7.resize(NumOfWireBDCUV);
       event.pos_l8.resize(NumOfWireBDCUV);

       event.res_l1.resize(NumOfWireBDCX);
       event.res_l2.resize(NumOfWireBDCX);
       event.res_l3.resize(NumOfWireBDCUV);
       event.res_l4.resize(NumOfWireBDCUV);
       event.res_l5.resize(NumOfWireBDCX);
       event.res_l6.resize(NumOfWireBDCX);
       event.res_l7.resize(NumOfWireBDCUV);
       event.res_l8.resize(NumOfWireBDCUV);

       {
          DCAna->TrackSearchBDC();
          int nt=DCAna->GetNtracksBDC();
          event.nt=nt;
          
          for( int it=0; it<nt; ++it ){
             DCLocalTrack *tp=DCAna->GetTrackBDC(it);
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

             double xltof=tp->GetX( 2195.);
             double ultof=u0;
             double yltof=tp->GetY( 2195.);
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
        
             for( int ih=0; ih<nh; ++ih ){
                DCLTrackHit *hit=tp->GetHit(ih);
                int layerId=hit->GetLayer()-PlOffsBDC; 
                double wire=hit->GetWire();
                event.layer.push_back(layerId);  
                double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
                double pos=hit->GetLocalHitPos(), res=hit->GetResidual();

                double wp=hit->GetWirePosition();
                double sign=1.;
                if( pos-wp<0. ) sign=-1.;

                HF1( 100*layerId+11, wire-0.5 );
                HF1( 100*layerId+12, dt );
                HF1( 100*layerId+13, dl );
                HF1( 100*layerId+14, pos );
                HF1( 100*layerId+15, res );
                HF2( 100*layerId+16, pos, res );
                HF2( 100*layerId+17, sign*dl, res );
                
                double xlcal=hit->GetLocalCalPos();
                HF2( 100*layerId+18, xlcal-wp, dt);
                HF2( 100*layerId+19, dt, xlcal-wp);

          //       if(layerId==1){
          //          event.dt_l1[int(wire)].push_back(dt);
          //          event.dl_l1[int(wire)].push_back(sign*dl);
          //          event.pos_l1[int(wire)].push_back(pos);
          //          event.res_l1[int(wire)].push_back(res);
          //       }
          //       if(layerId==2){
          //          event.dt_l2[int(wire)].push_back(dt);
          //          event.dl_l2[int(wire)].push_back(sign*dl);
          //          event.pos_l2[int(wire)].push_back(pos);
          //          event.res_l2[int(wire)].push_back(res);
          //       }
          //       if(layerId==3){
          //          event.dt_l3[int(wire)].push_back(dt);
          //          event.dl_l3[int(wire)].push_back(sign*dl);
          //          event.pos_l3[int(wire)].push_back(pos);
          //          event.res_l3[int(wire)].push_back(res);
          //       }
          //       if(layerId==4){
          //          event.dt_l4[int(wire)].push_back(dt);
          //          event.dl_l4[int(wire)].push_back(sign*dl);
          //          event.pos_l4[int(wire)].push_back(pos);
          //          event.res_l4[int(wire)].push_back(res);
          //       }
          //       if(layerId==5){
          //          event.dt_l5[int(wire)].push_back(dt);
          //          event.dl_l5[int(wire)].push_back(sign*dl);
          //          event.pos_l5[int(wire)].push_back(pos);
          //          event.res_l5[int(wire)].push_back(res);
          //       }
          //       if(layerId==6){
          //          event.dt_l6[int(wire)].push_back(dt);
          //          event.dl_l6[int(wire)].push_back(sign*dl);
          //          event.pos_l6[int(wire)].push_back(pos);
          //          event.res_l6[int(wire)].push_back(res);
          //       }
          //       if(layerId==7){
          //          event.dt_l7[int(wire)].push_back(dt);
          //          event.dl_l7[int(wire)].push_back(sign*dl);
          //          event.pos_l7[int(wire)].push_back(pos);
          //          event.res_l7[int(wire)].push_back(res);
          //       }
          //       if(layerId==8){
          //          event.dt_l8[int(wire)].push_back(dt);
          //          event.dl_l8[int(wire)].push_back(sign*dl);
          //          event.pos_l8[int(wire)].push_back(pos);
          //          event.res_l8[int(wire)].push_back(res);
          //       }
             }
          }
       }
          
       tree->Fill();
   }
   
   return true;
}

void EventBDCTracking::InitializeEvent( void )
{
   //////Timing counter
   event.nhits_utof = -1;
   event.nhits_ltof = -1;
   event.nhits_t1   = -1;
   
   event.nhits_beam_tof = -1;
   event.nhits_scat_tof = -1;
   event.nhits_scat_tof2 = -1;
   event.nhits_scat_tof3 = -1;
   
   for( int i=0; i<MaxHits; i++){
      event.id_utof[i] = -1;
      event.id_ltof[i] = -1;
      event.id_t1[i]   = -1;
   }
   
   for( int i=0; i<NumOfSegUTOF; i++){
      event.ltdc_utof_l[i] = -999.0;
      event.ltdc_utof_r[i] = -999.0;
      event.tot_utof_l[i] = -999.0;
      event.tot_utof_r[i] = -999.0;
      event.mt_utof[i]  = -999.0;
      event.de_utof[i]  = -999.0;
   }

   for( int i=0; i<NumOfSegLTOF; i++){
      event.ltdc_ltof_l[i] = -999.0;
      event.ltdc_ltof_r[i] = -999.0;
      event.tot_ltof_l[i] = -999.0;
      event.tot_ltof_r[i] = -999.0;
      event.mt_ltof[i]  = -999.0;
      event.de_ltof[i]  = -999.0;
   }
   
   for( int i=0; i<NumOfSegT1; i++){
      event.ltdc_t1_l[i] = -999.0;
      event.ltdc_t1_r[i] = -999.0;
      event.tot_t1_l[i] = -999.0;
      event.tot_t1_r[i] = -999.0;
      event.mt_t1[i]  = -999.0;
      event.de_t1[i]  = -999.0;
   }
   
   for( int i=0; i<NumOfSegBHT; i++){
      event.ltdc_bht[i] = -999.0;
      event.tot_bht[i] = -999.0;
   }
   
   for( int i=0; i<MaxHits; i++){
      event.beam_tof[i]   = -999.0;
      event.scat_tof[i]   = -999.0;
      event.scat_tof2[i]  = -999.0;
      event.scat_tof3[i]  = -999.0;
   }

   ////BDC
   event.id_bdc_l1.clear();
   event.id_bdc_l2.clear();
   event.id_bdc_l3.clear();
   event.id_bdc_l4.clear();
   event.id_bdc_l5.clear();
   event.id_bdc_l6.clear();
   event.id_bdc_l7.clear();
   event.id_bdc_l8.clear();

   for(int i=0; i<event.tot_bdc_l1.size(); i++) event.tot_bdc_l1[i].clear();
   for(int i=0; i<event.tot_bdc_l2.size(); i++) event.tot_bdc_l2[i].clear();
   for(int i=0; i<event.tot_bdc_l3.size(); i++) event.tot_bdc_l3[i].clear();
   for(int i=0; i<event.tot_bdc_l4.size(); i++) event.tot_bdc_l4[i].clear();
   for(int i=0; i<event.tot_bdc_l5.size(); i++) event.tot_bdc_l5[i].clear();
   for(int i=0; i<event.tot_bdc_l6.size(); i++) event.tot_bdc_l6[i].clear();
   for(int i=0; i<event.tot_bdc_l7.size(); i++) event.tot_bdc_l7[i].clear();
   for(int i=0; i<event.tot_bdc_l8.size(); i++) event.tot_bdc_l8[i].clear();
   for(int i=0; i<event.ltdc_bdc_l1.size(); i++) event.ltdc_bdc_l1[i].clear();
   for(int i=0; i<event.ltdc_bdc_l2.size(); i++) event.ltdc_bdc_l2[i].clear();
   for(int i=0; i<event.ltdc_bdc_l3.size(); i++) event.ltdc_bdc_l3[i].clear();
   for(int i=0; i<event.ltdc_bdc_l4.size(); i++) event.ltdc_bdc_l4[i].clear();
   for(int i=0; i<event.ltdc_bdc_l5.size(); i++) event.ltdc_bdc_l5[i].clear();
   for(int i=0; i<event.ltdc_bdc_l6.size(); i++) event.ltdc_bdc_l6[i].clear();
   for(int i=0; i<event.ltdc_bdc_l7.size(); i++) event.ltdc_bdc_l7[i].clear();
   for(int i=0; i<event.ltdc_bdc_l8.size(); i++) event.ltdc_bdc_l8[i].clear();

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

   for(int i=0; i<event.dt_l1.size(); i++) event.dt_l1[i].clear();
   for(int i=0; i<event.dt_l2.size(); i++) event.dt_l2[i].clear();
   for(int i=0; i<event.dt_l3.size(); i++) event.dt_l3[i].clear();
   for(int i=0; i<event.dt_l4.size(); i++) event.dt_l4[i].clear();
   for(int i=0; i<event.dt_l5.size(); i++) event.dt_l5[i].clear();
   for(int i=0; i<event.dt_l6.size(); i++) event.dt_l6[i].clear();
   for(int i=0; i<event.dt_l7.size(); i++) event.dt_l7[i].clear();
   for(int i=0; i<event.dt_l8.size(); i++) event.dt_l8[i].clear();

   for(int i=0; i<event.dl_l1.size(); i++) event.dl_l1[i].clear();
   for(int i=0; i<event.dl_l2.size(); i++) event.dl_l2[i].clear();
   for(int i=0; i<event.dl_l3.size(); i++) event.dl_l3[i].clear();
   for(int i=0; i<event.dl_l4.size(); i++) event.dl_l4[i].clear();
   for(int i=0; i<event.dl_l5.size(); i++) event.dl_l5[i].clear();
   for(int i=0; i<event.dl_l6.size(); i++) event.dl_l6[i].clear();
   for(int i=0; i<event.dl_l7.size(); i++) event.dl_l7[i].clear();
   for(int i=0; i<event.dl_l8.size(); i++) event.dl_l8[i].clear();

   for(int i=0; i<event.pos_l1.size(); i++) event.pos_l1[i].clear();
   for(int i=0; i<event.pos_l2.size(); i++) event.pos_l2[i].clear();
   for(int i=0; i<event.pos_l3.size(); i++) event.pos_l3[i].clear();
   for(int i=0; i<event.pos_l4.size(); i++) event.pos_l4[i].clear();
   for(int i=0; i<event.pos_l5.size(); i++) event.pos_l5[i].clear();
   for(int i=0; i<event.pos_l6.size(); i++) event.pos_l6[i].clear();
   for(int i=0; i<event.pos_l7.size(); i++) event.pos_l7[i].clear();
   for(int i=0; i<event.pos_l8.size(); i++) event.pos_l8[i].clear();

   for(int i=0; i<event.res_l1.size(); i++) event.res_l1[i].clear();
   for(int i=0; i<event.res_l2.size(); i++) event.res_l2[i].clear();
   for(int i=0; i<event.res_l3.size(); i++) event.res_l3[i].clear();
   for(int i=0; i<event.res_l4.size(); i++) event.res_l4[i].clear();
   for(int i=0; i<event.res_l5.size(); i++) event.res_l5[i].clear();
   for(int i=0; i<event.res_l6.size(); i++) event.res_l6[i].clear();
   for(int i=0; i<event.res_l7.size(); i++) event.res_l7[i].clear();
   for(int i=0; i<event.res_l8.size(); i++) event.res_l8[i].clear();
}

bool EventBDCTracking::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();

  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventBDCTracking;
}

bool ConfMan:: InitializeHistograms()
{  
   HBTree("tree","tree");
   TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
   
   //////Timing counter
   tree->Branch("nhits_utof", &event.nhits_utof,  "nhits_utof/I");
   tree->Branch("id_utof",     event.id_utof,     "id_utof[32]/I");
   tree->Branch("ltdc_utof_l", event.ltdc_utof_l, "ltdc_utof_l[1]/D");
   tree->Branch("ltdc_utof_r", event.ltdc_utof_r, "ltdc_utof_r[1]/D");
   tree->Branch("tot_utof_l",  event.tot_utof_l,  "tot_utof_l[1]/D");
   tree->Branch("tot_utof_r",  event.tot_utof_r,  "tot_utof_l[1]/D");
   tree->Branch("mt_utof",     event.mt_utof,     "mt_utof_l[1]/D");
   tree->Branch("de_utof",     event.de_utof,     "de_utof_l[1]/D");

   tree->Branch("nhits_ltof", &event.nhits_ltof,  "nhits_ltof/I");
   tree->Branch("id_ltof",     event.id_ltof,     "id_ltof[32]/I");
   tree->Branch("ltdc_ltof_l", event.ltdc_ltof_l, "ltdc_ltof_l[6]/D");
   tree->Branch("ltdc_ltof_r", event.ltdc_ltof_r, "ltdc_ltof_r[6]/D");
   tree->Branch("tot_ltof_l",  event.tot_ltof_l,  "tot_ltof_l[6]/D");
   tree->Branch("tot_ltof_r",  event.tot_ltof_r,  "tot_ltof_l[6]/D");
   tree->Branch("mt_ltof",     event.mt_ltof,     "mt_ltof_l[6]/D");
   tree->Branch("de_ltof",     event.de_ltof,     "de_ltof_l[6]/D");

   tree->Branch("nhits_t1", &event.nhits_t1,  "nhits_t1/I");
   tree->Branch("id_t1",     event.id_t1,     "id_t1[32]/I");
   tree->Branch("ltdc_t1_l", event.ltdc_t1_l, "ltdc_t1_l[1]/D");
   tree->Branch("ltdc_t1_r", event.ltdc_t1_r, "ltdc_t1_r[1]/D");
   tree->Branch("tot_t1_l",  event.tot_t1_l,  "tot_t1_l[1]/D");
   tree->Branch("tot_t1_r",  event.tot_t1_r,  "tot_t1_l[1]/D");
   tree->Branch("mt_t1",     event.mt_t1,     "mt_t1_l[1]/D");
   tree->Branch("de_t1",     event.de_t1,     "de_t1_l[1]/D");
   
   tree->Branch("nhits_beam_tof", &event.nhits_beam_tof, "nhits_beam_tof/I");
   tree->Branch("beam_tof",        event.beam_tof,       "beam_tof[32]/D");  
   
   tree->Branch("nhits_scat_tof", &event.nhits_scat_tof, "nhits_scat_tof/I");
   tree->Branch("scat_tof",        event.scat_tof,       "scat_tof[32]/D");  

   tree->Branch("nhits_scat_tof2", &event.nhits_scat_tof2, "nhits_scat_tof2/I");
   tree->Branch("scat_tof2",        event.scat_tof2,       "scat_tof2[32]/D");  

   tree->Branch("nhits_scat_tof3", &event.nhits_scat_tof3, "nhits_scat_tof3/I");
   tree->Branch("scat_tof3",        event.scat_tof3,       "scat_tof3[32]/D");  
   
   //BDC
   tree->Branch("id_bdc_l1", &event.id_bdc_l1);
   tree->Branch("id_bdc_l2", &event.id_bdc_l2);
   tree->Branch("id_bdc_l3", &event.id_bdc_l3);
   tree->Branch("id_bdc_l4", &event.id_bdc_l4);
   tree->Branch("id_bdc_l5", &event.id_bdc_l5);
   tree->Branch("id_bdc_l6", &event.id_bdc_l6);
   tree->Branch("id_bdc_l7", &event.id_bdc_l7);
   tree->Branch("id_bdc_l8", &event.id_bdc_l8);
   
   tree->Branch("tot_bdc_l1", &event.tot_bdc_l1);
   tree->Branch("tot_bdc_l2", &event.tot_bdc_l2);
   tree->Branch("tot_bdc_l3", &event.tot_bdc_l3);
   tree->Branch("tot_bdc_l4", &event.tot_bdc_l4);
   tree->Branch("tot_bdc_l5", &event.tot_bdc_l5);
   tree->Branch("tot_bdc_l6", &event.tot_bdc_l6);
   tree->Branch("tot_bdc_l7", &event.tot_bdc_l7);
   tree->Branch("tot_bdc_l8", &event.tot_bdc_l8);
   
   tree->Branch("ltdc_bdc_l1", &event.ltdc_bdc_l1);
   tree->Branch("ltdc_bdc_l2", &event.ltdc_bdc_l2);
   tree->Branch("ltdc_bdc_l3", &event.ltdc_bdc_l3);
   tree->Branch("ltdc_bdc_l4", &event.ltdc_bdc_l4);
   tree->Branch("ltdc_bdc_l5", &event.ltdc_bdc_l5);
   tree->Branch("ltdc_bdc_l6", &event.ltdc_bdc_l6);
   tree->Branch("ltdc_bdc_l7", &event.ltdc_bdc_l7);
   tree->Branch("ltdc_bdc_l8", &event.ltdc_bdc_l8);

   //Tracking
   tree->Branch("nt", &event.nt);
   tree->Branch("layer", &event.layer);
   tree->Branch("chisqr", &event.chisqr);
   tree->Branch("x0", &event.x0);
   tree->Branch("u0", &event.u0);
   tree->Branch("y0", &event.y0);
   tree->Branch("v0", &event.v0);
   tree->Branch("x1", &event.x1);
   tree->Branch("u1", &event.u1);
   tree->Branch("y1", &event.y1);
   tree->Branch("v1", &event.v1);
   
   tree->Branch("dt_l1", &event.dt_l1);
   tree->Branch("dt_l2", &event.dt_l2);
   tree->Branch("dt_l3", &event.dt_l3);
   tree->Branch("dt_l4", &event.dt_l4);
   tree->Branch("dt_l5", &event.dt_l5);
   tree->Branch("dt_l6", &event.dt_l6);
   tree->Branch("dt_l7", &event.dt_l7);
   tree->Branch("dt_l8", &event.dt_l8);

   tree->Branch("dl_l1", &event.dl_l1);
   tree->Branch("dl_l2", &event.dl_l2);
   tree->Branch("dl_l3", &event.dl_l3);
   tree->Branch("dl_l4", &event.dl_l4);
   tree->Branch("dl_l5", &event.dl_l5);
   tree->Branch("dl_l6", &event.dl_l6);
   tree->Branch("dl_l7", &event.dl_l7);
   tree->Branch("dl_l8", &event.dl_l8);

   tree->Branch("pos_l1", &event.pos_l1);
   tree->Branch("pos_l2", &event.pos_l2);
   tree->Branch("pos_l3", &event.pos_l3);
   tree->Branch("pos_l4", &event.pos_l4);
   tree->Branch("pos_l5", &event.pos_l5);
   tree->Branch("pos_l6", &event.pos_l6);
   tree->Branch("pos_l7", &event.pos_l7);
   tree->Branch("pos_l8", &event.pos_l8);

   tree->Branch("res_l1", &event.res_l1);
   tree->Branch("res_l2", &event.res_l2);
   tree->Branch("res_l3", &event.res_l3);
   tree->Branch("res_l4", &event.res_l4);
   tree->Branch("res_l5", &event.res_l5);
   tree->Branch("res_l6", &event.res_l6);
   tree->Branch("res_l7", &event.res_l7);
   tree->Branch("res_l8", &event.res_l8);

   //BDC Histograms
   for( int i=1; i<=NumOfLayersBDC; ++i ){
      std::ostringstream title1, title2, title3, title4, title5;
      title1 << "#Hits BDC#" << std::setw(2) << i;
      title2 << "Hitpat BDC#" << std::setw(2) << i;
      title3 << "Tdc BDC#" << std::setw(2) << i;
      title4 << "Drift Time BDC#" << std::setw(2) << i;
      title5 << "Drift Length BDC#" << std::setw(2) << i;

      HB1( 100*i+0, title1.str().c_str(), NumOfWireBDCX+1, 0., double(NumOfWireBDCX+1) );
      HB1( 100*i+1, title2.str().c_str(), NumOfWireBDCX+1, 0., double(NumOfWireBDCX+1) );
      HB1( 100*i+2, title3.str().c_str(), 800, 800, 1600 );
      HB1( 100*i+3, title4.str().c_str(), 750, -50., 700. );
      HB1( 100*i+4, title5.str().c_str(), 125, -0.5, 12.0 );
      
      std::ostringstream title10, title11, title12, title13, title14, title15, title16, title17;
      title10 << "wire for LayerId = " << std::setw(2) << i<< " [Track]";
      title11 << "drift time for LayerId = " << std::setw(2) << i<< " [Track]";
      title12 << "drift length for LayerId = " << std::setw(2) << i<< " [Track]";
      title13 << "Position " << std::setw(2) << i;
      title14 << "Residual " << std::setw(2) << i;
      title15 << "Resid%Pos " << std::setw(2) << i;
      title16 << "Resid%dl " << std::setw(2) << i;
      title17 << "Drift Length%Drift Time " << std::setw(2) << i;
      
      HB1( 100*i+11, title10.str().c_str(), NumOfWireBDCX+1, 0., double(NumOfWireBDCX+1) );
      HB1( 100*i+12, title11.str().c_str(), 310, -10, 300 );
      HB1( 100*i+13, title12.str().c_str(), 125, -0.5, 12.0);
      HB1( 100*i+14, title13.str().c_str(), 300, -1050., 1050. ); 
      HB1( 100*i+15, title14.str().c_str(), 200, -2.0, 2.0 );
      HB2( 100*i+16, title15.str().c_str(), 100, -1500., 1050., 100, -2.0, 2.0 );
      HB2( 100*i+17, title16.str().c_str(), 100, -15.0, 15.0, 100, -2.0, 2.0 );
      HB2( 100*i+18, title17.str().c_str(), 100, -15.0, 15.0, 100, -10, 350 );
      HB2( 100*i+19, title17.str().c_str(), 100, -10, 150., 100, -15.0, 15.0 );
      
      for (int wire=1; wire<=NumOfWireBDCX; wire++) {
         std::ostringstream title11, title12, title13;
         title11 << "Tdc DC#" << std::setw(2) << i << " Wire#" << wire;
         HB1( 10000*i+wire, title11.str().c_str(), 1000, 0, 1000 );
         title12 << "Drift Time #" << std::setw(2) << i << " Wire#" << wire;
         HB1( 10000*i+1000+wire, title12.str().c_str(), 500, -100., 500. );
         title13 << "Drift Length #" << std::setw(2) << i << " Wire#" << wire;
         HB1( 10000*i+2000+wire, title13.str().c_str(), 100, -5., 15. );
      }
   }
   
  return true;
}
