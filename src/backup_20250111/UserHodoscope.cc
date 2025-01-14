/*
  UserHodoscope.cc

  2024/02  K.Shirotori 
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

#include "TFile.h"
#include "TTree.h"

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventHodoscope : public VEvent
{

public:
  EventHodoscope();
  ~EventHodoscope();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( TFile*, int evnum );
  void InitializeEvent();

private:
  RawData *rawData;
  HodoAnalyzer *hodoAna;
};

EventHodoscope::EventHodoscope()
  : VEvent(),
    rawData(0),
    hodoAna(new HodoAnalyzer)
{
}

EventHodoscope::~EventHodoscope()
{
  delete hodoAna;
  if (rawData) delete rawData;
}

#ifndef MaxHits 
#define MaxHits 32
#endif

struct Event{
   int nhits_utof;
   int id_utof[MaxHits];
   double ltdc_utof_l[NumOfSegUTOF];
   double ltdc_utof_r[NumOfSegUTOF];
   double tot_utof_l[NumOfSegUTOF];
   double tot_utof_r[NumOfSegUTOF];
   double mt_utof[NumOfSegUTOF];
   double de_utof[NumOfSegUTOF];

   int nhits_dtof;
   int id_dtof[MaxHits];
   double ltdc_dtof_l[NumOfSegDTOF];
   double ltdc_dtof_r[NumOfSegDTOF];
   double tot_dtof_l[NumOfSegDTOF];
   double tot_dtof_r[NumOfSegDTOF];
   double mt_dtof[NumOfSegDTOF];
   double de_dtof[NumOfSegDTOF];

   int nhits_ltof;
   int id_ltof[MaxHits];
   double ltdc_ltof_l[NumOfSegLTOF];
   double ltdc_ltof_r[NumOfSegLTOF];
   double tot_ltof_l[NumOfSegLTOF];
   double tot_ltof_r[NumOfSegLTOF];
   double mt_ltof[NumOfSegLTOF];
   double de_ltof[NumOfSegLTOF];

   int nhits_t0;
   int id_t0[MaxHits];
   double ltdc_t0_l[NumOfSegT0];
   double ltdc_t0_r[NumOfSegT0];
   double tot_t0_l[NumOfSegT0];
   double tot_t0_r[NumOfSegT0];
   double mt_t0[NumOfSegT0];
   double de_t0[NumOfSegT0];

   int nhits_t0r;
   int id_t0r[MaxHits];
   double ltdc_t0r_l[NumOfSegT0r];
   double ltdc_t0r_r[NumOfSegT0r];
   double tot_t0r_l[NumOfSegT0r];
   double tot_t0r_r[NumOfSegT0r];
   double mt_t0r[NumOfSegT0r];
   double de_t0r[NumOfSegT0r];

   int nhits_bref;
   int id_bref[MaxHits];
   double ltdc_bref[NumOfSegBref];
   double tot_bref[NumOfSegBref];

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
   double de_tof3[MaxHits];
};
static Event event;

bool EventHodoscope::ProcessingBegin()
{
 return true;
}

bool EventHodoscope::ProcessingNormal( TFile *iFile, int evnum )
{
   const std::string funcname = "ProcessingNormal";
   
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
      
      TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
      
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
      //DTOF
      hodoAna->DecodeDTOFHits( rawData );
      {
         int nh=hodoAna->GetNHitsDTOF();
         int nh2=0;
         for( int i=0; i<nh; ++i ){
            Hodo2Hit *hit=hodoAna->GetHitDTOF(i);
            
            int seg=hit->SegmentId();
            event.id_dtof[i] = seg;
            event.ltdc_dtof_l[seg-1] = hit->GetTLeft();
            event.ltdc_dtof_r[seg-1] = hit->GetTRight();
            event.tot_dtof_l[seg-1] = hit->GetALeft();
            event.tot_dtof_r[seg-1] = hit->GetARight();
            event.mt_dtof[seg-1] = hit->MeanTime();
            event.de_dtof[seg-1] = hit->DeltaE();
            nh2++;
         }
         event.nhits_dtof = nh2;
         if(nh2>0) f_dtof_veto = true;
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
      //T0
      hodoAna->DecodeT0Hits( rawData );
      {
         int nh=hodoAna->GetNHitsT0();
         int nh2=0;
         for( int i=0; i<nh; ++i ){
            Hodo2Hit *hit=hodoAna->GetHitT0(i);
            
            int seg=hit->SegmentId();
            event.id_t0[i] = seg;
            event.ltdc_t0_l[seg-1] = hit->GetTLeft();
            event.ltdc_t0_r[seg-1] = hit->GetTRight();
            event.tot_t0_l[seg-1] = hit->GetALeft();
            event.tot_t0_r[seg-1] = hit->GetARight();
            event.mt_t0[seg-1] = hit->MeanTime();
            event.de_t0[seg-1] = hit->DeltaE();
            nh2++;
         }
         event.nhits_t0 = nh2;
      }
      //T0R
      hodoAna->DecodeT0rHits( rawData );
      {
         int nh=hodoAna->GetNHitsT0r();
         int nh2=0;
         for( int i=0; i<nh; ++i ){
            Hodo2Hit *hit=hodoAna->GetHitT0r(i);
            
            int seg=hit->SegmentId();
            event.id_t0r[i] = seg;
            event.ltdc_t0r_l[seg-1] = hit->GetTLeft();
            event.ltdc_t0r_r[seg-1] = hit->GetTRight();
            event.tot_t0r_l[seg-1] = hit->GetALeft();
            event.tot_t0r_r[seg-1] = hit->GetARight();
            event.mt_t0r[seg-1] = hit->MeanTime();
            event.de_t0r[seg-1] = hit->DeltaE();
            nh2++;
         }
         event.nhits_t0r = nh2;
      }
      //Bref
      hodoAna->DecodeBrefHits( rawData );
      {
         int nh=hodoAna->GetNHitsBref();
         int nh2=0;
         for( int i=0; i<nh; ++i ){
            Hodo1Hit *hit=hodoAna->GetHitBref(i);
            
            int seg=hit->SegmentId();
            event.id_bref[i] = seg;
            event.ltdc_bref[seg-1] = hit->GetT();
            event.tot_bref[seg-1] = hit->GetA();
            nh2++;
         }
         event.nhits_bref = nh2;
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
      //BHT
      hodoAna->DecodeBHTHits( rawData );
      {
         int nh=hodoAna->GetNHitsBHT();
         int nh2=0;
         for( int i=0; i<nh; ++i ){
            Hodo1Hit *hit=hodoAna->GetHitBHT(i);
            
            int seg=hit->SegmentId();
            event.id_bref[i] = seg;
            event.ltdc_bref[seg-1] = hit->GetT();
            event.tot_bref[seg-1] = hit->GetA();
            nh2++;
         }
         event.nhits_bht = nh2;
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

               event.de_tof3[nhits_s] = ltof->DeltaE();
               
               nhits_s++;
            }
         }
         event.nhits_scat_tof3 = nhits_s;
      }   
      
      tree->Fill();
   }

   return true;
}


void EventHodoscope::InitializeEvent( void )
{
   event.nhits_utof = -1;
   event.nhits_dtof = -1;
   event.nhits_ltof = -1;
   event.nhits_t0   = -1;
   event.nhits_t0r  = -1;
   event.nhits_bref = -1;
   event.nhits_t1   = -1;
   event.nhits_bht  = -1;

   event.nhits_beam_tof  = -1;
   event.nhits_scat_tof  = -1;
   event.nhits_scat_tof2 = -1;
   event.nhits_scat_tof3 = -1;
   
  for( int i=0; i<MaxHits; i++){
     event.id_utof[i] = -1;
     event.id_dtof[i] = -1;
     event.id_ltof[i] = -1;
     event.id_t0[i]   = -1;
     event.id_t0r[i]  = -1;
     event.id_bref[i] = -1;
     event.id_t1[i]   = -1;
     event.id_bht[i]  = -1;
  }
  
  for( int i=0; i<NumOfSegUTOF; i++){
     event.ltdc_utof_l[i] = -999.0;
     event.ltdc_utof_r[i] = -999.0;
     event.tot_utof_l[i] = -999.0;
     event.tot_utof_r[i] = -999.0;
     event.mt_utof[i]  = -999.0;
     event.de_utof[i]  = -999.0;
  }

  for( int i=0; i<NumOfSegDTOF; i++){
     event.ltdc_dtof_l[i] = -999.0;
     event.ltdc_dtof_r[i] = -999.0;
     event.tot_dtof_l[i] = -999.0;
     event.tot_dtof_r[i] = -999.0;
     event.mt_dtof[i]  = -999.0;
     event.de_dtof[i]  = -999.0;
  }

  for( int i=0; i<NumOfSegLTOF; i++){
     event.ltdc_ltof_l[i] = -999.0;
     event.ltdc_ltof_r[i] = -999.0;
     event.tot_ltof_l[i] = -999.0;
     event.tot_ltof_r[i] = -999.0;
     event.mt_ltof[i]  = -999.0;
     event.de_ltof[i]  = -999.0;
  }

  for( int i=0; i<NumOfSegT0; i++){
     event.ltdc_t0_l[i] = -999.0;
     event.ltdc_t0_r[i] = -999.0;
     event.tot_t0_l[i] = -999.0;
     event.tot_t0_r[i] = -999.0;
     event.mt_t0[i]  = -999.0;
     event.de_t0[i]  = -999.0;
  }

  for( int i=0; i<NumOfSegT0r; i++){
     event.ltdc_t0r_l[i] = -999.0;
     event.ltdc_t0r_r[i] = -999.0;
     event.tot_t0r_l[i] = -999.0;
     event.tot_t0r_r[i] = -999.0;
     event.mt_t0r[i]  = -999.0;
     event.de_t0r[i]  = -999.0;
  }

  for( int i=0; i<NumOfSegBref; i++){
     event.ltdc_bref[i] = -999.0;
     event.tot_bref[i] = -999.0;
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
     event.de_tof3[i]  = -999.0;
  }
}

bool EventHodoscope::ProcessingEnd()
{
   // gFile->Write();
   // gFile->Close();
   
   return true;
}

VEvent *ConfMan::EventAllocator()
{
   return new EventHodoscope;
}

bool ConfMan:: InitializeHistograms()
{
   HBTree("tree","tree");
   TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
   
   tree->Branch("nhits_utof", &event.nhits_utof,  "nhits_utof/I");
   tree->Branch("id_utof",     event.id_utof,     "id_utof[32]/I");
   tree->Branch("ltdc_utof_l", event.ltdc_utof_l, "ltdc_utof_l[1]/D");
   tree->Branch("ltdc_utof_r", event.ltdc_utof_r, "ltdc_utof_r[1]/D");
   tree->Branch("tot_utof_l",  event.tot_utof_l,  "tot_utof_l[1]/D");
   tree->Branch("tot_utof_r",  event.tot_utof_r,  "tot_utof_l[1]/D");
   tree->Branch("mt_utof",     event.mt_utof,     "mt_utof_l[1]/D");
   tree->Branch("de_utof",     event.de_utof,     "de_utof_l[1]/D");

   tree->Branch("nhits_dtof", &event.nhits_dtof,  "nhits_dtof/I");
   tree->Branch("id_dtof",     event.id_dtof,     "id_dtof[32]/I");
   tree->Branch("ltdc_dtof_l", event.ltdc_dtof_l, "ltdc_dtof_l[3]/D");
   tree->Branch("ltdc_dtof_r", event.ltdc_dtof_r, "ltdc_dtof_r[3]/D");
   tree->Branch("tot_dtof_l",  event.tot_dtof_l,  "tot_dtof_l[3]/D");
   tree->Branch("tot_dtof_r",  event.tot_dtof_r,  "tot_dtof_l[3]/D");
   tree->Branch("mt_dtof",     event.mt_dtof,     "mt_dtof_l[3]/D");
   tree->Branch("de_dtof",     event.de_dtof,     "de_dtof_l[3]/D");

   tree->Branch("nhits_ltof", &event.nhits_ltof,  "nhits_ltof/I");
   tree->Branch("id_ltof",     event.id_ltof,     "id_ltof[32]/I");
   tree->Branch("ltdc_ltof_l", event.ltdc_ltof_l, "ltdc_ltof_l[6]/D");
   tree->Branch("ltdc_ltof_r", event.ltdc_ltof_r, "ltdc_ltof_r[6]/D");
   tree->Branch("tot_ltof_l",  event.tot_ltof_l,  "tot_ltof_l[6]/D");
   tree->Branch("tot_ltof_r",  event.tot_ltof_r,  "tot_ltof_l[6]/D");
   tree->Branch("mt_ltof",     event.mt_ltof,     "mt_ltof_l[6]/D");
   tree->Branch("de_ltof",     event.de_ltof,     "de_ltof_l[6]/D");

   tree->Branch("nhits_t0", &event.nhits_t0,  "nhits_t0/I");
   tree->Branch("id_t0",     event.id_t0,     "id_t0[32]/I");
   tree->Branch("ltdc_t0_l", event.ltdc_t0_l, "ltdc_t0_l[8]/D");
   tree->Branch("ltdc_t0_r", event.ltdc_t0_r, "ltdc_t0_r[8]/D");
   tree->Branch("tot_t0_l",  event.tot_t0_l,  "tot_t0_l[8]/D");
   tree->Branch("tot_t0_r",  event.tot_t0_r,  "tot_t0_l[8]/D");
   tree->Branch("mt_t0",     event.mt_t0,     "mt_t0_l[8]/D");
   tree->Branch("de_t0",     event.de_t0,     "de_t0_l[8]/D");

   tree->Branch("nhits_t0r", &event.nhits_t0r,  "nhits_t0r/I");
   tree->Branch("id_t0r",     event.id_t0r,     "id_t0r[32]/I");
   tree->Branch("ltdc_t0r_l", event.ltdc_t0r_l, "ltdc_t0r_l[1]/D");
   tree->Branch("ltdc_t0r_r", event.ltdc_t0r_r, "ltdc_t0r_r[1]/D");
   tree->Branch("tot_t0r_l",  event.tot_t0r_l,  "tot_t0r_l[1]/D");
   tree->Branch("tot_t0r_r",  event.tot_t0r_r,  "tot_t0r_l[1]/D");
   tree->Branch("mt_t0r",     event.mt_t0r,     "mt_t0r_l[1]/D");
   tree->Branch("de_t0r",     event.de_t0r,     "de_t0r_l[1]/D");

   tree->Branch("nhits_bref", &event.nhits_bref, "nhits_bref/I");
   tree->Branch("id_bref",     event.id_bref,    "id_bref[32]/I");
   tree->Branch("ltdc_bref",   event.ltdc_bref,  "ltdc_bref[2]/D");
   tree->Branch("tot_bref",    event.tot_bref,   "tot_bref[2]/D");

   tree->Branch("nhits_t1", &event.nhits_t1,  "nhits_t1/I");
   tree->Branch("id_t1",     event.id_t1,     "id_t1[32]/I");
   tree->Branch("ltdc_t1_l", event.ltdc_t1_l, "ltdc_t1_l[1]/D");
   tree->Branch("ltdc_t1_r", event.ltdc_t1_r, "ltdc_t1_r[1]/D");
   tree->Branch("tot_t1_l",  event.tot_t1_l,  "tot_t1_l[1]/D");
   tree->Branch("tot_t1_r",  event.tot_t1_r,  "tot_t1_l[1]/D");
   tree->Branch("mt_t1",     event.mt_t1,     "mt_t1_l[1]/D");
   tree->Branch("de_t1",     event.de_t1,     "de_t1_l[1]/D");

   tree->Branch("nhits_bht", &event.nhits_bht, "nhits_bht/I");
   tree->Branch("id_bht",     event.id_bht,    "id_bht[32]/I");
   tree->Branch("ltdc_bht",   event.ltdc_bht,  "ltdc_bht[1]/D");
   tree->Branch("tot_bht",    event.tot_bht,   "tot_bht[1]/D");

   tree->Branch("nhits_beam_tof", &event.nhits_beam_tof, "nhits_beam_tof/I");
   tree->Branch("beam_tof",        event.beam_tof,       "beam_tof[32]/D");  
   
   tree->Branch("nhits_scat_tof", &event.nhits_scat_tof, "nhits_scat_tof/I");
   tree->Branch("scat_tof",        event.scat_tof,       "scat_tof[32]/D");  
   
   tree->Branch("nhits_scat_tof2", &event.nhits_scat_tof2, "nhits_scat_tof2/I");
   tree->Branch("scat_tof2",        event.scat_tof2,       "scat_tof2[32]/D");  

   tree->Branch("nhits_scat_tof3", &event.nhits_scat_tof3, "nhits_scat_tof3/I");
   tree->Branch("scat_tof3",        event.scat_tof3,       "scat_tof3[32]/D");  
   tree->Branch("de_tof3",          event.de_tof3,         "de_tof3[32]/D");  
   
  return true;
}
