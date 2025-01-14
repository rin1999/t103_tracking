//Brief Event Builder: 2024/04 K.Shirotori

#include "TTimingChargeData.h"
#include "TFile.h"
#include "TTree.h"

#include "RootHelper.hh"
#include "DetectorInfo.hh"

////////////////////////////////////////////////////////////////////////////////////
struct Event{
   //////Timing counter
   //UTOF
   std::vector<int> id_utof_l, id_utof_r;  
   std::vector<std::vector<double>> tot_utof_l, tot_utof_r;  
   std::vector<std::vector<double>> ltdc_utof_l, ltdc_utof_r;  

   //DTOF
   std::vector<int> id_dtof_l, id_dtof_r;  
   std::vector<std::vector<double>> tot_dtof_l, tot_dtof_r;  
   std::vector<std::vector<double>> ltdc_dtof_l, ltdc_dtof_r;  

   //LTOF
   std::vector<int> id_ltof_l, id_ltof_r;  
   std::vector<std::vector<double>> tot_ltof_l, tot_ltof_r;  
   std::vector<std::vector<double>> ltdc_ltof_l, ltdc_ltof_r;  

   //T0
   std::vector<int> id_t0_l, id_t0_r;  
   std::vector<std::vector<double>> tot_t0_l, tot_t0_r;  
   std::vector<std::vector<double>> ltdc_t0_l, ltdc_t0_r;  
   
   //T0r
   std::vector<int> id_t0r_l, id_t0r_r;  
   std::vector<std::vector<double>> tot_t0r_l, tot_t0r_r;  
   std::vector<std::vector<double>> ltdc_t0r_l, ltdc_t0r_r;  

   //Bref
   std::vector<int> id_bref_l, id_bref_r;  
   std::vector<std::vector<double>> tot_bref_l, tot_bref_r;  
   std::vector<std::vector<double>> ltdc_bref_l, ltdc_bref_r;  

   //T1 from K1.8BR   
   std::vector<int> id_t1_l, id_t1_r;  
   std::vector<std::vector<double>> tot_t1_l, tot_t1_r;  
   std::vector<std::vector<double>> ltdc_t1_l, ltdc_t1_r;  

   //BHT from K1.8BR
   std::vector<int> id_bht_l, id_bht_r;  
   std::vector<std::vector<double>> tot_bht_l, tot_bht_r;  
   std::vector<std::vector<double>> ltdc_bht_l, ltdc_bht_r;  

   //////Drift chamber
   //BDC
   std::vector<int> id_bdc_l1, id_bdc_l2, id_bdc_l3, id_bdc_l4;
   std::vector<int> id_bdc_l5, id_bdc_l6, id_bdc_l7, id_bdc_l8;
   std::vector<std::vector<double>> tot_bdc_l1, tot_bdc_l2, tot_bdc_l3, tot_bdc_l4;  
   std::vector<std::vector<double>> tot_bdc_l5, tot_bdc_l6, tot_bdc_l7, tot_bdc_l8;  
   std::vector<std::vector<double>> ltdc_bdc_l1, ltdc_bdc_l2, ltdc_bdc_l3, ltdc_bdc_l4;  
   std::vector<std::vector<double>> ltdc_bdc_l5, ltdc_bdc_l6, ltdc_bdc_l7, ltdc_bdc_l8;  

   //KLDC
   std::vector<int> id_kldc_l1, id_kldc_l2, id_kldc_l3, id_kldc_l4;
   std::vector<int> id_kldc_l5, id_kldc_l6, id_kldc_l7, id_kldc_l8;
   std::vector<std::vector<double>> tot_kldc_l1, tot_kldc_l2, tot_kldc_l3, tot_kldc_l4;  
   std::vector<std::vector<double>> tot_kldc_l5, tot_kldc_l6, tot_kldc_l7, tot_kldc_l8;  
   std::vector<std::vector<double>> ltdc_kldc_l1, ltdc_kldc_l2, ltdc_kldc_l3, ltdc_kldc_l4;  
   std::vector<std::vector<double>> ltdc_kldc_l5, ltdc_kldc_l6, ltdc_kldc_l7, ltdc_kldc_l8;  

   //////Fiber tracker
   //BFT
   std::vector<int> id_bft_l1, id_bft_l2, id_bft_l3, id_bft_l4, id_bft_l5, id_bft_l6;
   std::vector<std::vector<double>> tot_bft_l1, tot_bft_l2, tot_bft_l3, tot_bft_l4, tot_bft_l5, tot_bft_l6;
   std::vector<std::vector<double>> ltdc_bft_l1, ltdc_bft_l2, ltdc_bft_l3, ltdc_bft_l4, ltdc_bft_l5, ltdc_bft_l6;

   //SFT
   std::vector<int> id_sft_l1, id_sft_l2, id_sft_l3, id_sft_l4, id_sft_l5, id_sft_l6;
   std::vector<std::vector<double>> tot_sft_l1, tot_sft_l2, tot_sft_l3, tot_sft_l4, tot_sft_l5, tot_sft_l6;
   std::vector<std::vector<double>> ltdc_sft_l1, ltdc_sft_l2, ltdc_sft_l3, ltdc_sft_l4, ltdc_sft_l5, ltdc_sft_l6;
};
static Event event;

////////////////////////////////////////////////////////////////////////////////////
void InitializeEvent();

void DefineHistograms( const char *filename );

////////////////////////////////////////////////////////////////////////////////////
////Main
int EventBuilder00( int run, const char* filename, int nevent )
{
   char rootfile[100];
   sprintf(rootfile, "after_artemis_tfb00/run%06d.root",run);
   cout << rootfile << endl;

   TFile *fin1 = TFile::Open(rootfile);
   TTree *tree1 = (TTree *)gDirectory->Get("tree");
   
   DefineHistograms( filename );
   TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

   //TClonesArray definistions
   TClonesArray *utof_cal_r = 0;
   TClonesArray *utof_cal_l = 0;
   TClonesArray *dtof_cal_r = 0;
   TClonesArray *dtof_cal_l = 0;
   TClonesArray *ltof_cal_r = 0;
   TClonesArray *ltof_cal_l = 0;
   TClonesArray *t0_cal_l = 0;
   TClonesArray *t0_cal_r = 0;
   TClonesArray *t0r_cal_l = 0;
   TClonesArray *t0r_cal_r = 0;
   TClonesArray *bref_cal_l = 0;
   TClonesArray *bref_cal_r = 0;
   TClonesArray *t1_cal_l = 0;
   TClonesArray *t1_cal_r = 0;
   TClonesArray *bht_cal_l = 0;
   TClonesArray *bht_cal_r = 0;
   TClonesArray *bdc1_raw = 0;
   TClonesArray *bdc2_raw = 0;
   TClonesArray *kldc1_raw = 0;
   TClonesArray *kldc2_raw = 0;
   TClonesArray *bft1_raw = 0;
   TClonesArray *bft2_raw = 0;
   TClonesArray *sft1_raw = 0;
   TClonesArray *sft2_raw = 0;
   
   // TBranch *brutof_cal_r = 0;
   // TBranch *brutof_cal_l = 0;
   // TBranch *brdtof_cal_r = 0;
   // TBranch *brdtof_cal_l = 0;
   // TBranch *brltof_cal_r = 0;
   // TBranch *brltof_cal_l = 0;
   // TBranch *brt0_cal_l = 0;
   // TBranch *brt0_cal_r = 0;
   // TBranch *brt0r_cal_l = 0;
   // TBranch *brt0r_cal_r = 0;
   // TBranch *brbref_cal_l = 0;
   // TBranch *brbref_cal_r = 0;
   // TBranch *brt1_cal_l = 0;
   // TBranch *brt1_cal_r = 0;
   // // TBranch *bht_cal_l = 0;
   // // TBranch *bht_cal_r = 0;
   // TBranch *brbdc1_raw= 0;
   // TBranch *brbdc2_raw= 0;
   // TBranch *brkldc1_raw= 0;
   // TBranch *brkldc2_raw= 0;
   // TBranch *brbft1_raw= 0;
   // TBranch *brbft2_raw= 0;
   // TBranch *brsft1_raw= 0;
   // TBranch *brsft2_raw= 0;

   //Reading trees (tree1)
   tree1->SetBranchAddress("utof_cal_l", &utof_cal_l);
   tree1->SetBranchAddress("utof_cal_r", &utof_cal_r);
   tree1->SetBranchAddress("dtof_cal_l", &dtof_cal_l);
   tree1->SetBranchAddress("dtof_cal_r", &dtof_cal_r);
   tree1->SetBranchAddress("ltof_cal_l", &ltof_cal_l);
   tree1->SetBranchAddress("ltof_cal_r", &ltof_cal_r);
   tree1->SetBranchAddress("t0_cal_l", &t0_cal_l);
   tree1->SetBranchAddress("t0_cal_r", &t0_cal_r);
   tree1->SetBranchAddress("t0r_cal_l", &t0r_cal_l);
   tree1->SetBranchAddress("t0r_cal_r", &t0r_cal_r);
   tree1->SetBranchAddress("bref_cal_l", &bref_cal_l);
   tree1->SetBranchAddress("bref_cal_r", &bref_cal_r);
   tree1->SetBranchAddress("t1_cal_l", &t1_cal_l);
   tree1->SetBranchAddress("t1_cal_r", &t1_cal_r);
   tree1->SetBranchAddress("bht_cal_l", &bht_cal_l);
   tree1->SetBranchAddress("bht_cal_r", &bht_cal_r);
   tree1->SetBranchAddress("bdc1_raw", &bdc1_raw);
   tree1->SetBranchAddress("bdc2_raw", &bdc2_raw);
   tree1->SetBranchAddress("kldc1_raw", &kldc1_raw);
   tree1->SetBranchAddress("kldc2_raw", &kldc2_raw);
   tree1->SetBranchAddress("bft1_raw", &bft1_raw);
   tree1->SetBranchAddress("bft2_raw", &bft2_raw);
   tree1->SetBranchAddress("sft1_raw", &sft1_raw);
   tree1->SetBranchAddress("sft2_raw", &sft2_raw);

   //Event Tree Resize
   event.tot_utof_l.resize(NumOfSegUTOF);
   event.tot_utof_r.resize(NumOfSegUTOF);
   event.ltdc_utof_l.resize(NumOfSegUTOF);
   event.ltdc_utof_r.resize(NumOfSegUTOF);

   event.tot_dtof_l.resize(NumOfSegDTOF);
   event.tot_dtof_r.resize(NumOfSegDTOF);
   event.ltdc_dtof_l.resize(NumOfSegDTOF);
   event.ltdc_dtof_r.resize(NumOfSegDTOF);

   event.tot_ltof_l.resize(NumOfSegLTOF);
   event.tot_ltof_r.resize(NumOfSegLTOF);
   event.ltdc_ltof_l.resize(NumOfSegLTOF);
   event.ltdc_ltof_r.resize(NumOfSegLTOF);

   event.tot_t0_l.resize(NumOfSegT0);
   event.tot_t0_r.resize(NumOfSegT0);
   event.ltdc_t0_l.resize(NumOfSegT0);
   event.ltdc_t0_r.resize(NumOfSegT0);

   event.tot_t0r_l.resize(NumOfSegT0r);
   event.tot_t0r_r.resize(NumOfSegT0r);
   event.ltdc_t0r_l.resize(NumOfSegT0r);
   event.ltdc_t0r_r.resize(NumOfSegT0r);

   event.tot_bref_l.resize(NumOfSegBref);
   event.tot_bref_r.resize(NumOfSegBref);
   event.ltdc_bref_l.resize(NumOfSegBref);
   event.ltdc_bref_r.resize(NumOfSegBref);

   event.tot_t1_l.resize(NumOfSegT1);
   event.tot_t1_r.resize(NumOfSegT1);
   event.ltdc_t1_l.resize(NumOfSegT1);
   event.ltdc_t1_r.resize(NumOfSegT1);

   event.tot_bht_l.resize(NumOfSegBHT);
   event.tot_bht_r.resize(NumOfSegBHT);
   event.ltdc_bht_l.resize(NumOfSegBHT);
   event.ltdc_bht_r.resize(NumOfSegBHT);

   event.tot_bdc_l1.resize(NumOfWire);
   event.tot_bdc_l2.resize(NumOfWire);
   event.tot_bdc_l3.resize(NumOfWire);
   event.tot_bdc_l4.resize(NumOfWire);
   event.tot_bdc_l5.resize(NumOfWire);
   event.tot_bdc_l6.resize(NumOfWire);
   event.tot_bdc_l7.resize(NumOfWire);
   event.tot_bdc_l8.resize(NumOfWire);

   event.ltdc_bdc_l1.resize(NumOfWire);
   event.ltdc_bdc_l2.resize(NumOfWire);
   event.ltdc_bdc_l3.resize(NumOfWire);
   event.ltdc_bdc_l4.resize(NumOfWire);
   event.ltdc_bdc_l5.resize(NumOfWire);
   event.ltdc_bdc_l6.resize(NumOfWire);
   event.ltdc_bdc_l7.resize(NumOfWire);
   event.ltdc_bdc_l8.resize(NumOfWire);

   event.tot_kldc_l1.resize(NumOfWire);
   event.tot_kldc_l2.resize(NumOfWire);
   event.tot_kldc_l3.resize(NumOfWire);
   event.tot_kldc_l4.resize(NumOfWire);
   event.tot_kldc_l5.resize(NumOfWire);
   event.tot_kldc_l6.resize(NumOfWire);
   event.tot_kldc_l7.resize(NumOfWire);
   event.tot_kldc_l8.resize(NumOfWire);

   event.ltdc_kldc_l1.resize(NumOfWire);
   event.ltdc_kldc_l2.resize(NumOfWire);
   event.ltdc_kldc_l3.resize(NumOfWire);
   event.ltdc_kldc_l4.resize(NumOfWire);
   event.ltdc_kldc_l5.resize(NumOfWire);
   event.ltdc_kldc_l6.resize(NumOfWire);
   event.ltdc_kldc_l7.resize(NumOfWire);
   event.ltdc_kldc_l8.resize(NumOfWire);

   event.tot_bft_l1.resize(NumOfFiberBFT);
   event.tot_bft_l2.resize(NumOfFiberBFT);
   event.tot_bft_l3.resize(NumOfFiberBFT);
   event.tot_bft_l4.resize(NumOfFiberBFT);
   event.tot_bft_l5.resize(NumOfFiberBFT);
   event.tot_bft_l6.resize(NumOfFiberBFT);

   event.ltdc_bft_l1.resize(NumOfFiberBFT);
   event.ltdc_bft_l2.resize(NumOfFiberBFT);
   event.ltdc_bft_l3.resize(NumOfFiberBFT);
   event.ltdc_bft_l4.resize(NumOfFiberBFT);
   event.ltdc_bft_l5.resize(NumOfFiberBFT);
   event.ltdc_bft_l6.resize(NumOfFiberBFT);
   
   event.tot_sft_l1.resize(NumOfFiberSFT_X);
   event.tot_sft_l2.resize(NumOfFiberSFT_UV);
   event.tot_sft_l3.resize(NumOfFiberSFT_UV);
   event.tot_sft_l4.resize(NumOfFiberSFT_UV);
   event.tot_sft_l5.resize(NumOfFiberSFT_UV);
   event.tot_sft_l6.resize(NumOfFiberSFT_X);

   event.ltdc_sft_l1.resize(NumOfFiberSFT_X);
   event.ltdc_sft_l2.resize(NumOfFiberSFT_UV);
   event.ltdc_sft_l3.resize(NumOfFiberSFT_UV);
   event.ltdc_sft_l4.resize(NumOfFiberSFT_UV);
   event.ltdc_sft_l5.resize(NumOfFiberSFT_UV);
   event.ltdc_sft_l6.resize(NumOfFiberSFT_X);

   //Event building start
   cout<< "******************************" <<endl;
   const int n = tree1->GetEntries();

   int nhbf=0;
   if( nevent < n ) nhbf = nevent;
   if( nevent > n || nevent == -1 ) nhbf = n;

   cout<< "# of HBF = " << n << " : " << nhbf <<endl;
   
   for(int i=0; i<nhbf; i++){
      tree1->GetEntry(i);

      if(i%100==0) std::cout<< n << " : " << i <<std::endl;

      //Detector # of hits per HBF
      int nh_utof_l = utof_cal_l->GetEntries();
      int nh_utof_r = utof_cal_r->GetEntries();
      int nh_dtof_l = dtof_cal_l->GetEntries();
      int nh_dtof_r = dtof_cal_r->GetEntries();
      int nh_ltof_l = ltof_cal_l->GetEntries();
      int nh_ltof_r = ltof_cal_r->GetEntries();
      int nh_t0_l = t0_cal_l->GetEntries();
      int nh_t0_r = t0_cal_r->GetEntries();
      int nh_t0r_l = t0r_cal_l->GetEntries();
      int nh_t0r_r = t0r_cal_r->GetEntries();
      int nh_bref_l = bref_cal_l->GetEntries();
      int nh_bref_r = bref_cal_r->GetEntries();
      int nh_t1_l = t1_cal_l->GetEntries();
      int nh_t1_r = t1_cal_r->GetEntries();
      int nh_bht_l = bht_cal_l->GetEntries();
      int nh_bht_r = bht_cal_r->GetEntries();
      
      int nh_bdc1 = bdc1_raw->GetEntries();
      int nh_bdc2 = bdc2_raw->GetEntries();
      int nh_kldc1 = kldc1_raw->GetEntries();
      int nh_kldc2 = kldc2_raw->GetEntries();

      int nh_bft1 = bft1_raw->GetEntries();
      int nh_bft2 = bft2_raw->GetEntries();
      int nh_sft1 = sft1_raw->GetEntries();
      int nh_sft2 = sft2_raw->GetEntries();
      
      
      double ref_t_prev = -10000.0; //Previous reference timing
      
      //cout << "Events per HBF: " << nh_utof_l << endl;
      for( int j = 0; j<nh_utof_l; j++ ){
         InitializeEvent();
      
         double ref_t; //Reference timing (UTOF L timing)
         double offset = 1000.0; //TDC offset      

         //double cut  = 2000.0; //For Timing counter (including Bref and BFT)
         //double cut2 = 2000.0;//For Drift Chamber

         

         double cut  = 300.0; //For Timing counter
         double cut2 = 1000.0;//For Drift Chamber

         //Timing counter
         art::TTimingChargeData* data_utof_l = (art::TTimingChargeData*)utof_cal_l->UncheckedAt(j);
         int utof_l_id = data_utof_l->GetDetID();
         double utof_l_tot = data_utof_l->GetCharge();
         double utof_l_ltdc = data_utof_l->GetTiming();
         
         
         // ignore the UTOF hits when timing was close to the previous one
         if(fabs(utof_l_ltdc - ref_t_prev) < cut){ 
            continue;
         }
         else{
            ref_t = utof_l_ltdc;
            ref_t_prev = ref_t;
         }
         

         //UTOF L
         event.id_utof_l.push_back(utof_l_id+1);
         event.tot_utof_l[utof_l_id].push_back(utof_l_tot);
         event.ltdc_utof_l[utof_l_id].push_back(utof_l_ltdc-ref_t+offset);
         
         //UTOF R
         for( int k = 0; k<nh_utof_r; k++ ){
            art::TTimingChargeData* data_utof_r = (art::TTimingChargeData*)utof_cal_r->UncheckedAt(k);
            int utof_r_id = data_utof_r->GetDetID();
            double utof_r_tot = data_utof_r->GetCharge();
            double utof_r_ltdc = data_utof_r->GetTiming();
            if( fabs(utof_r_ltdc-ref_t)<cut ){
               event.id_utof_r.push_back(utof_r_id+1);
               event.tot_utof_r[utof_r_id].push_back(utof_r_tot);
               event.ltdc_utof_r[utof_r_id].push_back(utof_r_ltdc-ref_t+offset);
            }
         }

         //DTOF L
         for( int k = 0; k<nh_dtof_l; k++ ){
            art::TTimingChargeData* data_dtof_l = (art::TTimingChargeData*)dtof_cal_l->UncheckedAt(k);
            int dtof_l_id = data_dtof_l->GetDetID();
            double dtof_l_tot = data_dtof_l->GetCharge();
            double dtof_l_ltdc = data_dtof_l->GetTiming();
            if( fabs(dtof_l_ltdc-ref_t)<cut ){
               event.id_dtof_l.push_back(dtof_l_id+1);
               event.tot_dtof_l[dtof_l_id].push_back(dtof_l_tot);
               event.ltdc_dtof_l[dtof_l_id].push_back(dtof_l_ltdc-ref_t+offset);
            }
         }
         //DTOF R
         for( int k = 0; k<nh_dtof_r; k++ ){
            art::TTimingChargeData* data_dtof_r = (art::TTimingChargeData*)dtof_cal_r->UncheckedAt(k);
            int dtof_r_id = data_dtof_r->GetDetID();
            double dtof_r_tot = data_dtof_r->GetCharge();
            double dtof_r_ltdc = data_dtof_r->GetTiming();
            if( fabs(dtof_r_ltdc-ref_t)<cut ){
               event.id_dtof_r.push_back(dtof_r_id+1);
               event.tot_dtof_r[dtof_r_id].push_back(dtof_r_tot);
               event.ltdc_dtof_r[dtof_r_id].push_back(dtof_r_ltdc-ref_t+offset);
            }
         }

         //LTOF L
         for( int k = 0; k<nh_ltof_l; k++ ){
            art::TTimingChargeData* data_ltof_l = (art::TTimingChargeData*)ltof_cal_l->UncheckedAt(k);
            int ltof_l_id = data_ltof_l->GetDetID();
            double ltof_l_tot = data_ltof_l->GetCharge();
            double ltof_l_ltdc = data_ltof_l->GetTiming();
            if( fabs(ltof_l_ltdc-ref_t)<cut ){
               event.id_ltof_l.push_back(ltof_l_id+1);
               event.tot_ltof_l[ltof_l_id].push_back(ltof_l_tot);
               event.ltdc_ltof_l[ltof_l_id].push_back(ltof_l_ltdc-ref_t+offset);
            }
         }
         //LTOF R
         for( int k = 0; k<nh_ltof_r; k++ ){
            art::TTimingChargeData* data_ltof_r = (art::TTimingChargeData*)ltof_cal_r->UncheckedAt(k);
            int ltof_r_id = data_ltof_r->GetDetID();
            double ltof_r_tot = data_ltof_r->GetCharge();
            double ltof_r_ltdc = data_ltof_r->GetTiming();
            if( fabs(ltof_r_ltdc-ref_t)<cut ){
               event.id_ltof_r.push_back(ltof_r_id+1);
               event.tot_ltof_r[ltof_r_id].push_back(ltof_r_tot);
               event.ltdc_ltof_r[ltof_r_id].push_back(ltof_r_ltdc-ref_t+offset);
            }
         }

         //T0 L
         for( int k = 0; k<nh_t0_l; k++ ){
            art::TTimingChargeData* data_t0_l = (art::TTimingChargeData*)t0_cal_l->UncheckedAt(k);
            int t0_l_id = data_t0_l->GetDetID();
            double t0_l_tot = data_t0_l->GetCharge();
            double t0_l_ltdc = data_t0_l->GetTiming();
            if( fabs(t0_l_ltdc-ref_t)<cut ){
               event.id_t0_l.push_back(t0_l_id+1);
               event.tot_t0_l[t0_l_id].push_back(t0_l_tot);
               event.ltdc_t0_l[t0_l_id].push_back(t0_l_ltdc-ref_t+offset);
            }
         }
         //T0 R
         for( int k = 0; k<nh_t0_r; k++ ){
            art::TTimingChargeData* data_t0_r = (art::TTimingChargeData*)t0_cal_r->UncheckedAt(k);
            int t0_r_id = data_t0_r->GetDetID();
            double t0_r_tot = data_t0_r->GetCharge();
            double t0_r_ltdc = data_t0_r->GetTiming();
            if( fabs(t0_r_ltdc-ref_t)<cut ){
               event.id_t0_r.push_back(t0_r_id+1);
               event.tot_t0_r[t0_r_id].push_back(t0_r_tot);
               event.ltdc_t0_r[t0_r_id].push_back(t0_r_ltdc-ref_t+offset);
            }
         }

         //T0R L
         for( int k = 0; k<nh_t0r_l; k++ ){
            art::TTimingChargeData* data_t0r_l = (art::TTimingChargeData*)t0r_cal_l->UncheckedAt(k);
            int t0r_l_id = data_t0r_l->GetDetID();
            double t0r_l_tot = data_t0r_l->GetCharge();
            double t0r_l_ltdc = data_t0r_l->GetTiming();
            if( fabs(t0r_l_ltdc-ref_t)<cut ){
               event.id_t0r_l.push_back(t0r_l_id+1);
               event.tot_t0r_l[t0r_l_id].push_back(t0r_l_tot);
               event.ltdc_t0r_l[t0r_l_id].push_back(t0r_l_ltdc-ref_t+offset);
            }
         }
         //T0R R
         for( int k = 0; k<nh_t0r_r; k++ ){
            art::TTimingChargeData* data_t0r_r = (art::TTimingChargeData*)t0r_cal_r->UncheckedAt(k);
            int t0r_r_id = data_t0r_r->GetDetID();
            double t0r_r_tot = data_t0r_r->GetCharge();
            double t0r_r_ltdc = data_t0r_r->GetTiming();
            if( fabs(t0r_r_ltdc-ref_t)<cut ){
               event.id_t0r_r.push_back(t0r_r_id+1);
               event.tot_t0r_r[t0r_r_id].push_back(t0r_r_tot);
               event.ltdc_t0r_r[t0r_r_id].push_back(t0r_r_ltdc-ref_t+offset);
            }
         }
         
         //Bref L
         for( int k = 0; k<nh_bref_l; k++ ){
            art::TTimingChargeData* data_bref_l = (art::TTimingChargeData*)bref_cal_l->UncheckedAt(k);
            int bref_l_id = data_bref_l->GetDetID();
            double bref_l_tot = data_bref_l->GetCharge();
            double bref_l_ltdc = data_bref_l->GetTiming();
            if( fabs(bref_l_ltdc-ref_t)<cut ){
               event.id_bref_l.push_back(bref_l_id+1);
               event.tot_bref_l[bref_l_id].push_back(bref_l_tot);
               event.ltdc_bref_l[bref_l_id].push_back(bref_l_ltdc-ref_t+offset);
            }
         }
         //Bref R
         for( int k = 0; k<nh_bref_r; k++ ){
            art::TTimingChargeData* data_bref_r = (art::TTimingChargeData*)bref_cal_r->UncheckedAt(k);
            int bref_r_id = data_bref_r->GetDetID();
            double bref_r_tot = data_bref_r->GetCharge();
            double bref_r_ltdc = data_bref_r->GetTiming();
            if( fabs(bref_r_ltdc-ref_t)<cut ){
               event.id_bref_r.push_back(bref_r_id+1);
               event.tot_bref_r[bref_r_id].push_back(bref_r_tot);
               event.ltdc_bref_r[bref_r_id].push_back(bref_r_ltdc-ref_t+offset);
            }
         }

         //T1 L
         for( int k = 0; k<nh_t1_l; k++ ){
            art::TTimingChargeData* data_t1_l = (art::TTimingChargeData*)t1_cal_l->UncheckedAt(k);
            int t1_l_id = data_t1_l->GetDetID();
            double t1_l_tot = data_t1_l->GetCharge();
            double t1_l_ltdc = data_t1_l->GetTiming();
            if( fabs(t1_l_ltdc-ref_t)<cut ){
               event.id_t1_l.push_back(t1_l_id+1);
               event.tot_t1_l[t1_l_id].push_back(t1_l_tot);
               event.ltdc_t1_l[t1_l_id].push_back(t1_l_ltdc-ref_t+offset);
            }
         }
         //T1 R
         for( int k = 0; k<nh_t1_r; k++ ){
            art::TTimingChargeData* data_t1_r = (art::TTimingChargeData*)t1_cal_r->UncheckedAt(k);
            int t1_r_id = data_t1_r->GetDetID();
            double t1_r_tot = data_t1_r->GetCharge();
            double t1_r_ltdc = data_t1_r->GetTiming();
            if( fabs(t1_r_ltdc-ref_t)<cut ){
               event.id_t1_r.push_back(t1_r_id+1);
               event.tot_t1_r[t1_r_id].push_back(t1_r_tot);
               event.ltdc_t1_r[t1_r_id].push_back(t1_r_ltdc-ref_t+offset);
            }
         }

         //BHT L
         for( int k = 0; k<nh_bht_l; k++ ){
            art::TTimingChargeData* data_bht_l = (art::TTimingChargeData*)bht_cal_l->UncheckedAt(k);
            int bht_l_id = data_bht_l->GetDetID();
            double bht_l_tot = data_bht_l->GetCharge();
            double bht_l_ltdc = data_bht_l->GetTiming();
            if( fabs(bht_l_ltdc-ref_t)<cut ){
               event.id_bht_l.push_back(bht_l_id+1);
               event.tot_bht_l[bht_l_id].push_back(bht_l_tot);
               event.ltdc_bht_l[bht_l_id].push_back(bht_l_ltdc-ref_t+offset);
            }
         }
         //BHT R
         for( int k = 0; k<nh_bht_r; k++ ){
            art::TTimingChargeData* data_bht_r = (art::TTimingChargeData*)bht_cal_r->UncheckedAt(k);
            int bht_r_id = data_bht_r->GetDetID();
            double bht_r_tot = data_bht_r->GetCharge();
            double bht_r_ltdc = data_bht_r->GetTiming();
            if( fabs(bht_r_ltdc-ref_t)<cut ){
               event.id_bht_r.push_back(bht_r_id+1);
               event.tot_bht_r[bht_r_id].push_back(bht_r_tot);
               event.ltdc_bht_r[bht_r_id].push_back(bht_r_ltdc-ref_t+offset);
            }
         }

         //cout<< "*****" << endl;
         //BDC1
         for( int k = 0; k<nh_bdc1; k++ ){
            art::TTimingChargeData* data_bdc1 = (art::TTimingChargeData*)bdc1_raw->UncheckedAt(k);
            int bdc1_id = data_bdc1->GetDetID();
            double bdc1_tot = data_bdc1->GetCharge();
            double bdc1_ltdc = data_bdc1->GetTiming();

            //BDC ch offset
            int ofs1 = 112, ofs2 = 122;
            int ofs3 = 228, ofs4 = 237;
            int ofs5 = 342, ofs6 = 341;
            int ofs_wire1 = 116;
            int ofs_wire2 = 232;
            int ofs_wire3 = 341;

            if( fabs(bdc1_ltdc-ref_t)<cut2 ){
               if(bdc1_id<ofs1){
                  event.id_bdc_l1.push_back(bdc1_id);
                  event.tot_bdc_l1[bdc1_id].push_back(bdc1_tot);
                  event.ltdc_bdc_l1[bdc1_id].push_back(bdc1_ltdc-ref_t+offset);
               }
               if(ofs2<bdc1_id&&bdc1_id<ofs3){
                  event.id_bdc_l2.push_back(bdc1_id-ofs_wire1);
                  event.tot_bdc_l2[bdc1_id-ofs_wire1].push_back(bdc1_tot);
                  event.ltdc_bdc_l2[bdc1_id-ofs_wire1].push_back(bdc1_ltdc-ref_t+offset);
               }
               if(ofs4<bdc1_id&&bdc1_id<ofs5){
                  event.id_bdc_l3.push_back(bdc1_id-ofs_wire2);
                  event.tot_bdc_l3[bdc1_id-ofs_wire2].push_back(bdc1_tot);
                  event.ltdc_bdc_l3[bdc1_id-ofs_wire2].push_back(bdc1_ltdc-ref_t+offset);
               }
               if(ofs6<bdc1_id){
                  event.id_bdc_l4.push_back(bdc1_id-ofs_wire3);
                  event.tot_bdc_l4[bdc1_id-ofs_wire3].push_back(bdc1_tot);
                  event.ltdc_bdc_l4[bdc1_id-ofs_wire3].push_back(bdc1_ltdc-ref_t+offset);
               }
            }
         }
         //BDC2
         for( int k = 0; k<nh_bdc2; k++ ){
            art::TTimingChargeData* data_bdc2 = (art::TTimingChargeData*)bdc2_raw->UncheckedAt(k);
            int bdc2_id = data_bdc2->GetDetID();
            double bdc2_tot = data_bdc2->GetCharge();
            double bdc2_ltdc = data_bdc2->GetTiming();

            //BDC ch offset
            int ofs1 = 112, ofs2 = 122;
            int ofs3 = 228, ofs4 = 237;
            int ofs5 = 342, ofs6 = 341;
            int ofs_wire1 = 116;
            int ofs_wire2 = 232;
            int ofs_wire3 = 341;

            if( fabs(bdc2_ltdc-ref_t)<cut2 ){
               if(bdc2_id<ofs1){
                  event.id_bdc_l5.push_back(bdc2_id);
                  event.tot_bdc_l5[bdc2_id].push_back(bdc2_tot);
                  event.ltdc_bdc_l5[bdc2_id].push_back(bdc2_ltdc-ref_t+offset);
               }
               if(ofs2<bdc2_id&&bdc2_id<ofs3){
                  event.id_bdc_l6.push_back(bdc2_id-ofs_wire1);
                  event.tot_bdc_l6[bdc2_id-ofs_wire1].push_back(bdc2_tot);
                  event.ltdc_bdc_l6[bdc2_id-ofs_wire1].push_back(bdc2_ltdc-ref_t+offset);
               }
               if(ofs4<bdc2_id&&bdc2_id<ofs5){
                  event.id_bdc_l7.push_back(bdc2_id-ofs_wire2);
                  event.tot_bdc_l7[bdc2_id-ofs_wire2].push_back(bdc2_tot);
                  event.ltdc_bdc_l7[bdc2_id-ofs_wire2].push_back(bdc2_ltdc-ref_t+offset);
               }
               if(ofs6<bdc2_id){
                  event.id_bdc_l8.push_back(bdc2_id-ofs_wire3);
                  event.tot_bdc_l8[bdc2_id-ofs_wire3].push_back(bdc2_tot);
                  event.ltdc_bdc_l8[bdc2_id-ofs_wire3].push_back(bdc2_ltdc-ref_t+offset);
               }
            }
         }

         //cout<< "*****" << endl;
         //KLDC1
         for( int k = 0; k<nh_kldc1; k++ ){
            art::TTimingChargeData* data_kldc1 = (art::TTimingChargeData*)kldc1_raw->UncheckedAt(k);
            int kldc1_id = data_kldc1->GetDetID();
            double kldc1_tot = data_kldc1->GetCharge();
            double kldc1_ltdc = data_kldc1->GetTiming();

            //KLDC ch offset
            int ofs_wire1 = 128;
            int ofs_wire2 = 256;
            int ofs_wire3 = 384;

            if( fabs(kldc1_ltdc-ref_t)<cut2 ){
               if(kldc1_id<ofs_wire1){
                  event.id_kldc_l1.push_back(kldc1_id);
                  event.tot_kldc_l1[kldc1_id].push_back(kldc1_tot);
                  event.ltdc_kldc_l1[kldc1_id].push_back(kldc1_ltdc-ref_t+offset);
               }
               if(ofs_wire1<kldc1_id&&kldc1_id<ofs_wire2){
                  event.id_kldc_l2.push_back(kldc1_id-ofs_wire1);
                  event.tot_kldc_l2[kldc1_id-ofs_wire1].push_back(kldc1_tot);
                  event.ltdc_kldc_l2[kldc1_id-ofs_wire1].push_back(kldc1_ltdc-ref_t+offset);
               }
               if(ofs_wire2<kldc1_id&&kldc1_id<ofs_wire3){
                  event.id_kldc_l3.push_back(kldc1_id-ofs_wire2);
                  event.tot_kldc_l3[kldc1_id-ofs_wire2].push_back(kldc1_tot);
                  event.ltdc_kldc_l3[kldc1_id-ofs_wire2].push_back(kldc1_ltdc-ref_t+offset);
               }
               if(ofs_wire3<kldc1_id){
                  event.id_kldc_l4.push_back(kldc1_id-ofs_wire3);
                  event.tot_kldc_l4[kldc1_id-ofs_wire3].push_back(kldc1_tot);
                  event.ltdc_kldc_l4[kldc1_id-ofs_wire3].push_back(kldc1_ltdc-ref_t+offset);
               }
            }
         }
         //KLDC2
         for( int k = 0; k<nh_kldc2; k++ ){
            art::TTimingChargeData* data_kldc2 = (art::TTimingChargeData*)kldc2_raw->UncheckedAt(k);
            int kldc2_id = data_kldc2->GetDetID();
            double kldc2_tot = data_kldc2->GetCharge();
            double kldc2_ltdc = data_kldc2->GetTiming();

            //KLDC ch offset
            int ofs_wire1 = 128;
            int ofs_wire2 = 256;
            int ofs_wire3 = 384;

            if( fabs(kldc2_ltdc-ref_t)<cut2 ){
               if(kldc2_id<ofs_wire1){
                  event.id_kldc_l5.push_back(kldc2_id);
                  event.tot_kldc_l5[kldc2_id].push_back(kldc2_tot);
                  event.ltdc_kldc_l5[kldc2_id].push_back(kldc2_ltdc-ref_t+offset);
               }
               if(ofs_wire1<kldc2_id&&kldc2_id<ofs_wire2){
                  event.id_kldc_l6.push_back(kldc2_id-ofs_wire1);
                  event.tot_kldc_l6[kldc2_id-ofs_wire1].push_back(kldc2_tot);
                  event.ltdc_kldc_l6[kldc2_id-ofs_wire1].push_back(kldc2_ltdc-ref_t+offset);
               }
               if(ofs_wire2<kldc2_id&&kldc2_id<ofs_wire3){
                  event.id_kldc_l7.push_back(kldc2_id-ofs_wire2);
                  event.tot_kldc_l7[kldc2_id-ofs_wire2].push_back(kldc2_tot);
                  event.ltdc_kldc_l7[kldc2_id-ofs_wire2].push_back(kldc2_ltdc-ref_t+offset);
               }
               if(ofs_wire3<kldc2_id){
                  event.id_kldc_l8.push_back(kldc2_id-ofs_wire3);
                  event.tot_kldc_l8[kldc2_id-ofs_wire3].push_back(kldc2_tot);
                  event.ltdc_kldc_l8[kldc2_id-ofs_wire3].push_back(kldc2_ltdc-ref_t+offset);
               }
            }
         }

         //cout<< "*****" << endl;
         //BFT1
         for( int k = 0; k<nh_bft1; k++ ){
            art::TTimingChargeData* data_bft1 = (art::TTimingChargeData*)bft1_raw->UncheckedAt(k);
            int bft1_id = data_bft1->GetDetID();
            double bft1_tot = data_bft1->GetCharge();
            double bft1_ltdc = data_bft1->GetTiming();

            //BFT ch offset
            int ofs_fiber1 = 301;
            int ofs_fiber2 = 601;
            
            if( fabs(bft1_ltdc-ref_t)<cut ){
               //if( fabs(bft1_ltdc-ref_t)<cut2 ){
               if(bft1_id<ofs_fiber1){
                  event.id_bft_l1.push_back(bft1_id-1);
                  event.tot_bft_l1[bft1_id-1].push_back(bft1_tot);
                  event.ltdc_bft_l1[bft1_id-1].push_back(bft1_ltdc-ref_t+offset);
               }
               if(ofs_fiber1<bft1_id&&bft1_id<ofs_fiber2){
                  event.id_bft_l2.push_back(bft1_id-ofs_fiber1);
                  event.tot_bft_l2[bft1_id-ofs_fiber1].push_back(bft1_tot);
                  event.ltdc_bft_l2[bft1_id-ofs_fiber1].push_back(bft1_ltdc-ref_t+offset);
               }
               if(ofs_fiber2<bft1_id){
                  event.id_bft_l3.push_back(bft1_id-ofs_fiber2);
                  event.tot_bft_l3[bft1_id-ofs_fiber2].push_back(bft1_tot);
                  event.ltdc_bft_l3[bft1_id-ofs_fiber2].push_back(bft1_ltdc-ref_t+offset);
               }
            }
         }
         //BFT2
         for( int k = 0; k<nh_bft2; k++ ){
            art::TTimingChargeData* data_bft2 = (art::TTimingChargeData*)bft2_raw->UncheckedAt(k);
            int bft2_id = data_bft2->GetDetID();
            double bft2_tot = data_bft2->GetCharge();
            double bft2_ltdc = data_bft2->GetTiming();

            //BFT ch offset
            int ofs_fiber1 = 301;
            int ofs_fiber2 = 601;

            if( fabs(bft2_ltdc-ref_t)<cut2 ){
               if(bft2_id<ofs_fiber1){
                  event.id_bft_l4.push_back(bft2_id-1);
                  event.tot_bft_l4[bft2_id-1].push_back(bft2_tot);
                  event.ltdc_bft_l4[bft2_id-1].push_back(bft2_ltdc-ref_t+offset);
               }
               if(ofs_fiber1<bft2_id&&bft2_id<ofs_fiber2){
                  event.id_bft_l5.push_back(bft2_id-ofs_fiber1);
                  event.tot_bft_l5[bft2_id-ofs_fiber1].push_back(bft2_tot);
                  event.ltdc_bft_l5[bft2_id-ofs_fiber1].push_back(bft2_ltdc-ref_t+offset);
               }
               if(ofs_fiber2<bft2_id){
                  event.id_bft_l6.push_back(bft2_id-ofs_fiber2);
                  event.tot_bft_l6[bft2_id-ofs_fiber2].push_back(bft2_tot);
                  event.ltdc_bft_l6[bft2_id-ofs_fiber2].push_back(bft2_ltdc-ref_t+offset);
               }
            }
         }

         //cout<< "*****" << endl;
         //SFT1
         for( int k = 0; k<nh_sft1; k++ ){
            art::TTimingChargeData* data_sft1 = (art::TTimingChargeData*)sft1_raw->UncheckedAt(k);
            int sft1_id = data_sft1->GetDetID();
            double sft1_tot = data_sft1->GetCharge();
            double sft1_ltdc = data_sft1->GetTiming();

            //SFT ch offset
            int ofs1 =  850;
            int ofs2 =  900, ofs3 = 1050;
            int ofs4 = 1100;

            int ofs_fiber1 = 401;
            int ofs_fiber2 = 601;
            
            if( fabs(sft1_ltdc-ref_t)<cut ){
               if(sft1_id<ofs1){
                  event.id_sft_l1.push_back(sft1_id);
                  event.tot_sft_l1[sft1_id].push_back(sft1_tot);
                  event.ltdc_sft_l1[sft1_id].push_back(sft1_ltdc-ref_t+offset);
               }
               if(ofs2<sft1_id&&sft1_id<ofs3){
                  event.id_sft_l2.push_back(sft1_id-ofs_fiber1);
                  event.tot_sft_l2[sft1_id-ofs_fiber1].push_back(sft1_tot);
                  event.ltdc_sft_l2[sft1_id-ofs_fiber1].push_back(sft1_ltdc-ref_t+offset);
               }
               if(ofs4<sft1_id){
                  event.id_sft_l3.push_back(sft1_id-ofs_fiber2);
                  event.tot_sft_l3[sft1_id-ofs_fiber2].push_back(sft1_tot);
                  event.ltdc_sft_l3[sft1_id-ofs_fiber2].push_back(sft1_ltdc-ref_t+offset);
               }
            }
         }
         //SFT2
         for( int k = 0; k<nh_sft2; k++ ){
            art::TTimingChargeData* data_sft2 = (art::TTimingChargeData*)sft2_raw->UncheckedAt(k);
            int sft2_id = data_sft2->GetDetID();
            double sft2_tot = data_sft2->GetCharge();
            double sft2_ltdc = data_sft2->GetTiming();

            //SFT ch offset
            int ofs1 =  650;
            int ofs2 =  900, ofs3 = 1050;
            int ofs4 = 1100;

            int ofs_fiber1 = 401;
            int ofs_fiber2 = 601;

            if( fabs(sft2_ltdc-ref_t)<cut2 ){
               if(sft2_id<ofs1){
                  event.id_sft_l6.push_back(sft2_id);
                  event.tot_sft_l6[sft2_id].push_back(sft2_tot);
                  event.ltdc_sft_l6[sft2_id].push_back(sft2_ltdc-ref_t+offset);
               }
               if(ofs2<sft2_id&&sft2_id<ofs3){
                  event.id_sft_l5.push_back(sft2_id-ofs_fiber1);
                  event.tot_sft_l5[sft2_id-ofs_fiber1].push_back(sft2_tot);
                  event.ltdc_sft_l5[sft2_id-ofs_fiber1].push_back(sft2_ltdc-ref_t+offset);
               }
               if(ofs4<sft2_id){
                  event.id_sft_l4.push_back(sft2_id-ofs_fiber2);
                  event.tot_sft_l4[sft2_id-ofs_fiber2].push_back(sft2_tot);
                  event.ltdc_sft_l4[sft2_id-ofs_fiber2].push_back(sft2_ltdc-ref_t+offset);
               }
            }
         }

         tree->Fill();
      }
   }
   gFile->Write();

   return 0;
}

void InitializeEvent( void )
{
   //////Timing counter
   //UTOF
   event.id_utof_l.clear();
   event.id_utof_r.clear();
   for(int i=0; i<event.tot_utof_l.size(); i++) event.tot_utof_l[i].clear();
   for(int i=0; i<event.tot_utof_r.size(); i++) event.tot_utof_r[i].clear();
   for(int i=0; i<event.ltdc_utof_l.size(); i++) event.ltdc_utof_l[i].clear();
   for(int i=0; i<event.ltdc_utof_r.size(); i++) event.ltdc_utof_r[i].clear();

   //DTOF
   event.id_dtof_l.clear();
   event.id_dtof_r.clear();
   for(int i=0; i<event.tot_dtof_l.size(); i++) event.tot_dtof_l[i].clear();
   for(int i=0; i<event.tot_dtof_r.size(); i++) event.tot_dtof_r[i].clear();
   for(int i=0; i<event.ltdc_dtof_l.size(); i++) event.ltdc_dtof_l[i].clear();
   for(int i=0; i<event.ltdc_dtof_r.size(); i++) event.ltdc_dtof_r[i].clear();

   //LTOF
   event.id_ltof_l.clear();
   event.id_ltof_r.clear();
   for(int i=0; i<event.tot_ltof_l.size(); i++) event.tot_ltof_l[i].clear();
   for(int i=0; i<event.tot_ltof_r.size(); i++) event.tot_ltof_r[i].clear();
   for(int i=0; i<event.ltdc_ltof_l.size(); i++) event.ltdc_ltof_l[i].clear();
   for(int i=0; i<event.ltdc_ltof_r.size(); i++) event.ltdc_ltof_r[i].clear();

   //T0
   event.id_t0_l.clear();
   event.id_t0_r.clear();
   for(int i=0; i<event.tot_t0_l.size(); i++) event.tot_t0_l[i].clear();
   for(int i=0; i<event.tot_t0_r.size(); i++) event.tot_t0_r[i].clear();
   for(int i=0; i<event.ltdc_t0_l.size(); i++) event.ltdc_t0_l[i].clear();
   for(int i=0; i<event.ltdc_t0_r.size(); i++) event.ltdc_t0_r[i].clear();

   //T0R
   event.id_t0r_l.clear();
   event.id_t0r_r.clear();
   for(int i=0; i<event.tot_t0r_l.size(); i++) event.tot_t0r_l[i].clear();
   for(int i=0; i<event.tot_t0r_r.size(); i++) event.tot_t0r_r[i].clear();
   for(int i=0; i<event.ltdc_t0r_l.size(); i++) event.ltdc_t0r_l[i].clear();
   for(int i=0; i<event.ltdc_t0r_r.size(); i++) event.ltdc_t0r_r[i].clear();

   //Bref
   event.id_bref_l.clear();
   event.id_bref_r.clear();
   for(int i=0; i<event.tot_bref_l.size(); i++) event.tot_bref_l[i].clear();
   for(int i=0; i<event.tot_bref_r.size(); i++) event.tot_bref_r[i].clear();
   for(int i=0; i<event.ltdc_bref_l.size(); i++) event.ltdc_bref_l[i].clear();
   for(int i=0; i<event.ltdc_bref_r.size(); i++) event.ltdc_bref_r[i].clear();

   //T1 from K1.8BR   
   event.id_t1_l.clear();
   event.id_t1_r.clear();
   for(int i=0; i<event.tot_t1_l.size(); i++) event.tot_t1_l[i].clear();
   for(int i=0; i<event.tot_t1_r.size(); i++) event.tot_t1_r[i].clear();
   for(int i=0; i<event.ltdc_t1_l.size(); i++) event.ltdc_t1_l[i].clear();
   for(int i=0; i<event.ltdc_t1_r.size(); i++) event.ltdc_t1_r[i].clear();

   //BHT from K1.8BR   
   event.id_bht_l.clear();
   event.id_bht_r.clear();
   for(int i=0; i<event.tot_bht_l.size(); i++) event.tot_bht_l[i].clear();
   for(int i=0; i<event.tot_bht_r.size(); i++) event.tot_bht_r[i].clear();
   for(int i=0; i<event.ltdc_bht_l.size(); i++) event.ltdc_bht_l[i].clear();
   for(int i=0; i<event.ltdc_bht_r.size(); i++) event.ltdc_bht_r[i].clear();

   ////Drift chamber
   //BDC
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

   //KLDC
   event.id_kldc_l1.clear();
   event.id_kldc_l2.clear();
   event.id_kldc_l3.clear();
   event.id_kldc_l4.clear();
   event.id_kldc_l5.clear();
   event.id_kldc_l6.clear();
   event.id_kldc_l7.clear();
   event.id_kldc_l8.clear();
   for(int i=0; i<event.tot_kldc_l1.size(); i++) event.tot_kldc_l1[i].clear();
   for(int i=0; i<event.tot_kldc_l2.size(); i++) event.tot_kldc_l2[i].clear();
   for(int i=0; i<event.tot_kldc_l3.size(); i++) event.tot_kldc_l3[i].clear();
   for(int i=0; i<event.tot_kldc_l4.size(); i++) event.tot_kldc_l4[i].clear();
   for(int i=0; i<event.tot_kldc_l5.size(); i++) event.tot_kldc_l5[i].clear();
   for(int i=0; i<event.tot_kldc_l6.size(); i++) event.tot_kldc_l6[i].clear();
   for(int i=0; i<event.tot_kldc_l7.size(); i++) event.tot_kldc_l7[i].clear();
   for(int i=0; i<event.tot_kldc_l8.size(); i++) event.tot_kldc_l8[i].clear();
   for(int i=0; i<event.ltdc_kldc_l1.size(); i++) event.ltdc_kldc_l1[i].clear();
   for(int i=0; i<event.ltdc_kldc_l2.size(); i++) event.ltdc_kldc_l2[i].clear();
   for(int i=0; i<event.ltdc_kldc_l3.size(); i++) event.ltdc_kldc_l3[i].clear();
   for(int i=0; i<event.ltdc_kldc_l4.size(); i++) event.ltdc_kldc_l4[i].clear();
   for(int i=0; i<event.ltdc_kldc_l5.size(); i++) event.ltdc_kldc_l5[i].clear();
   for(int i=0; i<event.ltdc_kldc_l6.size(); i++) event.ltdc_kldc_l6[i].clear();
   for(int i=0; i<event.ltdc_kldc_l7.size(); i++) event.ltdc_kldc_l7[i].clear();
   for(int i=0; i<event.ltdc_kldc_l8.size(); i++) event.ltdc_kldc_l8[i].clear();

   ////Fiber tracker
   //BFT
   event.id_bft_l1.clear();
   event.id_bft_l2.clear();
   event.id_bft_l3.clear();
   event.id_bft_l4.clear();
   event.id_bft_l5.clear();
   event.id_bft_l6.clear();
   for(int i=0; i<event.tot_bft_l1.size(); i++) event.tot_bft_l1[i].clear();
   for(int i=0; i<event.tot_bft_l2.size(); i++) event.tot_bft_l2[i].clear();
   for(int i=0; i<event.tot_bft_l3.size(); i++) event.tot_bft_l3[i].clear();
   for(int i=0; i<event.tot_bft_l4.size(); i++) event.tot_bft_l4[i].clear();
   for(int i=0; i<event.tot_bft_l5.size(); i++) event.tot_bft_l5[i].clear();
   for(int i=0; i<event.tot_bft_l6.size(); i++) event.tot_bft_l6[i].clear();
   for(int i=0; i<event.ltdc_bft_l1.size(); i++) event.ltdc_bft_l1[i].clear();
   for(int i=0; i<event.ltdc_bft_l2.size(); i++) event.ltdc_bft_l2[i].clear();
   for(int i=0; i<event.ltdc_bft_l3.size(); i++) event.ltdc_bft_l3[i].clear();
   for(int i=0; i<event.ltdc_bft_l4.size(); i++) event.ltdc_bft_l4[i].clear();
   for(int i=0; i<event.ltdc_bft_l5.size(); i++) event.ltdc_bft_l5[i].clear();
   for(int i=0; i<event.ltdc_bft_l6.size(); i++) event.ltdc_bft_l6[i].clear();

   //SFT
   event.id_sft_l1.clear();
   event.id_sft_l2.clear();
   event.id_sft_l3.clear();
   event.id_sft_l4.clear();
   event.id_sft_l5.clear();
   event.id_sft_l6.clear();
   for(int i=0; i<event.tot_sft_l1.size(); i++) event.tot_sft_l1[i].clear();
   for(int i=0; i<event.tot_sft_l2.size(); i++) event.tot_sft_l2[i].clear();
   for(int i=0; i<event.tot_sft_l3.size(); i++) event.tot_sft_l3[i].clear();
   for(int i=0; i<event.tot_sft_l4.size(); i++) event.tot_sft_l4[i].clear();
   for(int i=0; i<event.tot_sft_l5.size(); i++) event.tot_sft_l5[i].clear();
   for(int i=0; i<event.tot_sft_l6.size(); i++) event.tot_sft_l6[i].clear();
   for(int i=0; i<event.ltdc_sft_l1.size(); i++) event.ltdc_sft_l1[i].clear();
   for(int i=0; i<event.ltdc_sft_l2.size(); i++) event.ltdc_sft_l2[i].clear();
   for(int i=0; i<event.ltdc_sft_l3.size(); i++) event.ltdc_sft_l3[i].clear();
   for(int i=0; i<event.ltdc_sft_l4.size(); i++) event.ltdc_sft_l4[i].clear();
   for(int i=0; i<event.ltdc_sft_l5.size(); i++) event.ltdc_sft_l5[i].clear();
   for(int i=0; i<event.ltdc_sft_l6.size(); i++) event.ltdc_sft_l6[i].clear();
}

void DefineHistograms( const char *filename )
{ 
  new TFile( filename, "recreate" );

  HBTree("tree","tree");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  
  //////Timing counter
  //UTOF
  tree->Branch("id_utof_l", &event.id_utof_l);
  tree->Branch("id_utof_r", &event.id_utof_r);
  tree->Branch("tot_utof_l", &event.tot_utof_l);
  tree->Branch("tot_utof_r", &event.tot_utof_r);
  tree->Branch("ltdc_utof_l", &event.ltdc_utof_l);
  tree->Branch("ltdc_utof_r", &event.ltdc_utof_r);

  //DTOF
  tree->Branch("id_dtof_l", &event.id_dtof_l);
  tree->Branch("id_dtof_r", &event.id_dtof_r);
  tree->Branch("tot_dtof_l", &event.tot_dtof_l);
  tree->Branch("tot_dtof_r", &event.tot_dtof_r);
  tree->Branch("ltdc_dtof_l", &event.ltdc_dtof_l);
  tree->Branch("ltdc_dtof_r", &event.ltdc_dtof_r);

  //LTOF
  tree->Branch("id_ltof_l", &event.id_ltof_l);
  tree->Branch("id_ltof_r", &event.id_ltof_r);
  tree->Branch("tot_ltof_l", &event.tot_ltof_l);
  tree->Branch("tot_ltof_r", &event.tot_ltof_r);
  tree->Branch("ltdc_ltof_l", &event.ltdc_ltof_l);
  tree->Branch("ltdc_ltof_r", &event.ltdc_ltof_r);

  //T0
  tree->Branch("id_t0_l", &event.id_t0_l);
  tree->Branch("id_t0_r", &event.id_t0_r);
  tree->Branch("tot_t0_l", &event.tot_t0_l);
  tree->Branch("tot_t0_r", &event.tot_t0_r);
  tree->Branch("ltdc_t0_l", &event.ltdc_t0_l);
  tree->Branch("ltdc_t0_r", &event.ltdc_t0_r);

  //T0R
  tree->Branch("id_t0r_l", &event.id_t0r_l);
  tree->Branch("id_t0r_r", &event.id_t0r_r);
  tree->Branch("tot_t0r_l", &event.tot_t0r_l);
  tree->Branch("tot_t0r_r", &event.tot_t0r_r);
  tree->Branch("ltdc_t0r_l", &event.ltdc_t0r_l);
  tree->Branch("ltdc_t0r_r", &event.ltdc_t0r_r);

  //Bref
  tree->Branch("id_bref_l", &event.id_bref_l);
  tree->Branch("id_bref_r", &event.id_bref_r);
  tree->Branch("tot_bref_l", &event.tot_bref_l);
  tree->Branch("tot_bref_r", &event.tot_bref_r);
  tree->Branch("ltdc_bref_l", &event.ltdc_bref_l);
  tree->Branch("ltdc_bref_r", &event.ltdc_bref_r);

  //T1 from K1.8BR
  tree->Branch("id_t1_l", &event.id_t1_l);
  tree->Branch("id_t1_r", &event.id_t1_r);
  tree->Branch("tot_t1_l", &event.tot_t1_l);
  tree->Branch("tot_t1_r", &event.tot_t1_r);
  tree->Branch("ltdc_t1_l", &event.ltdc_t1_l);
  tree->Branch("ltdc_t1_r", &event.ltdc_t1_r);

  //BHT from K1.8BR
  tree->Branch("id_bht_l", &event.id_bht_l);
  tree->Branch("id_bht_r", &event.id_bht_r);
  tree->Branch("tot_bht_l", &event.tot_bht_l);
  tree->Branch("tot_bht_r", &event.tot_bht_r);
  tree->Branch("ltdc_bht_l", &event.ltdc_bht_l);
  tree->Branch("ltdc_bht_r", &event.ltdc_bht_r);

  ////Drift chamber
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

  //KLDC
  tree->Branch("id_kldc_l1", &event.id_kldc_l1);
  tree->Branch("id_kldc_l2", &event.id_kldc_l2);
  tree->Branch("id_kldc_l3", &event.id_kldc_l3);
  tree->Branch("id_kldc_l4", &event.id_kldc_l4);
  tree->Branch("id_kldc_l5", &event.id_kldc_l5);
  tree->Branch("id_kldc_l6", &event.id_kldc_l6);
  tree->Branch("id_kldc_l7", &event.id_kldc_l7);
  tree->Branch("id_kldc_l8", &event.id_kldc_l8);

  tree->Branch("tot_kldc_l1", &event.tot_kldc_l1);
  tree->Branch("tot_kldc_l2", &event.tot_kldc_l2);
  tree->Branch("tot_kldc_l3", &event.tot_kldc_l3);
  tree->Branch("tot_kldc_l4", &event.tot_kldc_l4);
  tree->Branch("tot_kldc_l5", &event.tot_kldc_l5);
  tree->Branch("tot_kldc_l6", &event.tot_kldc_l6);
  tree->Branch("tot_kldc_l7", &event.tot_kldc_l7);
  tree->Branch("tot_kldc_l8", &event.tot_kldc_l8);

  tree->Branch("ltdc_kldc_l1", &event.ltdc_kldc_l1);
  tree->Branch("ltdc_kldc_l2", &event.ltdc_kldc_l2);
  tree->Branch("ltdc_kldc_l3", &event.ltdc_kldc_l3);
  tree->Branch("ltdc_kldc_l4", &event.ltdc_kldc_l4);
  tree->Branch("ltdc_kldc_l5", &event.ltdc_kldc_l5);
  tree->Branch("ltdc_kldc_l6", &event.ltdc_kldc_l6);
  tree->Branch("ltdc_kldc_l7", &event.ltdc_kldc_l7);
  tree->Branch("ltdc_kldc_l8", &event.ltdc_kldc_l8);

  ////Fiber tracker
  //BFT
  tree->Branch("id_bft_l1", &event.id_bft_l1);
  tree->Branch("id_bft_l2", &event.id_bft_l2);
  tree->Branch("id_bft_l3", &event.id_bft_l3);
  tree->Branch("id_bft_l4", &event.id_bft_l4);
  tree->Branch("id_bft_l5", &event.id_bft_l5);
  tree->Branch("id_bft_l6", &event.id_bft_l6);

  tree->Branch("tot_bft_l1", &event.tot_bft_l1);
  tree->Branch("tot_bft_l2", &event.tot_bft_l2);
  tree->Branch("tot_bft_l3", &event.tot_bft_l3);
  tree->Branch("tot_bft_l4", &event.tot_bft_l4);
  tree->Branch("tot_bft_l5", &event.tot_bft_l5);
  tree->Branch("tot_bft_l6", &event.tot_bft_l6);

  tree->Branch("ltdc_bft_l1", &event.ltdc_bft_l1);
  tree->Branch("ltdc_bft_l2", &event.ltdc_bft_l2);
  tree->Branch("ltdc_bft_l3", &event.ltdc_bft_l3);
  tree->Branch("ltdc_bft_l4", &event.ltdc_bft_l4);
  tree->Branch("ltdc_bft_l5", &event.ltdc_bft_l5);
  tree->Branch("ltdc_bft_l6", &event.ltdc_bft_l6);

  //SFT
  tree->Branch("id_sft_l1", &event.id_sft_l1);
  tree->Branch("id_sft_l2", &event.id_sft_l2);
  tree->Branch("id_sft_l3", &event.id_sft_l3);
  tree->Branch("id_sft_l4", &event.id_sft_l4);
  tree->Branch("id_sft_l5", &event.id_sft_l5);
  tree->Branch("id_sft_l6", &event.id_sft_l6);

  tree->Branch("tot_sft_l1", &event.tot_sft_l1);
  tree->Branch("tot_sft_l2", &event.tot_sft_l2);
  tree->Branch("tot_sft_l3", &event.tot_sft_l3);
  tree->Branch("tot_sft_l4", &event.tot_sft_l4);
  tree->Branch("tot_sft_l5", &event.tot_sft_l5);
  tree->Branch("tot_sft_l6", &event.tot_sft_l6);

  tree->Branch("ltdc_sft_l1", &event.ltdc_sft_l1);
  tree->Branch("ltdc_sft_l2", &event.ltdc_sft_l2);
  tree->Branch("ltdc_sft_l3", &event.ltdc_sft_l3);
  tree->Branch("ltdc_sft_l4", &event.ltdc_sft_l4);
  tree->Branch("ltdc_sft_l5", &event.ltdc_sft_l5);
  tree->Branch("ltdc_sft_l6", &event.ltdc_sft_l6);
}
