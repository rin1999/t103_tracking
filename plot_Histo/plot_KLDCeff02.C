#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"

void plot_KLDCeff02( const char* filename )
{
  TFile *fin = new TFile( filename );
  TTree *tree = (TTree*)fin->Get("tree");

  std::vector<std::vector<double>> *ltdc_utof_l=0, *ltdc_utof_r=0;  
  std::vector<std::vector<double>> *ltdc_ltof_l=0, *ltdc_ltof_r=0;  

  TBranch *Bltdc_utof_l=0, *Bltdc_utof_r=0; 
  TBranch *Bltdc_ltof_l=0, *Bltdc_ltof_r=0; 

  tree->SetBranchAddress("ltdc_utof_l", &ltdc_utof_l, &Bltdc_utof_l );
  tree->SetBranchAddress("ltdc_utof_r", &ltdc_utof_r, &Bltdc_utof_r );
  tree->SetBranchAddress("ltdc_ltof_l", &ltdc_ltof_l, &Bltdc_ltof_l );
  tree->SetBranchAddress("ltdc_ltof_r", &ltdc_ltof_r, &Bltdc_ltof_r );
  
  std::vector<std::vector<double>> *ltdc_kldc_l1 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l2 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l3 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l4 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l5 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l6 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l7 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l8 = 0;

  TBranch *Bltdc_kldc_l1 = 0;
  TBranch *Bltdc_kldc_l2 = 0;
  TBranch *Bltdc_kldc_l3 = 0;
  TBranch *Bltdc_kldc_l4 = 0;
  TBranch *Bltdc_kldc_l5 = 0;
  TBranch *Bltdc_kldc_l6 = 0;
  TBranch *Bltdc_kldc_l7 = 0;
  TBranch *Bltdc_kldc_l8 = 0;
  
  tree->SetBranchAddress("ltdc_kldc_l1", &ltdc_kldc_l1, &Bltdc_kldc_l1 );
  tree->SetBranchAddress("ltdc_kldc_l2", &ltdc_kldc_l2, &Bltdc_kldc_l2 );
  tree->SetBranchAddress("ltdc_kldc_l3", &ltdc_kldc_l3, &Bltdc_kldc_l3 );
  tree->SetBranchAddress("ltdc_kldc_l4", &ltdc_kldc_l4, &Bltdc_kldc_l4 );
  tree->SetBranchAddress("ltdc_kldc_l5", &ltdc_kldc_l5, &Bltdc_kldc_l5 );
  tree->SetBranchAddress("ltdc_kldc_l6", &ltdc_kldc_l6, &Bltdc_kldc_l6 );
  tree->SetBranchAddress("ltdc_kldc_l7", &ltdc_kldc_l7, &Bltdc_kldc_l7 );
  tree->SetBranchAddress("ltdc_kldc_l8", &ltdc_kldc_l8, &Bltdc_kldc_l8 );

  const int N = 8;
  const int MaxHits = 128;

  TH1F *histdc[N];
  histdc[0] = new TH1F("h11", "Nhits KLDC1 V1", 70, 0, 70);
  histdc[1] = new TH1F("h12", "Nhits KLDC1 V2", 70, 0, 70);
  histdc[2] = new TH1F("h13", "Nhits KLDC1 U1", 70, 0, 70);
  histdc[3] = new TH1F("h14", "Nhits KLDC1 U2", 70, 0, 70);
  histdc[4] = new TH1F("h15", "Nhits KLDC2 U1", 70, 0, 70);
  histdc[5] = new TH1F("h16", "Nhits KLDC2 U2", 70, 0, 70);
  histdc[6] = new TH1F("h17", "Nhits KLDC2 V1", 70, 0, 70);
  histdc[7] = new TH1F("h18", "Nhits KLDC2 V2", 70, 0, 70);

  TH1F *histdc2[N];
  histdc2[0] = new TH1F("h21", "Id KLDC1 V1", 120, 0, 120);
  histdc2[1] = new TH1F("h22", "Id KLDC1 V2", 120, 0, 120);
  histdc2[2] = new TH1F("h23", "Id KLDC1 U1", 120, 0, 120);
  histdc2[3] = new TH1F("h24", "Id KLDC1 U2", 120, 0, 120);
  histdc2[4] = new TH1F("h25", "Id KLDC2 U1", 120, 0, 120);
  histdc2[5] = new TH1F("h26", "Id KLDC2 U2", 120, 0, 120);
  histdc2[6] = new TH1F("h27", "Id KLDC2 V1", 120, 0, 120);
  histdc2[7] = new TH1F("h28", "Id KLDC2 V2", 120, 0, 120);

  TH1F *histdc3[N];
  histdc3[0] = new TH1F("h31", "Tdc 1st KLDC1 V1", 400, 800, 1300);
  histdc3[1] = new TH1F("h32", "Tdc 1st KLDC1 V2", 400, 800, 1300);
  histdc3[2] = new TH1F("h33", "Tdc 1st KLDC1 U1", 400, 800, 1300);
  histdc3[3] = new TH1F("h34", "Tdc 1st KLDC1 U2", 400, 800, 1300);
  histdc3[4] = new TH1F("h35", "Tdc 1st KLDC2 U1", 400, 800, 1300);
  histdc3[5] = new TH1F("h36", "Tdc 1st KLDC2 U2", 400, 800, 1300);
  histdc3[6] = new TH1F("h37", "Tdc 1st KLDC2 V1", 400, 800, 1300);
  histdc3[7] = new TH1F("h38", "Tdc 1st KLDC2 V2", 400, 800, 1300);

  const int dctdclow  =  900;
  const int dctdchigh = 1200;

  const int trigtdclow1  = 990;
  const int trigtdchigh1 = 1010;
  const int trigtdclow2  = 1040;
  const int trigtdchigh2 = 1060;
  
  int n = tree->GetEntries();
  for(int a = 0; a<n; ++a){
    tree->GetEntry(a);

    if(a%1000==0)
       std::cout<< n << " : " << a <<std::endl;

    bool trigflag = false;

    bool trigflag_utof = false;
    bool trigflag1 = false;
    bool trigflag2 = false;

    int size_utofl = ltdc_utof_l->at(0).size();
    int size_utofr = ltdc_utof_r->at(0).size();
    for(int j=0; j<size_utofl; ++j){
       double tdc1 = ltdc_utof_l->at(0).at(j);
       if( trigtdclow1<tdc1 && tdc1<trigtdchigh1 ) trigflag1=true;
    }
    for(int j=0; j<size_utofr; ++j){
       double tdc2 = ltdc_utof_r->at(0).at(j);
       if( trigtdclow1<tdc2 && tdc2<trigtdchigh1 ) trigflag2=true;
    }
    if( trigflag1 && trigflag2 ) trigflag_utof=true;
    
    bool trigflag_ltof = false;
    bool trigflag3[6];
    bool trigflag4[6];
    for( int i=0; i<6; i++ ){
          trigflag3[i]=false;
          trigflag4[i]=false;
    }

    for( int i=0; i<6; i++ ){
       int size_ltofl = ltdc_ltof_l->at(i).size();
       int size_ltofr = ltdc_ltof_r->at(i).size();

       for(int j=0; j<size_ltofl; ++j){
          double tdc3 = ltdc_ltof_l->at(i).at(j);
          if( trigtdclow2<tdc3 && tdc3<trigtdchigh2 ) trigflag3[i]=true;
       }
       for(int j=0; j<size_ltofr; ++j){
          double tdc4 = ltdc_ltof_r->at(i).at(j);
          if( trigtdclow2<tdc4 && tdc4<trigtdchigh2 ) trigflag4[i]=true;
       }
    }

    if( (trigflag3[0] && trigflag4[0]) || (trigflag3[1] && trigflag4[1]) ||
        (trigflag3[2] && trigflag4[2]) || (trigflag3[3] && trigflag4[3]) ||
        (trigflag3[4] && trigflag4[4]) || (trigflag3[5] && trigflag4[5]) ) trigflag_ltof=true;
    
    if(trigflag_utof && trigflag_ltof ) trigflag=true;

    //DC analysis
    int nhitsdc[N];
    bool flag_kldc_l1[MaxHits];
    bool flag_kldc_l2[MaxHits];
    bool flag_kldc_l3[MaxHits];
    bool flag_kldc_l4[MaxHits];
    bool flag_kldc_l5[MaxHits];
    bool flag_kldc_l6[MaxHits];
    bool flag_kldc_l7[MaxHits];
    bool flag_kldc_l8[MaxHits];

    double ltdc1st_kldc_l1[MaxHits];
    double ltdc1st_kldc_l2[MaxHits];
    double ltdc1st_kldc_l3[MaxHits];
    double ltdc1st_kldc_l4[MaxHits];
    double ltdc1st_kldc_l5[MaxHits];
    double ltdc1st_kldc_l6[MaxHits];
    double ltdc1st_kldc_l7[MaxHits];
    double ltdc1st_kldc_l8[MaxHits];
    
    if(trigflag){  
       for( int i=0; i<N; i++ ) nhitsdc[i]=0;
       for( int i=0; i<MaxHits; i++ ){
          flag_kldc_l1[i]=false;
          flag_kldc_l2[i]=false;
          flag_kldc_l3[i]=false;
          flag_kldc_l4[i]=false;
          flag_kldc_l5[i]=false;
          flag_kldc_l6[i]=false;
          flag_kldc_l7[i]=false;
          flag_kldc_l8[i]=false;

          ltdc1st_kldc_l1[i]=9999.;
          ltdc1st_kldc_l2[i]=9999.;
          ltdc1st_kldc_l3[i]=9999.;
          ltdc1st_kldc_l4[i]=9999.;
          ltdc1st_kldc_l5[i]=9999.;
          ltdc1st_kldc_l6[i]=9999.;
          ltdc1st_kldc_l7[i]=9999.;
          ltdc1st_kldc_l8[i]=9999.;
       }

       //cout<< "*****" << endl;
       for( int i=0; i<MaxHits; i++ ){
          int size1 = ltdc_kldc_l1->at(i).size();
          for(int j=0; j<size1; ++j){
             double ltdc = ltdc_kldc_l1->at(i).at(j);
             if( dctdclow<ltdc && ltdc<dctdchigh ){
                flag_kldc_l1[i]=true;
                if( ltdc<ltdc1st_kldc_l1[i] ) ltdc1st_kldc_l1[i]=ltdc;
             }
          }
          int size2 = ltdc_kldc_l2->at(i).size();
          for(int j=0; j<size2; ++j){
             double ltdc = ltdc_kldc_l2->at(i).at(j);
             if( dctdclow<ltdc && ltdc<dctdchigh ){
                flag_kldc_l2[i]=true;
                if( ltdc<ltdc1st_kldc_l2[i] ) ltdc1st_kldc_l2[i]=ltdc;
             }
          }
          int size3 = ltdc_kldc_l3->at(i).size();
          for(int j=0; j<size3; ++j){
             double ltdc = ltdc_kldc_l3->at(i).at(j);
             if( dctdclow<ltdc && ltdc<dctdchigh ){
                flag_kldc_l3[i]=true;
                if( ltdc<ltdc1st_kldc_l3[i] ) ltdc1st_kldc_l3[i]=ltdc;
             }
          }
          int size4 = ltdc_kldc_l4->at(i).size();
          for(int j=0; j<size4; ++j){
             double ltdc = ltdc_kldc_l4->at(i).at(j);
             if( dctdclow<ltdc && ltdc<dctdchigh ){
                flag_kldc_l4[i]=true;
                if( ltdc<ltdc1st_kldc_l4[i] ) ltdc1st_kldc_l4[i]=ltdc;
             }
          }
          int size5 = ltdc_kldc_l5->at(i).size();
          for(int j=0; j<size5; ++j){
             double ltdc = ltdc_kldc_l5->at(i).at(j);
             if( dctdclow<ltdc && ltdc<dctdchigh ){
                flag_kldc_l5[i]=true;
                if( ltdc<ltdc1st_kldc_l5[i] ) ltdc1st_kldc_l5[i]=ltdc;
             }
          }
          int size6 = ltdc_kldc_l6->at(i).size();
          for(int j=0; j<size6; ++j){
             double ltdc = ltdc_kldc_l6->at(i).at(j);
             if( dctdclow<ltdc && ltdc<dctdchigh ){
                flag_kldc_l6[i]=true;
                if( ltdc<ltdc1st_kldc_l6[i] ) ltdc1st_kldc_l6[i]=ltdc;
             }
          }
          int size7 = ltdc_kldc_l7->at(i).size();
          for(int j=0; j<size7; ++j){
             double ltdc = ltdc_kldc_l7->at(i).at(j);
             if( dctdclow<ltdc && ltdc<dctdchigh ){
                flag_kldc_l7[i]=true;
                if( ltdc<ltdc1st_kldc_l7[i] ) ltdc1st_kldc_l7[i]=ltdc;
             }
          }
          int size8 = ltdc_kldc_l8->at(i).size();
          for(int j=0; j<size8; ++j){
             double ltdc = ltdc_kldc_l8->at(i).at(j);
             if( dctdclow<ltdc && ltdc<dctdchigh ){
                flag_kldc_l8[i]=true;
                if( ltdc<ltdc1st_kldc_l8[i] ) ltdc1st_kldc_l8[i]=ltdc;
             }
          }
       }
       
       for( int i=0; i<MaxHits; i++ ){
          if( flag_kldc_l1[i] ){
             nhitsdc[0]++;
             histdc2[0]->Fill(i);
             histdc3[0]->Fill(ltdc1st_kldc_l1[i]);
          }
          if( flag_kldc_l2[i] ){
             nhitsdc[1]++;
             histdc2[1]->Fill(i);
             histdc3[1]->Fill(ltdc1st_kldc_l2[i]);
          }
          if( flag_kldc_l3[i] ){
             nhitsdc[2]++;
             histdc2[2]->Fill(i);
             histdc3[2]->Fill(ltdc1st_kldc_l3[i]);
          }
          if( flag_kldc_l4[i] ){
             nhitsdc[3]++;
             histdc2[3]->Fill(i);
             histdc3[3]->Fill(ltdc1st_kldc_l4[i]);
          }
          if( flag_kldc_l5[i] ){
             nhitsdc[4]++;
             histdc2[4]->Fill(i);
             histdc3[4]->Fill(ltdc1st_kldc_l5[i]);
          }
          if( flag_kldc_l6[i] ){
             nhitsdc[5]++;
             histdc2[5]->Fill(i);
             histdc3[5]->Fill(ltdc1st_kldc_l6[i]);
          }
          if( flag_kldc_l7[i] ){
             nhitsdc[6]++;
             histdc2[6]->Fill(i);
             histdc3[6]->Fill(ltdc1st_kldc_l7[i]);
          }
          if( flag_kldc_l8[i] ){
             nhitsdc[7]++;
             histdc2[7]->Fill(i);
             histdc3[7]->Fill(ltdc1st_kldc_l8[i]);
          }
       }
       for( int i=0; i<N; i++ ) histdc[i]->Fill(nhitsdc[i]);
    }
  }
  
  TCanvas *c1 = new TCanvas("c1","c1", 1200,600);
  c1->Divide(4,2);
  for( int i=0; i<N; i++ ){
     c1->cd(i+1);
     c1->cd(i+1)->SetLogy();
     histdc[i]->Draw();
  }
  
  TCanvas *c2 = new TCanvas("c2","c2", 1200,600);
  c2->Divide(4,2);
  for( int i=0; i<N; i++ ){
     c2->cd(i+1);
     histdc2[i]->Draw();
  }

  TCanvas *c3 = new TCanvas("c3","c3", 1200,600);
  c3->Divide(4,2);
  for( int i=0; i<N; i++ ){
     c3->cd(i+1);
     histdc3[i]->Draw();
  }

  double Tot_Num[N],NoHit[N],eff[N],effe[N];
  for( int i=0; i<N; i++ ){
     Tot_Num[i]=0.; NoHit[i]=0.;
     eff[i]=0.; effe[i]=0.;
  }

  for( int i=0; i<N; i++ ){
     Tot_Num[i]=double(histdc[i]->Integral());
     NoHit[i]=double(histdc[i]->Integral(1,1));
     eff[i] = 1.0-NoHit[i]/Tot_Num[i];
     effe[i] = pow(((NoHit[i]/Tot_Num[i])*(1.0-NoHit[i]/Tot_Num[i]))/Tot_Num[i], 0.5);
  }

  cout<< "***********************" << endl;
  cout<< "Eff.Lay KLDC1 X1 : "<< eff[0]*100.0 << " +- " << effe[0]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay KLDC1 X2 : "<< eff[1]*100.0 << " +- " << effe[1]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay KLDC1 U1 : "<< eff[2]*100.0 << " +- " << effe[2]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay KLDC1 V1 : "<< eff[3]*100.0 << " +- " << effe[3]*100.0 << " [%]" << endl;
  cout<< "***********************" << endl;
  cout<< "***********************" << endl;
  cout<< "Eff.Lay KLDC2 X1 : "<< eff[4]*100.0 << " +- " << effe[4]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay KLDC2 X2 : "<< eff[5]*100.0 << " +- " << effe[5]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay KLDC2 U1 : "<< eff[6]*100.0 << " +- " << effe[6]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay KLDC2 V1 : "<< eff[7]*100.0 << " +- " << effe[7]*100.0 << " [%]" << endl;
  cout<< "***********************" << endl;





  
}
