#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"

void plot_KLDCInfo02( const char* filename )
{
   TFile *fin = new TFile( filename );
   TTree *tree = (TTree*)fin->Get("tree");

   std::vector<double> *x0 = 0;
   std::vector<double> *y0 = 0;
   std::vector<double> *u0 = 0;
   std::vector<double> *v0 = 0;
   std::vector<double> *x1 = 0;
   std::vector<double> *y1 = 0;
   std::vector<double> *u1 = 0;
   std::vector<double> *v1 = 0;

   TBranch *Bx0 = 0;
   TBranch *By0 = 0;
   TBranch *Bu0 = 0;
   TBranch *Bv0 = 0;
   TBranch *Bx1 = 0;
   TBranch *By1 = 0;
   TBranch *Bu1 = 0;
   TBranch *Bv1 = 0;

   tree->SetBranchAddress("x0", &x0, &Bx0 );
   tree->SetBranchAddress("y0", &y0, &By0 );
   tree->SetBranchAddress("u0", &u0, &Bu0 );
   tree->SetBranchAddress("v0", &v0, &Bv0 );
   tree->SetBranchAddress("x1", &x1, &Bx1 );
   tree->SetBranchAddress("y1", &y1, &By1 );
   tree->SetBranchAddress("u1", &u1, &Bu1 );
   tree->SetBranchAddress("v1", &v1, &Bv1 );

   int nhits_scat_tof;
   int nhits_scat_tof2;

   tree->SetBranchAddress("nhits_scat_tof", &nhits_scat_tof );
   tree->SetBranchAddress("nhits_scat_tof2", &nhits_scat_tof2 );
   
   const int N = 8;

   TH2F *hist1[N];
   hist1[0] = new TH2F("h1", "KLDC Lyaer 1", 1000,-450,450,1000,-400,400);
   hist1[1] = new TH2F("h2", "KLDC Lyaer 2", 1000,-450,450,1000,-400,400);
   hist1[2] = new TH2F("h3", "KLDC Lyaer 3", 1000,-450,450,1000,-400,400);
   hist1[3] = new TH2F("h4", "KLDC Lyaer 4", 1000,-450,450,1000,-400,400);
   hist1[4] = new TH2F("h5", "KLDC Lyaer 5", 1000,-450,450,1000,-400,400);
   hist1[5] = new TH2F("h6", "KLDC Lyaer 6", 1000,-450,450,1000,-400,400);
   hist1[6] = new TH2F("h7", "KLDC Lyaer 7", 1000,-450,450,1000,-400,400);
   hist1[7] = new TH2F("h8", "KLDC Lyaer 8", 1000,-450,450,1000,-400,400);

   TH2F *hist2[N];
   hist2[0] = new TH2F("h11", "KLDC Lyaer 1", 1000,-450,450,1000,-400,400);
   hist2[1] = new TH2F("h12", "KLDC Lyaer 2", 1000,-450,450,1000,-400,400);
   hist2[2] = new TH2F("h13", "KLDC Lyaer 3", 1000,-450,450,1000,-400,400);
   hist2[3] = new TH2F("h14", "KLDC Lyaer 4", 1000,-450,450,1000,-400,400);
   hist2[4] = new TH2F("h15", "KLDC Lyaer 5", 1000,-450,450,1000,-400,400);
   hist2[5] = new TH2F("h16", "KLDC Lyaer 6", 1000,-450,450,1000,-400,400);
   hist2[6] = new TH2F("h17", "KLDC Lyaer 7", 1000,-450,450,1000,-400,400);
   hist2[7] = new TH2F("h18", "KLDC Lyaer 8", 1000,-450,450,1000,-400,400);

   double zpos[N] = {499.0, 506.8, 523.2, 531.0, 619.0, 626.8, 643.2, 651.0};

   TH1F *hist21 = new TH1F("h21", "LTOF Pos Y", 300,-400,400);
   TH1F *hist22 = new TH1F("h22", "LTOF Pos Y", 300,-400,400);

   TH2F *hist23 = new TH2F("h23", "LTOF Pos X vs Y", 300,-600,600,300,-400,400);
   TH2F *hist24 = new TH2F("h24", "LTOF Pos X vs Y", 300,-600,600,300,-400,400);
   
   int n = tree->GetEntries();
   for(int a = 0; a<100000; ++a){
      tree->GetEntry(a);
      
      bool f_scat = false;
      bool f_scat2 = false; 
      if( nhits_scat_tof>0 ) f_scat = true;
      if( nhits_scat_tof2>0 ) f_scat2 = true;
      
      // for(int i=0; i<N; ++i){
      //    int nh = x0->size();
      //    for(int j=0; j<nh; ++j){
      //       double X = x0->at(j)+u0->at(j)*zpos[i];
      //       double Y = y0->at(j)+v0->at(j)*zpos[i];
            
      //        if(f_scat) hist1[i]->Fill(X,Y,1.0);
      //        if(!f_scat2) hist2[i]->Fill(X,Y,1.0);
      //    }
      // }
      
      int nh2 = y1->size();
      for(int j=0; j<nh2; ++j){
         double X = x1->at(j);
         double Y = y1->at(j);
         if(nhits_scat_tof>0){
            hist21->Fill(Y);
            hist23->Fill(X,Y,1.0);
         }
         //if(nhits_scat_tof2==0){
         if(nhits_scat_tof>0&&nhits_scat_tof2==0){
             hist22->Fill(Y);
             hist24->Fill(X,Y,1.0);
         }
      }
   }
   
   // TCanvas *c1 = new TCanvas("c1","c1", 1600,900);
   // c1->Divide(4,2);
   // for( int i=0; i<N; ++i){
   //    c1->cd(i+1);
      
   //    hist1[i]->Draw("col");
   // }
   
   // TCanvas *c2 = new TCanvas("c2","c2", 1600,900);
   // c2->Divide(4,2);
   // for( int i=0; i<N; ++i){
   //    c2->cd(i+1);
      
   //    hist2[i]->Draw("col");
   // }

   TCanvas *c3 = new TCanvas("c3","c3", 1200,800);
   c3->Divide(2,2);
   c3->cd(1); hist21->Draw();
   c3->cd(2); hist22->Draw();
   c3->cd(3); hist23->Draw("col");
   c3->cd(4); hist24->Draw("col");

}
