#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"

void plot_KLDCInfo01( const char* filename )
{
   TFile *fin = new TFile( filename );
   TTree *tree = (TTree*)fin->Get("tree");

   std::vector<double> *x0 = 0;
   std::vector<double> *y0 = 0;
   std::vector<double> *u0 = 0;
   std::vector<double> *v0 = 0;

   TBranch *Bx0 = 0;
   TBranch *By0 = 0;
   TBranch *Bu0 = 0;
   TBranch *Bv0 = 0;

   tree->SetBranchAddress("x0", &x0, &Bx0 );
   tree->SetBranchAddress("y0", &y0, &By0 );
   tree->SetBranchAddress("u0", &u0, &Bu0 );
   tree->SetBranchAddress("v0", &v0, &Bv0 );

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

   double zpos[N] = {499.0, 506.8, 523.2, 531.0, 619.0, 626.8, 643.2, 651.0};
   
   int n = tree->GetEntries();
   for(int a = 0; a<3000000; ++a){
      tree->GetEntry(a);

       for(int i=0; i<N; ++i){
          int nh = x0->size();
          for(int j=0; j<nh; ++j){
             double X = x0->at(j)+u0->at(j)*zpos[i];
             double Y = y0->at(j)+v0->at(j)*zpos[i];

             hist1[i]->Fill(X,Y,1.0);
          }
       }
   }

   TCanvas *c1 = new TCanvas("c1","c1", 1600,900);
   c1->Divide(4,2);
   for( int i=0; i<N; ++i){
      c1->cd(i+1);
      
      hist1[i]->Draw("col");
   }
}
