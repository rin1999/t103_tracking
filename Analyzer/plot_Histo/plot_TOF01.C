void plot_TOF01( int run )
{
  TH1F *hist1 = new TH1F("h1", "Beam TOF1", 200,-5,20);
  TH1F *hist2 = new TH1F("h2", "Beam TOF1", 200,-5,20);
  TH1F *hist3 = new TH1F("h3", "Beam TOF1", 200,-5,20);

  char file1[100], file2[100];
  sprintf(file1, "run00%d_Hodo00.root", run);
  //sprintf(file2, "dcmrun%d_Hodo00.root", run);

  {
     TFile *fin = new TFile( "extrun1745_Hodo00.root" );
     TTree *tree = (TTree*)fin->Get("tree");

     double scat_tof3[32];
     double de_tof3[32];
     int nhits_scat_tof3;

     tree->SetBranchAddress("scat_tof3", &scat_tof3 );
     tree->SetBranchAddress("de_tof3", &de_tof3 );
     tree->SetBranchAddress("nhits_scat_tof3", &nhits_scat_tof3 );

     int n = tree->GetEntries();
     for(int a = 0; a<n; ++a){
        tree->GetEntry(a);

        for(int i=0; i<1; ++i){
           double tof = scat_tof3[i]+2.9*de_tof3[i]-3.2;
           hist1->Fill(tof);
        }
     }
  }
  
  {
     TFile *fin = new TFile( file1 );
     TTree *tree = (TTree*)fin->Get("tree");

     double scat_tof3[32];
     double de_tof3[32];
     int nhits_scat_tof3;

     tree->SetBranchAddress("scat_tof3", &scat_tof3 );
     tree->SetBranchAddress("de_tof3", &de_tof3 );
     tree->SetBranchAddress("nhits_scat_tof3", &nhits_scat_tof3 );

     int n = tree->GetEntries();
     for(int a = 0; a<n; ++a){
        tree->GetEntry(a);

        for(int i=0; i<1; ++i){
           double tof = scat_tof3[i]+2.9*de_tof3[i]-1.8;
           hist2->Fill(tof);
        }
     }
  }
  TCanvas *c1 = new TCanvas("c1","c1",900,800);

  c1->Divide(1,2);

  c1->cd(1);
  hist1->SetLineColor(1);
  hist1->Draw();

  c1->cd(2);
  hist2->SetLineColor(1);
  hist2->Draw();

  // hist2->SetLineColor(4);
  // hist2->Draw("same");

  // hist3->SetLineColor(2);
  // hist3->Draw("same");





}
