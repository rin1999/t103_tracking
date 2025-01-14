void plot_Filter01( int run, double factor )
{
  TH1F *hist1 = new TH1F("h1", "Beam TOF1", 100,-6,4);
  TH1F *hist2 = new TH1F("h2", "Beam TOF1", 100,-6,4);
  TH1F *hist3 = new TH1F("h3", "Beam TOF1", 100,-6,4);

  char file1[100], file2[100];
  sprintf(file1, "run%d_Hodo00.root", run);
  sprintf(file2, "dcmrun%d_Hodo00.root", run);

  {
     TFile *fin = new TFile( file1 );
     TTree *tree = (TTree*)fin->Get("tree");

     double beam_tof[32];
     int nhits_scat_tof;

     tree->SetBranchAddress("beam_tof", &beam_tof );
     tree->SetBranchAddress("nhits_scat_tof", &nhits_scat_tof );

     int n = tree->GetEntries();
     for(int a = 0; a<n; ++a){
        tree->GetEntry(a);

        if(nhits_scat_tof>0) hist1->Fill(beam_tof[0]+1.5);
     }
  }

  {
     TFile *fin = new TFile( file2 );
     TTree *tree = (TTree*)fin->Get("tree");

     double beam_tof[32];
     int nhits_scat_tof;

     tree->SetBranchAddress("beam_tof", &beam_tof );
     tree->SetBranchAddress("nhits_scat_tof", &nhits_scat_tof );

     int n = tree->GetEntries();
     for(int a = 0; a<n; ++a){
        tree->GetEntry(a);

        if(nhits_scat_tof>0) hist2->Fill(beam_tof[0]+1.5);
     }
  }

  {
     TFile *fin = new TFile( file2 );
     TTree *tree = (TTree*)fin->Get("tree");

     double beam_tof[32];
     int nhits_scat_tof;

     tree->SetBranchAddress("beam_tof", &beam_tof );
     tree->SetBranchAddress("nhits_scat_tof", &nhits_scat_tof );
     
     int n = tree->GetEntries();
     for(int a = 0; a<n; ++a){
        tree->GetEntry(a);

        if(nhits_scat_tof) hist3->Fill(beam_tof[0]+1.5, factor);
     }
  }

  TCanvas *c1 = new TCanvas("c1","c1",900,600);

  c1->Divide(1,1);

  c1->cd(1);
  hist1->SetLineColor(1);
  hist1->Draw();

  hist2->SetLineColor(4);
  hist2->Draw("same");

  hist3->SetLineColor(2);
  hist3->Draw("same");

  // c1->cd(2);
  // hist2->SetLineColor(1);
  // hist2->Draw();

  // hist1->SetLineColor(4);
  // hist1->Draw("same");

  // hist3->SetLineColor(2);
  // hist3->Draw("same");


  
  // hist3->SetLineColor(2);
  // hist3->Draw();

  // hist2->SetLineColor(1);
  // hist2->Draw("same");

  // hist1->SetLineColor(4);
  // hist1->Draw("same");





}
