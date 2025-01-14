void plot_Filter02( )
{
  gROOT->Reset();
  TCanvas *c1 = new TCanvas("c1","c1",900,900);

  c1->Divide(1,2);

  TH1F *hist1 = new TH1F("h1", "Beam TOF1", 1000,120,135);
  TH1F *hist2 = new TH1F("h2", "Beam TOF1", 1000,120,135);
  TH1F *hist3 = new TH1F("h3", "Beam TOF1", 1000,120,135);
 
  {
     TFile *fin = new TFile("rootfile/online_ana01/run1371_EB01.root");
     TTree *tree = (TTree*)fin->Get("tree");

     std::vector<std::vector<double>> *ltdc_utof_l=0, *ltdc_utof_r=0;  
     std::vector<std::vector<double>> *ltdc_t1_l=0, *ltdc_t1_r=0;  

     TBranch *Bltdc_utof_l=0, *Bltdc_utof_r=0;
     TBranch *Bltdc_t1_l=0, *Bltdc_t1_r=0; 

     tree->SetBranchAddress("ltdc_utof_l", &ltdc_utof_l, &Bltdc_utof_l );
     tree->SetBranchAddress("ltdc_utof_r", &ltdc_utof_r, &Bltdc_utof_r );
     tree->SetBranchAddress("ltdc_t1_l", &ltdc_t1_l, &Bltdc_t1_l );
     tree->SetBranchAddress("ltdc_t1_r", &ltdc_t1_r, &Bltdc_t1_r );
  
     int n = tree->GetEntries();
     for(int a = 0; a<n; ++a){
        tree->GetEntry(a);

        int nutof_l = ltdc_utof_l->at(0).size();
        int nutof_r = ltdc_utof_r->at(0).size();
        int nt1_l = ltdc_t1_l->at(0).size();
        int nt1_r = ltdc_t1_r->at(0).size();
        
        for(int i1=0; i1<nutof_l; i1++){
           for(int i2=i1; i2<nutof_r; i2++){
              for(int i3=i2; i3<nt1_l; i3++){
                 for(int i4=i3; i4<nt1_r; i4++){

                    double t1 = ltdc_utof_l->at(0).at(i1);
                    double t2 = ltdc_utof_r->at(0).at(i2);
                    double t3 = ltdc_t1_l->at(0).at(i3);
                    double t4 = ltdc_t1_r->at(0).at(i4);
                    
                    double tof = (t3+t4)/2.-(t1+t2)/2.;
                            
                    hist1->Fill(tof);
                 }
              }
           }
        }
     }
  }
  // {
  //    TFile *fin = new TFile("rootfile/online_ana01/dcmrun1371_EB01.root");
  //    TTree *tree = (TTree*)fin->Get("tree");

  //    std::vector<std::vector<double>> *ltdc_utof_l=0, *ltdc_utof_r=0;  
  //    std::vector<std::vector<double>> *ltdc_t1_l=0, *ltdc_t1_r=0;  

  //    TBranch *Bltdc_utof_l=0, *Bltdc_utof_r=0;
  //    TBranch *Bltdc_t1_l=0, *Bltdc_t1_r=0; 

  //    tree->SetBranchAddress("ltdc_utof_l", &ltdc_utof_l, &Bltdc_utof_l );
  //    tree->SetBranchAddress("ltdc_utof_r", &ltdc_utof_r, &Bltdc_utof_r );
  //    tree->SetBranchAddress("ltdc_t1_l", &ltdc_t1_l, &Bltdc_t1_l );
  //    tree->SetBranchAddress("ltdc_t1_r", &ltdc_t1_r, &Bltdc_t1_r );
  
  //    int n = tree->GetEntries();
  //    for(int a = 0; a<n; ++a){
  //       tree->GetEntry(a);

  //       int nutof_l = ltdc_utof_l->at(0).size();
  //       int nutof_r = ltdc_utof_r->at(0).size();
  //       int nt1_l = ltdc_t1_l->at(0).size();
  //       int nt1_r = ltdc_t1_r->at(0).size();
        
  //       for(int i1=0; i1<nutof_l; i1++){
  //          for(int i2=0; i2<nutof_r; i2++){
  //             for(int i3=0; i3<nt1_l; i3++){
  //                for(int i4=0; i4<nt1_r; i4++){

  //                   double t1 = ltdc_utof_l->at(0).at(i1);
  //                   double t2 = ltdc_utof_r->at(0).at(i2);
  //                   double t3 = ltdc_t1_l->at(0).at(i3);
  //                   double t4 = ltdc_t1_r->at(0).at(i4);
                    
  //                   double tof = (t3+t4)/2.-(t1+t2)/2.;
                            
  //                   hist2->Fill(tof);
  //                }
  //             }
  //          }
  //       }
  //    }
  // }
  // {
  //    TFile *fin = new TFile("rootfile/online_ana01/run1371_EB01.root");
  //    TTree *tree = (TTree*)fin->Get("tree");

  //    std::vector<std::vector<double>> *ltdc_utof_l=0, *ltdc_utof_r=0;  
  //    std::vector<std::vector<double>> *ltdc_t1_l=0, *ltdc_t1_r=0;  

  //    TBranch *Bltdc_utof_l=0, *Bltdc_utof_r=0;
  //    TBranch *Bltdc_t1_l=0, *Bltdc_t1_r=0; 

  //    tree->SetBranchAddress("ltdc_utof_l", &ltdc_utof_l, &Bltdc_utof_l );
  //    tree->SetBranchAddress("ltdc_utof_r", &ltdc_utof_r, &Bltdc_utof_r );
  //    tree->SetBranchAddress("ltdc_t1_l", &ltdc_t1_l, &Bltdc_t1_l );
  //    tree->SetBranchAddress("ltdc_t1_r", &ltdc_t1_r, &Bltdc_t1_r );
  
  //    int n = tree->GetEntries();
  //    for(int a = 0; a<n; ++a){
  //       tree->GetEntry(a);

  //       int nutof_l = ltdc_utof_l->at(0).size();
  //       int nutof_r = ltdc_utof_r->at(0).size();
  //       int nt1_l = ltdc_t1_l->at(0).size();
  //       int nt1_r = ltdc_t1_r->at(0).size();
        
  //       for(int i1=0; i1<nutof_l; i1++){
  //          for(int i2=0; i2<nutof_r; i2++){
  //             for(int i3=0; i3<nt1_l; i3++){
  //                for(int i4=0; i4<nt1_r; i4++){

  //                   double t1 = ltdc_utof_l->at(0).at(i1);
  //                   double t2 = ltdc_utof_r->at(0).at(i2);
  //                   double t3 = ltdc_t1_l->at(0).at(i3);
  //                   double t4 = ltdc_t1_r->at(0).at(i4);
                    
  //                   double tof = (t3+t4)/2.-(t1+t2)/2.;
                            
  //                   hist3->Fill(tof+1.5);
  //                }
  //             }
  //          }
  //       }
  //    }
  // }

  c1->cd(1);
  hist2->SetLineColor(1);
  hist2->Draw();

  // hist1->SetLineColor(4);
  // hist1->Draw("same");

  // hist3->SetLineColor(2);
  // hist3->Draw("same");

  // c1->cd(2);

  // hist3->SetLineColor(2);
  // hist3->Draw();

  // hist2->SetLineColor(1);
  // hist2->Draw("same");

  // hist1->SetLineColor(4);
  // hist1->Draw("same");





}
