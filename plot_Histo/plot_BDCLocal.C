void plot_BDCLocal( int num )
{
  gROOT->Reset();
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);

  const int n = 8;
  TCanvas* c11[n];

  c1->Divide(4,2);

  char histo[255];

  //TFile *f= new TFile("test.root");
  TFile *f= new TFile("extrun1745_BDC00.root");
  TH1F *h[n];

  for( int i=0; i<n; ++i){
    c11[i] = (TCanvas*) c1->GetPad(i+1);
    c1->cd(i+1);
    c11[i]->SetLogy(0); 

    sprintf(histo,"h%d",100*(i+1)+num);
    h[i] = (TH1F*) f->Get( histo );

    h[i]->Draw();
  }

}
