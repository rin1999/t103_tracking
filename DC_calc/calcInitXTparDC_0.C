void calcInitXTparDC_0(const char *filename)
{
  TFile *fin = new TFile(filename);
  fin->cd();
  fin->ReadAll();
  
  TCanvas *c1 = new TCanvas("c1","c1",900,600);

  gROOT->SetStyle("Plain");
  gROOT->SetStyle("Plain");

  for (int layer=2; layer<=2; layer++) {  
    char name[100];
    sprintf(name,"h%d",layer*100+12);
    //sprintf(name,"h%d",layer*100+3);
    TH1F *h1 = (TH1F *)fin->FindObject(name);
    //int bin1=h1->FindBin(-50);
    int bin1=h1->FindBin(-1);
    int bin2=h1->FindBin(300);

    double Entry = h1->Integral(bin1, bin2);

    double MaxLength=10.0;
    double x[1000], y[1000];
    int nIndex=0;
    for (int i=bin1; i<=bin2; i++) {
      double Integral = h1->Integral(bin1, i);
      x[nIndex] = h1->GetBinCenter(i);
      y[nIndex] = MaxLength*Integral/Entry;
      nIndex++;
    }

    TGraph *gr = new TGraph(nIndex, x, y);
    char gr_name[100];
    gr->SetName("gr");
    gr->SetMarkerColor(2);
    gr->SetMarkerStyle(2);
    gr->SetMarkerSize(1.0);
    //gr->Draw("ap");
    
    TF1 *func = new TF1("func", "pol5(0)", -5., 300.);
    
    func->SetParameter(0, 0);
    func->SetParameter(1, 1.0267e-1);
    func->SetParameter(2,-1.87541e-3);
    func->SetParameter(3, 2.0871e-5);
    func->SetParameter(4,-1.41421e-7);
    func->SetParameter(5, 4.25782e-10);
    
    gr->Fit("func","","same", -5., 250.);

    //c1->Clear();

    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.03);
    c1->SetBottomMargin(0.15);
    TH1F* Hc1=c1->DrawFrame( -5, 0, 250, 10.0 );
    Hc1->GetXaxis()->SetLabelSize(0.060);
    Hc1->GetYaxis()->SetLabelSize(0.060);
    Hc1->GetYaxis()->SetTitle("Drift length [mm]");
    Hc1->GetYaxis()->SetTitleSize(0.06);
    Hc1->GetXaxis()->SetTitle("Drift time [ns]");
    Hc1->GetXaxis()->SetTitleSize(0.06);
    c1->SetGridx();
    c1->SetGridy();

    sprintf(name,"h%d",layer*100+19);
    TH2F *h2 = (TH2F *)fin->FindObject(name);
    h2->RebinX(1.);
    h2->Draw("same");

    func->SetRange( -5, 300 );
    func->SetLineColor(2);
    func->Draw("same");

    gROOT->SetStyle("Plain");
    gROOT->SetStyle("Plain");
    //getchar();

//     char foutname[100];
//     sprintf(foutname, "DC_calc/xtInitpar_bc%d.d", layer+100+12);
//     fout = fopen(foutname, "w");

//     fprintf(fout, "%d %e %e %e %e %e\n", layer,
// 	    func->GetParameter(1), func->GetParameter(2),
// 	    func->GetParameter(3), func->GetParameter(4),
// 	    func->GetParameter(5));
    
//     fclose(fout);
//     func->Delete();
//     gr->Delete();
//     c1->Clear();
  }
  gROOT->SetStyle("Plain");
  gROOT->SetStyle("Plain");

}
