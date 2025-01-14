void calcInitXTparDC(const char *filename)
{

  TFile *fin = new TFile(filename);
  fin->cd();
  fin->ReadAll();
  
  TCanvas *c1;
  if (!(c1=(TCanvas *)gROOT->FindObject("c1")))
    c1 = new TCanvas("c1","c1",700,700);

  c1->SetFillColor(10);
  c1->SetGrid();

  FILE *fout;

  const int DT = 200;
  //const int DT = 100;
  
  for (int layer=1; layer<=8; layer++) {  
    char name[100];
    //sprintf(name,"h%d",layer*100+12);
    sprintf(name,"h%d",layer*100+3);
    TH1F *h1 = (TH1F *)fin->FindObject(name);
    //int bin1=h1->FindBin(-50);
    int bin1=h1->FindBin(-1);
    int bin2=h1->FindBin(DT);

//     /* for background study*/
//     int bin3=h1->FindBin(-50);
//     int bin4=h1->FindBin(-20);
//     double BG_bin =h1->Integral(bin3, bin4)/double(bin4-bin3+1);
//     double Entry = h1->Integral(bin1, bin2) - BG_bin*(double)(bin2-bin1+1);

    double Entry = h1->Integral(bin1, bin2);

    double MaxLength=10.41;
    //double MaxLength=4.5;
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
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(2);
    gr->SetMarkerSize(0.9);
    gr->Draw("ap");
    
    TF1 *func = new TF1("func", "pol5(0)", -5., DT );
    
    func->SetParameter(0, 0);
    func->SetParameter(1, 1.0267e-1);
    func->SetParameter(2,-1.87541e-3);
    func->SetParameter(3, 2.0871e-5);
    func->SetParameter(4,-1.41421e-7);
    func->SetParameter(5, 4.25782e-10);
    
    gr->Fit("func","","same",-1., DT );
    
    func->SetRange(-5, DT);
    func->SetLineColor(kBlue);
    func->Draw("same");

    sprintf(name,"h%d",layer*100+12);
    //sprintf(name,"h%d",layer*100+3);
    TH2F *h2 = (TH2F *)fin->FindObject(name);
    h2->RebinX(8.);
    h2->Draw("same");

    c1->Update();
    
    getchar();

    char foutname[100];
    sprintf(foutname, "DC_calc/xtInitpar_dc%d.d", layer+100);
    //sprintf(foutname, "DC_calc/xtInitpar_dc%d.d", layer+200);
    fout = fopen(foutname, "w");

    fprintf(fout, "%d %e %e %e %e %e\n", layer,
	    func->GetParameter(1), func->GetParameter(2),
	    func->GetParameter(3), func->GetParameter(4),
	    func->GetParameter(5));
    
    fclose(fout);
    
    func->Delete();
    gr->Delete();

    c1->Clear();
  }
}
