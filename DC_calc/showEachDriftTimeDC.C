#include <stdio.h>

const int NwireDC1  = 128;
const int NwireDC2  = 128;

void prompt()
{
  printf("Press return");
  getchar();
}

int getMaxWire(int layer) 
{
  switch (layer) {
  case 1:
    return NwireDC1;
    break;
  case 2:
    return NwireDC2;
    break;
  case 3:
    return NwireDC1;
    break;
  default:
    fprintf(stderr, "No such layer %d in DC\n", layer);
    exit(-1);
  }
}

void showEachDriftTimeDC(int layer, const char *filename)
{
  TFile *fin = new TFile(filename);
  fin->cd();
  fin->ReadAll();

  char name[100], buf[100];

  TCanvas *c1 = new TCanvas("c1","c1",700,700);
  c1->SetFillColor(10);
  c1->SetGrid();
  
  char fout_name[100];
  sprintf(fout_name, "DC_calc/dc_layer%d_offset_1.d",layer);
  FILE *fp = fopen(fout_name, "w");
  if (!fp) {
    fprintf(stderr, "cannot open %s\n", fout_name);
    exit(1);
  }

  for (int i=1; i<=getMaxWire(layer); i++) {
    sprintf(name, "h%d", 10000*layer+1000+i);
    TH1F *h1 = (TH1F *)fin->FindObject(name);

    int binx1, binx2;
    binx1 = h1->FindBin(-30);
    binx2 = h1->FindBin(50);
    h1->GetXaxis()->SetRange(binx1, binx2);

    h1->Draw();
    c1->Update();
    printf("Ok? 0 : ");
    double data;
    gets(buf);
    if (buf[0] == ' ' || buf[0] == '\n') {
      data = 0.0;
    }
    else 
      sscanf(buf, "%lf", &data);
    fprintf(fp, "%f\n", data);
  }
}
