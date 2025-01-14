void res1_2(){

//  TCanvas *c2 = new TCanvas("c2","",50,400,300,300);
//  c2->Divide(1,1);
//  TCanvas *c3 = new TCanvas("c3","",200,400,300,300);
//  c3->Divide(1,1);
  TCanvas *c5 = new TCanvas("c5","",2,1,1200,300);
  c5->Divide(4,1);
  TCanvas *c6 = new TCanvas("c6","",2,350,1200,300);
  c6->Divide(4,1);
  TCanvas *c7 = new TCanvas("c7","",4,450,300,300);
  c7->Divide(1,1);

  TTree *tr = (TTree*)_file0->Get("tree");

  unsigned short tdc[16];
  unsigned short adc[16];

  tr->SetBranchAddress("tdc",tdc);
  tr->SetBranchAddress("adc",adc);

  long N = tr->GetEntries();
//  long N = 1000;

//*/*/Change Here :)/*/*/*/*/*/*/*//
  int start = -75*12; 
//  int start = -75*12; 
//  int start = -75*12+55; 
//  int start = -75*20-35; 
  double fitMin = 80;
  double fitMax = 1000;

  int start2 = -75*12+50; 
//  int start2 = -75*10; 
//  int start2 = -75*9+35; 
//  int start2 = -75*19-80; 
  double fitMin2 = 80;
  double fitMax2 =1000;

  double tdc0M = 2000;
  double tdc1M = 1900;
//  double tdc0M = 2700;
//  double tdc1M = 2700;
  double sc0m = 1700;
  double sc1m = 1700;
  double sc3m = 1700;
  double sc0M = 1980;
  double sc1M = 1980;
  double sc3M = 1980;
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*//

  double STA, STO;
//  double d = 75.5;
  double d = 74.276;
  double STA2, STO2;
//  double d2 = 75.5;
  double d2 = 74.276;
  int J = 30;

  double p1[7],p2[7],p3[7];
  double P1,P2,P3;
  double res1[3],res2[3],res3[3],RES1,RES2,RES3;
  double rese1,rese2,rese3,RESE1,RESE2,RESE3;
  char r1[100],r2[100],r3[100];

  double q1[7];
  double Q1;
  double re[3],RE,re2[3],RE2;
  double rse,RSE,rse2,RSE2;
  char R[100],R2[100];

  TH2F * h1 = new TH2F ("h1","RPC-RF",2000,0,2000,500,-1800,-1300);

//  TH2F * h2 = new TH2F ("h2","Before Correction 1", 300,0,2000,75,0,75);
  TH2F * h2 = new TH2F ("h2","Before Correction 1", 450,50,500,75,0,75);
//  TH2F * H2 = new TH2F ("H2","After Correction 1",950,0,2000,75,-37,38);
  TH2F * H2 = new TH2F ("H2","After Correction 1",450,0,500,75,-37,38);

  TH1F * hh2 = new TH1F ("hh2","After Correction 1",75,-37,38);

  TH2F * f1 = new TH2F ("f1","RPC-RF",2000,0,2000,500,-1800,-1300);

//  TH2F * f2 = new TH2F ("f2","Before Correction 1", 300,0,2000,75,0,75);
  TH2F * f2 = new TH2F ("f2","Before Correction 1", 450,50,500,75,0,75);
//  TH2F * F2 = new TH2F ("F2","After Correction 1",950,50,2000,75,-37,38);
  TH2F * F2 = new TH2F ("F2","After Correction 1",450,50,500,75,-37,38);

  //TH2F * f2 = new TH2F ("f2","Before Correction 1", 200,50,2000,75,0,75);
  //TH2F * F2 = new TH2F ("F2","After Correction 1",950,50,2000,75,-37,38);

  TH1F * ff2 = new TH1F ("ff2","After Correction 1",75,-37,38);

//  TH1F * g = new TH1F ("g","After Correction 1",75,40,115);
//  TH1F * g = new TH1F ("g","After Correction 1",4000,-2000,2000);
  TH1F * g = new TH1F ("g","After Correction 1",75,-37,38);
  TH2F * g2 = new TH2F ("F2","After Correction 1",1950,50,2000,75,-37,38);

  for(long I=0; I<N; I++){

	tr->GetEntry(I);

	if(tdc[0]<tdc0M){

    h1->Fill(adc[0],tdc[6]-tdc[0]);

	for(int j=0; j<J; j++){

		STA = start + d*j;
	    STO = start + d*(j+1);

if (tdc[6]-tdc[0]>STA && tdc[6]-tdc[0]<STO){

	h2->Fill(adc[0],tdc[6]-tdc[0]-d*j-start);

}}


	}


	if(tdc[1]<tdc1M ){

    f1->Fill(adc[1],tdc[6]-tdc[1]);

	for(int j=0; j<J; j++){

		STA2 = start2 + d2*j;
	    STO2 = start2 + d2*(j+1);

if (tdc[6]-tdc[1]>STA2 && tdc[6]-tdc[1]<STO2){

	f2->Fill(adc[1],tdc[6]-tdc[1]-d2*j-start2);

}}


	}
	}



 c5->cd(1);
	h2->Draw();

 c6->cd(1);
	f2->Draw();
  
  h1->GetXaxis()->SetTitle("ADC[ch]"); 
  h1->GetYaxis()->SetTitle("TDC[ch]");
  h2->GetXaxis()->SetTitle("ADC[ch]"); 
  h2->GetYaxis()->SetTitle("TDC[ch]");
 
  f1->GetXaxis()->SetTitle("ADC[ch]"); 
  f1->GetYaxis()->SetTitle("TDC[ch]");
  f2->GetXaxis()->SetTitle("ADC[ch]"); 
  f2->GetYaxis()->SetTitle("TDC[ch]");

 c5->cd(2);
   h2->ProfileX();
   h2_pfx->Draw();
   h2_pfx->Fit("pol6","","",fitMin,fitMax); 
   pol6->GetParameters(p1);

 c6->cd(2);
   f2->ProfileX();
   f2_pfx->Draw();
   f2_pfx->Fit("pol6","","",fitMin2,fitMax2); 
   pol6->GetParameters(q1);


 for(long I=0; I<N; I++){

	tr->GetEntry(I);

	if(tdc[0]<tdc0M ){

    P1 =  p1[6]*pow(adc[0],6) + p1[5]*pow(adc[0],5) + p1[4]*pow(adc[0],4) 
	    + p1[3]*pow(adc[0],3) + p1[2]*pow(adc[0],2) + p1[1]*adc[0] + p1[0];

	for(int j=0; j<J; j ++){

		STA = start + d*j;
	    STO = start + d*(j+1);

if (tdc[6]-tdc[0]>STA && tdc[6]-tdc[0]<STO 
		&& tdc[0]<tdc0M/// && tdc[1]<tdc1M 
		&& ( (tdc[14]>sc0m && tdc[15]>sc1m) || (tdc[12]>sc0m && tdc[13]>sc1m) )    
		&& ( (tdc[14]<sc0M && tdc[15]<sc3M) || (tdc[12]<sc0M && tdc[13]<sc1M) )    
//	    && ( (adc[14]>150 && adc[15]>150 && adc[14]<250 && adc[15]<250) 
//	    || (adc[12]>150 && adc[13]>150 && adc[12]<250 && adc[13]<250) )
	    && ( (adc[14]>150 && adc[15]>150 && adc[14]<500 && adc[15]<500) 
	    || (adc[12]>150 && adc[13]>150 && adc[12]<500 && adc[13]<500) )
		){


	H2->Fill(adc[0],tdc[6]-tdc[0]-d*j-P1-start);
	hh2->Fill(tdc[6]-tdc[0]-d*j-P1-start);

}
	

	}

	}



	if(tdc[1]<tdc1M ){

	Q1 =  q1[6]*pow(adc[1],6) + q1[5]*pow(adc[1],5) + q1[4]*pow(adc[1],4) 
	    + q1[3]*pow(adc[1],3) + q1[2]*pow(adc[1],2) + q1[1]*adc[1] + q1[0];

		for(int j=0; j<J; j ++){

	
		STA2 = start2 + d2*j;
	    STO2 = start2 + d2*(j+1);

if (tdc[6]-tdc[1]>STA2 && tdc[6]-tdc[1]<STO2 
		&& tdc[0]<tdc0M && tdc[1]<tdc1M
		&& ( (tdc[14]>sc0m && tdc[15]>sc1m) || (tdc[12]>sc0m && tdc[13]>sc1m) )    
		&& ( (tdc[14]<sc0M && tdc[15]<sc3M) || (tdc[12]<sc0M && tdc[13]<sc1M) )    
//	    && ((adc[14]>200 && adc[15]>200 && adc[14]<300 && adc[15]<300) 
//	    || (adc[12]>200 && adc[13]>200 && adc[12]<300 && adc[13]<300) )
//	    && ( (adc[14]>150 && adc[15]>150 && adc[14]<250 && adc[15]<250) 
//	    || (adc[12]>150 && adc[13]>150 && adc[12]<250 && adc[13]<250) )
	    && ( (adc[14]>150 && adc[15]>150 && adc[14]<500 && adc[15]<500) 
	    || (adc[12]>150 && adc[13]>150 && adc[12]<500 && adc[13]<500) )
		){

	F2->Fill(adc[1],tdc[6]-tdc[1]-d2*j-Q1-start2);
	ff2->Fill(tdc[6]-tdc[1]-d2*j-Q1-start2);

}

		}}


		for(int j=0; j<J; j ++){

		STA = start + d*j;
	    STO = start + d*(j+1);
		STA2 = start2 + d2*j;
	    STO2 = start2 + d2*(j+1);

if (tdc[6]-tdc[0]>STA && tdc[6]-tdc[0]<STO && tdc[6]-tdc[1]>STA2 && tdc[6]-tdc[1]<STO2
		&& tdc[0]<tdc0M && tdc[1]<tdc1M 
		&& ( (tdc[14]>sc0m && tdc[15]>sc1m) || (tdc[12]>sc0m && tdc[13]>sc3m) )    
		&& ( (tdc[14]<sc0M && tdc[15]<sc1M) || (tdc[12]<sc0M && tdc[13]<sc3M) )    
//	    && ((adc[14]>200 && adc[15]>200 && adc[14]<300 && adc[15]<300) 
//	    || (adc[12]>200 && adc[13]>200 && adc[12]<300 && adc[13]<300) )
//	    && ( (adc[14]>150 && adc[15]>150 && adc[14]<250 && adc[15]<250) 
//	    || (adc[12]>150 && adc[13]>150 && adc[12]<250 && adc[13]<250) )
	    && ( (adc[14]>150 && adc[15]>150 && adc[14]<500 && adc[15]<500) 
	    || (adc[12]>150 && adc[13]>150 && adc[12]<500 && adc[13]<500) )
     ){

	
	g->Fill(tdc[6]-((tdc[0]+d*j+P1+start)+(tdc[1]+d2*j+Q1+start2))/2);
	
}


		}
	
 }




 c5->cd(3);
    H2->Draw();
 c6->cd(3);
    F2->Draw();

  H2->GetXaxis()->SetTitle("ADC[ch]"); 
  H2->GetYaxis()->SetTitle("TDC[ch]");
  F2->GetXaxis()->SetTitle("ADC[ch]"); 
  F2->GetYaxis()->SetTitle("TDC[ch]");

 c5->cd(4);
    hh2->Draw();
	hh2->Fit("gaus");
	hh2->GetXaxis()->SetTitle("TDCrpc-rf[ch]");
	hh2->GetYaxis()->SetTitle("counts");
    gaus->GetParameters(res1);
    rese1 = gaus->GetParError(2);
    RES1 = res1[2]*26;
    RESE1 = rese1*26;
	sprintf(r1,"Resolution = %.1f pm %.1f ps",RES1,RESE1);
    TText* tex1 = new TText(-35,2,r1); 
	tex1->SetTextSize(0.07);
	tex1->SetTextColor(8);
    tex1->Draw();  

 c6->cd(4);
    ff2->Draw();
	ff2->Fit("gaus");
	ff2->GetXaxis()->SetTitle("TDCrpc-rf[ch]");
	ff2->GetYaxis()->SetTitle("counts");
    gaus->GetParameters(re);
    rse = gaus->GetParError(2);
    RE = re[2]*26;
    RSE = rse*26;
	sprintf(R,"Resolution = %.1f pm %.1f ps",RE,RSE);
    TText* tex2 = new TText(-38,2,R); 
	tex2->SetTextSize(0.07);
	tex2->SetTextColor(8);
    tex2->Draw();  

 c7->cd(1);
    g->Draw();
	g->Fit("gaus");
	g->GetXaxis()->SetTitle("TDCrpc-rf[ch]");
	g->GetYaxis()->SetTitle("counts");
    gaus->GetParameters(re2);
    rse2 = gaus->GetParError(2);
    RE2 = re2[2]*26;
    RSE2 = rse2*26;
	sprintf(R2,"Resolution = %.1f pm %.1f ps",RE2,RSE2);
    TText* tex22 = new TText(-38,2,R2); 
	tex22->SetTextSize(0.07);
	tex22->SetTextColor(8);
    tex22->Draw();  
	}