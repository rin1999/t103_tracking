/*
  TrTrackSearch.cc

  2024/11  K.Shirotori
*/

#include "TrTrackSearch.hh"
#include "TrParameters.hh"
#include "TrLTrackHit.hh"
#include "TrPairHitCluster.hh"
#include "TrLocalTrack.hh"
#include "TemplateLib.hh"
#include "DetectorInfo.hh"
#include "GeomMan.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double MaxChisquare = 1000.;//30
const double MaxChisquareTr = 1000.;//30
const double MaxNumberOfClusters = 10.;//10.
const double MaxCombi = 1.0e6;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
int LocalTrackSearch( const TrHitContainer * HC, std::vector <TrLocalTrack *> &TrackCont, int NumOfLayers, int MinNumOfHits)
{
  static const std::string funcname = "[LocalTrackSearch]";

  std::vector < std::vector <TrPairHitCluster *> > CandCont;
  CandCont.resize(NumOfLayers);

  for( int i=0; i<NumOfLayers; ++i ){
    MakeUnPairPlaneHitCluster( HC[i], CandCont[i] );
  }
  
  std::vector <int> nCombi(NumOfLayers);
  for( int i=0; i<NumOfLayers; ++i ){ 
    nCombi[i]=(CandCont[i]).size();

    // If #Cluster>MaxNumerOfCluster,  error return

    if(nCombi[i]>MaxNumberOfClusters){
      for( int i=0; i<NumOfLayers; ++i )
	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
      return 0;
    } 
  }

#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)CandCont[i][j]->GetHit(0))->GetMeanFiber() << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector < std::vector <int> > 


    CombiIndex = makeindex( NumOfLayers, 
			    MinNumOfHits, 
			    NumOfLayers, 
			    &(nCombi[0]) );
  int nnCombi=CombiIndex.size();

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif
  if( nnCombi>MaxCombi )  return 0;

  for( int i=0; i<nnCombi; ++i ){
    TrLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHits && 
	    track->DoFit() &&
	    track->GetChiSquare()<MaxChisquare ){
      //hitlayersCont.push_back(track->GetHitLayers());
      TrackCont.push_back(track);
      double chisqr = track->GetChiSquare();
    }
    //else if(track->GetNHit()==MinNumOfHits-1 && 
	  //  track->DoFit() &&
	  //  track->GetChiSquare()<MaxChisquare){
    //    hitlayersCont.push_back(track->GetHitLayers());
    //    delete track;
    //}
    else{ // trackが見つからない場合
      //      std::cout << "No tracks available" << std::endl;
      delete track;
    }
  }

  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit(); 
    for( int j=0; j<nh; ++j ) {
      //if(tp->GetHit(j)->showFlags()) 
	tp->GetHit(j)->clearFlags();

    }
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), TrLTrackComp() );

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  // Delete Duplicated Tracks
  for( int i=0; i<int(TrackCont.size()); ++i ){
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      TrLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }

  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = GeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Deleting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  for( int i=0; i<NumOfLayers; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
  return TrackCont.size();
}

bool MakePairPlaneHitCluster( const TrHitContainer & HC1,
			      const TrHitContainer & HC2,
			      double CellSize,
			      std::vector <TrPairHitCluster *> & Cont )
{
  int nh1=HC1.size(), nh2=HC2.size();
  std::vector <int> UsedFlag(nh2,0);   

  for( int i1=0; i1<nh1; ++i1 ){
    TrHit *hit1=HC1[i1];

    double wp1=hit1->GetMPosition();
    bool flag=false;
    for( int i2=0; i2<nh2; ++i2 ){
      TrHit *hit2=HC2[i2];
      double wp2=hit2->GetMPosition();
      if( fabs(wp1-wp2)<CellSize ){

	int multi1 = hit1->GetDriftLengthSize();
	int multi2 = hit2->GetDriftLengthSize();
	for (int m1=0; m1<multi1; m1++) {
	  if( !(hit1->rangecheck(m1)) ){
	    continue;
	  }
	  for (int m2=0; m2<multi2; m2++) {
	    if( !(hit2->rangecheck(m2)) ){
	      continue;
	    }
	    double x1,x2;
	    if( wp1<wp2 ){
	      x1=wp1+hit1->GetDriftLength(m1);
	      x2=wp2-hit2->GetDriftLength(m2);
	    }
	    else {
	      x1=wp1-hit1->GetDriftLength(m1);
	      x2=wp2+hit2->GetDriftLength(m2);
	    }
	    Cont.push_back( new TrPairHitCluster( new TrLTrackHit(hit1,x1,m1),
						  new TrLTrackHit(hit2,x2,m2) ) );
	    flag=true; ++UsedFlag[i2];
	  }
	}
      }
    }      
#if 1
    if(!flag){
      int multi1 = hit1->GetDriftLengthSize();
      for (int m1=0; m1<multi1; m1++) {
	if( !(hit1->rangecheck(m1)) ) continue;
	double dl=hit1->GetDriftLength(m1);
	Cont.push_back( new TrPairHitCluster( new TrLTrackHit(hit1,wp1+dl,m1) ) );
	Cont.push_back( new TrPairHitCluster( new TrLTrackHit(hit1,wp1-dl,m1) ) );
      }
    }
#endif
  }
#if 1
  for( int i2=0; i2<nh2; ++i2 ){
    if( UsedFlag[i2]==0 ) {
      TrHit *hit2=HC2[i2];
      int multi2 = hit2->GetDriftLengthSize();
      for (int m2=0; m2<multi2; m2++) {
	if( !(hit2->rangecheck(m2)) ) continue;
	
	double wp=hit2->GetMPosition();
	double dl=hit2->GetDriftLength(m2);
	Cont.push_back( new TrPairHitCluster( new TrLTrackHit(hit2,wp+dl,m2) ) );
	Cont.push_back( new TrPairHitCluster( new TrLTrackHit(hit2,wp-dl,m2) ) );
      }
    }
  }
#endif
  return true;
}

bool MakeUnPairPlaneHitCluster( const TrHitContainer & HC,
				std::vector <TrPairHitCluster *> & Cont )
{
  int nh=HC.size();

  for( int i=0; i<nh; ++i ){
    TrHit *hit=HC[i];
    if( hit ){
      int multi = hit->GetDriftLengthSize();

      for (int m=0; m<multi; m++) {
	if( !(hit->rangecheck(m)) ) continue;

	double wp=hit->GetMPosition();
	double dl=hit->GetDriftLength(m);

	Cont.push_back( new TrPairHitCluster( new TrLTrackHit(hit,wp,m) ) );
	// Cont.push_back( new TrPairHitCluster( new TrLTrackHit(hit,wp+dl,m) ) );
	// Cont.push_back( new TrPairHitCluster( new TrLTrackHit(hit,wp-dl,m) ) );
      }
    }
  }

  return true;
}

std::vector< std::vector<int> > 
makeindex( int ndim_org, int minimumHit, int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex( ndim_org, minimumHit, ndim-1, index1+1 );
 
  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int validHitNum=0;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if (i != -1)
	validHitNum++;
      for( int k=0; k<n3; ++k ) {
        elem.push_back(index2[j][k]);
        if (index2[j][k] != -1)
          validHitNum++;
      }
      if (ndim==ndim_org) {
        if (validHitNum >= minimumHit)
          index.push_back(elem);
      } else {
        index.push_back(elem);
      }
      int size1=index.size();
    }
  }

  return index;
}

TrLocalTrack 
* MakeTrack(  const std::vector < std::vector <TrPairHitCluster *> > &CandCont,
	      const int *combination )
{
  static const std::string funcname = "[MakeTrack]";

  TrLocalTrack *tp=new TrLocalTrack();
  //TrLocalTrack *tp ;

  if(!tp){
    std::cerr << funcname << ": new fail" << std::endl;
    return 0;
  }    

  int n=CandCont.size();

  for( int i=0; i<n; ++i ){
    int m=combination[i];
    TrPairHitCluster *cluster=0;
    if(m>=0) cluster=CandCont[i][m];
#if 0
    std::cout << funcname << ":" << std::setw(3)
	      << i << std::setw(3) << m  << " "
	      << CandCont[i][m] << std::endl; 
#endif


#if 0

    /*****************************************************************/
    // Only SdcIn & BcOut Local Tracking !!                          
    /*****************************************************************/
    int mm = 0 ;
    //std::cout << "-------------------------------------" << std::endl;
    if(cluster){
     mm = cluster->NumberOfHits();
     if(mm>2) std::cout << "mm = " << mm << std::endl;

     double DL[mm] ;
     double DT[mm] ;
     int Layer[mm] ;

      for(int j=0; j<mm; ++j ){
	TrLTrackHit *hitp=cluster->GetHit(j);
#if 0
	std::cout << "mm = " << mm 
	          << " Layer = " << hitp->GetLayer() 
		  << " Drift Length = " << hitp->GetDriftLength()
		  << std::endl;
#endif
	DL[j] = hitp->GetDriftLength();
	DT[j] = hitp->GetDriftTime();
	Layer[j] = hitp->GetLayer();

	if(mm<2 || !((Layer[j]>0 && Layer[j]<11) || (Layer[j]>112 && Layer[j]<125))){
	  if((Layer[j]>=1 && Layer[j]<=4) && (DT[j]>-10 && DT[j]<40))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=5 && Layer[j]<=10) && (DT[j]>-10 && DT[j]<60))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=113 && Layer[j]<=118) && (DT[j]>-10 && DT[j]<40))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=119 && Layer[j]<=124) && (DT[j]>-10 && DT[j]<60))
	    if(hitp) tp->AddHit( hitp );
	}
      }
      //std::cout << "***********************************" << std::endl;
      if(mm==2 && ((Layer[0]>0 && Layer[0]<11) || (Layer[0]>112 && Layer[0]<125))){
	int Ok = 0 ;
	/*	
	std::cout << "DL[ " << Layer[0] << "] = " << DL[0] 
		  << " DL[ " << Layer[1] << "] = " << DL[1] << std::endl; 
	*/
	//std::cout << "***********************************" << std::endl;
	if(Layer[0] ==1 || Layer[0]==3 || Layer[0]==113 || Layer[0]==115 || Layer[0]==117){
	  if(((DL[0]+DL[1])>0.7) && ((DL[0]+DL[1])<2.3)) Ok = 1 ;
	}  
	if(Layer[0] ==5 || Layer[0]==7 || Layer[0]==9 || Layer[0]==119 || Layer[0]==121 || Layer[0]==123){
	  if(((DL[0]+DL[1])>1.5) && ((DL[0]+DL[1])<3.5)) Ok = 1 ;
	}
	if(Ok){
	  for(int j=0; j<mm; ++j ){
	    TrLTrackHit *hitp=cluster->GetHit(j);
	    if(hitp) tp->AddHit( hitp );
	  }
	}
      }
    }
#endif


#if 1
    if(cluster){
     int mm=cluster->NumberOfHits();

      for(int j=0; j<mm; ++j ){
	TrLTrackHit *hitp=cluster->GetHit(j);

	  if(hitp) tp->AddHit( hitp );
      }
    }
#endif

  }
  return tp;
}


