/*
  ScalerAna.cc

  2004/12/21

*/

#include "ScalerAna.hh"

#include <cstdio>
#include <iostream>
#include <sstream>
#include <iomanip>

const int BufSize = 144;

bool ScalerAna::Initialize( void )
{
  char str1[BufSize], str2[BufSize];
  FILE *fp;
  bool retval = true;
  int n, n1, n2, seg;

  if( ScalerDefinitionFileName_=="" )
    return false;

  if((fp=fopen(ScalerDefinitionFileName_.c_str(),"r"))){
    while( fgets(str1,BufSize,fp)!=0 ){
      if( str1[0]!='#' ){
	if((n=sscanf(str1,"%d %s",&n1,str2))>=1 ){
	  if( n1>0 && n1<=NumOfScalers ){
	    if(n==2) Name_[n1-1]=str2;
	    else     Name_[n1-1]="";
	  }
	}
      }
    }
    fclose(fp);
    retval=true;
  }
  else{ retval=false; }
  return retval;
}

void ScalerAna::ResetScaler( void )
{
  Counts_=0;
  for( int i=0; i<NumOfScalers; ++i )
    Sum_[i]=Mean_[i]=0.0;
}

bool ScalerAna::Processing( const unsigned int *buf )
{
  ++Counts_;
  for( int i=0; i<NumOfScalers; ++i ){
    Sum_[i] += double(buf[i+1]);
    Mean_[i] = Sum_[i]/double(Counts_);
  }
  return true;
}

void ScalerAna::PrintFile( int runnum ) const
{
  static const std::string funcname = "[ScalerAna::PrintFile(int)]";
  std::ostringstream ost;
  ost << "scaler/scaler" << std::setfill('0') << std::setw(7) << runnum 
      << ".out";
  PrintFile(ost.str());
  return;
}
    
void ScalerAna::PrintFile( const std::string & filename ) const
{
  static const std::string funcname =
    "[ScalerAna::PrintFile(const std::string &)]";
  FILE *fp;

  if((fp=fopen(filename.c_str(),"w"))==0){
    std::cerr << funcname << ": file open fail." << std::endl;
    return;
  }

  fprintf(fp,"Scaler Output :: Total #Scaler Trig. %d\n",Counts_);
  fprintf(fp,
          "========================== Scaler  ==========================\n");
  /*          1234 12345678901234567890 123456789012 1234567890123 */
  fprintf(fp,"Num  Name                 Total        Mean         \n");
  for(int i=1; i<=NumOfScalers; ++i){
    fprintf(fp,"%3d: %-20s %12.0lf %13.2lf\n",
            i,Name_[i-1].c_str(),Sum_[i-1],Mean_[i-1]);
  }
  fclose(fp);
}
