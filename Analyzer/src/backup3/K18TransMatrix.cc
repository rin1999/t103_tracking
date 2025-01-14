/*
  K18TransMatrix.cc
*/

#include "K18TransMatrix.hh"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <vector>
#include <cstdlib>

#include "lexical_cast.hh"

const double MMtoM = 1.E-3;
const double MtoMM = 1000.;

#define DEBUG 0

using namespace hddaq;

bool K18TransMatrix::Initialize( void )
{
  static const std::string funcname = "[K18TansMatrix::Initialize]";
  std::vector<std::vector<double> > parameters;
  parameters.resize(4);
  std::ifstream file( filename_.c_str() );

  if(file.fail()){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }

  while(file.good())
    {
      std::string buf;
      std::getline(file,buf);

      std::istringstream is(buf);
      std::istream_iterator<std::string> issBegin(is);
      std::istream_iterator<std::string> issEnd;
      std::vector<std::string> param(issBegin,issEnd);
      if(param.empty())continue;
      if(param[0].empty())continue;
      if(param[0]=="X" && param[1]!="A")
	{
	  for(int i=1;i<param.size();++i)
	    {
	      double p = a2d(param[i]);
	      parameters[0].push_back(p);
	    }
	}
      if(param[0]=="A")
	{
	  for(int i=1;i<param.size();++i)
	    {
	      double p = a2d(param[i]);
	      parameters[1].push_back(p);
	    }
	}
      if(param[0]=="Y" && param[1]!="B")
	{
	  for(int i=1;i<param.size();++i)
	    {
	      double p = a2d(param[i]);
	      parameters[2].push_back(p);
	    }
	}
      if(param[0]=="B")
	{
	  for(int i=1;i<param.size();++i)
	    {
	      double p = a2d(param[i]);
	      parameters[3].push_back(p);
	    }
	}
    }
  file.close();

  for(unsigned int i=0;i<parameters[0].size();++i)
    {
      Xpar[i] = parameters[0][i];
    }
  for(unsigned int i=0;i<parameters[1].size();++i)
    {
      Upar[i] = parameters[1][i];
    }
  for(unsigned int i=0;i<parameters[2].size();++i)
    {
      Ypar[i] = parameters[2][i];
    }
  for(unsigned int i=0;i<parameters[3].size();++i)
    {
      Vpar[i] = parameters[3][i];
    }
#if DEBUG
  for(unsigned int i=0;i<parameters[0].size();++i)
    {
      std::cout<<"Xpar[i]"<<" ";
    }
  std::cout<<std::endl;
  for(unsigned int i=0;i<parameters[1].size();++i)
    {
      std::cout<<"Upar[i]"<<" ";
    }
  std::cout<<std::endl;
  for(unsigned int i=0;i<parameters[2].size();++i)
    {
      std::cout<<"Ypar[i]"<<" ";
    }
  std::cout<<std::endl;
  for(unsigned int i=0;i<parameters[3].size();++i)
    {
      std::cout<<"Vpar[i]"<<" ";
    }
  std::cout<<std::endl;
#endif
  std::cout << funcname << ": Initialization finished" << std::endl;
  return true;
}

bool K18TransMatrix::
Transport( double xin, double yin, double uin, double vin, double delta,
	   double & xout, double & yout, double & uout, double & vout ) const
{
  xin *= -MMtoM; yin *= -MMtoM; uin = -uin; vin = -vin;
  xout=yout=uout=vout=0.0;

  double bj1[31], bj2[8];
  bj1[ 0] = xin;     
  bj1[ 1] = uin; 
  bj1[ 2] = delta; 
  bj1[ 3] = xin*xin;
  bj1[ 4] = xin*uin; 
  bj1[ 5] = xin*delta;
  bj1[ 6] = uin*uin;
  bj1[ 7] = uin*delta;
  bj1[ 8] = delta*delta;
  bj1[ 9] = yin*yin;
  bj1[10] = yin*vin;
  bj1[11] = vin*vin;
  bj1[12] = xin*xin*xin;
  bj1[13] = xin*xin*uin;
  bj1[14] = xin*xin*delta;
  bj1[15] = xin*uin*uin;
  bj1[16] = xin*uin*delta;
  bj1[17] = xin*delta*delta;
  bj1[18] = xin*yin*yin;
  bj1[19] = xin*yin*vin;
  bj1[20] = xin*vin*vin;
  bj1[21] = uin*uin*uin;
  bj1[22] = uin*uin*delta;
  bj1[23] = uin*delta*delta;
  bj1[24] = uin*yin*yin;
  bj1[25] = uin*yin*vin;
  bj1[26] = uin*vin*vin;
  bj1[27] = delta*delta*delta;
  bj1[28] = delta*yin*yin;
  bj1[29] = delta*yin*vin;
  bj1[30] = delta*vin*vin;

  bj2[ 0] = yin;
  bj2[ 1] = vin;
  bj2[ 2] = yin*xin;
  bj2[ 3] = yin*uin;
  bj2[ 4] = yin*delta;
  bj2[ 5] = vin*xin;
  bj2[ 6] = vin*uin;
  bj2[ 7] = vin*delta;

  for( int i=0; i<31; ++i ){
    xout += Xpar[i]*bj1[i];
    uout += Upar[i]*bj1[i];
  }

  for( int i=0; i<8; ++i ){
    yout += Ypar[i]*bj2[i];
    vout += Vpar[i]*bj2[i];
  }

  xout *= -MtoMM; yout *= -MtoMM; uout = -uout; vout = -vout;
  return true;
}

