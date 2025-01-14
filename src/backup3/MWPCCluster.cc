// -*- C++ -*-
/*
  MWPCCluster.cc
*/

#include <algorithm>
#include <set>
#include <cmath>
#include <limits>
#include <iostream>

#include "TemplateLib.hh"
#include "DCHit.hh"

#include "MWPCCluster.hh"

#ifdef MemoryLeak
debug::Counter MWPCCluster::sm_counter("MWPCCluster");
#endif

// #ifdef MemoryLeak
// debug::Counter MWPCCluster::sm_counter("MWPCCluster");
// #endif

//______________________________________________________________________________
template <typename T>
bool
equal(const T& lhs,
      const T& rhs)
{
  return ((std::fabs(lhs-rhs)/std::max(std::fabs(lhs), std::fabs(rhs)))
	  < std::numeric_limits<T>::epsilon());
}
//______________________________________________________________________________
void
calcFirst(const std::vector<DCHit*>& hits,
	  MWPCCluster::Statistics& first)
{
  int nHits = hits.size();
  if (nHits==0)
    return;

  for (int i=0; i<nHits; ++i)
    {
      const DCHit* h = hits[i];
      if (!h)
	continue;
      double le = h->GetDriftTime(0);
      if (i==0)
	{
	  first.m_wire     = h->GetWire();
	  first.m_wpos     = h->GetWirePosition();
	  first.m_leading  = le;
	  first.m_trailing = h->GetTrailingTime(0);
	  first.m_length   = first.m_trailing - le;
	}
      else if (le<first.m_leading) // update first hit 
	{
	  first.m_wire     = h->GetWire();
	  first.m_wpos     = h->GetWirePosition();
	  first.m_leading  = le;
	  first.m_trailing = h->GetTrailingTime(0);
	  first.m_length   = first.m_trailing - le;
	}
    }

  double wire     = 0;
  double wire_pos = 0;
  double trailing = 0;
  double length   = 0;
  
  std::set<int> wires;
  for (int i=0; i<nHits; ++i)
    {
      const DCHit* h = hits[i];
      if (!h)
	continue;
      double le = h->GetDriftTime(0);
      if (equal(le, first.m_leading))
	{
	  double tr   = h->GetTrailingTime();
	  double len  = tr - le;
	  double w    = h->GetWire();
	  wire       += w*len;
	  wire_pos   += h->GetWirePosition()*len;
	  trailing   += tr*len;
	  length     += len;
	  wires.insert(static_cast<int>(w));
	}
    }

  if (!wires.empty())
    {
      first.m_wire        = wire/length;
      first.m_wpos        = wire_pos/length;
      first.m_trailing    = trailing/length;
      first.m_length      = length/wires.size();
      first.m_totalLength = length;
      std::set<int>::iterator imin 
	= std::min_element(wires.begin(), wires.end());
      std::set<int>::iterator imax 
	= std::max_element(wires.begin(), wires.end());
      first.m_clusterSize = *imax - *imin + 1;
    }
  else
    {
      first.m_totalLength = first.m_length;
      first.m_clusterSize = 1;
    }
  
  return;
}

//______________________________________________________________________________
void
calcMean(const std::vector<DCHit*>& hits,
	 MWPCCluster::Statistics& mean)
{
  int nHits = hits.size();
  if (nHits==0)
    return;
  double wire     = 0;
  double wire_pos = 0;
  double leading  = 0;
  double trailing = 0;
  double length   = 0;

  std::set<int> wires;

  for (int i=0; i<nHits; ++i)
    {
      const DCHit* h = hits[i];
      if (!h)
	continue;
      double w    = h->GetWire();
      double wpos = h->GetWirePosition();
      double le   = h->GetDriftTime(0);
      double tr   = h->GetTrailingTime(0);
      double len  = tr - le;

      length   += len;
      leading  += le*len;
      trailing += tr*len;
      wire     += w*len;
      wire_pos += wpos*len;
      wires.insert(static_cast<int>(w));
    }
  
  mean.m_wire     = wire    /length;
  mean.m_wpos     = wire_pos/length;
  mean.m_leading  = leading /length;
  mean.m_trailing = trailing/length;
  mean.m_length   = length/wires.size();
  mean.m_totalLength = length;

  std::set<int>::iterator imin = std::min_element(wires.begin(), wires.end());
  std::set<int>::iterator imax = std::max_element(wires.begin(), wires.end());
  mean.m_clusterSize = *imax - *imin + 1;
  return;
}

//______________________________________________________________________________
// struct MWPCCluster::Statistics
//______________________________________________________________________________
MWPCCluster::
Statistics::Statistics()
  : m_wire(-0xffff),
    m_wpos(-0xffff),
    m_leading(-0xffff),
    m_trailing(-0xffff),
    m_length(-0xffff),
    m_totalLength(-0xffff),
    m_clusterSize(-0xffff)
{
}

//______________________________________________________________________________
MWPCCluster::
Statistics::~Statistics()
{
}

//______________________________________________________________________________
void
MWPCCluster::
Statistics::Print(const std::string& arg) const
{
  std::cout << "[MWPCCluster::Statistics::Print] " << arg
	    << "\n  wire   = " << m_wire
	    << "\n  wpos   = " << m_wpos        << " [mm]"
	    << "\n  dt     = " << m_leading     << " [nsec]" 
	    << "\n  trail  = " << m_trailing    << " [nsec]"
	    << "\n  siglen = " << m_length      << " [nsec]"
	    << "\n  total  = " << m_totalLength << " [nsec]"
	    << "\n  csize  = " << m_clusterSize
	    << std::endl;
  return;
}


//______________________________________________________________________________
// class MWPCCluster
//______________________________________________________________________________
MWPCCluster::MWPCCluster()
  : m_hits(0),
    m_mean(),
    m_first(),
    m_status(false)
{
#ifdef MemoryLeak
  ++sm_counter;
#endif
}

//______________________________________________________________________________
MWPCCluster::~MWPCCluster()
{
  std::for_each(m_hits.begin(), m_hits.end(), DeleteObject());
  m_hits.clear();
#ifdef MemoryLeak
  --sm_counter;
#endif
}

//______________________________________________________________________________
void
MWPCCluster::Add(DCHit* h)
{
  m_hits.push_back(h);
  return;
}

//______________________________________________________________________________
void
MWPCCluster::Calculate()
{
  if (m_hits.empty())
    return;

  calcMean(m_hits, m_mean);
  calcFirst(m_hits, m_first);

  return;
}

//______________________________________________________________________________
int
MWPCCluster::GetClusterSize() const
{
  return m_mean.m_clusterSize;
}

//______________________________________________________________________________
const MWPCCluster::Statistics&
MWPCCluster::GetFirst() const
{
  return m_first;
}

//______________________________________________________________________________
const std::vector<DCHit*>&
MWPCCluster::GetHits() const
{
  return m_hits;
}

//______________________________________________________________________________
const MWPCCluster::Statistics&
MWPCCluster::GetMean() const
{
  return m_mean;
}

//______________________________________________________________________________
double
MWPCCluster::GetMeanTime() const
{
  return m_mean.m_leading;
}

//______________________________________________________________________________
double
MWPCCluster::GetMeanWire() const
{
  return m_mean.m_wire;
}

//______________________________________________________________________________
double
MWPCCluster::GetMeanWirePos() const
{
  return m_mean.m_wpos;
}

//______________________________________________________________________________
int
MWPCCluster::GetNumOfHits() const
{
  return m_hits.size();
}

//______________________________________________________________________________
void
MWPCCluster::Print(const std::string& arg) const
{
  std::cout << "[MWPCCluster::Print] " << arg << std::endl;
  std::cout << " nhits = " << m_hits.size() << std::endl;
  m_mean.Print("mean");
  m_first.Print("first");
  return;
}

//______________________________________________________________________________
bool
MWPCCluster::IsGoodForAnalysis() const
{
  return m_status;
}

//______________________________________________________________________________
void
MWPCCluster::SetStatus(bool status)
{
  m_status = status;
  return;
}


// #include "MWPCCluster.hh"
// #include "DCHit.hh"
// //#include "SksObjectId.hh"
// #include "DCAnalyzer.hh"

// #include <cstring>
// #include <string>
// #include <iostream>
// #include <stdexcept>
// #include <sstream>

// MWPCCluster::MWPCCluster( DCHit *hitA, DCHit *hitB, DCHit *hitC,
// 			  DCHit *hitD, DCHit *hitE)
//   : hitA_(hitA), hitB_(hitB), hitC_(hitC), 
//     hitD_(hitD), hitE_(hitE), 
//     csize_(0), gfastatus_(true)
// {
//   if(hitA) ++csize_;
//   if(hitB) ++csize_;
//   if(hitC) ++csize_;
//   if(hitD) ++csize_;
//   if(hitE) ++csize_;

//   calculate();
// }

// DCHit * MWPCCluster::GetHit( int i ) const
// {
//   if( i==0 ) return hitA_;
//   else if( i==1 ) return hitB_;
//   else if( i==2 ) return hitC_;
//   else if( i==3 ) return hitD_;
//   else if( i==4 ) return hitE_;

//   else return 0;
// }

// void MWPCCluster::calculate( void )
// {
//   double mw=0., mwp=0., mt=0., dt=0.;
//   if( hitA_ ){
//     mw += hitA_->GetWire();
//     mwp+= hitA_->GetWirePosition();
//     //std::cout<<"CWireA="<< hitA_->GetWire() <<std::endl;
//     int multiA = hitA_->GetDriftTimeSize();
//     int nhitsA=0;
//     for (int mA=0; mA<multiA; mA++) {
//       if( !(hitA_->rangecheck(mA)) ) continue;
//       mt += hitA_->GetDriftTime(mA);
//       nhitsA++; if( nhitsA>1 ) continue;
//     }
//     mt /= double(nhitsA);
//   }
//   if( hitB_ ){
//     mw += hitB_->GetWire();
//     mwp+= hitB_->GetWirePosition();
//     //std::cout<<"CWireB="<< hitB_->GetWire() <<std::endl;
//     int multiB = hitB_->GetDriftTimeSize();
//     int nhitsB=0;
//     for (int mB=0; mB<multiB; mB++) {
//       if( !(hitB_->rangecheck(mB)) ) continue;
//       mt += hitB_->GetDriftTime(mB);
//       nhitsB++; if( nhitsB>1 ) continue;
//     }
//     mt /= double(nhitsB);
//   }
//   if( hitC_ ){
//     mw += hitC_->GetWire();
//     mwp+= hitC_->GetWirePosition();
//     //std::cout<<"CWireC="<< hitC_->GetWire() <<std::endl;
//     int multiC = hitC_->GetDriftTimeSize();
//     int nhitsC=0;
//     for (int mC=0; mC<multiC; mC++) {
//       if( !(hitC_->rangecheck(mC)) ) continue;
//       mt += hitC_->GetDriftTime(mC);
//       nhitsC++; if( nhitsC>1 ) continue;
//     }
//     mt /= double(nhitsC);
//   }
//   if( hitD_ ){
//     mw += hitD_->GetWire();
//     mwp+= hitD_->GetWirePosition();
//     //std::cout<<"CWireD="<< hitD_->GetWire() <<std::endl;
//     int multiD = hitD_->GetDriftTimeSize();
//     int nhitsD=0;
//     for (int mD=0; mD<multiD; mD++) {
//       if( !(hitD_->rangecheck(mD)) ) continue;
//       mt += hitD_->GetDriftTime(mD);
//       nhitsD++; if( nhitsD>1 ) continue;
//     }
//     mt /= double(nhitsD);
//   }
//   mw /= double(csize_);
//   mwp /= double(csize_);
//   mt /= double(csize_);

//   MeanWire_=mw;
//   MeanWirePos_=mwp;
//   MeanTime_=mt;

// //    std::cout<< "***************************"<<std::endl;
// //    std::cout<< "Clustersize=" << csize_ << std::endl;
// //   std::cout<< "MeanWire=" << mw << std::endl;
// //   std::cout<< "MeanWirePos=" << mwp << std::endl;
// //   std::cout<< "MeanTime=" << mt << std::endl;
// } 

// bool MWPCCluster::ReCalc( bool applyRecursively )
// {
//   if( applyRecursively ){
//     if(hitA_) hitA_->ReCalcMWPC(applyRecursively);
//     if(hitB_) hitB_->ReCalcMWPC(applyRecursively);
//     if(hitC_) hitC_->ReCalcMWPC(applyRecursively);
//   }
//   calculate();

//   return true;
// }
