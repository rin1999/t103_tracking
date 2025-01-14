// -*- C++ -*-

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <iterator>
#include <set>

#include "DetectorID.hh"
#include "HodoAnalyzer.hh"
#include "DCHit.hh" 
#include "DCAnalyzer.hh"
#include "Hodo2Hit.hh"
#include "HodoCluster.hh"
#include "BH1Filter.hh"

//______________________________________________________________________________
// struct BH1Filter::Param
//______________________________________________________________________________
BH1Filter::
Param::Param()
  : m_xmin(NumOfLayersBcIn),
    m_xmax(NumOfLayersBcIn)
{
}

//______________________________________________________________________________
BH1Filter::
Param::~Param()
{
}

//______________________________________________________________________________
void
BH1Filter::
Param::Print(const std::string& arg) const
{

  std::cout << "\n " << arg << " (xmin, xmax) = \n";
  for (int i=0, n=m_xmin.size(); i<n; ++i)
    {
      std::cout << " iplane " << std::setw(3) << i 
		<< " (" << std::setw(6) << m_xmin[i] << " "
		<< std::setw(6) << m_xmax[i] << ")\n";
    }
  std::cout << std::endl;
  return;
}

//______________________________________________________________________________
// class BH1Filter
//______________________________________________________________________________
BH1Filter::BH1Filter()
  : m_param(NumOfSegBH1)
{
}

//______________________________________________________________________________
BH1Filter::~BH1Filter()
{
}

//______________________________________________________________________________
void
BH1Filter::Apply(const HodoAnalyzer& hodo,
		 const DCAnalyzer& dc,
		 std::vector<std::vector<DCHitContainer> >& candidates)
{
  // +++++++++++++++++++++++++++++++++++++++++++
  // candidates [segment id] [plane id] [hit id]
  // +++++++++++++++++++++++++++++++++++++++++++

  m_dc   = &dc;
  m_hodo = &hodo;

  std::set<int> seg;
  for (int i=0, n=hodo.GetNHitsBH1(); i<n; ++i)
    {
      const Hodo2Hit* const h = hodo.GetHitBH1(i);
      if (!h)
	continue;
      seg.insert(h->SegmentId());
    }


  candidates.resize(seg.size());
  std::vector<std::vector<DCHitContainer> >::iterator itCont = candidates.begin();
  for (std::set<int>::const_iterator itSeg = seg.begin(), itSegEnd = seg.end();
       itSeg!=itSegEnd; ++itSeg, ++itCont)
    {
      const int iSeg = *itSeg;
      std::vector<DCHitContainer>& c = *itCont;
      c.resize(NumOfLayersBcIn);

//       std::cout << "  BH1 seg = " << iSeg << "\n";
      for (int iplane=0; iplane<NumOfLayersBcIn; ++iplane)
	{
	  DCHitContainer& after = c[iplane];
	  const double xmin = m_param[iSeg].m_xmin[iplane];
	  const double xmax = m_param[iSeg].m_xmax[iplane];
	  int iLayer = iplane + 1;
	  const DCHitContainer& before = dc.GetBcInHC(iLayer);
	  for (int ih=0, nh=before.size(); ih<nh; ++ih)
	    {
	      const DCHit* const h = before[ih];
	      if (!h)
		continue;
	      const double wpos = h->GetWirePosition();
	      const int layer   = h->GetLayer();

// 	      std::cout << " layer = " << iplane 
// 			<< "(" << layer << ") : " << wpos
// 			<< " (" << xmin << ", " << xmax << ")";
	      if (wpos<xmin || xmax<wpos)
		{
// 		  std::cout << std::endl;
		  continue;
		}
// 	      std::cout << " good " << std::endl;
	      after.push_back(const_cast<DCHit*>(h));
	    }
// 	  std::cout << __FILE__ << ":" << __LINE__
// 		    << " " << after.size() << std::endl;
	}
    }
//   std::cout << __FILE__ << ":" << __LINE__ 
// 	    << " " << candidates.size() << std::endl;
  
  return;
}

//______________________________________________________________________________
BH1Filter&
BH1Filter::GetInstance()
{
  static BH1Filter s_instance;
  return s_instance;
}


//______________________________________________________________________________
const std::vector<double>&
BH1Filter::GetXmax(int seg) const
{
  return m_param[seg].m_xmax;
}

//______________________________________________________________________________
const std::vector<double>&
BH1Filter::GetXmin(int seg) const
{
  return m_param[seg].m_xmin;
}

//______________________________________________________________________________
void
BH1Filter::Initialize(const std::string& filename)
{
  static const std::string funcname
    = std::string("[BH1Filter::") + __func__ + "]"; 

  std::ifstream f(filename.c_str());
  if (f.fail())
    {
      std::cerr << "#E " << funcname << " file open fail " 
		<<  filename << std::endl;
      std::exit(-1);
    }
  
  while (f.good())
    {
      std::string l;
      std::getline(f, l);
      if (l.empty() || l.find("#")!=std::string::npos)
	continue;
      std::istringstream iss(l);
      std::istream_iterator<double> issBegin(iss);
      std::istream_iterator<double> issEnd;
      std::vector<double> v(issBegin ,issEnd);
      if (v.size()<kNParam)
	{
// 	  std::cout << "#W " << funcname
// 		    << " number of parameters = " << v.size()
// 		    << ": required = " << kNParam
// 		    << std::endl;
	  continue;
	}
      const int bh1Seg  = static_cast<int>(v[kBH1Segment]);
      const int bcPlane = static_cast<int>(v[kLayerID]);
      const double xmin = v[kXMin];
      const double xmax = v[kXMax];

      int iplane = bcPlane - PlOffsBc -1;
//       std::cout << " seg = "    << std::setw(2) << bh1Seg
// 		<< ", plane = " << std::setw(3) << bcPlane
// 		<< "(" << iplane << ")"
// 		<< ", xmin = "  << std::setw(5) << xmin
// 		<< ", xmax = "  << std::setw(5) << xmax
// 		<< std::endl;
      m_param[bh1Seg].m_xmin[iplane] = xmin;
      m_param[bh1Seg].m_xmax[iplane] = xmax;

    }
//   Print();
  return;
}

//______________________________________________________________________________
void
BH1Filter::Print(const std::string& arg) const
{
  std::cout << "#D BH1Filter::print" << std::endl;
  for (int i=0, n=m_param.size(); i<n; ++i)
    {
      std::stringstream ss;
      ss << "isegment " << i;
      m_param[i].Print(ss.str());
    }
  return;
}
