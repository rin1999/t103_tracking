#ifndef DECODERH_
#define DECODERH_

#include<fstream>
#include<bitset>
#include<vector>
#include <boost/iostreams/filtering_stream.hpp>
#include<TGraphErrors.h>
#include"Datadrs.hh"
#include"MznHRTdcData.hh"

class Decoder{

public:
  bool Open(const char* name);
  bool isGood(int imod){return flag[imod][Good];};
  bool isEOF(int imod){return flag[imod][EoF];};

  bool getNextEvent();
  bool decode(DRS::dataDrs& obj);
  
  //Fiber
  void GetFiberData( int adch[128], int adcl[128], 
		     std::vector<std::vector<int> >& tdcl, 
		     std::vector<std::vector<int> >& tdct );
  
  //DRS4
  bool HasTDCData() const { return has_tdc_data; }
  void SetDecodeTDC(bool flag=true) { has_tdc_data = flag; }
  bool decodeTDC(DRS::dataDrs& obj);
  bool decodeADC(DRS::dataDrs& obj);

  //HR-TDC
  MZN_HRTDC::Data_t GetMznHRTdcData(){return sMZNData;}

  static Decoder& getInstance();

private:
  enum IndexBaseHeader{
    i_magic, i_datasize, i_tag, i_tic,
    sizeBaseHeader
  };

  enum IndexExHeader{
    i_drs0, i_drs1, i_drs2, i_drs3,
    sizeExHeader
  };

  typedef std::vector<unsigned int>           bufType;
  typedef std::vector<unsigned int>::iterator itrType;

  bufType head_buffer_[DRS::NofBoards];
  bufType exhead_buffer_[DRS::NofBoards];
  bufType body_buffer_[DRS::NofBoards];
  bufType tdc_buffer_[DRS::NofBoards];
  
  enum Flags{ Good, EoF, ExHead, FirstData, DoubleData,
	      sizeFlags};
  std::bitset<sizeFlags> flag[DRS::NofBoards];

  boost::iostreams::filtering_istream m_ifs;

  bool has_tdc_data;

  int read(bufType& buf, int nword);
  bool nReadDRS;

  Decoder();
  ~Decoder();
  Decoder(const Decoder& obj);
  Decoder& operator=(const Decoder& obj);

  // base line restorer
  struct dataTic{
    bool fl_read;
    int time_interval;
    bool fl_overflow;
  };

  typedef std::vector<dataTic> ticType;
  ticType contTic[DRS::NofDrs][DRS::NofBoards]; // for each DRS4 chip

  void   resetTicCont();
  void   setThisEvent(int tic, bool of);
  
  MZN_HRTDC::Data_t sMZNData;

  //functions for NIMEASIROC 
  static const int FiberNch = 128;
  unsigned int getBigEndian32(const char* b);
  unsigned int Decode32bitWord(unsigned int word32bit);
  bool isAdcHg(unsigned int);
  bool isAdcLg(unsigned int);
  bool isTdcLeading(unsigned int);
  bool isTdcTrailing(unsigned int);
  bool isScaler(unsigned int);
  
  //fiber data
  int fiberadch_[128];
  int fiberadcl_[128];
  std::vector<std::vector<int> > fibertdcl_;
  std::vector<std::vector<int> > fibertdct_;
};

// getInstance --------------------------------------------------------------
inline
Decoder&
Decoder::getInstance()
{
  static Decoder obj;
  return obj;
}

#endif
