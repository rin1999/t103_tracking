#include<ios>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<algorithm>
#include<functional>
#include<arpa/inet.h>
#include <cstdint>

#include "Decoder.hh"
#include "HexDump.hh"

#include "ConfMan.hh"

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#define BLCOR 0
#define VER2  0

static const std::string MyName = "Decoder";

const double kTimeUnit      = 1.0/960e6/1e-9; // nsec
const double kAmplitudeUnit = 0.445; // mV
const double kChargeUnit    = kTimeUnit*kAmplitudeUnit/50; // pC
const double kTDCTimeUnit   = 1.0;

struct FiberHit{
  int layer;
  int fiber;
  int easiroc;
  int type;
  double pos;
  int adchigh;
  bool otradchigh;
  int adclow;
  bool otradclow;
  std::vector<int> tdcleading;
  std::vector<int> tdctrailing;

  FiberHit(){
    layer = -1;
    fiber = -1;
    easiroc = -1;
    type = -1;
    pos = -99.;
    adchigh = -1;
    otradchigh = false;
    adclow = -1;
    otradclow = false;
    tdcleading.clear(); 
    tdctrailing.clear();
  }

  void clear(){
    layer = -1;
    fiber = -1;
    easiroc = -1;
    type = -1;
    pos = -99.;
    adchigh = -1;
    otradchigh = false;
    adclow = -1;
    otradclow = false;
    tdcleading.clear(); 
    tdctrailing.clear();
  }

};

//______________________________________________________________________________
unsigned int
decodeTdcWord(unsigned int v,
              unsigned int& ct,
              unsigned int& ft)
{
  unsigned int data_type = (v >> DRS::ktdc_data_type_shift) & DRS::ktdc_data_type_mask;
  ct = (v >> DRS::ktdc_coarse_time_shift) & DRS::ktdc_coarse_time_mask;
  ft = (v >> DRS::ktdc_fine_time_shift) & DRS::ktdc_fine_time_mask;
  return data_type;
}

//______________________________________________________________________________
unsigned int
decodeTdcWord(unsigned int v,
              unsigned int& ch,
              unsigned int& ct,
              unsigned int& ft)
{
  ch = (v >> DRS::ktdc_ch_id_shift) & DRS::ktdc_ch_id_mask;
  return decodeTdcWord(v, ct, ft);
}

//______________________________________________________________________________
unsigned int
getTDC(unsigned int traw)
{
  unsigned int ret = traw & 0x00ffffff;
  return ret;

  // following calculations are performed in FPGA
  unsigned int c1 = (ret>>11);
  c1 = c1 & 0x1fff;
  unsigned int c2 = (ret>>10);
  c2 = c2 & 0x1;
  c2 = 2 -c2;
  unsigned int ft = (ret & 0x3ff);
  ret = (c1<<11) - (c2<<10) - ft;
  return ret & 0x00ffffff;

}

//______________________________________________________________________________
unsigned int
calc_dt(unsigned int traw0,
        unsigned int traw1)
{
  unsigned int t0 = getTDC(traw0);
  unsigned int t1 = getTDC(traw1);

  return (t0 - t1) & 0x00ffffff;
}


// Open ----------------------------------------------------------------------
bool
Decoder::Open(const char* name)
{
  static const std::string MyFunc = "::Open ";

  {// local scope
    std::ifstream ifs_org(name, std::ios::binary);

    if(!ifs_org.is_open()){
      std::cerr << MyName << MyFunc
                << "File open error " << name << std::endl;
      return false;
    }
  }

  // add gzip compressor to filter chain
 // m_ifs.push(boost::iostreams::gzip_decompressor());
  m_ifs.push(boost::iostreams::file_source(name, std::ios::binary));
  m_ifs.auto_close();
  
  // Initialize process for DRS4
  for(int imod=0;imod<DRS::NofBoards;imod++){
    flag[imod].reset();
    flag[imod].set(Good);
    flag[imod].set(FirstData);
  }

  return true;
}

// getNextEvent --------------------------------------------------------------
bool
Decoder::getNextEvent()
{
  static const std::string MyFunc = "::getNextEvent ";
  // EOF ?
  if(m_ifs.eof()){
    std::cout << "#D : " << MyName << MyFunc
	      << "End Of File" << std::endl;
    return false;
  }

  const unsigned int DAQMW_HEADERSIZE = 8 ; //bytes
  const unsigned int DAQMW_FOOTERSIZE = 8 ;
  char buf[DAQMW_HEADERSIZE];
  m_ifs.read(buf,DAQMW_HEADERSIZE);

  //DAQMW header check 
  if( ((buf[0] & 0x000000ff)  != 0xe7 ) && ((buf[1] & 0x000000ff)  != 0xe7 )){
    std::cerr << "Invalid DAQMW header !"  << std::endl;
    for(unsigned int i=0;i<DAQMW_HEADERSIZE; i++){
      std::cout << std::setfill('0') << std::hex << std::setw(8) << (int) buf[i] << " " ;
    }
    std::cerr << " " << std::endl;

    return false;
  }

  //read Comp header
  char bufComp[2*sizeof(uint32_t)];
  m_ifs.read(bufComp,2*sizeof(uint32_t));

  //GO TO NIM-EASIROC 
  //std::cout << "start NIM-EASIROC" << std::endl;

  struct FiberHit sfiber[FiberNch];

  for(int ich=0; ich<FiberNch; ich++){
    sfiber[ich].tdcleading.clear();
    sfiber[ich].tdctrailing.clear();
  }

  for(int imod = 0 ;imod<6;++imod){
    //comp header
    char NERcompheader[2*sizeof(uint32_t)];
    m_ifs.read(NERcompheader,2*sizeof(uint32_t));

    char headerByte[4];
    m_ifs.read(headerByte, 4);
    //for(int i =0 ;i<4;i++){
    //  std::cout  << (int)headerByte[i] << std::endl;
    //}

    unsigned int header32 = getBigEndian32(headerByte);
    unsigned int header = Decode32bitWord(header32);
    bool isHeader = ((header >> 27) & 0x01) == 0x01;
    if(!isHeader) {
      //std::cerr << imod << std::endl;
      //std::cerr << "Frame Error of header data" << std::endl;
      std::cerr << "EOF " << std::endl;
      fprintf(stderr, "    %08X\n", header);
      //for(unsigned int i=0;i<DAQMW_HEADERSIZE; i++){
      //  std::cout << std::setfill('0') << std::hex << std::setw(8) << (int) buf[i] << " " ;
      //}
      //std::cout << std::dec << std::endl;
      //for(int i =0 ;i<4;i++){
      //  std::cout << std::setfill('0') << std::hex << std::setw(8) << (int) headerByte[i] << " " ;
      //}
      //std::cout << std::dec << std::endl;
      //std::exit(1);
      return false;
    }
    size_t dataSize = header & 0x0fff;
    //std::cout << "datasize " << dataSize << std::endl;
    unsigned int scalerValues[69];
    char* dataBytes = new char[dataSize * 4];
    m_ifs.read(dataBytes, dataSize * 4);  

    //std::cout << "Datasize=" << dataSize << std::endl;
    for(size_t i = 0; i < dataSize; ++i) {
      unsigned int data32 = getBigEndian32(dataBytes + 4 * i);
      unsigned int data = Decode32bitWord(data32); 
      if(!data){
        std::cout << "Invalid data event " <<  std::endl;
        std::cout << "data size " << dataSize << std::endl;
        break;
      }
      //int ch_test = (data >> 13) & 0x3f;
      //cout << "event " << events << " ch " << ch_test << endl;
      int ch0 = (data >> 13) & 0x3f;
      int ch = ch0 + 64*imod; 

      //std::cout << ch0 << ":" << ch << std::endl;

      if(isAdcHg(data)) {
        bool otr = ((data >> 12) & 0x01) != 0;
        int value = data & 0x0fff;
        sfiber[ch].otradchigh = otr;
        sfiber[ch].adchigh = value;
	// std::cout<< "***" << std::endl;
	// std::cout<< "ch= " << ch 
	// 	 << " Adc=" << value
	// 	 << std::endl;
        if(!otr) {
          //hit ?
          //if(value > thre[ch]){
          //cout << "HIT layer " << layer << "type" << type << " ch " << ch << " fiber "<< fiber << endl;
          //hitprofile[layer][type]->Fill(pos);
          //hitprofile_ch[layer][type]->Fill(fiber);
          // if(layer==0){
          //   hitpos[1][0][type][multi[layer][type]]=pos;//pat. 1 L0-L1
          //   hitpos[2][0][type][multi[layer][type]]=pos;//pat. 2 L0-L2
          //   hitpos[3][0][type][multi[layer][type]]=pos;//pat. 3 L0-L3
          //   fiberarr[1][0][type][multi[layer][type]]=fiber;//pat. 1 L0-L1
          //   fiberarr[2][0][type][multi[layer][type]]=fiber;//pat. 2 L0-L2
          //  fiberarr[3][0][type][multi[layer][type]]=fiber;//pat. 3 L0-L3
          // }else{
          //   hitpos[layer][1][type][multi[layer][type]]=pos;
          //   fiberarr[layer][1][type][multi[layer][type]]=fiber;
          // }
          // multi[layer][type]++;
          // if(multi[layer][type] > maxbuf) cout << "too many events" << endl;
          // }
          //adcHigh[ch]->Fill(value);
        }
      }else if(isAdcLg(data)) {
        bool otr = ((data >> 12) & 0x01) != 0;
        int value = data & 0x0fff;
        sfiber[ch].otradclow = otr;
        sfiber[ch].adclow = value;
	// std::cout<< "***" << std::endl;
	// std::cout<< "ch= " << ch 
	// 	 << " Adc=" << value
	// 	 << std::endl;
        if(!otr) {
          //adcLow[ch]->Fill(value);
        }
      }else if(isTdcLeading(data)) {
        int value = data & 0x0fff;
        sfiber[ch].tdcleading.push_back(value);
        //sfiber[ch].tdcleading = value;
        //tdcLeading[ch]->Fill(value);
	// std::cout<< "***" << std::endl;
	// std::cout<< "ch= " << ch 
	// 	 << " lTdc=" << value
	// 	 << std::endl;
      }else if(isTdcTrailing(data)) {
        int value = data & 0x0fff;
        sfiber[ch].tdctrailing.push_back(value);
        //sfiber[ch].tdctrailing = value;
        //tdcTrailing[ch]->Fill(value);
	// std::cout<< "***" << std::endl;
	// std::cout<< "ch= " << ch 
	// 	 << " tTdc=" << value
	// 	 << std::endl;
      }else if(isScaler(data)) {
        int value = data & 0x3fff;
        scalerValues[ch] = value;
        //cout << "event:"<<events<<"/scalerValues["<<ch<<"]:"<<scalerValues[ch] << endl; 
        if(ch == 68) {
          //int scalerValuesArrayIndex = events % 100;
          //memcpy(scalerValuesArray[scalerValuesArrayIndex], scalerValues,
          //    sizeof(scalerValues));
        }
      }else {
        int ch = (data >> 13) & 0x3f;
        int value = data & 0x0fff;
        std::cout << "adchg:"  << (data & 0x00680000);
        std::cout << "adclg:"  << (data & 0x00680000);
        std::cout << "tdcl:"   << (data & 0x00601000);
        std::cout << "tdct:"   << (data & 0x00601000);
        std::cout << "scaler:" << (data & 0x00600000);
        std::cout << "data:" << data << std::endl; 
        std::cout << "ch:" << ch << " value:" << value << std::endl;
        std::cerr << "Unknown data type" << std::endl;
      }

      //store ADC vs ToT correlation
      // if((sfiber[ch].tdctrailing != -1) && (sfiber[ch].tdcleading !=-1)){
      //double tot =  sfiber[ch].tdcleading - sfiber[ch].tdctrailing;
      //tdcToT[ch]->Fill(tot);
      //adcToTcorr[ch]->Fill(sfiber[ch].adchigh,tot);
      // }
    }//end of loop of dataSize 
    delete[] dataBytes;
  }//end of NIM EASIROC

  fibertdcl_.resize(FiberNch);
  fibertdct_.resize(FiberNch);
  for(int ich=0;ich<FiberNch;ich++){
    fiberadch_[ich] = sfiber[ich].adchigh;
    fiberadcl_[ich] = sfiber[ich].adclow;

    int sizeltdc = sfiber[ich].tdcleading.size();
    for(int ilt=0;ilt<sizeltdc;ilt++){
      fibertdcl_[ich].push_back(sfiber[ich].tdcleading[ilt]);
    }
    int sizettdc = sfiber[ich].tdctrailing.size();
    for(int itt=0;itt<sizettdc;itt++){
      fibertdct_[ich].push_back(sfiber[ich].tdctrailing[itt]);
    }
  }

  //std::cout<< "*********************" << std::endl;

  // //start DRS4 data decoding
  // for(int imod=0;imod<DRS::NofBoards;imod++){
  //   // read component header
  //   bufType compHeader;
  //   read(compHeader, 2);
  //   //std::cout << std::hex << compHeader[0] << " " << compHeader[1] << std::dec << std::endl;

  //   // read base header
  //   if (compHeader[0]!=0x504d4f43) {
  //     head_buffer_[imod].resize(2);
  //     head_buffer_[imod][0] = compHeader[0];
  //     head_buffer_[imod][1] = compHeader[1];
  //     compHeader.resize(sizeBaseHeader-2);
  //     read(compHeader, compHeader.size());
  //     std::copy(compHeader.begin(), compHeader.end(), std::back_inserter(head_buffer_[imod]));
  //   } else {
  //     read(head_buffer_[imod], (int)sizeBaseHeader);
  //   }

  //   //std::cout << "#D Base header " << std::endl;
  //   //std::for_each(head_buffer_[imod].begin(), head_buffer_[imod].end(), hddaq::HexDump());
  //   //std::cout<< "DRS4 Event Num = " << std::dec << ((head_buffer_[imod].at(2))&0x00ff) << std::endl;

  //   if(flag[imod][EoF]){
  //     std::cout <<"L." << __LINE__ << " EOF return" << std::endl;
  //     return false;
  //   }

  //   int n_word_body   = head_buffer_[imod][i_datasize] & DRS::data_size_mask;
  //   if(flag[imod][FirstData]){
  //     int n_word_header = (head_buffer_[imod][i_datasize] >> DRS::header_size_shift) & DRS::header_size_mask;
  //     if(n_word_header == sizeBaseHeader + sizeExHeader){
  //       flag[imod].set(ExHead);
  //     }// if(n_work)

  //     flag[imod].reset(FirstData);
  //   }//if(flag)

  //   // read drs header
  //   if(flag[imod][ExHead]) read(exhead_buffer_[imod], (int)sizeExHeader);

  //   // read body
  //   read(body_buffer_[imod], n_word_body);

  //   if (has_tdc_data)
  //   { // TDC data read header
  //     tdc_buffer_[imod].resize(DRS::sizeTdcHeader, 0);
  //     read(tdc_buffer_[imod], static_cast<int>(tdc_buffer_[imod].size()));
  //     // std::cout << "#D TDC header " << std::endl;
  //     // std::for_each(tdc_buffer_[imod].begin(), tdc_buffer_[imod].end(), hddaq::HexDump());

  //     for (;;) {
  //       //unsigned int d32 = tdc_header.back();
  //       unsigned int d32 = tdc_buffer_[imod][DRS::i_frame_size];
  //       unsigned int data_type = (d32 >> DRS::ktdc_data_type_shift) & DRS::ktdc_data_type_mask;
  //       int seq_id       = (d32 >> DRS::ktdc_num_frame_shift) & DRS::ktdc_num_frame_mask;
  //       int tdc_nword    = (d32 >> DRS::ktdc_frame_size_shift) & DRS::ktdc_frame_size_mask;
  //       // std::cout << "#D TDC num frae = " << seq_id
  //       //           << ", n word = " << tdc_nword
  //       //           << " data type = " << std::hex << data_type << std::dec
  //       //           << std::endl;

  //       bool isLastFrame =  (data_type == DRS::ktdc_type_last_frame_size);
  //       //        std::cout << "#D isLastFrame = " << isLastFrame << std::endl;

  //       if (seq_id==0) tdc_nword -= DRS::sizeTdcHeader;
  //       if (!isLastFrame) ++tdc_nword;

  //       //std::cout << "#D n word (mod) = " << tdc_nword << std::endl;

  //       bufType tdcBufTmp(tdc_nword);

  //       read(tdcBufTmp, tdc_nword);
  //       //        std::for_each(tdcBufTmp.begin(), tdcBufTmp.end(), hddaq::HexDump());
  //       tdc_buffer_[imod].insert(tdc_buffer_[imod].end(), tdcBufTmp.begin(), tdcBufTmp.end());
  //       //std::for_each(tdc_buffer_.begin(), tdc_buffer_.end(), hddaq::HexDump());
  //       if (isLastFrame) break;
  //     }
  //   }
  // }

  //MznHRTdc
  char mzncompheader[2*sizeof(uint32_t)];
  m_ifs.read(mzncompheader,2*sizeof(uint32_t));

  unsigned int mznheader1=0;
  unsigned int mznheader2=0;
  unsigned int mznheader3=0;
  int          mznnword=0;
  int          index_slot=0;

  sMZNData.clear();
  
  m_ifs.read((char*)&mznheader1, sizeof(unsigned int));
  if(mznheader1 !=MZN_HRTDC::magic) {
    std::cerr << "Not MZN_HRTDC magic word " << std::hex << std::setw(8) << mznheader1 << std::endl;
    return false;  
  }
  
  m_ifs.read((char*)&mznheader2, sizeof(unsigned int));
  mznnword = mznheader2 & MZN_HRTDC::mask_nword;
  m_ifs.read((char*)&mznheader3, sizeof(unsigned int));
  int eventNum = mznheader3 & 0xffff;
  
  //std::cout<< "MZN Event Num = " << eventNum << std::endl;
  sMZNData.EventNum = eventNum;

  for(int i=0;i<mznnword;i++){
    unsigned int rdata =0;
    m_ifs.read((char*)&rdata, sizeof(unsigned int));
    unsigned int reg_slot = rdata & MZN_HRTDC::mask_slot;
  
    if(MZN_HRTDC::magic_slot_u == reg_slot || MZN_HRTDC::magic_slot_d == reg_slot){
      // slot identification
      index_slot = MZN_HRTDC::magic_slot_u == reg_slot ? 0 : 1;
      sMZNData.OverFlow[index_slot] = (rdata >> MZN_HRTDC::shift_overflow) & 0x1;
      sMZNData.StopOut[index_slot]  = (rdata >> MZN_HRTDC::shift_stop_dout) & 0x1;
      sMZNData.Through[index_slot]  = (rdata >> MZN_HRTDC::shift_through) & 0x1;

      //std::cout << "**** " << index_slot << " " << reg_slot << std::endl;
    }
    // Data body
    unsigned int data_header = (rdata >> MZN_HRTDC::shift_data_head) & MZN_HRTDC::mask_data_head;
    // std::cout << "#D data type = " << data_header << std::endl;
    if(data_header == MZN_HRTDC::magic_leading){
      sMZNData.tdc_leading[sMZNData.Nhit] = rdata & MZN_HRTDC::mask_tdc;

      //std::cout << "**** " << sMZNData.tdc_leading[index_slot] << std::endl;

      sMZNData.Ch[sMZNData.Nhit]          = ((rdata >> MZN_HRTDC::shift_ch)  & MZN_HRTDC::mask_ch) + index_slot*32;
      unsigned int fc               = (rdata) & MZN_HRTDC::mask_fc;
      sMZNData.FineCount[sMZNData.Nhit]   = sMZNData.Through[index_slot] == 1 ? fc : -1;
      sMZNData.Estimator[sMZNData.Nhit]   = sMZNData.Through[index_slot] == 0 ? fc : -1;
      sMZNData.CoarseCount[sMZNData.Nhit] = (rdata >> MZN_HRTDC::shift_cc)  &  MZN_HRTDC::mask_cc;

      sMZNData.NhitCh[sMZNData.Ch[sMZNData.Nhit]]++;
      sMZNData.Nhit++;

    }// leading data
    
    if(data_header ==  MZN_HRTDC::magic_trailing){
      sMZNData.tdc_trailing[sMZNData.tNhit] = rdata & MZN_HRTDC::mask_tdc;

      //std::cout << "**** " << sMZNData.tdc_trailing[index_slot] << std::endl;

      sMZNData.tCh[sMZNData.tNhit]          = ((rdata >> MZN_HRTDC::shift_ch)  & MZN_HRTDC::mask_ch) + index_slot*32;
      unsigned int fc                 = (rdata)              & MZN_HRTDC::mask_fc;
      sMZNData.tFineCount[sMZNData.tNhit]   = sMZNData.Through[index_slot] == 1 ? fc : -1;
      sMZNData.tEstimator[sMZNData.tNhit]   = sMZNData.Through[index_slot] == 0 ? fc : -1;
      sMZNData.tCoarseCount[sMZNData.tNhit] = (rdata >> MZN_HRTDC::shift_cc)  & MZN_HRTDC::mask_cc;

      sMZNData.tNhitCh[sMZNData.tCh[sMZNData.tNhit]]++;
      sMZNData.tNhit++;
    }// trailing data

    if(data_header == MZN_HRTDC::magic_stop){
      sMZNData.common_stop     = rdata & MZN_HRTDC::mask_tdc;

      unsigned int fc       = (rdata) & MZN_HRTDC::mask_fc;
      sMZNData.StopFineCount   = sMZNData.Through[index_slot] == 1 ? fc : -1;
      sMZNData.StopEstimator   = sMZNData.Through[index_slot] == 0 ? fc : -1;
      sMZNData.StopCoarseCount = (rdata >> MZN_HRTDC::shift_cc)  & MZN_HRTDC::mask_cc;
    }// common stop data
  }//for mznnword

  char DAQMWfooter[DAQMW_FOOTERSIZE];
  m_ifs.read(DAQMWfooter,DAQMW_FOOTERSIZE);
  //DAQMW footer check 
  if( ((DAQMWfooter[0] & 0x000000ff)  != 0xcc ) && ((DAQMWfooter[1] & 0x000000ff)  != 0xcc )){
    std::cerr << "Invalid DAQMW footer !"  << std::endl;
    for(unsigned int i=0;i<DAQMW_FOOTERSIZE; i++){
      std::cout << std::setfill('0') << std::hex << std::setw(8) << (int) DAQMWfooter[i] << " " ;
      std::cout << std::dec  << std::endl;
    }
  }

  return true;
}

void Decoder::GetFiberData(int adch[FiberNch], int adcl[FiberNch], 
			   std::vector<std::vector<int> >& tdcl, 
			   std::vector<std::vector<int> >& tdct)
{
  for(int ich=0; ich<FiberNch; ich++){
    tdcl[ich].clear();
    tdct[ich].clear();
  }

  for(int ich=0;ich<FiberNch;ich++){
    adch[ich] = fiberadch_[ich];
    adcl[ich] = fiberadcl_[ich];
  }

  tdcl.resize(FiberNch);
  tdct.resize(FiberNch);
  for(int ich=0;ich<FiberNch;ich++){
    int sizeltdc = fibertdcl_[ich].size();
    for(int ilt=0;ilt<sizeltdc;ilt++){
      tdcl[ich].push_back(fibertdcl_[ich][ilt]);
    }
    int sizettdc = fibertdct_[ich].size();
    for(int itt=0;itt<sizettdc;itt++){
      tdct[ich].push_back(fibertdct_[ich][itt]);
    }
  }

  for(int ich=0; ich<FiberNch; ich++){
    fibertdcl_[ich].clear();
    fibertdct_[ich].clear();
  }
}

// decode --------------------------------------------------------------------
bool
Decoder::decode(DRS::dataDrs& cont)
{
  static const std::string MyFunc = "::decode ";

  // clear vector
  for(int i = 0; i<DRS::NofChModule*DRS::NofBoards; ++i){
    cont.data_wf[i].clear();
  }
  //std::cout << "clear" << std::endl;
  // decode base header
  for(int imod=0;imod<DRS::NofBoards;imod++){
    cont.nword_header[imod] = (head_buffer_[imod][i_datasize] >> DRS::header_size_shift) &DRS::header_size_mask;
    cont.nword_body[imod] =  head_buffer_[imod][i_datasize] & DRS::data_size_mask;
    cont.global_tag[imod]    = (head_buffer_[imod][i_tag] >> DRS::gtag_shift) & DRS::gtag_mask;
    cont.local_tag[imod]      =  head_buffer_[imod][i_tag] & DRS::ltag_mask;
    cont.tic_count[imod]    =  head_buffer_[imod][i_tic] & DRS::tic_mask;
    cont.fl_double_data[imod] = (head_buffer_[imod][i_tic] >> DRS::datastr_shift ) & DRS::datastr_mask;
    // std::cout << "#D " << std::hex << cont.nword_header
    //           << " " << cont.nword_body
    //           << " " << cont.local_tag
    //           << std::endl;

    //std::cout<< "DRS4 Event Num = " << std::dec << (head_buffer_[imod].at(2)&0x00ff)  << std::endl;
    cont.eventNum = (head_buffer_[imod].at(2)&0xffff);

    // decode ex header
    if(flag[imod][ExHead]){
      //std::cout << "ex" << std::endl;
      for(unsigned int i = 0; i<sizeExHeader; ++i){
        int icel = i + imod*DRS::NofDrs;
        cont.wsr[icel]     = (exhead_buffer_[imod][icel] >> DRS::wsr_shift) & DRS::wsr_mask;
        cont.cellnum[icel] =  exhead_buffer_[imod][icel] & DRS::cellnum_mask;
      }// for(i)
    }// if(flag)

    // decode body
    for(int i = 0; i<cont.nword_body[imod]; ++i){
      //std::cout << "body" << std::endl;
      unsigned int data_type = (body_buffer_[imod][i] >> DRS::type_shift) & DRS::type_mask;
      unsigned int ch        = (body_buffer_[imod][i] >> DRS::ch_shift) & DRS::ch_mask;
      int ich = ch + imod*DRS::NofChModule;
      if(data_type == DRS::type_wf){
        if(cont.fl_double_data){
          double data_high = (double)((body_buffer_[imod][i] >> DRS::wf_high_shift) & DRS::wf_mask_12bit);
          double data_low  = (double)(body_buffer_[imod][i] & DRS::wf_mask_12bit);
          cont.data_wf[ich].push_back(data_low);
          cont.data_wf[ich].push_back(data_high);
        }else{
          unsigned int data  = (body_buffer_[imod][i] & DRS::wf_mask_14bit);
          unsigned int extra_bits  = data & 0x3; // extra 2bits
          double wf_data = (double)(data >> 2) + extra_bits*0.25;
          cont.data_wf[ich].push_back(wf_data);
          //	if((ch == 0 || ch == 15 ) && cont.data_wf[ch].size()>511) std::cout << "ch:" << ch << " " << cont.data_wf[ch].size() << std::endl;
        }//Peak serch
      }else if(data_type == DRS::type_qdc){
        unsigned int data_qdc = body_buffer_[imod][i] & DRS::qdc_mask;
        cont.data_qdc[ich] = data_qdc;
      }else{
        std::cout << "#E : " << MyName << MyFunc
          << "Unknown data type " << std::hex << data_type 
          << std::dec << std::endl;
        std::cout << "mod id " << imod << std::endl;
        std::cout << "D ith = " << i
          << " nword header = " << cont.nword_header[imod] 
          << "\n nword body = " << cont.nword_body[imod]
          << "\n global tag = " << cont.global_tag[imod]
          << "\n local tag = "  << cont.local_tag[imod]
          << "\n tic count = " << cont.tic_count[imod]
          << "\n fl_double = " << cont.fl_double_data[imod]
          << std::endl;
        std::for_each(body_buffer_[imod].begin(),   body_buffer_[imod].end(),   hddaq::HexDump());
        return false;
      }// if(data_type)
    }// for(i)
  }//for (imod)

  return true;
}

//______________________________________________________________________________
bool
Decoder::decodeTDC(DRS::dataDrs& cont)
{
  if (!has_tdc_data) return true;
  
  for (unsigned int i=0; i<cont.data_tdc.size(); ++i) {
    cont.data_tdc[i].clear();
    cont.data_dt[i].clear();

    cont.data_tdc_2nd[i].clear();
    cont.data_dt_2nd[i].clear();

    cont.data_width[i].clear();
  }
  
  for(int imod = 0 ; imod <DRS::NofBoards ;imod++){
    int n = tdc_buffer_[imod].size();
    unsigned int l1_traw  = tdc_buffer_[imod][DRS::i_l1_tdc];
    unsigned int l1_t1raw = tdc_buffer_[imod][DRS::i_l1_tdc_2nd];
    {
      unsigned int ct;
      unsigned int ft;
      decodeTdcWord(l1_traw, ct, ft);
      cont.l1_tdc[imod] = getTDC(l1_traw);
    }
    {
      unsigned int ct;
      unsigned int ft;
      decodeTdcWord(l1_t1raw, ct, ft);
      cont.l1_tdc1[imod] = getTDC(l1_t1raw);
    }

    for (int i=DRS::sizeTdcHeader; i<n; ++i) {
      //std::cout << "#D i = " << i  << std::endl;
      unsigned int v = tdc_buffer_[imod][i];
      unsigned int ch;
      unsigned int ct;
      unsigned int ft;
      unsigned int data_type = decodeTdcWord(v, ch, ct, ft);
      unsigned int tdc = getTDC(v);
      unsigned int dt = calc_dt(l1_traw, v);
      unsigned int ich = ch + imod*DRS::NofChModule;
      if (data_type == DRS::ktdc_type_leading_edge){
        // std::cout << "data_type = " << std::hex << data_type << std::dec
        //           << ", ch = " << ch
        //           << ", ct = " << ct
        //           << ", ft = " << ft
        //           << std::endl;
        cont.data_tdc[ich].push_back(tdc*kTDCTimeUnit);
        cont.data_dt[ich].push_back(dt*kTDCTimeUnit);
      } else if ((data_type == DRS::ktdc_type_trailing_edge) || (data_type == DRS::ktdc_type_2nd_LE)){
        cont.data_tdc_2nd[ich].push_back(tdc*kTDCTimeUnit);
        cont.data_dt_2nd[ich].push_back(dt*kTDCTimeUnit);

	if(cont.data_tdc[ich].size()>0){
	  unsigned int width = (tdc - cont.data_tdc[ich].back())& 0xffffff;
	  cont.data_width[ich].push_back(width*kTDCTimeUnit);
	}
	else{
	  unsigned int width = 0;
	  cont.data_width[ich].push_back(width*kTDCTimeUnit);
	}
      }
    }
  }

  return true;

}

//______________________________________________________________________________
bool
Decoder::decodeADC(DRS::dataDrs& cont)
{
  ConfMan *confMan = ConfMan::GetConfManager();
  const int fIStart = confMan->PeakStart();
  const int fIEnd   = confMan->PeakEnd();
  const int BaseParam = confMan->BaseParam();
  const int RangeLow  = confMan->RangeLow();
  const int RangeHigh = confMan->RangeHigh();
  const int PileLow  = confMan->PileLow();
  const int PileHigh = confMan->PileHigh();

  for (unsigned int i=0; i<cont.data_adc.size(); ++i) {
    cont.data_adc[i].clear();
    cont.data_bl[i].clear();
    cont.data_amp[i].clear();
    cont.data_peakx[i].clear();
  }

  for(int ich = 0; ich<DRS::NofChModule*DRS::NofBoards; ++ich){
    int fPeakX  = -9999;
    int fPeakY  = -9999;
    double fBaseLine  = 0.0;
    double fAmplitude = 0.0;
    double fIntegral  = 0.0;
    double fWFx = 0.0;
    double fWFy = 0.0;
    double ph = 0.0;
    double x = 0.0;
    double max=-9999.;
    double xt = 0.0;

    //Serch peak 
    for (int j=fIStart; j<fIEnd; ++j){
      if (fPeakY < cont.data_wf[ich][j]){
	fPeakY = cont.data_wf[ich][j];
	fPeakX = j;
      }
    }

    //Calc Bese Line  
    int range1 = 1;
    int range2 = fPeakX - BaseParam;
    double sum=0.0;
    int    num=0;
    for (int j=range1; j<range2; ++j) {
      if (TMath::Abs(cont.data_wf[ich][j+1]-cont.data_wf[ich][j]) < 10) {
	++num;
	sum += cont.data_wf[ich][j];
      }
    }
    fBaseLine = sum/num;
    fAmplitude = fPeakY - fBaseLine;
	  
    //waveform
    const std::vector<double> wf = cont.data_wf[ich];
    int ncell = wf.size();
    
    for (int icell=fPeakX-PileLow; icell<fPeakX-PileHigh; ++icell){
      fWFx = icell;
      fWFy = wf[icell]-fBaseLine;
      double WFy = fWFy * kAmplitudeUnit;
      double WFx = fWFx * kTimeUnit;
      cont.data_bl[ich].push_back(WFy);
      cont.data_peakx[ich].push_back(WFx);
      if(max<WFy){
	max=WFy;
	x=WFx;
      }
      //std::cout << max << "   "<< x <<std::endl;
    }

    //Integral
    int n = cont.data_wf[ich].size();
    int imin = fPeakX - RangeLow;
    if (imin<0) imin = 0;
    int imax = fPeakX + RangeHigh;
    if (imax>n) imax = n;
    
    
    for (int j=imin; j<imax; ++j) {
      fIntegral += cont.data_wf[ich][j]-fBaseLine;
    }
    
    //Fill QDC information
    double Integral = fIntegral * kChargeUnit;
    double BaseLine = fBaseLine * kAmplitudeUnit;
    double Amplitude = fAmplitude * kAmplitudeUnit;
    double PeakX = fPeakX * kTimeUnit;
    
    cont.data_adc[ich].push_back(Integral);
    cont.data_bl[ich].push_back(BaseLine);
    cont.data_amp[ich].push_back(Amplitude);
    cont.data_peakx[ich].push_back(PeakX);
  }

  return true;
}

// read ----------------------------------------------------------------------
int 
Decoder::read(bufType& buf, int nword)
{
  static const std::string MyFunc = "::read ";

  buf.clear();
  buf.resize(nword);

  m_ifs.read((char*)&buf[0], nword*sizeof(unsigned int));
  int n_read_word = m_ifs.gcount()/(int)sizeof(unsigned int);
  if(n_read_word != nword){
    // std::cerr << "#E : " << MyName << MyFunc
    //           << "End of file" << std::endl;
    //std::cout << "#D hogee " << n_read_word << " " << nword << std::endl;
    for(int imod=0;imod<DRS::NofBoards;imod++){
      flag[imod].set(EoF);
    }
  }
  
  return n_read_word;  
}

// Constructor ---------------------------------------------------------------
Decoder::Decoder()
{
  for(int imod=0;imod<DRS::NofBoards;imod++){
    head_buffer_[imod].resize(sizeBaseHeader);
    exhead_buffer_[imod].resize(sizeExHeader);
    body_buffer_[imod].resize(1);
  }
}

// Destructor ----------------------------------------------------------------
Decoder::~Decoder()
{
  if (!m_ifs.empty()) m_ifs.reset();
}

// ---------------------------------------------------------------------------
// base line restorer
// ---------------------------------------------------------------------------
// resetTicCont --------------------------------------------------------------
void
Decoder::resetTicCont()
{
  for(int imod=0;imod<DRS::NofBoards;imod++){
    for(int i = 0; i<DRS::NofDrs; ++i){
      contTic[i][imod].clear();
      contTic[i][imod].resize(DRS::NofCell);

      for(int cell = 0; cell<DRS::NofCell; ++cell){
        contTic[i][imod][cell].fl_read       = false;
        contTic[i][imod][cell].time_interval = 0;
        contTic[i][imod][cell].fl_overflow   = true;
      }// for(cell)
    }//for NofDrs
  }// for(imod)

  return;
}// resetTicCont

// setThisEvent --------------------------------------------------------------
void
Decoder::setThisEvent(int tic, bool of)
{
  for(int imod=0;imod<DRS::NofBoards;imod++){
    for(int i = 0; i<DRS::NofDrs; ++i){
      for(int cell = 0; cell<DRS::NofCell; ++cell){
        if(of){
          contTic[i][imod][cell].time_interval = DRS::max_tic;
          contTic[i][imod][cell].fl_overflow   = true;	
        }else{
          if(!contTic[i][imod][cell].fl_overflow){
            contTic[i][imod][cell].time_interval += tic;
            if(contTic[i][imod][cell].time_interval > DRS::max_tic){
              contTic[i][imod][cell].time_interval = DRS::max_tic;
              contTic[i][imod][cell].fl_overflow   = true;
            }
          }// if(fl_overflow)
        }// if(of)
      }// for(cell)
    }// for(i)
  }

  return;
}// setThisEvent

//NIM-EASIROC decoding function
unsigned int Decoder::getBigEndian32(const char* b=NULL)
{
    //std::cout << "size of b " << sizeof(b) << std::endl;
    return ((b[0] << 24) & 0xff000000) |
           ((b[1] << 16) & 0x00ff0000) |
           ((b[2] <<  8) & 0x0000ff00) |
           ((b[3] <<  0) & 0x000000ff);
}

unsigned int Decoder::Decode32bitWord(unsigned int word32bit=0)
{
  //check data format
  unsigned int frame = word32bit & 0x80808080;
  if(frame != 0x80000000){
    //std::cerr << __FILE__ << " L." << __LINE__ << " Frame Error! " << std::endl;
    std::cerr << "32 bit word: " << std::hex << word32bit << std::dec << std::endl;
    return 0;
  }

  return ((word32bit & 0x7f000000) >> 3) | 
         ((word32bit & 0x007f0000) >> 2) |
         ((word32bit & 0x00007f00) >> 1) |
         ((word32bit & 0x0000007f) >> 0);
}


//ADC High Gain
bool Decoder::isAdcHg(unsigned int data=0)
{
    return (data & 0x00680000) == 0x00000000;
}

//ADC Low Gain
bool Decoder::isAdcLg(unsigned int data=0)
{
    return (data & 0x00680000) == 0x00080000;
}

//TDC Leading
bool Decoder::isTdcLeading(unsigned int data=0)
{
    return (data & 0x00601000) == 0x00201000;
}

//TDC Trailing
bool Decoder::isTdcTrailing(unsigned int data=0)
{
    return (data & 0x00601000) == 0x00200000;
}

//not impletented yet 
bool Decoder::isScaler(unsigned int data=0)
{
    return (data & 0x00600000) == 0x00400000;
}


