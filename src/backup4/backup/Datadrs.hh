#ifndef DATADRS_
#define DATADRS_

#include<vector>

namespace DRS{ 
  // constants for TDC
  enum IndexTdcHeader{
    i_tdc_magic, i_tdc_l1_cnt, i_tdc_l2_cnt, i_l1_tdc, i_l1_tdc_2nd,
    i_frame_size, sizeTdcHeader
  };

  static const unsigned int ktdc_data_type_shift  = 28;
  static const unsigned int ktdc_data_type_mask   = 0xf;
  enum ETdcDataType {
    ktdc_type_leading_edge    = 0x0,
    ktdc_type_trailing_edge   = 0x1,
    ktdc_type_2nd_LE          = 0x2,
    ktdc_type_magic           = 0x4,
    ktdc_type_l1_counter      = 0x5,
    ktdc_type_l2_counter      = 0x6,
    ktdc_type_l1_tdc          = 0x7,
    ktdc_type_frame_size      = 0xc,
    ktdc_type_last_frame_size = 0xd,
    ktdc_type_error           = 0x9,
    ktdc_type_num_hit         = 0xa,
    ktdc_type_trailer         = 0xb
  };
  static const unsigned int ktdc_num_frame_shift   = 16;
  static const unsigned int ktdc_num_frame_mask    = 0xfff;
  static const unsigned int ktdc_frame_size_shift  = 0;
  static const unsigned int ktdc_frame_size_mask   = 0xfff;
  static const unsigned int ktdc_ch_id_shift       = 24;
  static const unsigned int ktdc_ch_id_mask        = 0xf;
  static const unsigned int ktdc_coarse_time_shift = 10;
  static const unsigned int ktdc_coarse_time_mask  = 0x3fff; // 14 bits
  static const unsigned int ktdc_fine_time_shift   = 0;
  static const unsigned int ktdc_fine_time_mask    = 0x3ff; // 10 bits
  static const int NofCell     = 1024;
  static const int NofChModule = 16;
  static const int NofChDrs    = 16;
  static const int NofChInDrs  = 8;
  static const int NofDrs      = 4;
  static const int NofBoards   = 1;
  static const unsigned int data_size_mask    = 0xffff;
  static const unsigned int header_size_mask  = 0xf;
  static const unsigned int header_size_shift = 16;
  
  static const unsigned int gtag_shift        = 16;
  static const unsigned int gtag_mask         = 0xf;
  static const unsigned int ltag_mask         = 0xffff;

  static const unsigned int tic_mask          = 0xfffff;
  static const unsigned int datastr_mask      = 0x1;
  static const unsigned int datastr_shift     = 20;

  static const unsigned int wsr_shift         = 16;
  static const unsigned int wsr_mask          = 0xff;
  static const unsigned int cellnum_mask      = 0x3ff;

  static const unsigned int ch_mask           = 0xf;
  static const unsigned int ch_shift          = 24;

  static const unsigned int type_shift        = 28;
  static const unsigned int type_mask         = 0xf;
  static const unsigned int type_header       = 0xf;
  static const unsigned int type_wf           = 0xd;
  static const unsigned int type_qdc          = 0x9;

  // double data structure
  static const unsigned int wf_high_shift    = 12;
  static const unsigned int wf_mask_12bit    = 0xfff;
  // sigle data structure
  static const unsigned int wf_mask_14bit    = 0x3fff;

  static const unsigned int qdc_mask          = 0xffff;
  static const int non_cascade = 1;
  static const int cascade     = 2;

  static const int max_tic     = 13000;

  struct dataDrs{
    int nword_header[NofBoards];
    int nword_body[NofBoards];
    int global_tag[NofBoards];
    int local_tag[NofBoards];
    int tic_count[NofBoards];
    int fl_double_data[NofBoards]; // 1:12bit wf data x2, 0:14bit wf data x1

    int wsr[NofDrs*NofBoards];
    int cellnum[NofDrs*NofBoards];

    int eventNum;

    //module ID
    int modid;

    std::vector<std::vector<double> > data_wf;
    unsigned int data_qdc[NofChModule*NofBoards];

    std::vector<std::vector<double> > data_adc;
    std::vector<std::vector<double> > data_bl;
    std::vector<std::vector<double> > data_amp;
    std::vector<std::vector<double> > data_peakx;

    std::vector<std::vector<unsigned int> > data_tdc;
    std::vector<std::vector<unsigned int> > data_dt;
    std::vector<std::vector<unsigned int> > data_tdc_2nd;
    std::vector<std::vector<unsigned int> > data_dt_2nd;
    std::vector<std::vector<unsigned int> > data_width;
    //unsigned int l1_tdc;
    //unsigned int l1_tdc1;
    unsigned int l1_tdc[NofBoards];
    unsigned int l1_tdc1[NofBoards];

    dataDrs()
      : data_wf(NofChModule*NofBoards),
      data_adc(NofChModule*NofBoards),
      data_bl(NofChModule*NofBoards),
      data_amp(NofChModule*NofBoards),
      data_peakx(NofChModule*NofBoards),
      data_tdc(NofChModule*NofBoards),
      data_dt(NofChModule*NofBoards),
      data_tdc_2nd(NofChModule*NofBoards),
      data_dt_2nd(NofChModule*NofBoards),
      data_width(NofChModule*NofBoards)
    {}
  };

};

#endif
