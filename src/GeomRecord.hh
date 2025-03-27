/*
  GeomRecord.hh

  2024/04  K.Shirotori
*/

#ifndef GeomRecord_h
#define GeomRecord_h 1

#include "ThreeVector.hh"

#include <string>
#include <functional>

class GeomRecord
{
public:
   GeomRecord( int id, const char *name,
               double x, double y, double z, double ta,
               double ra1, double ra2, double length, double resol,
               double w0, double dd, double ofs )
      : id_(id), name_(name), pos_(x,y,z), tiltAngle_(ta),
        rotAngle1_(ra1), rotAngle2_(ra2),
        length_(length), resol_(resol), w0_(w0), dd_(dd), ofs_(ofs)
      { calcVectors(); }
   
   GeomRecord( int id, const std::string &name,
               double x, double y, double z, double ta,
               double ra1, double ra2, double length, double resol,
               double w0, double dd, double ofs )
      : id_(id), name_(name), pos_(x,y,z), tiltAngle_(ta),
        rotAngle1_(ra1), rotAngle2_(ra2),
        length_(length), resol_(resol), w0_(w0), dd_(dd), ofs_(ofs)
      { calcVectors(); }
   
   GeomRecord( int id, const char *name,
               const ThreeVector pos, double ta,
               double ra1, double ra2, double length, double resol,
               double w0, double dd, double ofs )
      : id_(id), name_(name), pos_(pos),  tiltAngle_(ta),
      rotAngle1_(ra1), rotAngle2_(ra2),
        length_(length), resol_(resol), w0_(w0), dd_(dd), ofs_(ofs)
      { calcVectors(); }
   
   GeomRecord( int id, const std::string &name,
               const ThreeVector pos, double ta,
               double ra1, double ra2, double length, double resol,
               double w0, double dd, double ofs )
      : id_(id), name_(name), pos_(pos),  tiltAngle_(ta),
        rotAngle1_(ra1), rotAngle2_(ra2),
        length_(length), resol_(resol), w0_(w0), dd_(dd), ofs_(ofs)
      { calcVectors(); }
   
   ~GeomRecord() {}
   GeomRecord( const GeomRecord & );
   GeomRecord & operator=( const GeomRecord );

public:
   const ThreeVector & Position( void ) const { return pos_; }
   ThreeVector NormalVector( void ) const 
      { return ThreeVector( dxdu_, dydu_, dzdu_ ); }
   ThreeVector UnitVector( void ) const
      { return ThreeVector( dxds_, dyds_, dzds_ ); }
   
   double dsdx( void ) const { return dsdx_; }
   double dsdy( void ) const { return dsdy_; }
   double dsdz( void ) const { return dsdz_; }
   double dtdx( void ) const { return dtdx_; }
   double dtdy( void ) const { return dtdy_; }
   double dtdz( void ) const { return dtdz_; }
   double dudx( void ) const { return dudx_; }
   double dudy( void ) const { return dudy_; }
   double dudz( void ) const { return dudz_; }
   
   double dxds( void ) const { return dxds_; }
   double dxdt( void ) const { return dxdt_; }
   double dxdu( void ) const { return dxdu_; }
   double dyds( void ) const { return dyds_; }
   double dydt( void ) const { return dydt_; }
   double dydu( void ) const { return dydu_; }
   double dzds( void ) const { return dzds_; }
   double dzdt( void ) const { return dzdt_; }
   double dzdu( void ) const { return dzdu_; }

   double CalcWirePos_BFT( int wire, int layer ) const {
      //layer1だけ特別処理を行う（layer4も追加した）

      if(layer==301){ //layer1 の時
         int split_start=119;
         int split_end  =136;
         if(wire < split_start){
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if (wire % 2 == 1){ // wireが奇数の場合
               position = dd_long * ((double)wire-1.0) / 2.0 + dd_short * ((double)wire-1.0) / 2.0;
            }
            else{ // wireが偶数の場合
               position = dd_long * ((double)wire-2.0) / 2.0 + dd_short * ((double)wire-2.0) / 2.0 + dd_short;
            }
            
            return position;
            //return position + 2.0*dd_;
         }
         else if (split_end < wire){
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if (wire % 2 == 1){ // wireが奇数の場合
               position = dd_long * ((double)wire-1.0) / 2.0 + dd_short * ((double)wire-1.0) / 2.0;
            }
            else{ // wireが偶数の場合
               position = dd_long * ((double)wire-2.0) / 2.0 + dd_short * ((double)wire-2.0) / 2.0 + dd_short;
            }
            return position + 1.0*dd_;
         }
         else{ // ここに分離部の処理
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if (wire % 2 == 1){ // wireが奇数の場合
               position = dd_long * ((double)wire-1.0) / 2.0 + dd_short * ((double)wire-1.0) / 2.0;
            }
            else{ // wireが偶数の場合
               position = dd_long * ((double)wire-2.0) / 2.0 + dd_short * ((double)wire-2.0) / 2.0 + dd_long;
            }
            //if(wire == 137.0){// fiberID 138 に存在すると考えているファイバー欠けの例外処理
            //   position = dd_long * (wire-2.0) / 2.0 + dd_short * (wire-2.0) / 2.0 + dd_long;
            //   position += 2.*dd_;
            //}
            return position;
            //return position + 2.0*dd_;
         }
      }
      else if(layer==304){ //layer4 の時
         int border1=130;
         int border2=160;
         int border3=170;
         if(wire < border1){
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if (wire % 2 == 1){ // wireが奇数の場合
               position = dd_long * ((double)wire-1.0) / 2.0 + dd_short * ((double)wire-1.0) / 2.0;
            }
            else{ // wireが偶数の場合
               position = dd_long * ((double)wire-2.0) / 2.0 + dd_short * ((double)wire-2.0) / 2.0 + dd_long;
            }
            return position + dd_;
         }
         else if (border1<=wire && wire<border2){
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if ((int)wire % 2 == 1){ // wireが奇数の場合
               position = dd_long * ((double)wire-1.0) / 2.0 + dd_short * ((double)wire-1.0) / 2.0;
            }
            else{ // wireが偶数の場合
               position = dd_long * ((double)wire-2.0) / 2.0 + dd_short * ((double)wire-2.0) / 2.0 + dd_long;
            }
            position = wire*dd_;
            return position;
         }
         else if (border2<=wire && wire<border3){
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if ((int)wire % 2 == 1){ // wireが奇数の場合
               position = dd_long * ((double)wire-1.0) / 2.0 + dd_short * ((double)wire-1.0) / 2.0;
            }
            else{ // wireが偶数の場合
               position = dd_long * ((double)wire-2.0) / 2.0 + dd_short * ((double)wire-2.0) / 2.0 + dd_long;
            }
            position = wire*dd_;
            return position;
         }
         else if (border3 <= wire){
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if ((int)wire % 2 == 1){ // wireが奇数の場合
               position = dd_long * ((double)wire-1.0) / 2.0 + dd_short * ((double)wire-1.0) / 2.0;
            }
            else{ // wireが偶数の場合
               position = dd_long * ((double)wire-2.0) / 2.0 + dd_short * ((double)wire-2.0) / 2.0 + dd_short;
            }
            return position + 1.0*dd_;
         }
         else{ // 現状使用されないはず
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if ((int)wire % 2 == 1){ // wireが奇数の場合
               position = dd_long * ((double)wire-1.0) / 2.0 + dd_short * ((double)wire-1.0) / 2.0;
            }
            else{ // wireが偶数の場合
               position = dd_long * ((double)wire-2.0) / 2.0 + dd_short * ((double)wire-2.0) / 2.0 + dd_long;
            }
            position = wire*dd_;
            return position;
         }
      }
      else{//layer1, 4 以外
         double position;
         double dd_long = dd_*1.5;
         double dd_short= dd_*0.5;
         if ((int)wire % 2 == 1){ // wireが奇数の場合
            position = dd_long * ((double)wire-1.0) / 2.0 + dd_short * ((double)wire-1.0) / 2.0;
         }
         else{ // wireが偶数の場合
            position = dd_long * ((double)wire-2.0) / 2.0 + dd_short * ((double)wire-2.0) / 2.0 + dd_short;
         }
         return position;
      }
   }
   

   double WirePos( double wire, int layer ) const { 
      
      double pos_wire = 0;
      double pos_w0   = 0;

      if(wire == std::floor(wire)){
         pos_wire = CalcWirePos_BFT((int)wire, layer);
      }else{
         int wire1 = (int) wire;
         int wire2 = wire1 + 1;
         double pos1 = CalcWirePos_BFT(wire1, layer);
         double pos2 = CalcWirePos_BFT(wire2, layer);
         pos_wire = (pos1 + pos2) / 2.0;
      }
      double wire0 = w0_;
      if(wire0 == std::floor(wire0)){
         pos_w0 = CalcWirePos_BFT((int)wire0, layer);
      }else{
         int wire1 = (int) wire0;
         int wire2 = wire1 + 1;
         double pos1 = CalcWirePos_BFT(wire1, layer);
         double pos2 = CalcWirePos_BFT(wire2, layer);
         pos_w0 = (pos1 + pos2) / 2.0;
      }

      return pos_wire - pos_w0 + ofs_;


      // 元のやつ
      //return dd_*(wire - w0_)+ofs_; 
   }   
   int WireNumber( double pos ) const; 
   
private:
   void calcVectors( void );
   
private:
   int id_;
   std::string name_;
   ThreeVector pos_;
   double tiltAngle_, rotAngle1_, rotAngle2_;
   double length_;
   double resol_;
   double w0_, dd_, ofs_;
   
   double dxds_, dxdt_, dxdu_;
   double dyds_, dydt_, dydu_;
   double dzds_, dzdt_, dzdu_;
   
   double dsdx_, dsdy_, dsdz_;
   double dtdx_, dtdy_, dtdz_;
   double dudx_, dudy_, dudz_;
   
   friend class GeomMan;
   friend class GeomRecordComp;
};

struct GeomRecordComp 
   : public std::binary_function <GeomRecord *, GeomRecord *, bool> 
{
   bool operator()( const GeomRecord * const p1,
                    const GeomRecord * const p2 ) const
      { return p1->id_ < p2->id_; }
};

#endif
