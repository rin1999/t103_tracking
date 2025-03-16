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

   double CalcWirePos_BFT( double wire, int layer ) const {
      //layer1だけ特別処理を行う（layer4も追加した）
      
      if(layer==301){ //layer1 の時
         double split_start=119.0;
         double split_end  =136.0;
         if(wire < split_start){
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if (std::floor(wire) == wire){ // wireが整数の場合
               if ((int)wire % 2 == 1){ // wireが奇数の場合
                  position = dd_long * (wire-1.0) / 2.0 + dd_short * (wire-1.0) / 2.0;
               }
               else{ // wireが偶数の場合
                  position = dd_long * (wire-2.0) / 2.0 + dd_short * (wire-2.0) / 2.0 + dd_short;
               }
            }
            else{ // wireが半整数の場合（一応wire=1.7とかの場合でも対応できるように書いておこう）
               if ((int)std::floor(wire) % 2 == 1){ // wireの整数部分が奇数の場合
                  position = dd_long * ((double)std::floor(wire)-1.0) / 2.0 + dd_short * ((double)std::floor(wire)-1.0) / 2.0 + dd_short / 2.0;
               }
               else{ // wireの整数部分が偶数の場合
                  position = dd_long * ((double)std::floor(wire)-2.0) / 2.0 + dd_short * ((double)std::floor(wire)-2.0) / 2.0 + dd_short + dd_long / 2.0;
               }
            
            }
            return position;
            //return position + 2.0*dd_;
         }
         else if (split_end < wire){
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if (std::floor(wire) == wire){ // wireが整数の場合
               if ((int)wire % 2 == 1){ // wireが奇数の場合
                  position = dd_long * (wire-1.0) / 2.0 + dd_short * (wire-1.0) / 2.0;
               }
               else{ // wireが偶数の場合
                  position = dd_long * (wire-2.0) / 2.0 + dd_short * (wire-2.0) / 2.0 + dd_short;
               }
            }
            else{ // wireが半整数の場合（一応wire=1.7とかの場合でも対応できるように書いておこう）
               if ((int)std::floor(wire) % 2 == 1){ // wireの整数部分が奇数の場合
                  position = dd_long * ((double)std::floor(wire)-1.0) / 2.0 + dd_short * ((double)std::floor(wire)-1.0) / 2.0 + dd_short / 2.0;
               }
               else{ // wireの整数部分が偶数の場合
                  position = dd_long * ((double)std::floor(wire)-2.0) / 2.0 + dd_short * ((double)std::floor(wire)-2.0) / 2.0 + dd_short + dd_long / 2.0;
               }
            }
            return position + 1.0*dd_;
            //return position + 3.0*dd_; // 後半は外れたファイバー分だけ位置がずれる
         }
         else{ // ここに分離部の処理
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if (std::floor(wire) == wire){ // wireが整数の場合
               if ((int)wire % 2 == 1){ // wireが奇数の場合
                  position = dd_long * (wire-1.0) / 2.0 + dd_short * (wire-1.0) / 2.0;
               }
               else{ // wireが偶数の場合
                  position = dd_long * (wire-2.0) / 2.0 + dd_short * (wire-2.0) / 2.0 + dd_long;
               }
            }
            else{ // wireが半整数の場合（一応wire=1.7とかの場合でも対応できるように書いておこう）
               if ((int)std::floor(wire) % 2 == 1){ // wireの整数部分が奇数の場合
                  position = dd_long * ((double)std::floor(wire)-1.0) / 2.0 + dd_short * ((double)std::floor(wire)-1.0) / 2.0 + dd_long / 2.0;
               }
               else{ // wireの整数部分が偶数の場合
                  position = dd_long * ((double)std::floor(wire)-2.0) / 2.0 + dd_short * ((double)std::floor(wire)-2.0) / 2.0 + dd_long + dd_short / 2.0;
               }
            }
            if(wire == 137.0){// fiberID 138 に存在すると考えているファイバー欠けの例外処理
               position = dd_long * (wire-2.0) / 2.0 + dd_short * (wire-2.0) / 2.0 + dd_long;
               position += 2.*dd_;
            }
            return position;
            //return position + 2.0*dd_;
         }
      }
      else if(layer==304){ //layer4 の時
         double split_start=130.0;
         double split_end  =170.0;
         if(wire < split_start){
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if (std::floor(wire) == wire){ // wireが整数の場合
               if ((int)wire % 2 == 1){ // wireが奇数の場合
                  position = dd_long * (wire-1.0) / 2.0 + dd_short * (wire-1.0) / 2.0;
               }
               else{ // wireが偶数の場合
                  position = dd_long * (wire-2.0) / 2.0 + dd_short * (wire-2.0) / 2.0 + dd_long;
               }
            }
            else{ // wireが半整数の場合（一応wire=1.7とかの場合でも対応できるように書いておこう）
               if ((int)std::floor(wire) % 2 == 1){ // wireの整数部分が奇数の場合
                  position = dd_long * ((double)std::floor(wire)-1.0) / 2.0 + dd_short * ((double)std::floor(wire)-1.0) / 2.0 + dd_long / 2.0;
               }
               else{ // wireの整数部分が偶数の場合
                  position = dd_long * ((double)std::floor(wire)-2.0) / 2.0 + dd_short * ((double)std::floor(wire)-2.0) / 2.0 + dd_long + dd_short / 2.0;
               }
            }
            return position;
         }
         else if (split_end < wire){
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if (std::floor(wire) == wire){ // wireが整数の場合
               if ((int)wire % 2 == 1){ // wireが奇数の場合
                  position = dd_long * (wire-1.0) / 2.0 + dd_short * (wire-1.0) / 2.0;
               }
               else{ // wireが偶数の場合
                  position = dd_long * (wire-2.0) / 2.0 + dd_short * (wire-2.0) / 2.0 + dd_short;
               }
            }
            else{ // wireが半整数の場合（一応wire=1.7とかの場合でも対応できるように書いておこう）
               if ((int)std::floor(wire) % 2 == 1){ // wireの整数部分が奇数の場合
                  position = dd_long * ((double)std::floor(wire)-1.0) / 2.0 + dd_short * ((double)std::floor(wire)-1.0) / 2.0 + dd_short / 2.0;
               }
               else{ // wireの整数部分が偶数の場合
                  position = dd_long * ((double)std::floor(wire)-2.0) / 2.0 + dd_short * ((double)std::floor(wire)-2.0) / 2.0 + dd_short + dd_long / 2.0;
               }
            }
            return position + 1.0*dd_;
            //return position + 3.0*dd_; // 後半は外れたファイバー分だけ位置がずれる
         }
         else{ // ここに分離部の処理
            double position;
            double dd_long = dd_*1.5;
            double dd_short= dd_*0.5;
            if (std::floor(wire) == wire){ // wireが整数の場合
               if ((int)wire % 2 == 1){ // wireが奇数の場合
                  position = dd_long * (wire-1.0) / 2.0 + dd_short * (wire-1.0) / 2.0;
               }
               else{ // wireが偶数の場合
                  position = dd_long * (wire-2.0) / 2.0 + dd_short * (wire-2.0) / 2.0 + dd_long;
               }
            }
            else{ // wireが半整数の場合（一応wire=1.7とかの場合でも対応できるように書いておこう）
               if ((int)std::floor(wire) % 2 == 1){ // wireの整数部分が奇数の場合
                  position = dd_long * ((double)std::floor(wire)-1.0) / 2.0 + dd_short * ((double)std::floor(wire)-1.0) / 2.0 + dd_long / 2.0;
               }
               else{ // wireの整数部分が偶数の場合
                  position = dd_long * ((double)std::floor(wire)-2.0) / 2.0 + dd_short * ((double)std::floor(wire)-2.0) / 2.0 + dd_long + dd_short / 2.0;
               }
            }
            position = (wire - 1.0)/2.0;
            return position;
            //return position + 2.0*dd_;
         }
      }
      else{//layer1, 4 以外
         double position;
         double dd_long = dd_*1.5;
         double dd_short= dd_*0.5;
         if (std::floor(wire) == wire){ // wireが整数の場合
            if ((int)wire % 2 == 1){ // wireが奇数の場合
               position = dd_long * (wire-1.0) / 2.0 + dd_short * (wire-1.0) / 2.0;
            }
            else{ // wireが偶数の場合
               position = dd_long * (wire-2.0) / 2.0 + dd_short * (wire-2.0) / 2.0 + dd_short;
            }
         }
         else{ // wireが半整数の場合（一応wire=1.7とかの場合でも対応できるように書いておこう）
            if ((int)std::floor(wire) % 2 == 1){ // wireの整数部分が奇数の場合
               position = dd_long * ((double)std::floor(wire)-1.0) / 2.0 + dd_short * ((double)std::floor(wire)-1.0) / 2.0 + dd_short / 2.0;
            }
            else{ // wireの整数部分が偶数の場合
               position = dd_long * ((double)std::floor(wire)-2.0) / 2.0 + dd_short * ((double)std::floor(wire)-2.0) / 2.0 + dd_short + dd_long / 2.0;
            }
         }
         return position;
      }


      /*
      double position;
      double dd_long = dd_*1.5;
      double dd_short= dd_*0.5;
      if (std::floor(wire) == wire){ // wireが整数の場合
         if ((int)wire % 2 == 1){ // wireが奇数の場合
            position = dd_long * (wire-1.0) / 2.0 + dd_short * (wire-1.0) / 2.0;
         }
         else{ // wireが偶数の場合
            position = dd_long * (wire-2.0) / 2.0 + dd_short * (wire-2.0) / 2.0 + dd_short;
         }
      }
      else{ // wireが半整数の場合（一応wire=1.7とかの場合でも対応できるように書いておこう）
         if ((int)std::floor(wire) % 2 == 1){ // wireの整数部分が奇数の場合
            position = dd_long * ((double)std::floor(wire)-1.0) / 2.0 + dd_short * ((double)std::floor(wire)-1.0) / 2.0 + dd_short / 2.0;
         }
         else{ // wireの整数部分が偶数の場合
            position = dd_long * ((double)std::floor(wire)-2.0) / 2.0 + dd_short * ((double)std::floor(wire)-2.0) / 2.0 + dd_short + dd_long / 2.0;
         }
      }
      return position;
      */

      //偶奇反転バージョン
      /*
      double position;
      double dd_long = dd_*1.5;
      double dd_short= dd_*0.5;
      if (std::floor(wire) == wire){ // wireが整数の場合
         if ((int)wire % 2 == 1){ // wireが奇数の場合
            position = dd_long * (wire-1.0) / 2.0 + dd_short * (wire-1.0) / 2.0;
         }
         else{ // wireが偶数の場合
            position = dd_long * (wire-2.0) / 2.0 + dd_short * (wire-2.0) / 2.0 + dd_long;
         }
      }
      else{ // wireが半整数の場合（一応wire=1.7とかの場合でも対応できるように書いておこう）
         if ((int)std::floor(wire) % 2 == 1){ // wireの整数部分が奇数の場合
            position = dd_long * ((double)std::floor(wire)-1.0) / 2.0 + dd_short * ((double)std::floor(wire)-1.0) / 2.0 + dd_long / 2.0;
         }
         else{ // wireの整数部分が偶数の場合
            position = dd_long * ((double)std::floor(wire)-2.0) / 2.0 + dd_short * ((double)std::floor(wire)-2.0) / 2.0 + dd_long + dd_short / 2.0;
         }
      }
      return position;
      */
      
   }
   
   double WirePos( double wire, int layer ) const { 

      // 改造案
      // 半整数ワイヤーの場合の処理をCalcWirePos_BFTに任せずにWirePos内で処理したい。
      // つまり、CalcWirePos_BFTは整数値のfiberIDを受け付けて座標を返すだけの関数にしたい。
      // （ちょっと現状だとCalcWirePos_BFTの中身が汚すぎる）
      
      double pos_wire = CalcWirePos_BFT(wire, layer);
      double pos_w0   = CalcWirePos_BFT(w0_, layer);

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
