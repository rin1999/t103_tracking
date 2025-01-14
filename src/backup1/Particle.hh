/*
  Particle.hh
  
  2015/12  K.Shirotori
*/

#ifndef Particle_h 
#define Particle_h

#include <cstddef>
#include <vector>

typedef std::vector <double> DoubleVec;
typedef std::vector <int> IntVec;

class Particle
{
private:
  IntVec  id_, pid_;
  DoubleVec momx_, momy_, momz_;
  IntVec fdet_;

public:
  Particle()
    : id_(0), pid_(0), 
      momx_(0), momy_(0), momz_(0),
      fdet_(0)
  {};
  ~Particle() {};

public:
  void SetId( int id ) { id_.push_back(id); }
  void SetPid( int pid ) { pid_.push_back(pid); }
  void SetMomX( double momx ) { momx_.push_back(momx); }
  void SetMomY( double momy ) { momy_.push_back(momy); }
  void SetMomZ( double momz ) { momz_.push_back(momz); }
  void SetDetFlag( int fdet ) { fdet_.push_back(fdet); }

  int    GetSize( void ) const { return id_.size(); };
  int    GetId( int n ) const { return id_[n]; };
  int    GetPid( int n ) const { return pid_[n]; };
  double GetMomX( int n ) const { return momx_[n]; };
  double GetMomY( int n ) const { return momy_[n]; };
  double GetMomZ( int n ) const { return momz_[n]; };
  int    GetDetFlag( int n ) const { return fdet_[n]; };
};
#endif
