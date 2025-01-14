/*
  FieldMan.hh

  2012/5  K.Shirotori
*/

#ifndef FieldMan_h
#define FieldMan_h 1

#include "ThreeVector.hh"
#include <string>
#include <vector>

class SpecFieldMap;
class FieldElements;

class FieldMan
{
private:
  FieldMan();
  FieldMan( const FieldMan & );
  FieldMan & operator = ( const FieldMan & );
public:
  ~FieldMan(); 

public:
  void SetFileName( const char *filename ) { filename_=filename; }
  void SetFileName( const std::string &filename ) { filename_=filename; }

  bool Initialize( void );
  bool Initialize( const char *filename )
  { filename_=filename; Initialize(); }
  bool Initialize( const std::string &filename )
  { filename_=filename; Initialize(); }

  static FieldMan & GetInstance( void );
  ThreeVector GetField( const ThreeVector & position ) const;
  ThreeVector GetdBdX( const ThreeVector & position ) const;
  ThreeVector GetdBdY( const ThreeVector & position ) const;
  ThreeVector GetdBdZ( const ThreeVector & position ) const;

  void cleanupElementsList( void );
  void AddElement( FieldElements *elem );

  double StepSize( const ThreeVector & position, 
                   double defStepSize, 
                   double MinStepSize ) const; 

private:
  static FieldMan *fMan_;
  std::string filename_;
  SpecFieldMap *Specmap_;

  typedef std::vector <FieldElements *> FEContainer;
  typedef std::vector <FieldElements *>::const_iterator FEIterator; 

  FEContainer elemList_;
};
    
#endif
