/*
  FieldMan.cc

  2012/5  K.Shirotori
*/

#include "FieldMan.hh"
#include "SpecFieldMap.hh"
#include "FieldElements.hh"

FieldMan *FieldMan::fMan_ = 0;

const double Delta = 0.1;

FieldMan::FieldMan()
  : Specmap_(0)
{}

FieldMan::~FieldMan()
{
  delete Specmap_;
}

FieldMan & FieldMan::GetInstance( void )
{
  if( !fMan_ ){
    fMan_ = new FieldMan();
  }
  return *fMan_;
}

bool FieldMan::Initialize( void )
{
  SpecFieldMap *pSpec = new SpecFieldMap( filename_.c_str() );
  if( pSpec ){
    delete Specmap_;
    Specmap_ = pSpec;
  }
  if( Specmap_ )
    return Specmap_->Initialize();
  else
    return false;
}

ThreeVector FieldMan::GetField( const ThreeVector & position ) const
{
  ThreeVector field( 0.0, 0.0, 0.0 );
  if( Specmap_ ){
    double p[3], b[3];
    p[0]=position.x()*0.1; p[1]=position.y()*0.1; p[2]=position.z()*0.1; 
    if( Specmap_->GetFieldValue( p, b ) ){
      field.setX(b[0]); field.setY(b[1]); field.setZ(b[2]); 
    }
  }

#if 1
  FEIterator end = elemList_.end();
  for( FEIterator itr=elemList_.begin(); itr!=end; ++itr ){
    if( (*itr)->ExistField( position ) )
      field += (*itr)->GetField( position );
  }   
#endif

  return field;
}

ThreeVector FieldMan::GetdBdX( const ThreeVector & position ) const
{
  ThreeVector p1 = position + ThreeVector( Delta, 0., 0. );
  ThreeVector p2 = position - ThreeVector( Delta, 0., 0. );

  ThreeVector B1 = GetField( p1 );
  ThreeVector B2 = GetField( p2 );
 
  return 0.5*(B1-B2)/Delta;
} 

ThreeVector FieldMan::GetdBdY( const ThreeVector & position ) const
{
  ThreeVector p1 = position + ThreeVector( 0., Delta,  0. );
  ThreeVector p2 = position - ThreeVector( 0., Delta, 0. );

  ThreeVector B1 = GetField( p1 );
  ThreeVector B2 = GetField( p2 );

  return 0.5*(B1-B2)/Delta;
} 

ThreeVector FieldMan::GetdBdZ( const ThreeVector & position ) const
{
  ThreeVector p1 = position + ThreeVector( 0., 0., Delta );
  ThreeVector p2 = position - ThreeVector( 0., 0., Delta );

  ThreeVector B1 = GetField( p1 );
  ThreeVector B2 = GetField( p2 );

  return 0.5*(B1-B2)/Delta;
} 

void FieldMan::cleanupElementsList( void )
{
  elemList_.clear();
}

void FieldMan::AddElement( FieldElements *elem )
{
  elemList_.push_back( elem );
}

double FieldMan::StepSize( const ThreeVector &position,
                           double defStepSize,
                           double MinStepSize ) const
{
  double d=fabs(defStepSize);
  double s=defStepSize/d;
  MinStepSize = fabs(MinStepSize);

  bool flag = true;
  FEIterator end = elemList_.end();
  while ( flag && d>MinStepSize ){
    for( FEIterator itr=elemList_.begin(); itr!=end; ++itr ){
      if( (*itr)->checkRegion( position, d) != FERSurface ){
        flag=false;
        break;
      }
    }
    d *= 0.5;
  }
  if( flag ) return s*MinStepSize;
  else       return s*d;
}
