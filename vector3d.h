// This is 3d vector class
// VECTOR3D is a vector in cartesian coordinate system

#ifndef _VECTOR3D_H
#define _VECTOR3D_H

class VECTOR3D 
{
  public:
    long double x, y, z;									// component along each axis (cartesian)
    VECTOR3D(long double xx = 0.0, long double yy = 0.0, long double zz = 0.0) : x(xx), y(yy), z(zz) 	// make a 3d vector
    {
    } 
    long double GetMagnitude()								// magnitude of the vector
    {
      return sqrt(x*x + y*y +z*z);							
    }
    VECTOR3D operator-(const VECTOR3D& vec)						// subtract two vectors
    {
      return VECTOR3D(x - vec.x, y - vec.y, z - vec.z);					
    }
    long double operator*(const VECTOR3D& vec)						// dot product of two vectors
    {
      return x*vec.x + y*vec.y + z*vec.z;						
    }
    VECTOR3D operator^(long double scalar)							// product of a vector and a scalar
    {
      return VECTOR3D(x*scalar, y*scalar, z*scalar);					
    }
    VECTOR3D operator+(const VECTOR3D& vec)						// add two vectors
    {
      return VECTOR3D(x + vec.x, y + vec.y, z + vec.z);					
    }
    bool operator==(const VECTOR3D& vec)						// compare two vectors
    {
      if (x==vec.x && y==vec.y && z==vec.z) 
	return true;
      return false;
    } 
};

#endif
