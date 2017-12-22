// This is a header file for the INTERFACE class.  

#ifndef _INTERFACE_H
#define _INTERFACE_H

#include "utility.h"
#include "vector3d.h"
#include "particle.h"
#include "vertex.h"

class INTERFACE
{
  public:

  VECTOR3D posvec;		// position vector of the inteface, for sphere its center
  double radius;		// radius of the dielectric sphere, interface geometry variable
  double ein; 			// permittivity of inside medium
  double eout; 			// permittivity of outside medium
  double em;			// permittivity at the interface, mean
  double ed;			// permittivity change at the inteface, difference scaled by 4 pi
  double box_radius;		// radius of the spherical box in which all the action occurs
  
  double lB_in;			// Bjerrum length inside
  double lB_out;		// Bjerrum length outside
  double inv_kappa_in;		// debye length inside
  double inv_kappa_out;		// debye length outside	
  double mean_sep_in;		// mean separation inside
  double mean_sep_out;		// mean separation outside
  
  int number_of_vertices;	// number of points used to discretize the interface
  
  bool POLARIZED;		// is the nanoparticle polarized; depends on ein, eout
  
  void set_up(double, double, double, double, int, double);
  void put_counterions(PARTICLE&, vector<PARTICLE>&, int, double, vector<PARTICLE>&);
  void put_saltions_inside(vector<PARTICLE>&, int, double, double, vector<PARTICLE>&);
  void put_saltions_outside(vector<PARTICLE>&, int, double, double, vector<PARTICLE>&);
  void discretize(vector<VERTEX>&);
  
  INTERFACE(VECTOR3D posvec = VECTOR3D(0,0,0), double radius = 0, double ein = 1, double eout = 1) : posvec(posvec), radius(radius), ein(ein), eout(eout)
  {
  }
  
  // total charge inside
  double total_charge_inside(vector<PARTICLE>& ion)
  {
    double charge = 0;
    for (unsigned int i = 0; i < ion.size(); i++)
      if (ion[i].posvec.GetMagnitude() < radius) charge += ion[i].q;	
    return charge;
  }
  
  // total induced charge
  double total_induced_charge(vector<VERTEX>& s)
  {
    double charge = 0;
    for (unsigned int k = 0; k < s.size(); k++)
      charge += s[k].w * s[k].a;
    return charge;
  }
};

#endif

