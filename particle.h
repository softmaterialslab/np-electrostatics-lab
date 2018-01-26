// This is particle class

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "utility.h"
#include "vector3d.h"
#include "thermostat.h"

class PARTICLE 
{
  public:

  // members
  int id;		// id of the particle
  double diameter;	// diameter of the particle
  int valency;		// valency of the ion
  double q;		// charge of the particle
  double m; 		// mass of the particle
  double epsilon;	// dielectric constant of the medium
  VECTOR3D posvec;	// position vector of the particle
  VECTOR3D velvec;	// velocity vector of the particle
  VECTOR3D forvec;	// force vector on the particle
  double pe;		// potential energy
  long double ke;	// kinetic energy
  double energy;	// energy
  
  // member functions
  
  // make a particle
  PARTICLE(int get_id = 0, double get_diameter = 0, int get_valency = 0, double get_charge = 0, double get_mass = 0, double get_medium_dielectric_constant = 0, VECTOR3D get_position = VECTOR3D(0,0,0))
  {
    id = get_id;
    diameter = get_diameter;
    valency = get_valency;
    q = get_charge;
    m = get_mass;
    epsilon = get_medium_dielectric_constant;
    posvec = get_position;
  }
  
  // update position of the particle
  void update_position(double dt)			
  {
    posvec = ( posvec + (velvec ^ dt) );
    return;
  }
  
  // update velocity of the particle
  void update_velocity(double dt)	
  {
    velvec = ( velvec + ( forvec ^ ( 0.5 * dt ) ) );
    return;
  }
  
  void new_update_velocity(double dt, THERMOSTAT main_bath, long double expfac)
  {
    velvec = ( ( velvec ^ (expfac)  ) + ( forvec ^ (0.5 * dt * sqrt(expfac)) ) );
    return;
  }
  
  // calculate kinetic energy of a particle
  void kinetic_energy()				
  {
    ke = 0.5 * m * velvec.GetMagnitude() * velvec.GetMagnitude();
    return;
  }
};

#endif
