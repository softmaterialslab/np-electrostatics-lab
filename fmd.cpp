// This is fictitious molecular dynamics
// This program is used to estimate the correct w(k)'s on the interface

#include "fmd.h"

void fmd(vector<VERTEX>& s, vector<PARTICLE>& ion, INTERFACE& dsphere, PARTICLE& colloid, CONTROL& fmdremote, CONTROL& cpmdremote) 
{
  // Part I : Initialize
  for (unsigned int k = 0; k < s.size(); k++) 	
  {										
    s[k].mu = fmdremote.fakemass * s[k].a * s[k].a;						// Assign mass to the fake degree
    s[k].w = 0.0;								// Initialize fake degree value		(unconstrained)
    s[k].vw = 0.0;								// Initialize fake degree velocity	(unconstrained)
  } 
  long double sigma = constraint(s, ion, dsphere);
  for (unsigned int k = 0; k < s.size(); k++)
    s[k].w = s[k].w - sigma / (s[k].a * s.size());				// Satisfy constraint 
  long double sigmadot = dotconstraint(s);  
  for (unsigned int k = 0; k < s.size(); k++)
    s[k].vw = s[k].vw - sigmadot / (s[k].a * s.size());				// Satisfy time derivative of the constraint
  for_fmd_calculate_force(s, ion, dsphere, colloid);				// Compute initial force
  long double kinetic_energy = fake_kinetic_energy(s);				// Compute initial kinetic energy
  double potential_energy = energy_functional(s, ion, dsphere, colloid);	// Compute initial potential energy
  fmdremote.annealfreq = 1000;
  
  // Output fmd essentials
  cout << "\n";
  cout << "F M D" << " at " << fmdremote.verify << endl;
  cout << "Mass assigned to the fake degrees " << s[0].mu << endl;
  cout << "Total induced charge on the interface " << dsphere.total_induced_charge(s) << endl;
  cout << "Constraint is (zero if satisfied) " << constraint(s, ion, dsphere) << endl;
  cout << "Time derivative of the constraint is " << dotconstraint(s) << endl;
  cout << "Initial force on fake degree at vertex 0 " << s[0].fw << endl;
  cout << "Initial fake kinetic energy " << kinetic_energy << endl;
  cout << "Inital potential energy " << potential_energy << endl;
  cout << "Initial total energy " << kinetic_energy + potential_energy << endl;
  cout << "Time step " << fmdremote.timestep << endl;
  cout << "Number of steps " << fmdremote.steps << endl;
  
  // create fmd output files
  ofstream fmdv, fmde, fmdtic;		
  char data[200];
  sprintf(data, "outfiles/fmdv_%.06d.dat", fmdremote.verify);
  fmdv.open(data);
  sprintf(data, "outfiles/fmde_%.06d.dat", fmdremote.verify);
  fmde.open(data);
  sprintf(data, "outfiles/fmdtic_%.06d.dat", fmdremote.verify);
  fmdtic.open(data);
  
  double average_total_induced_charge = 0.0;					// initialize useful averages
  for (unsigned int k = 0; k < s.size(); k++)
    s[k].wmean = 0;
  int samples = 0;								// number of samples

  fmdremote.extra_compute = 100;						// new addition, if one wants to have longer fmd since the startind density is not coming out right
  
  // PART II : Propagate
  /*.......................................Fictitious molecular dynamics.......................................*/ 
  for (int num = 1; num <= fmdremote.steps; num++) 
  {
    fmdv << num << "  " << s[0].w << "  " << s[0].vw << "  " << s[0].fw << "  " << dsphere.total_induced_charge(s) << endl;
    
    // total induced charge
    fmdtic << num << "  " << dsphere.total_induced_charge(s) << endl;
    
    // INTEGRATOR
    //! begins
    for (unsigned int k = 0; k < s.size(); k++)
      s[k].update_velocity(fmdremote.timestep);		// update velocity half time step
    for (unsigned int k = 0; k < s.size(); k++)
      s[k].update_position(fmdremote.timestep);		// update position full time step 
    SHAKE(s, ion, dsphere, fmdremote);			// shake to ensure constraint  
    for_fmd_calculate_force(s, ion, dsphere, colloid);	// calculate forces using the new positions
    for (unsigned int k = 0; k < s.size(); k++)
      s[k].update_velocity(fmdremote.timestep);		// update velocity full time step
    RATTLE(s);						// rattle to ensure time derivative of the constraint
    //! ends
    
    // anneal (only if necessary)
    if (fmdremote.anneal == 'y' && num >= fmdremote.hiteqm &&  num % fmdremote.annealfreq == 0)
    {
      for (unsigned int k = 0; k < s.size(); k++)
	s[k].vw = 0.9 * s[k].vw;
    }
    
    // sample collection for averages
    if (num > fmdremote.hiteqm && (num % fmdremote.freq) == 0)
    {
      samples = samples + 1;
      for (unsigned int k = 0; k < s.size(); k++)
	s[k].wmean += s[k].w;
      average_total_induced_charge += dsphere.total_induced_charge(s);
    }
    
    
    // additional computations (turn off for faster simulations)
    if (num % fmdremote.extra_compute == 0)
    {
      long double kinetic_energy = fake_kinetic_energy(s);
      double potential_energy = energy_functional(s, ion, dsphere, colloid);
      double extended_energy = kinetic_energy + potential_energy;
      fmde << num << "  " << extended_energy << "  " << kinetic_energy << "  " <<  potential_energy << endl;
    }
  }					
  /*.................................................fmd ends..................................................*/

  // Part III: Compute averages
  for (unsigned int k = 0; k < s.size(); k++) 
    if (samples != 0) s[k].wmean /= samples;				// average induced charge at each vertex
  average_total_induced_charge /= samples;		// average total induced charge
  
  cout << "Number of samples used to estimate induced charges " << samples << endl;
  cout << "Total induced charge on average " << average_total_induced_charge << endl;      
  return;
}

