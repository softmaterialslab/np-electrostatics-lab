// This is main.
// This is Car Parrinello molecular dynamics for simulating dynamics of ions near a nanoparticle (NP) surface
// Problem : Compute the density profile of ions around a NP
/* Useful studies :	
		     1. Role of dielectric contrast
		     2. Role of valency of ions
		     3. Role of varying salt concentration
		     4. Role of NP charge
*/

#include <boost/program_options.hpp>
#include "utility.h"
#include "interface.h"
#include "particle.h"
#include "vertex.h"
#include "BIN.h"
#include "control.h"
#include "functions.h"
#include "precalculations.h"
#include "thermostat.h"
#include "fmd.h"

void cpmd(vector<PARTICLE>&, vector<VERTEX>&, INTERFACE&, vector<THERMOSTAT>&, vector<THERMOSTAT>&, vector<BIN>&, PARTICLE&, CONTROL&, CONTROL&);

using namespace boost::program_options;

int main(int argc, char* argv[]) 
{
  // Electrostatic system variables
  double radius;		// radius of the dielectric sphere
  double ein; 			// permittivity of inside medium
  double eout; 			// permittivity of outside medium
  int counterion_valency;	// counterion valency (positive by convention)
  double counterion_diameter;	// counterion diameter
  int colloid_valency;		// valency of the colloid (negative by convention in case of colloid problem). for interface problem set this variable to zero
  int salt_valency_in;		// salt valency inside
  int salt_valency_out;		// salt valency outside
  double salt_conc_in;		// salt concentration outside	(enter in M)
  double salt_conc_out;		// salt concentration outside	(enter in M)
  double saltion_diameter_in;	// inside salt ion diameter	(positive and negative ions assumed to have same diameter at this point)
  double saltion_diameter_out;  // outside salt ion diameter	(positive and negative ions assumed to have same diameter at this point)
  double T;			// temperature at which the system of ions is
  
  // Simulation related variables
  int total_gridpoints;		// total number of grid points that discretize the continuous interface
  double Q;			// thermostat mass required to generate canonical ensemble
  double fake_T;			// fake temperature useful for Car Parrinello dynamics
  double fake_Q;			// fake thermostat mass required to generate canonical ensemble for fake degrees
  unsigned int chain_length_real;
  unsigned int chain_length_fake;
  double box_radius;		// simulation box size, measured as radius in case of a sphere
  double bin_width;		// width of the bins used to compute density profiles
  CONTROL fmdremote;		// remote control for fmd
  CONTROL cpmdremote;		// remote control for cpmd

  // Different parts of the system
  INTERFACE dsphere;		// interface(s)
  PARTICLE colloid;		// nano object (acts as interface)
  vector<PARTICLE> counterion;	// counterions
  vector<PARTICLE> saltion_in;	// salt ions inside
  vector<PARTICLE> saltion_out;	// salt ions outside
  vector<PARTICLE> ion;		// all ions in the system
  vector<VERTEX> s;		// all vertices
  
  // Analysis
  vector<BIN> bin;		// bins		
  char short_run;		// run length if short invokes post analysis
  
  // Get input values from the user
  options_description desc("Usage:\nrandom_mesh <options>");
  desc.add_options()
      ("help,h", "print usage message")
      ("short_run,r", value<char>(&short_run)->default_value('n'), "short run?")
      ("radius,a", value<double>(&radius)->default_value(3.57), "sphere radius")				// enter in nanometers
      ("epsilon_in,e", value<double>(&ein)->default_value(80), "dielectric const inside")
      ("epsilon_out,E", value<double>(&eout)->default_value(80), "dielectric const outside")
      ("counterion_valency,v", value<int>(&counterion_valency)->default_value(1), "counterion valency")
      ("colloid_valency,V", value<int>(&colloid_valency)->default_value(-60), "colloid valency")
      ("salt_valency_in,z", value<int>(&salt_valency_in)->default_value(1), "salt valency inside")
      ("salt_valency_out,Z", value<int>(&salt_valency_out)->default_value(1), "salt valency outside")
      ("salt_conc_in,c", value<double>(&salt_conc_in)->default_value(0.0), "salt concentration inside")
      ("salt_conc_out,C", value<double>(&salt_conc_out)->default_value(0.0), "salt concentration outside")
      ("counterion_diameter,x", value<double>(&counterion_diameter)->default_value(0.357), "counterion diameter")			// enter in nanometers
      ("saltion_diameter_in,d", value<double>(&saltion_diameter_in)->default_value(0.357), "salt ion diameter inside")		// enter in nanometers
      ("saltion_diameter_out,D", value<double>(&saltion_diameter_out)->default_value(0.357), "salt ion diameter outside")		// enter in nanometers
      ("total_gridpoints,g", value<int>(&total_gridpoints)->default_value(132), "gridpoints")
      ("thermostat_mass,Q", value<double>(&Q)->default_value(1.0), "thermostat mass")
      ("chain_length_real,L", value<unsigned int>(&chain_length_real)->default_value(5), "chain length for real system: enter L+1 if you want L thermostats")
      ("fake_temperature,k", value<double>(&fake_T)->default_value(0.01), "fake temperature")
      ("fake_thermostat_mass,q", value<double>(&fake_Q)->default_value(1.0), "fake thermostat mass")
      ("chain_length_fake,l", value<unsigned int>(&chain_length_fake)->default_value(5), "chain length for fake system: enter L+1 if you want L thermostats")
      ("box_radius,b", value<double>(&box_radius)->default_value(14.28), "simulation box radius")		// enter in nanometers
      ("bin_width,B", value<double>(&bin_width)->default_value(0.1), "bin width")
      ("anneal_fmd,A", value<char>(&fmdremote.anneal)->default_value('n'), "anneal in fmd on?")
      ("fmd_fake_mass,m", value<double>(&fmdremote.fakemass)->default_value(1.0), "fmd fake mass")
      ("cpmd_fake_mass,M", value<double>(&cpmdremote.fakemass)->default_value(1.0), "cpmd fake mass")
      ("fmd_timestep,t", value<double>(&fmdremote.timestep)->default_value(0.001), "time step used in fmd")
      ("cpmd_timestep,T", value<double>(&cpmdremote.timestep)->default_value(0.001), "time step used in cpmd")
      ("fmd_steps,s", value<int>(&fmdremote.steps)->default_value(10000), "steps used in fmd")
      ("cpmd_steps,S", value<int>(&cpmdremote.steps)->default_value(50000), "steps used in cpmd")
      ("fmd_eqm,p", value<int>(&fmdremote.hiteqm)->default_value(100), "production begin (fmd)")
      ("cpmd_eqm,P", value<int>(&cpmdremote.hiteqm)->default_value(10000), "production begin (cpmd)")
      ("fmd_freq,f", value<int>(&fmdremote.freq)->default_value(10), "sample frequency (fmd)")
      ("cpmd_freq,F", value<int>(&cpmdremote.freq)->default_value(100), "sample frequency (cpmd)")
      ("fmd_verify,y", value<int>(&fmdremote.verify)->default_value(0), "verify (fmd)")
      ("cpmd_verify,Y", value<int>(&cpmdremote.verify)->default_value(10000), "verify (cpmd)")
      ("cpmd_writedata,U", value<int>(&cpmdremote.writedata)->default_value(1000), "write data files")
      ("cpmd_extra_compute,X", value<int>(&cpmdremote.extra_compute)->default_value(1000), "compute additional (cpmd)")
      ("cpmd_writedensity,W", value<int>(&cpmdremote.writedensity)->default_value(10000), "write density files");
      
  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);
  if (vm.count("help"))
  {
    std::cout << desc << "\n";
    return 0;
  }  
  
  cout << "\nProgram starts\n";
  cout << "Number of processors used  " << THREADSIZE << endl;
  cout << "Make sure that number of grid points / ions is greater than  " << THREADSIZE << endl;
  
  // Set up the system
  T = 1;															// set temperature
  dsphere = INTERFACE(VECTOR3D(0,0,0), radius/unitlength, ein, eout);								// make interface
  dsphere.set_up(salt_conc_in, salt_conc_out, salt_valency_in, salt_valency_out, total_gridpoints, box_radius/unitlength);	// set up properties inside and outside the interface
  colloid = PARTICLE(0, 2*radius/unitlength, colloid_valency, 1.0*colloid_valency, 1.0, dsphere.ein, VECTOR3D(0,0,0));		// make the nano object. note the diameter is supplied, not the radius
  colloid.velvec = VECTOR3D(0,0,0);												// nano object is stationary
  dsphere.put_counterions(colloid, counterion, counterion_valency, counterion_diameter, ion);					// put counterions	Note: ion contains all ions
  dsphere.put_saltions_inside(saltion_in, salt_valency_in, salt_conc_in, saltion_diameter_in, ion);				// put salt ions inside
  dsphere.put_saltions_outside(saltion_out, salt_valency_out, salt_conc_out, saltion_diameter_out, ion); 			// put salt ions outside
  dsphere.discretize(s);								// discretize interface								
  // if dielectric environment inside and outside NP are different, NPs get polarized
  if (dsphere.ein == dsphere.eout)
    dsphere.POLARIZED = false;
  else
    dsphere.POLARIZED = true;
  
  cout << "np is polarized " << dsphere.POLARIZED << endl;
  
  // NOTE: sizing the arrays employed in precalculate functions
  for (unsigned int k = 0; k < s.size(); k++)
  {
    s[k].presumgwEw.resize(s.size());
    s[k].presumgEwEq.resize(s.size());
    s[k].presumgEwEw.resize(s.size());
    s[k].presumfwEw.resize(s.size());
    s[k].presumfEwEq.resize(s.size());
    s[k].presumhEqEw.resize(s.size());
  }
  
  // could only do precalculate if CPMD
  precalculate(s, dsphere);								// precalculate 
  
  for (unsigned int k = 0; k < s.size(); k++)						// get polar coordinates for the vertices
    s[k].get_polar();
  make_bins(bin, dsphere, bin_width);							// set up bins to be used for computing density profiles
  vector<double> initial_density;
  bin_ions(ion, dsphere, initial_density, bin);						// bin the ions to get initial density profile
  
  // output to screen the parameters of the problem
  cout << "\n";
  if (colloid.q == 0 ) 
    cout << "Interface problem " << endl;
  if (colloid.q != 0) 
    cout << "Colloid problem " << endl;
  cout << "Short run " << short_run << endl;
  cout << "Reduced units: scalefactor entering in Coloumb interaction is " << scalefactor << endl;
  cout << "Other units : length (cms)" << setw(5) << unitlength*pow(10.0,-7) << setw(10) << "energy(ergs)" << setw(5) << unitenergy << setw(10) << "mass(g)" << setw(5) << unitmass << setw(10) << "time(s)" << setw(5) << unittime << endl;
  cout << "Radius of the dielectric sphere (interface) " << dsphere.radius << endl;
  cout << "Colloid charge " << colloid.q << endl;
  cout << "Permittivity inside " << dsphere.ein << endl;
  cout << "Permittivity outside " << dsphere.eout << endl;
  cout << "Contrast strength " << 2*(dsphere.eout - dsphere.ein)/(dsphere.eout+dsphere.ein) << endl;
  cout << "Counterion valency " << counterion_valency << endl;
  cout << "Salt ion valency inside " << salt_valency_in << endl;
  cout << "Salt ion valency outside " << salt_valency_out << endl;
  cout << "Counterion diameter " << counterion_diameter/unitlength << endl;
  cout << "Salt ion diameter inside " << saltion_diameter_in/unitlength << endl;
  cout << "Salt ion diameter outside " << saltion_diameter_out/unitlength << endl;
  cout << "Salt concentration inside " << salt_conc_in << endl;
  cout << "Salt concentration outside " << salt_conc_out << endl;
  cout << "Debye length inside " << dsphere.inv_kappa_in << endl;
  cout << "Debye length outside " << dsphere.inv_kappa_out << endl;
  cout << "Mean separation inside " << dsphere.mean_sep_in << endl;
  cout << "Mean separation outside " << dsphere.mean_sep_out << endl;
  cout << "Box radius " << dsphere.box_radius << endl;
  cout << "Number of counterions " << counterion.size() << endl;
  cout << "Number of salt ions inside " << saltion_in.size() << endl;
  cout << "Number of salt ions outside " << saltion_out.size() << endl;
  cout << "Temperature " << T << endl;
  cout << "Number of points discretizing the interface " << s.size() << endl;
  cout << "Binning width (uniform) " << bin[0].width << endl;
  
  // write to files
  // initial configuration
  ofstream initial_configuration("outfiles/initialconfig.dat");
  for (unsigned int i = 0; i < ion.size(); i++)
    initial_configuration << "ion" << setw(5) << ion[i].id << setw(15) << "charge" << setw(5) << ion[i].q << setw(15) << "position" << setw(15) << ion[i].posvec << endl;
  initial_configuration.close();
  // initial density
  ofstream density_profile("outfiles/initial_density_profile.dat", ios::out);
  for (unsigned int b = 0; b < initial_density.size(); b++)
    density_profile << b*bin[b].width << setw(15) << initial_density.at(b) << endl;
  density_profile.close();
  
  // some calculations before simulation begins
  cout << "Total charge inside the sphere " << dsphere.total_charge_inside(ion) << endl;
  
  // NEW NOTE : resizing the member arrays Gion and gradGion to store dynamic precalculations in fmd and cpmd force routines
  for (unsigned int k = 0; k < s.size(); k++)
  {
    s[k].Gion.resize(ion.size());
    s[k].gradGion.resize(ion.size());
  }
  
  // Fictitious molecular dynamics
  fmd(s, ion, dsphere, colloid, fmdremote, cpmdremote);
  
  // result of fmd
  ofstream induced_density ("outfiles/induced_density.dat");
  for (unsigned int k = 0; k < s.size(); k++) 
    induced_density << k+1 << setw(15) << s[k].theta  << setw(15) << s[k].phi << setw(15) << s[k].w << setw(15) << s[k].wmean << endl;
  induced_density.close();
  
  // prepare for cpmd : make real and fake baths
  
  vector<THERMOSTAT> real_bath;
  if (chain_length_real == 1)
    real_bath.push_back(THERMOSTAT(0,T,3*ion.size(),0.0,0,0));		
  else
  {
    real_bath.push_back(THERMOSTAT(Q, T, 3*ion.size(), 0, 0, 0));
    while (real_bath.size() != chain_length_real - 1)
      real_bath.push_back(THERMOSTAT(Q/(3*ion.size()), T, 1, 0, 0, 0));				
    real_bath.push_back(THERMOSTAT(0,T,3*ion.size(),0.0,0,0));			// finally, the coding trick: dummy bath (dummy bath always has zero mass)
  }
  
  vector<THERMOSTAT> fake_bath;
  if (chain_length_fake == 1)
    fake_bath.push_back(THERMOSTAT(0,fake_T,s.size(),0.0,0,0));
  else
  {
    fake_bath.push_back( THERMOSTAT(fake_Q, fake_T, s.size(), 0, 0, 0) );
    while (fake_bath.size() != chain_length_fake - 1)
      fake_bath.push_back(THERMOSTAT(fake_Q/s.size(), fake_T, 1, 0, 0, 0));
    fake_bath.push_back(THERMOSTAT(0,fake_T,s.size(),0.0,0,0));			// finally, the coding trick: dummy bath (dummy bath always has zero mass)
  }
  
  cout << "Number of chains for real system" << setw(3) << real_bath.size() - 1 << endl;
  cout << "Number of chains for fake system" << setw(3) << fake_bath.size() - 1 << endl;
  
  // Car-Parrinello Molecular Dynamics
  cpmd(ion, s, dsphere, real_bath, fake_bath, bin, colloid, fmdremote, cpmdremote);
  
  // Post simulation analysis (useful for short runs, but performed otherwise too)
  cout << "MD trust factor R (should be < 0.05) is " << compute_MD_trust_factor_R(cpmdremote.hiteqm) << endl;
  auto_correlation_function();
  
  cout << "Program ends" << endl;
  cout << endl;
  return 0;
} 
// End of main
