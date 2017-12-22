// This file contains the routine to compute the total potential energy of the system
// Electrostatic and Excluded volume contributions

#include "energies.h"

// Potential energy
double energy_functional(vector<VERTEX>& s, vector<PARTICLE>& ion, INTERFACE& dsphere, PARTICLE& colloid)
{
  // Electrostatic interaction
  // fww : potential energy due to induced charge and induced charge interaction
  // fwEw : potential energy due to induced charge and E of induced charge interaction
  // fEwEw : potential due to E of induced charge and E of induced charge interaction
  // fqqc : potential energy due to ion and central charge 
  // fqq : potential energy due to ion and ion interaction
  // fwq : potential energy due to induced charge and ion interaction
  // fqEq : potential due to ion and E due to ion interaction
  // fEwq : potential energy due to E of induced charge and ion interaction
  // fwEq : potential energy due to induced charge and E of ion interaction
  // fEqEq : potential due to E of ion and E of ion interaction
  // fEwEq : potential due to E of induced charge and E of ion interaction
  
  double potential;
  
  if (dsphere.POLARIZED)
  {
    vector<double> saveinner1(s.size(),0.0);
    vector<double> inner2(s.size(),0.0);
    vector<double> inner3(s.size(),0.0);
    vector<double> inner4(s.size(),0.0);
    vector<double> inner5(s.size(),0.0);

    unsigned int k, l;
    unsigned int i, j;
    double fqq, fwq, fqEq_qEw, fwEq_EqEq_EwEq;
    double insum;
    double ind_ind;				
    double ion_ion = 0;
    vector<double> ind_energy(s.size(),0.0);
    vector<double> ion_energy(ion.size(), 0.0);

    #pragma omp parallel default(shared) private(k, l, i, j, insum, ind_ind, fqq, fwq, fqEq_qEw, fwEq_EqEq_EwEq) num_threads(THREADSIZE)
    {
      #pragma omp for schedule(dynamic,CHUNKSIZE)
      for (k = 0; k < s.size(); k++)
      {
	insum = 0;
	for (unsigned int j = 0; j < ion.size(); j++)
	  insum += (s[k].normalvec*Grad(s[k].posvec,ion[j].posvec))*ion[j].q/ion[j].epsilon;
	saveinner1[k] = insum;
      }
	
      #pragma omp for schedule(dynamic,CHUNKSIZE)
      for (k = 0; k < s.size(); k++)
      {
	insum = 0;
	for(l = 0; l < s.size(); l++)
	  insum += s[k].ndotGradGreens[l] * s[l].w * s[l].a;
	inner2[k] = insum;
	
	insum = 0;
	for(l = 0; l < s.size(); l++)
	  insum += s[k].Greens[l] * s[l].w * s[l].a;
	inner3[k] = insum;
	
	insum = 0;
	for (l = 0; l < s.size(); l++)
	  insum += s[k].Greens[l] * saveinner1[l] * s[l].a;
	inner4[k] = insum;
	
	insum = 0;
	for (l = 0; l < s.size(); l++)
	  insum += s[k].presumfEwEq[l] * s[l].w * s[l].a;
	inner5[k] = insum;
      }

      #pragma omp for schedule(dynamic,CHUNKSIZE) nowait
      for (k = 0; k < s.size(); k++)
      {
	ind_ind = 0;
	for (l = 0; l < s.size(); l++)
	  ind_ind += s[k].w * s[k].a * ( (-0.5)*dsphere.ed*(2*dsphere.em - 1)*s[k].presumfwEw[l] +  0.5*dsphere.em*(dsphere.em - 1)*s[k].Greens[l] + 0.5*dsphere.ed*dsphere.ed*s[k].presumgEwEw[l] ) * s[l].w * s[l].a;
	ind_energy[k] = ind_ind;
      }

      #pragma omp for schedule(dynamic,CHUNKSIZE) nowait
      for (i = 0; i < ion.size(); i++)
      {
	fqq = 0;
	for (j = 0; j < ion.size(); j++)
	{
	  if (i == j) continue;
	  fqq += 0.5 * ion[i].q * ion[j].q * (1.0 / ion[i].epsilon) / ((ion[i].posvec - ion[j].posvec).GetMagnitude());
	}
      
	fwq = 0;
	for (k = 0; k < s.size(); k++)
	  fwq += ion[i].q * (0.5 - 0.5 * dsphere.em / ion[i].epsilon) * (1/((s[k].posvec - ion[i].posvec).GetMagnitude())) * s[k].w * s[k].a;

	insum = 0;
	for (k = 0; k < s.size(); k++)
	  insum += (1.0/((ion[i].posvec-s[k].posvec).GetMagnitude())) * (saveinner1[k] + inner2[k]) * s[k].a;
	fqEq_qEw = 0.5 * dsphere.ed * ion[i].q * insum / ion[i].epsilon;

	insum = 0;
	for(k = 0; k < s.size(); k++)
	  insum += (s[k].normalvec * Grad(s[k].posvec,ion[i].posvec)) * ( (-1) * 0.5 * dsphere.ed * (2*dsphere.em - 1) * inner3[k] + 0.5 * dsphere.ed * dsphere.ed * inner4[k] + dsphere.ed * dsphere.ed * inner5[k] ) * s[k].a;
	fwEq_EqEq_EwEq = (ion[i].q / ion[i].epsilon) * insum; 

	ion_energy[i] = fqq + fwq + fqEq_qEw + fwEq_EqEq_EwEq + colloid.q * ion[i].q * (1.0 / ion[i].epsilon) / ( (ion[i].posvec - colloid.posvec).GetMagnitude() );
      }
    }

    // ion-ion + ion-induced charge energy
    ion_ion = 0;
    for (unsigned int i = 0; i < ion.size(); i++)
      ion_ion = ion_ion + ion_energy[i];
    ind_ind = 0;
    for (unsigned int k = 0; k < s.size(); k++)
      ind_ind = ind_ind + ind_energy[k];
    
    // electrostatic potential energy
    potential = (ion_ion + ind_ind) * scalefactor;
  }
  
  else	// if not POLARIZED
  {
    unsigned int i, j;
    double fqq;
    vector<double> ion_energy(ion.size(), 0.0);
    #pragma omp parallel default(shared) private(i, j, fqq) num_threads(THREADSIZE)
    {
      #pragma omp for schedule(dynamic,CHUNKSIZE) nowait
      for (i = 0; i < ion.size(); i++)
      {
	fqq = 0;
	for (j = 0; j < ion.size(); j++)
	{
	  if (i == j) continue;
	  fqq += 0.5 * ion[i].q * ion[j].q * (1.0 / ion[i].epsilon) / ((ion[i].posvec - ion[j].posvec).GetMagnitude());
	}
	ion_energy[i] = fqq + colloid.q * ion[i].q * (1.0 / ion[i].epsilon) / ( (ion[i].posvec - colloid.posvec).GetMagnitude() );
      }
    }

    double ion_ion = 0;
    for (unsigned int i = 0; i < ion.size(); i++)
      ion_ion = ion_ion + ion_energy[i];

    // electrostatic potential energy
    potential = (ion_ion) * scalefactor;
  }

  // Excluded volume interaction energy given by purely repulsive LJ
  
  // ion-sphere (ions outsdie)
  // make a dummy particle with the same diameter as the ion just beneath the interface
  vector<double> lj1;
  for (unsigned int i = 0; i < ion.size(); i++)
  {
    double uljcs = 0.0;
    if (ion[i].posvec.GetMagnitude() < dsphere.radius)
    {
      lj1.push_back(uljcs);
      continue;
    }
    PARTICLE dummy = PARTICLE(0,ion[i].diameter,0,0,0,dsphere.ein,ion[i].posvec^((dsphere.radius - 0.5*ion[i].diameter)/ion[i].posvec.GetMagnitude()));
    VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
    double r = r_vec.GetMagnitude();
    double d = 0.5 * (ion[i].diameter + dummy.diameter);
    double elj = 1.0;
    if (r < dcut * d)
    {
      double r2 = r * r;
      double r6 = r2 * r2 * r2;
      double d2 = d * d;
      double d6 = d2 * d2 * d2;
      uljcs = 4 * elj *  (d6 / r6) * ( ( d6 / r6 ) -  1 )  + elj;
    }
    lj1.push_back(uljcs);
  }
  
  // ion-ion
  vector<double> lj2;
  for (unsigned int i = 0; i < ion.size(); i++)
  {
    double uljcc = 0.0;
    for (unsigned int j = 0; j < ion.size(); j++)
    {
      if (j == i) continue;
      VECTOR3D r_vec = ion[i].posvec - ion[j].posvec;
      double r = r_vec.GetMagnitude();
      double d = 0.5 * (ion[i].diameter + ion[j].diameter);
      double elj = 1.0;
      if (r < dcut * d)
      {
	double r2 = r * r;
	double r6 = r2 * r2 * r2;
	double d2 = d * d;
	double d6 = d2 * d2 * d2;
	uljcc = uljcc +  4 * elj * (d6 / r6) * ( ( d6 / r6 ) - 1 ) + elj;
      }
      else
	uljcc = uljcc + 0.0;
    }
    lj2.push_back(uljcc);
  }
  
  // ion-box
  // make a dummy particle with the same diameter as the ion just above the simulation box
  vector<double> lj3;
  for (unsigned int i = 0; i < ion.size(); i++)
  {
    double uljcb = 0.0;
    PARTICLE dummy = PARTICLE(0,ion[i].diameter,0,0,0,dsphere.eout,ion[i].posvec^((dsphere.box_radius + 0.5*ion[i].diameter)/ion[i].posvec.GetMagnitude()));
    VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
    double r = r_vec.GetMagnitude();
    double d = 0.5 * (ion[i].diameter + dummy.diameter);
    double elj = 1.0;
    if (r < dcut * d)
    {
      double r2 = r * r;
      double r6 = r2 * r2 * r2;
      double d2 = d * d;
      double d6 = d2 * d2 * d2;
      uljcb = 4 * elj * (d6 / r6) * ( (d6 / r6) - 1 ) + elj;
    }
    lj3.push_back(uljcb);
  }
  
  // ion-sphere (ions inside)
  // make a dummy particle with the same diameter as the ion just above the interface
  vector<double> lj4;
  for (unsigned int i = 0; i < ion.size(); i++)
  {
    double uljcs = 0.0;					
    if (ion[i].posvec.GetMagnitude() > dsphere.radius)
    {
      lj4.push_back(uljcs);
      continue;
    }
    PARTICLE dummy = PARTICLE(0,ion[i].diameter,0,0,0,dsphere.eout,ion[i].posvec^((dsphere.radius + 0.5*ion[i].diameter)/ion[i].posvec.GetMagnitude()));
    VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
    double r = r_vec.GetMagnitude();
    double d = 0.5 * (ion[i].diameter + dummy.diameter);
    double elj = 1.0;
    if (r < dcut * d)
    {
      double r2 = r * r;
      double r6 = r2 * r2 * r2;
      double d2 = d * d;
      double d6 = d2 * d2 * d2;
      uljcs = 4 * elj * (d6 / r6) * ( (d6 / r6) - 1 ) + elj;
    }
    lj4.push_back(uljcs);
  }
  
  double lj1energy = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    lj1energy += lj1[i];
  
  double lj2energy = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    lj2energy += lj2[i];
  lj2energy = 0.5 * lj2energy;
  
  double lj3energy = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    lj3energy += lj3[i];
  
  double lj4energy = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    lj4energy += lj4[i];
  
  potential = potential + lj1energy + lj2energy + lj3energy + lj4energy;
  
  return potential;
}

