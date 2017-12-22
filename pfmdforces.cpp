// This file contains the routine that computes the force
// on the induced charge at vertex k and the force on the particle i
// for all k and i

#include "forces.h"

// Total Force on all degrees of freedom
void for_fmd_calculate_force(vector<VERTEX>& s, vector<PARTICLE>& ion, INTERFACE& dsphere, PARTICLE& colloid) 
{
  
  // force calculation for fake degrees of freedom
  // gwq : force due to induced charge (w) - ion (q) interaction
  // gww : force due to induced charge (w) - induced charge (w) interaction
  // gwEw : force due to induced charge (w) - Electric field due to induced charge (Ew) interaction
  // gEwEw : force due to Electric field due to induced charge (Ew) - Electric field due to induced charge (Ew) interaction
  // gEwq : force due to electric field of induced charge (Ew) - ion (q) interaction
  // gwEq : force due to induced charge (w) - Electric field due to ion (Eq) interaction
  // gEwEq : force due to Electric field due to induced charge (Ew) - Electric field due to ion (Eq) interaction
  
  if (dsphere.POLARIZED)
  {
  
    // declarations (necessary beforehand for parallel implementation)
    long double gwq, gww_wEw_EwEw, gEwq, gwEq, gwEq_EwEq;
    unsigned int kloop, l1, i1;
    vector<long double> innerg3(s.size(),0.0);
    vector<long double> innerg4(s.size(),0.0);
    
    // some pre-summations (Green's function, gradient of Green's function, gEwq, gwEq)
    #pragma omp parallel default(shared) private(kloop, l1, i1, gEwq, gwEq, gwq, gwEq_EwEq, gww_wEw_EwEw) num_threads(THREADSIZE)
    {
      #pragma omp for schedule(dynamic,CHUNKSIZE)
      for (kloop = 0; kloop < s.size(); kloop++)
      {
	for (i1 = 0; i1 < ion.size(); i1++)
	{
	  s[kloop].Gion[i1] = (1.0/((s[kloop].posvec - ion[i1].posvec).GetMagnitude()));	// push_back avoided (with new code in main)
	  s[kloop].gradGion[i1] = (Grad(s[kloop].posvec, ion[i1].posvec));		// push_back avoided similarly
	}
	
	gEwq = 0;
	for (i1 = 0; i1 < ion.size(); i1++)
	  gEwq += s[kloop].Gion[i1] * (ion[i1].q / ion[i1].epsilon);
	innerg3[kloop] = gEwq;
      
	gwEq = 0;
	for (i1 = 0; i1 < ion.size(); i1++) 
	  gwEq += (s[kloop].normalvec * s[kloop].gradGion[i1]) * (ion[i1].q/ion[i1].epsilon);
	innerg4[kloop] = gwEq;
      }
      
      // calculate force
      #pragma omp for schedule(dynamic,CHUNKSIZE) nowait
      for (kloop = 0; kloop < s.size(); kloop++)
      {
	gwq = 0;
	for (i1 = 0; i1 < ion.size(); i1++) 
	  gwq += (-1.0) * (0.5 - 0.5 * dsphere.em / ion[i1].epsilon) * ion[i1].q * s[kloop].Gion[i1];
	
	gww_wEw_EwEw = 0;
	for (l1 = 0; l1 < s.size(); l1++) 
	  gww_wEw_EwEw +=  ( (-1.0)*dsphere.em*(dsphere.em - 1)*s[kloop].Greens[l1] + 0.5*dsphere.ed*(2*dsphere.em - 1)*s[kloop].presumgwEw[l1] + (-1.0)*dsphere.ed*dsphere.ed*s[kloop].presumgEwEw[l1] ) * s[l1].w * s[l1].a;

	gEwq = 0;
	for (l1 = 0; l1 < s.size(); l1++)
	  gEwq += (-1.0) * 0.5 * dsphere.ed * s[kloop].ndotGradGreens[l1] * innerg3[l1] * s[l1].a;
	
	gwEq_EwEq = 0;
	for (l1 = 0; l1 < s.size(); l1++)
	  gwEq_EwEq += ( 0.5 * dsphere.ed * (2*dsphere.em - 1) * s[kloop].Greens[l1] + (-1.0) * dsphere.ed * dsphere.ed * s[kloop].presumgEwEq[l1] ) * innerg4[l1] * s[l1].a;
	
	s[kloop].fw = gwq + gww_wEw_EwEw + gEwq + gwEq_EwEq;
      }
    }
    
    // force on the fake degrees of freedom
    for (unsigned int k = 0; k < s.size(); k++)
      s[k].fw = s[k].a * s[k].fw * scalefactor;

    innerg3.clear();
    innerg4.clear();
  }
  
  else //if not POLARIZED, that is no induced charges due to difference in dielectric properties of NP and its environment
  {
  // force on the fake degrees of freedom
  for (unsigned int k = 0; k < s.size(); k++)
    s[k].fw = 0.0;
  }
  
  return;  
}
