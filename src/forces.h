// This is a header file containing the forces routines on induced charges 
// and particles respectively

#ifndef _FORCES_H
#define _FORCES_H

#include "vertex.h"
#include "particle.h"
#include "functions.h"

void for_fmd_calculate_force(vector<VERTEX> &, vector<PARTICLE> &, NanoParticle *);

void
for_cpmd_calculate_force(vector<VERTEX> &, vector<PARTICLE> &, NanoParticle *);

#endif 
