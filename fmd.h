#ifndef _FMD_H
#define _FMD_H

#include "vertex.h"
#include "particle.h"
#include "interface.h"
#include "control.h"
#include "forces.h"
#include "energies.h"
#include "functions.h"

void fmd(vector<VERTEX>&, vector<PARTICLE>&, INTERFACE&, PARTICLE&, CONTROL&, CONTROL&);

#endif

