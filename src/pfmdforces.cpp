// This file contains the routine that computes the force
// on the induced charge at vertex k and the force on the particle i
// for all k and i

#include "forces.h"

// Total Force on all degrees of freedom
void for_fmd_calculate_force(vector<VERTEX> &s, vector<PARTICLE> &ion, NanoParticle *nanoParticle) {

    // force calculation for fake degrees of freedom
    // gwq : force due to induced charge (w) - ion (q) interaction
    // gww : force due to induced charge (w) - induced charge (w) interaction
    // gwEw : force due to induced charge (w) - Electric field due to induced charge (Ew) interaction
    // gEwEw : force due to Electric field due to induced charge (Ew) - Electric field due to induced charge (Ew) interaction
    // gEwq : force due to electric field of induced charge (Ew) - ion (q) interaction
    // gwEq : force due to induced charge (w) - Electric field due to ion (Eq) interaction
    // gEwEq : force due to Electric field due to induced charge (Ew) - Electric field due to ion (Eq) interaction

    if (nanoParticle->POLARIZED) {

        // declarations (necessary beforehand for parallel implementation)
        long double gwq, gww_wEw_EwEw, gEwq, gwEq, gwEq_EwEq;
        unsigned int kloop, l1, i1;

        
        
        /////////////POLARIZED only MPI Message objects
        vector<long double> innerg3(sizFVecMesh, 0.0);
        vector<long double> innerg3Gather(s.size() + extraElementsMesh, 0.0);
        vector<long double> innerg4(sizFVecMesh, 0.0);
        vector<long double> innerg4Gather(s.size() + extraElementsMesh, 0.0);
        vector<long double> fw(sizFVecMesh, 0.0);
        vector<long double> fwGather(s.size() + extraElementsMesh, 0.0);

        // some pre-summations (Green's function, gradient of Green's function, gEwq, gwEq)

#pragma parallel omp for schedule(dynamic) default(shared) private(kloop, i1, gEwq, gwEq)
            for (kloop = lowerBoundMesh; kloop <= upperBoundMesh; kloop++) {
                for (i1 = 0; i1 < ion.size(); i1++) {
                    s[kloop].Gion[i1] = (1.0 / ((s[kloop].posvec -
                                                 ion[i1].posvec).GetMagnitude()));    // push_back avoided (with new code in main)
                    s[kloop].gradGion[i1] = (Grad(s[kloop].posvec,
                                                  ion[i1].posvec));        // push_back avoided similarly
                }

                gEwq = 0;
                for (i1 = 0; i1 < ion.size(); i1++)
                    gEwq += s[kloop].Gion[i1] * (ion[i1].q / ion[i1].epsilon);
                innerg3[kloop - lowerBoundMesh] = gEwq;

                gwEq = 0;
                for (i1 = 0; i1 < ion.size(); i1++)
                    gwEq += (s[kloop].normalvec * s[kloop].gradGion[i1]) * (ion[i1].q / ion[i1].epsilon);
                innerg4[kloop - lowerBoundMesh] = gwEq;
            }


        //innerg3,innerg4 broadcasting using all gather = gather + broadcast
        if (world.size() > 1) {

            all_gather(world, &innerg3[0], innerg3.size(), innerg3Gather);
            all_gather(world, &innerg4[0], innerg4.size(), innerg4Gather);

        } else {
            for (kloop = lowerBoundMesh; kloop <= upperBoundMesh; kloop++) {

                innerg3Gather[kloop] = innerg3[kloop - lowerBoundMesh];
                innerg4Gather[kloop] = innerg4[kloop - lowerBoundMesh];

            }
        }

            // calculate force
#pragma parallel omp for schedule(dynamic) default(shared) private(kloop, l1, gEwq, gwq, gwEq_EwEq, gww_wEw_EwEw)
            for (kloop = lowerBoundMesh; kloop <= upperBoundMesh; kloop++) {
                gwq = 0;
                for (l1 = 0; l1 < ion.size(); l1++)
                    gwq += (-1.0) * (0.5 - 0.5 * nanoParticle->em / ion[l1].epsilon) * ion[l1].q * s[kloop].Gion[l1];

                gww_wEw_EwEw = 0;
                for (l1 = 0; l1 < s.size(); l1++)
                    gww_wEw_EwEw += ((-1.0) * nanoParticle->em * (nanoParticle->em - 1) * s[kloop].Greens[l1] +
                                     0.5 * nanoParticle->ed * (2 * nanoParticle->em - 1) * s[kloop].presumgwEw[l1] +
                                     (-1.0) * nanoParticle->ed * nanoParticle->ed * s[kloop].presumgEwEw[l1]) * s[l1].w *
                                    s[l1].a;

                gEwq = 0;
                for (l1 = 0; l1 < s.size(); l1++)
                    gEwq += (-1.0) * 0.5 * nanoParticle->ed * s[kloop].ndotGradGreens[l1] * innerg3Gather[l1] * s[l1].a;

                gwEq_EwEq = 0;
                for (l1 = 0; l1 < s.size(); l1++)
                    gwEq_EwEq += (0.5 * nanoParticle->ed * (2 * nanoParticle->em - 1) * s[kloop].Greens[l1] +
                                  (-1.0) * nanoParticle->ed * nanoParticle->ed * s[kloop].presumgEwEq[l1]) * innerg4Gather[l1] *
                                 s[l1].a;

                fw[kloop - lowerBoundMesh] = gwq + gww_wEw_EwEw + gEwq + gwEq_EwEq;
            }

        //fw broadcasting using all gather = gather + broadcast
        if (world.size() > 1)
            all_gather(world, &fw[0], fw.size(), fwGather);
        else
            for (kloop = lowerBoundMesh; kloop <= upperBoundMesh; kloop++)
                fwGather[kloop] = fw[kloop - lowerBoundMesh];

        // force on the fake degrees of freedom
        for (unsigned int k = 0; k < s.size(); k++)
            s[k].fw = s[k].a * fwGather[k] * scalefactor;

        innerg3.clear();
        innerg3Gather.clear();
        innerg4.clear();
        innerg4Gather.clear();

    } else //if not POLARIZED, that is no induced charges due to difference in dielectric properties of NP and its environment
    {
        // force on the fake degrees of freedom
        for (unsigned int k = 0; k < s.size(); k++)
            s[k].fw = 0.0;
    }

    return;
}
