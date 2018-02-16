// This file contains the routine that computes the force
// on the induced charge at vertex k and the force on the particle i
// for all k and i

#include "forces.h"

// Total Force on all degrees of freedom
void for_cpmd_calculate_force(vector<VERTEX> &s, vector<PARTICLE> &ion, INTERFACE &nanoparticle) {

    // force calculation for fake degrees of freedom
    // gwq : force due to induced charge (w) - ion (q) interaction
    // gww : force due to induced charge (w) - induced charge (w) interaction
    // gwEw : force due to induced charge (w) - Electric field due to induced charge (Ew) interaction
    // gEwEw : force due to Electric field due to induced charge (Ew) - Electric field due to induced charge (Ew) interaction
    // gEwq : force due to electric field of induced charge (Ew) - ion (q) interaction
    // gwEq : force due to induced charge (w) - Electric field due to ion (Eq) interaction
    // gEwEq : force due to Electric field due to induced charge (Ew) - Electric field due to ion (Eq) interaction

    // force calculation for ions (electrostatic)
    // hqqc : force due to central charge qc
    // hqq : force due to ion (q) - ion (q) interaction
    // hqw : force due to ion (q) - induced charge (w) interaction
    // part of hqEq : force due to ion (q) - electric field of ion (Eq) interaction
    // hqEw : force due to ion (q) - electric field of induced charge (Ew) interaction
    // other part of hqEq : force due to ion (q) - electric field of ion (Eq) interaction
    // hEqw : force due to electric field of ion (Eq) - induced charge (w) interaction
    // hEqEq : force due to electric field of ion (Eq) - electric field  of ion (Eq) interaction
    // hEqEw : force due to electric field of ion (Eq) - electric field of induced degree (Ew) interaction

    if (nanoparticle.POLARIZED) {
        // declarations (necessary beforehand for parallel implementation)
        long double gwq, gww_wEw_EwEw, gEwq, gwEq, gwEq_EwEq;
        long double hqEw, hqEq, hEqw, hEqEq, hEqEw;
        unsigned int kloop, l1, i1;
        unsigned int iloop, j1;
        long double insum;
        vector<long double> innerg3(s.size(), 0.0);
        vector<long double> innerg4(s.size(), 0.0);
        vector<long double> saveinsum(s.size(), 0.0);        // this is for h force calculation
        vector<long double> innerh2(s.size(), 0.0);
        vector<long double> innerh4(s.size(), 0.0);
        VECTOR3D h0, h1, h2, h3;

        // parallel calculation of fake and real forces
#pragma omp parallel default(shared) private(kloop, iloop, l1, i1, j1, gEwq, gwEq, gwEq_EwEq, gwq, gww_wEw_EwEw, hqEw, hqEq, hEqw, hEqEq, hEqEw, h0, h1, h2, h3)
        {
            // inner loop calculations for fake forces and one inner loop for real force
#pragma omp for schedule(dynamic)
            for (kloop = 0; kloop < s.size(); kloop++) {
                for (i1 = 0; i1 < ion.size(); i1++) {
                    s[kloop].Gion[i1] = (1.0 / ((s[kloop].posvec -
                                                 ion[i1].posvec).GetMagnitude()));    // push_back avoided (with new code in main)
                    s[kloop].gradGion[i1] = (Grad(s[kloop].posvec,
                                                  ion[i1].posvec));        // push_back avoided similarly
                }

                gEwq = 0;
                for (i1 = 0; i1 < ion.size(); i1++)
                    gEwq += s[kloop].Gion[i1] * (ion[i1].q / ion[i1].epsilon);
                innerg3[kloop] = gEwq;

                gwEq = 0;
                for (i1 = 0; i1 < ion.size(); i1++)
                    gwEq += (s[kloop].normalvec * s[kloop].gradGion[i1]) * (ion[i1].q / ion[i1].epsilon);
                innerg4[kloop] = gwEq;

                insum = 0;
                for (i1 = 0; i1 < ion.size(); i1++)
                    insum = insum + (s[kloop].normalvec * s[kloop].gradGion[i1]) * (ion[i1].q / ion[i1].epsilon);
                saveinsum[kloop] = insum;
            }

            // continuing with inner loop calculations for force on real ions
#pragma omp for schedule(dynamic)
            for (kloop = 0; kloop < s.size(); kloop++) {
                hqEw = saveinsum[kloop];
                for (l1 = 0; l1 < s.size(); l1++)
                    hqEw = hqEw + s[kloop].ndotGradGreens[l1] * s[l1].w * s[l1].a;

                innerh2[kloop] = hqEw;

                hqEq = 0;
                for (i1 = 0; i1 < ion.size(); i1++)
                    hqEq = hqEq + (s[kloop].Gion[i1] * (ion[i1].q / ion[i1].epsilon));
                hqEq = hqEq * (-1.0 * 0.5 * nanoparticle.ed);

                hEqw = 0;
                for (l1 = 0; l1 < s.size(); l1++)
                    hEqw = hEqw + s[kloop].Greens[l1] * s[l1].w * s[l1].a;
                hEqw = hEqw * (-1.0 * (-0.5) * nanoparticle.ed * (2 * nanoparticle.em - 1));

                hEqEq = 0;
                for (l1 = 0; l1 < s.size(); l1++)
                    hEqEq = hEqEq + s[kloop].Greens[l1] * saveinsum[l1] * s[l1].a;
                hEqEq = hEqEq * (-1.0 * nanoparticle.ed * nanoparticle.ed);

                hEqEw = 0;
                for (l1 = 0; l1 < s.size(); l1++)
                    hEqEw = hEqEw + s[kloop].presumhEqEw[l1] * s[l1].w * s[l1].a;
                hEqEw = hEqEw * (-1.0 * nanoparticle.ed * nanoparticle.ed);

                innerh4[kloop] = (hqEq + hEqw + hEqEq + hEqEw);
            }

            // fake force computation
#pragma omp for schedule(dynamic) nowait
            for (kloop = 0; kloop < s.size(); kloop++) {
                gwq = 0;
                for (i1 = 0; i1 < ion.size(); i1++)
                    gwq += (-1.0) * (0.5 - 0.5 * nanoparticle.em / ion[i1].epsilon) * ion[i1].q * s[kloop].Gion[i1];

                gww_wEw_EwEw = 0;
                for (l1 = 0; l1 < s.size(); l1++)
                    gww_wEw_EwEw += ((-1.0) * nanoparticle.em * (nanoparticle.em - 1) * s[kloop].Greens[l1] +
                                     0.5 * nanoparticle.ed * (2 * nanoparticle.em - 1) * s[kloop].presumgwEw[l1] +
                                     (-1.0) * nanoparticle.ed * nanoparticle.ed * s[kloop].presumgEwEw[l1]) * s[l1].w *
                                    s[l1].a;

                gEwq = 0;
                for (l1 = 0; l1 < s.size(); l1++)
                    gEwq += (-1.0) * 0.5 * nanoparticle.ed * s[kloop].ndotGradGreens[l1] * innerg3[l1] * s[l1].a;

                gwEq_EwEq = 0;
                for (l1 = 0; l1 < s.size(); l1++)
                    gwEq_EwEq += (0.5 * nanoparticle.ed * (2 * nanoparticle.em - 1) * s[kloop].Greens[l1] +
                                  (-1.0) * nanoparticle.ed * nanoparticle.ed * s[kloop].presumgEwEq[l1]) * innerg4[l1] *
                                 s[l1].a;

                s[kloop].fw = gwq + gww_wEw_EwEw + gEwq + gwEq_EwEq;
            }

            // force calculation for real ions (this was in parallel with the previous for loop)
#pragma omp for schedule(dynamic) nowait
            for (iloop = 0; iloop < ion.size(); iloop++) {
                h0 = ((Grad(ion[iloop].posvec, nanoparticle.posvec)) ^
                      ((-1.0) * nanoparticle.bare_charge * ion[iloop].q * 1.0 / ion[iloop].epsilon));

                h1 = VECTOR3D(0, 0, 0);
                for (j1 = 0; j1 < ion.size(); j1++) {
                    if (j1 == iloop) continue;
                    h1 = h1 + ((Grad(ion[iloop].posvec, ion[j1].posvec)) ^
                               ((-0.5) * ion[iloop].q * ion[j1].q * (1 / ion[iloop].epsilon + 1 / ion[j1].epsilon)));
                }

                double aqw = -1.0 * ion[iloop].q * (0.5 - nanoparticle.em / (2.0 * ion[iloop].epsilon));
                double bqEqw = -1.0 * 0.5 * nanoparticle.ed * ion[iloop].q / ion[iloop].epsilon;
                h2 = VECTOR3D(0, 0, 0);
                for (l1 = 0; l1 < s.size(); l1++)
                    h2 = h2 + (((s[l1].gradGion[iloop]) ^ (-1.0)) ^ ((aqw * s[l1].w + bqEqw * innerh2[l1]) * s[l1].a));

                h3 = VECTOR3D(0, 0, 0);
                for (l1 = 0; l1 < s.size(); l1++)
                    h3 = h3 +
                         (GradndotGrad(s[l1].posvec, ion[iloop].posvec, s[l1].normalvec) ^ (innerh4[l1] * s[l1].a));
                h3 = (h3 ^ (ion[iloop].q / ion[iloop].epsilon));

                ion[iloop].forvec = (h0 + h1 + h2 + h3);
            }
        }

        innerg3.clear();
        innerg4.clear();
        innerh2.clear();
        innerh4.clear();

        // force on the fake degrees of freedom
        for (unsigned int k = 0; k < s.size(); k++)
            s[k].fw = s[k].a * s[k].fw * scalefactor;

        // force on the particles (electrostatic)
        for (unsigned int i = 0; i < ion.size(); i++)
            ion[i].forvec = ((ion[i].forvec) ^ (scalefactor));

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    } else //if not POLARIZED, that is no induced charges due to difference in dielectric properties of NP and its environment
    {
        unsigned int iloop, j1;
        VECTOR3D h0, h1, h2, h3;
        // parallel calculation of real forces (uniform case)
#pragma omp parallel default(shared) private(iloop, j1, h0, h1)
        {
#pragma omp for schedule(dynamic) nowait
            for (iloop = 0; iloop < ion.size(); iloop++) {
                h0 = ((Grad(ion[iloop].posvec, nanoparticle.posvec)) ^
                      ((-1.0) * nanoparticle.bare_charge * ion[iloop].q * 1.0 / ion[iloop].epsilon));

                h1 = VECTOR3D(0, 0, 0);
                for (j1 = 0; j1 < ion.size(); j1++) {
                    if (j1 == iloop) continue;
                    h1 = h1 + ((Grad(ion[iloop].posvec, ion[j1].posvec)) ^
                               ((-0.5) * ion[iloop].q * ion[j1].q * (1 / ion[iloop].epsilon + 1 / ion[j1].epsilon)));
                }
                ion[iloop].forvec = (h0 + h1);
            }
        }

        // force on the fake degrees of freedom
        for (unsigned int k = 0; k < s.size(); k++)
            s[k].fw = 0.0;

        // force on the particles (electrostatic)
        for (unsigned int i = 0; i < ion.size(); i++)
            ion[i].forvec = ((ion[i].forvec) ^ (scalefactor));
    }
    // Excluded volume interactions given by purely repulsive LJ

    // ion-sphere (when ion is outside)
    // make a dummy particle with the same diameter as the ion just beneath the interface
    vector<VECTOR3D> lj1;
    for (unsigned int i = 0; i < ion.size(); i++) {
        VECTOR3D fljcs = VECTOR3D(0, 0, 0);
        if (ion[i].posvec.GetMagnitude() < nanoparticle.radius) {
            lj1.push_back(fljcs);
            continue;
        }
        PARTICLE dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, nanoparticle.ein, ion[i].posvec ^ ((nanoparticle.radius -
                                                                                                   0.5 *
                                                                                                   ion[i].diameter) /
                                                                                                  ion[i].posvec.GetMagnitude()));
        VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
        double r = r_vec.GetMagnitude();
        double d = 0.5 * (ion[i].diameter + dummy.diameter);
        double elj = 1.0;
        if (r < dcut * d) {
            double r2 = r * r;
            double r6 = r2 * r2 * r2;
            double r12 = r6 * r6;
            double d2 = d * d;
            double d6 = d2 * d2 * d2;
            double d12 = d6 * d6;
            fljcs = r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
        }
        lj1.push_back(fljcs);
    }

    // ion-ion
    vector<VECTOR3D> lj2;
    for (unsigned int i = 0; i < ion.size(); i++) {
        VECTOR3D fljcc = VECTOR3D(0, 0, 0);
        for (unsigned int j = 0; j < ion.size(); j++) {
            if (j == i) continue;
            VECTOR3D r_vec = ion[i].posvec - ion[j].posvec;
            double r = r_vec.GetMagnitude();
            double d = 0.5 * (ion[i].diameter + ion[j].diameter);
            double elj = 1.0;
            if (r < dcut * d) {
                double r2 = r * r;
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d2 = d * d;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                fljcc = fljcc + (r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            } else
                fljcc = fljcc + VECTOR3D(0, 0, 0);
        }
        lj2.push_back(fljcc);
    }

    // ion-box
    // make a dummy particle with the same diameter as the ion just above the simulation box
    vector<VECTOR3D> lj3;
    for (unsigned int i = 0; i < ion.size(); i++) {
        VECTOR3D fljcb = VECTOR3D(0, 0, 0);
        PARTICLE dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, nanoparticle.eout, ion[i].posvec ^
                                                                                  ((nanoparticle.box_radius +
                                                                                    0.5 * ion[i].diameter) /
                                                                                   ion[i].posvec.GetMagnitude()));
        VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
        double r = r_vec.GetMagnitude();
        double d = 0.5 * (ion[i].diameter + dummy.diameter);
        double elj = 1.0;
        if (r < dcut * d) {
            double r2 = r * r;
            double r6 = r2 * r2 * r2;
            double r12 = r6 * r6;
            double d2 = d * d;
            double d6 = d2 * d2 * d2;
            double d12 = d6 * d6;
            fljcb = r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
        }
        lj3.push_back(fljcb);
    }

    // ion-sphere (works when ions are inside)
    // make a dummy particle with the same diameter as the ion just above the interface
    vector<VECTOR3D> lj4;
    for (unsigned int i = 0; i < ion.size(); i++) {
        VECTOR3D fljcs = VECTOR3D(0, 0, 0);
        if (ion[i].posvec.GetMagnitude() > nanoparticle.radius) {
            lj4.push_back(fljcs);
            continue;
        }
        PARTICLE dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, nanoparticle.eout, ion[i].posvec ^
                                                                                  ((nanoparticle.radius +
                                                                                    0.5 * ion[i].diameter) /
                                                                                   ion[i].posvec.GetMagnitude()));
        VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
        double r = r_vec.GetMagnitude();
        double d = 0.5 * (ion[i].diameter + dummy.diameter);
        double elj = 1.0;
        if (r < dcut * d) {
            double r2 = r * r;
            double r6 = r2 * r2 * r2;
            double r12 = r6 * r6;
            double d2 = d * d;
            double d6 = d2 * d2 * d2;
            double d12 = d6 * d6;
            fljcs = r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
        }
        lj4.push_back(fljcs);
    }

    // total force on the particle = the electrostatic force + the Lennard-Jones force
    for (unsigned int i = 0; i < ion.size(); i++)
        ion[i].forvec = ((ion[i].forvec) + lj1[i] + lj2[i] + lj3[i] + lj4[i]);

    return;
}
