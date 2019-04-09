// This file contains the routine that computes the force
// on the induced charge at vertex k and the force on the particle i
// for all k and i

#include "forces.h"

// Total Force on all degrees of freedom
void
for_cpmd_calculate_force(vector<VERTEX> &s, vector<PARTICLE> &ion, NanoParticle *nanoParticle) {

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

    
    

    //Common MPI Message objects
    vector<VECTOR3D> forvec(sizFVecIons, VECTOR3D(0, 0, 0));
    vector<VECTOR3D> lj1(sizFVecIons, VECTOR3D(0, 0, 0));
    vector<VECTOR3D> lj2(sizFVecIons, VECTOR3D(0, 0, 0));
    vector<VECTOR3D> lj3(sizFVecIons, VECTOR3D(0, 0, 0));
    vector<VECTOR3D> lj4(sizFVecIons, VECTOR3D(0, 0, 0));
    vector<VECTOR3D> forvecGather(ion.size() + extraElementsIons, VECTOR3D(0, 0, 0));

    unsigned int iloop, j1;


    if (nanoParticle->POLARIZED) {
        // declarations (necessary beforehand for parallel implementation)
        long double gwq, gww_wEw_EwEw, gEwq, gwEq, gwEq_EwEq;
        long double hqEw, hqEq, hEqw, hEqEq, hEqEw;
        unsigned int kloop, l1, i1;

        long double insum;

        VECTOR3D h0, h1, h2, h3;

        /////////////POLARIZED only MPI Message objects
        vector<long double> saveinsum(sizFVecMesh, 0.0);
        vector<long double> saveinsumGather(s.size() + extraElementsMesh, 0.0);
        vector<long double> innerg3(sizFVecMesh, 0.0);
        vector<long double> innerg3Gather(s.size() + extraElementsMesh, 0.0);
        vector<long double> innerg4(sizFVecMesh, 0.0);
        vector<long double> innerg4Gather(s.size() + extraElementsMesh, 0.0);
        vector<long double> innerh2(sizFVecMesh, 0.0);
        vector<long double> innerh2Gather(s.size() + extraElementsMesh, 0.0);
        vector<long double> innerh4(sizFVecMesh, 0.0);
        vector<long double> innerh4Gather(s.size() + extraElementsMesh, 0.0);
        vector<long double> fw(sizFVecMesh, 0.0);
        vector<long double> fwGather(s.size() + extraElementsMesh, 0.0);
        //////////

        // parallel calculation of fake and real forces
        // inner loop calculations for fake forces and one inner loop for real force : V1
#pragma omp parallel for schedule(dynamic) private(kloop, i1)
        for (kloop = 0; kloop < s.size(); kloop++) {
            for (i1 = 0; i1 < ion.size(); i1++)
                s[kloop].gradGion[i1] = (Grad(s[kloop].posvec,
                                              ion[i1].posvec));        // push_back avoided similarly
        }

        // inner loop calculations for fake forces and one inner loop for real force: V2
#pragma omp parallel for schedule(dynamic) default(shared) private(kloop, i1, gEwq, gwEq, insum)
        for (kloop = lowerBoundMesh; kloop <= upperBoundMesh; kloop++) {

            for (i1 = 0; i1 < ion.size(); i1++)
                s[kloop].Gion[i1] = (1.0 / ((s[kloop].posvec -
                                             ion[i1].posvec).GetMagnitude()));    // push_back avoided (with new code in main)

            gEwq = 0;
            for (i1 = 0; i1 < ion.size(); i1++)
                gEwq += s[kloop].Gion[i1] * (ion[i1].q / ion[i1].epsilon);
            innerg3[kloop - lowerBoundMesh] = gEwq;

            gwEq = 0;
            for (i1 = 0; i1 < ion.size(); i1++)
                gwEq += (s[kloop].normalvec * s[kloop].gradGion[i1]) * (ion[i1].q / ion[i1].epsilon);
            innerg4[kloop - lowerBoundMesh] = gwEq;

            insum = 0;
            for (i1 = 0; i1 < ion.size(); i1++)
                insum = insum + (s[kloop].normalvec * s[kloop].gradGion[i1]) * (ion[i1].q / ion[i1].epsilon);

            saveinsum[kloop - lowerBoundMesh] = insum;
        }

        //saveinsum,innerg3,innerg4 broadcasting using all gather = gather + broadcast
        if (world.size() > 1) {

            //cout <<"This is proc : " << world.rank() << ", LowerBound"<<lowerBoundMesh<< ", UpperBound"<<upperBoundMesh<<endl;
            //cout  << world.rank() <<" : Size of the data trying to send : " << saveinsum.size() <<"saveinsumGather size : " << saveinsumGather.size() <<endl;
            all_gather(world, &saveinsum[0], saveinsum.size(), saveinsumGather);
            all_gather(world, &innerg3[0], innerg3.size(), innerg3Gather);
            all_gather(world, &innerg4[0], innerg4.size(), innerg4Gather);

            //cout <<"Done proc : " << world.rank() <<endl;


        } else {
            for (kloop = lowerBoundMesh; kloop <= upperBoundMesh; kloop++) {
                saveinsumGather[kloop] = saveinsum[kloop - lowerBoundMesh];
                innerg3Gather[kloop] = innerg3[kloop - lowerBoundMesh];
                innerg4Gather[kloop] = innerg4[kloop - lowerBoundMesh];
            }
        }

        // continuing with inner loop calculations for force on real ions
#pragma omp parallel for schedule(dynamic) default(shared) private(kloop, l1, hqEw, hqEq, hEqw, hEqEq, hEqEw)
        for (kloop = lowerBoundMesh; kloop <= upperBoundMesh; kloop++) {
            hqEw = saveinsumGather[kloop];
            for (l1 = 0; l1 < s.size(); l1++)
                hqEw = hqEw + s[kloop].ndotGradGreens[l1] * s[l1].w * s[l1].a;

            innerh2[kloop - lowerBoundMesh] = hqEw;

            hqEq = 0;
            for (l1 = 0; l1 < ion.size(); l1++)
                hqEq = hqEq + (s[kloop].Gion[l1] * (ion[l1].q / ion[l1].epsilon));
            hqEq = hqEq * (-1.0 * 0.5 * nanoParticle->ed);

            hEqw = 0;
            for (l1 = 0; l1 < s.size(); l1++)
                hEqw = hEqw + s[kloop].Greens[l1] * s[l1].w * s[l1].a;
            hEqw = hEqw * (-1.0 * (-0.5) * nanoParticle->ed * (2 * nanoParticle->em - 1));

            hEqEq = 0;
            for (l1 = 0; l1 < s.size(); l1++)
                hEqEq = hEqEq + s[kloop].Greens[l1] * saveinsumGather[l1] * s[l1].a;
            hEqEq = hEqEq * (-1.0 * nanoParticle->ed * nanoParticle->ed);

            hEqEw = 0;
            for (l1 = 0; l1 < s.size(); l1++)
                hEqEw = hEqEw + s[kloop].presumhEqEw[l1] * s[l1].w * s[l1].a;
            hEqEw = hEqEw * (-1.0 * nanoParticle->ed * nanoParticle->ed);

            innerh4[kloop - lowerBoundMesh] = (hqEq + hEqw + hEqEq + hEqEw);
        }


        //innerh2,innerh4 broadcasting using all gather = gather + broadcast
        if (world.size() > 1) {

            all_gather(world, &innerh2[0], innerh2.size(), innerh2Gather);
            all_gather(world, &innerh4[0], innerh4.size(), innerh4Gather);


        } else {
            for (kloop = lowerBoundMesh; kloop <= upperBoundMesh; kloop++) {

                innerh2Gather[kloop] = innerh2[kloop - lowerBoundMesh];
                innerh4Gather[kloop] = innerh4[kloop - lowerBoundMesh];
            }
        }

        // fake force computation
#pragma omp parallel for schedule(dynamic) default(shared) private(kloop, l1, gEwq, gwEq_EwEq, gwq, gww_wEw_EwEw)
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
                              (-1.0) * nanoParticle->ed * nanoParticle->ed * s[kloop].presumgEwEq[l1]) *
                             innerg4Gather[l1] *
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

        // force calculation for real ions (this was in parallel with the previous for loop)
#pragma omp parallel for schedule(dynamic) default(shared) private(iloop, l1, h0, h1, h2, h3)
        for (iloop = lowerBoundIons; iloop <= upperBoundIons; iloop++) {
            //h0 = ((Grad(ion[iloop].posvec, nanoParticle->posvec)) ^
            //      ((-1.0) * nanoParticle->bare_charge * ion[iloop].q * 1.0 / ion[iloop].epsilon));

            h0 = VECTOR3D(0, 0, 0);
            for (int k = 0; k < s.size(); k++)
                h0 = h0 + ((Grad(ion[iloop].posvec, s[k].posvec)) ^
                           ((-1.0) * s[k].realQ * ion[iloop].q * 1.0 / ion[iloop].epsilon));


            h1 = VECTOR3D(0, 0, 0);
            for (l1 = 0; l1 < ion.size(); l1++) {
                if (l1 == iloop) continue;
                h1 = h1 + ((Grad(ion[iloop].posvec, ion[l1].posvec)) ^
                           ((-0.5) * ion[iloop].q * ion[l1].q * (1 / ion[iloop].epsilon + 1 / ion[l1].epsilon)));
            }

            double aqw = -1.0 * ion[iloop].q * (0.5 - nanoParticle->em / (2.0 * ion[iloop].epsilon));
            double bqEqw = -1.0 * 0.5 * nanoParticle->ed * ion[iloop].q / ion[iloop].epsilon;
            h2 = VECTOR3D(0, 0, 0);
            for (l1 = 0; l1 < s.size(); l1++)
                h2 = h2 +
                     (((s[l1].gradGion[iloop]) ^ (-1.0)) ^ ((aqw * s[l1].w + bqEqw * innerh2Gather[l1]) * s[l1].a));


            h3 = VECTOR3D(0, 0, 0);
            for (l1 = 0; l1 < s.size(); l1++)
                h3 = h3 +
                     (GradndotGrad(s[l1].posvec, ion[iloop].posvec, s[l1].normalvec) ^ (innerh4Gather[l1] * s[l1].a));
            h3 = (h3 ^ (ion[iloop].q / ion[iloop].epsilon));

            forvec[iloop - lowerBoundIons] = (h0 + h1 + h2 + h3);

        }

        saveinsum.clear();
        saveinsumGather.clear();
        innerg3.clear();
        innerg3Gather.clear();
        innerg4.clear();
        innerg4Gather.clear();
        innerh2.clear();
        innerh2Gather.clear();
        innerh4.clear();
        innerh4Gather.clear();
        fw.clear();
        fwGather.clear();



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    } else //if not POLARIZED, that is no induced charges due to difference in dielectric properties of NP and its environment
    {

        VECTOR3D h0, h1, h2, h3;
        // parallel calculation of real forces (uniform case)

#pragma omp parallel for schedule(dynamic) default(shared) private(iloop, j1, h0, h1)
        for (iloop = lowerBoundIons; iloop <= upperBoundIons; iloop++) {
            //h0 = ((Grad(ion[iloop].posvec, nanoParticle->posvec)) ^
             //     ((-1.0) * nanoParticle->bare_charge * ion[iloop].q * 1.0 / ion[iloop].epsilon));

            h0 = VECTOR3D(0, 0, 0);
            for (int k = 0; k < s.size(); k++)
            h0 = h0 + ((Grad(ion[iloop].posvec, s[k].posvec)) ^
                  ((-1.0) * s[k].realQ * ion[iloop].q * 1.0 / ion[iloop].epsilon));

            //if (iloop == 0)
            //cout << iloop << " : " << h0.GetMagnitude() << endl;

            h1 = VECTOR3D(0, 0, 0);
            for (j1 = 0; j1 < ion.size(); j1++) {
                if (j1 == iloop) continue;
                h1 = h1 + ((Grad(ion[iloop].posvec, ion[j1].posvec)) ^
                           ((-0.5) * ion[iloop].q * ion[j1].q * (1 / ion[iloop].epsilon + 1 / ion[j1].epsilon)));
            }
            forvec[iloop - lowerBoundIons] = (h0 + h1);
        }


        // force on the fake degrees of freedom
        for (unsigned int k = 0; k < s.size(); k++)
            s[k].fw = 0.0;

    }

    /////////////////////// Not POLARIZED over

    // Excluded volume interactions given by purely repulsive LJ

    // ion-sphere (when ion is outside)
    // make a dummy particle with the same diameter as the ion just beneath the interface

    for (iloop = lowerBoundIons; iloop <= upperBoundIons; iloop++) {
        VECTOR3D fljcs = VECTOR3D(0, 0, 0);
        if (ion[iloop].posvec.GetMagnitude() < nanoParticle->radius)
            continue;

        PARTICLE dummy = PARTICLE(0, ion[iloop].diameter, 0, 0, 0, nanoParticle->ein,
                                  ion[iloop].posvec ^ ((nanoParticle->radius -
                                                        0.5 *
                                                        ion[iloop].diameter) /
                                                       ion[iloop].posvec.GetMagnitude()));
        VECTOR3D r_vec = ion[iloop].posvec - dummy.posvec;
        double r = r_vec.GetMagnitude();
        double d = 0.5 * (ion[iloop].diameter + dummy.diameter);
        double elj = 1.0;
        if (r < dcut * d) {
            double r2 = r * r;
            double r6 = r2 * r2 * r2;
            double r12 = r6 * r6;
            double d2 = d * d;
            double d6 = d2 * d2 * d2;
            double d12 = d6 * d6;
            fljcs = r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
            lj1[iloop - lowerBoundIons] = fljcs;
        }

    }

    // ion-ion
#pragma omp parallel for schedule(dynamic) default(shared) private(iloop, j1)
    for (iloop = lowerBoundIons; iloop <= upperBoundIons; iloop++) {
        VECTOR3D fljcc = VECTOR3D(0, 0, 0);
        for (j1 = 0; j1 < ion.size(); j1++) {
            if (j1 == iloop) continue;
            VECTOR3D r_vec = ion[iloop].posvec - ion[j1].posvec;
            double r = r_vec.GetMagnitude();
            double d = 0.5 * (ion[iloop].diameter + ion[j1].diameter);
            double elj = 1.0;
            if (r < dcut * d) {
                double r2 = r * r;
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d2 = d * d;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                fljcc = fljcc + (r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            }
        }
        lj2[iloop - lowerBoundIons] = fljcc;
    }

    // ion-box
    // make a dummy particle with the same diameter as the ion just above the simulation box
    for (iloop = lowerBoundIons; iloop <= upperBoundIons; iloop++) {
        PARTICLE dummy = PARTICLE(0, ion[iloop].diameter, 0, 0, 0, nanoParticle->eout, ion[iloop].posvec ^
                                                                                      ((nanoParticle->box_radius +
                                                                                        0.5 * ion[iloop].diameter) /
                                                                                       ion[iloop].posvec.GetMagnitude()));
        VECTOR3D r_vec = ion[iloop].posvec - dummy.posvec;
        double r = r_vec.GetMagnitude();
        double d = 0.5 * (ion[iloop].diameter + dummy.diameter);
        double elj = 1.0;
        if (r < dcut * d) {
            double r2 = r * r;
            double r6 = r2 * r2 * r2;
            double r12 = r6 * r6;
            double d2 = d * d;
            double d6 = d2 * d2 * d2;
            double d12 = d6 * d6;
            lj3[iloop - lowerBoundIons] = r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
        }

    }

    // ion-sphere (works when ions are inside)
    // make a dummy particle with the same diameter as the ion just above the interface
    for (iloop = lowerBoundIons; iloop <= upperBoundIons; iloop++) {

        if (ion[iloop].posvec.GetMagnitude() > nanoParticle->radius)
            continue;

        PARTICLE dummy = PARTICLE(0, ion[iloop].diameter, 0, 0, 0, nanoParticle->eout, ion[iloop].posvec ^
                                                                                      ((nanoParticle->radius +
                                                                                        0.5 * ion[iloop].diameter) /
                                                                                       ion[iloop].posvec.GetMagnitude()));
        VECTOR3D r_vec = ion[iloop].posvec - dummy.posvec;
        double r = r_vec.GetMagnitude();
        double d = 0.5 * (ion[iloop].diameter + dummy.diameter);
        double elj = 1.0;
        if (r < dcut * d) {
            double r2 = r * r;
            double r6 = r2 * r2 * r2;
            double r12 = r6 * r6;
            double d2 = d * d;
            double d6 = d2 * d2 * d2;
            double d12 = d6 * d6;
            lj4[iloop - lowerBoundIons] = r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
        }

    }

    // Total force on the particle = the electrostatic force + the Lennard-Jones force
    for (iloop = 0; iloop < forvec.size(); iloop++)
        forvec[iloop] = ((forvec[iloop]) ^ (scalefactor)) + lj1[iloop] + lj2[iloop] + lj3[iloop] + lj4[iloop];

    //forvec broadcasting using all gather = gather + broadcast
    if (world.size() > 1)
        all_gather(world, &forvec[0], forvec.size(), forvecGather);
    else
        for (iloop = lowerBoundIons; iloop <= upperBoundIons; iloop++)
            forvecGather[iloop] = forvec[iloop - lowerBoundIons];

    // force on the particles (electrostatic)
    for (iloop = 0; iloop < ion.size(); iloop++)
        ion[iloop].forvec = forvecGather[iloop];

    forvec.clear();
    forvecGather.clear();
    lj1.clear();
    lj2.clear();
    lj3.clear();
    lj4.clear();

    return;
}
