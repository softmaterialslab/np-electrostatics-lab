// This file contains the routine to compute the total potential energy of the system
// Electrostatic and Excluded volume contributions

#include "energies.h"

// Potential energy
double energy_functional(vector<VERTEX> &s, vector<PARTICLE> &ion, NanoParticle *nanoParticle) {
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

    //Common MPI Message objects
    vector<double> ion_energy(sizFVecIons, 0.0);
    vector<double> lj1(sizFVecIons, 0.0);
    vector<double> lj2(sizFVecIons, 0.0);
    vector<double> lj3(sizFVecIons, 0.0);
    vector<double> lj4(sizFVecIons, 0.0);

    double potential,totalPotential;
    unsigned int i, j;

    if (nanoParticle->POLARIZED) {

        /////////////POLARIZED only MPI Message objects
        vector<long double> saveinner1(sizFVecMesh, 0.0);
        vector<long double> saveinner1Gather(s.size() + extraElementsMesh, 0.0);
        vector<long double> inner2(sizFVecMesh, 0.0);
        vector<long double> inner2Gather(s.size() + extraElementsMesh, 0.0);
        vector<long double> inner3(sizFVecMesh, 0.0);
        vector<long double> inner3Gather(s.size() + extraElementsMesh, 0.0);
        vector<long double> inner4(sizFVecMesh, 0.0);
        vector<long double> inner4Gather(s.size() + extraElementsMesh, 0.0);
        vector<long double> inner5(sizFVecMesh, 0.0);
        vector<long double> inner5Gather(s.size() + extraElementsMesh, 0.0);
        vector<double> ind_energy(sizFVecMesh, 0.0);


        unsigned int k, l;
        double fqq, fwq, fqEq_qEw, fwEq_EqEq_EwEq;
        double insum;
        double ind_ind;
        double ion_ion = 0;



#pragma omp parallel for schedule(dynamic) default(shared) private(k, l, insum)
        for (k = lowerBoundMesh; k <= upperBoundMesh; k++) {
            insum = 0;
            for (unsigned int l = 0; l < ion.size(); l++)
                insum += (s[k].normalvec * Grad(s[k].posvec, ion[l].posvec)) * ion[l].q / ion[l].epsilon;
            saveinner1[k - lowerBoundMesh] = insum;
        }

        //saveinner1 broadcasting using all gather = gather + broadcast
        if (world.size() > 1)
            all_gather(world, &saveinner1[0], saveinner1.size(), saveinner1Gather);
        else
            for (k = lowerBoundMesh; k <= upperBoundMesh; k++)
                saveinner1Gather[k] = saveinner1[k - lowerBoundMesh];


#pragma omp parallel for schedule(dynamic) default(shared) private(k, l, insum)
        for (k = lowerBoundMesh; k <= upperBoundMesh; k++) {
            insum = 0;
            for (l = 0; l < s.size(); l++)
                insum += s[k].ndotGradGreens[l] * s[l].w * s[l].a;
            inner2[k - lowerBoundMesh] = insum;

            insum = 0;
            for (l = 0; l < s.size(); l++)
                insum += s[k].Greens[l] * s[l].w * s[l].a;
            inner3[k - lowerBoundMesh] = insum;

            insum = 0;
            for (l = 0; l < s.size(); l++)
                insum += s[k].Greens[l] * saveinner1Gather[l] * s[l].a;
            inner4[k - lowerBoundMesh] = insum;

            insum = 0;
            for (l = 0; l < s.size(); l++)
                insum += s[k].presumfEwEq[l] * s[l].w * s[l].a;
            inner5[k - lowerBoundMesh] = insum;
        }
        //inner2,inner3,inner4 broadcasting using all gather = gather + broadcast
        if (world.size() > 1) {

            all_gather(world, &inner2[0], inner2.size(), inner2Gather);
            all_gather(world, &inner3[0], inner3.size(), inner3Gather);
            all_gather(world, &inner4[0], inner4.size(), inner4Gather);
            all_gather(world, &inner5[0], inner5.size(), inner5Gather);

        } else {
            for (k = lowerBoundMesh; k <= upperBoundMesh; k++) {
                inner2Gather[k] = inner2[k - lowerBoundMesh];
                inner3Gather[k] = inner3[k - lowerBoundMesh];
                inner4Gather[k] = inner4[k - lowerBoundMesh];
                inner5Gather[k] = inner5[k - lowerBoundMesh];
            }
        }

#pragma omp parallel for schedule(dynamic) default(shared) private(k, l, ind_ind)
        for (k = lowerBoundMesh; k <= upperBoundMesh; k++) {
            ind_ind = 0;
            for (l = 0; l < s.size(); l++)
                ind_ind += s[k].w * s[k].a *
                           ((-0.5) * nanoParticle->ed * (2 * nanoParticle->em - 1) * s[k].presumfwEw[l] +
                            0.5 * nanoParticle->em * (nanoParticle->em - 1) * s[k].Greens[l] +
                            0.5 * nanoParticle->ed * nanoParticle->ed * s[k].presumgEwEw[l]) * s[l].w * s[l].a;
            ind_energy[k - lowerBoundMesh] = ind_ind;
        }

#pragma omp parallel for schedule(dynamic) default(shared) private(k, i, insum, fqq, fwq, fqEq_qEw, fwEq_EqEq_EwEq)
        for (i = lowerBoundIons; i <= upperBoundIons; i++) {
            fqq = 0;
            for (k = 0; k < ion.size(); k++) {
                if (i == k) continue;
                fqq += 0.5 * ion[i].q * ion[k].q * (1.0 / ion[i].epsilon) /
                       ((ion[i].posvec - ion[k].posvec).GetMagnitude());
            }

            fwq = 0;
            for (k = 0; k < s.size(); k++)
                fwq += ion[i].q * (0.5 - 0.5 * nanoParticle->em / ion[i].epsilon) *
                       (1 / ((s[k].posvec - ion[i].posvec).GetMagnitude())) * s[k].w * s[k].a;

            insum = 0;
            for (k = 0; k < s.size(); k++)
                insum += (1.0 / ((ion[i].posvec - s[k].posvec).GetMagnitude())) * (saveinner1Gather[k] + inner2Gather[k]) *
                         s[k].a;
            fqEq_qEw = 0.5 * nanoParticle->ed * ion[i].q * insum / ion[i].epsilon;

            insum = 0;
            for (k = 0; k < s.size(); k++)
                insum += (s[k].normalvec * Grad(s[k].posvec, ion[i].posvec)) *
                         ((-1) * 0.5 * nanoParticle->ed * (2 * nanoParticle->em - 1) * inner3Gather[k] +
                          0.5 * nanoParticle->ed * nanoParticle->ed * inner4Gather[k] +
                          nanoParticle->ed * nanoParticle->ed * inner5Gather[k]) * s[k].a;
            fwEq_EqEq_EwEq = (ion[i].q / ion[i].epsilon) * insum;

            insum = 0;
            for (k = 0; k < s.size(); k++)
                insum += s[k].realQ * ion[i].q * (1.0 / ion[i].epsilon) /
                ((ion[i].posvec - s[k].posvec).GetMagnitude());
            /*
            ion_energy[i-lowerBoundIons] = fqq + fwq + fqEq_qEw + fwEq_EqEq_EwEq +
                            nanoParticle->bare_charge * ion[i].q * (1.0 / ion[i].epsilon) /
                            ((ion[i].posvec - nanoParticle->posvec).GetMagnitude());*/

            ion_energy[i-lowerBoundIons] = fqq + fwq + fqEq_qEw + fwEq_EqEq_EwEq + insum;
        }


        // ion-ion + ion-induced charge energy
        ion_ion = 0;
        for (i = 0; i < ion_energy.size(); i++)
            ion_ion = ion_ion + ion_energy[i];
        ind_ind = 0;
        for (k = 0; k < ind_energy.size(); k++)
            ind_ind = ind_ind + ind_energy[k];

        // electrostatic potential energy
        potential = (ion_ion + ind_ind) * scalefactor;
    } else    // if not POLARIZED
    {
        double fqq;

#pragma omp parallel for schedule(dynamic) default(shared) private(i, j, fqq)
            for (i = lowerBoundIons; i <= upperBoundIons; i++) {
                fqq = 0;
                for (j = 0; j < ion.size(); j++) {
                    if (i == j) continue;
                    fqq += 0.5 * ion[i].q * ion[j].q * (1.0 / ion[i].epsilon) /
                           ((ion[i].posvec - ion[j].posvec).GetMagnitude());
                }

                double insum = 0;
                for (int k = 0; k < s.size(); k++)
                    insum += s[k].realQ * ion[i].q * (1.0 / ion[i].epsilon) /
                             ((ion[i].posvec - s[k].posvec).GetMagnitude());

                //ion_energy[i-lowerBoundIons] = fqq + nanoParticle->bare_charge * ion[i].q * (1.0 / ion[i].epsilon) /
                //                      ((ion[i].posvec - nanoParticle->posvec).GetMagnitude());

                ion_energy[i-lowerBoundIons] = fqq + insum;
            }


        double ion_ion = 0;
        for (unsigned int i = 0; i < ion_energy.size(); i++)
            ion_ion = ion_ion + ion_energy[i];

        // electrostatic potential energy
        potential = (ion_ion) * scalefactor;
    }

    // Excluded volume interaction energy given by purely repulsive LJ

    // ion-sphere (ions outsdie)
    // make a dummy particle with the same diameter as the ion just beneath the interface

    for (i = lowerBoundIons; i <= upperBoundIons; i++) {
        if (ion[i].posvec.GetMagnitude() < nanoParticle->radius)
            continue;

        PARTICLE dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, nanoParticle->ein, ion[i].posvec ^ ((nanoParticle->radius -
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
            double d2 = d * d;
            double d6 = d2 * d2 * d2;
            lj1[i-lowerBoundIons] = 4 * elj * (d6 / r6) * ((d6 / r6) - 1) + elj;
        }
    }

    // ion-ion
#pragma omp parallel for schedule(dynamic) default(shared) private(i, j)
    for (i = lowerBoundIons; i <= upperBoundIons; i++) {
        double uljcc = 0.0;
        for (j = 0; j < ion.size(); j++) {
            if (j == i) continue;
            VECTOR3D r_vec = ion[i].posvec - ion[j].posvec;
            double r = r_vec.GetMagnitude();
            double d = 0.5 * (ion[i].diameter + ion[j].diameter);
            double elj = 1.0;
            if (r < dcut * d) {
                double r2 = r * r;
                double r6 = r2 * r2 * r2;
                double d2 = d * d;
                double d6 = d2 * d2 * d2;
                uljcc = uljcc + 4 * elj * (d6 / r6) * ((d6 / r6) - 1) + elj;
            }
        }
        lj2[i-lowerBoundIons]=uljcc;
    }

    // ion-box
    // make a dummy particle with the same diameter as the ion just above the simulation box

    for (i = lowerBoundIons; i <= upperBoundIons; i++) {
        PARTICLE dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, nanoParticle->eout, ion[i].posvec ^
                                                                                  ((nanoParticle->box_radius +
                                                                                    0.5 * ion[i].diameter) /
                                                                                   ion[i].posvec.GetMagnitude()));
        VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
        double r = r_vec.GetMagnitude();
        double d = 0.5 * (ion[i].diameter + dummy.diameter);
        double elj = 1.0;
        if (r < dcut * d) {
            double r2 = r * r;
            double r6 = r2 * r2 * r2;
            double d2 = d * d;
            double d6 = d2 * d2 * d2;
            lj3[i-lowerBoundIons] = 4 * elj * (d6 / r6) * ((d6 / r6) - 1) + elj;
        }

    }

    // ion-sphere (ions inside)
    // make a dummy particle with the same diameter as the ion just above the interface

    for (i = lowerBoundIons; i <= upperBoundIons; i++) {
        if (ion[i].posvec.GetMagnitude() > nanoParticle->radius)
            continue;

        PARTICLE dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, nanoParticle->eout, ion[i].posvec ^
                                                                                  ((nanoParticle->radius +
                                                                                    0.5 * ion[i].diameter) /
                                                                                   ion[i].posvec.GetMagnitude()));
        VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
        double r = r_vec.GetMagnitude();
        double d = 0.5 * (ion[i].diameter + dummy.diameter);
        double elj = 1.0;
        if (r < dcut * d) {
            double r2 = r * r;
            double r6 = r2 * r2 * r2;
            double d2 = d * d;
            double d6 = d2 * d2 * d2;
            lj4[i-lowerBoundIons] = 4 * elj * (d6 / r6) * ((d6 / r6) - 1) + elj;
        }
    }

    double lj1energy = 0,lj2energy = 0,lj3energy = 0,lj4energy = 0;

    for (unsigned int i = 0; i < lj1.size(); i++) {
        lj1energy += lj1[i];
        lj2energy += lj2[i];
        lj3energy += lj3[i];
        lj4energy += lj4[i];
    }

    lj2energy = 0.5 * lj2energy;

    potential = potential + lj1energy + lj2energy + lj3energy + lj4energy;

    //MPI Operations
    if (world.size() > 1) {
        //broadcasting using all_reduce = reduce + broadcast
        //void all_reduce(const communicator & comm, const T & in_value, T & out_value, Op op);
        all_reduce(world, potential, totalPotential,std::plus<double>());
    }else{

        totalPotential=potential;
    }

    return totalPotential;
}

