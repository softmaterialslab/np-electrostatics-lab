// This routine sets up and performs precalculations

#include "precalculations.h"

void precalculate(vector<VERTEX> &s, NanoParticle *nanoParticle) {
    
    

    if (world.rank() == 0)
        cout << "Precalculations..." << endl;

    for (unsigned int k = 0; k < s.size(); k++) {
        for (unsigned int l = 0; l < s.size(); l++)
            s[k].Greens.push_back(G(s, k, l));
    }
    for (unsigned int k = 0; k < s.size(); k++) {
        for (unsigned int l = 0; l < s.size(); l++)
            s[k].ndotGradGreens.push_back(H(s, l, k, nanoParticle->radius));
    }
    //if (world.rank() == 0)
    //    cout << "Next, one inner loop computation for gEwEw" << endl;

    // this will be a serial computation
    vector<vector<long double> > inner(s.size() * s.size());
    vector<long double> row(s.size(), 0.0);
    unsigned int k, m, l, n;
    long double gwEw, gEwEq, gEwEw, fwEw, fEwEq, pre_gEwEw;

    for (l = 0; l < s.size(); l++) {
        for (n = 0; n < s.size(); n++) {
            pre_gEwEw = 0;
            for (m = 0; m < s.size(); m++)
                pre_gEwEw += G(s, l, m) * H(s, n, m, nanoParticle->radius) * s[m].a;
            row[n] = pre_gEwEw;
        }
        inner[l] = row;
    }
    //if (world.rank() == 0)
    //    cout << "Precalculate gwEw, gEwEq, gEwEw, fwEw, fEwEq in parallel" << endl;

    // parallel computation
#pragma omp parallel default(shared) private(k, m, l, gwEw, gEwEq, gEwEw, fwEw, fEwEq)
    {
#pragma omp for schedule(dynamic) nowait
        for (k = 0; k < s.size(); k++) {
            for (m = 0; m < s.size(); m++) {
                // gwEw
                gwEw = 0;
                for (l = 0; l < s.size(); l++)
                    gwEw += (H(s, k, l, nanoParticle->radius) * G(s, m, l) +
                             G(s, k, l) * H(s, m, l, nanoParticle->radius)) * s[l].a;
                s[k].presumgwEw[m] = gwEw;

                // gEwEq
                gEwEq = 0;
                for (l = 0; l < s.size(); l++)
                    gEwEq += (H(s, k, l, nanoParticle->radius) * G(s, m, l)) * s[l].a;
                s[k].presumgEwEq[m] = gEwEq;

                // gEwEw
                gEwEw = 0;
                for (l = 0; l < s.size(); l++)
                    gEwEw += H(s, k, l, nanoParticle->radius) * inner[l][m] *
                             s[l].a;        // NOTE the switch from n to m
                s[k].presumgEwEw[m] = gEwEw;

                // fwEw
                fwEw = 0;
                for (l = 0; l < s.size(); l++)
                    fwEw += G(s, k, l) * H(s, m, l, nanoParticle->radius) * s[l].a;
                s[k].presumfwEw[m] = fwEw;

                // fEwEq
                fEwEq = 0;
                for (l = 0; l < s.size(); l++)
                    fEwEq += G(s, k, l) * H(s, m, l, nanoParticle->radius) * s[l].a;
                s[k].presumfEwEq[m] = fEwEq;
            }
        }
    }

    // hEqEw is the same as fEwEq
    for (unsigned int k = 0; k < s.size(); k++) {
        for (unsigned int m = 0; m < s.size(); m++) {
            s[k].presumhEqEw[m] = s[k].presumfEwEq[m];
        }
    }
    inner.clear();

    return;
}

