// This is particle class

#include "NanoParticle.h"

class NanoParticleDisk : public NanoParticle {

private:
    // members
    vector<vector<BinRing> > bin_pos;
    vector<vector<BinRing> > bin_neg;
    double bin_width_R;
    double bin_width_Z;
    vector<PARTICLE> *ion;
    vector<vector<double> > density_pos;
    vector<vector<double> > density_neg;
    int cpmdstep;
    double density_profile_samples;
    vector<vector<double> > mean_density_pos;                // average density profile of 2D sampling
    vector<vector<double> > mean_sq_density_pos;            // average of square of density of 2D sampling
    vector<vector<double> > mean_density_neg;
    vector<vector<double> > mean_sq_density_neg;
    CONTROL cpmdremote;
    string np_shape;

public:

    NanoParticleDisk(string , vector<vector<BinRing> > &, vector<vector<BinRing> > &, double , double , vector<PARTICLE> &,
                     vector<vector<double> > &, vector<vector<double> > &, int , double , CONTROL &, VECTOR3D &, double ,
                     double , double , double );


// make bins for disk
    void make_bins();

    // bin ions to get density profile for disk
    void bin_ions();

    // compute initial density profile
    void compute_initial_density_profile();

// compute density profile disk
    void compute_density_profile();

    // compute final density profile
    void compute_final_density_profile();

    void updateStep(int );

    void updateSamples(double );

    string getType();

    void printType();

    //print number of bins used
    void printBinSize();

};

