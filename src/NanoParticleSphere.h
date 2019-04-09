// This is particle class


#include "NanoParticle.h"


class NanoParticleSphere : public NanoParticle {

private:
    vector<BinShell> bin_pos;
    vector<BinShell> bin_neg;
    double bin_width;
    vector<PARTICLE> *ion;
    vector<double> density_pos;
    vector<double> density_neg;
    int cpmdstep;
    double density_profile_samples;
    vector<double> mean_pos_density;                // average density profile
    vector<double> mean_pos_sq_density;            // average of square of density
    vector<double> mean_neg_density;                // average density profile
    vector<double> mean_neg_sq_density;            // average of square of density
    CONTROL cpmdremote;
    string np_shape;

public:

    NanoParticleSphere(string , vector<BinShell> &,vector<BinShell> &, double, vector<PARTICLE> &, vector<double> &, vector<double> &, int, double,
                       CONTROL &,VECTOR3D &, double, double, double, double);

    // members
    // make bins for disk
    void make_bins() ;

// bin ions to get density profile
    void bin_ions() ;

    // compute initial density profile
    void compute_initial_density_profile() ;

// compute density profile
    void compute_density_profile() ;

    // compute final density profile
    void compute_final_density_profile() ;

    void updateStep(int ) ;

    void updateSamples(double density_profile_samplesL) ;

    string getType() ;

    void printType();

    //print number of bins used
    void printBinSize() ;

};

