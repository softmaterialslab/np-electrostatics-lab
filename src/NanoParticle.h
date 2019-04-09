// This is particle class

#ifndef _NANOPARTICLE_H
#define _NANOPARTICLE_H

#include <iostream>
#include<string>
#include "particle.h"
#include "vector3d.h"
#include "vertex.h"
#include "control.h"
#include "BinShell.h"
#include "BinRing.h"
#include "thermostat.h"
#include "mpi_utility.h"



using namespace std;

class NanoParticle {

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & posvec;
        ar & radius;
        ar & area_np;
        ar & ein;
        ar & eout;
        ar & em;
        ar & ed;
        ar & box_radius;
        ar & lB_in;
        ar & lB_out;
        ar & inv_kappa_in;
        ar & inv_kappa_out;
        ar & mean_sep_in;
        ar & mean_sep_out;
        ar & number_of_vertices;
        ar & bare_charge;
        ar & POLARIZED;
        ar & RANDOMIZE_ION_FEATURES;

    }



public:

    VECTOR3D posvec;        // position vector of the inteface, for sphere its center
    double radius;        // radius of the dielectric sphere, interface geometry variable
    double ein;            // permittivity of inside medium
    double eout;            // permittivity of outside medium
    double em;            // permittivity at the interface, mean
    double ed;            // permittivity change at the inteface, difference scaled by 4 pi
    double box_radius;        // radius of the spherical box in which all the action occurs
    double lB_in;            // Bjerrum length inside
    double lB_out;        // Bjerrum length outside
    double inv_kappa_in;        // debye length inside
    double inv_kappa_out;        // debye length outside
    double mean_sep_in;        // mean separation inside
    double mean_sep_out;        // mean separation outside
    int number_of_vertices;    // number of points used to discretize the interface
    double bare_charge;        // bare charge on the membrane (interface)
    double area_np;        // area of the nano particle
    bool POLARIZED;        // is the nanoparticle polarized; depends on ein, eout
    bool RANDOMIZE_ION_FEATURES;    // are selections randomized
    int shape_id = -1;                   //Shape id number -> initialized to non type

    // make a particle constructor
    NanoParticle();

    // member functions definitions

    void set_up(double, double, double, double, int, double);

    void put_counterions(vector<PARTICLE> &, int, double, vector<PARTICLE> &);

    void put_saltions_inside(vector<PARTICLE> &, int, double, double, vector<PARTICLE> &);

    void put_saltions_outside(vector<PARTICLE> &, int, double, double, vector<PARTICLE> &);

    void discretize(vector<VERTEX> &);

    // total charge inside
    double total_charge_inside(vector<PARTICLE> &);

    // total induced charge
    double total_induced_charge(vector<VERTEX> &);

    // calculate the area of np
    double set_surface_area_np(double);

    // member functions definitions for child classes
    // make bins for disk
    virtual void make_bins();

    // bin ions to get density profile for disk
    virtual void bin_ions();

    // compute initial density profile
    virtual void compute_initial_density_profile();

    // compute density profile
    virtual void compute_density_profile();

    // compute final density profile
    virtual void compute_final_density_profile();

    //update time step
    virtual void updateStep(int );

    //update the number of samples used for density profile
    virtual void updateSamples(double );

    //get NP type
    virtual string getType();

    //print number of bins used
    virtual void printBinSize();

    //get NP type
    virtual void printType();

};

#endif