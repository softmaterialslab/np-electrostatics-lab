// This is a header file for the INTERFACE class.  

#ifndef _INTERFACE_H
#define _INTERFACE_H

#include "utility.h"
#include "vector3d.h"
#include "particle.h"
#include "vertex.h"

class INTERFACE {

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

    double area_np;		// area of the nano particle

    bool POLARIZED;        // is the nanoparticle polarized; depends on ein, eout

    bool RANDOMIZE_ION_FEATURES;    // are selections randomized

    void set_up(double, double, double, double, int, double);

    void put_counterions(vector<PARTICLE> &, int, double, vector<PARTICLE> &);

    void put_saltions_inside(vector<PARTICLE> &, int, double, double, vector<PARTICLE> &);

    void put_saltions_outside(vector<PARTICLE> &, int, double, double, vector<PARTICLE> &);

    void discretize(vector<VERTEX> &);

    INTERFACE(VECTOR3D get_posvec = VECTOR3D(0, 0, 0), double get_radius = 0, double get_ein = 1, double get_eout = 1,
              double get_bare_charge = 0) {
        posvec = get_posvec;
        radius = get_radius;
        ein = get_ein;
        eout = get_eout;
        bare_charge = get_bare_charge;
        area_np = set_surface_area_np(get_radius);
    }

    // total charge inside
    double total_charge_inside(vector<PARTICLE> &ion) {
        double charge = 0;
        for (unsigned int i = 0; i < ion.size(); i++)
            if (ion[i].posvec.GetMagnitude() < radius) charge += ion[i].q;
        return charge;
    }

    // total induced charge
    double total_induced_charge(vector<VERTEX> &s) {
        double charge = 0;
        for (unsigned int k = 0; k < s.size(); k++)
            charge += s[k].w * s[k].a;
        return charge;
    }
    // calculate the area of np
    double set_surface_area_np(double radius){

        return (4.0) * pi * pow(radius, 2);

    }
};

#endif

