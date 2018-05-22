// This file contains member functions for interface class

#include "interface.h"
#include "functions.h"

void
INTERFACE::set_up(double salt_conc_in, double salt_conc_out, double salt_valency_in, double salt_valency_out, int N,
                  double b) {
    
    
    // useful combinations of different dielectric constants (inside and outside)
    em = 0.5 * (ein + eout);
    ed = (eout - ein) / (4 * pi);

    // useful length scales signifying competition between electrostatics and entropy
    lB_in = (lB_water * epsilon_water / ein) / unitlength;
    lB_out = (lB_water * epsilon_water / eout) / unitlength;
    if (salt_conc_in != 0) {
        inv_kappa_in = (0.257 / (salt_valency_in * sqrt(lB_in * unitlength * salt_conc_in))) / unitlength;
        mean_sep_in = pow(1.2 * salt_conc_in, -1.0 / 3.0) / unitlength;
    } else {
        inv_kappa_in = 0;
        mean_sep_in = 0;
    }
    if (salt_conc_out != 0) {
        inv_kappa_out = (0.257 / (salt_valency_out * sqrt(lB_out * unitlength * salt_conc_out))) / unitlength;
        mean_sep_out = pow(1.2 * salt_conc_out, -1.0 / 3.0) / unitlength;
    } else {
        inv_kappa_out = 0;
        mean_sep_out = 0;
    }

    // simulation box size
    if (salt_conc_out == 0 && salt_conc_in == 0)
        box_radius = b;
    else
        box_radius = b;
//    box_radius = radius + 0.5 + 5 * inv_kappa_out;				// NOTE change this number for testing purposes...
    
    // discretization parameters
    number_of_vertices = N;

    return;
}

void
INTERFACE::put_counterions(vector<PARTICLE> &counterion, int ion_valency, double ion_diameter, vector<PARTICLE> &ion) {
    
    
    // establish the number of counterions first
    unsigned int total_counterions = int(abs(bare_charge / ion_valency));
    // express diameter in consistent units
    ion_diameter = ion_diameter / unitlength;
    // distance of closest approach between counterion and the nanoparticle
    double r0 = radius + 0.5 * ion_diameter;
    // distance of closest approach between counterion and the simulation box boundary (treated as a sphere)
    double r0_box = box_radius - 0.5 * ion_diameter;

    UTILITY ugsl;

    // generate counterions in the box
    while (counterion.size() != total_counterions) {
        double x = gsl_rng_uniform(ugsl.r);
        x = (1 - x) * (-r0_box) + x * (r0_box);
        double y = gsl_rng_uniform(ugsl.r);
        y = (1 - y) * (-r0_box) + y * (r0_box);
        double z = gsl_rng_uniform(ugsl.r);
        z = (1 - z) * (-r0_box) + z * (r0_box);
        VECTOR3D posvec = VECTOR3D(x, y, z);
        if (posvec.GetMagnitude() < r0 +
                                    ion_diameter)                        // giving an extra ion diameter length to be safe with the initial generation
            continue;
        if (posvec.GetMagnitude() >= box_radius - ion_diameter)
            continue;
        bool continuewhile = false;
        for (unsigned int i = 0; i < ion.size() && !continuewhile; i++)
            if ((posvec - ion[i].posvec).GetMagnitude() <= ion_diameter)
                continuewhile = true;        // avoid overlap of initial ion positions
        if (continuewhile)
            continue;
        counterion.push_back(PARTICLE(int(ion.size()) + 1, ion_diameter, ion_valency, ion_valency * 1.0, 1.0, eout,
                                      posvec));        // create a counterion
        ion.push_back(PARTICLE(int(ion.size()) + 1, ion_diameter, ion_valency, ion_valency * 1.0, 1.0, eout,
                               posvec));            // copy the counterion as a general ion
    }
    if (world.rank() == 0) {
        ofstream listcounterions("outfiles/counterions.xyz", ios::out);
        listcounterions << counterion.size() << endl;
        listcounterions << "counterions" << endl;

        for (unsigned int i = 0; i < counterion.size(); i++)
            listcounterions << "C" << setw(15) << counterion[i].posvec << setw(15)
                            << counterion[i].posvec.GetMagnitude() << endl;
        listcounterions.close();
    }
    return;
}

void INTERFACE::put_saltions_inside(vector<PARTICLE> &saltion_in, int valency, double concentration, double diameter,
                                    vector<PARTICLE> &ion) {
    
    
    // establish the number of inside salt ions first
    // Note: salt concentration is the concentration of one kind of ions, so for total ions a factor of 2 needs to be multiplied. also some factors appear to be consistent with units.
    double volume_sphere = (4.0 / 3.0) * pi * radius * radius * radius;
    unsigned int total_saltions_inside = int(2 * (concentration * 0.6) * (volume_sphere * unitlength * unitlength *
                                                                          unitlength));        // NOTE there is a change in the factor, i think it was wrong before
    if (total_saltions_inside % 2 != 0)
        total_saltions_inside = total_saltions_inside + 1;
    // express diameter in consistent units
    diameter = diameter / unitlength;
    // distance of closest approach between the ion and the interface
    double r0 = radius - 0.5 * diameter;

    UTILITY ugsl;                                                    // utility used for making initial configuration

    // generate salt ions inside
    while (saltion_in.size() != total_saltions_inside) {
        double x = gsl_rng_uniform(ugsl.r);
        x = (1 - x) * (-r0) + x * (r0);
        double y = gsl_rng_uniform(ugsl.r);
        y = (1 - y) * (-r0) + y * (r0);
        double z = gsl_rng_uniform(ugsl.r);
        z = (1 - z) * (-r0) + z * (r0);
        VECTOR3D posvec = VECTOR3D(x, y, z);
        if (posvec.GetMagnitude() > r0 -
                                    diameter)                                    // putting an extra ion diameter length away from interface
            continue;
        bool continuewhile = false;
        for (unsigned int i = 0; i < ion.size() && continuewhile == false; i++)
            if ((posvec - ion[i].posvec).GetMagnitude() <= (0.5 * diameter + 0.5 * ion[i].diameter))
                continuewhile = true;
        if (continuewhile == true)
            continue;
        saltion_in.push_back(PARTICLE(int(ion.size()) + 1, diameter, valency, valency * 1.0, 1.0, ein,
                                      posvec));        // create a salt ion
        ion.push_back(PARTICLE(int(ion.size()) + 1, diameter, valency, valency * 1.0, 1.0, ein,
                               posvec));            // copy the salt ion to the stack of all ions
        valency = (-1) *
                  valency;                                            // switch between creating positive and negative ion
    }
    if (world.rank() == 0) {
        ofstream list_salt_ions_inside("outfiles/salt_ions_inside.xyz", ios::out);
        list_salt_ions_inside << saltion_in.size() << endl;
        list_salt_ions_inside << "salt ions inside" << endl;
        for (unsigned int i = 0; i < saltion_in.size(); i++)
            list_salt_ions_inside << "Si" << setw(15) << saltion_in[i].posvec << endl;
        list_salt_ions_inside.close();
    }
    return;
}

void INTERFACE::put_saltions_outside(vector<PARTICLE> &saltion_out, int valency, double concentration, double diameter,
                                     vector<PARTICLE> &ion) {
    
    
    // establish the number of outside salt ions first
    // Note: salt concentration is the concentration of one kind of ions, so for total ions a factor of 2 needs to be multiplied. also some factors appear to be consistent with units.
    double volume_box = (4.0 / 3.0) * pi * (box_radius * box_radius * box_radius - radius * radius * radius);
    unsigned int total_saltions_outside = int(2 * (concentration * 0.6) * (volume_box * unitlength * unitlength *
                                                                           unitlength));        // NOTE there is a change in the factor, i think it was wrong before
    if (total_saltions_outside % 2 != 0)
        total_saltions_outside = total_saltions_outside + 1;
    // express diameter in consistent units
    diameter = diameter / unitlength;
    // distance of closest approach between the ion and the interface
    double r0 = radius + 0.5 * diameter;
    // distance of closest approach between the ion and the box
    double r0_box = box_radius - 0.5 * diameter;

    UTILITY ugsl;                                                        // utility used for making initial configuration

    // generate salt ions outside
    while (saltion_out.size() != total_saltions_outside) {
        double x = gsl_rng_uniform(ugsl.r);
        x = (1 - x) * (-r0_box) + x * (r0_box);
        double y = gsl_rng_uniform(ugsl.r);
        y = (1 - y) * (-r0_box) + y * (r0_box);
        double z = gsl_rng_uniform(ugsl.r);
        z = (1 - z) * (-r0_box) + z * (r0_box);
        VECTOR3D posvec = VECTOR3D(x, y, z);
        if (posvec.GetMagnitude() <
            r0 + diameter)                                        // giving an extra ion diameter length
            continue;
        if (posvec.GetMagnitude() >= r0_box)
            continue;
        bool continuewhile = false;
        for (unsigned int i = 0; i < ion.size() && continuewhile == false; i++)
            if ((posvec - ion[i].posvec).GetMagnitude() <= (0.5 * diameter + 0.5 * ion[i].diameter))
                continuewhile = true;        // avoid overlap of initial ion positions
        if (continuewhile == true)
            continue;
        saltion_out.push_back(PARTICLE(int(ion.size()) + 1, diameter, valency, valency * 1.0, 1.0, eout,
                                       posvec));            // create a salt ion
        ion.push_back(PARTICLE(int(ion.size()) + 1, diameter, valency, valency * 1.0, 1.0, eout,
                               posvec));                // copy the salt ion to the stack of all ions
        valency = (-1) *
                  valency;                                                // switch between creating positive and negative ion
    }
    if (world.rank() == 0) {
        ofstream list_salt_ions_outside("outfiles/salt_ions_outside.xyz", ios::out);
        list_salt_ions_outside << saltion_out.size() << endl;
        list_salt_ions_outside << "salt ions outside" << endl;
        for (unsigned int i = 0; i < saltion_out.size(); i++)
            list_salt_ions_outside << "So" << setw(15) << saltion_out[i].posvec << endl;
        list_salt_ions_outside.close();
    }
    return;
}

// discretize interface
void INTERFACE::discretize(vector<VERTEX> &s) {
    
    
    //for the disk temp change. : infiles_a10_disk
    char filename[200];
    sprintf(filename, "infiles_a7.5/grid%d.dat",
            number_of_vertices); // change infiles folder if nanoparticle radius changes; for a = 2.67m nm = 7.5 sigma in reduced units, infiles_a7.5 is the folder
    ifstream in(filename, ios::in);
    if (!in) {
        if (world.rank() == 0)
            cout << "File could not be opened" << endl;
        exit(1);
    }

    unsigned int col1;
    double col2, col3, col4, col5, col6, col7, col8;
    while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8)
        s.push_back(VERTEX(VECTOR3D(col2, col3, col4), col5, VECTOR3D(col6, col7, col8), area_np, bare_charge));

    if (world.rank() == 0) {
        ofstream listvertices("outfiles/interface.xyz", ios::out);
        listvertices << number_of_vertices << endl;
        listvertices << "interface" << endl;
        for (unsigned int k = 0; k < s.size(); k++)
            listvertices << "I" << setw(15) << s[k].posvec << endl;
        listvertices.close();
    }
    return;
}

