// This is particle class

#include "NanoParticleDisk.h"


NanoParticleDisk::NanoParticleDisk(string shape, vector<vector<BinRing> > &bin_posL, vector<vector<BinRing> > &bin_negL, double bin_width_RL,
                                   double bin_width_ZL, vector<PARTICLE> &ionL, vector<vector<double> > &density_posL, vector<vector<double> > &density_negL,
                                   int cpmdstepL, double density_profile_samplesL, CONTROL &cpmdremoteL,
                                   VECTOR3D &get_posvec,
                                   double get_radius, double get_ein, double get_eout,
                                   double get_bare_charge) {

    np_shape = shape;
    shape_id = 1;
    bin_posL = bin_posL;
    bin_negL = bin_negL;
    bin_width_R = bin_width_RL;
    bin_width_Z = bin_width_ZL;
    ion = &ionL;
    density_pos = density_posL;
    density_neg = density_negL;
    cpmdstep = cpmdstepL;
    density_profile_samples = density_profile_samplesL;
    cpmdremote = cpmdremoteL;
    posvec = get_posvec;
    radius = get_radius;
    ein = get_ein;
    eout = get_eout;
    bare_charge = get_bare_charge;
    area_np = set_surface_area_np(get_radius);

}


// make bins for disk
void NanoParticleDisk::make_bins() {

    unsigned int number_of_bins_R = int(box_radius / bin_width_R);
    unsigned int number_of_bins_Z = int(box_radius / bin_width_Z);
    bin_pos.resize(number_of_bins_Z, vector<BinRing>(number_of_bins_R));
    bin_neg.resize(number_of_bins_Z, vector<BinRing>(number_of_bins_R));

    for (unsigned int bin_num_Z = 0; bin_num_Z < bin_pos.size(); bin_num_Z++)
        for (unsigned int bin_num_R = 0; bin_num_R < bin_pos[bin_num_Z].size(); bin_num_R++) {
            bin_pos[bin_num_Z][bin_num_R].set_up(bin_num_Z, bin_num_R, bin_width_Z, bin_width_R);
            bin_neg[bin_num_Z][bin_num_R].set_up(bin_num_Z, bin_num_R, bin_width_Z, bin_width_R);
        }


    // This is only done for positive ions.
    if (world.rank() == 0) {
        ofstream listbin("outfiles/listbin.dat");
        for (unsigned int bin_num_Z = 0; bin_num_Z < bin_pos.size(); bin_num_Z++)
            for (unsigned int bin_num_R = 0; bin_num_R < bin_pos[bin_num_Z].size(); bin_num_R++)
                listbin << bin_pos[bin_num_Z][bin_num_R].n << setw(15) << bin_pos[bin_num_Z][bin_num_R].width_R << setw(15)
                        << setw(15) << bin_pos[bin_num_Z][bin_num_R].width_Z << setw(15)
                        << bin_pos[bin_num_Z][bin_num_R].volume
                        << setw(15) << bin_pos[bin_num_Z][bin_num_R].lower_R << setw(15)
                        << bin_pos[bin_num_Z][bin_num_R].higher_R << setw(15)
                        << bin_pos[bin_num_Z][bin_num_R].lower_Z << setw(15) << bin_pos[bin_num_Z][bin_num_R].higher_Z
                        << endl;
        listbin.close();
    }

    mean_density_pos.resize(bin_pos.size(), vector<double>(bin_pos[0].size()));
    mean_sq_density_pos.resize(bin_pos.size(), vector<double>(bin_pos[0].size()));
    mean_density_neg.resize(bin_neg.size(), vector<double>(bin_neg[0].size()));
    mean_sq_density_neg.resize(bin_neg.size(), vector<double>(bin_neg[0].size()));

    return;
}

// bin ions to get density profile for disk
void NanoParticleDisk::bin_ions() {

    for (unsigned int bin_num_Z = 0; bin_num_Z < bin_pos.size(); bin_num_Z++)
        for (unsigned int bin_num_R = 0; bin_num_R < bin_pos[bin_num_Z].size(); bin_num_R++) {
            bin_pos[bin_num_Z][bin_num_R].n = 0;
            bin_neg[bin_num_Z][bin_num_R].n = 0;
        }
    for (unsigned int i = 0; i < (*ion).size(); i++) {

        int binNumberZ = int(fabs((*ion)[i].posvec.z) / bin_pos[0][0].width_Z);
        int binNumberR = int(sqrt(pow((*ion)[i].posvec.x, 2) + pow((*ion)[i].posvec.y, 2)) / bin_pos[0][0].width_R);

        if ((*ion)[i].valency > 0) {
            bin_pos[binNumberZ][binNumberR].n = bin_pos[binNumberZ][binNumberR].n + 1;
        }else{
            //Assuming 0 valency never exists
            bin_neg[binNumberZ][binNumberR].n = bin_neg[binNumberZ][binNumberR].n + 1;
        }
    }

    density_pos.resize(bin_pos.size(), vector<double>(bin_pos[0].size()));
    density_neg.resize(bin_neg.size(), vector<double>(bin_neg[0].size()));

    for (unsigned int bin_num_Z = 0; bin_num_Z < bin_pos.size(); bin_num_Z++)
        for (unsigned int bin_num_R = 0; bin_num_R < bin_pos[bin_num_Z].size(); bin_num_R++) {
            density_pos[bin_num_Z][bin_num_R] = bin_pos[bin_num_Z][bin_num_R].n /
                                            bin_pos[bin_num_Z][bin_num_R].volume;            // push_back is the culprit, array goes out of bound
            density_neg[bin_num_Z][bin_num_R] = bin_neg[bin_num_Z][bin_num_R].n /
                                            bin_neg[bin_num_Z][bin_num_R].volume;
        }
    return;
}

// compute initial density profile
void NanoParticleDisk::compute_initial_density_profile() {

    density_pos.clear();
    density_neg.clear();

    if (world.rank() == 0) {
        bin_ions();
        ofstream density_pos_profile("outfiles/initial_positive_density_profile.dat", ios::out);
        ofstream density_neg_profile("outfiles/initial_negative_density_profile.dat", ios::out);

        // initial density for disk
        for (unsigned int b = 0; b < density_pos.size(); b++)
            for (unsigned int c = 0; c < density_pos[b].size(); c++) {
                density_pos_profile << b * bin_pos[b][c].width_Z << setw(15) << c * bin_pos[b][c].width_R << setw(15)
                                << density_pos[b][c] << endl;
                density_neg_profile << b * bin_neg[b][c].width_Z << setw(15) << c * bin_neg[b][c].width_R << setw(15)
                                << density_neg[b][c] << endl;
            }
        density_pos_profile.close();
        density_neg_profile.close();
    }

}

// compute density profile disk
void NanoParticleDisk::compute_density_profile() {

    density_pos.clear();
    density_neg.clear();

    if (world.rank() == 0) {
        ofstream file_for_auto_corr("outfiles/for_auto_corr.dat", ios::app);

        bin_ions();

        for (unsigned int b = 0; b < mean_density_pos.size(); b++)
            for (unsigned int c = 0; c < mean_density_pos[b].size(); c++) {
                mean_density_pos[b][c] = mean_density_pos[b][c] + density_pos[b][c];
                mean_density_neg[b][c] = mean_density_neg[b][c] + density_neg[b][c];
            }

        for (unsigned int b = 0; b < density_pos.size(); b++)
            for (unsigned int c = 0; c < density_pos[b].size(); c++) {
                mean_sq_density_pos[b][c] = mean_sq_density_pos[b][c] + (density_pos[b][c] * density_pos[b][c]);
                mean_sq_density_neg[b][c] = mean_sq_density_neg[b][c] + (density_neg[b][c] * density_neg[b][c]);
            }

        // write a file for post analysis to get auto correlation time		// NOTE this is assuming ions do not cross the interface
        //Only used positive ion densities here
        if ((*ion)[0].posvec.GetMagnitude() > radius)
            file_for_auto_corr << cpmdstep - cpmdremote.hiteqm << "\t"
                               << density_pos[int((radius + 2) / bin_pos[0][0].width_Z)][int(
                                       (radius + 2) / bin_pos[0][0].width_R)] << endl;
        else
            file_for_auto_corr << cpmdstep - cpmdremote.hiteqm << "\t"
                               << density_pos[int((radius - 2) / bin_pos[0][0].width_Z)][int(
                                       (radius - 2) / bin_pos[0][0].width_R)] << endl;

        // write files
        if (cpmdstep % cpmdremote.writedensity == 0) {
            char data_pos[200], data_neg[200];
            sprintf(data_pos, "datafiles/_den_pos_%.06d.dat", cpmdstep);
            sprintf(data_neg, "datafiles/_den_neg_%.06d.dat", cpmdstep);
            ofstream outden_pos, outden_neg;
            outden_pos.open(data_pos);
            outden_neg.open(data_neg);
            for (unsigned int b = 0; b < mean_density_pos.size(); b++)
                for (unsigned int c = 0; c < mean_density_pos[b].size(); c++) {
                    outden_pos << b * bin_pos[b][c].width_Z << setw(15) << c * bin_pos[b][c].width_R << setw(15)
                           << mean_density_pos[b][c] / density_profile_samples << endl;
                    outden_neg << b * bin_neg[b][c].width_Z << setw(15) << c * bin_neg[b][c].width_R << setw(15)
                           << mean_density_neg[b][c] / density_profile_samples << endl;

                }
            outden_pos.close();
            outden_neg.close();
        }
    }
    return;
}

// compute final density profile
void NanoParticleDisk::compute_final_density_profile() {
    if (world.rank() == 0) {
        // 1. density profile
        vector<vector<double> > density_profile_pos;
        vector<vector<double> > density_profile_neg;
        density_profile_pos.resize(bin_pos.size(), vector<double>(bin_pos[0].size()));
        density_profile_neg.resize(bin_neg.size(), vector<double>(bin_neg[0].size()));
        for (unsigned int b = 0; b < mean_density_pos.size(); b++)
            for (unsigned int c = 0; c < mean_density_pos[b].size(); c++) {
                density_profile_pos[b][c] = (mean_density_pos[b][c] / density_profile_samples);
                density_profile_neg[b][c] = (mean_density_neg[b][c] / density_profile_samples);
            }

        // 2. error bars
        vector<vector<double> > error_bar_disk_pos;
        vector<vector<double> > error_bar_disk_neg;
        error_bar_disk_pos.resize(bin_pos.size(), vector<double>(bin_pos[0].size()));
        error_bar_disk_neg.resize(bin_neg.size(), vector<double>(bin_neg[0].size()));
        for (unsigned int b = 0; b < density_profile_pos.size(); b++)
            for (unsigned int c = 0; c < density_profile_pos[b].size(); c++) {
                error_bar_disk_pos[b][c] =
                        sqrt(1.0 / density_profile_samples) * sqrt(mean_sq_density_pos[b][c] / density_profile_samples -
                                                                   density_profile_pos[b][c] * density_profile_pos[b][c]);
                error_bar_disk_neg[b][c] =
                        sqrt(1.0 / density_profile_samples) * sqrt(mean_sq_density_neg[b][c] / density_profile_samples -
                                                                   density_profile_neg[b][c] * density_profile_neg[b][c]);
            }

        // 3. write results
        ofstream list_profile_pos("outfiles/positive_density_profile.dat", ios::out);
        ofstream list_profile_neg("outfiles/negative_density_profile.dat", ios::out);
        for (unsigned int b = 0; b < density_profile_pos.size(); b++)
            for (unsigned int c = 0; c < density_profile_pos[b].size(); c++) {
                list_profile_pos << b * bin_pos[b][c].width_Z << setw(15) << c * bin_pos[b][c].width_R << setw(15)
                             << density_profile_pos[b][c] << setw(15) << error_bar_disk_pos[b][c] << endl;
                list_profile_neg << b * bin_neg[b][c].width_Z << setw(15) << c * bin_neg[b][c].width_R << setw(15)
                             << density_profile_neg[b][c] << setw(15) << error_bar_disk_neg[b][c] << endl;
            }

        list_profile_pos.close();
        list_profile_neg.close();

    }
}

void NanoParticleDisk::updateStep(int cpmdstepL) {

    cpmdstep = cpmdstepL;

}

void NanoParticleDisk::updateSamples(double density_profile_samplesL) {

    density_profile_samples = density_profile_samplesL;

}

string NanoParticleDisk::getType() {

    return np_shape;

}

void NanoParticleDisk::printType() {

    cout << np_shape << endl;
}

//print number of bins used
void NanoParticleDisk::printBinSize() {

    if (world.rank() == 0) {
        cout << "Number of bins used (Z direction) for computing density profiles " << bin_pos.size() << endl;
        cout << "Number of bins used (R direction) for computing density profiles " << bin_pos[0].size() << endl;
    }
}



