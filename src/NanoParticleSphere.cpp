#include "NanoParticleSphere.h"

NanoParticleSphere::NanoParticleSphere(string shape, vector<BinShell> &bin_pos_L, vector<BinShell> &bin_neg_L,
                                       double bin_widthL,
                                       vector<PARTICLE> &ionL, vector<double> &density_posL, vector<double> &density_negL, int cpmdstepL,
                                       double density_profile_samplesL, CONTROL &cpmdremoteL, VECTOR3D &get_posvec,
                                       double get_radius = 0, double get_ein = 1, double get_eout = 1,
                                       double get_bare_charge = 0) {

    np_shape = shape;
    shape_id = 0;
    bin_pos = bin_pos_L;
    bin_pos = bin_neg_L;
    bin_width = bin_widthL;
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
void NanoParticleSphere::make_bins() {


    unsigned int number_of_bins = int(box_radius / bin_width);
    bin_pos.resize(number_of_bins);
    bin_neg.resize(number_of_bins);
    for (unsigned int bin_num = 0; bin_num < bin_pos.size(); bin_num++) {
        bin_pos[bin_num].set_up(bin_num, bin_width);
        bin_neg[bin_num].set_up(bin_num, bin_width);
    }

    // This is only done for positive ions.
    if (world.rank() == 0) {
        ofstream listbin("outfiles/listbin.dat");
        for (unsigned int num = 0; num < bin_pos.size(); num++)
            listbin << bin_pos[num].n << setw(15) << bin_pos[num].width << setw(15) << bin_pos[num].volume << setw(15)
                    << bin_pos[num].lower
                    << setw(15) << bin_pos[num].higher << endl;
        listbin.close();
    }

    //set mean and mean sq vectors
    for (unsigned int b = 0; b < bin_pos.size(); b++) {
        mean_pos_density.push_back(0.0);
        mean_pos_sq_density.push_back(0.0);
        mean_neg_density.push_back(0.0);
        mean_neg_sq_density.push_back(0.0);
    }

    return;
}

// bin ions to get density profile
void NanoParticleSphere::bin_ions() {
    double r;
    int bin_number;
    for (unsigned int bin_num = 0; bin_num < bin_pos.size(); bin_num++) {
        bin_pos[bin_num].n = 0;
        bin_neg[bin_num].n = 0;
    }
    for (unsigned int i = 0; i < (*ion).size(); i++) {
        r = (*ion)[i].posvec.GetMagnitude();
        bin_number = int(r / bin_pos[0].width);
        if ((*ion)[i].valency > 0) {
            bin_pos[bin_number].n = bin_pos[bin_number].n + 1;
        } else {
            //Assuming 0 valency never exists
            bin_neg[bin_number].n = bin_neg[bin_number].n + 1;
        }
    }
    for (unsigned int bin_num = 0; bin_num < bin_pos.size(); bin_num++) {
        density_pos.push_back(
                bin_pos[bin_num].n /
                bin_pos[bin_num].volume);            // push_back is the culprit, array goes out of bound

        density_neg.push_back(
                bin_neg[bin_num].n /
                bin_neg[bin_num].volume);            // push_back is the culprit, array goes out of bound
    }
    return;
}


// compute initial density profile
void NanoParticleSphere::compute_initial_density_profile() {

    density_pos.clear();
    density_neg.clear();
    if (world.rank() == 0) {
        bin_ions();
        ofstream density_profile_pos("outfiles/initial_positive_density_profile.dat", ios::out);
        ofstream density_profile_neg("outfiles/initial_negative_density_profile.dat", ios::out);

        for (unsigned int b = 0; b < density_pos.size(); b++) {
            density_profile_pos << b * bin_pos[b].width << setw(15) << density_pos.at(b) << endl;
            density_profile_neg << b * bin_neg[b].width << setw(15) << density_neg.at(b) << endl;
        }
        density_profile_pos.close();
        density_profile_neg.close();
    }

}

// compute density profile
void NanoParticleSphere::compute_density_profile() {

    density_pos.clear();
    density_neg.clear();

    if (world.rank() == 0) {
        ofstream file_for_auto_corr("outfiles/for_auto_corr.dat", ios::app);

        bin_ions();

        for (unsigned int b = 0; b < mean_pos_density.size(); b++)
            mean_pos_density.at(b) = mean_pos_density.at(b) + density_pos.at(b);
        for (unsigned int b = 0; b < density_pos.size(); b++)
            mean_pos_sq_density.at(b) = mean_pos_sq_density.at(b) + density_pos.at(b) * density_pos.at(b);
        for (unsigned int b = 0; b < mean_neg_density.size(); b++)
            mean_neg_density.at(b) = mean_neg_density.at(b) + density_neg.at(b);
        for (unsigned int b = 0; b < density_neg.size(); b++)
            mean_neg_sq_density.at(b) = mean_neg_sq_density.at(b) + density_neg.at(b) * density_neg.at(b);

        // write a file for post analysis to get auto correlation time		// NOTE this is assuming ions do not cross the interface
        //Only used positive ion densities here
        if ((*ion)[0].posvec.GetMagnitude() > radius)
            file_for_auto_corr << cpmdstep - cpmdremote.hiteqm << "\t"
                               << density_pos[int((radius + 2) / bin_pos[0].width)] << endl;
        else
            file_for_auto_corr << cpmdstep - cpmdremote.hiteqm << "\t"
                               << density_pos[int((radius - 2) / bin_pos[0].width)] << endl;

        // write files
        if (cpmdstep % cpmdremote.writedensity == 0) {
            char data_pos[200], data_neg[200];
            sprintf(data_pos, "datafiles/_den_pos_%.06d.dat", cpmdstep);
            sprintf(data_neg, "datafiles/_den_neg_%.06d.dat", cpmdstep);
            ofstream outden_pos, outden_neg;
            outden_pos.open(data_pos);
            outden_neg.open(data_neg);
            for (unsigned int b = 0; b < mean_pos_density.size(); b++) {
                outden_pos << b * bin_pos[b].width << setw(15) << mean_pos_density.at(b) / density_profile_samples
                           << endl;
                outden_neg << b * bin_neg[b].width << setw(15) << mean_neg_density.at(b) / density_profile_samples
                           << endl;
            }
            outden_pos.close();
            outden_neg.close();
        }
    }
    return;
}


// compute final density profile
void NanoParticleSphere::compute_final_density_profile() {

    // 1. density profile
    if (world.rank() == 0) {

        vector<double> density_profile_pos;
        vector<double> density_profile_neg;

        for (unsigned int b = 0; b < mean_pos_density.size(); b++) {
            density_profile_pos.push_back(mean_pos_density.at(b) / density_profile_samples);
            density_profile_neg.push_back(mean_neg_density.at(b) / density_profile_samples);
        }

        // 2. error bars
        vector<double> error_bar_pos;
        vector<double> error_bar_neg;
        for (unsigned int b = 0; b < density_profile_pos.size(); b++) {
            error_bar_pos.push_back(
                    sqrt(1.0 / density_profile_samples) * sqrt(mean_pos_sq_density.at(b) / density_profile_samples -
                                                               density_profile_pos.at(b) * density_profile_pos.at(b)));
            error_bar_neg.push_back(
                    sqrt(1.0 / density_profile_samples) * sqrt(mean_neg_sq_density.at(b) / density_profile_samples -
                                                               density_profile_neg.at(b) * density_profile_neg.at(b)));
        }

        // 3. write results
        ofstream list_profile_pos("outfiles/positive_density_profile.dat", ios::out);
        ofstream list_profile_neg("outfiles/negative_density_profile.dat", ios::out);

        for (unsigned int b = 0; b < density_profile_pos.size(); b++) {
            list_profile_pos << b * bin_pos[b].width << setw(15) << density_profile_pos.at(b) << setw(15)
                             << error_bar_pos.at(b)
                             << endl;
            list_profile_neg << b * bin_neg[b].width << setw(15) << density_profile_neg.at(b) << setw(15)
                             << error_bar_neg.at(b)
                             << endl;
        }

        list_profile_pos.close();
        list_profile_neg.close();
    }
}

void NanoParticleSphere::updateStep(int cpmdstepL) {

    cpmdstep = cpmdstepL;

}

void NanoParticleSphere::updateSamples(double density_profile_samplesL) {

    density_profile_samples = density_profile_samplesL;

}

string NanoParticleSphere::getType() {

    return np_shape;

}

void NanoParticleSphere::printType() {

    cout << np_shape << endl;
}

//print number of bins used
void NanoParticleSphere::printBinSize() {

    if (world.rank() == 0)
        cout << "Number of bins used for computing density profiles " << bin_pos.size() << endl;
}



