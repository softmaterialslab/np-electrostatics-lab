// This file contains the routines 

#include "functions.h"

// overload out
ostream &operator<<(ostream &os, VECTOR3D vec) {
    os << vec.x << setw(15) << vec.y << setw(15) << vec.z;
    return os;
}

// make bins
void make_bins(vector<BIN> &bin, INTERFACE &nanoparticle, double bin_width) {


    unsigned int number_of_bins = int(nanoparticle.box_radius / bin_width);
    bin.resize(number_of_bins);
    for (unsigned int bin_num = 0; bin_num < bin.size(); bin_num++)
        bin[bin_num].set_up(bin_num, bin_width);
    if (world.rank() == 0) {
        ofstream listbin("outfiles/listbin.dat");
        for (unsigned int num = 0; num < bin.size(); num++)
            listbin << bin[num].n << setw(15) << bin[num].width << setw(15) << bin[num].volume << setw(15)
                    << bin[num].lower
                    << setw(15) << bin[num].higher << endl;
        listbin.close();
    }
    return;
}

// bin ions to get density profile
void bin_ions(vector<PARTICLE> &ion, INTERFACE &nanoparticle, vector<double> &density, vector<BIN> &bin) {
    double r;
    int bin_number;
    for (unsigned int bin_num = 0; bin_num < bin.size(); bin_num++)
        bin[bin_num].n = 0;
    for (unsigned int i = 0; i < ion.size(); i++) {
        r = ion[i].posvec.GetMagnitude();
        bin_number = int(r / bin[0].width);
        bin[bin_number].n = bin[bin_number].n + 1;
    }
    for (unsigned int bin_num = 0; bin_num < bin.size(); bin_num++)
        density.push_back(
                bin[bin_num].n / bin[bin_num].volume);            // push_back is the culprit, array goes out of bound
    return;
}

// initialize velocities of particles to start simulation
void initialize_particle_velocities(vector<PARTICLE> &ion, vector<THERMOSTAT> &bath, INTERFACE &nanoparticle) {


    if (bath.size() == 1) {
        for (unsigned int i = 0; i < ion.size(); i++)
            ion[i].velvec = VECTOR3D(0, 0, 0);                    // initialized velocities
        if (world.rank() == 0)
            cout << "No thermostat for real system" << endl;
        return;
    }

    if (nanoparticle.RANDOMIZE_ION_FEATURES) {
        double p_sigma = sqrt(kB * bath[0].T / (2.0 * ion[0].m));        // Maxwell distribution width
        UTILITY ugsl;
        for (unsigned int i = 0; i < ion.size(); i++)
            ion[i].velvec = VECTOR3D(gsl_ran_gaussian(ugsl.r, p_sigma), gsl_ran_gaussian(ugsl.r, p_sigma),
                                     gsl_ran_gaussian(ugsl.r, p_sigma));    // initialized velocities
    } else {
        for (unsigned int i = 0; i < ion.size(); i++)
            ion[i].velvec = VECTOR3D(0, 0, 0);
    }

    VECTOR3D average_velocity_vector = VECTOR3D(0, 0, 0);
    for (unsigned int i = 0; i < ion.size(); i++)
        average_velocity_vector = average_velocity_vector + ion[i].velvec;
    average_velocity_vector = average_velocity_vector ^ (1.0 / ion.size());
    for (unsigned int i = 0; i < ion.size(); i++)
        ion[i].velvec = ion[i].velvec - average_velocity_vector;
    return;
}

// initialize velocities of fake degrees to start simulation
void initialize_fake_velocities(vector<VERTEX> &s, vector<THERMOSTAT> &fake_bath, INTERFACE &nanoparticle) {


    if (fake_bath.size() == 1)    // only the dummy is there
    {
        for (unsigned int k = 0; k < s.size(); k++)
            s[k].vw = 0.0;
        if (world.rank() == 0)
            cout << "No thermostat for fake system" << endl;
        return;
    }

    if (nanoparticle.RANDOMIZE_ION_FEATURES) {
        double w_sigma = sqrt(kB * fake_bath[0].T / (2.0 * s[0].mu));        // Maxwell distribution width
        UTILITY ugsl;
        for (unsigned int k = 0; k < s.size(); k++) {
            s[k].vw = gsl_ran_gaussian(ugsl.r, w_sigma);    // initialized velocities
        }
    } else {
        for (unsigned int k = 0; k < s.size(); k++)
            s[k].vw = 0.0;
    }

    double average_fake_momentum = 0.0;
    for (unsigned int k = 0; k < s.size(); k++)
        average_fake_momentum = average_fake_momentum + s[k].mu * s[k].vw;
    average_fake_momentum = average_fake_momentum * (1.0 / s.size());
    for (unsigned int k = 0; k < s.size(); k++)
        s[k].vw = s[k].vw - average_fake_momentum / (s[k].mu);
    return;
}

// compute additional quantities
void compute_n_write_useful_data(int cpmdstep, vector<PARTICLE> &ion, vector<VERTEX> &s, vector<THERMOSTAT> &real_bath,
                                 vector<THERMOSTAT> &fake_bath, INTERFACE &nanoparticle, CONTROL &cpmdremote) {


    double potential_energy = energy_functional(s, ion, nanoparticle);

    if (world.rank() == 0) {
        ofstream list_tic("outfiles/total_induced_charge.dat", ios::app);
        ofstream list_temperature("outfiles/temperature.dat", ios::app);
        ofstream list_energy("outfiles/energy.dat", ios::app);
        list_temperature << cpmdstep << setw(15) << 2 * particle_kinetic_energy(ion) / (real_bath[0].dof * kB)
                         << setw(15)
                         << real_bath[0].T << setw(15) << 2 * fake_kinetic_energy(s) / (fake_bath[0].dof * kB)
                         << setw(15)
                         << fake_bath[0].T << endl;
        list_tic << cpmdstep << setw(15) << nanoparticle.total_induced_charge(s) << endl;
        double fake_ke = fake_kinetic_energy(s);
        double particle_ke = particle_kinetic_energy(ion);
        double real_bath_ke = bath_kinetic_energy(real_bath);
        double real_bath_pe = bath_potential_energy(real_bath);
        double fake_bath_ke = bath_kinetic_energy(fake_bath);
        double fake_bath_pe = bath_potential_energy(fake_bath);
        double extenergy =
                fake_ke + particle_ke + potential_energy + real_bath_ke + real_bath_pe + fake_bath_ke + fake_bath_pe;
        list_energy << cpmdstep << setw(15) << extenergy << setw(15) << particle_ke << setw(15) << potential_energy
                    << setw(15) << particle_ke + potential_energy + real_bath_ke + real_bath_pe << setw(15) << fake_ke
                    << setw(15) << fake_ke + fake_bath_ke + fake_bath_pe << setw(15) << real_bath_ke << setw(15)
                    << real_bath_pe << setw(15) << fake_bath_ke << setw(15) << fake_bath_pe << endl;

    }
}

// verify on the fly properties with exact
double
verify_with_FMD(int cpmdstep, vector<VERTEX> s, vector<PARTICLE> &ion, INTERFACE &nanoparticle, CONTROL &fmdremote,
                CONTROL &cpmdremote) {


    vector<VERTEX> exact_s;
    exact_s = s;
    fmdremote.verify = cpmdstep;
    fmd(exact_s, ion, nanoparticle, fmdremote, cpmdremote);
    for (unsigned int k = 0; k < s.size(); k++)
        exact_s[k].w = exact_s[k].wmean;
    double exact_functional = energy_functional(exact_s, ion, nanoparticle);
    double on_the_fly_functional = energy_functional(s, ion, nanoparticle);
    double functional_deviation = 0;
    for (unsigned int k = 0; k < s.size(); k++)
        functional_deviation =
                functional_deviation + 100 * (on_the_fly_functional - exact_functional) / exact_functional;
    functional_deviation = functional_deviation / s.size();
    if (world.rank() == 0) {
        ofstream track_density("outfiles/track_density.dat", ios::app);
        ofstream track_functional("outfiles/track_functional.dat", ios::app);
        ofstream track_functional_deviation("outfiles/track_deviation.dat", ios::app);
        track_density << cpmdstep << setw(15) << s[0].w << setw(15) << exact_s[0].w << endl;
        track_functional << cpmdstep << setw(15) << on_the_fly_functional << setw(15) << exact_functional << endl;
        track_functional_deviation << cpmdstep << setw(15) << functional_deviation << endl;

        // write exact induced density
        char data[200];
        sprintf(data, "verifiles/_ind_%.06d.dat", cpmdstep);
        ofstream out_correct_ind;
        out_correct_ind.open(data);
        for (unsigned int k = 0; k < s.size(); k++)
            out_correct_ind << k + 1 << " " << exact_s[k].theta << " " << exact_s[k].phi << " " << exact_s[k].wmean
                            << endl;

        // write cpmd computed induced density
        sprintf(data, "computedfiles/_cpmdind_%.06d.dat", cpmdstep);
        ofstream out_cpmd_ind;
        out_cpmd_ind.open(data);
        for (unsigned int k = 0; k < s.size(); k++)
            out_cpmd_ind << k + 1 << " " << s[k].theta << " " << s[k].phi << " " << s[k].w << endl;
    }
    return functional_deviation;
}

// make movie
void make_movie(int num, vector<PARTICLE> &ion, INTERFACE &nanoparticle, CONTROL &cpmdremote) {

    if (world.rank() == 0) {

        std::string ions_pos_str = "\n####_Ions_Position_Wrapper_####\n";

        ofstream outdump("outfiles/p.lammpstrj", ios::app);
        outdump << "ITEM: TIMESTEP" << endl;
        outdump << num - 1 << endl;
        outdump << "ITEM: NUMBER OF ATOMS" << endl;
        outdump << ion.size() << endl;
        outdump << "ITEM: BOX BOUNDS" << endl;
        outdump << -nanoparticle.box_radius << "\t" << nanoparticle.box_radius << endl;
        outdump << -nanoparticle.box_radius << "\t" << nanoparticle.box_radius << endl;
        outdump << -nanoparticle.box_radius << "\t" << nanoparticle.box_radius << endl;
        outdump << "ITEM: ATOMS index type x y z" << endl;
        string type;
        for (unsigned int i = 0; i < ion.size(); i++) {
            if (ion[i].valency > 0)
                type = "1";
            else
                type = "-1";
            outdump << setw(6) << i << "\t" << type << "\t" << setw(8) << ion[i].posvec.x << "\t" << setw(8)
                    << ion[i].posvec.y << "\t" << setw(8) << ion[i].posvec.z << endl;
            if (!cpmdremote.verbose)
                ions_pos_str =
                        ions_pos_str + std::to_string(i) + "," + type + "," + std::to_string(ion[i].posvec.x) +
                        "," + std::to_string(ion[i].posvec.y) + "," + std::to_string(ion[i].posvec.z) + "\n";
        }
        outdump.close();
        ions_pos_str = ions_pos_str + "####_Ions_Position_Wrapper__Over_####";
        if (!cpmdremote.verbose)
            cout << ions_pos_str << "\n";

    }
    return;
}

// compute density profile
void compute_density_profile(int cpmdstep, double density_profile_samples, vector<double> &mean_density,
                             vector<double> &mean_sq_density, vector<PARTICLE> &ion, INTERFACE &nanoparticle,
                             vector<BIN> &bin, CONTROL &cpmdremote) {


    vector<double> sample_density;


    if (world.rank() == 0) {
        ofstream file_for_auto_corr("outfiles/for_auto_corr.dat", ios::app);

        bin_ions(ion, nanoparticle, sample_density, bin);
        for (unsigned int b = 0; b < mean_density.size(); b++)
            mean_density.at(b) = mean_density.at(b) + sample_density.at(b);
        for (unsigned int b = 0; b < sample_density.size(); b++)
            mean_sq_density.at(b) = mean_sq_density.at(b) + sample_density.at(b) * sample_density.at(b);

        // write a file for post analysis to get auto correlation time		// NOTE this is assuming ions do not cross the interface
        if (ion[0].posvec.GetMagnitude() > nanoparticle.radius)
            file_for_auto_corr << cpmdstep - cpmdremote.hiteqm << "\t"
                               << sample_density[int((nanoparticle.radius + 2) / bin[0].width)] << endl;
        else
            file_for_auto_corr << cpmdstep - cpmdremote.hiteqm << "\t"
                               << sample_density[int((nanoparticle.radius - 2) / bin[0].width)] << endl;

        // write files
        if (cpmdstep % cpmdremote.writedensity == 0) {

            std::string desnity_pr_str = "\n####_Density_Profile_Wrapper_####";
            desnity_pr_str = desnity_pr_str + "#Stepsize:" + std::to_string(cpmdstep) + "\n";

            char data[200];
            sprintf(data, "datafiles/_den_%.06d.dat", cpmdstep);
            ofstream outden;
            outden.open(data);
            for (unsigned int b = 0; b < mean_density.size(); b++) {
                outden << b * bin[b].width << setw(15) << mean_density.at(b) / density_profile_samples << endl;
                if (!cpmdremote.verbose)
                    desnity_pr_str = desnity_pr_str + std::to_string(b * bin[b].width) + "," +
                                     std::to_string(mean_density.at(b) / density_profile_samples) + "\n";
            }
            outden.close();
            desnity_pr_str = desnity_pr_str + "####_Density_Profile_Wrapper_Over_####";
            if (!cpmdremote.verbose) {
                cout << desnity_pr_str << "\n";
                cout_energy_data();
            }
        }


    }
    return;
}

void cout_energy_data() {

    char filename[200];
    sprintf(filename, "outfiles/energy.dat");
    ifstream in(filename, ios::in);
    if (!in) {
        if (world.rank() == 0)
            cout << "File could not be opened" << endl;
        return;
    }
    if (world.rank() == 0) {

        std::string energy_pr_str = "\n####_Energy_Profile_Wrapper_####\n";

        int col1;
        double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
        while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11) {
            energy_pr_str =
                    energy_pr_str + std::to_string(col1) + "," + std::to_string(col2) + "," + std::to_string(col3) +
                    "," + std::to_string(col4) + "," + std::to_string(col6) + "\n";
        }
        energy_pr_str = energy_pr_str + "####_Energy_Profile_Wrapper_Over_####";
        cout << energy_pr_str << "\n";
    }
}


// compute MD trust factor R
double compute_MD_trust_factor_R(int hiteqm) {


    char filename[200];
    sprintf(filename, "outfiles/energy.dat");
    ifstream in(filename, ios::in);
    if (!in) {
        if (world.rank() == 0)
            cout << "File could not be opened" << endl;
        return 0;
    }

    int col1;
    double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
    vector<double> ext, ke, pe, fake, real;
    while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11) {
//     if (col1 < hiteqm) continue;
        ext.push_back(col2);
        ke.push_back(col3);
        pe.push_back(col4);
        fake.push_back(col6);
        real.push_back(col5);
    }
//   cout << "Sizes of ext and ke arrays" << setw(10) << ext.size() << setw(10) << ke.size() << endl;

    double ext_mean = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
        ext_mean += ext[i];
    ext_mean = ext_mean / ext.size();
    double ke_mean = 0;
    for (unsigned int i = 0; i < ke.size(); i++)
        ke_mean += ke[i];
    ke_mean = ke_mean / ke.size();

//   cout << "Mean ext and ke" << setw(10) << ext_mean << setw(10) << ke_mean << endl;

    double ext_sd = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
        ext_sd += (ext[i] - ext_mean) * (ext[i] - ext_mean);
    ext_sd = ext_sd / ext.size();
    ext_sd = sqrt(ext_sd);

    double ke_sd = 0;
    for (unsigned int i = 0; i < ke.size(); i++)
        ke_sd += (ke[i] - ke_mean) * (ke[i] - ke_mean);
    ke_sd = ke_sd / ke.size();
    ke_sd = sqrt(ke_sd);


    double R = ext_sd / ke_sd;
//   cout << "R" << setw(15) <<  R << endl;

    if (world.rank() == 0) {
        ofstream out("outfiles/R.dat");
        out << "Sample size " << ext.size() << endl;
        out << "Sd: ext, kinetic energy and R" << endl;
        out << ext_sd << setw(15) << ke_sd << setw(15) << R << endl;
    }
    return R;
}


// compute MD trust factor R_v
double compute_MD_trust_factor_R_v(int hiteqm) {


    char filename[200];
    sprintf(filename, "outfiles/energy.dat");
    ifstream in(filename, ios::in);
    if (!in) {
        if (world.rank() == 0)
            cout << "File could not be opened" << endl;
        return 0;
    }

    int col1;
    double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
    vector<double> ext, fake;
    while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11) {

        ext.push_back(col2);
        fake.push_back(col6);

    }
    double ext_mean = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
        ext_mean += ext[i];
    ext_mean = ext_mean / ext.size();
    double ke_mean_fake = 0;
    for (unsigned int i = 0; i < fake.size(); i++)
        ke_mean_fake += fake[i];
    ke_mean_fake = ke_mean_fake / fake.size();

//   cout << "Mean ext and ke" << setw(10) << ext_mean << setw(10) << ke_mean << endl;

    double ext_sd = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
        ext_sd += (ext[i] - ext_mean) * (ext[i] - ext_mean);
    ext_sd = ext_sd / ext.size();
    ext_sd = sqrt(ext_sd);

    double ke_sd_fake = 0;
    for (unsigned int i = 0; i < fake.size(); i++)
        ke_sd_fake += (fake[i] - ke_mean_fake) * (fake[i] - ke_mean_fake);
    ke_sd_fake = ke_sd_fake / fake.size();
    ke_sd_fake = sqrt(ke_sd_fake);

//   cout << "Standard deviations in ext and ke" << setw(10) << ext_sd << setw(10) << ke_sd << endl;

    double RV = ext_sd / ke_sd_fake;
//   cout << "RV" << setw(15) <<  RV << endl;

    if (world.rank() == 0) {
        ofstream out("outfiles/RV.dat");
        out << "Sample size " << ext.size() << endl;
        out << "Sd: ext, kinetic energy_fake and RV" << endl;
        out << ext_sd << setw(15) << ke_sd_fake << setw(15) << RV << endl;
    }
    return RV;
}


void progressBar(double fraction_completed) {


    if (world.rank() == 0) {
        int val = (int) (fraction_completed * 100);
        int lpad = (int) (fraction_completed * PBWIDTH);
        int rpad = PBWIDTH - lpad;
        printf("\r%3d%% |%.*s%*s|", val, lpad, PBSTR, rpad, "");
        fflush(stdout);
    }
}

/*
// auto correlation function
void auto_correlation_function()
{
  char filename[200];
  sprintf(filename, "outfiles/for_auto_corr.dat");
  ifstream in(filename, ios::in);
  if (!in) 
  {
    cout << "File could not be opened" << endl; 
    return;
  }

  double col1, col2;
  vector<double> n, autocorr;
  while (in >> col1 >> col2)
    n.push_back(col2);
//   cout << "Number of samples after eqm" << setw(10) << n.size() << endl;

  double avg = 0;
  for (unsigned int j = 0; j< n.size(); j++)
    avg = avg + n[j];
  avg = avg / n.size();

  int ntau = 5000;			// time to which the auto correlation function is computed 
  
  for (int i = 0; i < ntau; i++)
  {
    double A = 0;
    for (unsigned int j = 0; j< n.size(); j++)
      A = A + n[j+i]*n[j];
    A = A / n.size();
    autocorr.push_back(A - avg*avg);
  }

  ofstream out ("outfiles/auto_correlation.dat");
  for (int i = 0; i < ntau; i++)
    out << i << setw(15) << autocorr[i]/autocorr[0] << endl;
  
  cout << "Auto correlation function generated" << endl;
  return;
}
*/
