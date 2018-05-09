// This is a header file to set control

#ifndef _CONTROL_H
#define _CONTROL_H

class CONTROL
{
  private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & verbose;
        ar & fakemass;
        ar & timestep;
        ar & steps;
        ar & hiteqm;
        ar & freq;
        ar & on;
        ar & anneal;
        ar & annealfreq;
        ar & extra_compute;
        ar & verify;
        ar & writeverify;
        ar & writedensity;
        ar & writedata;
    }

  public:
    bool verbose;   //more control for nanoHUB deployment
    double fakemass;		// fictitious mass used in molecular dynamics
    double timestep;		// timestep used in molecular dynamics
    int steps;			// number of steps in molecular dynamics
    int hiteqm; 		// wait till this step
    int freq; 			// frequency of sampling
    char on;			// y if fmd needs to turn on
    char anneal;		// if you want to anneal, applies only to fmd
    int annealfreq;		// frequency of annealing
    int extra_compute;		// energy computed after these many steps
    int verify;			// verification with correct answer for ind. density done after these many steps
    int writeverify;		// write the verification files
    int writedensity; 		// write the density files
    int writedata; 		// write the data files
};

#endif
