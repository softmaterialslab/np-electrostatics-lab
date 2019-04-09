// This is bin class

#ifndef _Disk_BIN_H
#define _Disk_BIN_H

#include "utility.h"

class BinRing {

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & n;
        ar & bin_num_R;
        ar & bin_num_Z;
        ar & width_R;
        ar & width_Z;
        ar & volume;
        ar & lower_R;
        ar & higher_R;
        ar & lower_Z;
        ar & higher_Z;
    }

public:

    // members
    int n;              // number of ions in the bin
    int bin_num_R;    // bin index for R direction
    int bin_num_Z;    //  bin index for Z direction
    double width_R;    // width of the bin in R direction
    double width_Z;    // width of the bin in Z direction
    double volume;    // volume of the bin
    double lower_R;    // lower value of bin in R direction
    double higher_R;    // higher value of bin in R direction
    double lower_Z;    // lower value of bin in Z direction
    double higher_Z;    // higher value of bin in Z direction

    // member functions
    // make a bin
    BinRing();

    BinRing(int, int, int, double, double, double, double, double, double, double, double, double);

    // set up bin in 2D
    void set_up(int , int , double , double  );

};

#endif

