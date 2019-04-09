// This is bin class

#ifndef _BIN_H
#define _BIN_H

#include "utility.h"

class BinShell
{

  private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & n;
        ar & width;
        ar & volume;
        ar & lower;
        ar & higher;
    }

  public:
    
    // members
    int n;		// number of ions in the bin
    double width;	// width of the bin
    double volume;	// volume of the bin
    double lower;	// lower value of bin
    double higher;	// higher value of bin
    
    // member functions
    // make a bin
    BinShell();

    BinShell(int, double, double, double, double);
    
    // set up bin
    void set_up(int, double);

};

#endif

