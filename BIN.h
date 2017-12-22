// This is bin class

#ifndef _BIN_H
#define _BIN_H

#include "utility.h"

class BIN
{
  public:
    
    // members
    int n;		// number of ions in the bin
    double width;	// width of the bin
    double volume;	// volume of the bin
    double lower;	// lower value of bin
    double higher;	// higher value of bin
    
    // member functions
    
    // make a bin
    BIN(int n = 0, double width = 0, double volume = 0, double lower = 0, double higher = 0 ) : n(n), width(width), volume(volume), lower(lower), higher(higher)
    {
    }
    
    // set up bin
    void set_up(int bin_num, double bin_width)
    {
      n = 0;
      width = bin_width;
      if (bin_num == 0)
	volume = (4.0/3.0) * pi * pow(bin_width, 3);
      else
	volume = 4.0 * pi * 0.25* (bin_num * bin_width + (bin_num+1) * bin_width) * (bin_num * bin_width + (bin_num+1)*bin_width) * bin_width;
// 	volume = 4.0 * pi * (bin_num * bin_width) * (bin_num * bin_width) * bin_width;		NOTE change of the definition of the volume of a bin
      lower = bin_num * bin_width;
      higher = (bin_num + 1) * bin_width;
    }
};

#endif

