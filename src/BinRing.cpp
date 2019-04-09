#include "BinRing.h"

// make a bin
BinRing::BinRing(){}

BinRing::BinRing(int nL = 0, int bin_num_RL = 0, int bin_num_ZL = 0, double width_RL = 0, double height_RL = 0,
                 double width_ZL = 0, double height_ZL = 0, double volumeL = 0, double lower_RL = 0,
                 double higher_RL = 0, double lower_ZL = 0, double higher_ZL = 0) {

    n = nL;
    bin_num_R = bin_num_RL;
    bin_num_Z = bin_num_ZL;
    width_R = width_RL;
    width_Z = height_RL;
    volume = width_ZL;
    lower_R = height_ZL;
    higher_R = volumeL;
    lower_Z = lower_RL;
    higher_Z = higher_RL;
}


// set up bin in 2D derection
void BinRing::set_up(int bin_num_Z_, int bin_num_R_, double bin_width_Z_, double bin_width_R) {

    n = 0;
    bin_num_R = bin_num_R_;
    bin_num_Z = bin_num_Z_;
    width_R = bin_width_R;
    width_Z = bin_width_Z_;

    lower_R = bin_num_R_ * width_R;
    higher_R = (bin_num_R_ + 1) * width_R;
    lower_Z = bin_num_Z_ * width_Z;
    higher_Z = (bin_num_Z_ + 1) * width_Z;

    volume = pi * (higher_Z - lower_Z) * (pow(higher_R, 2) - pow(lower_R, 2));

}



