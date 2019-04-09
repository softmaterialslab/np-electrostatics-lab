#include "BinShell.h"

// make a bin
BinShell::BinShell(){}

BinShell::BinShell(int nL = 0, double widthL = 0, double volumeL = 0, double lowerL = 0, double higherL = 0) {
    n = nL;
    width = widthL;
    volume = volumeL;
    lower = lowerL;
    higher = higherL;
}

// set up bin
void BinShell::set_up(int bin_num, double bin_width) {
    n = 0;
    width = bin_width;
    if (bin_num == 0)
        volume = (4.0 / 3.0) * pi * pow(bin_width, 3);
    else
        volume = 4.0 * pi * 0.25 * (bin_num * bin_width + (bin_num + 1) * bin_width) *
                 (bin_num * bin_width + (bin_num + 1) * bin_width) * bin_width;
// 	volume = 4.0 * pi * (bin_num * bin_width) * (bin_num * bin_width) * bin_width;		NOTE change of the definition of the volume of a bin
    lower = bin_num * bin_width;
    higher = (bin_num + 1) * bin_width;
}


