#include "fftw++.h"

unsigned int fftw::effort=FFTW_PATIENT;
bool fftw::Wise=false;
const char *fftw::WisdomName="wisdom3.txt";
ifstream fftw::ifWisdom;
ofstream fftw::ofWisdom;
