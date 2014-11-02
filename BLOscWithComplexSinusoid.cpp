#include "SC_PlugIn.h"
#include <tgmath.h>
#include <complex>

/////////////////////////////////////////////////////////////////

// Band Limited Oscillator made with complex sinusoid, without wave table.

/////////////////////////////////////////////////////////////////

// BASIC ADMINISTRATION

static InterfaceTable *ft;

struct BLOscWithComplexSinusoid : public Unit
{
    float partialTheta; // 2pi / sr
	float currentphase;
};

extern "C"
{
	void BLOscWithComplexSinusoid_next(BLOscWithComplexSinusoid *unit, int inNumSamples);
	void BLOscWithComplexSinusoid_Ctor(BLOscWithComplexSinusoid* unit);
};

//////////////////////////////////////////////////////////////////

// CONSTRUCTOR

void BLOscWithComplexSinusoid_Ctor(BLOscWithComplexSinusoid* unit)
{

	SETCALC(BLOscWithComplexSinusoid_next);

	unit->partialTheta = twopi_f / SAMPLERATE;
	unit->currentphase = twopi_f * 0.75; // output is cosine wave. Set the initial phase to 0.75pi to make the beginning displacement 0 and avoid clicking noise

	BLOscWithComplexSinusoid_next(unit, 1);
}

//////////////////////////////////////////////////////////////////

// UGEN CALCULATION

void BLOscWithComplexSinusoid_next(BLOscWithComplexSinusoid *unit, int inNumSamples)
{
    
    float *out = ZOUT(0);
    float freqin = ZIN0(0);
    int32 loHarmonics = ZIN0(1);
    int32 numHarmonics = ZIN0(2);
    float slope = ZIN0(3);
    float evenOddRatio = ZIN0(4);

    int32 hiHarmonics = loHarmonics + numHarmonics - 1; // The highest harmonic index
    int32 hiHarmonicsPlusOne = hiHarmonics + 1;
    int32 loEvenHarmonics = loHarmonics%2 == 0? loHarmonics : loHarmonics + 1; // The lowest even harmonic index
    int32 hiEvenHarmonics = hiHarmonics%2 == 0? hiHarmonics : hiHarmonics - 1; // The highest even harmonic index
    int32 hiEvenHarmonicsPlusTwo = hiEvenHarmonics + 2;
    int32 numEvenHarmonics = (hiEvenHarmonics - loEvenHarmonics) / 2 + 1; //The total number of even harmonics
    float evenOddFactor = 1 - evenOddRatio;
    std::complex<double> evenOddFactorC(evenOddFactor, 0);
    float ampFactor = 0.99<slope&&slope<1.01? numHarmonics - evenOddFactor * numEvenHarmonics:((pow(slope,loHarmonics) - pow(slope,hiHarmonicsPlusOne)) / (1 - slope)) - (evenOddFactor * (pow(slope, loEvenHarmonics) - pow(slope, hiEvenHarmonicsPlusTwo)) / (1 - pow(slope, 2)));  //ampFactor will be used to normalize the output amplitude. To avoid the denominator of this calculation to be 0 when slope = 1, the different formula is used when slope falls between 0.99 and 1.01.
    
    float phaseinc = ZIN0(0) * unit->partialTheta;
    float currentphase = unit->currentphase;
    std::complex<double> baseOsc;
    std::complex<double> r;
    std::complex<double> signalC; // signal as complex number
    float signalR; // signal as real number
    std::complex<double> oneC(1,0);
    std::complex<double> slopeC(slope, 0);
    float z;
    
    LOOP(inNumSamples,
         baseOsc = std::exp(std::complex<double>(0, currentphase));

         r = slopeC * baseOsc;
         signalC = (std::pow(r, loHarmonics) - std::pow(r, hiHarmonicsPlusOne))/(oneC - r) - evenOddFactorC * (std::pow(r, loEvenHarmonics) - std::pow(r, hiEvenHarmonicsPlusTwo))/(oneC - std::pow(r, 2));
         signalR = signalC.real();
         z = signalR/ampFactor;
         
         currentphase += phaseinc;
         while (currentphase >= twopi_f)
         currentphase -= twopi_f;
         
         ZXP(out) = z;
         )
    unit->currentphase = currentphase;
}

////////////////////////////////////////////////////////////////////

// LOAD FUNCTION

PluginLoad(BLOscWithComplexSinusoid)
{
	ft = inTable;

	DefineSimpleUnit(BLOscWithComplexSinusoid);
}

////////////////////////////////////////////////////////////////////

