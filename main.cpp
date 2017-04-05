#include <iostream>
#include <fstream>
#include <cmath>
#include "FourVector.hpp"
#include "LorentzBoost.hpp"
#include "QuadratureRule.hpp"
#include "ProbabilityDistribution.hpp"
#include "TabulatedFunction.hpp"
#include "CubicInterpolant.hpp"
#include "Distributions.hpp"
#include "RandomVariable.hpp"




class Electron
{
public:
    Electron () : momentum (1, 0, 0, 0) {}
    Electron (FourVector momentum) : momentum (momentum) {}

    FourVector momentum;
};




class Photon
{
public:
    Photon (FourVector momentum) : momentum (momentum) {}

    FourVector momentum;
};




void doComptonScattering (Photon& photon, Electron& electron)
{
    RandomVariable crossSection = RandomVariable::uniformOver (-1, 1);
    UnitVector finalPhotonDir = photon.momentum.getUnitThreeVector().sampleAxisymmetric (crossSection);

    double cosTheta = finalPhotonDir.pitchAngleMu;
    double e0 = photon.momentum.getTimeComponent();
    double e1 = e0 / (1 + e0 * (1 - cosTheta));

    FourVector p1 = FourVector::nullWithUnitVector (finalPhotonDir) * e1;
    FourVector q1 = electron.momentum + photon.momentum - p1;

    photon.momentum = p1;
    electron.momentum = q1;
}


Electron generateElectronGivenPhotonWasScattered (const Photon& photon, RandomVariable& electronGammaBeta)
{
    // Sample the electron speed (Note: convert from gammaBeta if MJ).
    double electronU = electronGammaBeta.sample();
    double electronV = FourVector::betaFromGammaBeta (electronU);

    // Sample the electron propagation direction. The distribution of electron
    // velocities, given that scattering has occurred, is axisymmetric around
    // the photon propagation vector k, but heavier around -k due to the
    // relative motion.
    RandomVariable relativeMu (Distributions::makePitchAngle (electronV, Distributions::Qnt));
    UnitVector nhat = photon.momentum.getUnitThreeVector().sampleAxisymmetric (relativeMu);

    return FourVector::fromBetaAndUnitVector (electronV, nhat);
}




int main (int argc, char **argv)
{
    double kT = 1.0;
    auto maxwellPdf = Distributions::makeMaxwellBoltzmann (kT, Distributions::Pdf);

    RandomVariable electronBetaRv (RandomVariable::fromPdf (maxwellPdf, 0, 5));


    return 0;
}








int testHistogram()
{


    //auto muQnt = Distributions::makePitchAngle (0.5, Distributions::Qnt);

    // auto pdf = Distributions::makePitchAngle (0.4, Distributions::Pdf);
    // auto qnt = Distributions::makePitchAngle (0.4, Distributions::Qnt);

    // DistributionSampler sampler;
    // //sampler.setQnt (qnt);
    // sampler.computeQntFromDensity (pdf, -1, 1);

    // std::vector<double> samples = sampler.generateSamples (1 << 20);

    // TabulatedFunction hist = TabulatedFunction::makeHistogram (samples, 128,
    //     TabulatedFunction::useEqualBinWidthsLinear, true, true, true);

    // std::cout.flags (std::ios::fixed | std::ios::showpos);
    // std::cout.precision (10);
    // hist.outputTable (std::cout);

    return 0;
}


