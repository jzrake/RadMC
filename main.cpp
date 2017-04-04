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
    template <class RandomVariableType, class EngineType>
    static Electron generateIsotropic (RandomVariableType& randomGammaBeta, EngineType& engine)
    {
        return generateIsotropic (randomGammaBeta (engine), engine);
    }

    template <class EngineType>
     static Electron generateIsotropic (double gammaBeta, EngineType& engine)
    {
        const UnitVector nhat = UnitVector::generateIsotropic (engine);
        return FourVector::fromGammaBetaAndUnitVector (gammaBeta, nhat);
    }

    Electron () : momentum (1, 0, 0, 0) {}
    Electron (FourVector momentum) : momentum (momentum) {}

    FourVector momentum;
};



class Photon
{
public:
    template <class RandomVariableType, class EngineType>
    static Photon generateIsotropic (RandomVariableType& randomEnergy, EngineType& engine)
    {
        return generateIsotropic (randomEnergy (engine), engine);
    }

    template <class EngineType>
     static Photon generateIsotropic (double energy, EngineType& engine)
    {
        const UnitVector nhat = UnitVector::generateIsotropic (engine);
        return FourVector::nullWithUnitVector (nhat) * energy;
    }

    Photon (FourVector momentum) : momentum (momentum) {}

    FourVector momentum;
};



void doComptonScattering (Photon& photon, Electron& electron)
{
    std::mt19937 engine;

    UnitVector finalPhotonDir = UnitVector::generateIsotropic (engine);

    double cosTheta = finalPhotonDir.pitchAngleMu;
    double e0 = photon.momentum.getTimeComponent();
    double e1 = e0 / (1 + e0 * (1 - cosTheta));

    FourVector p1 = FourVector::nullWithUnitVector (finalPhotonDir) * e1;
    FourVector q1 = electron.momentum + photon.momentum - p1;

    photon.momentum = p1;
    electron.momentum = q1;
}



Electron generateElectronGivenPhotonWasScattered (const Photon& photon)
{
    std::mt19937 engine;

    // 0. The photon is measured in the lab frame.

    // 1. Choose a temperature for the electron distribution.
    double kT = 1.0;
    auto maxwellPdf = Distributions::makeMaxwellBoltzmann (kT, Distributions::probabilityDensityFunction);
    RandomVariable::FromProbabilityDensityFunction electronBetaRv (maxwellPdf, 0, 5);

    // 2. Sample the electron speed (Note: convert from gammaBeta if MJ).
    double electronBeta = electronBetaRv (engine);

    // 3. Create a distribution of pitch angles over which the scattering
    //    might have occurred.
    auto scatteringMuQnt = Distributions::makePitchAngleGivenScattered (electronBeta, Distributions::quantileFunction);
    RandomVariable::FromQuantileFunction scatteringMuRv (scatteringMuQnt);
    std::uniform_real_distribution<double> scatteringPhiRv (0, 2 * M_PI);

    // 4. Sample that distribution. This mu is the angle between the photon
    //    and electron propagation, measured in the lab frame. Also generate a
    //    uniform azimuthal angle.
    double scatteringMu = scatteringMuRv (engine);
    double scatteringPhi = scatteringPhiRv (engine);

    // 5. Form the unit vector nhatPrime of the electron, that is the
    //    electron's propagation in a coordinate system where zhat is aligned
    //    with the photon momentum.
    UnitVector nhatPrime (scatteringMu, scatteringPhi);

    // 6. Rotate the electron propagation direction so that its pitch angle
    //    relative to the photon's is mu.
    UnitVector nhat = nhatPrime.withPolarAxis (photon.momentum.getUnitThreeVector());

    std::cout << scatteringMu <<std::endl;
    std::cout << nhat.getPitchAngleWithRespectTo (photon.momentum.getUnitThreeVector()) << std::endl;

    std::cout << scatteringPhi <<std::endl;
    std::cout << nhat.azimuthalAnglePhi << std::endl;

    return FourVector::fromBetaAndUnitVector (electronBeta, nhat);
}


int main (int argc, char **argv)
{

    std::mt19937 engine;

    Photon p = Photon::generateIsotropic (1.0, engine);
    Electron e = Electron::generateIsotropic (2.0, engine);

    generateElectronGivenPhotonWasScattered (p);
    return 0;

    std::cout << p.momentum + e.momentum << std::endl;

    doComptonScattering (p, e);

    std::cout << p.momentum + e.momentum << std::endl;
    
    return 0;
}








int testHistogram()
{
    // auto pdf = Distributions::makePitchAngleGivenScattered (0.4, Distributions::probabilityDensityFunction);
    // auto qnt = Distributions::makePitchAngleGivenScattered (0.4, Distributions::quantileFunction);

    // DistributionSampler sampler;
    // //sampler.setQuantileFunction (qnt);
    // sampler.computeQuantileFunctionFromDensity (pdf, -1, 1);

    // std::vector<double> samples = sampler.generateSamples (1 << 20);

    // TabulatedFunction hist = TabulatedFunction::makeHistogram (samples, 128,
    //     TabulatedFunction::useEqualBinWidthsLinear, true, true, true);

    // std::cout.flags (std::ios::fixed | std::ios::showpos);
    // std::cout.precision (10);
    // hist.outputTable (std::cout);

    return 0;
}


