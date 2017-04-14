#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <sys/stat.h> // mkdir
#include "FourVector.hpp"
#include "LorentzBoost.hpp"
#include "QuadratureRule.hpp"
#include "ProbabilityDistribution.hpp"
#include "TabulatedFunction.hpp"
#include "CubicInterpolant.hpp"
#include "Distributions.hpp"
#include "RandomVariable.hpp"
#include "RichardsonCascade.hpp"
#include "Variant.hpp"




class Electron
{
public:
    Electron () : momentum (1, 0, 0, 0) {}
    Electron (FourVector momentum) : momentum (momentum) {}

    /**
    Return the electron four-velocity (this is the same as its momentum since
    electron mass is 1).
    */
    FourVector getFourVelocity()
    {
        return momentum;
    }

    FourVector momentum;
};




class Photon
{
public:
    static Photon sampleIsotropic (RandomVariable& photonEnergy)
    {
        double E = photonEnergy.sample();
        return FourVector::nullWithUnitVector (UnitVector::sampleIsotropic()) * E;
    }
    Photon (FourVector momentum) : momentum (momentum) {}

    FourVector momentum;
};




void doComptonScattering (Photon& photon, Electron& electron)
{
    // Photon four-momentum in the electron rest frame
    LorentzBoost L (electron.getFourVelocity());
    FourVector p0 = photon.momentum.transformedBy (L);

    // Photon scattering angle and direction (k1) in the electron rest frame
    RandomVariable crossSection = RandomVariable::uniformOver (-1, 1); // Thomson is ~ isotropic
    double cosTheta = crossSection.sample();
    UnitVector k1 = p0.getUnitThreeVector().sampleAxisymmetric (cosTheta);

    // Photon energies in the electron frame, before and after scattering
    double e0 = p0.getTimeComponent();
    double e1 = e0 / (1 + e0 * (1 - cosTheta));

    // Photon four-momentum seen in the electron's initial rest frame (after
    // scattering) and the associated impulse - dp - seen in the lab frame.
    FourVector p1 = FourVector::nullWithUnitVector (k1) * e1;
    FourVector dp = (p1 - p0).transformedBy (L.inverted());

    // Change the photon and electron energy-momentum vector.
    photon.momentum   += dp;
    electron.momentum -= dp;
}




Electron sampleElectronForScattering (const Photon& photon, RandomVariable& electronGammaBeta)
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




int comptonize (int argc, const char *argv[])
{
    Variant::NamedValues userParams;

    userParams["kT"] = 0.01;
    userParams["Ephot"] = 0.1;
    userParams["nphot"] = 1000000;
    userParams["iter"] = 10;

    Variant::updateFromCommandLine (userParams, argc - 1, argv + 1);

    std::cout << userParams << std::endl;

    double kT = userParams["kT"];
    int nphot = userParams["nphot"];
    double uphot = userParams["Ephot"];

    auto maxwellPdf = Distributions::makeMaxwellBoltzmann (kT, Distributions::Pdf);

    RandomVariable photonEnergy = RandomVariable::diracDelta (uphot);
    RandomVariable electronGammaBeta (RandomVariable::fromPdf (maxwellPdf, 0, 5 * kT));

    std::vector<Photon> photons;

    for (int i = 0; i < nphot; ++i)
    {
        photons.push_back (Photon::sampleIsotropic (photonEnergy));
    }

    for (int n = 0; n < int (userParams["iter"]); ++n)
    {
        for (auto it = photons.begin(); it != photons.end(); ++it)
        {
            if (RandomVariable::sampleUniform() < 1.0)
            {
                Electron e = sampleElectronForScattering (*it, electronGammaBeta);
                doComptonScattering (*it, e);
            }
        }
        
        std::cout << n << std::endl;

        if (n % 1 == 0)
        {
            std::vector<double> samples;

            std::for_each (photons.begin(), photons.end(), [&] (Photon& p)
            {
                samples.push_back (p.momentum.getTimeComponent());
            });

            TabulatedFunction hist = TabulatedFunction::makeHistogram (samples, 256,
                TabulatedFunction::useEqualBinWidthsLinear, true, true, false);

            std::stringstream stream;
            stream << "compton" << n << ".dat";
            std::ofstream out (stream.str());
            hist.outputTable (out);
        }
    }

    return 0;
}





// ============================================================================
class TurbulenceCascadeDriver
{
public:
    TurbulenceCascadeDriver()
    {

    }

    void run (int argc, const char *argv[])
    {
        makeUserParameters();

        try
        {
            Variant::updateFromCommandLine (userParams, argc - 1, argv + 1);
            std::cout << "=====================================================\n";
            std::cout << userParams;
            std::cout << "=====================================================\n";
        }
        catch (std::runtime_error& error)
        {
            std::cout << error.what() << std::endl;
            return;
        }

        while (shouldContinue())
        {
            advance();
        }
    }

    bool shouldContinue()
    {
        return simulationTime < double (userParams["tmax"]);
    }

    void advance()
    {
        if (simulationIter % 100 == 0)
        {
            std::ofstream out = makeOutputStream();
            cascade.spectralEnergy.outputTable (out);
        }

        double dt = 0.5 * cascade.getShortestTimeScale();
        cascade.advance (dt);
        simulationTime += dt;
        simulationIter += 1;
    }

    void makeUserParameters()
    {
        userParams["tmax"] = 1.0;
        userParams["outdir"] = ".";
        userParams["Re"] = 1e6;
    }

private:
    std::ofstream makeOutputStream()
    {
        std::string outdir = userParams["outdir"];
        mkdir (outdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::string filename = outdir + "/cascade" + std::to_string (simulationIter / 100) + ".dat";
        std::ofstream out (filename);
        return out;
    }

    Variant::NamedValues userParams;
    RichardsonCascade cascade;
    double simulationTime;
    int simulationIter;
};




int main (int argc, const char *argv[])
{
    TurbulenceCascadeDriver driver;
    driver.run (argc, argv);

    return 0;
}
