#include "StructuredJetModel.hpp"
#include "LorentzBoost.hpp"
#include "RelativisticWind.hpp"
#include "TabulatedFunction.hpp"
#include "Distributions.hpp"




// ========================================================================
class SamplingScheme : public RandomVariable::SamplingScheme
{
public:
    SamplingScheme (std::function<double (double)> densityFunction, double x0, double x1)
    {
        const int numberOfTableEntries = 256;
        const double accuracyParameter = 1e-10;
        const GaussianQuadrature gauss (8);

        tabulatedCDF = TabulatedFunction::createTabulatedIntegral (
            densityFunction, x0, x1, numberOfTableEntries,
            TabulatedFunction::useEqualBinWidthsLogarithmic, gauss,
            accuracyParameter, true);
    }

    double generate (double F) override
    {
        return tabulatedCDF.lookupArgumentValue (F);
    }

private:
    TabulatedFunction tabulatedCDF;
};




// ========================================================================
StructuredJetModel::StructuredJetModel()
{

}

std::vector<StructuredJetModel::Photon> StructuredJetModel::generatePhotonPath (double innerRadius, double theta)
{
    auto photon = generatePhoton (innerRadius, theta);
    auto path = std::vector<StructuredJetModel::Photon>();

    for (int n = 0; n < 100; ++n)
    {
        // Electron electron;
        // doComptonScattering (photon, electron);
        // path.push_back (photon);
    }
    return path;
}

StructuredJetModel::Photon StructuredJetModel::generatePhoton (double r, double theta) const
{
    auto photon = Photon();
    const double phi = RandomVariable::sampleUniformAzimuth();
    const double ur = sampleWind (r, theta).u;
    const auto radiusUnitVector = UnitVector (std::cos (theta), phi);
    const auto u = FourVector::fromGammaBetaAndUnitVector (ur, radiusUnitVector);

    photon.position = FourVector::spaceLikeInDirection (r, radiusUnitVector);
    photon.momentum = FourVector::nullWithUnitVector (UnitVector::sampleIsotropic()).transformedBy(u);

    const double kT = sampleWind (r, theta).temperature();
    auto photonEnergy = RandomVariable::diracDelta (kT);

    photon.momentum *= photonEnergy.sample();

    return photon;
}

StructuredJetModel::Electron StructuredJetModel::generateElectron (double r, double theta) const
{

    const double kT = sampleWind (r, theta).temperature();
    auto electronGammaBeta = RandomVariable::diracDelta (kT);

    auto electron = Electron();
    const double phi = RandomVariable::sampleUniformAzimuth();
    const double ur = sampleWind (r, theta).u;
    const auto radiusUnitVector = UnitVector (std::cos (theta), phi);

    const auto ubulk = FourVector::fromGammaBetaAndUnitVector (ur, radiusUnitVector);
    const auto uther = electronGammaBeta.sample();

    electron.momentum = FourVector::fromGammaBetaAndUnitVector (uther, UnitVector::sampleIsotropic());

    return electron.transformedBy(ubulk);
}

RelativisticWind::WindState StructuredJetModel::sampleWind (double r, double theta) const
{
    auto wind = RelativisticWind().setSpecificWindPower (10.0).setInitialFourVelocity (2.0);
    return wind.integrate(r);
}
