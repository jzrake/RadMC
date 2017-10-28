#include "StructuredJetModel.hpp"
#include "LorentzBoost.hpp"
#include "RelativisticWind.hpp"
#include "TabulatedFunction.hpp"
#include "Distributions.hpp"
#include "ScatteringOperations.hpp"




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
    // auto photon = generatePhoton (innerRadius, theta);
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
    const auto rhat = UnitVector (std::cos (theta), phi); // wind propagation angle
    const auto u = FourVector::fromGammaBetaAndUnitVector (ur, rhat);

    photon.position = FourVector::spaceLikeInDirection (r, rhat);
    photon.momentum = FourVector::nullWithUnitVector (UnitVector::sampleIsotropic()).transformedBy(u);

    const double kT = sampleWind (r, theta).temperature();
    auto photonEnergy = RandomVariable::diracDelta (kT);

    photon.momentum *= photonEnergy.sample();

    return photon;
}

StructuredJetModel::Electron StructuredJetModel::generateElectron (const Photon& photon) const
{
    auto sops = ScatteringOperations();
    auto state = sampleWind (photon.position);
    auto electron = Electron();
    double uth = std::sqrt (3. * state.temperature()); // TODO: properly sample a Maxwellian here
    electron.momentum = sops.sampleScatteredParticles (state.fourVelocity(), photon.momentum, uth);
    return electron;
}

RelativisticWind::WindState StructuredJetModel::sampleWind (double r, double theta) const
{
    const auto x = FourVector::spaceLikeInDirection (r, UnitVector (std::cos (theta), 0.0));
    return sampleWind(x);
}

RelativisticWind::WindState StructuredJetModel::sampleWind (const FourVector& position) const
{
    auto wind = RelativisticWind();
    wind.setSpecificWindPower (10.0);
    wind.setInitialFourVelocity (2.0); // TODO: set wind parameters and structure here
    const double r = position.radius();
    const double t = position.theta();
    auto state = sampleWind (r, t);
    return state.setPropagationAngle (position.getUnitThreeVector());
}
