#include <iostream>
#include <cassert>
#include "StructuredJetModel.hpp"
#include "LorentzBoost.hpp"
#include "RelativisticWind.hpp"
#include "TabulatedFunction.hpp"
#include "Distributions.hpp"
#include "ScatteringOperations.hpp"
#include "PhysicsConstants.hpp"




// ========================================================================
class SamplingScheme : public RandomVariable::SamplingScheme
{
public:
    SamplingScheme (std::function<double (double)> densityFunction, double x0, double x1)
    {
        const int numberOfTableEntries = 256;
        const double accuracyParameter = 1e-8;
        const GaussianQuadrature gauss(8);

        tabulatedCDF = TabulatedFunction::createTabulatedIntegral (
            densityFunction, x0, x1, numberOfTableEntries,
            TabulatedFunction::useEqualBinWidthsLinear, gauss,
            accuracyParameter, true);
    }
    double generate (double F) override
    {
        return tabulatedCDF.lookupArgumentValue(F);
    }
private:
    TabulatedFunction tabulatedCDF;
};




// ========================================================================
StructuredJetModel::StructuredJetModel (Config config) : config (config)
{
    tabulateWindAllAngles (config.outermostRadius);
}

double StructuredJetModel::sampleTheta() const
{
    auto theta = RandomVariable::uniformOver (0.0, config.jetOpeningAngle * 4);
    return theta.sample();
}

double StructuredJetModel::approximatePhotosphere (double theta) const
{
    auto physics = PhysicsConstants();
    const double G = jetStructureEtaOfTheta (theta);
    const double F = config.luminosityPerSteradian / config.specificWindPower / physics.gramToErg (physics.mp);
    return 0.5 * F * physics.st / physics.c / G / G;
}

std::vector<StructuredJetModel::Photon> StructuredJetModel::generatePhotonPath (double innerRadius, double theta)
{
    auto photon = generatePhoton (innerRadius, theta);
    auto path = std::vector<StructuredJetModel::Photon>();

    while (photon.position.radius() < 1e4)
    {
        std::cout << "[StructuredJetModel] " << ": " << "r = " << photon.position.radius() << std::endl;
        path.push_back (photon = stepPhoton (photon));
    }
    return path;
}

StructuredJetModel::Photon StructuredJetModel::stepPhoton (const StructuredJetModel::Photon& photon0) const
{
    auto scatteringOps = ScatteringOperations();
    auto state = sampleWind (photon0.position);
    auto photon = photon0;

    Electron electron = generateElectron (photon, state);
    scatteringOps.comptonScatter (photon.momentum, electron.momentum);

    // TODO: sample scattering length from exponential
    const double dr = state.thomsonMeanFreePath (photon.momentum.getUnitThreeVector());
    const double dt = dr / config.innerRadiusCm;

    photon.position += photon.momentum / photon.momentum[0] * dt;
    return photon;
}

StructuredJetModel::Photon StructuredJetModel::generatePhoton (double r, double theta) const
{
    auto photon = Photon();
    const double phi = RandomVariable::sampleUniformAzimuth();
    const auto state = sampleWindSpherical (r, theta);
    const auto rhat = UnitVector (std::cos (theta), phi); // wind propagation angle
    const auto u = FourVector::fromGammaBetaAndUnitVector (state.u, rhat);

    const double kT = state.temperature();
    auto photonEnergy = RandomVariable::diracDelta (3 * kT);

    photon.position = FourVector::spaceLikeInDirection (r, rhat);
    photon.momentum = FourVector::nullWithUnitVector (UnitVector::sampleIsotropic());
    photon.momentum *= photonEnergy.sample();
    photon.momentum.transformBy(u);

    return photon;
}

StructuredJetModel::Electron StructuredJetModel::generateElectron (const Photon& photon,
    RelativisticWind::WindState state) const
{
    const double uth = sampleElectronGammaBeta(state.temperature());

    auto sops = ScatteringOperations();
    return sops.sampleScatteredParticlesInFrame (state.fourVelocity(), photon.momentum, uth);
}

RelativisticWind::WindState StructuredJetModel::sampleWindSpherical (double r, double theta) const
{
    const auto x = FourVector::spaceLikeInDirection (r, UnitVector (std::cos (theta), 0.0));
    return sampleWind(x);
}

RelativisticWind::WindState StructuredJetModel::sampleWind (const FourVector& position) const
{
    const double eta = jetStructureEtaOfTheta (position.theta());
    const double r = position.radius();
    const double t = position.theta();

    auto wind = RelativisticWind();
    wind.setSpecificWindPower (eta);
    wind.setInitialFourVelocity (1.0);

    try
    {
        const double u = getTableForTheta(t).lookupFunctionValue(r);
        auto state = RelativisticWind::WindState (wind, r, u);
        return configureWindState (state, position);
    }
    catch (const std::exception& e)
    {
        std::cout
        << "Warning: sampleWind is out of tabulated bounds at r = "
        << r
        << " theta = "
        << t << ". Falling back to exact integration, which is slow."
        << std::endl;

        return configureWindState (wind.integrate(r), position);
    }
}

double StructuredJetModel::jetStructureEtaOfTheta (double theta) const
{
    const double Q = std::pow (theta / config.jetOpeningAngle, config.jetStructureExponent);
    return config.specificWindPower * std::exp (-Q);
}

double StructuredJetModel::sampleElectronGammaBeta (double kT) const
{
    double urms = std::sqrt (kT < 1 ? 3 * kT : 12 * kT * kT); // approximate RMS four-velocity

    if (false)
    {
        auto pdf = Distributions::makeMaxwellian (kT, Distributions::Pdf);
        auto electronGammaBeta = RandomVariable (new SamplingScheme (pdf, urms * 0.01, urms * 20));
        return electronGammaBeta.sample();
    }
    return urms;
}

const TabulatedFunction& StructuredJetModel::getTableForTheta (double theta) const
{
    if (theta >= tableOfThetas.back())
    {
        throw std::runtime_error ("theta value " + std::to_string (theta) + " is outside tabulated solution");
    }
    const int thetaBin = tableOfThetas.size() * theta / tableOfThetas.back();
    assert (0 <= thetaBin && thetaBin < tableOfSolutions.size());

    return tableOfSolutions[thetaBin];
}

TabulatedFunction StructuredJetModel::tabulateWindSolution (double rmax, double theta) const
{
    const int N = config.tableResolutionRadius;
    const double eta = jetStructureEtaOfTheta (theta);
    auto wind = RelativisticWind();
    wind.setSpecificWindPower (eta);
    wind.setInitialMachNumber (1.5);
    // wind.setInitialFourVelocity (2.0);

    auto table = TabulatedFunction (1.0, rmax, N, TabulatedFunction::useEqualBinWidthsLogarithmic);
    auto solution = wind.integrate (table.getDataX());

    for (int i = 0; i < table.size(); ++i)
    {
        table[i] = solution[i].u;
    }
    return table;
}

void StructuredJetModel::tabulateWindAllAngles (double rmax)
{
    const int N = config.tableResolutionTheta;
    const double thetaMax = config.jetOpeningAngle * 4;

    for (int n = 0; n < N; ++n)
    {
        const double theta = double(n) / (N - 1) * thetaMax;

        tableOfSolutions.push_back (tabulateWindSolution (rmax, theta));
        tableOfThetas.push_back (theta);
    }
}

RelativisticWind::WindState StructuredJetModel::configureWindState (RelativisticWind::WindState state, FourVector position) const
{
    const double eta = jetStructureEtaOfTheta (position.theta());
    state.setLuminosityPerSteradian (config.luminosityPerSteradian * eta);
    state.setInnerRadiusCm (config.innerRadiusCm);
    state.setLeptonsPerBaryon (config.leptonsPerBaryon);
    state.setPhotonsPerBaryon (config.photonsPerBaryon);
    state.setPropagationAngle (position.getUnitThreeVector());
    return state;
}
