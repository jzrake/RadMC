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
double StructuredJetModel::Photon::lagTime (double lengthUnits) const
{
    PhysicsConstants physics;
    auto x = position;
    auto k = momentum.getUnitThreeVector();
    return (x[0] - x.projectedAlong(k).radius()) * lengthUnits / physics.c;
}




// ========================================================================
StructuredJetModel::StructuredJetModel (Config config) : config (config)
{
    tabulateWindAllAngles (config.outermostRadius);
}

double StructuredJetModel::sampleTheta (double fraction) const
{
    auto mu = RandomVariable::uniformOver (std::cos (fraction * config.jetPolarBoundary), 1.0);
    return std::acos (mu.sample());
}

double StructuredJetModel::approximatePhotosphere (double theta) const
{
    const double eta = jetStructureEtaOfTheta (theta);
    const double eff = jetStructureEffOfTheta (theta);
    const double Ndot = eff / physics.mp / physics.c / physics.c;
    return 0.5 * Ndot * physics.st / physics.c * std::pow (eta, -2);
}

double StructuredJetModel::approximateLagTime (double theta) const
{
    const double G = jetStructureEtaOfTheta (theta);
    const double r = approximatePhotosphere (theta);
    return 0.5 * r / G / G / physics.c;
}

double StructuredJetModel::angularLuminosity (double theta) const
{
    const double eta = jetStructureEtaOfTheta (theta);
    const double eff = jetStructureEffOfTheta (theta);
    return eta * eff;
}

double StructuredJetModel::totalLuminosity() const
{
    auto L = [this] (double theta)
    {
        const double eta = jetStructureEtaOfTheta (theta);
        const double eff = jetStructureEffOfTheta (theta);
        return 2 * M_PI * eta * eff * std::sin (theta);
    };
    const GaussianQuadrature gauss(8);
    return gauss.computeDefiniteIntegral (L, 0.0, M_PI / 2, 1e-8);
}

double StructuredJetModel::fluidPropagationTimeToRadius (double r, double theta) const
{
    // TODO: Improve estimate here (account for acceleration)

    if (r < jetStructureEtaOfTheta (theta))
    {
        throw std::runtime_error ("The radius is less than the saturation radius; "
            "estimate of propagation time is wrong");
    }
    const auto& table = getTableForTheta (theta);
    const double u = table.lookupFunctionValue(r);
    const double v = u / std::sqrt (1 + u * u);
    const double t = r / v * config.innerRadiusCm / physics.c;
    return t;
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
    const double phi = RandomVariable::sampleUniformAzimuth();
    const auto state = sampleWindSpherical (r, theta);
    const auto rhat = UnitVector (std::cos (theta), phi); // wind propagation angle
    const auto u = FourVector::fromGammaBetaAndUnitVector (state.u, rhat);

    const double kT = state.temperature();
    const auto photonEnergy = RandomVariable::diracDelta (3 * kT);

    auto photon = Photon();
    photon.position = FourVector::spaceLikeInDirection (r, rhat);
    photon.position[0] = fluidPropagationTimeToRadius (r, theta) * physics.c / config.innerRadiusCm;
    photon.momentum = FourVector::nullWithUnitVector (UnitVector::sampleIsotropic());
    photon.momentum *= photonEnergy.sample();
    photon.momentum.transformBy(-u);

    return photon;
}

StructuredJetModel::Electron StructuredJetModel::generateElectron (const Photon& photon,
    RelativisticWind::WindState state) const
{
    const double uth = sampleElectronGammaBeta(state.temperature());
    const auto sops = ScatteringOperations();
    return sops.sampleScatteredParticlesInFrame (state.fourVelocity(), photon.momentum, uth);
}

RelativisticWind::WindState StructuredJetModel::sampleWindSpherical (double r, double theta) const
{
    const auto x = FourVector::spaceLikeInDirection (r, UnitVector (std::cos (theta), 0.0));
    return sampleWind(x);
}

RelativisticWind::WindState StructuredJetModel::sampleWind (const FourVector& position) const
{
    const double r = position.radius();
    const double t = position.theta();
    const double u = getTableForTheta(t).lookupFunctionValue(r);
    auto wind = makeWindSolver(t);
    auto state = RelativisticWind::WindState (wind, r, u);
    return configureWindState (state, position);
}

double StructuredJetModel::jetStructureEtaOfTheta (double theta) const
{
    if (config.jetStructureExponent == 0) return config.specificWindPower;
    const double Q = std::pow (theta / config.jetOpeningAngle, config.jetStructureExponent);
    return std::exp (-Q) * config.specificWindPower;
}

double StructuredJetModel::jetStructureEffOfTheta (double theta) const
{
    if (config.jetStructureExponent == 0) return config.luminosityPerSteradian / config.specificWindPower;
    const double Q = std::pow (theta / config.jetOpeningAngle, config.jetStructureExponent);
    return std::exp (-Q) * config.luminosityPerSteradian / config.specificWindPower;
    // return config.luminosityPerSteradian / config.specificWindPower;
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
        throw std::runtime_error ("theta = " + std::to_string (theta) + " is outside tabulated solution");
    }
    const int thetaBin = tableOfThetas.size() * theta / tableOfThetas.back();
    assert (0 <= thetaBin && thetaBin < tableOfSolutions.size());

    return tableOfSolutions[thetaBin];
}

TabulatedFunction StructuredJetModel::tabulateWindSolution (double rmax, double theta) const
{
    const int N = config.tableResolutionRadius;
    auto table = TabulatedFunction (1.0, rmax, N, TabulatedFunction::useEqualBinWidthsLogarithmic);
    auto wind = makeWindSolver (theta);
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

    for (int n = 0; n < N; ++n)
    {
        const double theta = double(n) / (N - 1) * config.jetPolarBoundary;
        tableOfSolutions.push_back (tabulateWindSolution (rmax, theta));
        tableOfThetas.push_back (theta);
    }
}

RelativisticWind::WindState StructuredJetModel::configureWindState (RelativisticWind::WindState state, FourVector position) const
{
    const double eta = jetStructureEtaOfTheta (position.theta());
    const double eff = jetStructureEffOfTheta (position.theta());
    state.setLuminosityPerSteradian (eta * eff);
    state.setInnerRadiusCm (config.innerRadiusCm);
    state.setLeptonsPerBaryon (config.leptonsPerBaryon);
    state.setPhotonsPerBaryon (config.photonsPerBaryon);
    state.setPropagationAngle (position.getUnitThreeVector());
    return state;
}

RelativisticWind StructuredJetModel::makeWindSolver (double theta) const
{
    const double eta = jetStructureEtaOfTheta (theta);
    auto wind = RelativisticWind();
    wind.setSpecificWindPower (eta);
    wind.setInitialFourVelocity (1.0);
    return wind;
}
