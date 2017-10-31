#include <iostream>
#include <cassert>
#include "StructuredJetModel.hpp"
#include "LorentzBoost.hpp"
#include "RelativisticWind.hpp"
#include "TabulatedFunction.hpp"
#include "Distributions.hpp"
#include "ScatteringOperations.hpp"





// ========================================================================
StructuredJetModel::StructuredJetModel (Config config) : config (config)
{
    tabulateWindAllAngles (config.outermostRadius);
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
    auto photonEnergy = RandomVariable::diracDelta (kT);

    photon.position = FourVector::spaceLikeInDirection (r, rhat);
    photon.momentum = FourVector::nullWithUnitVector (UnitVector::sampleIsotropic());
    photon.momentum *= 3 * photonEnergy.sample();
    photon.momentum.transformBy(u);

    return photon;
}

StructuredJetModel::Electron StructuredJetModel::generateElectron (const Photon& photon,
    RelativisticWind::WindState state) const
{
    auto sops = ScatteringOperations();
    double uth = std::sqrt (3. * state.temperature()); // TODO: properly sample a Maxwellian here
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
    wind.setInitialFourVelocity (2.0);

    try {
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
    wind.setInitialFourVelocity (2.0);

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
