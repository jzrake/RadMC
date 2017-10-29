#include "StructuredJetModel.hpp"
#include "LorentzBoost.hpp"
#include "RelativisticWind.hpp"
#include "TabulatedFunction.hpp"
#include "Distributions.hpp"
#include "ScatteringOperations.hpp"





// ========================================================================
StructuredJetModel::StructuredJetModel()
{

}

std::vector<StructuredJetModel::Photon> StructuredJetModel::generatePhotonPath (double innerRadius, double theta)
{
    auto scatteringOps = ScatteringOperations();
    auto photon = generatePhoton (innerRadius, theta);
    auto path = std::vector<StructuredJetModel::Photon>();

    for (int n = 0; n < 100; ++n)
    {
        auto state = sampleWind (photon.position);

        // TODO: sample scattering length from exponential
        const double dr = state.thomsonMeanFreePath (photon.momentum.getUnitThreeVector());
        const double dt = dr / innerRadiusCm;

        photon.position += photon.momentum / photon.momentum[0] * dt;

        // TODO: improve efficiency by passing wind stat to generateElectron
        Electron electron = generateElectron (photon);
        scatteringOps.comptonScatter (photon.momentum, electron.momentum);
        path.push_back (photon);
    }
    return path;
}

StructuredJetModel::Photon StructuredJetModel::generatePhoton (double r, double theta) const
{
    auto photon = Photon();
    const double phi = RandomVariable::sampleUniformAzimuth();
    const auto state = sampleWind (r, theta);
    const auto rhat = UnitVector (std::cos (theta), phi); // wind propagation angle
    const auto u = FourVector::fromGammaBetaAndUnitVector (state.u, rhat);

    const double kT = sampleWind (r, theta).temperature();
    auto photonEnergy = RandomVariable::diracDelta (kT);

    photon.position = FourVector::spaceLikeInDirection (r, rhat);
    photon.momentum = FourVector::nullWithUnitVector (UnitVector::sampleIsotropic());
    photon.momentum *= photonEnergy.sample();
    photon.momentum.transformBy(u);

    return photon;
}

StructuredJetModel::Electron StructuredJetModel::generateElectron (const Photon& photon) const
{
    auto sops = ScatteringOperations();
    auto state = sampleWind (photon.position);
    double uth = std::sqrt (3. * state.temperature()); // TODO: properly sample a Maxwellian here
    return sops.sampleScatteredParticles (state.fourVelocity(), photon.momentum, uth);
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
    state.setLuminosityPerSteradian (luminosityPerSteradian);
    state.setInnerRadiusCm (innerRadiusCm);
    state.setLeptonsPerBaryon (leptonsPerBaryon);
    state.setPhotonsPerBaryon (photonsPerBaryon);
    state.setPropagationAngle (position.getUnitThreeVector());

    return state;
}
