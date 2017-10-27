#include "StructuredJetModel.hpp"
#include "LorentzBoost.hpp"
#include "RelativisticWind.hpp"




StructuredJetModel::StructuredJetModel()
{

}

std::vector<StructuredJetModel::Photon> StructuredJetModel::generatePhotonPath (double innerRadius, double theta)
{
    auto photon = generatePhoton (innerRadius, theta);

    auto path = std::vector<StructuredJetModel::Photon>();

    for (int n = 0; n < 100; ++n)
    {
        Electron electron;
        doComptonScattering (photon, electron);
        path.push_back (photon);
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

    // TODO: initialize photon energy

    return photon;
}

StructuredJetModel::Electron StructuredJetModel::generateElectron (double r, double theta) const
{
    auto electron = Electron();
    const double phi = RandomVariable::sampleUniformAzimuth();
    const double ur = sampleWind (r, theta).u;
    const auto radiusUnitVector = UnitVector (std::cos (theta), phi);
    const auto u = FourVector::fromGammaBetaAndUnitVector (ur, radiusUnitVector);

    electron.momentum = FourVector::nullWithUnitVector (UnitVector::sampleIsotropic()).transformedBy(u);

    return electron;
}

RelativisticWind::WindState StructuredJetModel::sampleWind (double r, double theta) const
{
    auto wind = RelativisticWind().setSpecificWindPower (10.0).setInitialFourVelocity (2.0);
    return wind.integrate(r);
}

FourVector StructuredJetModel::doComptonScattering (Photon& photon, Electron& electron) const
{
    // Photon four-momentum in the electron rest frame
    LorentzBoost L (electron.getFourVelocity());
    FourVector p0 = photon.momentum.transformedBy(L);

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

    return dp;
}
