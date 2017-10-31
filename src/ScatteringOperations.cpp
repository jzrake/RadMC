#include "ScatteringOperations.hpp"
#include "Distributions.hpp"
#include "LorentzBoost.hpp"




ScatteringOperations::ScatteringOperations()
{

}

FourVector ScatteringOperations::sampleScatteredParticles (FourVector k, double u) const
{
    double v = FourVector::betaFromGammaBeta(u);

    // Sample the electron propagation direction. The distribution of electron
    // velocities, given that scattering has occurred, is axisymmetric around
    // the photon propagation vector k, but heavier at -k than +k due to the
    // relative motion.
    RandomVariable relativeMu (Distributions::makePitchAngle (v, Distributions::Qnt));
    UnitVector nhat = k.getUnitThreeVector().sampleAxisymmetric (relativeMu);

    return FourVector::fromBetaAndUnitVector (v, nhat);
}

FourVector ScatteringOperations::sampleScatteredParticlesInFrame (FourVector restFrame, FourVector k, double u) const
{
    LorentzBoost L (restFrame);

    // Get the photon in the fluid rest frame
    FourVector q = k.transformedBy(L);
    FourVector f = sampleScatteredParticles (q, u);
    FourVector e = f.transformedBy (L.inverted());
    return e;
}

FourVector ScatteringOperations::comptonScatter (FourVector& photon, FourVector& electron) const
{
    // Photon four-momentum in the electron rest frame
    LorentzBoost L (electron);
    FourVector p0 = photon.transformedBy(L);

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
    photon   += dp;
    electron -= dp;

    return dp;
}

std::vector<FourVector> ScatteringOperations::comptonize (double temperature, int scatterings) const
{
    auto path = std::vector<FourVector>();
    auto photon = FourVector::nullWithUnitVector (UnitVector::sampleIsotropic());

    double kT = temperature;
    double urms = std::sqrt (kT < 1 ? 3 * kT : 12 * kT * kT); // approximate RMS four-velocity
    auto pdf = Distributions::makeMaxwellian (kT, Distributions::Pdf);
    auto electronGammaBeta = RandomVariable (RandomVariable::fromPdf (pdf, urms * 0.01, urms * 10));        

    for (int n = 0; n < scatterings; ++n)
    {
        const double uth = electronGammaBeta.sample();
        FourVector electron = sampleScatteredParticles (photon, uth);
        comptonScatter (photon, electron);
        path.push_back (photon);
    }
    return path;
}

