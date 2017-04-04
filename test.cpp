#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "FourVector.hpp"
#include "LorentzBoost.hpp"
#include "QuadratureRule.hpp"
#include "NewtonRaphesonSolver.hpp"
#include "RootBracketingSolver.hpp"
#include "RotationMatrix.hpp"


TEST_CASE ("FourVector")
{
    SECTION ("u = [1 0 0 0]")
    {
        FourVector u (1, 0, 0, 0);
        REQUIRE (u.isNull() == 0);
        REQUIRE (u.isSpacelike() == 1);
        REQUIRE (u.isTimelike() == 0);
    }
    SECTION ("u = [0 1 1 1]")
    {
        FourVector u (0, 1, 1, 1);
        REQUIRE (u.isNull() == 0);
        REQUIRE (u.isSpacelike() == 0);
        REQUIRE (u.isTimelike() == 1);
    }
    SECTION ("u = [-1 1 0 0]")
    {
        FourVector u (-1, 1, 0, 0);
        REQUIRE (u.isNull() == 1);
        REQUIRE (u.isSpacelike() == 0.0);
        REQUIRE (u.isTimelike() == 0);
    }
}

TEST_CASE ("LorentzBoost")
{
    FourVector u = FourVector::fromThreeVelocity (0.5, 0.0, 0.0);
    LorentzBoost L (u);
    REQUIRE (u.getTimeComponent() == 1.0 / std::sqrt (1.0 - 0.5 * 0.5));
    REQUIRE ((L * u).getTimeComponent() == Approx (1.0));
    REQUIRE ((L *(-u)).isFourVelocity() == true);
    REQUIRE ((L *(-u)).getThreeVelocityMagnitude() == Approx ((0.5 + 0.5) / (1 + 0.5 * 0.5)));
}

SCENARIO ("Numerical integration", "[QuadratureRule]" )
{
    GIVEN ("The function y = sin(x)^2")
    {
        auto f = [] (double x) { return sin (x) * sin (x) / (2 * M_PI); };
        GaussianQuadrature gaussianQuadrature;

        WHEN ("q = 5")
        {
            gaussianQuadrature.setPolynomialDegree (5);
            THEN ("<y> = 1/2")
            {
                REQUIRE (gaussianQuadrature.integrate (f, 0, 2 * M_PI) == Approx (0.5).epsilon (1e-2));
            }
        }
        WHEN ("q = 9")
        {
            gaussianQuadrature.setPolynomialDegree (9);
            THEN ("<y> = 1/2")
            {
                REQUIRE (gaussianQuadrature.integrate (f, 0, 2 * M_PI) == Approx (0.5).epsilon (1e-7));
            }
        }
        WHEN ("q = 12")
        {
            gaussianQuadrature.setPolynomialDegree (12);
            THEN ("<y> = 1/2")
            {
                REQUIRE (gaussianQuadrature.integrate (f, 0, 2 * M_PI) == Approx (0.5).epsilon (1e-12));
            }
        }
    }


    GIVEN ("The function y = 1 / sqrt (M_PI) * exp (-x * x)")
    {
        double sigma = 1 / std::sqrt (2);
        auto f = [] (double x) { return 1 / std::sqrt (M_PI) * std::exp (-x * x); };

        WHEN ("q = 24 with a single bin")
        {
            GaussianQuadrature gaussianQuadrature;
            gaussianQuadrature.setPolynomialDegree (24);
            THEN ("The integral is exact out to 6 sigma")
            {
                for (int n = 0; n < 4; ++n)
                {
                    double x = 6 * sigma * double (n) / 4;
                    double F = gaussianQuadrature.integrate (f, -x, x);
                    REQUIRE (F == Approx (std::erf (x)).epsilon (1e-12));
                }
            }
        }

        WHEN ("q = 8 with adaptive bins")
        {
            THEN ("The integral is exact out to 6 sigma, and converges at 32 bins")
            {
                GaussianQuadrature gaussianQuadrature;
                gaussianQuadrature.setPolynomialDegree (8);
                QuadratureRule::EvaluationDetails details;
                double F = gaussianQuadrature.computeDefiniteIntegral (f, -6, 6, 1e-12, details);
                REQUIRE (F == Approx (std::erf (6)).epsilon (1e-12));
                REQUIRE (details.numberOfBinsUsed == 32);
            }
        }

        WHEN ("We use Simpson's rule with adaptive bins")
        {
            THEN ("The integral is exact out to 6 sigma, and converges at 64 bins")
            {
                SimpsonRule simpson;
                QuadratureRule::EvaluationDetails details;
                double F = simpson.computeDefiniteIntegral (f, -6, 6, 1e-12, details);
                REQUIRE (F == Approx (std::erf (6)).epsilon (1e-12));
                REQUIRE (details.numberOfBinsUsed == 64);
            }
        }
    }
}

SCENARIO ("RootFinding1", "[NewtonRaphesonSolver]")
{
    GIVEN ("The function F = (x + 2) * (x + 4) * (x - 3)")
    {
        auto F = [] (double x) { return (x + 2) * (x + 4) * (x - 3); };
        auto f = [] (double x) { return -10 + 3 * x * (2 + x); };

        NewtonRaphesonSolver solver (F, f);
        NewtonRaphesonSolver::Result result;
        double tolerance = 1e-14;

        THEN ("The solver finds the root to machine accuracy in 5 iterations")
        {
            REQUIRE (solver.solve (0, 3.5, tolerance, result) == Approx (3.0));
            REQUIRE (result.numIterations <= 5);
            REQUIRE (solver.solve (0, 2.5, tolerance, result) == Approx (3.0));
            REQUIRE (result.numIterations <= 5);
        }
    }
}

SCENARIO ("RootFinding2", "[RootBracketingSolver]")
{
    GIVEN ("The function F = (x + 2) * (x + 4) * (x - 3) + 2")
    {
        auto F = [] (double x) { return (x + 2) * (x + 4) * (x - 3) + 2; };

        RootBracketingSolver solver (F);
        RootBracketingSolver::Result result;
        double tolerance = 1e-14;

        THEN ("The solver finds the root to machine accuracy in 48 iterations")
        {
            REQUIRE (solver.solve (2, 2.71, 3.58, tolerance, result) == Approx (3.0));
            REQUIRE (result.numIterations <= 48);
        }
    }
}

SCENARIO ("RotationMatrix", "[]")
{
    GIVEN ("The unit vector zhat")
    {
        UnitVector xhat = UnitVector::normalizeFrom (1, 0, 0);
        UnitVector yhat = UnitVector::normalizeFrom (0, 1, 0);
        UnitVector zhat = UnitVector::normalizeFrom (0, 0, 1);

        THEN ("Y (pi / 2) zhat = xhat")
        {
            REQUIRE ((RotationMatrix::aboutY (M_PI / 2) * zhat).pitchAngleMu == Approx (xhat.pitchAngleMu));
            REQUIRE ((RotationMatrix::aboutY (M_PI / 2) * zhat).azimuthalAnglePhi == Approx (xhat.azimuthalAnglePhi));
        }

        THEN ("X (-pi / 2) zhat = yhat")
        {
            REQUIRE ((RotationMatrix::aboutX (-M_PI / 2) * zhat).pitchAngleMu == Approx (yhat.pitchAngleMu));
            REQUIRE ((RotationMatrix::aboutX (-M_PI / 2) * zhat).azimuthalAnglePhi == Approx (yhat.azimuthalAnglePhi));
        }

        THEN ("Z (pi / 2) xhat = yhat")
        {
            REQUIRE ((RotationMatrix::aboutZ (M_PI / 2) * xhat).pitchAngleMu == Approx (yhat.pitchAngleMu));
            REQUIRE ((RotationMatrix::aboutZ (M_PI / 2) * xhat).azimuthalAnglePhi == Approx (yhat.azimuthalAnglePhi));
        }

        THEN ("zhat.withPolarAxis (xhat) = xhat")
        {
            REQUIRE (zhat.withPolarAxis (xhat).pitchAngleMu == Approx (xhat.pitchAngleMu));
            REQUIRE (zhat.withPolarAxis (xhat).azimuthalAnglePhi == Approx (xhat.azimuthalAnglePhi));
        }

        THEN ("zhat.withPolarAxis (yhat) = yhat")
        {
            REQUIRE (zhat.withPolarAxis (yhat).pitchAngleMu == Approx (yhat.pitchAngleMu));
            REQUIRE (zhat.withPolarAxis (yhat).azimuthalAnglePhi == Approx (yhat.azimuthalAnglePhi));
        }

        THEN ("zhat.withPolarAxis (zhat) = zhat")
        {
            REQUIRE (zhat.withPolarAxis (zhat).pitchAngleMu == Approx (zhat.pitchAngleMu));
            REQUIRE (zhat.withPolarAxis (zhat).azimuthalAnglePhi == Approx (zhat.azimuthalAnglePhi));
        }
    }
}



