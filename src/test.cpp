#define COMPILE_TESTS
#ifdef COMPILE_TESTS

#define CATCH_CONFIG_FAST_COMPILE
#include "../lib/catch.hpp"
#include "FourVector.hpp"
#include "LorentzBoost.hpp"
#include "QuadratureRule.hpp"
#include "NewtonRaphesonSolver.hpp"
#include "RootBracketingSolver.hpp"
#include "RotationMatrix.hpp"
#include "RungeKutta.hpp"


SCENARIO ("FourVector")
{
    GIVEN ("u = [1 0 0 0]")
    {
        THEN ("Spacelike, timelike, and null methods work")
        {
            FourVector u (1, 0, 0, 0);
            REQUIRE (u.isNull() == 0);
            REQUIRE (u.isSpacelike() == 1);
            REQUIRE (u.isTimelike() == 0);
        }
    }
    GIVEN ("u = [0 1 1 1]")
    {
        THEN ("Spacelike, timelike, and null methods work")
        {
            FourVector u (0, 1, 1, 1);
            REQUIRE (u.isNull() == 0);
            REQUIRE (u.isSpacelike() == 0);
            REQUIRE (u.isTimelike() == 1);
        }
    }
    GIVEN ("u = [-1 1 0 0]")
    {
        THEN ("Spacelike, timelike, and null methods work")
        {
            FourVector u (-1, 1, 0, 0);
            REQUIRE (u.isNull() == 1);
            REQUIRE (u.isSpacelike() == 0.0);
            REQUIRE (u.isTimelike() == 0);
        }
    }
}

SCENARIO ("LorentzBoost")
{
    GIVEN ("Lorentz boost with v=0.5 along x")
    {
        THEN ("Four vectors are transformed correctly")
        {
            FourVector u = FourVector::fromThreeVelocity (0.5, 0.0, 0.0);
            LorentzBoost L (u);
            REQUIRE (u.getTimeComponent() == 1.0 / std::sqrt (1.0 - 0.5 * 0.5));
            REQUIRE ((L * u).getTimeComponent() == Approx (1.0));
            REQUIRE ((L *(-u)).isFourVelocity() == true);
            REQUIRE ((L *(-u)).getThreeVelocityMagnitude() == Approx ((0.5 + 0.5) / (1 + 0.5 * 0.5)));
        }
    }
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

SCENARIO ("Root finding 1", "[NewtonRaphesonSolver]")
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

SCENARIO ("Root finding 2", "[RootBracketingSolver]")
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

SCENARIO ("Rotation matrices", "[RotationMatrix]")
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

SCENARIO ("Sampling unit vectors works OK")
{
    GIVEN ("Pitch angles are distributed around -1")
    {
        RandomVariable delta = RandomVariable::diracDelta (-1);

        WHEN ("The target vector is xhat")
        {
            UnitVector n = UnitVector::xhat;
            THEN ("Sampled unit vectors are -xhat")
            {
                REQUIRE (n.sampleAxisymmetric (delta).getX() == -1);
            }
        }
        WHEN ("The target vector is yhat")
        {
            UnitVector n = UnitVector::yhat;
            THEN ("Sampled unit vectors are -yhat")
            {
                REQUIRE (n.sampleAxisymmetric (delta).getY() == -1);
            }
        }
        WHEN ("The target vector is zhat")
        {
            UnitVector n = UnitVector::zhat;
            THEN ("Sampled unit vectors are -zhat")
            {
                REQUIRE (n.sampleAxisymmetric (delta).getZ() == -1);
            }
        }
    }
}

SCENARIO ("RungeKutta integration works", "[RungeKutta]")
{
    RungeKutta solver;

    GIVEN ("A the function y'=1")
    {
        solver.setFunction ([] (double y, double t) { return 1.0;} );
        solver.setValues (0.0, 0.0);
        solver.integrate (1.0);

        THEN ("Integrating to t=1 gives y=1")
        {
            REQUIRE (solver.getT() == 1.0);
            REQUIRE (solver.getY() == 1.0);
            REQUIRE (solver.getDivisionsForLast() == 1);
        }
    }
    GIVEN ("A the function y'=y")
    {
        solver.setFunction ([] (double y, double t) { return y;} );
        solver.setValues (1.0, 0.0);
        solver.integrate (10.0);
        solver.setTolerance (1e-13);

        THEN ("Integrating to t=1 gives y=e")
        {
            REQUIRE (solver.getT() == 10.0);
            REQUIRE (solver.getY() == Approx (std::exp (10.0)).epsilon (1e-13));
            REQUIRE (solver.getDivisionsForLast() == 8192);
        }
    }
}
#endif // COMPILE_TESTS

