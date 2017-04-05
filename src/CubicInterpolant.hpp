#ifndef CubicInterpolant_hpp
#define CubicInterpolant_hpp



class CubicInterpolant
{
public:
    /**
    Returns a cubic interpolating function f(x) which passes through the
    coordinates (xa, fa) and (xb, fb), and has derivatitves f'(xa) = ga and
    f'(xb) = gb.
    */
    static CubicInterpolant fromValuesAndDerivatives (double xa, double xb, double fa, double fb, double ga, double gb);

    CubicInterpolant();

    CubicInterpolant (double c0, double c1, double c2, double c3);

    double evaluate (double x);

private:
    double c0;
    double c1;
    double c2;
    double c3;
};



#endif
