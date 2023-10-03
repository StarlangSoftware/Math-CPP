//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#include <cmath>
#include "Distribution.h"

/**
 * The ex method takes a double x as an input, if x is less than -BIGX it returns 0, otherwise it returns Euler's number
 * <i>e</i> raised to the power of x.
 *
 * @param x double input.
 * @return 0 if input is less than -BIGX, Euler's number <i>e</i> raised to the power of x otherwise.
 */
double Distribution::ex(double x) {
    if (x < -BIGX) {
        return 0;
    }
    return exp(x);
}

/**
 * The beta method takes a double array x as an input. It loops through x and accumulates
 * the value of gammaLn(x), also it sums up the items of x and returns (accumulated result - gammaLn of this summation).
 *
 * @param x double array input.
 * @param size size of x.
 * @return gammaLn(sum).
 */
double Distribution::beta(double *x, int size) {
    double sum = 0.0, result = 0.0;
    int i;
    for (i = 0; i < size; i++) {
        result += gammaLn(x[i]);
        sum += x[i];
    }
    result -= gammaLn(sum);
    return result;
}

/**
 * The gammaLn method takes a double x as an input and returns the logarithmic result of the gamma distribution at point x.
 *
 * @param x double input.
 * @return the logarithmic result of the gamma distribution at point x.
 */
double Distribution::gammaLn(double x) {
    double y, tmp, ser;
    double cof[] = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
    int j;
    y = x;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++) {
        ser += cof[j] / ++y;
    }
    return -tmp + log(2.5066282746310005 * ser / x);
}

/**
 * The zNormal method performs the Z-Normalization. It ensures, that all elements of the input vector are transformed
 * into the output vector whose mean is approximately 0 while the standard deviation is in a range close to 1.
 *
 * @param z double input.
 * @return normalized value of given input.
 */
double Distribution::zNormal(double z) {
    double y, x, w;
    if (z == 0.0) {
        x = 0.0;
    } else {
        y = 0.5 * abs(z);
        if (y >= (Z_MAX * 0.5)) {
            x = 1.0;
        } else {
            if (y < 1.0) {
                w = y * y;
                x = ((((((((0.000124818987 * w - 0.001075204047) * w + 0.005198775019) * w - 0.019198292004) * w + 0.059054035642) * w - 0.151968751364) * w + 0.319152932694) * w - 0.531923007300) * w + 0.797884560593) * y * 2.0;
            } else {
                y -= 2.0;
                x = (((((((((((((-0.000045255659 * y + 0.000152529290) * y - 0.000019538132) * y - 0.000676904986) * y + 0.001390604284) * y - 0.000794620820) * y - 0.002034254874) * y + 0.006549791214) * y - 0.010557625006) * y + 0.011630447319) * y - 0.009279453341) * y + 0.005353579108) * y - 0.002141268741) * y + 0.000535310849) * y + 0.999936657524;
            }
        }
    }
    if (z > 0.0) {
        return ((x + 1.0) * 0.5);
    } else {
        return ((1.0 - x) * 0.5);
    }
}

/**
 * The zInverse method returns the Z-Inverse of given probability.
 *
 * @param p double probability.
 * @return the Z-Inverse of given probability.
 */
double Distribution::zInverse(double p) {
    double minz = -Z_MAX;
    double maxz = Z_MAX;
    double zval = 0.0;
    double pval;
    if (p <= 0.0 || p >= 1.0) {
        return (0.0);
    }
    while (maxz - minz > Z_EPSILON) {
        pval = zNormal(zval);
        if (pval > p) {
            maxz = zval;
        } else {
            minz = zval;
        }
        zval = (maxz + minz) * 0.5;
    }
    return zval;
}

/**
 * The chiSquare method is used to determine whether there is a significant difference between the expected
 * frequencies and the observed frequencies in one or more categories. It takes a double input x and an integer freedom
 * for degrees of freedom as inputs. It returns the Chi Squared result.
 *
 * @param x       double input.
 * @param freedom integer input for degrees of freedom.
 * @return the Chi Squared result.
 */
double Distribution::chiSquare(double x, int freedom) {
    double a, y = 0, s;
    double e, c, z;
    bool even;
    if (x <= 0.0 || freedom < 1) {
        return (1.0);
    }
    a = 0.5 * x;
    even = (freedom % 2 == 0);
    if (freedom > 1) {
        y = ex(-a);
    }
    if (even) {
        s = y;
    } else {
        s = (2.0 * zNormal(-sqrt(x)));
    }
    if (freedom > 2) {
        x = 0.5 * (freedom - 1.0);
        if (even) {
            z = 1.0;
        } else {
            z = 0.5;
        }
        if (a > BIGX) {
            if (even) {
                e = 0.0;
            } else {
                e = LOG_SQRT_PI;
            }
            c = log(a);
            while (z <= x) {
                e = log(z) + e;
                s += ex(c * z - a - e);
                z += 1.0;
            }
            return (s);
        } else {
            if (even) {
                e = 1.0;
            } else {
                e = (I_SQRT_PI / sqrt(a));
            }
            c = 0.0;
            while (z <= x) {
                e = e * (a / z);
                c = c + e;
                z += 1.0;
            }
            return (c * y + s);
        }
    } else {
        return s;
    }
}

/**
 * The chiSquareInverse method returns the Chi Square-Inverse of given probability with given degree of freedom.
 *
 * @param p       double probability.
 * @param freedom integer input for degrees of freedom.
 * @return chi square-Inverse of given probability.
 */
double Distribution::chiSquareInverse(double p, int freedom) {
    double minchisq = 0.0;
    double maxchisq = CHI_MAX;
    double chisqval;
    if (p <= 0.0) {
        return maxchisq;
    } else {
        if (p >= 1.0) {
            return 0.0;
        }
    }
    chisqval = freedom / sqrt(p);
    while (maxchisq - minchisq > CHI_EPSILON) {
        if (chiSquare(chisqval, freedom) < p) {
            maxchisq = chisqval;
        } else {
            minchisq = chisqval;
        }
        chisqval = (maxchisq + minchisq) * 0.5;
    }
    return chisqval;
}

/**
 * The fDistribution method is used to observe whether two samples have the same variance. It takes a double input F
 * and two integer freedom1 and freedom2 for degrees of freedom as inputs. It returns the F-Distribution result.
 *
 * @param F        double input.
 * @param freedom1 integer input for degrees of freedom.
 * @param freedom2 integer input for degrees of freedom.
 * @return the F-Distribution result.
 */
double Distribution::fDistribution(double F, int freedom1, int freedom2) {
    int i, j;
    int a, b;
    double w, y, z, d, p;
    if (F < F_EPSILON || freedom1 < 1 || freedom2 < 1) {
        return (1.0);
    }
    if (freedom1 % 2 != 0) {
        a = 1;
    } else {
        a = 2;
    }
    if (freedom2 % 2 != 0) {
        b = 1;
    } else {
        b = 2;
    }
    w = (F * freedom1) / freedom2;
    z = 1.0 / (1.0 + w);
    if (a == 1) {
        if (b == 1) {
            p = sqrt(w);
            y = I_PI;
            d = y * z / p;
            p = 2.0 * y * atan(p);
        } else {
            p = sqrt(w * z);
            d = 0.5 * p * z / w;
        }
    } else {
        if (b == 1) {
            p = sqrt(z);
            d = 0.5 * z * p;
            p = 1.0 - p;
        } else {
            d = z * z;
            p = w * z;
        }
    }
    y = 2.0 * w / z;
    for (j = b + 2; j <= freedom2; j += 2) {
        d *= (1.0 + a / (j - 2.0)) * z;
        if (a == 1) {
            p = p + d * y / (j - 1.0);
        } else {
            p = (p + w) * z;
        }
    }
    y = w * z;
    z = 2.0 / z;
    b = freedom2 - 2;
    for (i = a + 2; i <= freedom1; i += 2) {
        j = i + b;
        d *= y * j / (i - 2.0);
        p -= z * d / j;
    }
    if (p < 0.0) {
        p = 0.0;
    } else {
        if (p > 1.0) {
            p = 1.0;
        }
    }
    return 1.0 - p;
}

/**
 * The fDistributionInverse method returns the F-Distribution Inverse of given probability.
 *
 * @param p        double probability.
 * @param freedom1 integer input for degrees of freedom.
 * @param freedom2 integer input for degrees of freedom.
 * @return the F-Distribution Inverse of given probability.
 */
double Distribution::fDistributionInverse(double p, int freedom1, int freedom2) {
    double fval;
    double maxf = F_MAX;
    double minf = 0.0;
    if (p <= 0.0 || p >= 1.0) {
        return (0.0);
    }
    if (freedom1 == freedom2 && freedom1 > 2500) {
        return 1 + 4.0 / freedom1;
    }
    fval = 1.0 / p;
    while (abs(maxf - minf) > F_EPSILON) {
        if (fDistribution(fval, freedom1, freedom2) < p) {
            maxf = fval;
        } else {
            minf = fval;
        }
        fval = (maxf + minf) * 0.5;
    }
    return fval;
}

/**
 * The tDistribution method is used instead of the normal distribution when there is small samples. It takes a double input T
 * and an integer freedom for degree of freedom as inputs. It returns the T-Distribution result by using F-Distribution method.
 *
 * @param T       double input.
 * @param freedom integer input for degrees of freedom.
 * @return the T-Distribution result.
 */
double Distribution::tDistribution(double T, int freedom) {
    if (T >= 0) {
        return fDistribution(T * T, 1, freedom) / 2;
    } else {
        return 1 - fDistribution(T * T, 1, freedom) / 2;
    }
}

/**
 * The tDistributionInverse method returns the T-Distribution Inverse of given probability.
 *
 * @param p       double probability.
 * @param freedom integer input for degrees of freedom.
 * @return the T-Distribution Inverse of given probability.
 */
double Distribution::tDistributionInverse(double p, int freedom) {
    if (p < 0.5) {
        return sqrt(fDistributionInverse(p * 2, 1, freedom));
    } else {
        return -sqrt(fDistributionInverse((1 - p) * 2, 1, freedom));
    }
}
