#include <math.h>

/* (350.0)^-4 */
double sigmoid_algebraic_4(double x){
    const double POW350_MINUS4 = 6.663890045814244e-11;
    return x/sqrt(sqrt(1 + POW350_MINUS4*x*x*x*x)); /* limit=350 */
}
