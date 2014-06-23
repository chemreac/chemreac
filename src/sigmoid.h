// C++11 templated source code for algebraic sigmoids of type:
// s(x) := x*((x/lim)**n + 1)**(-1/n)

template <int n, int lim>
double sigmoid(double x){
    return x*pow(pow(x/static_cast<double>(lim), n) + 1, -1.0/n);
}

template <int n, int lim>
double asigmoid(double x){
    return x*pow(-pow(x/static_cast<double>(lim), n) + 1, -1.0/n);
}

template <int n, int lim>
double exp_sigm(double x){
    return exp(sigmoid<n, lim>(x)); // exp of sigmoid
}

template <int n, int lim>
double log_asigm(double x){
    return asigmoid<n, lim>(log(x)) // anti function to exp_sigm
}

template <int n, int lim>
double Dsigmoid(double x){
    return pow(pow(x/static_cast<double>(lim), n) + 1, -1 - 1.0/n);
}

template <int n, int lim>
double D2sigmoid(double x){
    return -pow(x/static_cast<double>(lim), n)*(n + 1)*\
        pow(pow(x/static_cast<double>(lim), n) + 1, -(2*n + 1.0)/n)/x;
}

template <int n, int lim>
double Dexps(double x){
    return exps<n, lim>(x)*Dsigmoid<n, lim>(x);
}
