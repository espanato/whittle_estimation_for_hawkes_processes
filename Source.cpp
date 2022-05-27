
#include <iostream>
#include <armadillo>
#include <complex>
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include <random>
#include <list>
#include <chrono>


using namespace alglib;
using namespace std;


// Cardinal sine: sinc(x) = sin(x) / x

arma::vec _sinc(arma::vec x) {
    arma::vec y(x.n_elem);

    // Iterators
    arma::vec::iterator it_x = x.begin();
    arma::vec::iterator it_x_end = x.end();
    arma::vec::iterator it_y = y.begin();

    // Loop on x
    for (; it_x != it_x_end; ++it_x, ++it_y) {
        if (*it_x == 0.0)
            *it_y = 1.0;
        else
            *it_y = sin(*it_x) / *it_x;
    }
    return y;
};


// Model class

class Model {

protected:
    arma::vec param;
    double binsize;

public:
    Model() {};
    Model(arma::vec param) : param(param), binsize(1.0) {};
    Model(double binsize) : param(arma::zeros<arma::vec>(1)), binsize(binsize) {};
    Model(arma::vec param, double binsize) : param(param), binsize(binsize) {};
    virtual ~Model() {};

    // Methods for long term mean and its derivatives
    double mean();


    arma::vec G(arma::vec xi);    // G(w) = |1-H(w)|^{-2}
    arma::vec f(arma::vec xi);    // f(w) = m * binsize * sinc²(w/2) * G(w/binsize)
    arma::vec f1(arma::vec xi, int trunc);  // f1(w) = sum_{k=-trunc}^{+trunc} f(w + 2*k*pi)
    arma::cx_vec H(arma::vec xi);

    // Whittle likelihood estimation methods
    double whittle(arma::vec& I, int trunc);
};

double Model::mean() {
    return param(0) / (1.0 - param(1));
}

arma::cx_vec Model::H(arma::vec xi) {
    arma::vec factor = param(1) * param(2) / (param(2) * param(2) + xi % xi);
    arma::cx_vec zeta = arma::cx_vec(factor * param(2), -factor % xi);
    return zeta;
}

// G(w) = |1-H(w)|^{-2}

arma::vec Model::G(arma::vec xi) {
    arma::cx_vec temp = arma::cx_double(1.0, 0.0) - H(xi);
    return 1.0 / arma::conv_to<arma::vec>::from(temp % arma::conj(temp));
}

// f(w) = m * binsize * sinc²(w/2) * G(w/binsize)

arma::vec Model::f(arma::vec xi) {
    arma::vec term1 = _sinc(.5 * xi);
    return mean() * binsize * term1 % term1 % G(xi / binsize);
}

// f1(w) = sum_{k=-trunc}^{+trunc} f(w + 2*k*pi)

arma::vec Model::f1(arma::vec xi, int trunc) {
    arma::vec omega = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(-trunc, trunc);
    arma::vec y(xi.n_elem);
    for (arma::uword k = 0; k < xi.n_elem; k++) {
        y(k) = arma::sum(f(xi(k) + omega));
    }
    return y;
}


// Whittle likelihood estimation methods
double Model::whittle(arma::vec& periodogram, int trunc) {
    arma::uword n = periodogram.n_elem + 1;     // + 1 because we removed element 0 in 'whittle.R'
    arma::uword n2 = n / 2;
    arma::vec omega = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(1, n2) / (double)n;
    arma::vec spectrum = f1(omega, trunc);
    return arma::sum(arma::log(spectrum) + periodogram.subvec(0, n2 - 1) / spectrum);
}

// The function to minimize
double nlopt_fn(arma::vec param, arma::vec I, double binsize, int trunc = 5) {
    Model Expo(param, binsize);
    return Expo.whittle(I, trunc);
}

arma::vec periodogram(arma::vec counts, double binsize, int trunc = 5) {
    int n = counts.n_elem;
    arma::cx_vec dft = arma::fft(counts - arma::mean(counts));
    arma::vec I = arma::abs(dft);
    I = I % I / (double)n;
    int M = I.n_elem;
    M = M - 1;
    I = I.rows(1, M); 
    return I;
}

arma::vec init_bfgs(arma::vec counts, double binsize, int trunc = 5) {
    double ymean = arma::mean(counts);
    double wmin = std::numeric_limits<double>::infinity();
    double mu = 0.25;
    double rate = 1;
    arma::vec param = arma::vec({ 0,0,0 });
    for (double mu_ = 0.25; mu_ <= 0.75; mu_ += 0.25) {
        for (double rate_ = 1; rate_ < 6; rate_ += 1) {
            param[0] = ymean * (1 - mu_);
            param[1] = mu_;
            param[2] = rate_;
            Model Expo(param, binsize);
            arma::vec I = periodogram(counts, binsize, 5);
            double whit = Expo.whittle(I, trunc);
            if (whit < wmin) {
                mu = mu_;
                rate = rate_;
                wmin = whit;
            }
        }
    }
    double p1 = ymean * (1 - mu);
    arma::vec x0 = arma::vec({ p1,mu, rate });

    return x0;
}


arma::vec gen_counts(){
    arma::vec counts = arma::vec({});

    ifstream ip("data.csv");

    if (!ip.is_open()) std::cout << "ERROR: File Open";

    string value;

    while (ip.good()) {

        getline(ip, value, ',');
        arma::vec next = arma::vec({ value });
        counts = arma::join_vert(counts, next);


    }
    return counts ;
}

arma::vec counts = gen_counts();


double binsize = 0.25;

arma::vec I = periodogram(::counts, ::binsize, 5);


void function1_func(const real_1d_array& x, double& func, void* ptr) {

    double x0 = x[0];
    double x1 = x[1];
    double x2 = x[2];
    arma::vec param = arma::vec({ x0, x1, x2 });
    Model Expo(param, ::binsize);
    func = Expo.whittle(::I, 5);
}




int main(int argc, char** argv)
{


    auto start = std::chrono::high_resolution_clock::now();

    arma::vec init = init_bfgs(::counts, ::binsize, 5);

    double x0 = init[0];
    double x1 = init[1];
    double x2 = init[2];
    real_1d_array x;
    double content[] = { x0, x1, x2 };
    x.setcontent(3, content);

    real_1d_array s = "[1,1,1]";
    real_1d_array bndx0 = "[0.01,0.01, 0.01]";
    real_1d_array bndx1 = "[100,0.99, 100]";
    minbleicstate state;
    double epsg = 0;
    double epsf = 0;
    double epsx = 0.000001;
    ae_int_t maxits = 0;
    double diffstep = 1.0e-6;

    minbleiccreatef(x, diffstep, state);
    minbleicsetbc(state, bndx0, bndx1);
    minbleicsetscale(state, s);
    minbleicsetcond(state, epsg, epsf, epsx, maxits);


    //
    minbleicoptguardsmoothness(state);

    //
    // Optimize and evaluate results
    //
    minbleicreport rep;
    minbleicoptimize(state, function1_func);
    minbleicresults(state, x, rep);
    printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
    printf("%s\n", x.tostring(3).c_str()); // EXPECTED: [-1,1]


    optguardreport ogrep;
    minbleicoptguardresults(state, ogrep);
    printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
    printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false





    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_ms = end - start;
    cout << "elapsed time is " << float_ms.count() << " milliseconds" << std::endl;

}


