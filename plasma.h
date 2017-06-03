#ifndef PLASMA_H
#define PLASMA_H

#include <fftw3.h>
#include <ctime>
#include <eigen3/Eigen/Dense>
#include <cmath>

#include <QDebug>
using namespace Eigen;

class plasma{
public:
    explicit plasma(float l, int j, int n);
	VectorXf eval(float tmax, float dt);
    void vDist(float vb);
    void rDist();
    void setN(int l);
    void setL(float l);
    void setJ(int j);

private:
    void density(VectorXf r0, VectorXf& n);
    void fftw_forward(VectorXf f, VectorXcf& F);
    void fftw_backward(VectorXcf F, VectorXf& f);
    void poisson1d(VectorXf v, VectorXf& u, float kappa);
    void electric(VectorXf phi, VectorXf& E);
    void rk4_fixed_vector(float& t, VectorXf& y, float h);
    void function(float t, VectorXf y, VectorXf& dydt);
    int N, J;
    float L;
    VectorXf v, r;
    MatrixXf plasVals;
};

#endif // PLASMA_H
