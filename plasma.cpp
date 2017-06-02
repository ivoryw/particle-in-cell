#include "plasma.h"

const float PI = 3.141592865358979;

plasma::plasma(float l, int j, int n){
	L = l;
	J = j;
	N = n;
	v.resize(N);
	r.resize(N);
}

void plasma::vDist(float vb){
	int i = 0;
	float vmin, vmax, vi, Pv, PvMax, x;
	while(i<N){
		static int flag = 0;
		if(flag==0){
			srand(time(NULL));
			flag = 1;
		}
		vmin = -5 *vb;
		vmax = 5 * vb;
		vi = vmin + (vmax-vmin) * float(rand()) / float(RAND_MAX);

		Pv = 0.5*(exp(-1*(vi-vb)*(vi-vb)/2) + exp(-1*(vi+vb)*(vi+vb)/2));
		PvMax = (1+exp(-2*vb*vb))/2;
		x = PvMax * float(rand()) / float(RAND_MAX);
		if(Pv < x)
			continue;
		else
			v(i) = vi;
		i++;
	}
}

void plasma::rDist(){
    r = (VectorXf::Random(N) + VectorXf::Constant(N,1)) * L/2;
}
void plasma::setJ(int j){
	J = j;
}
void plasma::setL(float l){
	L =l;
}
void plasma::setN(int n){
	N = n;
	v.resize(N);
	r.resize(N);
}

VectorXf plasma::eval(float tmax, float dt){
	float t = 0;
    VectorXf y(2*N);
//    plasVals.resize(20,2*N);
	y << r,v;
    do{
		rk4_fixed_vector(t, y, dt);
        for(int i=0; i<N; i++){
            y(i) = fmod(y(i), L);
            if(y(i) < 0)
                y(i) += L;
            // This catches cases where r0(i) = L
            else if(y(i) == L)
                y(i) -=L;
        }
//        if(fmod((t-dt)/dt,0.1) == 0)
//            plasVals.col((t-dt)/dt) = y;
	}while(t<=tmax);
//    for(int i = 0; i<2*N;i++)
//        qDebug() << plasVals(i);
    return y;
}

void plasma::rk4_fixed_vector(float&t, VectorXf&y, float dt){
	int I = N*2;
	VectorXf k1(I), k2(I), k3(I), k4(I), dydt(I);

	function(t, y, dydt);
	k1 = dt * dydt;
	dydt = VectorXf::Zero(I);
	function(t + dt/2,y + k1/2, dydt);
	k2 = dt * dydt;
	dydt = VectorXf::Zero(I);
	function(t + dt/2, y +k2/2, dydt);
	k3 = dt * dydt;
	dydt = VectorXf::Zero(I);
	function(t + dt, y + k3, dydt);
	k4 = dt * dydt;
	y += (k1+k4)/6 + (k2+k3)/3;
	t += dt;
}

void plasma::function(float t, VectorXf y, VectorXf& dydt){
	VectorXf Efield(N), dvdt(N), drdt(N), r0(N), v0(N);
	VectorXf n(J), rho(J), E(J), phi(J);
	float n0, kappa,dx;

	r0 = y.head(N);
	v0 = y.tail(N);

    for(int i=0; i<N; i++){
        r0(i) = std::fmod(r0(i),L);
        if(r0(i) < 0)
            r0(i) += L;
    }

	density(r0, n);

	n0 = float(N)/L;
	kappa = 2 * PI / L;
	rho = n / (n0 - 1);

	poisson1d(rho, phi, kappa);

	electric(phi, E);

	dx = L/float(J);
	for(int i = 0; i<N; i++){
		int j = r0(i)/dx;
        if(j == J)
            j = 0;
        float val = r0(i)/dx - j;

		if(j == J-1)
            Efield(i) = E(j)*(1-val) + E(0) *val;
		else
            Efield(i) = E(j)*(1-val) + E(j+1) * val;
	}
	drdt = v0;
	dvdt = -1 * Efield;
	dydt << drdt, dvdt;
}

void plasma::density(VectorXf r0, VectorXf& n){
	float dx = L/J;
	int j;
    n = VectorXf::Zero(J);

	for(int i = 0; i < N; i++){
        j = r0(i)/dx;
        if(j == J)
            j = 0;
        n(j) += (dx-r0(i) + float(j)*dx)/ (dx*dx);
        if(j == J-1)
            n(0) += (r0(i) - float(j)*dx)/(dx*dx);
        else
            n(j+1)+=(r0(i) - float(j)*dx)/(dx*dx);
    }
}

void plasma::poisson1d(VectorXf rho, VectorXf& phi, float kappa){
	VectorXcf V(J), U(J);
    fftw_forward(rho, V);

	U(0) = 0;
	for(int j = 1; j <= J/2; j++)
        U(j) = float(-1.) * V(j) / (j*j*kappa*kappa);
	for(int j = J/2; j<J; j++)
		U(j) = conj(U(J-j));

    fftw_backward(U,phi);
}

void plasma::fftw_forward(VectorXf f, VectorXcf& F){
	fftw_complex in[J], out[J];

	for(int j = 0; j<J; j++){
		in[j][0] = f(j); in[j][1] = 0;
	}

	fftw_plan p = fftw_plan_dft_1d(J, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);

	for(int j=0; j<J; j++){
		F(j).real(out[j][0]); F(j).imag(out[j][1]);
	}
	F /= J;
}

void plasma::fftw_backward(VectorXcf F, VectorXf& f){
	fftw_complex in[J], out[J];

    for(int j=0; j<J; j++){
		in[j][0] = real(F(j)); in[j][1] = imag(F(j));
	}

	fftw_plan p = fftw_plan_dft_1d(J, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);

	for(int j=0; j<J; j++)
		f(j) = out[j][0];	
}

void plasma::electric(VectorXf phi, VectorXf& E){
	float dx = L/J;
	for(int j=1; j<J-1; j++)
		E(j) = (phi(j-1) -phi(j+1))/(2*dx);
	E(0) = (phi(J-1) - phi(1))/(2*dx);
	E(J-1) = (phi(J-2) - phi(0))/(2*dx);
}
