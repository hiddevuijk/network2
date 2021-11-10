#ifndef GUARD_FIRE_H
#define GUARD_FIRE_H

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>


// the ebject E must have a method dE( vector<double>& F, const vector<double> &x)
// where F is the container for the calculated force
// and x is the vector describing the state of the system.
template<typename E>
class Fire {
public:
	Fire() {}
	Fire( int N, E *energy );

	void initialize( const std::vector<double>&);
	void minimizeVV( const std::vector<double>&);
	void minimizeVV2( const std::vector<double>&);
	void minimizeEE( const std::vector<double>&);
	void minimizeSIE( const std::vector<double>&);
	void make_VV_step(); 
	void make_VV2_step(); 
	void make_EE_step(); 
	void make_SIE_step(); 
	void make_FIRE_step();

	int N; // number of variables
	E *energy;

	int Nmin;
	double finc;
	double fdec;
	double alpha0;
	double falpha;
	double m;

	double dtmax;
	double dt0;
	double error;

	// changes	
	std::vector<double> x; // state 
	std::vector<double> v; // "speed"
	std::vector<double> F; // - nabla energy

	double alpha;
	double dt;
	double NPneg;
	double Fnorm;

};

template<typename E>
Fire<E>::Fire( int N, E *energy )
: N(N), energy(energy), x(N), v(N), F(N)
{

	Nmin = 5;
	finc = 1.1;
	fdec = 0.5;
	alpha0 = 0.25;
	falpha = 0.99;
	m = 1e-1;

	dtmax = 0.1;
	dt0 = 0.03;
	error = 1e-6;


}

template<typename E>
void Fire<E>::initialize( const std::vector<double>& xInit)
{
	for(int i=0; i< N; ++i) x[i] = xInit[i];

	// calculate F(t)
	energy->dE(F,x);	
	// set v[i] = 0
	std::fill(v.begin(), v.end(), 0.0);
	dt = dt0;
	alpha = alpha0;
	NPneg = 0;
}

template<typename E>
void Fire<E>::minimizeVV( const std::vector<double>& xInit)
{
	initialize(xInit);

	int it  =0 ;
	do {
		make_FIRE_step();
		make_VV_step();
		it++;
	} while(Fnorm > error );
	//std::cout << "\t\t" << it << std::endl;
}

template<typename E>
void Fire<E>::minimizeVV2( const std::vector<double>& xInit)
{
	initialize(xInit);

	int it  =0 ;
	do {
		make_FIRE_step();
		make_VV_step();
		it++;
	} while(Fnorm > error );
	//std::cout << "\t\t" << it << std::endl;
}

template<typename E>
void Fire<E>::minimizeEE( const std::vector<double>& xInit)
{
	initialize(xInit);

	int it  =0 ;
	do {
		make_FIRE_step();
		make_EE_step();
		it++;
	} while(Fnorm > error );
	//std::cout << "\t\t" << it << std::endl;
}

template<typename E>
void Fire<E>::minimizeSIE( const std::vector<double>& xInit)
{
	initialize(xInit);

	int it  =0 ;
	do {
		make_FIRE_step();
		make_EE_step();
		it++;
	} while(Fnorm > error );
	//std::cout << "\t\t" << it << std::endl;
}



template<typename E>
void Fire<E>::make_VV_step()
{
	// assumes F(t) is correct

	// set x(t+dt)
	for(int i=0; i<N; ++i) {
		x[i] += v[i]*dt + 0.5*F[i]*dt*dt/m;
	}

	// set v(t+dt) with F(t)
	for(int i=0; i<N; ++i) {
		v[i] += 0.5*F[i]*dt/m;
	}


	// calculate F(t+dt)
	energy->dE(F,x);

	// set v(t+dt) with F(t+dt)
	for(int i=0; i<N; ++i) {
		v[i] += 0.5*F[i]*dt/m;
	}

}

template<typename E>
void Fire<E>::make_VV2_step()
{
	// assumes F(t) is correct

	// set v(t+dt/2) with F(t)
	for(int i=0; i<N; ++i) {
		v[i] += 0.5*F[i]*dt/m;
		x[i] += dt*v[i];
	}

	// calculate F(t+dt)
	energy->dE(F,x);

	// set v(t+dt) with F(t+dt)
	for(int i=0; i<N; ++i) {
		v[i] += 0.5*F[i]*dt/m;
	}

}


template<typename E>
void Fire<E>::make_EE_step()
{
	// assumes F(t) is correct

	// set x(t+dt)
	for(int i=0; i<N; ++i) {
		x[i] += v[i]*dt;
	}

	// calculate F(t+dt)
	energy->dE(F,x);

	// set v(t+dt) with F(t+dt)
	for(int i=0; i<N; ++i) {
		v[i] += F[i]*dt/m;
	}

}


template<typename E>
void Fire<E>::make_SIE_step()
{
	// assumes F(t) is correct

	// set x(t+dt)
	for(int i=0; i<N; ++i) {
		v[i] += dt*F[i]/m;
		x[i] += v[i]*dt;
	}

	// calculate F(t+dt)
	energy->dE(F,x);
}



template<typename E>
void Fire<E>::make_FIRE_step()
{

	// step F1
	Fnorm = 0;
	double vnorm = 0;
	double P = 0;	
	for(int i=0;i<N; ++i) {
		Fnorm += F[i]*F[i];
		vnorm += v[i]*v[i];
		P += v[i]*F[i];
	}
	Fnorm = std::sqrt(Fnorm);
	vnorm = std::sqrt(vnorm);


	

	// step F3
	if( P <= 0 ) {
		dt *= fdec;
		alpha = alpha0;
		// v[i] = 0
		std::fill( v.begin(), v.end(), 0.);
		NPneg = 0;
	} else {
		NPneg++;
		for(int i=0;i<N; ++i) {
			v[i] = (1-alpha)*v[i] + alpha*F[i]*vnorm/Fnorm;
		}

		if( NPneg > Nmin ) {
			dt = std::min(dt*finc, dtmax);
			alpha *= falpha;
		}
	} 
}


#endif
