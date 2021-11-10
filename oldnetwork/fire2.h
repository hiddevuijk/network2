#ifndef GUARD_FIRE2_H
#define GUARD_FIRE2_H

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
	void minimizeEE( const std::vector<double>&);
	void minimizeSIE( const std::vector<double>&);

	void make_VV_step(); 
	void make_EE_step(); 
	void make_SIE_step(); 

	void make_FIRE_step();
    void set_Fnorm() {
        Fnorm = 0;
        for(int i=0; i<N; ++i) Fnorm += F[i]*F[i];
        Fnorm = std::sqrt(Fnorm);
    }

    double Fmax() {
        double max = 0;
        for(int i=0;i<N; ++i) {
             max = std::fabs(F[i]) > max ? std::fabs(F[i]) : max;
        }
        return max;
    }
	int N; // number of variables
	E *energy;

	int Nmin;
	double finc;
	double fdec;

	double alpha0;
	double falpha;

	double m;

	double dtmax;
    double dtmin;

	double dt0;
	double error;
    double error_max;

	// changes	
	std::vector<double> x; // state 
	std::vector<double> v; // "speed"
	std::vector<double> F; // - nabla energy

	double alpha;
	double dt;
	int NPneg;
	double Fnorm;
    int NF;
    

};

template<typename E>
Fire<E>::Fire( int N, E *energy )
: N(N), energy(energy), x(N), v(N), F(N) { }

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
    NF = 0;
}

template<typename E>
void Fire<E>::minimizeVV( const std::vector<double>& xInit)
{
	initialize(xInit);

	int it  =0 ;
	do {
		make_FIRE_step();
		make_VV_step();
        set_Fnorm();
		it++;
        NF++;
	} while(Fnorm > N*error or Fmax() > error_max );
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
        set_Fnorm();
		it++;
        NF++;
	} while(Fnorm > N*error or Fmax() > error_max );
	//std::cout << "\t\t" << it << std::endl;
}

template<typename E>
void Fire<E>::minimizeSIE( const std::vector<double>& xInit)
{
	initialize(xInit);

	int it  =0 ;
	do {
        //std::cout << it << '\t' <<  dt << std::endl;
		make_FIRE_step();
		make_EE_step();
        set_Fnorm();
		it++;
        NF++;
	} while(Fnorm > N*error or Fmax() > error_max );
    //} while( Fmax() > error );
}



template<typename E>
void Fire<E>::make_VV_step()
{
	double Fnorm = 0;
	double vnorm = 0;
	for(int i=0;i<N;++i) {
		v[i] += 0.5*dt*F[i]/m;
		Fnorm += F[i]*F[i];
		vnorm += v[i]*v[i];
	}
	Fnorm = std::sqrt(Fnorm);
	vnorm = std::sqrt(vnorm);

	for(int i=0; i<N; ++i) {
		v[i] += -alpha*v[i] + alpha*F[i]*vnorm/Fnorm;
		x[i] += dt*v[i];
	}

	energy->dE(F,x);

	for(int i=0;i<N; ++i) {
		v[i] += 0.5*dt*F[i]/m;
	}
}


template<typename E>
void Fire<E>::make_EE_step()
{
	double Fnorm = 0;	
	double vnorm = 0;	
	for(int i=0; i<N; ++i) {
		Fnorm += F[i]*F[i];
		vnorm += v[i]*v[i];
	}
	Fnorm = std::sqrt(Fnorm);
	vnorm = std::sqrt(vnorm);

	for(int i=0;i<N; ++i) {
		v[i] += -alpha*v[i] + alpha*F[i]*vnorm/Fnorm;
		x[i] += dt*v[i];
	}
	energy->dE(F,x);
	for(int i=0;i<N; ++i) {
		v[i] += dt*F[i]/m;	
	}

}


template<typename E>
void Fire<E>::make_SIE_step()
{
	double Fnorm = 0;
	double vnorm = 0;
	for(int i=0;i<N; ++i) {
		v[i] += dt*F[i]/m;
		Fnorm += F[i]*F[i];
		vnorm += v[i]*v[i];
	}
	Fnorm = std::sqrt(Fnorm);
	vnorm = std::sqrt(vnorm);

	for(int i=0; i<N; ++i) {
		v[i] += -alpha*v[i] + alpha*F[i]*vnorm/Fnorm;
		x[i] += dt*v[i];
	}

	energy->dE(F,x);

}



template<typename E>
void Fire<E>::make_FIRE_step()
{

	double P = 0;	
	for(int i=0;i<N; ++i) {
		P += v[i]*F[i];
	}

	if( P > 0 ) {
		NPneg++;
		if(NPneg > Nmin ) {
			dt = std::min(dt*finc, dtmax);
			alpha *= falpha;
		}

	} else { // P <= 0
		NPneg = 0;

		for(int i=0;i<N; ++i) {
			x[i] -= 0.5*dt*v[i];
		}

		//dt *= fdec;
        dt = std::max( dt*fdec, dtmin);
		alpha = alpha0;

		std::fill( v.begin(), v.end(), 0.0);
	}
	
}


#endif
