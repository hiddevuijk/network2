#ifndef GUARD_GENERATE_GRAPH_H
#define GUARD_GENERATE_GRAPH_H

#include "graph.h"

#include "boost/random.hpp"

#include <math.h>
#include <iostream>



// calculate vertex index from xi, yi
int xy2n(int xi, int yi, int Nx, int Ny)
{ return yi*Nx + xi; }

int n2x(int n, int Nx)
{   return n%Nx; }

int n2y( int n, int Nx)
{ return (n - n%Nx)/Nx; }

// vertex index from West neighbor of vertex at (xi, yi)
int niNeighborW(int xi, int yi, int Nx, int Ny)
{ return xy2n( (xi + Nx - 1) % Nx, yi, Nx, Ny); }

int niNeighborE(int xi,int yi, int Nx, int Ny)
{ return xy2n( (xi + 1) % Nx, yi, Nx, Ny); }

int niNeighborNE(int xi, int yi, int Nx, int Ny)
{ return xy2n( (xi + yi%2)%Nx, (yi+1)%Ny, Nx,Ny); }

int niNeighborNW(int xi, int yi, int Nx, int Ny)
{ return xy2n( (xi  +Nx -1 + yi%2)%Nx, (yi+1)%Ny , Nx, Ny);}

int niNeighborSE(int xi, int yi, int Nx, int Ny)
{ return xy2n( (xi + yi%2)%Nx  , (yi + Ny -1)%Ny  ,Nx, Ny); }

int niNeighborSW(int xi, int yi, int Nx, int Ny)
{ return xy2n( (xi + Nx -1 + yi%2)%Nx  , (yi+Ny-1)%Ny, Nx, Ny); }

// neighbor vertex index from vertex index
int niNeighborW(int ni, int Nx, int Ny)
{ return niNeighborW( n2x(ni,Nx), n2y(ni,Nx), Nx, Ny); }

int niNeighborE(int ni, int Nx, int Ny)
{ return niNeighborE( n2x(ni,Nx), n2y(ni,Nx), Nx, Ny); }

int niNeighborNE(int ni, int Nx, int Ny)
{ return niNeighborNE( n2x(ni,Nx), n2y(ni,Nx), Nx, Ny); }

int niNeighborNW(int ni, int Nx, int Ny)
{ return niNeighborNW( n2x(ni,Nx), n2y(ni,Nx), Nx, Ny); }

int niNeighborSE(int ni, int Nx, int Ny)
{ return niNeighborSE( n2x(ni,Nx), n2y(ni,Nx), Nx, Ny); }

int niNeighborSW(int ni, int Nx, int Ny)
{ return niNeighborSW( n2x(ni,Nx), n2y(ni,Nx), Nx, Ny); }



template<class V>
void shuffle( V& v, boost::random::mt19937 & rng)
{
    for( long unsigned int i=0;i<v.size(); ++i) {
        boost::random::uniform_int_distribution<int> rand(i, v.size()-1 );
        int j = rand(rng);
        auto temp = v[i];
        v[i] = v[j];
        v[j] = temp;
    }
}




Graph generateGraph(int Nx, int Ny, double Lx, double z, long int seed = 123456789, double sigma = 0)
{
    boost::random::mt19937 rng(seed);
    boost::random::normal_distribution<double> ndist(0,1);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > rndist(rng,ndist);

    boost::random::uniform_int_distribution<int> randX(0,Nx-1);
    boost::random::uniform_int_distribution<int> randY(0,Ny-1);

    // Ny must be even!!

    int Nnodes = Nx*Ny; // number of vertices
    Graph graph(Nnodes);
     
    double dx = Lx/Nx;
    double dy = dx*std::sqrt(3./4.);

    int ni, niNeighbor;
    int xb, yb;
    for(int xi=0; xi<Nx; ++xi) {
    for(int yi=0; yi<Ny; ++yi) {
        ni = xy2n(xi,yi,Nx,Ny);
        graph.setNodePosition( ni, xi*dx + (yi%2)*0.5*dx, yi*dy);

        // add East edge
        niNeighbor = niNeighborE(xi,yi,Nx,Ny);
        if( xi == Nx -1 ) xb = 1;
        else xb = 0;
        yb = 0;
        
        graph.addBond( ni, niNeighbor, xb, yb);

        // add South-East edge
        niNeighbor = niNeighborSE(xi,yi,Nx,Ny);
        if( xi == Nx -1 and yi % 2 == 1) xb = 1; 
        else xb = 0;
        if( yi == 0 ) yb = -1;
        else yb = 0;
        
        graph.addBond( ni, niNeighbor, xb, yb);

        // add South-West edge
        niNeighbor = niNeighborSW(xi,yi,Nx,Ny);
        if( xi == 0 and yi % 2 == 0 ) xb = -1; 
        else xb = 0;

        if( yi ==0 ) yb = -1;
        else yb = 0;
        
        graph.addBond( ni, niNeighbor, xb, yb);


    }}

    // add Bends
    int niPrev, niNext;

    // add W-E bends
    for(int yi=0; yi<Ny; ++yi ){
        ni = xy2n(0,yi,Nx,Ny);
        for(int xi = 0; xi < Nx; ++xi ) {
            ni = xy2n(xi,yi,Nx,Ny); 
            niPrev = niNeighborW(xi,yi,Nx,Ny);
            niNext = niNeighborE(xi,yi,Nx,Ny);
            graph.addBend(ni, niPrev, niNext);
        } 
        for(int xi=0; xi<Nx; ++xi ) {
            ni = xy2n(xi,yi,Nx,Ny); 
            niNext = niNeighborE(xi,yi,Nx,Ny);
            graph.polymerize(ni, niNext);
        }

    }

    // add SW-NE bends
    for(int yi=0; yi<Ny; ++yi ){
    for(int xi=0; xi<Nx; ++xi ) {
        ni = xy2n(xi,yi,Nx,Ny);
        niPrev = niNeighborSW(ni,Nx,Ny);  
        niNext = niNeighborNE(ni,Nx,Ny);  
        graph.addBend(ni, niPrev, niNext);

    }}

    //polymerize SW-NE bends
    for(int yi=0; yi<Ny; ++yi ){
    for(int xi=0; xi<Nx; ++xi ) {
        ni = xy2n(xi,yi,Nx,Ny);
        niNext = niNeighborNE(ni,Nx,Ny);  
        graph.polymerize(ni, niNext);
    }}

    // add SE-NW bends
    for(int yi=0; yi<Ny; ++yi ){
    for(int xi=0; xi<Nx; ++xi ) {
        ni = xy2n(xi,yi,Nx,Ny);
        niPrev = niNeighborSE(ni,Nx,Ny);  
        niNext = niNeighborNW(ni,Nx,Ny);  
        graph.addBend(ni, niPrev, niNext);

    }}
    //polymerize SE-NW bends
    for(int yi=0; yi<Ny; ++yi ){
    for(int xi=0; xi<Nx; ++xi ) {
        ni = xy2n(xi,yi,Nx,Ny);
        niNext = niNeighborNW(ni,Nx,Ny);  
        graph.polymerize(ni, niNext);
    }}


    // cut every filament once
    // cut WE filament
    for(int yi=0; yi<Ny; ++yi) {
        int xi = randX(rng);
        int ni = xy2n(xi,yi,Nx,Ny);
        int niNext = niNeighborE(ni,Nx,Ny);
        graph.delBond(ni,niNext);
    }

    // cut SE-NW filament
    for(int xi=0; xi<Nx; ++xi) {
        int yi = randY(rng);
        int ni = xy2n(xi,0, Nx,Ny);
        while( yi > 0 ) {
            ni = niNeighborNW(ni,Nx,Ny);
            --yi;
        }
        int niNext = niNeighborNW(ni,Nx,Ny);
        graph.delBond(ni,niNext);
    }

    // cut SW-NE filament
    for(int xi=0; xi<Nx; ++xi) {
        int yi = randY(rng);
        int ni = xy2n(xi,0, Nx,Ny);
        while( yi > 0 ) {
            ni = niNeighborNE(ni,Nx,Ny);
            --yi;
        }
        int niNext = niNeighborNE(ni,Nx,Ny);
        graph.delBond(ni,niNext);
    }

    std::vector<std::vector<int> > bonds = graph.getBonds();
    shuffle(bonds, rng);



    int i=0;
    // edges removed by pruning are not removed from edges list!!
    while( graph.getConnectivity2() > z ) {
        graph.delBond( bonds[i][0], bonds[i][1] );
		graph.prune( bonds[i][0]);
		graph.prune( bonds[i][1]);
        ++i;
    }
    graph.removeUnconnectedNodes();

    // shake vertices
    //vec2 r;
    //for(std::vector<Graph::Node*>::size_type ni = 0; ni < graph.Nnodes(); ++ni) {
    //    r = graph.getNodePosition(ni); 
    //    r.x += sigma*rndist();
    //    r.y += sigma*rndist();
    //    graph.setNodePosition(ni, r.x, r.y);
    //}

    return graph;

}




#endif
