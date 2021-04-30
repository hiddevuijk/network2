#ifndef GUARD_GENERATE_NETWORK_H
#define GUARD_GENERATE_NETWORK_H



#include "graph.h"

#include "boost/random.hpp"

#include <math.h>
#include <iostream>

// calculate vertex index from xi, yi
int xy2v(int xi, int yi, int Nx, int Ny)
{ return yi*Nx + xi; }

int v2x(int v, int Nx)
{   return v%Nx; }

int v2y( int v, int Nx)
{ return (v - v%Nx)/Nx; }

// vertex index from West neighbor of vertex at (xi, yi)
int viNeighborW(int xi, int yi, int Nx, int Ny)
{ return xy2v( (xi + Nx - 1) % Nx, yi, Nx, Ny); }

int viNeighborE(int xi,int yi, int Nx, int Ny)
{ return xy2v( (xi + 1) % Nx, yi, Nx, Ny); }

int viNeighborNE(int xi, int yi, int Nx, int Ny)
{ return xy2v( (xi + yi%2)%Nx, (yi+1)%Ny, Nx,Ny); }

int viNeighborNW(int xi, int yi, int Nx, int Ny)
{ return xy2v( (xi  +Nx -1 + yi%2)%Nx, (yi+1)%Ny , Nx, Ny);}

int viNeighborSE(int xi, int yi, int Nx, int Ny)
{ return xy2v( (xi + yi%2)%Nx  , (yi + Ny -1)%Ny  ,Nx, Ny); }

int viNeighborSW(int xi, int yi, int Nx, int Ny)
{ return xy2v( (xi + Nx -1 + yi%2)%Nx  , (yi+Ny-1)%Ny, Nx, Ny); }

// neighbor vertex index from vertex index
int viNeighborW(int vi, int Nx, int Ny)
{ return viNeighborW( v2x(vi,Nx), v2y(vi,Nx), Nx, Ny); }

int viNeighborE(int vi, int Nx, int Ny)
{ return viNeighborE( v2x(vi,Nx), v2y(vi,Nx), Nx, Ny); }

int viNeighborNE(int vi, int Nx, int Ny)
{ return viNeighborNE( v2x(vi,Nx), v2y(vi,Nx), Nx, Ny); }

int viNeighborNW(int vi, int Nx, int Ny)
{ return viNeighborNW( v2x(vi,Nx), v2y(vi,Nx), Nx, Ny); }

int viNeighborSE(int vi, int Nx, int Ny)
{ return viNeighborSE( v2x(vi,Nx), v2y(vi,Nx), Nx, Ny); }

int viNeighborSW(int vi, int Nx, int Ny)
{ return viNeighborSW( v2x(vi,Nx), v2y(vi,Nx), Nx, Ny); }





Graph generateNetwork(int Nx, int Ny, double Lx)
{
    long int seed = 123456789;
    boost::random::mt19937 rng(seed);
    boost::random::uniform_int_distribution<int> randX(1,Nx-1);
    boost::random::uniform_int_distribution<int> randY(1,Ny-1);

    // Ny must be even!!

    int Nv = Nx*Ny; // number of vertices

    Graph graph(Nv);

    double dx = Lx/Nx;
    double dy = dx*std::sqrt(3./4.);

    int vi; 
    for(int xi=0; xi<Nx; ++xi) {
    for(int yi=0; yi<Ny; ++yi) {
        vi = xy2v(xi,yi,Nx,Ny);
        graph.set_position( vi, xi*dx, yi*dy);
        graph.addEdge( vi, viNeighborE(xi,yi, Nx, Ny));
        graph.addEdge( vi, viNeighborSE(xi,yi, Nx, Ny));
        graph.addEdge( vi, viNeighborSW(xi,yi, Nx, Ny));
    }}


    // add bends 
    int v0;
    int viPrev;
    int viNext;

    // add W-E bends
    for(int yi =0; yi < Ny; ++yi ){
        // filament starts at vi
        vi = xy2v(0,yi,Nx,Ny); 
    

        for(int xi = 0; xi < Nx; ++xi ) {
            vi = xy2v(xi,yi,Nx,Ny); 
            viPrev = viNeighborW(xi,yi,Nx,Ny);
            viNext = viNeighborE(xi,yi,Nx,Ny);
            graph.addBend(vi, viPrev, viNext);
        } 
        for(int xi=0; xi<Nx; ++xi ) {
            vi = xy2v(xi,yi,Nx,Ny); 
            viNext = viNeighborE(xi,yi,Nx,Ny);
            graph.polymerize(vi, viNext);
        }
    }

    // add SW-NE bends
    for(int xi=0; xi<Nx; ++xi ) {
        v0 = xy2v(xi,0,Nx,Ny);
        vi = v0;
        viPrev = viNeighborSW(vi,Nx,Ny);  
        viNext = viNeighborNE(vi,Nx,Ny);  
        // check if it already exists 
        // if it exists it is the  last one added
        if( graph.vertices[vi].bends.size() > 0 and graph.vertices[vi].bends.back().next == viNext ) break;
       
        graph.addBend(vi, viPrev, viNext);
        while( viNext != v0 ) {
            viPrev = vi;
            vi = viNext;
            viNext = viNeighborNE(vi,Nx,Ny);  
            //graph.addBend(vi, viPrev, viNext);
        } 
    }

    //// add SE-NW bends    
    //for(int xi=0; xi<Nx; ++xi) {
    //    v0 = xy2v(xi,0,Nx,Ny);
    //    vi = v0;
    //    viPrev = viNeighborSE(vi,Nx,Ny);
    //    viNext = viNeighborNW(vi,Nx,Ny);

    //    // check if it already exists 
    //    // if it exists it is the  last one added
    //    if( graph.vertices[vi].bends.size() > 0 and graph.vertices[vi].bends.back().next == viNext ) break;

    //    // add a filament that starts at vi to the filament list
    //    graph.addBend(vi, viPrev, viNext);
    //    while(viNext != v0) {
    //        viPrev = vi;
    //        vi = viNext;
    //        viNext = viNeighborNW(vi,Nx,Ny);
    //        graph.addBend(vi, viPrev, viNext);
    //    }
    //}
       

   return graph;

}




#endif
