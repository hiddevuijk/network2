#ifndef GUARD_NETWORK_H
#define GUARD_NETWORK_H


/*
    to do:
            -- destructors
            -- throw error when polymerizing nodes without a bend
*/

#include "vec2.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

class Network {
  public:

    class Vertex;
    class Edge;
    class Bend;


    Network();
    Network(int Nv);
    // Network ( infile) read network from file
    ~Network();

    void write( std::ofstream& out);
    int addVertex();
    int addVertex(double x, double y);

    void setVertexPosition(int i, double x, double y);

    void addEdge(int i, int j);
    void deleteEdge(int i, int j);

    void addBend(int mid, int prev, int next);
    void deleteBend(int mid, int prev, int next); // --
    void deleteBend(int vi, int vj); // delete all Bends with edge vi-vj

    void polymerize(int i, int j);
    void depolymerize(int i, int j); // test 
    void depolymerize(int i); // test // remove from all polymers

    // this changes the indices!!!
    void removeUnconectedVertices(); // --
    // ....

    std::vector<std::vector<int> > getEdges() const;
    std::vector<std::vector<int> > getBends() const; // ---
    std::vector<std::vector<int> > getPolymers() const; // ---

    void showAdj() const;
    void showBends() const;
    // print polymer starting at i through j
    void showPolymer(int i, int j) const;


    class Vertex {
      public:
        Vertex(): index(-1) {}
        Vertex(int i) : index(i) {}
        Vertex(int i, double x, double y) : index(i), r(x,y) {}
        // destructor ----
        int index;
        std::vector<Vertex*> adj;
        std::vector<Edge*> edges;
        std::vector<Bend*> bends;

        // add location info
        vec2 r;
    };

    class Edge {
      public:
        Vertex *from, *to;
        int xBoundary, yBoundary;
        double l0; 
    };


    class Bend {
      public:
        Bend()
            : mid(nullptr), prev(nullptr), next(nullptr),
              prevBend(nullptr), nextBend(nullptr),
              filament(-1) {}

        Bend(Vertex *mid, Vertex *prev, Vertex *next)
            : mid(mid), prev(prev), next(next),
              prevBend(nullptr), nextBend(nullptr),
              filament(-1) {}
        Vertex *mid;
        Vertex *prev;
        Vertex *next;
        Bend *prevBend;
        Bend *nextBend;
        int filament;
        double theta0;
    };

    private:
        std::vector<Vertex*> vertices;

        void deleteVertex(int i); // test
        void deleteVertex(Vertex *vi); // test
        void exchangeVertex(int i, int j); // test
        void exchangeVertex(Vertex *vi, Vertex *vj); // test

        void addEdge(Vertex *vi, Vertex *vj);
        void deleteEdge(Vertex *vi, Vertex *vj);

        void addBend(Vertex *mid, Vertex *prev, Vertex *next);
        void deleteBend(Vertex *mid, Vertex *prev, Vertex *next); // ---
        void deleteBend(Vertex *vi, Vertex *vj); // --
        void deleteBend(Bend *bend); // --
        void deleteBend( std::vector<Bend*>::iterator it_bend);

        void polymerize(Vertex *vi, Vertex *vj);
        void polymerize(Bend *bi, Bend *bj); 

        void depolymerize( Vertex *vi, Vertex *vj); // --
        void depolymerize( Vertex *vi); // test // remove from all polymers
        void depolymerize( Bend *bi, Bend *bj); // test
        void depolymerize( Bend *bi); // test // remove from all polymers

};


Network::Network() {};

Network::Network( int Nv)
: vertices( std::vector<Vertex*>(Nv, nullptr) )
{
    for( int vi=0; vi<Nv; ++vi) {
        vertices[vi] = new Vertex(vi);
    }
}


Network::~Network()
{ 
    for( std::vector<Vertex*>::size_type vi=0; vi < vertices.size(); ++vi) {
        delete vertices[vi];
        vertices[vi] = nullptr;
    }
}

int Network::addVertex()
{
    std::vector<Vertex*>::size_type  Nv = vertices.size();
    vertices.push_back( new Vertex(Nv) );
    return Nv;
}

int Network::addVertex(double x, double y)
{
    std::vector<Vertex*>::size_type  Nv = vertices.size();
    vertices.push_back( new Vertex(Nv, x, y) );
    return Nv;
}

void Network::setVertexPosition(int i, double x, double y)
{
    vertices[i]->r.x = x;
    vertices[i]->r.y = y;
}

void Network::deleteVertex(int i)
{
    // remove all edges
    while( vertices[i]->adj.size() > 0 ) {
        deleteEdge( vertices[i], vertices[i]->adj.back() ); 
        vertices[i]->adj.pop_back();
    }

    std::vector<Vertex*>::size_type Nv = vertices.size();
    exchangeVertex(i,Nv);
    // delete vertex

    vertices.pop_back();
}

void Network::deleteVertex(Vertex *vi)
{ deleteVertex( vi->index ); }

void Network::exchangeVertex(int i, int j)
{
    Vertex *temp = vertices[i];
    vertices[i] = vertices[j];
    vertices[i]->index = i;
    vertices[j] = temp;
    vertices[j]->index = j;
}

void Network::exchangeVertex(Vertex *vi, Vertex *vj )
{ exchangeVertex( vi->index, vj->index ); }


void Network::addEdge(int i, int j)
{ addEdge(vertices[i], vertices[j] ); }

void Network::addEdge(Vertex *vi, Vertex *vj)
{
    // check if already exists??
    vi->adj.push_back( vj ); 
    vj->adj.push_back( vi ); 
}

void Network::deleteEdge(int i, int j)
{ deleteEdge( vertices[i], vertices[j]); }

void Network::deleteEdge( Vertex *vi, Vertex *vj)
{
    // delete bends with edge vi-vj
    deleteBend(vi,vj);

    //in adj of vi, find vj
    std::vector<Vertex*>::iterator it = vi->adj.begin();
    while( it != vi->adj.end() ) {
        if( *it == vj ){
            vi->adj.erase(it);
            break;
        }
        ++it;
    }

    // in adj of vj, find vi
    it = vj->adj.begin();
    while( it != vi->adj.end() ) {
        if( *it == vi ) {
            vj->adj.erase(it);
            break;
        }
        ++it;
    }

}

void Network::deleteBend(Bend *b)
{
    if( b == nullptr ) return;

    // remove pointers to this bend
    depolymerize(b);

    // exchange b with the last bend
    // in the bend list of b->mid
    std::vector<Bend*>::iterator it_b = std::find( b->mid->bends.begin(), b->mid->bends.end(), b);
    if( it_b == b->mid->bends.end() ) return;

    *it_b = b->mid->bends.back();
    b->mid->bends.pop_back();
    delete b;

}


// no guaranty that it is a valid pointer
void Network::deleteBend( std::vector<Bend*>::iterator it_bend)
{
    depolymerize(*it_bend);
    *it_bend = (*it_bend)->mid->bends.back();
    (*it_bend)->mid->bends.pop_back();
    //delete *it_bend;
    ERROR
}


void Network::deleteBend(int vi, int vj)
{ deleteBend(vertices[vi], vertices[vj] ); }

void Network::deleteBend(Vertex *vi, Vertex *vj )
{
    std::vector<Bend*>::iterator it_bend = vi->bends.begin();
    while( it_bend != vi->bends.end() ) {
        if( (*it_bend)->next == vj or (*it_bend)->prev == vj ) {
            deleteBend( it_bend );
        } 
        ++it_bend;
    }
}

void Network::deleteBend(int mid, int prev, int next)
{ deleteBend(vertices[mid], vertices[prev], vertices[next]); }

void Network::deleteBend(Vertex *mid, Vertex *prev, Vertex *next)
{
    // find this bend in bends of mid
    std::vector<Bend*>::iterator it_bend = mid->bends.begin();
    while( it_bend != mid->bends.end() ) {
        if( (*it_bend)->next == next and (*it_bend)->prev == prev ) {
            deleteBend( it_bend );
            break;
        }
        ++it_bend;
    }
}


void Network::addBend(int mid, int prev, int next)
{ addBend( vertices[mid], vertices[prev], vertices[next] ); }

void Network::addBend(Vertex *mid, Vertex *prev, Vertex *next)
{ mid->bends.push_back( new Bend(mid,prev,next) ); }

void Network::polymerize( Bend *bi, Bend *bj)
{
    bi->nextBend = bj;
    bj->prevBend = bi;
}

void Network::polymerize( Vertex *vi, Vertex *vj )
{

    //find the bend with ? - vi - vj
    std::vector<Bend*>::iterator it_bi = vi->bends.begin(); 
    while( it_bi != vi->bends.end() ) {
        if( (*it_bi)->next == vj ) break;
        ++it_bi;
    }
    //if( it_bi == vi->bends.end() ) return;

    //find the bend with  vi - vj - ?
    std::vector<Bend*>::iterator it_bj = vj->bends.begin(); 
    while( it_bj != vj->bends.end() ) {
        if( (*it_bj)->prev == vi ) break;
        ++it_bj;
    }
    //if( it_bj == vj->bends.end() ) return;

    // polymerize bends
    polymerize( *it_bi, *it_bj);
    
}

void Network::polymerize( int i, int j)
{ polymerize(vertices[i], vertices[j] ); }


void Network::depolymerize(int i, int j)
{ depolymerize( vertices[i], vertices[j]); }

void Network::depolymerize(int i)
{ depolymerize(vertices[i]); }


void Network::depolymerize( Vertex *vi, Vertex *vj)
{

    std::vector<Bend*>::iterator it_bi = vi->bends.begin();
    while( it_bi != vi->bends.end() ){
        if( (*it_bi)->next == vj ){
            // ???? remore if ???
            if( (*it_bi)->nextBend != nullptr) {
                depolymerize( *it_bi, (*it_bi)->nextBend);
            }
            break;
        }
        if( (*it_bi)->prev == vj ) {
            if( (*it_bi)->prevBend != nullptr ) {
                depolymerize( (*it_bi)->prevBend, *it_bi );
            }
            break;
        }
        ++it_bi;
    }
   
}

void Network::depolymerize( Vertex *vi)
{
    std::vector<Bend*>::iterator it_b = vi->bends.begin();
    while( it_b != vi->bends.end() ) {
        depolymerize( *it_b);
        ++it_b; 
    }
}

// remove polymer connection between bend bi and bj
void Network::depolymerize( Bend *bi, Bend *bj)
{
    bi->nextBend = nullptr;
    bj->prevBend = nullptr;
}

void Network::depolymerize( Bend *b) 
{
    if( b->nextBend != nullptr ) depolymerize(b, b->nextBend);
    if( b->prevBend != nullptr ) depolymerize(b, b->prevBend);
}


std::vector<std::vector<int> > Network::getEdges() const
{
    // list with all edges
    std::vector<std::vector<int> > edges;

    // a single edge to be added to edges
    std::vector<int> edge(2);

    std::vector<Vertex*>::const_iterator it_edge;
    std::vector<Vertex*>::const_iterator it_vertex = vertices.begin();
    while( it_vertex != vertices.end() ) { 
        int index_from = (*it_vertex)->index;
        it_edge = (*it_vertex)->adj.begin(); 
        // loop over adjacency list
        while( it_edge != (*it_vertex)->adj.end() ) {
            int index_to = (*it_edge)->index;
            // add all edges only once
            if( index_to > index_from) {
                edge[0] = index_from; 
                edge[1] = index_to; 
                edges.push_back(edge);             
            }
            ++it_edge;
        }
        ++it_vertex;
    }


    return edges;
}

std::vector<std::vector<int> > Network::getBends() const
{
    std::vector<std::vector<int> > bends;
    std::vector<int> bend(3);
    std::vector<Vertex*>::const_iterator it_v = vertices.begin();
    std::vector<Bend*>::const_iterator it_b;
    while( it_v != vertices.end() ) {
        it_b = (*it_v)->bends.begin();
        while( it_b != (*it_v)->bends.end() ) {
            bend[0] = (*it_b)->prev->index; 
            bend[1] = (*it_b)->mid->index; 
            bend[2] = (*it_b)->next->index; 
            bends.push_back( bend );
            ++it_b;
        }
        ++it_v;
    }
    return bends;
}
void Network::showAdj() const
{
    std::vector<Vertex*>::const_iterator it_v = vertices.begin();
    std::vector<Vertex*>::const_iterator it_adj;
    while( it_v != vertices.end() ) {
        std::cout << (*it_v)->index << ":";
        it_adj = (*it_v)->adj.begin();
        while( it_adj != (*it_v)->adj.end() ){
            std::cout << '\t' << (*it_adj)->index;
            ++it_adj;
        }
        std::cout << '\n';
        ++it_v;
    }
}

void Network::showBends() const
{

    std::vector<Vertex*>::const_iterator it_v = vertices.begin();
    std::vector<Bend*>::const_iterator it_b;
    while( it_v != vertices.end() ) {
        std::cout << (*it_v)->index << ":\n";
        it_b = (*it_v)->bends.begin();
        while( it_b != (*it_v)->bends.end() ) {
            std::cout << '\t' << (*it_b)->prev->index;
            std::cout << '\t' << (*it_b)->mid->index;
            std::cout << '\t' << (*it_b)->next->index;
            std::cout << '\n';
            ++it_b;
        }
        ++it_v;
    }

}

void Network::showPolymer(int i, int j) const
{
    std::cout<< i << '\t';
   // find the first bend object  
    std::vector<Bend*>::iterator it_b = vertices[i]->bends.begin();
    Bend *firstBend;
    while( it_b != vertices[i]->bends.end() ){
        if( (*it_b)->next == vertices[j] ){
            firstBend = *it_b;
            break;
        }
        ++it_b;
    }
    Bend *nextBend = firstBend->nextBend;
    while( nextBend != nullptr) {
        std::cout<< nextBend->mid->index << '\t';
        if( nextBend == firstBend) break;
        nextBend = nextBend->nextBend;
    }
}

void Network::write( std::ofstream& out)
{

    std::vector<std::vector<int> > edges = getEdges();
    std::vector<std::vector<int> > bends = getBends();
    out << vertices.size() << '\n';
    out << edges.size() << '\n';
    out << bends.size() << '\n';
    

    for( std::vector<Vertex*>::size_type vi = 0; vi< vertices.size() ; ++ vi) {
        out << vi << '\t' << vertices[vi]->r.x << '\t' << vertices[vi]->r.y << '\n';
    }
   
    for( std::vector<std::vector<int> >::size_type ei =0; ei<edges.size(); ++ei){
        out << edges[ei][0] << '\t' << edges[ei][1] << '\n';
    }

    for(std::vector<std::vector<int> >::size_type bi=0; bi<bends.size(); ++bi) {
        out << bends[bi][0] << '\t' << bends[bi][1] << '\t' << bends[bi][2] << '\n';
    }

    

}

#endif
