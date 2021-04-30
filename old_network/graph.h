#ifndef GUARD_GRAPH_H
#define GUARD_GRAPH_H


/*
    to do:
            -- throw error when adding bend without edge
            -- throw error when polymerizing nodes without a bend
            -- throw other errors

			-- copy constructor
*/

#include "vec2.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>

class Graph {
  private:
    class Vertex;
    class Edge;
    class Bend;

  public:

    Graph();
    Graph(int Nv);
    Graph ( std::ifstream& in);
    ~Graph();

    void write( std::ofstream& out);

    std::vector<Vertex*>::size_type Nvertices() const { return vertices.size(); };

    int addVertex();
    int addVertex(double x, double y);

    void setVertexPosition(int i, double x, double y);
    vec2 getVertexPosition(int i) const { return vertices[i]->r; }

    void addEdge(int i, int j, int xb=0, int yb=0, double l0=0);
    void deleteEdge(int i, int j);

    void prune(int i);

    void addBend(int mid, int va, int vb);
    void deleteBend(int mid, int va);
    void deleteAllBend(int vi, int vj); // delete all Bends with edge vi-vj

    void polymerize(int i, int j);
    void depolymerize(int i, int j); // test 
    void depolymerize(int i); // test // remove from all polymers

    // this changes the indices!!!
    void removeUnconnectedVertices(); 

    std::vector<vec2> getPositions() const;
    std::vector<std::vector<int> > getEdges() const;
    std::vector<std::vector<int> > getBends() const;
    std::vector<std::vector<int> > getPolymers();

    double averageConnectivity() const;

    void showAdj() const;
    void showBends() const;
    void showPolymers();

  private:
    class Vertex {
      public:
        Vertex(): index(-1) {}
        Vertex(int i) : index(i) {}
        Vertex(int i, double x, double y) : index(i), r(x,y) {}
        ~Vertex();
        int index;
        std::vector<Edge*> edges;
        std::vector<Bend*> bends;

        // add location info
        vec2 r;
    };

    class Edge {
      public:

        Edge(): from(nullptr), to(nullptr) {}

        Edge( Vertex *from, Vertex *to, int xb=0, int yb=0, double l0=0)
            : from(from), to(to), xb(xb), yb(yb), l0(l0) {}

        Vertex *from, *to;
        int xb, yb;
        double l0; 
    };


    class Bend {
      public:
        Bend()
            : mid(nullptr), a(nullptr), b(nullptr),
              prevBend(nullptr), nextBend(nullptr),
              index(-1), filament(-1) {}

        Bend( Vertex *mid, Edge* a, Edge *b )
            : mid(mid), a(a), b(b),
              prevBend(nullptr), nextBend(nullptr),
              filament(-1) {}

        Bend( Vertex *mid, Edge* a, Edge *b, int index)
            : mid(mid), a(a), b(b),
              prevBend(nullptr), nextBend(nullptr),
              index(index), filament(-1) {}
        Vertex* previousVertex() const;
        Vertex* nextVertex() const;

        Vertex *mid;
        Edge *a, *b;
        Bend *prevBend, *nextBend;
        int index;
        int filament;
        double theta0;
    };

        std::vector<Vertex*> vertices;

        void deleteVertex(int i);
        void deleteVertex(Vertex *vi);
        void prune( Vertex *vi);

        void exchangeVertex(int i, int j);
        void exchangeVertex(Vertex *vi, Vertex *vj);

        void addEdge(Vertex *vi, Vertex *vj, int xb=0, int yb=0, double l0=0);
        void deleteEdge(Vertex *vi, Vertex *vj);
        void deleteEdge( std::vector<Edge*>::iterator it);
        void deleteAllEdges( Vertex *vi );

        void addBend(Vertex *mid, Vertex *prev, Vertex *next);
        void addBend(Vertex *mid, Edge *a, Edge *b);
        void deleteBend(Vertex *mid, Vertex *va);
        void deleteAllBend(Vertex *vi, Vertex *vj);
        void deleteBend(Bend *bend); 
        void deleteBend( std::vector<Bend*>::iterator it_bend);

        void polymerize(Vertex *vi, Vertex *vj);
        void polymerize(Bend *bi, Bend *bj); 

        void depolymerize( Vertex *vi, Vertex *vj); 
        void depolymerize( Vertex *vi);
        void depolymerize( Bend *bi, Bend *bj);
        void depolymerize( Bend *bi); 

        void resetFilamentIndex();
        Bend* firstBend( Bend *bend) const;
};

Graph::Graph() {};

Graph::Graph( int Nv)
: vertices( std::vector<Vertex*>(Nv, nullptr) )
{
    for( int vi=0; vi<Nv; ++vi) {
        vertices[vi] = new Vertex(vi);
    }
}


Graph::~Graph()
{ 
    for( std::vector<Vertex*>::size_type vi=0; vi < vertices.size(); ++vi) {
        delete vertices[vi];
        vertices[vi] = nullptr;
    }
}

int Graph::addVertex()
{
    std::vector<Vertex*>::size_type  Nv = vertices.size();
    vertices.push_back( new Vertex(Nv) );
    return Nv;
}

int Graph::addVertex(double x, double y)
{
    std::vector<Vertex*>::size_type  Nv = vertices.size();
    vertices.push_back( new Vertex(Nv, x, y) );
    return Nv;
}

void Graph::deleteVertex(int i)
{ deleteVertex( vertices[i] ); }

void Graph::deleteVertex(Vertex *vi)
{
    deleteAllEdges(vi); 
    exchangeVertex( vi, vertices.back() ); 
    delete vertices.back();
    vertices.pop_back();
}
void Graph::exchangeVertex(int i, int j)
{ exchangeVertex( vertices[i], vertices[j] ); }

void Graph::exchangeVertex(Vertex *vi, Vertex *vj)
{
    int i = vi->index;
    int j = vj->index;
    vertices[i] = vj;
    vertices[i]->index = i;
    vertices[j] = vi;
    vertices[j]->index = j;
}


void Graph::setVertexPosition(int i, double x, double y)
{
    vertices[i]->r.x = x;
    vertices[i]->r.y = y;
}

void Graph::polymerize( Bend *bi, Bend *bj)
{
    if( bi == nullptr or bj == nullptr ) return;
    // ERROR if bi->nextBend != nullptr or bj->prevBend != nullptr

    bi->nextBend = bj;
    bj->prevBend = bi;
}

void Graph::polymerize( Vertex *vi, Vertex *vj )
{
    // find the bend in vi that goes to vj
    std::vector<Bend*>::iterator it_bi = vi->bends.begin();
    while( it_bi != vi->bends.end() ) {
        if( (*it_bi)->a->to == vj or (*it_bi)->b->to == vj ) break;
        ++it_bi;
    }

    // do nothing if bend not found
    if( it_bi == vi->bends.end() ) return;

    std::vector<Bend*>::iterator it_bj = vj->bends.begin();
    while( it_bj != vj->bends.end() ) {
        if( (*it_bj)->a->to == vi or (*it_bj)->b->to == vi ) break;
        ++it_bj;
    }
    if( it_bj == vj->bends.end() ) return;

    polymerize(*it_bi, *it_bj);
}

void Graph::polymerize( int i, int j)
{ polymerize(vertices[i], vertices[j] ); }

void Graph::depolymerize( Bend *bi, Bend *bj)
{ 
    if( bi == nullptr or bj == nullptr ) return;
    if( bi->nextBend != bj or bj->prevBend != bi) return;
    bi->nextBend = nullptr;
    bj->prevBend = nullptr;
} 

void Graph::depolymerize( Bend *b)
{ 
    if( b == nullptr ) return;
    if( b->nextBend != nullptr ) {
        b->nextBend->prevBend = nullptr;
        b->nextBend = nullptr;
    }
    if( b->prevBend != nullptr ) {
        b->prevBend->nextBend = nullptr;
        b->prevBend = nullptr;
    }
}

void Graph::depolymerize( Vertex *vi)
{
    std::vector<Bend*>::iterator it_b = vi->bends.begin();
    while( it_b != vi->bends.end() ) {
        depolymerize(*it_b);
        ++it_b;
    }
}

void Graph::depolymerize(int i)
{ depolymerize(vertices[i]); }


void Graph::depolymerize( Vertex *vi, Vertex *vj)
{
    // find bend from vi to vj
    std::vector<Bend*>::iterator it_b = vi->bends.begin();
    while( it_b != vi->bends.end() ) {
        if( (*it_b)->a->to == vj or (*it_b)->b->to == vj ) break; 
        ++it_b;
    }

    // do nothing if bend not found
    if( it_b == vi->bends.end() ) return;

    if( (*it_b)->nextBend != nullptr and (*it_b)->nextBend->mid == vj ) {
        depolymerize( *it_b, (*it_b)->nextBend );
    } else if( (*it_b)->prevBend != nullptr and (*it_b)->prevBend->mid == vi ) {
        depolymerize( (*it_b)->prevBend, *it_b );
    }

}

void Graph::depolymerize(int i, int j)
{ depolymerize( vertices[i], vertices[j]); }



void Graph::addBend(Vertex *mid, Edge *a, Edge *b)
{
    Bend* bend = new Bend( mid,a,b,mid->bends.size() );
    mid->bends.push_back( bend );
}

void Graph::addBend(int mid, int va, int vb)
{ addBend( vertices[mid], vertices[va], vertices[vb]); }

void Graph::addBend(Vertex *mid, Vertex *va, Vertex *vb)
{

    // in edges of mid, find mid-va and mid-vb 
    std::vector<Edge*>::iterator it_e = mid->edges.begin();
    Edge *a = nullptr;
    Edge *b = nullptr;
    while( it_e != mid->edges.end() ) {
        if( (*it_e)->to == va ) a = *it_e;
        if( (*it_e)->to == vb ) b = *it_e;
        if( a != nullptr and b != nullptr ) {
            addBend(mid, a, b);
            break;
        }
        ++it_e;
    }
    // ERROR: NO EDGE ADDED
}

void Graph::deleteBend( std::vector<Bend*>::iterator it_bend)
{
    depolymerize( *it_bend);
    // move last bend pointer to current location in the bend list
    Bend *temp = *it_bend;
    *it_bend = temp->mid->bends.back();
    (*it_bend)->index = temp->index;
    temp->mid->bends.pop_back();
    delete temp;
    
}

void Graph::deleteBend(Bend *bend)
{
    if( bend == nullptr ) return;

    depolymerize(bend); 

    deleteBend( bend->mid->bends.begin() + bend->index ); 

}


void Graph::deleteBend(Vertex *mid, Vertex *va)
{
    std::vector<Bend*>::iterator it_b = mid->bends.begin();
    while( it_b != mid->bends.end() ) {
        if( (*it_b)->a->to == va or (*it_b)->b->to == va ) {
            deleteBend( it_b);
            break;
        }
        ++it_b;
    }
    
}

void Graph::deleteBend(int mid, int va)
{ deleteBend( vertices[mid], vertices[va]); }


void Graph::deleteAllBend(int vi, int vj)
{ deleteAllBend( vertices[vi], vertices[vj]); }

void Graph::deleteAllBend(Vertex *vi, Vertex *vj)
{
    deleteBend(vi,vj);
    deleteBend(vj,vi);
}

void Graph::addEdge(int i, int j, int xb, int yb, double l0)
{ addEdge( vertices[i], vertices[j], xb, yb, l0); }

void Graph::addEdge(Vertex *vi, Vertex *vj, int xb, int yb, double l0)
{
    // check if already exists??
    vi->edges.push_back( new Edge(vi,vj,xb,yb,l0) );
    vj->edges.push_back( new Edge(vj,vi,-xb,-yb,l0) );

}

void Graph::deleteEdge(int vi, int vj)
{ deleteEdge(vertices[vi], vertices[vj]); }

void Graph::deleteEdge(Vertex *vi, Vertex *vj)
{
    // find the edge from vi to vj
    std::vector<Edge*>::iterator it_e = vi->edges.begin();
    while( it_e != vi->edges.end() ) {
        if( (*it_e)->to == vj ) {
            deleteEdge(it_e);
            break;
        }
        ++it_e;
    }
    it_e = vj->edges.begin(); 
    while( it_e != vj->edges.end() ) {
        if( (*it_e)->to == vi ) {
            deleteEdge(it_e);
            break;
        }
        ++it_e;
    }
    
}

void Graph::deleteEdge( std::vector<Edge*>::iterator it)
{

    // delete the corresponding bend
    deleteAllBend( (*it)->from, (*it)->to );

    Edge *temp = *it;
    *it = temp->from->edges.back();
    temp->from->edges.pop_back();
    delete temp;
}


void Graph::deleteAllEdges( Vertex *vi )
{
    while( vi->edges.size() > 0 ) {
        deleteEdge( vi->edges.end() - 1 );
    }
}

Graph::Vertex::~Vertex()
{
    //delete edges
    for( std::vector<Edge*>::size_type ei=0;ei<edges.size(); ++ei) {
        delete edges[ei];
        edges[ei] = nullptr;
    }
    //delete bends
    for( std::vector<Bend*>::size_type bi=0;bi<bends.size(); ++bi) {
        delete bends[bi];
        bends[bi] = nullptr;
    }
}

void Graph::showAdj() const
{
    std::vector<Vertex*>::const_iterator it_v = vertices.begin();
    std::vector<Edge*>::const_iterator it_e;
    while( it_v != vertices.end() ) {
        std::cout << (*it_v)->index << ":\t";
        it_e = (*it_v)->edges.begin();
        while( it_e != (*it_v)->edges.end() ) {
            std::cout << (*it_e)->to->index << '\t';
            ++it_e;
        } 
        std::cout << '\n';
        ++it_v;
    }
}

void Graph::showBends() const
{
    std::vector<Vertex*>::const_iterator it_v = vertices.begin();
    std::vector<Bend*>::const_iterator it_b;
    while( it_v != vertices.end() ) {
        std::cout << (*it_v)->index << ":\n";
        it_b = (*it_v)->bends.begin();
        while( it_b != (*it_v)->bends.end() ) {
            std::cout << '\t';
            std::cout << (*it_b)->a->to->index << '\t';
            std::cout << (*it_b)->mid->index << '\t';
            std::cout << (*it_b)->b->to->index << '\n';
            ++it_b;
        }
        ++it_v;
    }
}


std::vector<vec2> Graph::getPositions() const
{
    std::vector<vec2> positions(vertices.size() );
    for(unsigned long int vi=0;vi<vertices.size(); ++vi) {
        positions[vi] = vertices[vi]->r;
    }

    return positions;
}

std::vector<std::vector<int> > Graph::getEdges() const
{
    std::vector<std::vector<int> > edges;
    std::vector<int> edge(5);

    std::vector<Vertex*>::const_iterator it_v = vertices.begin();
    std::vector<Edge*>::const_iterator it_e;
    while( it_v != vertices.end() ) {
        it_e = (*it_v)->edges.begin(); 
        while( it_e != (*it_v)->edges.end() ) {
            if( (*it_e)->to->index > (*it_e)->from->index ) {
                edge[0] = (*it_e)->from->index; 
                edge[1] = (*it_e)->to->index; 
                edge[2] = (*it_e)->xb;
                edge[3] = (*it_e)->yb;
                edge[4] = (*it_e)->l0;
                edges.push_back(edge);
            }
            ++it_e;
        }
        ++it_v;
    }
    return edges;
}

std::vector<std::vector<int> > Graph::getBends() const
{

    std::vector<std::vector<int> > bends;
    std::vector<int> bend(8);

    std::vector<Vertex*>::const_iterator it_v = vertices.begin();
    std::vector<Bend*>::const_iterator it_b;
    while( it_v != vertices.end() ) {
        it_b = (*it_v)->bends.begin();
        while( it_b != (*it_v)->bends.end() ) {
            bend[0] = (*it_b)->mid->index;
            bend[1] = (*it_b)->a->to->index;
            bend[2] = (*it_b)->a->xb;
            bend[3] = (*it_b)->a->yb;
            bend[4] = (*it_b)->b->to->index;
            bend[5] = (*it_b)->b->xb;
            bend[6] = (*it_b)->b->yb;
            bend[7] = (*it_b)->theta0;
            bends.push_back(bend); 
            
            ++it_b;
        }
        ++it_v;
    }
    return bends; 
}

double Graph::averageConnectivity() const
{
    int Nv = 0; // number of vertices with at least one bond
    double avgc = 0;
    
    std::vector<Vertex*>::const_iterator it_v = vertices.begin();
    while( it_v != vertices.end() ) {
        avgc += (*it_v)->edges.size();
        if( (*it_v)->edges.size() > 0 ) Nv+=1;
        ++it_v;
    }

    return avgc/Nv;

}


Graph::Bend* Graph::firstBend( Bend *bend) const
{
    if( bend == nullptr ) return nullptr;
    Bend *first = bend;
    while( bend->prevBend != nullptr) {
        bend = bend->prevBend;
        if( bend == first ) break;
    }
    return bend;
}

std::vector<std::vector<int> > Graph::getPolymers() 
{
    std::vector<std::vector<int> > polymers;
    polymers.push_back( std::vector<int>() );

    resetFilamentIndex();

    int filament_index = 0;
    std::vector<Vertex*>::const_iterator it_v = vertices.begin();
    std::vector<Bend*>::const_iterator it_b;
    while( it_v != vertices.end() ) {
        it_b = (*it_v)->bends.begin(); 
        while( it_b != (*it_v)->bends.end() ) {

            // if the bend is part of a not yet added filament
            if( (*it_b)->filament == -1 ){
                Bend *bend = *it_b;            
                Bend *first = firstBend(bend);
                bend = first;
                while( bend != nullptr) {
                    // set filament index
                    bend->filament = filament_index;


                    // include the first vertex
                    if( bend == first ) {
                        polymers[filament_index].push_back( bend->previousVertex()->index );
                    }

                    // add vertex index
                    polymers[filament_index].push_back(bend->mid->index);

                    // if it is the last bend, also include the next vertex index
                    if( bend->nextBend == nullptr ) {
                        polymers[filament_index].push_back( bend->nextVertex()->index );
                    }

                    if(  bend->nextBend == first ) {
                        polymers[filament_index].push_back( bend->nextVertex()->index );
                        polymers[filament_index].push_back( first->nextVertex()->index );
                        break;
                    }

                    bend = bend->nextBend;
                }
                polymers.push_back( std::vector<int>() );
                ++filament_index;
            }
            ++it_b;
        }
        ++it_v;
    }
    polymers.pop_back();
    return polymers;
}

void Graph::resetFilamentIndex()
{
    std::vector<Vertex*>::const_iterator it_v = vertices.begin();
    std::vector<Bend*>::const_iterator it_b;
    while( it_v != vertices.end() ) {
        it_b = (*it_v)->bends.begin(); 
        while( it_b != (*it_v)->bends.end() ) {
            (*it_b)->filament = -1;
            ++it_b;
        }
        ++it_v;
    }   
}


void Graph::showPolymers()
{
    std::vector<std::vector<int> > p = getPolymers();
    std::vector<std::vector<int> >::iterator it_p = p.begin();
    while( it_p != p.end() ) {
        std::vector<int>::iterator it = it_p->begin();
        while( it != it_p->end() ) {
            std::cout << *it << '\t';
            ++it;
        }
        std::cout << '\n';
        ++it_p;
    }

}


Graph::Vertex* Graph::Bend::previousVertex() const
{
    if( prevBend != nullptr ) return prevBend->mid;

    // if bend is not part of filament, then a->to is called the next vertex
    if( nextBend == nullptr ) return a->to;
    if( a->to == nextBend->mid ) return b->to;
    return a->to;
}

Graph::Vertex* Graph::Bend::nextVertex() const
{
    if( nextBend != nullptr ) return nextBend->mid;

    // if bend is not part of filament, then b->to is called the vertex
    if( prevBend == nullptr ) return b->to;

    if( a->to == prevBend->mid ) return b->to;
    return a->to;
}

void Graph::write( std::ofstream& out) 
{
    std::vector<std::vector<int> > edges = getEdges();
    std::vector<std::vector<int> > bends = getBends();
    std::vector<std::vector<int> > polymers = getPolymers();
    int Nv = vertices.size();
    int Ne = edges.size();
    int Nb = bends.size();
    int Np = polymers.size();

    out << Nv << '\n';
    out << Ne << '\n';
    out << Nb << '\n';
    out << Np << '\n';

    std::vector<Vertex*>::iterator it_v = vertices.begin();
    while( it_v != vertices.end() ) {
        out << (*it_v)->index << '\t';
        out << (*it_v)->r.x << '\t';
        out << (*it_v)->r.y << '\n';
        ++it_v;
    }

    std::vector<std::vector<int> >::iterator it_e = edges.begin();
    while( it_e != edges.end() ) {
        out << (*it_e)[0] << '\t';
        out << (*it_e)[1] << '\t';
        out << (*it_e)[2] << '\t';
        out << (*it_e)[3] << '\t';
        out << (*it_e)[4] << '\n';
        ++it_e;
    }

    std::vector<std::vector<int> >::iterator it_b = bends.begin();
    while( it_b != bends.end() ) {
        out << (*it_b)[0] << '\t';
        out << (*it_b)[1] << '\t';
        out << (*it_b)[2] << '\t';
        out << (*it_b)[3] << '\t';
        out << (*it_b)[4] << '\t';
        out << (*it_b)[5] << '\t';
        out << (*it_b)[6] << '\t';
        out << (*it_b)[7] << '\n';
        ++it_b;
    }

    std::vector<std::vector<int> >::iterator it_p = polymers.begin();
    std::vector<int>::iterator it_pi;
    while( it_p != polymers.end() ) {
        it_pi = it_p->begin();
        while( it_pi != it_p->end() ) {
            out << *it_pi;
            if( it_pi == it_p->end() - 1) out << '\n';
            else out << '\t';
            ++it_pi;
        }    
        ++it_p;
    }
 
}

Graph::Graph( std::ifstream& in) 
{
    int Nv, Ne, Nb, Np; 
    in >> Nv;
    in >> Ne;
    in >> Nb;
    in >> Np;
    std::vector<std::vector<int> > edges(Ne, std::vector<int>(5) );
    std::vector<std::vector<int> > bends(Nb, std::vector<int>(8) );
    std::vector<std::vector<int> > polymers(Np);

    for( int li = 0; li < Ne; ++li ) {
        for( int i =0;i<5; ++i) in >> edges[li][i];    
    } 

    for( int li=0; li < Nb; ++li ) {
        for(int i=0;i<8;++i) in >> bends[li][i];
    }

   std::string line; 
    int pi;
    std::getline( in, line);
    for( int li=0; li < Np; ++li ){
        std::getline( in, line);
        std::stringstream ss(line);
        while( ss >> pi ) polymers[li].push_back(pi);
    } 

    for(int vi=0;vi<Nv;++vi) addVertex();

    for(int ei=0;ei<Ne;++ei) {
        addEdge( edges[ei][0], edges[ei][1],edges[ei][2], edges[ei][3], edges[ei][4]);
    }
    for(int bi=0; bi<Nb;++bi) {
        addBend( bends[bi][0], bends[bi][1], bends[bi][4]);
    }
   
    for(int pi=0; pi<Np;++pi) {

        // do not polymerize if polymer is single bend
        if( polymers[pi].size() < 4 ) continue; 
        int Npi = polymers[pi].size();
        // if polymer is a loop
        if( polymers[pi][0] == polymers[pi][Npi-3] ) {
            for( int bi = 1; bi < Npi-2; ++bi) {
                int a = polymers[pi][bi];
                int b = polymers[pi][bi+1];
                polymerize(a,b); 
            }
        } else {
            for( int bi = 1; bi < Npi-2; ++bi) {
                int a = polymers[pi][bi];
                int b = polymers[pi][bi+1];
                polymerize(a,b); 
            }

        }
    } 
}

void Graph::prune( int i)
{ prune( vertices[i] ); }

void Graph::prune(Vertex *vi)
{
    if( vi->edges.size() == 1 ){
        deleteEdge( vi, vi->edges[0]->to );
        prune( vi->edges[0]->to );
    }
}

void Graph::removeUnconnectedVertices()
{
    std::vector<Vertex*>::iterator it_v = vertices.begin();
    while( it_v != vertices.end() ) {
        if( (*it_v)->edges.size() == 0 ) {
            deleteVertex(*it_v);
        } else {
            ++it_v;
        }
    }


}
#endif
