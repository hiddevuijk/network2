#ifndef GUARD_GRAPH_H
#define GUARD_GRAPH_H

#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>



class Graph
{
  public:
    Graph(): Nv(0) {};
    Graph(int Nv); 
    Graph( std::ifstream& in);

    void showAdj() const;
    void showBends() const;
      
    int V() const { return Nv; }

    void set_position(int vi, double x, double y)
    { vertices[vi].x = x; vertices[vi].y = y; }

    void addEdge(int i, int j)
    {
        vertices[i].addAdj(j);
        vertices[j].addAdj(i);
    }

    bool checkEdge(int i, int j);
    void deleteEdge(int i, int j);
    void exchangeVertices(int i,int j);
    void addVertex() {
        vertices.push_back( Vertex(Nv) );
        ++Nv;
    }


    void deleteVertex(int i);
    int Nneighbours (int i) const { return vertices[i].Nneighbours(); } 
  //private:


    class Bend {
        // bend triplet prev-I-next
      public:
        Bend(): I(-1), prev(-1), next(-1), prevBend(nullptr), nextBend(nullptr) {}

        Bend(int I, int prev, int next) 
            : I(I), prev(prev), next(next), prevBend(nullptr), nextBend(nullptr) {}

        Bend( int I, int prev, int next, Bend *prevBend, Bend *nextBend)
            : I(I), prev(prev), next(next), prevBend(prevBend), nextBend(nextBend) {}


        bool operator==( const Bend& rhs) { return this == &rhs; }
        int I;
        int prev, next;
        Bend *prevBend;
        Bend *nextBend;
        double theta0;
    };

    void addBend(int vi, int viPrev, int viNext, Graph::Bend *prevBend, Graph::Bend *nextBend);
    void addBend(int vi, int viPrev, int viNext);
    
    void deleteBend( Bend& bend);
    void polymerize( Graph::Bend& b1, Graph::Bend& b2);
    void polymerize( int i, int j );

    class Vertex {
      public:
        Vertex(): index(-1) {} 
        Vertex(int index): index(index) {} 

        int index;
        double x,y;



        void set_index(int i) { index = i; }
        void addAdj(int vi);
        int Nneighbours() const { return adj.size(); }

        std::vector<int> adj;
        std::vector<Bend> bends;
    };
  
    Graph::Bend* nextBend( Graph::Bend bend); // return next bend
    Graph::Bend* prevBend( Graph::Bend bend); // return previous bend

    int Nv; // number ov vertices
    std::vector<Vertex> vertices;
};


Graph::Bend* Graph::nextBend( Graph::Bend bend)
{ return bend.nextBend; }

Graph::Bend* Graph::prevBend( Graph::Bend bend)
{ return bend.prevBend; }


void Graph::Vertex::addAdj(int vi)
{
    // if vp not in adj, add vp to adj
    if( std::find(adj.begin(), adj.end(), vi) == adj.end() ) {
        adj.push_back( vi );
    }
}



Graph::Graph(int Nv)
: Nv(Nv), vertices(Nv)
{
    for(int i=0; i<Nv; ++i) {
        vertices[i].index = i;
    }
}

Graph::Graph( std::ifstream& in)
{
    in >> Nv;
    vertices = std::vector<Vertex>(Nv);
    for(int i=0; i<Nv; ++i) {
        vertices[i].index = i;
    }

    int v,w;
    while( ( in >> v ) and ( in >> w) ){
        addEdge(v,w);
    }

}


bool Graph::checkEdge(int i, int j) 
{ return std::find( vertices[i].adj.begin(), vertices[i].adj.end(), j  ) != vertices[i].adj.end(); }

void Graph::deleteEdge(int i, int j)
{
    std::vector<int>::iterator it;
    std::vector<Bend>::iterator it_b;

    // remove bends
    it_b = vertices[i].bends.begin(); 
    while( it_b != vertices[i].bends.end() ) {
        if( it_b->prev == j or it_b->next == j) deleteBend( *it_b );
        if( it_b == vertices[i].bends.end() ) break;
        ++it_b;
    }
    it_b = vertices[j].bends.begin();
    while( it_b != vertices[j].bends.end() ) {
        if( it_b->prev == i or it_b->next == i ) deleteBend( *it_b);
        if( it_b == vertices[j].bends.end() ) break;
        ++it_b;
    }


    // if edge exist, remove it from adj[i] and adj[j]
    it = std::find( vertices[i].adj.begin(), vertices[i].adj.end(), j);
    if( it != vertices[i].adj.end() ) {
        vertices[i].adj.erase(it);
        it = std::find( vertices[j].adj.begin(), vertices[j].adj.end(), i);
        vertices[j].adj.erase(it);
    }

    

}

void Graph::exchangeVertices(int i, int j)
{
    if( i == j ) return;
    bool connected = checkEdge(i,j);
    
    Vertex vi = vertices[i]; 
    Vertex vj = vertices[j];

    vertices[i] = vj;
    vertices[i].index = i;
    vertices[j] = vi;
    vertices[j].index = j;

    // fix adjacencies in adj of the adjacent vertices
    std::vector<int>::iterator it = vi.adj.begin();
    while( it != vi.adj.end() ) {
        *std::find( vertices[*it].adj.begin(), vertices[*it].adj.end(), i) = j;
        ++it;
    }
    
    it = vj.adj.begin();
    while( it != vj.adj.end() ) {
        *std::find( vertices[*it].adj.begin(), vertices[*it].adj.end(), j) = i;
        ++it;
    }

    // if i and j connected they have self connections
    // change i -> i to i -> j and j->j to j->i
    if(connected) {
        *std::find( vertices[i].adj.begin(), vertices[i].adj.end(), i ) = j;
        *std::find( vertices[j].adj.begin(), vertices[j].adj.end(), j ) = i;
    }

    // fix bends
    

    // fix bends
    // loop over all edges of i, if the connected vertex has a bend with prev or next is j
    // change to i 
    std::vector<Bend>::iterator it_bend;
    it = vertices[i].adj.begin();

    while( it != vertices[i].adj.end() ){
        it_bend = vertices[*it].bends.begin();
        while( it_bend != vertices[*it].bends.end() ) {
            if( it_bend->prev == i ) it_bend->prev = j;
            if( it_bend->next == i ) it_bend->next = j;
            if( it_bend->prev == j ) it_bend->prev = i;
            if( it_bend->next == j ) it_bend->next = i;
            ++it_bend;
        } 
        ++it;
    }
    it = vertices[j].adj.begin();
    while( it != vertices[j].adj.end() ){
        it_bend = vertices[*it].bends.begin();
        // if thie vertex has already been visited, skip it
        if( !checkEdge(i,*it) ) {
            while( it_bend != vertices[*it].bends.end() ) {
                if( it_bend->prev == j ) it_bend->prev = i;
                if( it_bend->next == j ) it_bend->next = i;
                if( it_bend->prev == i ) it_bend->prev = j;
                if( it_bend->next == i ) it_bend->next = j;
                ++it_bend;
            } 
        }
        ++it;
    }

    // fix pointers to prev and next bends
    it_bend = vertices[i].bends.begin(); 
    while( it_bend != vertices[i].bends.end() ) {
        if( it_bend->nextBend != nullptr ) {
            it_bend->nextBend->prevBend->nextBend = nullptr;
            it_bend->nextBend->prevBend = &(*it_bend);
        }
        if( it_bend->prevBend != nullptr ) {
            it_bend->prevBend->nextBend->prevBend = nullptr;
            it_bend->prevBend->nextBend = &(*it_bend);
        }
        ++it_bend;
    }

    it_bend = vertices[j].bends.begin(); 
    while( it_bend != vertices[j].bends.end() ) {
        if( it_bend->nextBend != nullptr ) {
            it_bend->nextBend->prevBend = &(*it_bend);
        }
        if( it_bend->prevBend != nullptr ) {
            it_bend->prevBend->nextBend = &(*it_bend);
        }
        ++it_bend;
    }
}

void Graph::polymerize( Graph::Bend& b1, Graph::Bend& b2)
{
    b1.nextBend = &b2;
    b2.prevBend = &b1;
}

void Graph::polymerize( int i, int j)
{

    // find bend * - i - j
    std::vector<Bend>::iterator it_b1 = vertices[i].bends.begin();
    while( it_b1->next != j ){
         ++it_b1;
        if( it_b1 == vertices[i].bends.end() ) std::cout << "error 1\n";
    }

    // find bned i-j-*
    std::vector<Bend>::iterator it_b2 = vertices[j].bends.begin();
    while( it_b2->prev != i ) {
        ++it_b2;
        if( it_b2 == vertices[j].bends.end() ) std::cout << "error 2\n";
    }

    polymerize( *it_b1, *it_b2);

}

void Graph::deleteBend( Graph::Bend& bend)
{

    if( bend.prevBend != nullptr) bend.prevBend->nextBend = nullptr;
    bend.prevBend = nullptr;
    if( bend.nextBend != nullptr) bend.nextBend->prevBend = nullptr;
    bend.nextBend = nullptr;

    if( vertices[bend.I].bends.size() > 1 ) {
        std::vector<Bend>::iterator it = std::find( vertices[bend.I].bends.begin(), vertices[bend.I].bends.end(), bend );
        *it = vertices[bend.I].bends.back();
        if( it->prevBend != nullptr ) it->prevBend->nextBend = &(*it);     
        if( it->nextBend != nullptr ) it->nextBend->prevBend = &(*it);     

    }

    
    vertices[bend.I].bends.pop_back();

}


void Graph::deleteVertex(int i)
{

    // delete all bends with * - i - *
    while( vertices[i].bends.size() > 0 ) deleteBend( vertices[i].bends.back() );
  
    // delete all edge connected to vertex i
    while( vertices[i].adj.size() > 0 ) {
        int j = vertices[i].adj[0];
        deleteEdge(i,j);
    }

    int viLast = vertices.size() - 1; // index of the last vertex

    // exchange the indeces of vertex i and the last vertex
    exchangeVertices(i, viLast);   
    
    // remove last vertex
    vertices.pop_back();
    Nv -= 1;
}

void Graph::showAdj() const
{
    
    for(int iv=0; iv < Nv; ++iv) {
        std::cout << iv << " : ";
        for( std::vector<int>::size_type ie=0; ie < vertices[iv].adj.size(); ++ie ) {
            std::cout << vertices[iv].adj[ie] << " ";
        }
        std::cout << std::endl;
    }
}

void Graph::showBends() const
{
    
    for(int iv=0; iv < Nv; ++iv) {
        std::cout << iv << " :\n";
        for( std::vector<int>::size_type ib=0; ib < vertices[iv].bends.size(); ++ib ) {
            std::cout << "\t " << vertices[iv].bends[ib].prev << '\t';
            std::cout << "\t " << vertices[iv].bends[ib].next << '\n';
        }
        std::cout << std::endl;
    }
}

void Graph::addBend(int vi, int viPrev, int viNext, Graph::Bend *prevBend, Graph::Bend *nextBend) 
{ vertices[vi].bends.push_back( Graph::Bend(vi, viPrev,  viNext, prevBend, nextBend) ); }

void Graph::addBend(int vi, int viPrev, int viNext) 
{ vertices[vi].bends.push_back( Graph::Bend(vi, viPrev,  viNext) ); }

#endif
