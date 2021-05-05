#ifndef GUARD_GRAPH_H
#define GUARD_GRAPH_H

/*
    copy constructor ...
*/

#include "vec2.h"
#include <vector>

#include <iostream>

class Graph {
public:
    Graph();
    ~Graph();

private:
    class Node;
    class Bond;
    class Bend;


    class Node {
      public:
        Node(): i(-1) {};
        Node(int index): i(index) {};
        Node(int index, double x, double y): i(index), r(x,y) {}
        ~Node(); // destroy all objects that rely on this Node (all bonds and bends)

        std::vector<Bond*>::size_type Nbonds() { return bonds.size(); }
        std::vector<Bond*>::size_type Nbends() { return bends.size(); }

        int i;
        vec2 r;

        // all bond that ara going out of this node
        std::vector<Bond*> bonds;
        // all bends of which this is the mid point
        std::vector<Bend*> bends;
    };


   
  public:
    // Member functions
    int addNode();  // add a node
    int addNode(double x, double y); // ad a node with position (x,y)
    void delNode(int i);

    bool addBond(int i, int j); // add a bond between node i and j
    void delBond(int i, int j);
    bool checkBond(int i,int j);
    

    bool addBend(int mid, int i, int j); // add a bend i-mid-j, only if bonds already exist
    void delBend(int mid, int i, int j);
    bool checkBend(int mid, int i, int j);

    // nodes i and j are next to each other on a polymer, only if bends have bonds between these nodes
    bool polymerize(int i, int j);
    // node i and j are no longer part of the same polymer
    void depolymerize(int i, int j);

    // true if there is a bond between node i and j
    bool isBonded(int i,int j);

    // set all polymer_index of the bends to the same value if they belong to the same polymer
    void setPolymers();
    // returns a vector with the length of each polymer, .size() gives the number of polymers.
    std::vector<int> getPolymerLength();

    void showBonds();
    void showBends();
    std::vector<Node*>::size_type Nnodes() const {return nodes.size(); }

  private:


    bool addBond(Node *ni, Node *nj);
    void delBond(Node *ni, Node *nj);
    bool checkBond( Node *ni, Node *nj);

    bool addBend(Node *mid, Node *a, Node *b);
    void addBend(Node *mid, Bond *a, Bond *b);
    void delBend(Node *mid, Node  *a, Node *b);
    bool checkBend(Node *mid, Node *a, Node *b);

    bool polymerize(Node *ni_ptr, Node *nj_ptr);
    bool polymerize(Bend *bi_ptr, Bend *bj_ptr);
    void depolymerize(Node *ni_ptr, Node *nj_ptr);
    void depolymerize(Bend *bi_ptr, Bend *bj_ptr);


    bool isBonded(Node *ni, Node *nj);

    std::vector<Node*> nodes;   // all nodes
    std::vector<Bend*> polymers; // contains one bend from each polymer

};


Graph::Graph()  {}

Graph::~Graph()
{
    for(std::vector<Node*>::size_type i = 0; i < nodes.size(); ++i ) {
        delete nodes[i];
    }
}

class Graph::Bond {
public:
    Bond(): from_ptr(nullptr), to_ptr(nullptr) {}
    Bond(Node *from_ptr, Node *to_ptr)
        : from_ptr(from_ptr), to_ptr(to_ptr) {}


    // this is a bond between node from_ptr to node to_ptr
    // this bond is in the vector from_ptr->bonds
    Node *from_ptr, *to_ptr;

    // info about boundary crossing
    int xb, yb;
    // rest length
    double l0;



};

class Graph::Bend {
public:
    Bend(): mid(nullptr), bond_a(nullptr), bond_b(nullptr),
        prevBend(nullptr), nextBend(nullptr), index(-1) {}
    Bend(Node *mid, Bond *a, Bond *b, int index);

    void free();

    // returns the first bend of the polymer that this is part of
    Bend* findFirst();

    // mid points to the mid point node
    Node *mid;
    // the two bonds associated with this bend
    Bond *bond_a, *bond_b;
    // previous and next bend of the polymer this bend is part of
    Bend *prevBend, *nextBend;

    int index; // index in min->bends
    int polymer_index; // index of the polymer
    double theta; // rest angle
};

// Member functions
int Graph::addNode()
{
    std::vector<Node*>::size_type Nv = nodes.size();
    nodes.push_back( new Node(Nv) );
    return Nv;
}


int Graph::addNode(double x, double y)
{
    std::vector<Node*>::size_type Nv = nodes.size();
    nodes.push_back( new Node(Nv, x, y) );
    return Nv;
}

void Graph::delNode(int i)
{
    Node *temp = nodes[i];
    nodes[i] = nodes.back();
    nodes.pop_back();
    delete temp;
}

bool Graph::addBond(int ni, int nj)
{ return addBond(nodes[ni], nodes[nj]); }

bool Graph::addBond(Node *ni_ptr, Node *nj_ptr)
{
    // if exists, do nothing
    if( checkBond(ni_ptr, nj_ptr) == true ) return false;

    ni_ptr->bonds.push_back(new Bond(ni_ptr, nj_ptr) );
    nj_ptr->bonds.push_back(new Bond(nj_ptr, ni_ptr) );
    return true;
}

bool Graph::checkBond(int i, int j)
{ return checkBond( nodes[i], nodes[j]); }

bool Graph::checkBond(Node *ni_ptr, Node *nj_ptr)
{
    // Only check if bond ni to nj exists because they always come in pairs
    std::vector<Bond*>::iterator bond_i_iter = ni_ptr->bonds.begin();
    while( bond_i_iter != ni_ptr->bonds.end() ) {
       if( (*bond_i_iter)->to_ptr == nj_ptr) return true; 
        ++bond_i_iter;
    }
    
    return false;
}

void Graph::delBond(int ni, int nj)
{ delBond(nodes[ni], nodes[nj]); }

void Graph::delBond(Node *ni_ptr, Node *nj_ptr)
{


    // check for bends that depend on this bond
    std::vector<Bend*>::iterator bend_i_iter = ni_ptr->bends.begin();
    while( bend_i_iter != ni_ptr->bends.end() ) {
        if( (*bend_i_iter)->bond_a->to_ptr == nj_ptr or
            (*bend_i_iter)->bond_b->to_ptr == nj_ptr ) {

            
            Bend *temp = *bend_i_iter;    
            // remove bend from polymer
            temp->free();
            // remove bend from list with bends, and destroy
            *bend_i_iter = ni_ptr->bends.back();
            ni_ptr->bends.pop_back();
            delete temp;

        } else {
            ++bend_i_iter;
        }
    }

    std::vector<Bend*>::iterator bend_j_iter = nj_ptr->bends.begin();
    while( bend_j_iter != nj_ptr->bends.end() ) {
        if( (*bend_j_iter)->bond_a->to_ptr == ni_ptr or
            (*bend_j_iter)->bond_b->to_ptr == ni_ptr ) {

            
            Bend *temp = *bend_j_iter;    
            // remove bend from polymer
            temp->free();
            // remove bend from list with bends, and destroy
            *bend_j_iter = nj_ptr->bends.back();
            nj_ptr->bends.pop_back();
            delete temp;

        } else {
            ++bend_j_iter;
        }
    }


    // find the bond object in ni->ptr->bonds
    std::vector<Bond*>::iterator bond_i_iter = ni_ptr->bonds.begin();
    while( bond_i_iter != ni_ptr->bonds.end() ) {
        if( (*bond_i_iter)->to_ptr == nj_ptr ) {
            Bond *temp = *bond_i_iter;
            *bond_i_iter = ni_ptr->bonds.back();
            ni_ptr->bonds.pop_back();
            delete temp;
            break;
        }
        ++bond_i_iter;
    }

    // find the bond object in ni->ptr->bonds and delete it
    std::vector<Bond*>::iterator bond_j_iter = nj_ptr->bonds.begin();
    while( bond_j_iter != nj_ptr->bonds.end() ) {
        if( (*bond_j_iter)->to_ptr == ni_ptr ) {
            Bond *temp = *bond_j_iter;
            *bond_j_iter = nj_ptr->bonds.back();
            nj_ptr->bonds.pop_back();
            delete temp;
            break;
        }
        ++bond_j_iter;
    }


}


bool Graph::polymerize(int i, int j)
{ return polymerize(nodes[i], nodes[j]); }

bool Graph::polymerize(Node *ni_ptr, Node *nj_ptr)
{

    // find the bend in ni->bends that connects ni to nj
    std::vector<Bend*>::iterator b_ni_iter = ni_ptr->bends.begin();
    while(b_ni_iter != ni_ptr->bends.end() ){
        if( (*b_ni_iter)->bond_a->to_ptr == nj_ptr 
            or (*b_ni_iter)->bond_b->to_ptr == nj_ptr ) break;
        ++b_ni_iter;
    }
    //if it is not found, return false
    if( b_ni_iter == ni_ptr->bends.end() ) return false;

    // find the bend in nj->bends that connects nj to ni    
    std::vector<Bend*>::iterator b_nj_iter = nj_ptr->bends.begin();
    while( b_nj_iter != nj_ptr->bends.end() ) {
        if( (*b_nj_iter)->bond_a->to_ptr == ni_ptr 
            or (*b_nj_iter)->bond_b->to_ptr == ni_ptr ) break;
        ++b_nj_iter;
    }

    // if it is not found, return false
    if( b_nj_iter == nj_ptr->bends.end() ) return false;

   return polymerize(*b_ni_iter, *b_nj_iter);
}

bool Graph::polymerize(Bend *bi_ptr, Bend *bj_ptr)
{
    if( bi_ptr->nextBend != nullptr ) return false; 
    if( bj_ptr->prevBend != nullptr ) return false; 

    bi_ptr->nextBend = bj_ptr;
    bj_ptr->prevBend = bi_ptr;

    return true;
}


void Graph::depolymerize(int ni, int nj)
{ depolymerize(nodes[ni], nodes[nj]); }

void Graph::depolymerize(Node *ni_ptr, Node *nj_ptr)
{
    // find the two bend objects
    Bend *bi_ptr = nullptr;
    for(std::vector<Bend*>::size_type bi = 0; bi < ni_ptr->bends.size(); ++bi ){
        if( ni_ptr->bends[bi]->bond_a->to_ptr == nj_ptr or
            ni_ptr->bends[bi]->bond_b->to_ptr == nj_ptr) {
        
            bi_ptr = ni_ptr->bends[bi];
        }
    }
    Bend *bj_ptr = nullptr;
    for(std::vector<Bend*>::size_type bj = 0; bj < nj_ptr->bends.size(); ++bj ){
        if( nj_ptr->bends[bj]->bond_a->to_ptr == ni_ptr or
            nj_ptr->bends[bj]->bond_b->to_ptr == ni_ptr) {
        
            bj_ptr = nj_ptr->bends[bj];
        }
    }
    if( bi_ptr == nullptr or bj_ptr == nullptr) return;

    depolymerize( bi_ptr, bj_ptr);
}

void Graph::depolymerize(Bend *bi_ptr, Bend *bj_ptr)
{
    if( bi_ptr == nullptr or bj_ptr == nullptr) return;

    if( bi_ptr->nextBend == bj_ptr and bj_ptr->prevBend == bi_ptr) {
        bi_ptr->nextBend = nullptr;
        bj_ptr->prevBend = nullptr;
    } else if(  bi_ptr->prevBend == bj_ptr and bj_ptr->nextBend == bi_ptr) {
        bi_ptr->prevBend = nullptr;
        bj_ptr->nextBend = nullptr;
    }
}

bool Graph::isBonded(int ni, int nj)
{ return isBonded( nodes[ni], nodes[nj]); }


bool Graph::isBonded(Node *ni_ptr, Node *nj_ptr)
{
    
    std::vector<Bond*>::iterator bond_iter = ni_ptr->bonds.begin();
    while( bond_iter != ni_ptr->bonds.end() ) {
        if( (*bond_iter)->to_ptr == nj_ptr) return true;
        ++bond_iter;
    }
    return false;

}

bool Graph::addBend(int mid, int i, int j)
{ return addBend(nodes[mid], nodes[i], nodes[j]); }

bool Graph::addBend( Node *mid, Node *a, Node *b)
{
    if( checkBend(mid,a,b) ) return false;

    Bond *bond_a = nullptr;
    Bond *bond_b = nullptr;

    std::vector<Bond*>::iterator b_iter = mid->bonds.begin();
    while( b_iter != mid->bonds.end() ) {
        if( (*b_iter)->to_ptr == a) {
            bond_a = *b_iter; 
            break;
        }
        ++b_iter;
    }

    // look for bond mid->b
    b_iter = mid->bonds.begin();
    while( b_iter != mid->bonds.end() ) {
        if( (*b_iter)->to_ptr == b) {
            bond_b = *b_iter; 
            break;
        }
        ++b_iter;
    }

    if( (bond_a != nullptr) and (bond_b != nullptr ) ) {
        addBend(mid,bond_a, bond_b);
        return  true;
    } else {
        return false;
    }

}

void Graph::addBend(Node *mid, Bond *ba, Bond *bb)
{
    std::vector<Bend*>::size_type Nb = mid->bends.size();
    mid->bends.push_back( new Bend(mid,ba,bb, Nb) );
}

bool Graph::checkBend(int mid, int a, int b)
{ return checkBend(nodes[mid], nodes[a], nodes[b]); }

bool Graph::checkBend(Node *mid_ptr, Node *ni_ptr, Node *nj_ptr )
{
    std::vector<Bend*>::iterator bend_iter = mid_ptr->bends.begin();
    while( bend_iter != mid_ptr->bends.end() ) {
        if( ( (*bend_iter)->bond_a->to_ptr == ni_ptr and (*bend_iter)->bond_b->to_ptr == nj_ptr ) or
            ( (*bend_iter)->bond_a->to_ptr == nj_ptr and (*bend_iter)->bond_b->to_ptr == ni_ptr ) )
        {
            return true;
        }
        ++bend_iter;
    }
    return false;
}

void Graph::delBend(int mid, int i, int j)
{ delBend(nodes[mid], nodes[i], nodes[j]); }

void Graph::delBend(Node *mid_ptr, Node *ni_ptr, Node *nj_ptr)
{
    // finde the pointer to this bend objedt 
    std::vector<Bend*>::iterator bend_iter = mid_ptr->bends.begin();
    while( bend_iter != mid_ptr->bends.end() ){
        if( (*bend_iter)->bond_a->to_ptr == ni_ptr and 
            (*bend_iter)->bond_b->to_ptr == nj_ptr ) {

            Bend* temp = *bend_iter;
            *bend_iter = mid_ptr->bends.back();
            mid_ptr->bends.pop_back();
            temp->free(); 
            delete temp;
            break;
        }
        if( (*bend_iter)->bond_a->to_ptr == nj_ptr and 
            (*bend_iter)->bond_b->to_ptr == ni_ptr ) {

            Bend* temp = *bend_iter;
            *bend_iter = mid_ptr->bends.back();
            mid_ptr->bends.pop_back();
            temp->free();
            delete temp;
            break;
        }

    }

}

void Graph::showBonds()
{

    std::vector<Node*>::iterator node_iter = nodes.begin();
    while( node_iter != nodes.end() ) {
        //std::cout << (*node_iter)->i << "\n";
        std::vector<Bond*>::iterator bond_iter = (*node_iter)->bonds.begin();
        while( bond_iter != (*node_iter)->bonds.end() ) {
            std::cout << (*bond_iter)->from_ptr->i << "\t" 
                    << (*bond_iter)->to_ptr->i << "\n";
            ++bond_iter;
        }
        ++node_iter;
    }
}
void Graph::showBends()
{
    setPolymers();

    for(std::vector<Node*>::size_type ni = 0; ni<nodes.size(); ++ni ){
        std::cout <<"Node: " <<  ni << "\t";
        std::cout <<"Total bends: " <<  nodes[ni]->bends.size() << std::endl;
        for(std::vector<Bend*>::size_type bi=0; bi < nodes[ni]->bends.size(); ++bi) {
            std::cout << "pi\t from \t mid \t to \n";
            std::cout << nodes[ni]->bends[bi]->polymer_index << "\t";
            std::cout << nodes[ni]->bends[bi]->bond_a->to_ptr->i << "\t";
            std::cout << nodes[ni]->i << "\t";
            std::cout << nodes[ni]->bends[bi]->bond_b->to_ptr->i << "\n";
        }
    std::cout << "_________________________________\n";
    }

}
void Graph::setPolymers()
{
    // delete previous info
    polymers.clear();
    // set all polymer indices of bends to 0
    for( std::vector<Node*>::size_type ni = 0; ni < nodes.size(); ++ni) {
        for( std::vector<Bend*>::size_type bi=0; bi < nodes[ni]->bends.size(); ++bi) {
            nodes[ni]->bends[bi]->polymer_index = -1;
        }
    }

    // for each node go through all bends
    // for each bend that has no valid index (>=0), 
    // find the first bend by using bendPrev, (if loop pick any bend)
    // set each index of the polymer
    // add pointer to first bend to polymers

    //loop over each node
     for( std::vector<Node*>::size_type ni = 0; ni < nodes.size(); ++ni) {
        // loop over each bend of node ni
        for( std::vector<Bend*>::size_type bi=0; bi < nodes[ni]->bends.size(); ++bi) {

            // if the polymer index is not set
            if( nodes[ni]->bends[bi]->polymer_index == -1) {
                // find the first bend (if loop first is this bend)
                Bend *first = nodes[ni]->bends[bi]->findFirst();
                // add first bend of polymer to polymers
                polymers.push_back(first);
                int index = polymers.size() - 1;
                first->polymer_index = index;
                // walk foreward over the polymer and set all the polymer_index
                while( first->nextBend != nullptr and first->nextBend != nodes[ni]->bends[bi] ) {
                    first = first->nextBend;
                    first->polymer_index = index;
                }
            } 
        }
    }
}

std::vector<int> Graph::getPolymerLength()
{
    setPolymers();

    std::vector<int> length(polymers.size(), 0);
    Bend *bend;
    for( std::vector<Bend*>::size_type pi = 0; pi < polymers.size(); ++pi ){
        bend = polymers[pi];
        do {
            ++length[pi];
            bend = bend->nextBend;
        }while( bend != nullptr and bend != polymers[pi] );
    }
    return length;

}

// Node
Graph::Node::~Node()
{
    // detach the bend objects of this node from the nodes, and destroy them
    for(std::vector<Bend*>::size_type i = 0; i < bends.size(); ++i ) {
        bends[i]->free();
        delete bends[i];
    }

    //  destroy all bond objects of this node
    for(std::vector<Bond*>::size_type i = 0; i < bonds.size(); ++i ) {
        // also remove the other from bond[i]->to_ptr->bonds
        std::vector<Bond*>::iterator bond_iter = bonds[i]->to_ptr->bonds.begin();
        while( bond_iter != bonds[i]->to_ptr->bonds.end() ) {
            if( (*bond_iter)->to_ptr == this ) {
                Bond *temp = *bond_iter;
                *bond_iter = bonds[i]->to_ptr->bonds.back();
                bonds[i]->to_ptr->bonds.pop_back();
                delete temp;
                break;
            }
            ++bond_iter;
        }
        delete bonds[i];
    }

}



// Bend
Graph::Bend::Bend(Node *m, Bond *ba, Bond *bb, int i)
{
    index = i;
    nextBend = nullptr;
    prevBend = nullptr;
    mid = m;
    bond_a = ba;
    bond_b = bb;

}


void Graph::Bend::free()
{

    // remove from polymer
    if( prevBend != nullptr ) {
        prevBend->nextBend = nullptr;
    }
    if( nextBend != nullptr ) {
        nextBend->prevBend = nullptr;
    }
}

Graph::Bend *Graph::Bend::findFirst()
{
    Bend *first = this;
    while( first->prevBend != nullptr ){
        first = first->prevBend;
        // check if polymer is not a loop, if so, stop
        if(first == this ) break;
    }

    return first;
}

#endif
