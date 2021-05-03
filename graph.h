#ifndef GUARD_GRAPH_H
#define GUARD_GRAPH_H

/* 
    add const functions 
*/


#include "vec2.h"
#include <vector>

#include <iostream>

class Graph {
public:
    class Node;
    class Bond;
    class Bend;

    Graph();
    ~Graph();

    class Bond {
    public:
        Bond(): from_ptr(nullptr), to_ptr(nullptr) {}
        Bond(Node *from_ptr, Node *to_ptr)
            : from_ptr(from_ptr), to_ptr(to_ptr) {}

        Node *from_ptr, *to_ptr;

        // all bends that rely on this bond
        std::vector<Bend*> bends;

        int xb, yb;
        double l0;
    };
    
    class Bend {
    public:
        Bend(): mid(nullptr), bond_a(nullptr), bond_b(nullptr),
            prevBend(nullptr), nextBend(nullptr), index(-1) {}

        Bend(Node *mid, Bond *a, Bond *b, int index);

        // removes all connections of bend, still exists in bends of Node
        void remove();

        Node *mid;
        Bond *bond_a, *bond_b;
        Bend *prevBend, *nextBend;

        int index;
        int polymer;
        double theta;
    };

    class Node {
    public:
        Node(): i(-1) {};
        Node(int index): i(index) {};
        Node(int index, double x, double y): i(index), r(x,y) {}
        ~Node();

        std::vector<Bond*>::size_type Nbonds() { return bonds.size(); }
        std::vector<Bond*>::size_type Nbends() { return bends.size(); }

        int i;
        vec2 r;

        std::vector<Bond*> bonds;
        std::vector<Bend*> bends;
    };

    // Member functions
    int addNode();
    int addNode(double x, double y);
    void delNode(int i);

    void addBond(int i, int j);
    void delBond(int i, int j);

    bool addBend(int mid, int i, int j);

    // true if there is a bond between node i and j
    bool isBonded(int i,int j);
    //private

    void delNode(Node *ni);

    void addBond(Node *ni, Node *nj);
    void delBond(Node *ni, Node *nj);
    void delBond(std::vector<Bond*>::iterator);

    bool addBend(Node *mid, Node *a, Node *b);
    void addBend(Node *mid, Bond *a, Bond *b);
    void delBend(Node *mid, Node *a, Node *b);

    bool isBonded(Node *ni, Node *nj);

    void showBonds();
    std::vector<Node*>::size_type Nnodes() const {return nodes.size(); }
//private:
    std::vector<Node*> nodes;   // all nodes
    std::vector<Node*> polymers; // contains one bend from each polymer

};


Graph::Graph()  {}

Graph::~Graph()
{
    for(std::vector<Node*>::size_type i = 0; i < nodes.size(); ++i ) {
        delete nodes[i];
    }
}


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
    delete nodes[i];
    nodes[i] = nodes.back();
    nodes.pop_back();
}

void Graph::delNode(Node *ni)
{ delNode(ni->i); }

void Graph::addBond(int ni, int nj)
{ addBond(nodes[ni], nodes[nj]); }

void Graph::addBond(Node *ni_ptr, Node *nj_ptr)
{
    ni_ptr->bonds.push_back(new Bond(ni_ptr, nj_ptr) );
    nj_ptr->bonds.push_back(new Bond(nj_ptr, ni_ptr) );
}

void Graph::delBond(int ni, int nj)
{ delBond(nodes[ni], nodes[nj]); }

void Graph::delBond(Node *ni, Node *nj )
{
    // find the bond from ni to nj
    std::vector<Bond*>::iterator bond_iter = ni->bonds.begin();
    while( bond_iter != ni->bonds.end() ) {
        if( (*bond_iter)->to_ptr == nj ) {
            delBond( bond_iter);
            break;
        }
        ++bond_iter;
    }

    // find the bond from nj to ni
    bond_iter = nj->bonds.begin();
    while( bond_iter != nj->bonds.end() ) {
        if( (*bond_iter)->to_ptr == ni ) {
            delBond( bond_iter);
            break;
        }
        ++bond_iter;
    }
        

}

void Graph::delBond( std::vector<Bond*>::iterator bond_iter)
{
    // delete all Bends that rely on the bond *bond_iter points to
    while( (*bond_iter)->bends.end() != (*bond_iter)->bends.begin() ) {
        Bend *temp = (*bond_iter)->bends[0];
        (*bond_iter)->bends[0]->remove();
        delete temp;
    }

    Bond *temp = *bond_iter;
    // put last bond in place of *bond_iter's place
    *bond_iter = temp->from_ptr->bonds.back();

    // remove last element 
    temp->from_ptr->bonds.pop_back();
    // delete temp (the object that *bond_iter pointed to)
    delete temp;

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

void Graph::delBend(Node *mid, Node *a, Node *b)
{
    // find the pointer in the list of bends of node mid that
    // correspond to the bend a-mid-b 
    std::vector<Bend*>::iterator bend_iter = mid->bends.begin();
    while( bend_iter != mid->bends.end() ) {
        if( (*bend_iter)->bond_a->to_ptr == a 
                and (*bend_iter)->bond_b->to_ptr == b ){
            break;
        }
        if( (*bend_iter)->bond_a->to_ptr == b 
                and (*bend_iter)->bond_b->to_ptr == a ){
            break;
        }

        ++bend_iter;
    }
    
    if(bend_iter == mid->bends.end() ) return; 
    Bend *bend_ptr = *bend_iter;

    // delete the bend pointers from the list in the bonds
    std::vector<Bend*>::iterator bi = bend_ptr->bond_a->bends.begin();
    while( bi != bend_ptr->bond_a->bends.end() ){
        if( *bi == bend_ptr ) break;
        ++bi;
    }
    *bi = bend_ptr->bond_a->bends.back();
    bend_ptr->bond_a->bends.pop_back();

    bi = bend_ptr->bond_b->bends.begin();
    while( bi != bend_ptr->bond_b->bends.end() ){
        if( *bi == bend_ptr ) break;
        ++bi;
    }
    *bi = bend_ptr->bond_b->bends.back();
    bend_ptr->bond_b->bends.pop_back();


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

// Node
Graph::Node::~Node()
{
    // delete all bend objects
    std::vector<Bend*>::iterator bend_iter = bends.begin();
    while(bend_iter != bends.end() ) {
        (*bend_iter)->remove();
        delete *bend_iter;
        ++bend_iter;
    }

    // delete all bond objects
    std::vector<Bond*>::iterator bond_iter = bonds.begin();
    while(bond_iter != bonds.end() ) {
        delete *bond_iter;
        ++bond_iter;
    }
}

// Bend
Graph::Bend::Bend(Node *mid, Bond *ba, Bond *bb, int i)
{
    index = i;
    nextBend = nullptr;
    prevBend = nullptr;
    bond_a = ba;
    bond_b = bb;

    // add this bend to the list of bends in bonds bond_a and bond_b    
    bond_a->bends.push_back(this);
    bond_b->bends.push_back(this);


}

void Graph::Bend::remove()
{

    // remove from Bond list bond_a
    std::vector<Bend*>::iterator bend_iter = bond_a->bends.begin();
    while(bend_iter != bond_a->bends.end() ) {
        if( (*bend_iter) == this ) break;
        ++bend_iter;
    }
    *bend_iter = bond_a->bends.back();
    bond_a->bends.pop_back();

    // remove from Bond list bond_b
    bend_iter = bond_b->bends.begin();
    while(bend_iter != bond_b->bends.end() ) {
        if( (*bend_iter) == this ) break;
        ++bend_iter;
    }
    *bend_iter = bond_b->bends.back();
    bond_b->bends.pop_back();
    

    // set nextBend of prevBend
    if( nextBend != nullptr ) nextBend->prevBend = nullptr;
    // set prevBend of nextBend
    if( prevBend != nullptr ) prevBend->nextBend = nullptr;
    
    // set values 
    nextBend = nullptr;
    prevBend = nullptr;
    bond_a = nullptr;
    bond_b = nullptr;
}


#endif
