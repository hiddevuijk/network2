#ifndef GUARD_GRAPH_H
#define GUARD_GRAPH_H

#include "vec2.h"
#include <vector>


class Graph {
public:
    Graph();
    ~Graph();

    class Bond {
    public:
        Node *n1, *n2;
    };
    
    class Bend {
    public:
        Node *n;
        Bond *b1, *b2;
        Bend *next;
        int polymer;
    };

    class Node {
    public:
        Node(): i(-1) {};

        int i;
        vec2 r;

        std::vector<Bond*> bond;
        std::vector<Bend*> bends;
    };

    // Member functions
    int addNode(); // adds a node, returns index of new node
    void addBond(int i, int j);

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
    Node *new_node = new Node();
    nodes.push_back(new_node);

    return nodes.size()-1;
}

void Graph::addBond(int ni, int nj)
{
    Bond *new_bond = new Bond(); 
    new_bond->n1 = ni; 
    new_bond->n2 = nj; 
   
     

}


#endif
