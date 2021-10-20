#ifndef GRAPH_H
#define GRAPH_H

#include "vec2.h"

#include <vector>
#include <fstream>

#include <iostream> // for debug only


class Graph {

 public:
  // types
  // constructor, assignment, destructor
  Graph(int number_of_nodes = 0);
  ~Graph();

  // methods
  unsigned int NumberOfNodes() const { return (unsigned int) nodes_.size(); }

  // node manipolation
  int AddNode();
  int AddNode(double x, double y);
  void DelNode(int i);

  // bond manipulation
  bool AddBond(int i, int j, int xboundary = 0, int yboundary = 0);
  bool HasBond(int i, int j);
  void DelBond(int i, int j);

  // data members

 //private:
  // private types
  class Node;
  class Bond;
  //class Bend;

  // private constants 
  // private methods

  // private data members
  std::vector<Node*> nodes_;     // all nodes in the network

};


// private types

// The "index" of a node correponds to the indexin the "nodes"
//    vector of the Graph ( graph.nodes[i]->index = i )
//    all uninitialized nodes have index = -1
// The "positions" is a Vec2 with the x-y location of the node
// The "bonds" ("bends") vector has pointers to all bonds (bends) that are
//    attached to this node.
// The Node object is responsible for the memory menagement of the Bond and
//    Bend objects.

class Graph::Node {
 public:
  Node() : index(-1) {} // uninitialized node
  Node(int index) : index(index) {}
  Node(int index, double x, double y)
      : index(index), position(x, y) {}
  ~Node();

  // methods
  unsigned int NBonds() const { return (unsigned int) bonds.size(); }
  //unsigned int NBends() const { return (unsigned int) bends.size(); }

  // data members
  int index;
  Vec2 position;

  std::vector<Bond*> bonds; // all bonds that depend on this node
  //std::vector<Bend*> bends; // all bends that depend on this node
 
};


class Graph::Bond {
 public:
  Bond(Node *from_ptr = nullptr, Node *to_ptr = nullptr)
    : from_ptr(from_ptr), to_ptr(to_ptr) {}

  Bond(Node *from_ptr, Node *to_ptr, int xb, int yb)
    : from_ptr(from_ptr), to_ptr(to_ptr), xb(xb), yb(yb) {}
    
  //void free(); 

  Node *from_ptr;
  Node *to_ptr;

  int xb, yb;

  //double rest_length;

  void free(); 
};


// This object does not manage memory, but the destructor makes sure that
// no other objects depend on this bend.
// The index of a bend is equal to the correponding index of the pointer to
//    this bend in "bends" of graph.nodes[ mid->index ].bends[index] = index
// The pointer mid points to the middle node of the bend
// The pointers bond_a and bond_b point to bond objects that are the two
//    arms of the bend.
// The next_bend and prev_bend pointers point to bend ojects along the polymer
//    that this bend belongs to.
/*
class Graph::Bend {
 public:
  Bend(): index(-1), mid(nullptr), bond_a(nullptr), bond_b(nullptr),
          prev_bend(nullptr), next_bend(nullptr) {} 

  Bend(int index, Node* mid, Bond* bond_a, Bond* bond_b)
      : index(index), mid(mid), bond_a(bond_a), bond_b(bond_b),
        prev_bend(nullptr), next_bend(nullptr) {}

  ~Bend() { free(); }

  void free(); // removes all dependencies on this bend object

 
  // returns true if there is a bend along the next/prev direction of 
  // the polymer that this bend belongs to. 
  bool hasNext() const { return next_bend != nullptr; }
  bool hasPrev() const { return prev_bend != nullptr; }

  // returns the index of the node  TO DO
  //int prevNodeIndex() const;
  //int nextNodeIndex() const;

  int nodeIndex() const { return mid->index; }

  int index; 
  Node* mid;
  Bond* bond_a;
  Bond* bond_b;
  Bend* prev_bend; 
  Bend* next_bend;

  int polymer_index;
  //double theta; // rest angle NEEDED??



};


void Graph::Bend::free()
{
  // if this bend is part of a polymer, remove it from the polymer 
  if (prev_bend != nullptr) {
    prev_bend->next_bend = nullptr;
  }
  prev_bend = nullptr;

  if (next_bend != nullptr) {
    next_bend->prev_bend = nullptr;
  }
  next_bend = nullptr;

}
*/

// the Node destructor deletes a node, and removes and deletes
// all bond and bend objects in the graph that are dependent on this node
//  TO DO: discribe free 
Graph::Node::~Node()
{
    //delete bends

  for(std::vector<Bond*>::size_type i = 0; i < bonds.size(); ++i) {
    //bonds[i]->free();
    //delete bonds[i];
  }
  bonds.clear(); 

}

// constructor, assigmnent, destructor
Graph::Graph(int number_of_nodes)
: nodes_(number_of_nodes, nullptr)
{
  for (int i = 0; i < number_of_nodes; ++i)
    nodes_[i] = new Node(i);
}

Graph::~Graph()
{
  // delete all node objects
  for (std::vector<Node*>::size_type i = 0; i < nodes_.size(); ++i) {
    delete nodes_[i];
  }
}


// Graph member functions

int Graph::AddNode()
{
  nodes_.push_back(new Node(nodes_.size() ) );
  return nodes_.size(); 
}

int Graph::AddNode(double x, double y)
{
  nodes_.push_back(new Node(nodes_.size(), x, y) );
  return nodes_.size();
}

// Delete node[i] by putting the last node in position i
// and set the index
void Graph::DelNode(int i)
{
  Node *temp = nodes_[i];  
  nodes_[i] = nodes_.back();
  nodes_[i]->index = i;
  nodes_.pop_back();
  delete temp;
}

bool Graph::AddBond(int i, int j, int xb, int yb)
{ 
  if (HasBond(i, j) ) return false;

  nodes_[i]->bonds.push_back(new Bond(nodes_[i], nodes_[j],  xb,  yb) );
  nodes_[j]->bonds.push_back(new Bond(nodes_[j], nodes_[i], -xb, -yb) );
  return true;
}

bool Graph::HasBond(int i, int j)
{
  Node *ni_ptr = nodes_[i];
  Node *nj_ptr = nodes_[j];

  // check if there is a bond in the bond list of node i
  // that goes from node i to node j
  for (unsigned int bi = 0; bi < ni_ptr->bonds.size(); ++bi ) {
    if (ni_ptr->bonds[bi]->to_ptr == nj_ptr) return true;
  }

  return false;
}

// delete the bond between node i and j
// also delete all bending objects that depend on this bond
void Graph::DelBond(int i, int j)
{
  
}

#endif // GRAPH_H
