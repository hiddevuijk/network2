#ifndef GUARD_GRAPH_H
#define GUARD_GRAPH_H

#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

class Graph
{
  public:
    // create a graph with V vertices
    Graph(int V);
    // create a graph from input stream
    Graph( std::ifstream& in);

    // add edge between v and w
    void addEdge(int v, int w);
    int addVertex();
    void deleteEdge(int v, int w);
    void deleteVertex(int v);

    // return # vertices
    int V() const { return NV; }
    // return # edges
    int E() const { return NE; }

    void showAdj() const;

    const std::vector<int>::const_iterator get_iterator(int v)
            const { return adj[v].begin(); }
    int Nneighbours(int vi) const { return adj[vi].size(); }

    int neighbourIndex(int vi, int ni) const { return adj[vi][ni]; }
  private:
    int NV; 
    int NE;

    std::vector<std::vector<int> > adj;

    void exchangeVertexIndex(int v, int w);
    void freeVertex(int v);

};

void Graph::showAdj() const
{
    
    for(int iv=0; iv < NV; ++iv) {
        std::cout << iv << " : ";
        for(std::vector<int>::size_type ie=0;ie< adj[iv].size(); ++ie) {
            std::cout << adj[iv][ie] << " ";
        }
        std::cout << std::endl;

    }
}

void Graph::addEdge(int v, int w)
{
    NE += 1;
    adj[v].push_back(w);
    adj[w].push_back(v);

}

int Graph::addVertex() 
{
    NV += 1;
    adj.push_back( std::vector<int>() );

    //return the index of the created vertex
    return NV-1;
}
void Graph::deleteEdge(int v, int w)
{
    std::vector<int>::iterator it = std::find(adj[v].begin(), adj[v].end(), w);
    if( it != adj[v].end() ) {
        adj[v].erase(it);
        it = std::find(adj[w].begin(), adj[w].end(), v);
        adj[w].erase(it);
    }

}

void Graph::exchangeVertexIndex(int v, int w)
{
    // exchange location in adj
    std::vector<int> adjv = adj[v];
    adj[v] = adj[w];
    adj[w] = adjv;

    // in adj exchange v <--> w
    for(int i=0; i < NV; ++i) {
        for(std::vector<int>::size_type j=0;j< adj[i].size(); ++j) {
            if( adj[i][j] == v ) adj[i][j] = w;
            else if( adj[i][j] == w ) adj[i][j] = v;
        }
    } 
}

void Graph::freeVertex(int v)
{
    std::vector<int>::iterator it;
    for(int i=0; i<NV; ++i) {
        it = std::find(adj[i].begin(), adj[i].end(), v);
        while( it != adj[i].end() ) {
            adj[i].erase(it);
            it = std::find(adj[i].begin(), adj[i].end(), v);
        }
    }
}

void Graph::deleteVertex(int v)
{
    // remove all connection toand from vertex v
    freeVertex(v);
    // exchange v with the last index
    exchangeVertexIndex(v, NV-1);
    // remove the vertex from adj list
    adj.erase( adj.end()-1 );
    NV -= 1;
}

Graph::Graph( std::ifstream& in)
{
    NE = 0;
    in >> NV;
    adj = std::vector<std::vector<int> >(NV);
    int v,w;
    while( ( in >> v ) and ( in >> w) ){
        addEdge(v,w);
    }

}

#endif
