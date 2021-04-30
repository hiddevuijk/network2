#ifndef GUARD_PATHS_BFS_H
#define GUARD_PATHS_BFS_H

#include "queue_ll.h"

template<class G>
class PathsBFS
{
  public:
    PathsBFS( G &g, int s);

    bool hasPathTo(int v) const {
        return marked[v];
    }

    int distance(int v) const {
        return distTo[v];
    }

  private:
    int s;
    G& graph;

    std::vector<bool> marked;
    std::vector<int> distTo;
    std::vector<int> edgeTo;
    void mark_neighbours( int vi);
};


template<class G>
PathsBFS<G>::PathsBFS( G& g, int s)
: s(s), graph(g), marked(g.V() ), distTo(g.V(), -1), edgeTo(g.V())
{ mark_neighbours(s); }

template<class G>
void PathsBFS<G>::mark_neighbours( int s)
{
    Queue<int> q;
    q.enqueue(s); 
    marked[s] = true;
    int dist = 0;
    distTo[s] = dist;
    int vi, Nneighbours, nIndex;
    while( q.size() > 0 ) {
        vi = q.dequeue();
        Nneighbours = graph.Nneighbours(vi);
        //dist += 1;
        dist = distTo[vi] + 1;
        // ni = neighbour number
        // nIndex = vertex index corresponding to ni
        for(int ni =0; ni<Nneighbours; ++ni ) {
            nIndex = graph.neighbourIndex(vi,ni);
            if(  !marked[ nIndex ] ) {
                marked[ nIndex ] = true;
                edgeTo[ nIndex ] = vi;
                distTo[ nIndex ] = dist;
                q.enqueue(nIndex);
            }
        }
    }
}

#endif
