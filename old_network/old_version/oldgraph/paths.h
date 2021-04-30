#ifndef GUARD_PATHS_H
#define GUARD_PATHS_H


template<class G>
class Paths
{
  public:
    Paths( G &g, int s);

    bool hasPathTo(int v) const {
        return marked[v];
    }

  private:
    int s;
    G& graph;

    std::vector<bool> marked;
    std::vector<int> edgeTo;
    void mark_neighbours( int vi);
};


template<class G>
Paths<G>::Paths( G& g, int s)
: s(s), graph(g), marked(g.V()), edgeTo(g.V())
{ mark_neighbours(s); }

template<class G>
void Paths<G>::mark_neighbours( int vi)
{
    marked[vi] = true;
    int Nneighbours = graph.Nneighbours(vi);
    int nIndex;
    for(int ni =0; ni<Nneighbours; ++ni ) {
        nIndex = graph.neighbourIndex(vi,ni);
        if( marked[ nIndex ] == false ) {
            edgeTo[ nIndex ] = vi;
            mark_neighbours( nIndex );
        }
    }

}

#endif
