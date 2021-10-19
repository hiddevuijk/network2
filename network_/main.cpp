
#include "graph.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


int main()
{
  Graph graph(10);
  graph.AddNode();
  graph.DelNode(0);

  cout << graph.NumberOfNodes() << endl;
  cout << graph.nodes_[0]->index << endl;

  return 0;
}


