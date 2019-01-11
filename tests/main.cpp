#include <iostream>

#include "graph.hpp"
#include "executor.hpp"

int main(){
    const int nbThreads = 2;
    std::vector<std::pair<int,int>> someDeps{{0,1}, {0,2}, {1,2}};
    Graph aGraph(someDeps);
    aGraph.partition(1,1,nbThreads);
    aGraph.saveToDot("/tmp/agraph.dot");

    Graph depGraph = aGraph.getPartitionGraph();
    aGraph.saveToDot("/tmp/depgraph.dot");

    int duration;
    std::vector<Executor::Event> events;
    std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
    Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace.svg", events, nbThreads);

    return 0;
}
