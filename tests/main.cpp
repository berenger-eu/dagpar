#include <iostream>

#include "graph.hpp"
#include "executor.hpp"

int main(){    
    std::vector<std::pair<int,int>> someDeps;
    for(int idxTask = 1 ; idxTask < 10 ; ++idxTask){
        someDeps.push_back(std::pair<int,int>{(idxTask-1)/2, idxTask});
    }

    const int nbThreads = 2;
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
