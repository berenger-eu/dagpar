#include <iostream>

#include "graph.hpp"
#include "executor.hpp"

std::vector<std::pair<int,int>> GenerateBinaryTreeTasks(const int inHeight){
    std::vector<std::pair<int,int>> someDeps;
    int cellsOffset = 0;
    for(int level = 1 ; level < inHeight ; ++level){
        const int nbCellsAtPreviousLevel = (1 << (level-1));
        const int nbCellsAtLevel = (1 << level);
        for(int idxCell = 0 ; idxCell < nbCellsAtLevel ; ++idxCell){
            someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell/2, cellsOffset + nbCellsAtPreviousLevel + idxCell});
        }
        cellsOffset += nbCellsAtPreviousLevel;
    }
    return someDeps;
}

std::vector<std::pair<int,int>> GenerateDepTreeTasks(const int inHeight){
    std::vector<std::pair<int,int>> someDeps;
    int cellsOffset = 0;
    for(int level = 1 ; level < inHeight ; ++level){
        const int nbCellsAtPreviousLevel = level;
        const int nbCellsAtLevel = (level+1);

        someDeps.push_back(std::pair<int,int>{cellsOffset, cellsOffset + nbCellsAtPreviousLevel});

        for(int idxCell = 1 ; idxCell < nbCellsAtLevel-1 ; ++idxCell){
            someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell-1, cellsOffset + nbCellsAtPreviousLevel + idxCell});
            someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell, cellsOffset + nbCellsAtPreviousLevel + idxCell});
        }

        someDeps.push_back(std::pair<int,int>{cellsOffset + nbCellsAtPreviousLevel-1, cellsOffset + nbCellsAtPreviousLevel + nbCellsAtLevel-1});

        cellsOffset += nbCellsAtPreviousLevel;
    }
    return someDeps;
}

int main(){    
    //std::vector<std::pair<int,int>> someDeps = GenerateBinaryTreeTasks(5);
    std::vector<std::pair<int,int>> someDeps = GenerateDepTreeTasks(8);

    const int nbThreads = 2;
    const int partMinSize = 2;
    const int partMaxSize = 2;
    Graph aGraph(someDeps);
    aGraph.partition(partMinSize,partMaxSize,nbThreads);
    aGraph.saveToDot("/tmp/agraph.dot");

    Graph depGraph = aGraph.getPartitionGraph();
    depGraph.saveToDot("/tmp/depgraph.dot");

    int duration;
    std::vector<Executor::Event> events;
    std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
    Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace.svg", events, nbThreads);

    return 0;
}
