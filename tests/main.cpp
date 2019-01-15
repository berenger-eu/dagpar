#include <iostream>

#include "graph.hpp"
#include "executor.hpp"

std::pair<int, std::vector<std::pair<int,int>>> GenerateBinaryTreeTasks(const int inHeight){
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
    return std::pair<int, std::vector<std::pair<int,int>>>(cellsOffset,someDeps);
}

std::pair<int, std::vector<std::pair<int,int>>> GenerateDepTreeTasks(const int inHeight){
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
    return std::pair<int, std::vector<std::pair<int,int>>>(cellsOffset,someDeps);
}

std::pair<int, std::vector<std::pair<int,int>>> GenerateDoubleDepTreeTasks(const int inHeight){
    std::vector<std::pair<int,int>> someDeps;
    int cellsOffset = 0;
    for(int level = 1 ; level < (inHeight+1)/2 ; ++level){
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

    if((inHeight & 1) == 0){
        const int nbCellsAtPreviousLevel = inHeight/2;
        const int nbCellsAtLevel = inHeight/2;

        for(int idxCell = 0 ; idxCell < nbCellsAtLevel ; ++idxCell){
            someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell, cellsOffset + nbCellsAtPreviousLevel + idxCell});
        }

        cellsOffset += nbCellsAtPreviousLevel;
    }

    for(int level = (inHeight+2)/2 ; level < inHeight ; ++level){
        const int nbCellsAtPreviousLevel = (inHeight+1) - level;
        const int nbCellsAtLevel = (inHeight+1) - level - 1;

        for(int idxCell = 0 ; idxCell < nbCellsAtLevel ; ++idxCell){
            someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell, cellsOffset + nbCellsAtPreviousLevel + idxCell});
            someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell+1, cellsOffset + nbCellsAtPreviousLevel + idxCell});
        }

        cellsOffset += nbCellsAtPreviousLevel;
    }

    return std::pair<int, std::vector<std::pair<int,int>>>(cellsOffset,someDeps);
}

std::pair<int, std::vector<std::pair<int,int>>> Generate2DGrid(const int inGridDim){
    std::vector<std::pair<int,int>> someDeps;
    for(int idxRow = 1 ; idxRow < inGridDim ; ++idxRow){
        for(int idxCol = 1 ; idxCol < inGridDim ; ++idxCol){
            someDeps.push_back(std::pair<int,int>{(idxCol*inGridDim)+idxRow-1, (idxCol*inGridDim)+idxRow});
            someDeps.push_back(std::pair<int,int>{((idxCol-1)*inGridDim)+idxRow, (idxCol*inGridDim)+idxRow});
        }
    }
    return std::pair<int, std::vector<std::pair<int,int>>>(inGridDim*inGridDim,someDeps);
}

int main(){    
    //std::pair<int, std::vector<std::pair<int,int>>> someDeps = GenerateBinaryTreeTasks(6);
    //std::pair<int, std::vector<std::pair<int,int>>> someDeps = GenerateDepTreeTasks(8);
    std::pair<int, std::vector<std::pair<int,int>>> someDeps = Generate2DGrid(8);
    //std::pair<int, std::vector<std::pair<int,int>>> someDeps = GenerateDoubleDepTreeTasks(16);

    const int nbThreads = 2;
    const int partMinSize = 2;
    const int partMaxSize = 4;
    std::cout << "nbThreads : " << nbThreads << " / partMinSize : " << partMinSize << " / partMaxSize : " << partMaxSize << "\n";
    {
        Graph aGraph(someDeps.first, someDeps.second);
        std::pair<int,double> degGraph = aGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism one the original graph : " << degGraph.first << "  " << degGraph.second << "\n";

        aGraph.partitionRandom(partMaxSize);
        aGraph.saveToDot("/tmp/agraph-rand.dot");

        Graph depGraph = aGraph.getPartitionGraph();
        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism after random partitioning : " << degPar.first << "  " << degPar.second << "\n";

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace-rand.svg", events, nbThreads);
    }
    {
        Graph aGraph(someDeps.first, someDeps.second);
        aGraph.partition(partMinSize,partMaxSize,nbThreads);
        aGraph.saveToDot("/tmp/agraph.dot");

        Graph depGraph = aGraph.getPartitionGraph();
        depGraph.saveToDot("/tmp/depgraph.dot");

        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism after depth partitioning : " << degPar.first << "  " << degPar.second << "\n";

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace.svg", events, nbThreads);
    }

    return 0;
}
