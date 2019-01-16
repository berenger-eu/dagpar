#include <iostream>

#include "graph.hpp"
#include "executor.hpp"

///////////////////////////////////////////////////////////
///
/// Methods to generate some graphs
///
///////////////////////////////////////////////////////////

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
    return std::pair<int, std::vector<std::pair<int,int>>>(someDeps.back().second+1,someDeps);
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
    return std::pair<int, std::vector<std::pair<int,int>>>(someDeps.back().second+1,someDeps);
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

    return std::pair<int, std::vector<std::pair<int,int>>>(someDeps.back().second+1,someDeps);
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

///////////////////////////////////////////////////////////
///
/// Main:
/// Generate a graph
/// Test the random method
/// Test the greedy method
/// Test the advanced method
///
///////////////////////////////////////////////////////////

int main(int argc, char** argv){
    std::vector<std::string> params;
    params.insert(params.end(), argv, argv+argc);

    const std::string helpContent = "[HELP] You can only pass the graph generation.\n"
                                    "[HELP] $ ./main [generation method]\n"
                                    "[HELP] Where generation method is among: tree, deptree, 2dgrid, doubletree\n";

    if(params.size() > 2){
        std::cout << "[ERROR] Invalid number of parameters.\n" << helpContent;
        return 1;
    }
    if(params.size() == 2 && params[1] == "--help"){
        std::cout << "[HELP] Asked for help.\n" << helpContent;
        return 1;
    }

    //////////////////////////////////////////////////////////////////////////

    std::pair<int, std::vector<std::pair<int,int>>> someDeps;
    if(params.size() < 2){
        std::cout << "[INFO] use doubletree\n";
        someDeps = GenerateDoubleDepTreeTasks(16);
    }
    else{
        const std::vector<std::string> methodNames{"tree", "deptree", "2dgrid", "doubletree"};
        const size_t choice = std::distance(methodNames.begin(), std::find(methodNames.begin(), methodNames.end(), params[1]));
        switch(choice){
        case 0:
            std::cout << "[INFO] use tree\n";
            someDeps = GenerateDepTreeTasks(8);
            break;
        case 1:
            std::cout << "[INFO] use deptree\n";
            someDeps = GenerateDoubleDepTreeTasks(16);
            break;
        case 2:
            std::cout << "[INFO] use 2dgrid\n";
            someDeps = Generate2DGrid(8);
            break;
        case 3:
        default:
            std::cout << "[INFO] use doubletree\n";
            someDeps = GenerateDoubleDepTreeTasks(16);
            break;
        }
    }

    //////////////////////////////////////////////////////////////////////////

    const int nbThreads = 2;
    const int partMinSize = 2;
    const int partMaxSize = 4;
    std::cout << "nbThreads : " << nbThreads << " / partMinSize : " << partMinSize << " / partMaxSize : " << partMaxSize << "\n";
    {
        Graph aGraph(someDeps.first, someDeps.second);
        std::pair<int,double> degGraph = aGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism one the original graph : " << degGraph.first << "  " << degGraph.second << "\n";
        aGraph.saveToDot("/tmp/agraph-original.dot");
        std::cout << "Generate pdf of the original graph with: dot -Tpdf /tmp/agraph-original.dot -o /tmp/agraph-original.pdf\n";

        aGraph.partitionRandom(partMaxSize);
        aGraph.saveToDot("/tmp/agraph-rand.dot");
        std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/agraph-rand.dot -o /tmp/agraph-rand.pdf\n";

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
        aGraph.partitionGreedy(partMaxSize);
        aGraph.saveToDot("/tmp/agraph-greedy.dot");
        std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/agraph-greedy.dot -o /tmp/agraph-greedy.pdf\n";

        Graph depGraph = aGraph.getPartitionGraph();
        depGraph.saveToDot("/tmp/depgraph-greedy.dot");
        std::cout << "Generate pdf of the final partition graph with greedy method with: dot -Tpdf /tmp/depgraph-greedy.dot -o /tmp/depgraph-greedy.pdf\n";

        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism after greedy partitioning : " << degPar.first << "  " << degPar.second << "\n";

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace-greedy.svg", events, nbThreads);
    }
    {
        Graph aGraph(someDeps.first, someDeps.second);
        aGraph.partition(partMinSize,partMaxSize,nbThreads);
        aGraph.saveToDot("/tmp/agraph.dot");
        std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/agraph.dot -o /tmp/agraph.pdf\n";

        Graph depGraph = aGraph.getPartitionGraph();
        depGraph.saveToDot("/tmp/depgraph.dot");
        std::cout << "Generate pdf of the final partition graph with: dot -Tpdf /tmp/depgraph.dot -o /tmp/depgraph.pdf\n";

        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism after depth partitioning : " << degPar.first << "  " << degPar.second << "\n";

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace.svg", events, nbThreads);
    }

    return 0;
}
