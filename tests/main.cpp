#include <iostream>

#include "graph.hpp"
#include "executor.hpp"
#include "generator.hpp"


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

    const std::string helpContent = "[HELP] You can only pass the graph generation or a filename.\n"
                                    "[HELP] $ ./main [generation method]\n"
                                    "[HELP] Where generation method is among: tree, deptree, 2dgrid, doubletree, filename\n";

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
        someDeps = Generator::GenerateDoubleDepTreeTasks(16);
    }
    else{
        const std::vector<std::string> methodNames{"tree", "deptree", "2dgrid", "doubletree"};
        const size_t choice = std::distance(methodNames.begin(), std::find(methodNames.begin(), methodNames.end(), params[1]));
        switch(choice){
        case 0:
            std::cout << "[INFO] use tree\n";
            someDeps = Generator::GenerateDepTreeTasks(8);
            break;
        case 1:
            std::cout << "[INFO] use deptree\n";
            someDeps = Generator::GenerateDoubleDepTreeTasks(16);
            break;
        case 2:
            std::cout << "[INFO] use 2dgrid\n";
            someDeps = Generator::Generate2DGrid(8);
            break;
        case 3:
            std::cout << "[INFO] use doubletree\n";
            someDeps = Generator::GenerateDoubleDepTreeTasks(16);
            break;
        default:
            std::cout << "[INFO] load " << params[1] << "\n";
            someDeps = Generator::LoadEdgeFile(params[1]);
            break;
        }
    }

    {
        const int cutGraphLimit = 3000;
        if(cutGraphLimit < someDeps.first){
            int idxEdge = 0;
            while(idxEdge < int(someDeps.second.size())){
                if(cutGraphLimit <= someDeps.second[idxEdge].first
                        || cutGraphLimit <= someDeps.second[idxEdge].second){
                    someDeps.second[idxEdge] = someDeps.second.back();
                    someDeps.second.pop_back();
                }
                else{
                    idxEdge += 1;
                }
            }
            someDeps.first = cutGraphLimit;
        }
    }

    //////////////////////////////////////////////////////////////////////////

    const int nbThreads = 2;
    const int partMaxSize = 10;
    const int partMinSize = partMaxSize;
    std::cout << "nbThreads : " << nbThreads << " / partMinSize : " << partMinSize << " / partMaxSize : " << partMaxSize << "\n";
    {
        Graph aGraph(someDeps.first, someDeps.second);
        std::cout << "Number of nodes : " << aGraph.getNbNodes() << "\n";
        assert(aGraph.isDag());
        std::pair<int,double> degGraph = aGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism one the original graph : " << degGraph.first << "  " << degGraph.second << "\n";
        aGraph.saveToDot("/tmp/agraph-original.dot");
        std::cout << "Generate pdf of the original graph with: dot -Tpdf /tmp/agraph-original.dot -o /tmp/agraph-original.pdf\n";

        aGraph.partitionRandom(partMaxSize);
        aGraph.saveToDot("/tmp/agraph-rand.dot");
        std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/agraph-rand.dot -o /tmp/agraph-rand.pdf\n";

        Graph depGraph = aGraph.getPartitionGraph();
        assert(depGraph.isDag());
        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism after random partitioning : " << degPar.first << "  " << degPar.second << "\n";
        std::cout << "Number of partitions : " << depGraph.getNbNodes() << " -- avg part size : " << double(aGraph.getNbNodes())/double(depGraph.getNbNodes()) << "\n";

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace-rand.svg", events, nbThreads);
    }
    {
        Graph aGraph(someDeps.first, someDeps.second);
        assert(aGraph.isDag());
        aGraph.partitionGreedy(partMaxSize);
        aGraph.saveToDot("/tmp/agraph-greedy.dot");
        std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/agraph-greedy.dot -o /tmp/agraph-greedy.pdf\n";

        Graph depGraph = aGraph.getPartitionGraph();
        assert(depGraph.isDag());
        depGraph.saveToDot("/tmp/depgraph-greedy.dot");
        std::cout << "Generate pdf of the final partition graph with greedy method with: dot -Tpdf /tmp/depgraph-greedy.dot -o /tmp/depgraph-greedy.pdf\n";

        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism after greedy partitioning : " << degPar.first << "  " << degPar.second << "\n";
        std::cout << "Number of partitions : " << depGraph.getNbNodes() << " -- avg part size : " << double(aGraph.getNbNodes())/double(depGraph.getNbNodes()) << "\n";

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace-greedy.svg", events, nbThreads);
    }
    {
        Graph aGraph(someDeps.first, someDeps.second);
        assert(aGraph.isDag());
        aGraph.partitionBacktrack(partMaxSize);
        aGraph.saveToDot("/tmp/agraph-backtrack.dot");
        std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/agraph-backtrack.dot -o /tmp/agraph-backtrack.pdf\n";

        Graph depGraph = aGraph.getPartitionGraph();
        assert(depGraph.isDag());
        depGraph.saveToDot("/tmp/depgraph-backtrack.dot");
        std::cout << "Generate pdf of the final partition graph with: dot -Tpdf /tmp/depgraph-backtrack.dot -o /tmp/depgraph-backtrack.pdf\n";

        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism after backtrack partitioning : " << degPar.first << "  " << degPar.second << "\n";
        std::cout << "Number of partitions : " << depGraph.getNbNodes() << " -- avg part size : " << double(aGraph.getNbNodes())/double(depGraph.getNbNodes()) << "\n";

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace-backtrack.svg", events, nbThreads);
    }
    {
        Graph aGraph(someDeps.first, someDeps.second);
        assert(aGraph.isDag());
        aGraph.partition(partMinSize,partMaxSize);
        aGraph.saveToDot("/tmp/agraph.dot");
        std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/agraph.dot -o /tmp/agraph.pdf\n";

        Graph depGraph = aGraph.getPartitionGraph();
        assert(depGraph.isDag());
        depGraph.saveToDot("/tmp/depgraph.dot");
        std::cout << "Generate pdf of the final partition graph with: dot -Tpdf /tmp/depgraph.dot -o /tmp/depgraph.pdf\n";

        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism after advanced partitioning : " << degPar.first << "  " << degPar.second << "\n";
        std::cout << "Number of partitions : " << depGraph.getNbNodes() << " -- avg part size : " << double(aGraph.getNbNodes())/double(depGraph.getNbNodes()) << "\n";

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace.svg", events, nbThreads);
    }
    {
        Graph aGraph(someDeps.first, someDeps.second);
        assert(aGraph.isDag());
        aGraph.partitionHorizontal(partMaxSize);
        aGraph.saveToDot("/tmp/agraph-horizontal.dot");
        std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/agraph-horizontal.dot -o /tmp/agraph-horizontal.pdf\n";

        Graph depGraph = aGraph.getPartitionGraph();
        assert(depGraph.isDag());
        depGraph.saveToDot("/tmp/depgraph-horizontal.dot");
        std::cout << "Generate pdf of the horizontal partition graph with: dot -Tpdf /tmp/depgraph-horizontal.dot -o /tmp/depgraph-horizontal.pdf\n";

        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism after horizontal partitioning : " << degPar.first << "  " << degPar.second << "\n";
        std::cout << "Number of partitions : " << depGraph.getNbNodes() << " -- avg part size : " << double(aGraph.getNbNodes())/double(depGraph.getNbNodes()) << "\n";

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace.svg", events, nbThreads);
    }
    {
        Graph aGraph(someDeps.first, someDeps.second);
        assert(aGraph.isDag());
        aGraph.partitionDiamond(/*partMaxSize*/4);
        aGraph.saveToDot("/tmp/agraph-diamond.dot");
        std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/agraph-diamond.dot -o /tmp/agraph-diamond.pdf\n";

        Graph depGraph = aGraph.getPartitionGraph();
        assert(depGraph.isDag());
        depGraph.saveToDot("/tmp/depgraph-diamond.dot");
        std::cout << "Generate pdf of the diamond partition graph with: dot -Tpdf /tmp/depgraph-diamond.dot -o /tmp/depgraph-diamond.pdf\n";

        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << "Degree of parallelism after diamond partitioning : " << degPar.first << "  " << degPar.second << "\n";
        std::cout << "Number of partitions : " << depGraph.getNbNodes() << " -- avg part size : " << double(aGraph.getNbNodes())/double(depGraph.getNbNodes()) << "\n";

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace.svg", events, nbThreads);
    }

    return 0;
}
