#include <iostream>
#include <functional>

#include "graph.hpp"
#include "executor.hpp"
#include "generator.hpp"
#include "utils.hpp"

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
    const std::string helpContent = "[HELP] $ ./executeone [options]\n"
                                    "[HELP] Where options are:\n"
                                    "[HELP] --help to get help\n"
                                    "[HELP] Or\n"
                                    "[HELP] --filename [filename]\n"
                                    "[HELP] --nbt [number of threads to execute the dag]\n"
                                    "[HELP] --opt [overhead per task, % of the total execution]\n"
                                    "[HELP] --oppu [overhead per push, % of the total execution]\n"
                                    "[HELP] --oppo [overhead per pop, % of the total execution]\n"
                                    "[HELP] --cs [cluster size]\n"
                                    "[HELP] --h1 [height one for the cluster method]\n"
                                    "[HELP] --h2 [height two for the cluster method]\n"
                                    "[HELP] --export-dot (optional, to export the dot files)\n"
                                    "[HELP] Example:\n"
                                    "[HELP] ./executeone --filename ../data/chukrut_disque.dot --nbt 5 --opt 0 --oppu 0 --oppo 0 --cs 10 --h1 3 --h2 2\n";

    Utils::ParamHelper params(argc, argv);

    if(params.paramExist({"--help", "--h", "-h", "-help"})){
        std::cout << "[HELP] Asked for help.\n" << helpContent;
        return 1;
    }

    const std::string filename = params.getStr({"--filename", "--f", "-filename", "-f"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at filename.\n" << helpContent;
        return 1;
    }

    const int nbThreads = params.getValue<int>({"--nbt", "-nbt"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at nbt.\n" << helpContent;
        return 1;
    }

    const double overheadPerTask = params.getValue<int>({"--opt", "-opt"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at opt.\n" << helpContent;
        return 1;
    }

    const double overheadPerPush = params.getValue<int>({"--oppu", "-oppu"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at oppu.\n" << helpContent;
        return 1;
    }

    const double overheadPerPop = params.getValue<int>({"--oppo", "-oppo"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at oppo.\n" << helpContent;
        return 1;
    }

    const int maxSize = params.getValue<int>({"--cs", "-cs"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at cs.\n" << helpContent;
        return 1;
    }

    const int h1 = params.getValue<int>({"--h1", "-h1"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at h1.\n" << helpContent;
        return 1;
    }

    const int h2 = params.getValue<int>({"--h2", "-h2"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at h2.\n" << helpContent;
        return 1;
    }

    std::cout << "filename = " << filename << std::endl;
    std::cout << "nbThreads = " << nbThreads << std::endl;
    std::cout << "overheadPerTask = " << overheadPerTask << std::endl;
    std::cout << "overheadPerPush = " << overheadPerPush << std::endl;
    std::cout << "overheadPerPop = " << overheadPerPop << std::endl;
    std::cout << "maxSize = " << maxSize << std::endl;
    std::cout << "h1 = " << h1 << std::endl;
    std::cout << "h2 = " << h2 << std::endl;


    const bool exportDot = params.paramExist({"--export-dot", "-export-dot"});

    //////////////////////////////////////////////////////////////////////////

    std::pair<int, std::vector<std::pair<int,int>>> someDeps = Generator::LoadEdgeFile(filename);

    if(someDeps.second.size() == 0){
        std::cout << "[INFO] file is empty, exit...\n";
        return 1;
    }

    std::vector<double> costs = Generator::GetCostIfExist(filename);

    //////////////////////////////////////////////////////////////////////////

    {
        std::cout << "Original graph:\n";
        Graph aGraph(someDeps.first, someDeps.second);
        std::cout << " - Number of nodes : " << aGraph.getNbNodes() << "\n";
        assert(aGraph.isDag());
        std::pair<int,double> degGraph = aGraph.estimateDegreeOfParallelism();
        std::cout << " - Degree of parallelism one the sequential graph : " << degGraph.first << "  " << degGraph.second << "\n";

        if(costs.size()){
            for(int idxNode = 0 ; idxNode < aGraph.getNbNodes() ; ++idxNode){
                auto node = aGraph.getNode(idxNode);
                assert(node->getId() < int(costs.size()));
                node->setCost(costs[node->getId()]);
            }
        }

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(aGraph, nbThreads, overheadPerTask, overheadPerPop, overheadPerPush);
        std::cout << " - Without clustering duration = " << duration << "\n";
    }

    //////////////////////////////////////////////////////////////////////////

    const std::vector<std::pair<std::string,std::function<void(Graph&,int)>>> allPartitionMethods= {
                                    {"diamond", [](Graph& graph, const int clusterSize){
                                        graph.partitionDiamond(clusterSize);
                                     }},
                                    {"final", [h1,h2](Graph& graph, const int clusterSize){
                                        graph.partitionFinal(clusterSize, h1,h2);
                                     }},
                                    {"final-with-neighbor-rafinement", [h1,h2](Graph& graph, const int clusterSize){
                                        graph.partitionFinalWithNeighborRefinement(clusterSize, h1, h2, clusterSize);
                                     }},
                                    {"final-with-emulated-rafinement", [h1, h2, overheadPerTask, nbThreads, overheadPerPop, overheadPerPush](Graph& graph, const int clusterSize){
                                        graph.partitionFinalWithEmulationRefinement(clusterSize, h1, h2, clusterSize, overheadPerTask, nbThreads, overheadPerPop, overheadPerPush);
                                     }},
                                    {"final-with-neighbor-rafinement-2", [h1, h2](Graph& graph, const int clusterSize){
                                        graph.partitionFinalWithNeighborRefinement(clusterSize/2, h1, h2, clusterSize);
                                     }},
                                    {"final-with-emulated-rafinement-2", [h1, h2, overheadPerTask, nbThreads, overheadPerPop, overheadPerPush](Graph& graph, const int clusterSize){
                                        graph.partitionFinalWithEmulationRefinement(clusterSize/2, h1, h2, clusterSize, overheadPerTask, nbThreads, overheadPerPop, overheadPerPush);
                                     }}
                                    };


    for(const auto& method : allPartitionMethods){
        std::cout << "=======================[" << method.first <<  "]======================================\n";
        Graph aGraph(someDeps.first, someDeps.second);
        std::cout << " - Number of nodes : " << aGraph.getNbNodes() << "\n";
        assert(aGraph.isDag());
        std::pair<int,double> degGraph = aGraph.estimateDegreeOfParallelism();
        std::cout << " - Degree of parallelism one the " << method.first << " graph : " << degGraph.first << "  " << degGraph.second << "\n";

        if(costs.size()){
            for(int idxNode = 0 ; idxNode < aGraph.getNbNodes() ; ++idxNode){
                auto node = aGraph.getNode(idxNode);
                assert(node->getId() < int(costs.size()));
                node->setCost(costs[node->getId()]);
            }
        }

        method.second(aGraph, maxSize);
        if(exportDot){
            aGraph.saveToDot("/tmp/agraph-" + method.first + ".dot");
        }

        Graph depGraph = aGraph.getPartitionGraph();
        assert(depGraph.isDag());
        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << " - Degree of parallelism after " << method.first << " partitioning : " << degPar.first << "  " << degPar.second << "\n";
        std::cout << " - Number of partitions : " << depGraph.getNbNodes() << " -- avg part size : " << double(aGraph.getNbNodes())/double(depGraph.getNbNodes()) << "\n";

        if(exportDot){
            depGraph.saveToDot("/tmp/depgraph-" + method.first + ".dot");
        }

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads, overheadPerTask, overheadPerPop, overheadPerPush);
        std::cout << " - with " << method.first << " clustering duration = " << duration << "\n";
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace-" + method.first + ".svg", events, nbThreads);
        std::cout << "=============================================================\n";
    }

    return 0;
}