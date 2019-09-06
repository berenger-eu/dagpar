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
    const std::string helpContent = "[HELP] You can only pass the graph generation or a filename.\n"
                                    "[HELP] $ ./main [generation method]\n"
                                    "[HELP] Where generation method is among: tree, deptree, 2dgrid, doubletree, filename\n";

    Utils::ParamHelper params(argc, argv);

    if(params.getNbParams() == 1 || params.getNbParams() > 2){
        std::cout << "[ERROR] Invalid number of parameters.\n" << helpContent;
        return 1;
    }
    if(params.paramExist({"--help", "--h", "-h", "-help"})){
        std::cout << "[HELP] Asked for help.\n" << helpContent;
        return 1;
    }

    //////////////////////////////////////////////////////////////////////////

    std::pair<int, std::vector<std::pair<int,int>>> someDeps;
    std::vector<double> costs;

    if(params.getNbParams() != 2){
        std::cout << "[INFO] use doubletree\n";
        someDeps = Generator::GenerateDoubleDepTreeTasks(16);
    }
    else{
        const std::vector<std::string> methodNames{"tree", "deptree", "2dgrid", "doubletree"};
        const size_t choice = std::distance(methodNames.begin(), std::find(methodNames.begin(), methodNames.end(), argv[1]));
        switch(choice){
        case 0:
            std::cout << "[INFO] use tree\n";
            someDeps = Generator::GenerateDepTreeTasks(8);
            break;
        case 1:
            std::cout << "[INFO] use deptree\n";
            someDeps = Generator::GenerateDoubleDepTreeTasks(32);
            break;
        case 2:
            std::cout << "[INFO] use 2dgrid\n";
            someDeps = Generator::Generate2DGrid(5);
            break;
        case 3:
            std::cout << "[INFO] use doubletree\n";
            someDeps = Generator::GenerateDoubleDepTreeTasks(18);
            break;
        default:
            std::cout << "[INFO] load " << argv[1] << "\n";
            someDeps = Generator::LoadEdgeFile(argv[1]);
            costs = Generator::GetCostIfExist(argv[1]);
            if(someDeps.second.size() == 0){
                std::cout << "[INFO] file is empty, exit...\n";
                return 1;
            }
            break;
        }
    }

    //////////////////////////////////////////////////////////////////////////

    const int nbThreads = 40;

    const std::vector<std::pair<std::string,std::function<void(Graph&,int)>>> allPartitionMethods= {
                                    {"G", [](Graph& graph, const int clusterSize){
                                        graph.G(clusterSize);
                                     }},
                                    {"GUpdate", [](Graph& graph, const int clusterSize){
                                        graph.GUpdate(clusterSize);
                                     }},
                                    {"GStop", [](Graph& graph, const int clusterSize){
                                        graph.GStop(clusterSize);
                                     }}
                                    };

    for(int partMaxSize = 6 ; partMaxSize < 7 ; ++partMaxSize){

        for(const auto& method : allPartitionMethods){
            std::cout << "=======================[" << method.first <<  "]======================================\n";
            Graph aGraph(someDeps.first, someDeps.second);
            std::cout << "Number of nodes : " << aGraph.getNbNodes() << "\n";
            assert(aGraph.isDag());
            std::pair<int,double> degGraph = aGraph.estimateDegreeOfParallelism();
            std::cout << "Degree of parallelism one the " << method.first << " graph : " << degGraph.first << "  " << degGraph.second << "\n";

            if(costs.size()){
                for(int idxNode = 0 ; idxNode < aGraph.getNbNodes() ; ++idxNode){
                    auto node = aGraph.getNode(idxNode);
                    assert(node->getId() < int(costs.size()));
                    node->setCost(costs[node->getId()]);
                }
            }

            method.second(aGraph, partMaxSize);
            aGraph.saveToDot("/tmp/agraph-" + method.first + ".dot");
            std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/agraph-" << method.first << ".dot -o /tmp/agraph-" << method.first << ".pdf\n";

            Graph depGraph = aGraph.getPartitionGraph();
            assert(depGraph.isDag());
            std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
            std::cout << "Degree of parallelism after " << method.first << " partitioning : " << degPar.first << "  " << degPar.second << "\n";
            std::cout << "Number of partitions : " << depGraph.getNbNodes() << " -- avg part size : " << double(aGraph.getNbNodes())/double(depGraph.getNbNodes()) << "\n";

            depGraph.saveToDot("/tmp/depgraph-" + method.first + ".dot");
            std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/depgraph-" << method.first << ".dot -o /tmp/depgraph-" << method.first << ".pdf\n";

            double duration;
            std::vector<Executor::Event> events;
            std::tie(duration, events) = Executor::Execute(depGraph, nbThreads, 0.1, 0.2, 0.2);
            Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace-" + method.first + ".svg", events, nbThreads);
            std::cout << "=============================================================\n";
        }
    }

    return 0;
}
