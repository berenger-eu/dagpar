#include <iostream>
#include <functional>

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
    std::vector<double> costs;

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
            costs = Generator::GetCostIfExist(params[1]);
            if(someDeps.second.size() == 0){
                std::cout << "[INFO] file is empty, exit...\n";
                return 1;
            }
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

    const int nbThreads = 8;
    const int partMaxSize = 16;
    const int partMinSize = partMaxSize;
    std::cout << "nbThreads : " << nbThreads << " / partMinSize : " << partMinSize << " / partMaxSize : " << partMaxSize << "\n";

    const std::vector<std::pair<std::string,std::function<void(Graph&,int)>>> allPartitionMethods= {
                                    {"random", [](Graph& graph, const int clusterSize){
                                        graph.partitionRandom(clusterSize);
                                     }},
                                     {"greedy", [](Graph& graph, const int clusterSize){
                                         graph.partitionGreedy(clusterSize);
                                      }},
                                     {"backtrack", [](Graph& graph, const int clusterSize){
                                         graph.partitionBacktrack(clusterSize);
                                      }},
                                     {"advanced", [](Graph& graph, const int clusterSize){
                                         graph.partition(clusterSize, clusterSize);
                                      }},
                                    {"horizontal", [](Graph& graph, const int clusterSize){
                                        graph.partitionHorizontal(clusterSize);
                                     }},
                                    {"diamond", [](Graph& graph, const int clusterSize){
                                        graph.partitionDiamond(clusterSize);
                                     }},
                                    {"temporal", [](Graph& graph, const int clusterSize){
                                        graph.partitionTemporalPart(clusterSize);
                                     }},
                                    {"final", [](Graph& graph, const int clusterSize){
                                        graph.partitionFinal(clusterSize, 5, 1, false);
                                     }},
                                    {"final-with-rafinement", [](Graph& graph, const int clusterSize){
                                        graph.partitionFinal(clusterSize, 5, 1, true);
                                     }}
#ifdef USE_ACYCLIC
                                    ,
                                    {"acyclic", [](Graph& graph, const int clusterSize){
                                        graph.partitionAcyclic(clusterSize);
                                     }}
#endif
                                    };


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

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads, 0, 0.1, 0.2);
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace-" + method.first + ".svg", events, nbThreads);
        std::cout << "=============================================================\n";
    }
#ifdef USE_METIS
    {
        Graph aGraph(someDeps.first, someDeps.second);
        assert(aGraph.isDag());
        aGraph.partitionMetis(partMaxSize);
        aGraph.saveToDot("/tmp/agraph-metis.dot");
        std::cout << "Generate pdf of the graph with: dot -Tpdf /tmp/agraph-metis.dot -o /tmp/agraph-metis.pdf\n";
    }
#endif

    return 0;
}
