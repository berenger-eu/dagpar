#include <iostream>
#include <functional>

#include "graph.hpp"
#include "executor.hpp"
#include "generator.hpp"
#include "utils.hpp"
#include "timer.hpp"

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

    const double overheadPerTask = params.getValue<double>({"--opt", "-opt"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at opt.\n" << helpContent;
        return 1;
    }

    const double overheadPerPush = params.getValue<double>({"--oppu", "-oppu"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at oppu.\n" << helpContent;
        return 1;
    }

    const double overheadPerPop = params.getValue<double>({"--oppo", "-oppo"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at oppo.\n" << helpContent;
        return 1;
    }

    std::cout << "filename = " << filename << std::endl;
    std::cout << "nbThreads = " << nbThreads << std::endl;
    std::cout << "overheadPerTask = " << overheadPerTask << std::endl;
    std::cout << "overheadPerPush = " << overheadPerPush << std::endl;
    std::cout << "overheadPerPop = " << overheadPerPop << std::endl;


    const bool exportDot = params.paramExist({"--export-dot", "-export-dot"});

    //////////////////////////////////////////////////////////////////////////

    std::pair<int, std::vector<std::pair<int,int>>> someDeps = Generator::LoadEdgeFile(filename);

    if(someDeps.second.size() == 0){
        std::cout << "[INFO] file is empty, exit...\n";
        return 1;
    }

    std::vector<double> costs = Generator::GetCostIfExist(filename);

    //////////////////////////////////////////////////////////////////////////

    double totalcost = 0;
    double overheadPerTaskOne = 0;
    double overheadPerPushOne = 0;
    double overheadPerPopOne = 0;

    bool isBig;
    int bestGranularity = 1;
    double bestExecTime;
    int nbNodes;

    {
        std::cout << "Original graph:\n";
        Graph aGraph(someDeps.first, someDeps.second);
        std::cout << " - Number of nodes : " << aGraph.getNbNodes() << "\n";
        nbNodes = aGraph.getNbNodes();
        assert(aGraph.isDag());
        isBig = (aGraph.getNbNodes() > 20000);
        std::pair<int,double> degGraph = aGraph.estimateDegreeOfParallelism();
        std::cout << " - Degree of parallelism one the sequential graph : " << degGraph.first << "  " << degGraph.second << "\n";

        if(costs.size()){
            for(int idxNode = 0 ; idxNode < aGraph.getNbNodes() ; ++idxNode){
                auto node = aGraph.getNode(idxNode);
                assert(node->getId() < int(costs.size()));
                node->setCost(costs[node->getId()]);
            }
        }

        for(int idxNode = 0 ; idxNode < aGraph.getNbNodes() ; ++idxNode){
            auto node = aGraph.getNode(idxNode);
            totalcost += node->getCost();
        }

        overheadPerTaskOne = (totalcost*overheadPerTask)/double(aGraph.getNbNodes());
        overheadPerPushOne = (totalcost*overheadPerPush)/double(aGraph.getNbNodes());
        overheadPerPopOne = (totalcost*overheadPerPop)/double(aGraph.getNbNodes());

        std::cout << "totalcost = " << totalcost << std::endl;
        std::cout << "overheadPerTaskOne = " << overheadPerTaskOne << std::endl;
        std::cout << "overheadPerPushOne = " << overheadPerPushOne << std::endl;
        std::cout << "overheadPerPopOne = " << overheadPerPopOne << std::endl;

        int duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(aGraph, nbThreads, overheadPerTaskOne, overheadPerPopOne, overheadPerPushOne);
        std::cout << " - Without clustering duration = " << duration << "\n";
        bestExecTime = duration;
    }

    //////////////////////////////////////////////////////////////////////////

    auto doItFunc = [&someDeps, &costs, exportDot, nbThreads, overheadPerTaskOne, overheadPerPopOne, overheadPerPushOne](const int idxGranularity, const std::string& methodName, auto method) -> double {
        std::cout << "=======================[" << methodName <<  "]======================================\n";

        std::cout << " - Granularity : " << idxGranularity << std::endl;

        Graph aGraph(someDeps.first, someDeps.second);
        std::cout << " - Number of nodes " << methodName << " : " << aGraph.getNbNodes() << "\n";
        assert(aGraph.isDag());
        std::pair<int,double> degGraph = aGraph.estimateDegreeOfParallelism();
        std::cout << " - Degree of parallelism one the " << methodName << " graph : " << degGraph.first << "  " << degGraph.second << "\n";

        if(costs.size()){
            for(int idxNode = 0 ; idxNode < aGraph.getNbNodes() ; ++idxNode){
                auto node = aGraph.getNode(idxNode);
                assert(node->getId() < int(costs.size()));
                node->setCost(costs[node->getId()]);
            }
        }

        Timer timer;
        method(aGraph, idxGranularity);
        timer.stop();

        if(exportDot){
            aGraph.saveToDot("/tmp/agraph-" + methodName + ".dot");
        }

        Graph depGraph = aGraph.getPartitionGraph();
        if(!depGraph.isDag()){
            std::cout << "[Error] not a DAG, will not exit...." << std::endl;
            exit(-1);
        }
        std::pair<int,double> degPar = depGraph.estimateDegreeOfParallelism();
        std::cout << " - Degree of parallelism after " << methodName << " partitioning : " << degPar.first << "  " << degPar.second << "\n";
        std::cout << " - Number of partitions " << methodName << " : " << depGraph.getNbNodes() << "\n";
        std::cout << " - Avg part size " << methodName << " : " << double(aGraph.getNbNodes())/double(depGraph.getNbNodes()) << "\n";
        std::cout << " - Time to partition with " << methodName << " : " << timer.getElapsed() << "\n";

        if(exportDot){
            depGraph.saveToDot("/tmp/depgraph-" + methodName + ".dot");
        }

        double duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(depGraph, nbThreads, overheadPerTaskOne, overheadPerPopOne, overheadPerPushOne);
        std::cout << " - with " << methodName << " clustering duration = " << duration << "\n";
        Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace-" + methodName + ".svg", events, nbThreads);
        std::cout << "=============================================================\n";
        return duration;
    };

    int bestH1 = 0;
    int bestH2 = 0;

    {
        int idxGranularity = 2;
        while(idxGranularity <= (bestGranularity+1)*2 && idxGranularity <= nbNodes){
            int h1 = std::max(1,int(idxGranularity*3./4.));
            int h2 = std::max(1,int(idxGranularity/4.));
            double execFinal = doItFunc(idxGranularity, "final", [h1,h2](Graph& graph, const int clusterSize){
                std::cout << " - h1 : " << h1 << std::endl;
                std::cout << " - h2 : " << h2 << std::endl;
                graph.partitionFinal(clusterSize, h1,h2);
             });

            if(execFinal < bestExecTime){
                bestExecTime = execFinal;
                bestGranularity = idxGranularity;
                bestH1 = h1;
                bestH2 = h2;
            }

            h1 = std::max(1,int((idxGranularity+1)/2));
            h2 = std::max(1,int((idxGranularity+1)/2));

            execFinal = doItFunc(idxGranularity, "final", [h1,h2](Graph& graph, const int clusterSize){
                std::cout << " - h1 : " << h1 << std::endl;
                std::cout << " - h2 : " << h2 << std::endl;
                graph.partitionFinal(clusterSize, h1,h2);
             });

            if(execFinal < bestExecTime){
                bestExecTime = execFinal;
                bestGranularity = idxGranularity;
                bestH1 = h1;
                bestH2 = h2;
            }

            h1 = std::max(1,int(idxGranularity/4.));
            h2 = std::max(1,int(idxGranularity*3./4.));

            execFinal = doItFunc(idxGranularity, "final", [h1,h2](Graph& graph, const int clusterSize){
                std::cout << " - h1 : " << h1 << std::endl;
                std::cout << " - h2 : " << h2 << std::endl;
                graph.partitionFinal(clusterSize, h1,h2);
             });

            if(execFinal < bestExecTime){
                bestExecTime = execFinal;
                bestGranularity = idxGranularity;
                bestH1 = h1;
                bestH2 = h2;
            }

            doItFunc(idxGranularity, "diamond", [h1,h2](Graph& graph, const int clusterSize){
                graph.partitionDiamond(clusterSize);
            });

            idxGranularity +=1;
        }
    }
//    {
//            int testGranularities[5];
//            double testExecutionTimes[5];
//        {
//            testGranularities[0] = 1;
//            const int idxGranularity = testGranularities[0];
//            int h1 = std::max(1,int(idxGranularity*3./4.));
//            int h2 = std::max(1,int(idxGranularity/4.));
//            testExecutionTimes[0] = doItFunc(idxGranularity, "final", [h1,h2](Graph& graph, const int clusterSize){
//                std::cout << " - h1 : " << h1 << std::endl;
//                std::cout << " - h2 : " << h2 << std::endl;
//                graph.partitionFinal(clusterSize, h1,h2);
//             });

//            if(testExecutionTimes[0] < bestExecTime){
//                bestExecTime = testExecutionTimes[0];
//                bestGranularity = testGranularities[0];
//                bestH1 = h1;
//                bestH2 = h2;
//            }
//        }
//        {
//            testGranularities[4] = nbNodes;
//            const int idxGranularity = testGranularities[4];
//            int h1 = std::max(1,int(idxGranularity*3./4.));
//            int h2 = std::max(1,int(idxGranularity/4.));
//            testExecutionTimes[4] = doItFunc(idxGranularity, "final", [h1,h2](Graph& graph, const int clusterSize){
//                std::cout << " - h1 : " << h1 << std::endl;
//                std::cout << " - h2 : " << h2 << std::endl;
//                graph.partitionFinal(clusterSize, h1,h2);
//             });

//            if(testExecutionTimes[4] < bestExecTime){
//                bestExecTime = testExecutionTimes[4];
//                bestGranularity = testGranularities[4];
//                bestH1 = h1;
//                bestH2 = h2;
//            }
//        }
//        {
//            testGranularities[2] = (nbNodes+1)/2;
//            const int idxGranularity = testGranularities[2];
//            int h1 = std::max(1,int(idxGranularity*3./4.));
//            int h2 = std::max(1,int(idxGranularity/4.));
//            testExecutionTimes[2] = doItFunc(idxGranularity, "final", [h1,h2](Graph& graph, const int clusterSize){
//                std::cout << " - h1 : " << h1 << std::endl;
//                std::cout << " - h2 : " << h2 << std::endl;
//                graph.partitionFinal(clusterSize, h1,h2);
//             });

//            if(testExecutionTimes[2] < bestExecTime){
//                bestExecTime = testExecutionTimes[2];
//                bestGranularity = testGranularities[2];
//                bestH1 = h1;
//                bestH2 = h2;
//            }
//        }
//        {
//            testGranularities[1] = (nbNodes+1)/4;
//            const int idxGranularity = testGranularities[1];
//            int h1 = std::max(1,int(idxGranularity*3./4.));
//            int h2 = std::max(1,int(idxGranularity/4.));
//            testExecutionTimes[1] = doItFunc(idxGranularity, "final", [h1,h2](Graph& graph, const int clusterSize){
//                std::cout << " - h1 : " << h1 << std::endl;
//                std::cout << " - h2 : " << h2 << std::endl;
//                graph.partitionFinal(clusterSize, h1,h2);
//             });

//            if(testExecutionTimes[1] < bestExecTime){
//                bestExecTime = testExecutionTimes[1];
//                bestGranularity = testGranularities[1];
//                bestH1 = h1;
//                bestH2 = h2;
//            }
//        }
//        {
//            testGranularities[3] = (3*(nbNodes+1))/4;
//            const int idxGranularity = testGranularities[3];
//            int h1 = std::max(1,int(idxGranularity*3./4.));
//            int h2 = std::max(1,int(idxGranularity/4.));
//            testExecutionTimes[3] = doItFunc(idxGranularity, "final", [h1,h2](Graph& graph, const int clusterSize){
//                std::cout << " - h1 : " << h1 << std::endl;
//                std::cout << " - h2 : " << h2 << std::endl;
//                graph.partitionFinal(clusterSize, h1,h2);
//             });

//            if(testExecutionTimes[3] < bestExecTime){
//                bestExecTime = testExecutionTimes[3];
//                bestGranularity = testGranularities[3];
//                bestH1 = h1;
//                bestH2 = h2;
//            }
//        }
//        while(testGranularities[0] + 4 <= testGranularities[4]){
//            std::cout << "[" << testGranularities[0] << "|" << testExecutionTimes[0] << "] "
//                                      << "[" << testGranularities[1] << "|" << testExecutionTimes[1] << "] "
//                                      << "[" << testGranularities[2] << "|" << testExecutionTimes[2] << "] "
//                                      << "[" << testGranularities[3] << "|" << testExecutionTimes[3] << "] "
//                                      << "[" << testGranularities[4] << "|" << testExecutionTimes[4] << "] " << std::endl;

//            if(testExecutionTimes[1] < testExecutionTimes[3]){
//                testGranularities[4] = testGranularities[2];
//                testExecutionTimes[4] = testExecutionTimes[2];

//                testGranularities[2] = testGranularities[1];
//                testExecutionTimes[2] = testExecutionTimes[1];
//            }
//            else{
//                testGranularities[0] = testGranularities[2];
//                testExecutionTimes[0] = testExecutionTimes[2];

//                testGranularities[2] = testGranularities[3];
//                testExecutionTimes[2] = testExecutionTimes[3];
//            }
//            {
//                testGranularities[1] = (testGranularities[2]-testGranularities[0]+1)/2 + testGranularities[0];
//                const int idxGranularity = testGranularities[1];
//                int h1 = std::max(1,int(idxGranularity*3./4.));
//                int h2 = std::max(1,int(idxGranularity/4.));
//                testExecutionTimes[1] = doItFunc(idxGranularity, "final", [h1,h2](Graph& graph, const int clusterSize){
//                    std::cout << " - h1 : " << h1 << std::endl;
//                    std::cout << " - h2 : " << h2 << std::endl;
//                    graph.partitionFinal(clusterSize, h1,h2);
//                 });

//                if(testExecutionTimes[1] < bestExecTime){
//                    bestExecTime = testExecutionTimes[1];
//                    bestGranularity = testGranularities[1];
//                    bestH1 = h1;
//                    bestH2 = h2;
//                }
//            }
//            {
//                testGranularities[3] = (testGranularities[4]-testGranularities[2]+1)/2 + testGranularities[2];
//                const int idxGranularity = testGranularities[3];
//                int h1 = std::max(1,int(idxGranularity*3./4.));
//                int h2 = std::max(1,int(idxGranularity/4.));
//                testExecutionTimes[3] = doItFunc(idxGranularity, "final", [h1,h2](Graph& graph, const int clusterSize){
//                    std::cout << " - h1 : " << h1 << std::endl;
//                    std::cout << " - h2 : " << h2 << std::endl;
//                    graph.partitionFinal(clusterSize, h1,h2);
//                 });

//                if(testExecutionTimes[3] < bestExecTime){
//                    bestExecTime = testExecutionTimes[3];
//                    bestGranularity = testGranularities[3];
//                    bestH1 = h1;
//                    bestH2 = h2;
//                }
//            }
//        }
//    }

    std::cout << " - Best granularity : " << bestGranularity << std::endl;
    std::cout << " - Best duration : " << bestExecTime << std::endl;
    std::cout << " - Best h1 : " << bestH1 << std::endl;
    std::cout << " - Best h2 : " << bestH2 << std::endl;

    doItFunc(bestGranularity, "final-with-neighbor-rafinement", [bestH1,bestH2](Graph& graph, const int clusterSize){
        std::cout << " - h1 : " << bestH1 << std::endl;
        std::cout << " - h2 : " << bestH2 << std::endl;
        graph.partitionFinalWithNeighborRefinement(clusterSize, bestH1, bestH2, clusterSize);
    });

    doItFunc(bestGranularity, "final-with-neighbor-rafinement-2", [bestH1,bestH2](Graph& graph, const int clusterSize){
        std::cout << " - h1 : " << std::max(1,bestH1/2) << std::endl;
        std::cout << " - h2 : " << bestH2 << std::endl;
        graph.partitionFinalWithNeighborRefinement(std::max(1,clusterSize/2), std::max(1,bestH1/2), bestH2, clusterSize);
    });

    if(!isBig){
        doItFunc(bestGranularity, "final-with-emulated-rafinement", [bestH1,bestH2, overheadPerTaskOne, nbThreads, overheadPerPopOne, overheadPerPushOne](Graph& graph, const int clusterSize){
            std::cout << " - h1 : " << bestH1 << std::endl;
            std::cout << " - h2 : " << bestH2 << std::endl;
            graph.partitionFinalWithEmulationRefinement(clusterSize, bestH1,bestH2, clusterSize, overheadPerTaskOne, nbThreads, overheadPerPopOne, overheadPerPushOne);
         });

        doItFunc(bestGranularity, "final-with-emulated-rafinement-2", [bestH1, bestH2, overheadPerTaskOne, nbThreads, overheadPerPopOne, overheadPerPushOne](Graph& graph, const int clusterSize){
            std::cout << " - h1 : " << std::max(1,bestH1/2) << std::endl;
            std::cout << " - h2 : " << bestH2 << std::endl;
            graph.partitionFinalWithEmulationRefinement(std::max(1,clusterSize/2), std::max(1, bestH1/2), bestH2, clusterSize, overheadPerTaskOne, nbThreads, overheadPerPopOne, overheadPerPushOne);
         });
    }


    return 0;
}
