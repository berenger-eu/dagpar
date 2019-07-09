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
                                    "[HELP] --strategy [none,stop,update]"
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

    const std::string strategy = params.getStr({"--strategy", "--s", "-strategy", "-s"});
    if(params.parseHasFailed()){
        std::cout << "[HELP] Invalid command at strategy.\n" << helpContent;
        return 1;
    }

    std::cout << "strategy = " << strategy << std::endl;

    if(strategy != "none" && strategy != "stop" && strategy != "update"){
        std::cout << "[HELP] Invalid strategy " << strategy << ".\n" << helpContent;
        return 1;
    }

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

    int bestGranularity = 1;
    double bestExecTime;
    double bestAvgClusterSize = 0;

    int nbNodes;

    {
        std::cout << "Original graph:\n";
        Graph aGraph(someDeps.first, someDeps.second);
        std::cout << " - Number of nodes : " << aGraph.getNbNodes() << "\n";
        nbNodes = aGraph.getNbNodes();
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

        double duration;
        std::vector<Executor::Event> events;
        std::tie(duration, events) = Executor::Execute(aGraph, nbThreads, overheadPerTaskOne, overheadPerPopOne, overheadPerPushOne);
        std::cout << " - Without clustering duration = " << duration << "\n";
        bestExecTime = duration;
    }

    //////////////////////////////////////////////////////////////////////////

    auto doItFunc = [&someDeps, &costs, exportDot, nbThreads, overheadPerTaskOne, overheadPerPopOne, overheadPerPushOne](const int idxGranularity, const std::string& methodName, auto method) -> std::tuple<double,double> {
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
            aGraph.saveToDot("/tmp/agraph-" + std::to_string(idxGranularity) + "-" + methodName + ".dot");
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
        //Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace-" + methodName + ".svg", events, nbThreads);
        std::cout << "=============================================================\n";
        return std::make_tuple(duration, double(aGraph.getNbNodes())/double(depGraph.getNbNodes()));
    };

    {
        int idxGranularity = 2;
        while(idxGranularity <= (bestGranularity+1)*2 && idxGranularity <= nbNodes){
            double execFinal;
            double avgClusterSize;
            if(strategy == "update"){
                std::tie(execFinal, avgClusterSize) = doItFunc(idxGranularity, "update", [](Graph& graph, const int clusterSize){
                    graph.GUpdate(clusterSize);
                 });
            }
            else if(strategy == "stop"){
                std::tie(execFinal, avgClusterSize) = doItFunc(idxGranularity, "stop", [](Graph& graph, const int clusterSize){
                    graph.GStop(clusterSize);
                 });
            }
            else{
                assert(strategy == "none");
                std::tie(execFinal, avgClusterSize) = doItFunc(idxGranularity, "none", [](Graph& graph, const int clusterSize){
                    graph.G(clusterSize);
                 });
            }

            if(execFinal < bestExecTime){
                bestExecTime = execFinal;
                bestGranularity = idxGranularity;
                bestAvgClusterSize = avgClusterSize;
            }

            idxGranularity += 1;
        }
    }

    std::cout << " - Best granularity : " << bestGranularity << std::endl;
    std::cout << " - Best duration : " << bestExecTime << std::endl;
    std::cout << " - Best avg cluster size : " << bestAvgClusterSize << std::endl;

    if(nbNodes/bestAvgClusterSize <= 8000){
        if(strategy == "update"){
            doItFunc(bestGranularity, "update-with-emulated-rafinement-", [overheadPerTaskOne, nbThreads, overheadPerPopOne, overheadPerPushOne](Graph& graph, const int clusterSize){
                graph.GUpdatePartitionWithEmulationRefinement(clusterSize, clusterSize*2, overheadPerTaskOne, nbThreads, overheadPerPopOne, overheadPerPushOne);
             });
        }
        else if(strategy == "stop"){
            doItFunc(bestGranularity, "stop-with-emulated-rafinement", [overheadPerTaskOne, nbThreads, overheadPerPopOne, overheadPerPushOne](Graph& graph, const int clusterSize){
                graph.GStopPartitionWithEmulationRefinement(clusterSize, clusterSize*2, overheadPerTaskOne, nbThreads, overheadPerPopOne, overheadPerPushOne);
             });
        }
    }

    return 0;
}
