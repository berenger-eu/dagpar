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
    const std::string helpContent = "[HELP] $ ./graphstat [options]\n"
                                    "[HELP] Where options are:\n"
                                    "[HELP] --help to get help\n"
                                    "[HELP] Or\n"
                                    "[HELP] --filename [filename]\n"
                                    "[HELP] Example:\n"
                                    "[HELP] ./graphstat --filename ../data/chukrut_disque.dot\n";

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

    std::cout << "filename = " << filename << std::endl;


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

        double totalcost = 0;
        double mincost = std::numeric_limits<double>::max();
        double maxcost = std::numeric_limits<double>::min();

        for(int idxNode = 0 ; idxNode < aGraph.getNbNodes() ; ++idxNode){
            auto node = aGraph.getNode(idxNode);
            totalcost += node->getCost();
            mincost = std::min(mincost, node->getCost());
            maxcost = std::max(maxcost, node->getCost());
        }

        std::cout << " - totalcost : " << totalcost << "\n";
        std::cout << " - mincost : " << mincost << "\n";
        std::cout << " - maxcost : " << maxcost << "\n";
        std::cout << " - avgcost : " << totalcost/double(aGraph.getNbNodes()) << "\n";

        int totaledges = 0;
        int maxedges = std::numeric_limits<int>::min();

        for(int idxNode = 0 ; idxNode < aGraph.getNbNodes() ; ++idxNode){
            auto node = aGraph.getNode(idxNode);
            totaledges += int(node->getSuccessors().size());
            maxedges = std::max(maxedges, int(node->getSuccessors().size()));
        }

        std::cout << " - totaledges : " << totaledges << "\n";
        std::cout << " - maxedges : " << maxedges << "\n";

        const auto distHist = aGraph.getDistHistogram();
        int maxHist = std::numeric_limits<int>::min();
        std::cout << " - dist histogram : ";
        for(const auto& dist : distHist){
            std::cout << dist << " ";
            maxHist = std::max(maxHist, dist);
        }
        std::cout << "\n";
        std::cout << " - maxHist : " << maxHist << "\n";
        std::cout << " - avgHist : " << double(aGraph.getNbNodes())/double(distHist.size()) << "\n";
    }

    return 0;
}
