#include "generator.hpp"
#include "graph.hpp"

#include "UTester.hpp"

#include <vector>
#include <functional>
#include <set>
#include <unordered_map>

class TestAll : public UTester< TestAll > {
    using Parent = UTester< TestAll >;

    static const int NbGenerators = 4;

    std::pair<int, std::vector<std::pair<int,int>>> GenerateDeps(const int choice, const int size){
        switch(choice){
        case 0:
            return Generator::GenerateDepTreeTasks(size);
            break;
        case 1:
            return Generator::GenerateDoubleDepTreeTasks(size);
            break;
        case 2:
            return Generator::Generate2DGrid(size);
            break;
        case 3:
            return Generator::GenerateDoubleDepTreeTasks(size);
            break;
        default:
            std::cout << "[ERROR] invalid choice\n";
            exit(1);
            break;
        }
        return std::pair<int, std::vector<std::pair<int,int>>>();
    }

    void TestDeps(const std::pair<int, std::vector<std::pair<int,int>>>& inDeps){
        const int nbNodes = inDeps.first;

        for(const auto& dep : inDeps.second){
            UASSERTETRUE(dep.first < nbNodes);
            UASSERTETRUE(dep.second < nbNodes);
        }
    }

    void TestGraph(const std::pair<int, std::vector<std::pair<int,int>>>& inDeps){
        Graph aGraph(inDeps.first, inDeps.second);
        UASSERTETRUE(aGraph.isDag());

        int nbDepsInGraph = 0;
        for(int idxNode = 0 ; idxNode < aGraph.getNbNodes() ; ++idxNode){
            const auto& node = aGraph.getNodeFromId(idxNode);
            UASSERTETRUE(node->getId() == idxNode);

            for(const auto& pred : node->getPredecessors()){
                UASSERTETRUE(std::find(pred->getSuccessors().begin(), pred->getSuccessors().end(),
                                      node) != pred->getSuccessors().end());
            }

            for(const auto& next : node->getSuccessors()){
                UASSERTETRUE(std::find(next->getPredecessors().begin(), next->getPredecessors().end(),
                                      node) != next->getPredecessors().end());
                nbDepsInGraph += 1;
            }
        }

        UASSERTETRUE(nbDepsInGraph == int(inDeps.second.size()));

        for(const auto& dep : inDeps.second){
            const auto srcNode = aGraph.getNodeFromId(dep.first);
            const auto dstNode = aGraph.getNodeFromId(dep.second);

            UASSERTETRUE(std::find(srcNode->getSuccessors().begin(), srcNode->getSuccessors().end(),
                                  dstNode) != srcNode->getSuccessors().end());
            UASSERTETRUE(std::find(dstNode->getPredecessors().begin(), dstNode->getPredecessors().end(),
                                  srcNode) != dstNode->getPredecessors().end());
        }
    }

    template <class FuncType>
    void TestClusteringGraph(const std::pair<int, std::vector<std::pair<int,int>>>& inDeps,
                             const int inClusterSize,
                             FuncType&& func){
        Graph aGraph(inDeps.first, inDeps.second);

        func(aGraph, inClusterSize);

        Graph depGraph = aGraph.getPartitionGraph();

        UASSERTETRUE(depGraph.isDag());

        {
            int nbDepsInGraph = 0;
            double totalCost = 0;
            for(int idxNode = 0 ; idxNode < depGraph.getNbNodes() ; ++idxNode){
                const auto& node = depGraph.getNodeFromId(idxNode);
                UASSERTETRUE(node->getId() == idxNode);

                totalCost += node->getCost();

                for(const auto& pred : node->getPredecessors()){
                    UASSERTETRUE(std::find(pred->getSuccessors().begin(), pred->getSuccessors().end(),
                                          node) != pred->getSuccessors().end());
                }

                for(const auto& next : node->getSuccessors()){
                    UASSERTETRUE(std::find(next->getPredecessors().begin(), next->getPredecessors().end(),
                                          node) != next->getPredecessors().end());
                    nbDepsInGraph += 1;
                }
            }

            UASSERTETRUE(std::abs(totalCost - aGraph.getNbNodes())/totalCost < 1E-14);
        }
        {
            std::unordered_map<int,int> costPerPartition;
            std::unordered_map<int,int> clusterIdPerNode;
            for(int idxNode = 0 ; idxNode < aGraph.getNbNodes() ; ++idxNode){
                const auto& node = aGraph.getNodeFromId(idxNode);
                clusterIdPerNode[node->getId()] = node->getPartitionId();
                costPerPartition[node->getPartitionId()] += 1;
            }

            for(const auto& dep : inDeps.second){
                const int idxSrcCluster = clusterIdPerNode[dep.first];
                const int idxDestCluster = clusterIdPerNode[dep.second];

                if(idxSrcCluster != idxDestCluster){
                    const auto& srcNode = depGraph.getNodeFromId(idxSrcCluster);
                    const auto& dstNode = depGraph.getNodeFromId(idxDestCluster);

                    UASSERTETRUE(std::find(srcNode->getSuccessors().begin(), srcNode->getSuccessors().end(),
                                          dstNode) != srcNode->getSuccessors().end());
                    UASSERTETRUE(std::find(dstNode->getPredecessors().begin(), dstNode->getPredecessors().end(),
                                          srcNode) != dstNode->getPredecessors().end());
                }
            }


            for(int idxNode = 0 ; idxNode < depGraph.getNbNodes() ; ++idxNode){
                const auto& node = depGraph.getNodeFromId(idxNode);
                UASSERTETRUE(node->getCost() == costPerPartition[node->getId()]);
            }
        }
        {
            std::unordered_map<int,std::vector<int>> depsPerSrcCluster;
            for(const auto& dep : inDeps.second){
                depsPerSrcCluster[aGraph.getNodeFromId(dep.first)->getPartitionId()].push_back(aGraph.getNodeFromId(dep.second)->getPartitionId());
            }

            for(int idxNode = 0 ; idxNode < depGraph.getNbNodes() ; ++idxNode){
                const auto& node = depGraph.getNodeFromId(idxNode);

                for(const auto& succ : node->getSuccessors()){
                    UASSERTETRUE(std::find(depsPerSrcCluster[node->getId()].begin(), depsPerSrcCluster[node->getId()].end(),
                                           succ->getId()) != depsPerSrcCluster[node->getId()].end());
                }
            }
        }


        if(inClusterSize == 1){
            UASSERTETRUE(aGraph.getNbNodes() == depGraph.getNbNodes());

            std::unordered_map<int,int> clusterIdPerNode;
            for(int idxNode = 0 ; idxNode < aGraph.getNbNodes() ; ++idxNode){
                const auto& node = aGraph.getNodeFromId(idxNode);
                const auto& clusterNode = depGraph.getNodeFromId(node->getPartitionId());

                UASSERTETRUE(node->getPredecessors().size() == clusterNode->getPredecessors().size());
                UASSERTETRUE(node->getSuccessors().size() == clusterNode->getSuccessors().size());
            }
        }
    }

    void TestConfig(){
        const std::vector<int> testSizes{4, 10, 40};
        const std::vector<int> clusterSizes{1, 10, 100};

        const std::vector<std::function<void(Graph&,int)>> allPartitionMethods= {
                                        [](Graph& graph, const int clusterSize){
                                            graph.G(clusterSize);
                                         },
                                        [](Graph& graph, const int clusterSize){
                                            graph.GUpdate(clusterSize);
                                         },
                                            [](Graph& graph, const int clusterSize){
                                                graph.GStop(clusterSize);
                                             },
                                        [](Graph& graph, const int clusterSize){
                                            graph.GStopPartitionWithEmulationRefinement(clusterSize, clusterSize, 0, 10, 0, 0);
                                         },
                                        [](Graph& graph, const int clusterSize){
                                            graph.GUpdatePartitionWithEmulationRefinement(clusterSize, clusterSize, 0, 10, 0, 0);
                                         }
                                        };

        for(const int testSize : testSizes){
            std::cout << "testSize = " << testSize << std::endl;
            for(int idxGenerator = 0 ; idxGenerator < NbGenerators ; ++idxGenerator){
                std::cout << "idxGenerator = " << idxGenerator << std::endl;
                const auto deps = GenerateDeps(idxGenerator, testSize);

                UASSERTETRUE(deps.second.size() != 0);

                TestDeps(deps);
                TestGraph(deps);

                for(const int clusterSize : clusterSizes){
                    std::cout << "clusterSize = " << clusterSize << std::endl;
                    for(const auto& partMethod : allPartitionMethods){
                        TestClusteringGraph(deps, clusterSize, partMethod);
                    }
                }
            }
        }

        const std::vector<std::string> testFilenames{"../data/chukrut_disque.txt",
                                                     "../data/chukrut_quadcurve.txt",
                                                     "../data/chukrut_quad.txt"};

        for(const auto& filename : testFilenames){
            std::cout << "filename = " << filename << std::endl;
            const auto deps = Generator::LoadEdgeFile(filename);

            UASSERTETRUE(deps.second.size() != 0);

            TestDeps(deps);
            TestGraph(deps);

            for(const int clusterSize : clusterSizes){
                std::cout << "clusterSize = " << clusterSize << std::endl;
                for(const auto& partMethod : allPartitionMethods){
                    TestClusteringGraph(deps, clusterSize, partMethod);
                }
            }
        }
    }

    void SetTests() {
        Parent::AddTest(&TestAll::TestConfig, "Test everything");
    }
};

// You must do this
TestClass(TestAll)

