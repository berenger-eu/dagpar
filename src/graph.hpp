// Please read the corresponding licence file
#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <cassert>
#include <utility>
#include <vector>
#include <memory>
#include <deque>
#include <set>

// For the export to dot file
#include <string>
#include <fstream>
#include <unordered_map>
#include <array>
// For colors
#include <random>
// For warning messages
#include <iostream>
#include <queue>

#ifdef USE_METIS
#include <metis.h>
#endif

#include "node.hpp"

class Graph{    

    struct InfoPartition{
        int startingLevel;
        int limiteLevel;
        std::set<int> idsParentPartitionsSuccessors;
        std::set<int> idsParentPartitionsPredecessors;
        std::set<Node*> idsNodes;
    };

    std::vector<std::unique_ptr<Node>> nodesNotTopological;
    std::vector<Node*> nodes;

    static bool CoreIsDag(const std::vector<std::unique_ptr<Node>>& inNodes){
        std::deque<Node*> sources;
        std::set<Node*> alreadyVisistedNodes;
        for(auto& node : inNodes){
            if(node->getPredecessors().size() == 0){
                sources.push_back(node.get());
                alreadyVisistedNodes.insert(node.get());
            }
        }

        std::vector<int> counterRelease(inNodes.size(), 0);

        while(sources.size()){
            Node* selectedNode = sources.front();
            sources.pop_front();


            for(const auto& otherNode : selectedNode->getSuccessors()){
                counterRelease[otherNode->getId()] += 1;

                if(counterRelease[otherNode->getId()] > int(otherNode->getPredecessors().size())){
                    return false;
                }

                if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                    if(alreadyVisistedNodes.find(otherNode) != alreadyVisistedNodes.end()){
                        return false;
                    }
                    sources.push_back(otherNode);
                    alreadyVisistedNodes.insert(otherNode);
                }
            }
        }

        return alreadyVisistedNodes.size() == inNodes.size();
    }


    static double ExecPartions(const std::vector<InfoPartition>& partitions, const double overheadPerTask,
            const int nbWorkers, const int nbNodes, const double popOverhead, const double pushOverhead) {

         struct Worker{
             double scheduledAvailable;

             int currentPartitionId;

             int currentWorkerId;
             bool operator<(const Worker& other) const{
                 return scheduledAvailable > other.scheduledAvailable;
             }
         };

         std::vector<int> idleWorkerCount;

         idleWorkerCount.resize(nbWorkers);
         std::iota(idleWorkerCount.begin(), idleWorkerCount.end(), 0);

         std::vector<int> readyPartitions;
         double currentTime = 0;

         for(int idxPartition = 0 ; idxPartition < int(partitions.size()) ; ++idxPartition){
             if(partitions[idxPartition].idsNodes.size() && partitions[idxPartition].idsParentPartitionsPredecessors.size() == 0){
                 currentTime += pushOverhead;
                 readyPartitions.emplace_back(idxPartition);
             }
         }

         assert(readyPartitions.size());

         std::priority_queue<Worker> workers;
         int nbComputedTask = 0;

         while(readyPartitions.size() && idleWorkerCount.size()){
             const int readyPartitionId = readyPartitions.back();
             readyPartitions.pop_back();

             const int workerId = idleWorkerCount.back();
             idleWorkerCount.pop_back();

             double totalDuration = 0;
             for(const auto& node : partitions[readyPartitionId].idsNodes){
                 totalDuration += node->getCost();
                 // Not true since we did not update the nodes assert(node->getPartitionId() == readyPartitionId);
             }

             currentTime += popOverhead;

             Worker wk{currentTime + totalDuration + overheadPerTask,
                      readyPartitionId,
                      workerId};
             workers.push(wk);

             nbComputedTask += int(partitions[readyPartitionId].idsNodes.size());
         }

         assert(workers.size() != 0);

         std::vector<int> countPredecessorsOverPartitions(partitions.size(), 0);

         while(nbComputedTask != nbNodes){
             {
                 assert(workers.size());
                 Worker worker = workers.top();
                 workers.pop();

                 const int currentPartitionId = worker.currentPartitionId;

                 currentTime = std::max(currentTime, worker.scheduledAvailable);

                 // release dependencies
                 for(const auto& successorPartition : partitions[currentPartitionId].idsParentPartitionsSuccessors){
                     countPredecessorsOverPartitions[successorPartition] += 1;
                     if(countPredecessorsOverPartitions[successorPartition] == int(partitions[successorPartition].idsParentPartitionsPredecessors.size())){
                         readyPartitions.emplace_back(successorPartition);
                         currentTime+= pushOverhead;
                     }
                 }

                 // Make worker available again
                 idleWorkerCount.push_back(worker.currentWorkerId);
             }

             while(readyPartitions.size() && idleWorkerCount.size()){
                 const int readyPartitionId = readyPartitions.back();
                 readyPartitions.pop_back();

                 const int workerId = idleWorkerCount.back();
                 idleWorkerCount.pop_back();

                 double totalDuration = 0;
                 for(const auto& node : partitions[readyPartitionId].idsNodes){
                     totalDuration += node->getCost();
                     // Not true anymore since node has not been updated: assert(node->getPartitionId() == readyPartitionId);
                 }

                 currentTime += popOverhead;

                 Worker wk{currentTime + totalDuration + overheadPerTask,
                          readyPartitionId,
                          workerId};
                 workers.push(wk);

                 nbComputedTask += int(partitions[readyPartitionId].idsNodes.size());
             }

             assert(workers.size() != 0);
         }

         assert(workers.size() != 0);
         while(workers.size() != 0){
             assert(workers.size());
             Worker worker = workers.top();
             workers.pop();

             assert(currentTime <= worker.scheduledAvailable);
             currentTime = worker.scheduledAvailable;

             assert(partitions[worker.currentPartitionId].idsParentPartitionsSuccessors.size() == 0);
         }

         return currentTime;
    }


    std::vector<InfoPartition> partitionCore(const int maxSize,
                        const int inPhase1Height, const int inPhase2Height){
        std::vector<Node*> originalSources;
        for(auto& node : nodes){
            node->setPartitionId(-1);
            if(node->getPredecessors().size() == 0){
                originalSources.push_back(node);
            }
        }

        std::vector<int> maxDistFromTop(nodes.size(), -1);
        {
            std::vector<Node*> sources = originalSources;

            for(auto& node : sources){
                maxDistFromTop[node->getId()] = 0;
            }

            std::vector<int> counterRelease(nodes.size(), 0);
            while(sources.size()){
                Node* selectedNode = sources.back();
                sources.pop_back();

                // Add deps if released
                for(const auto& otherNode : selectedNode->getSuccessors()){
                    if(maxDistFromTop[otherNode->getId()] == -1){
                        maxDistFromTop[otherNode->getId()] = maxDistFromTop[selectedNode->getId()] + 1;
                    }
                    else{
                        maxDistFromTop[otherNode->getId()] = std::max(maxDistFromTop[selectedNode->getId()] + 1, maxDistFromTop[otherNode->getId()]);
                    }

                    counterRelease[otherNode->getId()] += 1;
                    assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                    if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                        sources.push_back(otherNode);
                    }
                }
            }
        }

        std::vector<std::pair<int,Node*>> maxDistFromTopWithNode(nodes.size());
        for(auto& node : nodes){
            maxDistFromTopWithNode[node->getId()].first = maxDistFromTop[node->getId()];
            maxDistFromTopWithNode[node->getId()].second = node;
        }

        std::sort(maxDistFromTopWithNode.begin(), maxDistFromTopWithNode.end(),
                  [](const std::pair<int,Node*>& p1, const std::pair<int,Node*>& p2){
            return p1.first < p2.first;
        });

        std::vector<std::pair<int,int>> consecuviteNodesSameDist;
        consecuviteNodesSameDist.reserve(maxDistFromTopWithNode.size()/1000);
        if(maxDistFromTopWithNode.size()){
            int indexStart = 0;
            int indexEnd = 1;
            while(indexEnd != int(maxDistFromTopWithNode.size())){
                if(maxDistFromTopWithNode[indexStart].first != maxDistFromTopWithNode[indexEnd].first){
                    consecuviteNodesSameDist.emplace_back(indexStart,indexEnd);
                    indexStart = indexEnd;
                }
                indexEnd += 1;
            }
            consecuviteNodesSameDist.emplace_back(indexStart,indexEnd);
        }

        std::vector<InfoPartition> proceedPartitionsInfo;
        proceedPartitionsInfo.reserve(nodes.size()/maxSize);

        for(const auto& consecutiveNodes : consecuviteNodesSameDist){
            std::vector<std::tuple<Node*,int,int>> nodeRelations(consecutiveNodes.second - consecutiveNodes.first);

            for(int idxNode = consecutiveNodes.first ; idxNode < consecutiveNodes.second ; ++idxNode){
                Node* selectedNode = maxDistFromTopWithNode[idxNode].second;

                const int nbPrevious = int(selectedNode->getPredecessors().size());
                std::set<int> predecessorPartitions;
                for(const auto& selectedNodeParent : selectedNode->getPredecessors()){
                    predecessorPartitions.insert(selectedNodeParent->getPartitionId());
                }
                const int nbPartitionPredecessors = int(predecessorPartitions.size());

                nodeRelations[idxNode-consecutiveNodes.first] = std::make_tuple(selectedNode, nbPrevious, nbPartitionPredecessors);
            }

            std::sort(nodeRelations.begin(), nodeRelations.end(),
                    [](const std::tuple<Node*,int,int>& n1, const std::tuple<Node*,int,int>& n2){
//                return std::get<2>(n1) > std::get<2>(n2)
//                        || (std::get<2>(n1) == std::get<2>(n2) && std::get<1>(n1) > std::get<1>(n2));
                return std::get<1>(n1) > std::get<1>(n2)
                        || (std::get<1>(n1) == std::get<1>(n2) && std::get<2>(n1) < std::get<2>(n2));
            });


            for(auto& iter : nodeRelations){
                Node* selectedNode = std::get<0>(iter);
                const int selectedNodeDistFromTop = maxDistFromTop[selectedNode->getId()];

                std::set<int> selectedNodeParentPartitionIds;
                for(const auto& selectedNodeParent : selectedNode->getPredecessors()){
                    assert(selectedNodeParent->getPartitionId() != -1);
                    selectedNodeParentPartitionIds.insert(selectedNodeParent->getPartitionId());
                }

                if(selectedNodeParentPartitionIds.size() == 0){
                    const int currentPartitionId = int(proceedPartitionsInfo.size());
                    selectedNode->setPartitionId(currentPartitionId);

                    proceedPartitionsInfo.resize(proceedPartitionsInfo.size() + 1);
                    proceedPartitionsInfo.back().startingLevel = selectedNodeDistFromTop;
                    proceedPartitionsInfo.back().limiteLevel = selectedNodeDistFromTop + 1;
                    proceedPartitionsInfo.back().idsNodes.insert(selectedNode);
                }
                else if(selectedNodeParentPartitionIds.size() == 1){
                    const int uniqueParentPartitionId = (*selectedNodeParentPartitionIds.begin());
                    assert(uniqueParentPartitionId < int(proceedPartitionsInfo.size()));
                    const auto& uniqueParentPartitionInfo = proceedPartitionsInfo[uniqueParentPartitionId];

                    if(selectedNodeDistFromTop <= uniqueParentPartitionInfo.startingLevel + inPhase1Height + inPhase2Height
                            && int(uniqueParentPartitionInfo.idsNodes.size()) < maxSize){
                        selectedNode->setPartitionId(uniqueParentPartitionId);
                        assert(proceedPartitionsInfo[uniqueParentPartitionId].limiteLevel <= selectedNodeDistFromTop+1);
                        proceedPartitionsInfo[uniqueParentPartitionId].limiteLevel = selectedNodeDistFromTop + 1;
                        proceedPartitionsInfo[uniqueParentPartitionId].idsNodes.insert(selectedNode);
                    }
                    else{
                        const int currentPartitionId = int(proceedPartitionsInfo.size());
                        selectedNode->setPartitionId(currentPartitionId);

                        proceedPartitionsInfo.resize(proceedPartitionsInfo.size() + 1);
                        proceedPartitionsInfo.back().startingLevel = selectedNodeDistFromTop;
                        proceedPartitionsInfo.back().limiteLevel = selectedNodeDistFromTop + 1;
                        proceedPartitionsInfo.back().idsParentPartitionsPredecessors.insert(uniqueParentPartitionId);
                        proceedPartitionsInfo.back().idsNodes.insert(selectedNode);
                        proceedPartitionsInfo[uniqueParentPartitionId].idsParentPartitionsSuccessors.insert(currentPartitionId);
                    }
                }
                else{
                    std::vector<std::tuple<int,int,int>> possibleCurrentPartitionsWithDistAndNb;
                    for(const auto& parentPartitionId : selectedNodeParentPartitionIds){
                        const auto& parentPartitionInfo = proceedPartitionsInfo[parentPartitionId];
                        if(selectedNodeDistFromTop <= parentPartitionInfo.startingLevel + inPhase1Height
                                && int(parentPartitionInfo.idsNodes.size()) < maxSize){
                            possibleCurrentPartitionsWithDistAndNb.emplace_back(parentPartitionId,
                                                                                parentPartitionInfo.startingLevel,
                                                                                int(parentPartitionInfo.idsNodes.size()));
                        }
                    }

                    std::sort(possibleCurrentPartitionsWithDistAndNb.begin(), possibleCurrentPartitionsWithDistAndNb.end(),
                              [](const std::tuple<int,int,int>& p1, const std::tuple<int,int,int>& p2){
                        return std::get<1>(p1) > std::get<1>(p2)
                                || (std::get<1>(p1) == std::get<1>(p2)
                                    && std::get<2>(p1) < std::get<2>(p2));
                    });

                    bool nodeHasBeenInserted = false;

                    for(const auto& testPartitionIdWithDist : possibleCurrentPartitionsWithDistAndNb){
                        const int testPartitionId = std::get<0>(testPartitionIdWithDist);
                        bool testPartitionIdIsLinkedToAParent = false;

                        std::deque<int> children;
                        std::set<int> childrenProceed;
                        for(const auto& idxChild : proceedPartitionsInfo[testPartitionId].idsParentPartitionsSuccessors){
                            if(selectedNodeParentPartitionIds.find(idxChild) != selectedNodeParentPartitionIds.end()){
                                testPartitionIdIsLinkedToAParent = true;
                                break;
                            }
                            if(childrenProceed.find(idxChild) == childrenProceed.end()){
                                children.push_back(idxChild);
                                childrenProceed.insert(idxChild);
                            }
                        }

                        while(testPartitionIdIsLinkedToAParent == false && children.size()){
                            const int testPartitionIter = children.front();
                            children.pop_front();

                            for(const auto& idxChild : proceedPartitionsInfo[testPartitionIter].idsParentPartitionsSuccessors){
                                if(selectedNodeParentPartitionIds.find(idxChild) != selectedNodeParentPartitionIds.end()){
                                    testPartitionIdIsLinkedToAParent = true;
                                    break;
                                }
                                if(childrenProceed.find(idxChild) == childrenProceed.end()){
                                    children.push_back(idxChild);
                                    childrenProceed.insert(idxChild);
                                }
                            }
                        }

                        if(testPartitionIdIsLinkedToAParent == false){
                            selectedNode->setPartitionId(testPartitionId);
                            assert(proceedPartitionsInfo[testPartitionId].limiteLevel <= selectedNodeDistFromTop+1);
                            proceedPartitionsInfo[testPartitionId].limiteLevel = selectedNodeDistFromTop + 1;
                            proceedPartitionsInfo[testPartitionId].idsNodes.insert(selectedNode);

                            selectedNodeParentPartitionIds.erase(testPartitionId);
                            proceedPartitionsInfo[testPartitionId].idsParentPartitionsPredecessors.insert(selectedNodeParentPartitionIds.begin(),
                                                                                                 selectedNodeParentPartitionIds.end());

                            for(const auto& newParent : selectedNodeParentPartitionIds){
                                proceedPartitionsInfo[newParent].idsParentPartitionsSuccessors.insert(testPartitionId);
                            }

                            nodeHasBeenInserted = true;
                            break;
                        }
                    }

                    if(nodeHasBeenInserted == false){
                        const int currentPartitionId = int(proceedPartitionsInfo.size());
                        selectedNode->setPartitionId(currentPartitionId);

                        proceedPartitionsInfo.resize(proceedPartitionsInfo.size() + 1);
                        proceedPartitionsInfo.back().startingLevel = selectedNodeDistFromTop;
                        proceedPartitionsInfo.back().limiteLevel = selectedNodeDistFromTop + 1;
                        proceedPartitionsInfo.back().idsNodes.insert(selectedNode);
                        proceedPartitionsInfo.back().idsParentPartitionsPredecessors = selectedNodeParentPartitionIds;

                        for(const auto& parentPartitionId : selectedNodeParentPartitionIds){
                            proceedPartitionsInfo[parentPartitionId].idsParentPartitionsSuccessors.insert(currentPartitionId);
                        }
                    }
                }
            }
        }
        return proceedPartitionsInfo;
    }


public:
    explicit Graph(const int inNodes, const std::vector<std::pair<int,int>>& inDependencyList){
        nodesNotTopological.resize(inNodes);

        for(int idxNode = 0 ; idxNode < int(nodesNotTopological.size()) ; ++idxNode){
            nodesNotTopological[idxNode].reset(new Node(idxNode));
        }

        for(const auto& dep : inDependencyList){
            //Is not true when index are not in order: assert(dep.first < dep.second);
            assert(dep.second < inNodes);
            nodesNotTopological[dep.first]->addSuccessor(nodesNotTopological[dep.second].get());
            nodesNotTopological[dep.second]->addPredecessor(nodesNotTopological[dep.first].get());
        }

        if(CoreIsDag(nodesNotTopological)){
            std::vector<Node*> sources;
            for(auto& node : nodesNotTopological){
                node->setPartitionId(-1);
                if(node->getPredecessors().size() == 0){
                    sources.push_back(node.get());
                }
            }
            std::vector<int> counterRelease(nodesNotTopological.size(), 0);
            nodes.reserve(inNodes);
            while(sources.size()){
                Node* selectedNode = sources.back();
                sources.pop_back();

                nodes.push_back(selectedNode);

                for(const auto& otherNode : selectedNode->getSuccessors()){
                    counterRelease[otherNode->getId()] += 1;
                    assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                    if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                        sources.push_back(otherNode);
                    }
                }
            }
        }
        else{
            nodes.reserve(nodesNotTopological.size());
            for(auto& node : nodesNotTopological){
                nodes.push_back(node.get());
            }
        }
    }

    explicit Graph(const int inDependencyList[], const int inNbDep){
        for(int idxDep = 0 ; idxDep < inNbDep ; ++idxDep){
            const int depSrc = inDependencyList[idxDep*2+0];
            const int depDest = inDependencyList[idxDep*2+1];

            if(int(nodesNotTopological.size()) <= std::max(depSrc,depDest)){
                const int currentNbNodes = std::max(depSrc,depDest) + 1;
                const int pastNbNodes = int(nodesNotTopological.size());
                nodesNotTopological.resize(currentNbNodes);
                for(int idxNode = pastNbNodes ; idxNode < currentNbNodes ; ++idxNode){
                    nodesNotTopological[idxNode].reset(new Node(idxNode));
                }
            }

            nodesNotTopological[depSrc]->addSuccessor(nodesNotTopological[depDest].get());
            nodesNotTopological[depDest]->addPredecessor(nodesNotTopological[depSrc].get());
        }

        if(CoreIsDag(nodesNotTopological)){
            std::vector<Node*> sources;
            for(auto& node : nodesNotTopological){
                node->setPartitionId(-1);
                if(node->getPredecessors().size() == 0){
                    sources.push_back(node.get());
                }
            }

            std::vector<int> counterRelease(nodesNotTopological.size(), 0);
            nodes.reserve(nodesNotTopological.size());
            while(sources.size()){
                Node* selectedNode = sources.back();
                sources.pop_back();

                nodes.push_back(selectedNode);

                for(const auto& otherNode : selectedNode->getSuccessors()){
                    counterRelease[otherNode->getId()] += 1;
                    assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                    if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                        sources.push_back(otherNode);
                    }
                }
            }
        }
        else{
            nodes.reserve(nodesNotTopological.size());
            for(auto& node : nodesNotTopological){
                nodes.push_back(node.get());
            }
        }
    }


    Graph(const int inNodes,
          const std::vector<std::pair<int,int>>& inDependencyList,
          const std::vector<int>& inPartitionPerNode,
          const std::vector<double>& inCostPerNode) : Graph(inNodes, inDependencyList){
        assert(nodes.size() == inPartitionPerNode.size());
        assert(nodes.size() == inCostPerNode.size());

        for(int idxNode = 0 ; idxNode < int(nodesNotTopological.size()) ; ++idxNode){
            nodesNotTopological[idxNode]->setPartitionId(inPartitionPerNode[idxNode]);
            nodesNotTopological[idxNode]->setCost(inCostPerNode[idxNode]);
        }
    }

    Graph(const Graph&) = delete;
    Graph& operator=(const Graph&) = delete;

    Graph(Graph&&) = default;
    Graph& operator=(Graph&&) = default;

    int getNbNodes() const{
        return int(nodes.size());
    }

    const Node* getNode(const int inIdxNode) const{
        assert(inIdxNode < int(nodes.size()));
        return nodes[inIdxNode];
    }

    const Node* getNodeFromId(const int inIdNode) const{
        assert(inIdNode < int(nodesNotTopological.size()));
        return nodesNotTopological[inIdNode].get();
    }

    Node* getNode(const int inIdxNode) {
        assert(inIdxNode < int(nodes.size()));
        return nodes[inIdxNode];
    }

    Node* getNodeFromId(const int inIdNode) {
        assert(inIdNode < int(nodesNotTopological.size()));
        return nodesNotTopological[inIdNode].get();
    }

    void saveToDot(const std::string& inFilename){
        std::unordered_map<int, std::array<double, 3>> partitionColors;

        std::ofstream dotFile(inFilename);
        dotFile << "digraph G {\n";

        for(const auto& node : nodes){
            for(const auto& otherNode : node->getSuccessors()){
                dotFile << node->getId() << " -> " << otherNode->getId() << ";\n";
            }

            if(partitionColors.find(node->getPartitionId()) == partitionColors.end()){
                std::mt19937 rng(node->getPartitionId());
                std::uniform_real_distribution<double> uni(0,1);
                partitionColors[node->getPartitionId()] = std::array<double, 3>{uni(rng), uni(rng), uni(rng)};
            }

            const auto& color = partitionColors[node->getPartitionId()];
            dotFile << node->getId() << " [label=\"" << node->getId() << " [" << node->getCost() << "/" << node->getPartitionId() << "]\", style=filled,color=\"" << color[0] << " " << color[1] << " " << color[2] << "\"]\n";
        }


        dotFile << "}\n";
        dotFile.close();
    }

    void partitionDiamond(const int maxSize, const bool warnIfInvalid = false){
        // This algorithm will work only if:
        // - there is one root
        // - each node has at most 3 pred/succ dependencies
        std::vector<Node*> originalSources;
        for(auto& node : nodes){
            node->setPartitionId(-1);
            if(node->getPredecessors().size() == 0){
                originalSources.push_back(node);
            }
            if(warnIfInvalid){
                const bool depsAreValid = (int(node->getPredecessors().size() + node->getSuccessors().size()) <= 4)
                        && int(node->getPredecessors().size()) >= 0 && int(node->getSuccessors().size()) >= 0
                        && int(node->getPredecessors().size()) < 4 && int(node->getSuccessors().size()) < 4;
                if(!depsAreValid){
                    std::cerr << "[GRAPH][ERROR] The number of dependencies are invalid for node " << node->getId() << "\n";
                    std::cerr << "[GRAPH][ERROR] - nb predecessors " << node->getPredecessors().size() << "\n";
                    std::cerr << "[GRAPH][ERROR] - nb successors " << node->getSuccessors().size() << "\n";
                    std::cerr << "[GRAPH][ERROR] - the code will continue anyway..." << std::endl;
                }
            }
        }

        std::vector<int> maxDistFromTop(nodes.size(), -1);
        {
            std::vector<Node*> sources = originalSources;

            for(auto& node : sources){
                maxDistFromTop[node->getId()] = 0;
            }

            std::vector<int> counterRelease(nodes.size(), 0);
            while(sources.size()){
                Node* selectedNode = sources.back();
                sources.pop_back();

                // Add deps if released
                for(const auto& otherNode : selectedNode->getSuccessors()){
                    if(maxDistFromTop[otherNode->getId()] == -1){
                        maxDistFromTop[otherNode->getId()] = maxDistFromTop[selectedNode->getId()] + 1;
                    }
                    else{
                        maxDistFromTop[otherNode->getId()] = std::max(maxDistFromTop[selectedNode->getId()] + 1, maxDistFromTop[otherNode->getId()]);
                    }

                    counterRelease[otherNode->getId()] += 1;
                    assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                    if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                        sources.push_back(otherNode);
                    }
                }
            }
        }

        std::vector<std::pair<int,Node*>> maxDistFromTopWithNode(nodes.size());
        for(auto& node : nodes){
            maxDistFromTopWithNode[node->getId()].first = maxDistFromTop[node->getId()];
            maxDistFromTopWithNode[node->getId()].second = node;
        }

        std::sort(maxDistFromTopWithNode.begin(), maxDistFromTopWithNode.end(),
                  [](const std::pair<int,Node*>& p1, const std::pair<int,Node*>& p2){
            return p1.first < p2.first;
        });

        int nbLevelsInHalfDiamond = 0;
        while(((nbLevelsInHalfDiamond+2)*(nbLevelsInHalfDiamond+1))/2 < maxSize/2){
            nbLevelsInHalfDiamond += 1;
        }

        struct InfoPartition{
            int nbNnodesInPartition;
            int startingLevel;
            std::set<int> idsParentPartitionsSuccessors;
            std::set<int> idsParentPartitionsPredecessors;
        };

        std::vector<InfoPartition> proceedPartitionsInfo;
        proceedPartitionsInfo.reserve(nodes.size()/maxSize);

        for(int idxNode = 0 ; idxNode < int(maxDistFromTopWithNode.size()) ; ++idxNode){
            Node* selectedNode = maxDistFromTopWithNode[idxNode].second;
            const int selectedNodeDistFromTop = maxDistFromTop[selectedNode->getId()];

            std::set<int> selectedNodeParentPartitionIds;
            for(const auto& selectedNodeParent : selectedNode->getPredecessors()){
                assert(selectedNodeParent->getPartitionId() != -1);
                selectedNodeParentPartitionIds.insert(selectedNodeParent->getPartitionId());
            }

            if(selectedNodeParentPartitionIds.size() == 0){
                const int currentPartitionId = int(proceedPartitionsInfo.size());
                selectedNode->setPartitionId(currentPartitionId);

                proceedPartitionsInfo.resize(proceedPartitionsInfo.size() + 1);
                proceedPartitionsInfo.back().startingLevel = selectedNodeDistFromTop;
                proceedPartitionsInfo.back().nbNnodesInPartition = 1;
            }
            else if(selectedNodeParentPartitionIds.size() == 1){
                const int uniqueParentPartitionId = (*selectedNodeParentPartitionIds.begin());
                assert(uniqueParentPartitionId < int(proceedPartitionsInfo.size()));
                const auto& uniqueParentPartitionInfo = proceedPartitionsInfo[uniqueParentPartitionId];

                if(selectedNodeDistFromTop < uniqueParentPartitionInfo.startingLevel + 2 * nbLevelsInHalfDiamond + 1
                        && uniqueParentPartitionInfo.nbNnodesInPartition < maxSize){
                    selectedNode->setPartitionId(uniqueParentPartitionId);
                    proceedPartitionsInfo[uniqueParentPartitionId].nbNnodesInPartition += 1;
                }
                else{
                    const int currentPartitionId = int(proceedPartitionsInfo.size());
                    selectedNode->setPartitionId(currentPartitionId);

                    proceedPartitionsInfo.resize(proceedPartitionsInfo.size() + 1);
                    proceedPartitionsInfo.back().startingLevel = selectedNodeDistFromTop;
                    proceedPartitionsInfo.back().nbNnodesInPartition = 1;
                    proceedPartitionsInfo.back().idsParentPartitionsPredecessors.insert(uniqueParentPartitionId);
                    proceedPartitionsInfo[uniqueParentPartitionId].idsParentPartitionsSuccessors.insert(currentPartitionId);
                }
            }
            else{
                std::vector<std::tuple<int,int,int>> possibleCurrentPartitionsWithDistAndNb;
                for(const auto& parentPartitionId : selectedNodeParentPartitionIds){
                    const auto& parentPartitionInfo = proceedPartitionsInfo[parentPartitionId];
                    if(selectedNodeDistFromTop < parentPartitionInfo.startingLevel + nbLevelsInHalfDiamond + 1
                            && parentPartitionInfo.nbNnodesInPartition < maxSize){
                        possibleCurrentPartitionsWithDistAndNb.emplace_back(parentPartitionId,
                                                                            parentPartitionInfo.startingLevel,
                                                                            parentPartitionInfo.nbNnodesInPartition);
                    }
                }

                std::sort(possibleCurrentPartitionsWithDistAndNb.begin(), possibleCurrentPartitionsWithDistAndNb.end(),
                          [](const std::tuple<int,int,int>& p1, const std::tuple<int,int,int>& p2){
                    return std::get<1>(p1) > std::get<1>(p2)
                            || (std::get<1>(p1) == std::get<1>(p2)
                                && std::get<2>(p1) < std::get<2>(p2));
                });

                bool nodeHasBeenInserted = false;

                for(const auto& testPartitionIdWithDist : possibleCurrentPartitionsWithDistAndNb){
                    const int testPartitionId = std::get<0>(testPartitionIdWithDist);
                    bool testPartitionIdIsLinkedToAParent = false;

                    std::deque<int> children;
                    std::set<int> childrenProceed;
                    for(const auto& idxChild : proceedPartitionsInfo[testPartitionId].idsParentPartitionsSuccessors){
                        if(selectedNodeParentPartitionIds.find(idxChild) != selectedNodeParentPartitionIds.end()){
                            testPartitionIdIsLinkedToAParent = true;
                            break;
                        }
                        if(childrenProceed.find(idxChild) == childrenProceed.end()){
                            children.push_back(idxChild);
                            childrenProceed.insert(idxChild);
                        }
                    }

                    while(testPartitionIdIsLinkedToAParent == false && children.size()){
                        const int testPartitionIter = children.front();
                        children.pop_front();

                        for(const auto& idxChild : proceedPartitionsInfo[testPartitionIter].idsParentPartitionsSuccessors){
                            if(selectedNodeParentPartitionIds.find(idxChild) != selectedNodeParentPartitionIds.end()){
                                testPartitionIdIsLinkedToAParent = true;
                                break;
                            }
                            if(childrenProceed.find(idxChild) == childrenProceed.end()){
                                children.push_back(idxChild);
                                childrenProceed.insert(idxChild);
                            }
                        }
                    }

                    if(testPartitionIdIsLinkedToAParent == false){
                        selectedNode->setPartitionId(testPartitionId);
                        proceedPartitionsInfo[testPartitionId].nbNnodesInPartition += 1;

                        selectedNodeParentPartitionIds.erase(testPartitionId);
                        proceedPartitionsInfo[testPartitionId].idsParentPartitionsPredecessors.insert(selectedNodeParentPartitionIds.begin(),
                                                                                             selectedNodeParentPartitionIds.end());

                        for(const auto& newParent : selectedNodeParentPartitionIds){
                            proceedPartitionsInfo[newParent].idsParentPartitionsSuccessors.insert(testPartitionId);
                        }

                        nodeHasBeenInserted = true;
                        break;
                    }
                }

                if(nodeHasBeenInserted == false){
                    const int currentPartitionId = int(proceedPartitionsInfo.size());
                    selectedNode->setPartitionId(currentPartitionId);

                    proceedPartitionsInfo.resize(proceedPartitionsInfo.size() + 1);
                    proceedPartitionsInfo.back().startingLevel = selectedNodeDistFromTop;
                    proceedPartitionsInfo.back().nbNnodesInPartition = 1;
                    proceedPartitionsInfo.back().idsParentPartitionsPredecessors = selectedNodeParentPartitionIds;

                    for(const auto& parentPartitionId : selectedNodeParentPartitionIds){
                        proceedPartitionsInfo[parentPartitionId].idsParentPartitionsSuccessors.insert(currentPartitionId);
                    }
                }
            }
        }
    }

    std::pair<int,double> estimateDegreeOfParallelism() const{
        int maxSourcesSize = 0;
        int sumSourcesSize = 0;

        std::deque<Node*> sources;
        for(auto& node : nodes){
            if(node->getPredecessors().size() == 0){
                sources.push_back(node);
            }
        }

        std::vector<int> counterRelease(nodes.size(), 0);

        while(sources.size()){
            maxSourcesSize = std::max(maxSourcesSize, int(sources.size()));
            sumSourcesSize += int(sources.size());

            Node* selectedNode = sources.front();
            sources.pop_front();
            for(const auto& otherNode : selectedNode->getSuccessors()){
                counterRelease[otherNode->getId()] += 1;
                assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                    sources.push_back(otherNode);
                }
            }
        }

        return std::pair<int,double>(maxSourcesSize, double(sumSourcesSize)/double(nodes.size()));
    }

    bool isDag() const {
        return CoreIsDag(nodesNotTopological);
    }

    Graph getPartitionGraph() const{
        std::vector<std::pair<int,int>> dependencyBetweenPartitions;
        std::vector<double> partitionCosts;

        for(const auto& node : nodes){
            for(const auto& otherNode : node->getSuccessors()){
                if(node->getPartitionId() != otherNode->getPartitionId()){
                    dependencyBetweenPartitions.emplace_back(std::pair<int,int>{node->getPartitionId(), otherNode->getPartitionId()});
                }
            }

            assert(node->getPartitionId() != -1);

            if(int(partitionCosts.size()) <= node->getPartitionId()){
                partitionCosts.resize(node->getPartitionId()+1);
            }
            partitionCosts[node->getPartitionId()] += node->getCost();
        }

        std::sort(dependencyBetweenPartitions.begin(), dependencyBetweenPartitions.end(),
                  [](const std::pair<int,int>& dep1, const std::pair<int,int>& dep2){
            return dep1.first < dep2.first || (dep1.first == dep2.first && dep1.second < dep2.second);
        });

        auto last = std::unique(dependencyBetweenPartitions.begin(), dependencyBetweenPartitions.end());
        dependencyBetweenPartitions.erase(last, dependencyBetweenPartitions.end());

        std::vector<int> partitionIds(partitionCosts.size());
        std::iota(partitionIds.begin(), partitionIds.end(), 0);

        return Graph(int(partitionCosts.size()), dependencyBetweenPartitions, partitionIds, partitionCosts);
    }

    bool isPartitionGraphDag() const{
        return getPartitionGraph().isDag();
    }

    void partitionFinal(const int maxSize,
                        const int inPhase1Height, const int inPhase2Height){
        partitionCore(maxSize, inPhase1Height, inPhase2Height);
    }

    void partitionFinalWithNeighborRefinement(const int maxSize,
                        const int inPhase1Height, const int inPhase2Height,
                                              const int maxSizeAfterRefinement){
        std::vector<InfoPartition> proceedPartitionsInfo = partitionCore(maxSize, inPhase1Height, inPhase2Height);

        bool hasChanged = true;
        while(hasChanged){
            hasChanged = false;

            for(int idxPart = 0 ; idxPart < int(proceedPartitionsInfo.size()) ; ++idxPart){
                auto& part = proceedPartitionsInfo[idxPart];
                if(part.idsNodes.size() && int(part.idsNodes.size()) < maxSizeAfterRefinement){
                    std::set<int> possibleOtherParts;

                    for(const Node* selectedNode : part.idsNodes){
                        for(const auto& otherNode : selectedNode->getSuccessors()){

                            if(otherNode->getPartitionId() != idxPart // TODO
                                    && int(proceedPartitionsInfo[otherNode->getPartitionId()].idsNodes.size() + part.idsNodes.size()) <= maxSizeAfterRefinement){
                                possibleOtherParts.insert(otherNode->getPartitionId());
                            }

                            for(const auto& otherOtherNode : otherNode->getPredecessors()){
                                if(otherOtherNode->getPartitionId() != idxPart
                                        && int(proceedPartitionsInfo[otherOtherNode->getPartitionId()].idsNodes.size() + part.idsNodes.size()) <= maxSizeAfterRefinement){
                                    possibleOtherParts.insert(otherOtherNode->getPartitionId());
                                }
                            }
                        }

                        for(const auto& otherNode : selectedNode->getPredecessors()){

                            if(otherNode->getPartitionId() != idxPart // TODO
                                    && int(proceedPartitionsInfo[otherNode->getPartitionId()].idsNodes.size() + part.idsNodes.size()) <= maxSizeAfterRefinement){
                                possibleOtherParts.insert(otherNode->getPartitionId());
                            }

                            for(const auto& otherOtherNode : otherNode->getSuccessors()){
                                if(otherOtherNode->getPartitionId() != idxPart
                                        && int(proceedPartitionsInfo[otherOtherNode->getPartitionId()].idsNodes.size() + part.idsNodes.size()) <= maxSize*2){
                                    possibleOtherParts.insert(otherOtherNode->getPartitionId());
                                }
                            }
                        }
                    }

                    std::vector<std::pair<int,int>> potentialPerParts;
                    potentialPerParts.reserve(possibleOtherParts.size());
                    for(auto idxOtherPart : possibleOtherParts){
                        std::set<int> intersectionPredecessors;
                        std::set_intersection(part.idsParentPartitionsPredecessors.begin(),
                                              part.idsParentPartitionsPredecessors.end(),
                                              proceedPartitionsInfo[idxOtherPart].idsParentPartitionsPredecessors.begin(),
                                              proceedPartitionsInfo[idxOtherPart].idsParentPartitionsPredecessors.end(),
                                          std::inserter(intersectionPredecessors,intersectionPredecessors.begin()));

                        std::set<int> intersectionSuccessors;
                        std::set_intersection(part.idsParentPartitionsSuccessors.begin(),
                                              part.idsParentPartitionsSuccessors.end(),
                                              proceedPartitionsInfo[idxOtherPart].idsParentPartitionsSuccessors.begin(),
                                              proceedPartitionsInfo[idxOtherPart].idsParentPartitionsSuccessors.end(),
                                          std::inserter(intersectionSuccessors,intersectionSuccessors.begin()));

                        const int score = int(intersectionPredecessors.size() + intersectionSuccessors.size());
                        potentialPerParts.emplace_back(idxOtherPart, score);
                    }

                    std::sort(potentialPerParts.begin(), potentialPerParts.end(),
                              [](const std::pair<int,int>& p1, const std::pair<int,int>& p2){
                        return p1.second > p2.second;
                    });

                    for(int idxPotential = 0 ; idxPotential < int(potentialPerParts.size()) ; ++idxPotential){
                        const int idxOtherPart = potentialPerParts[idxPotential].first;

                        bool testPartitionIdIsLinkedToAParent = false;

                        for(const auto& pair : {std::pair<int,int>{idxPart, idxOtherPart}, std::pair<int,int>{idxOtherPart, idxPart}}){
                            const int src = pair.first;
                            const int dest = pair.second;

                            std::deque<int> children;
                            std::set<int> childrenProceed;
                            for(const auto& idxChild : proceedPartitionsInfo[src].idsParentPartitionsSuccessors){
                                if(idxChild != dest && childrenProceed.find(idxChild) == childrenProceed.end()){
                                    children.push_back(idxChild);
                                    childrenProceed.insert(idxChild);
                                }
                            }

                            while(testPartitionIdIsLinkedToAParent == false && children.size()){
                                const int testPartitionIter = children.front();
                                children.pop_front();

                                for(const auto& idxChild : proceedPartitionsInfo[testPartitionIter].idsParentPartitionsSuccessors){
                                    if(idxChild == dest){
                                        assert(proceedPartitionsInfo[testPartitionIter].startingLevel < proceedPartitionsInfo[dest].limiteLevel);
                                        testPartitionIdIsLinkedToAParent = true;
                                        break;
                                    }
                                    if(proceedPartitionsInfo[idxChild].startingLevel < proceedPartitionsInfo[dest].limiteLevel
                                             && childrenProceed.find(idxChild) == childrenProceed.end()){
                                        children.push_back(idxChild);
                                        childrenProceed.insert(idxChild);
                                    }
                                }
                            }
                        }

                        if(testPartitionIdIsLinkedToAParent == false){
                            for(auto idxNext : proceedPartitionsInfo[idxOtherPart].idsParentPartitionsSuccessors){
                                proceedPartitionsInfo[idxNext].idsParentPartitionsPredecessors.insert(idxPart);
                                proceedPartitionsInfo[idxNext].idsParentPartitionsPredecessors.erase(idxOtherPart);
                            }
                            for(auto idxPrev : proceedPartitionsInfo[idxOtherPart].idsParentPartitionsPredecessors){
                                proceedPartitionsInfo[idxPrev].idsParentPartitionsSuccessors.insert(idxPart);
                                proceedPartitionsInfo[idxPrev].idsParentPartitionsSuccessors.erase(idxOtherPart);
                            }

                            if(part.startingLevel != proceedPartitionsInfo[idxOtherPart].startingLevel){
                                std::deque<int> parents;
                                std::set<int> parentsProceed;
                                if(part.startingLevel > proceedPartitionsInfo[idxOtherPart].startingLevel){
                                    part.startingLevel = proceedPartitionsInfo[idxOtherPart].startingLevel;
                                    for(const auto& idxParent : proceedPartitionsInfo[idxPart].idsParentPartitionsPredecessors){
                                        if(part.startingLevel < proceedPartitionsInfo[idxParent].startingLevel && parentsProceed.find(idxParent) == parentsProceed.end()){
                                            parents.push_back(idxParent);
                                            parentsProceed.insert(idxParent);
                                        }
                                    }
                                }
                                else{
                                    for(const auto& idxParent : proceedPartitionsInfo[idxOtherPart].idsParentPartitionsPredecessors){
                                        if(part.startingLevel < proceedPartitionsInfo[idxParent].startingLevel && parentsProceed.find(idxParent) == parentsProceed.end()){
                                            parents.push_back(idxParent);
                                            parentsProceed.insert(idxParent);
                                        }
                                    }
                                }

                                while(parents.size()){
                                    const int testPartitionIter = parents.front();
                                    parents.pop_front();

                                    proceedPartitionsInfo[testPartitionIter].startingLevel = part.startingLevel;

                                    for(const auto& idxParent : proceedPartitionsInfo[testPartitionIter].idsParentPartitionsPredecessors){
                                        if(part.startingLevel < proceedPartitionsInfo[idxParent].startingLevel && parentsProceed.find(idxParent) == parentsProceed.end()){
                                            parents.push_back(idxParent);
                                            parentsProceed.insert(idxParent);
                                        }
                                    }
                                }
                            }

                            part.limiteLevel = std::max(part.limiteLevel, proceedPartitionsInfo[idxOtherPart].limiteLevel);
                            part.idsParentPartitionsSuccessors.insert(proceedPartitionsInfo[idxOtherPart].idsParentPartitionsSuccessors.begin(),
                                                                      proceedPartitionsInfo[idxOtherPart].idsParentPartitionsSuccessors.end());
                            part.idsParentPartitionsSuccessors.erase(idxPart);
                            part.idsParentPartitionsSuccessors.erase(idxOtherPart);
                            part.idsParentPartitionsPredecessors.insert(proceedPartitionsInfo[idxOtherPart].idsParentPartitionsPredecessors.begin(),
                                                                        proceedPartitionsInfo[idxOtherPart].idsParentPartitionsPredecessors.end());
                            part.idsParentPartitionsPredecessors.erase(idxPart);
                            part.idsParentPartitionsPredecessors.erase(idxOtherPart);
                            part.idsNodes.insert(proceedPartitionsInfo[idxOtherPart].idsNodes.begin(),
                                                 proceedPartitionsInfo[idxOtherPart].idsNodes.end());

                            for(Node* otherNode : proceedPartitionsInfo[idxOtherPart].idsNodes){
                                otherNode->setPartitionId(idxPart);
                            }

                            proceedPartitionsInfo[idxOtherPart].idsParentPartitionsSuccessors.clear();
                            proceedPartitionsInfo[idxOtherPart].idsParentPartitionsPredecessors.clear();
                            proceedPartitionsInfo[idxOtherPart].idsNodes.clear();

                            hasChanged = true;
                            break;
                        }
                    }
                }
            }
        }

        // Reset partition number
        std::vector<int> partitionsPermut(proceedPartitionsInfo.size(), -1);
        int permutPartitionCounter = 0;
        for(Node* node : nodes){
            if(partitionsPermut[node->getPartitionId()] == -1){
                partitionsPermut[node->getPartitionId()] = (permutPartitionCounter++);
            }
            node->setPartitionId(partitionsPermut[node->getPartitionId()]);
        }
    }

    void partitionFinalWithEmulationRefinement(const int maxSize,
                        const int inPhase1Height, const int inPhase2Height,
                        const int maxSizeAfterRefinement, const double inOverheadPerTask,
                        const int inNbWorkers, const double inPopOverhead, const double inPushOverhead){
        std::vector<InfoPartition> proceedPartitionsInfo = partitionCore(maxSize, inPhase1Height, inPhase2Height);

        double currentExecutionTime = ExecPartions(proceedPartitionsInfo, inOverheadPerTask, inNbWorkers, int(nodes.size()),
                                                   inPopOverhead, inPushOverhead);

        std::vector<InfoPartition> copyProceedPartitionsInfo = proceedPartitionsInfo;

        bool hasChanged = true;
        while(hasChanged){
            hasChanged = false;

            for(int idxPart = 0 ; idxPart < int(proceedPartitionsInfo.size()) ; ++idxPart){
                auto& part = proceedPartitionsInfo[idxPart];
                if(part.idsNodes.size() && int(part.idsNodes.size()) < maxSizeAfterRefinement){
                    std::set<int> possibleOtherParts;

                    for(const Node* selectedNode : part.idsNodes){
                        for(const auto& otherNode : selectedNode->getSuccessors()){

                            if(otherNode->getPartitionId() != idxPart // TODO
                                    && int(proceedPartitionsInfo[otherNode->getPartitionId()].idsNodes.size() + part.idsNodes.size()) <= maxSizeAfterRefinement){
                                possibleOtherParts.insert(otherNode->getPartitionId());
                            }

                            for(const auto& otherOtherNode : otherNode->getPredecessors()){
                                if(otherOtherNode->getPartitionId() != idxPart
                                        && int(proceedPartitionsInfo[otherOtherNode->getPartitionId()].idsNodes.size() + part.idsNodes.size()) <= maxSizeAfterRefinement){
                                    possibleOtherParts.insert(otherOtherNode->getPartitionId());
                                }
                            }
                        }

                        for(const auto& otherNode : selectedNode->getPredecessors()){

                            if(otherNode->getPartitionId() != idxPart // TODO
                                    && int(proceedPartitionsInfo[otherNode->getPartitionId()].idsNodes.size() + part.idsNodes.size()) <= maxSizeAfterRefinement){
                                possibleOtherParts.insert(otherNode->getPartitionId());
                            }

                            for(const auto& otherOtherNode : otherNode->getSuccessors()){
                                if(otherOtherNode->getPartitionId() != idxPart
                                        && int(proceedPartitionsInfo[otherOtherNode->getPartitionId()].idsNodes.size() + part.idsNodes.size()) <= maxSizeAfterRefinement){
                                    possibleOtherParts.insert(otherOtherNode->getPartitionId());
                                }
                            }
                        }
                    }

                    int bestPartitionId = -1;
                    double bestExecTime = currentExecutionTime;

                    for(const auto& parentPartitionId : possibleOtherParts){
                        bool testPartitionIdIsLinkedToAParent = false;

                        for(const auto& pair : {std::pair<int,int>{idxPart, parentPartitionId}, std::pair<int,int>{parentPartitionId, idxPart}}){
                            const int src = pair.first;
                            const int dest = pair.second;

                            std::deque<int> children;
                            std::set<int> childrenProceed;
                            for(const auto& idxChild : proceedPartitionsInfo[src].idsParentPartitionsSuccessors){
                                if(idxChild != dest && childrenProceed.find(idxChild) == childrenProceed.end()){
                                    children.push_back(idxChild);
                                    childrenProceed.insert(idxChild);
                                }
                            }

                            while(testPartitionIdIsLinkedToAParent == false && children.size()){
                                const int testPartitionIter = children.front();
                                children.pop_front();

                                for(const auto& idxChild : proceedPartitionsInfo[testPartitionIter].idsParentPartitionsSuccessors){
                                    if(idxChild == dest){
                                        assert(proceedPartitionsInfo[testPartitionIter].startingLevel < proceedPartitionsInfo[dest].limiteLevel);
                                        testPartitionIdIsLinkedToAParent = true;
                                        break;
                                    }
                                    if(proceedPartitionsInfo[idxChild].startingLevel < proceedPartitionsInfo[dest].limiteLevel
                                             && childrenProceed.find(idxChild) == childrenProceed.end()){
                                        children.push_back(idxChild);
                                        childrenProceed.insert(idxChild);
                                    }
                                }
                            }
                        }

                        if(testPartitionIdIsLinkedToAParent == false){
                            for(auto idxNext : copyProceedPartitionsInfo[parentPartitionId].idsParentPartitionsSuccessors){
                                copyProceedPartitionsInfo[idxNext].idsParentPartitionsPredecessors.insert(idxPart);
                                copyProceedPartitionsInfo[idxNext].idsParentPartitionsPredecessors.erase(parentPartitionId);
                            }
                            for(auto idxPrev : copyProceedPartitionsInfo[parentPartitionId].idsParentPartitionsPredecessors){
                                copyProceedPartitionsInfo[idxPrev].idsParentPartitionsSuccessors.insert(idxPart);
                                copyProceedPartitionsInfo[idxPrev].idsParentPartitionsSuccessors.erase(parentPartitionId);
                            }

                            auto& copyPart = copyProceedPartitionsInfo[idxPart];

                            copyPart.limiteLevel = std::max(copyPart.limiteLevel, copyProceedPartitionsInfo[parentPartitionId].limiteLevel);
                            copyPart.idsParentPartitionsSuccessors.insert(copyProceedPartitionsInfo[parentPartitionId].idsParentPartitionsSuccessors.begin(),
                                                                      copyProceedPartitionsInfo[parentPartitionId].idsParentPartitionsSuccessors.end());
                            copyPart.idsParentPartitionsSuccessors.erase(idxPart);
                            copyPart.idsParentPartitionsSuccessors.erase(parentPartitionId);
                            copyPart.idsParentPartitionsPredecessors.insert(copyProceedPartitionsInfo[parentPartitionId].idsParentPartitionsPredecessors.begin(),
                                                                        copyProceedPartitionsInfo[parentPartitionId].idsParentPartitionsPredecessors.end());
                            copyPart.idsParentPartitionsPredecessors.erase(idxPart);
                            copyPart.idsParentPartitionsPredecessors.erase(parentPartitionId);
                            copyPart.idsNodes.insert(copyProceedPartitionsInfo[parentPartitionId].idsNodes.begin(),
                                                 copyProceedPartitionsInfo[parentPartitionId].idsNodes.end());

                            copyProceedPartitionsInfo[parentPartitionId].idsParentPartitionsSuccessors.clear();
                            copyProceedPartitionsInfo[parentPartitionId].idsParentPartitionsPredecessors.clear();
                            copyProceedPartitionsInfo[parentPartitionId].idsNodes.clear();


                            double testExecutionTime = ExecPartions(copyProceedPartitionsInfo, inOverheadPerTask, inNbWorkers, int(nodes.size()),
                                                                    inPopOverhead, inPushOverhead);

                            if(testExecutionTime < bestExecTime){
                                bestPartitionId = parentPartitionId;
                                bestExecTime = bestExecTime;
                            }

                            copyProceedPartitionsInfo[idxPart] = proceedPartitionsInfo[idxPart];
                            copyProceedPartitionsInfo[parentPartitionId] = proceedPartitionsInfo[parentPartitionId];

                            for(auto idxNext : copyProceedPartitionsInfo[parentPartitionId].idsParentPartitionsSuccessors){
                                copyProceedPartitionsInfo[idxNext] = proceedPartitionsInfo[idxNext];
                            }
                            for(auto idxPrev : copyProceedPartitionsInfo[parentPartitionId].idsParentPartitionsPredecessors){
                                copyProceedPartitionsInfo[idxPrev] = proceedPartitionsInfo[idxPrev];
                            }
                        }
                    }

                    if(bestPartitionId != -1){
                        for(auto idxNext : proceedPartitionsInfo[bestPartitionId].idsParentPartitionsSuccessors){
                            proceedPartitionsInfo[idxNext].idsParentPartitionsPredecessors.insert(idxPart);
                            proceedPartitionsInfo[idxNext].idsParentPartitionsPredecessors.erase(bestPartitionId);
                        }
                        for(auto idxPrev : proceedPartitionsInfo[bestPartitionId].idsParentPartitionsPredecessors){
                            proceedPartitionsInfo[idxPrev].idsParentPartitionsSuccessors.insert(idxPart);
                            proceedPartitionsInfo[idxPrev].idsParentPartitionsSuccessors.erase(bestPartitionId);
                        }

                        if(part.startingLevel != proceedPartitionsInfo[bestPartitionId].startingLevel){
                            std::deque<int> parents;
                            std::set<int> parentsProceed;
                            if(part.startingLevel > proceedPartitionsInfo[bestPartitionId].startingLevel){
                                part.startingLevel = proceedPartitionsInfo[bestPartitionId].startingLevel;
                                for(const auto& idxParent : proceedPartitionsInfo[idxPart].idsParentPartitionsPredecessors){
                                    if(part.startingLevel < proceedPartitionsInfo[idxParent].startingLevel && parentsProceed.find(idxParent) == parentsProceed.end()){
                                        parents.push_back(idxParent);
                                        parentsProceed.insert(idxParent);
                                    }
                                }
                            }
                            else{
                                for(const auto& idxParent : proceedPartitionsInfo[bestPartitionId].idsParentPartitionsPredecessors){
                                    if(part.startingLevel < proceedPartitionsInfo[idxParent].startingLevel && parentsProceed.find(idxParent) == parentsProceed.end()){
                                        parents.push_back(idxParent);
                                        parentsProceed.insert(idxParent);
                                    }
                                }
                            }

                            while(parents.size()){
                                const int testPartitionIter = parents.front();
                                parents.pop_front();

                                proceedPartitionsInfo[testPartitionIter].startingLevel = part.startingLevel;
                                copyProceedPartitionsInfo[testPartitionIter].startingLevel = part.startingLevel;

                                for(const auto& idxParent : proceedPartitionsInfo[testPartitionIter].idsParentPartitionsPredecessors){
                                    if(part.startingLevel < proceedPartitionsInfo[idxParent].startingLevel && parentsProceed.find(idxParent) == parentsProceed.end()){
                                        parents.push_back(idxParent);
                                        parentsProceed.insert(idxParent);
                                    }
                                }
                            }
                        }

                        part.limiteLevel = std::max(part.limiteLevel, proceedPartitionsInfo[bestPartitionId].limiteLevel);
                        part.idsParentPartitionsSuccessors.insert(proceedPartitionsInfo[bestPartitionId].idsParentPartitionsSuccessors.begin(),
                                                                  proceedPartitionsInfo[bestPartitionId].idsParentPartitionsSuccessors.end());
                        part.idsParentPartitionsSuccessors.erase(idxPart);
                        part.idsParentPartitionsSuccessors.erase(bestPartitionId);
                        part.idsParentPartitionsPredecessors.insert(proceedPartitionsInfo[bestPartitionId].idsParentPartitionsPredecessors.begin(),
                                                                    proceedPartitionsInfo[bestPartitionId].idsParentPartitionsPredecessors.end());
                        part.idsParentPartitionsPredecessors.erase(idxPart);
                        part.idsParentPartitionsPredecessors.erase(bestPartitionId);
                        part.idsNodes.insert(proceedPartitionsInfo[bestPartitionId].idsNodes.begin(),
                                             proceedPartitionsInfo[bestPartitionId].idsNodes.end());

                        for(Node* otherNode : proceedPartitionsInfo[bestPartitionId].idsNodes){
                            otherNode->setPartitionId(idxPart);
                        }

                        proceedPartitionsInfo[bestPartitionId].idsParentPartitionsSuccessors.clear();
                        proceedPartitionsInfo[bestPartitionId].idsParentPartitionsPredecessors.clear();
                        proceedPartitionsInfo[bestPartitionId].idsNodes.clear();


                        for(auto idxNext : copyProceedPartitionsInfo[bestPartitionId].idsParentPartitionsSuccessors){
                            copyProceedPartitionsInfo[idxNext].idsParentPartitionsPredecessors.insert(idxPart);
                            copyProceedPartitionsInfo[idxNext].idsParentPartitionsPredecessors.erase(bestPartitionId);
                        }
                        for(auto idxPrev : copyProceedPartitionsInfo[bestPartitionId].idsParentPartitionsPredecessors){
                            copyProceedPartitionsInfo[idxPrev].idsParentPartitionsSuccessors.insert(idxPart);
                            copyProceedPartitionsInfo[idxPrev].idsParentPartitionsSuccessors.erase(bestPartitionId);
                        }

                        copyProceedPartitionsInfo[idxPart] = proceedPartitionsInfo[idxPart];
                        copyProceedPartitionsInfo[bestPartitionId] = proceedPartitionsInfo[bestPartitionId];

                        hasChanged = true;
                        // Do not restart from here break;
                    }
                }
            }
        }

        // Reset partition number
        std::vector<int> partitionsPermut(proceedPartitionsInfo.size(), -1);
        int permutPartitionCounter = 0;
        for(Node* node : nodes){
            if(partitionsPermut[node->getPartitionId()] == -1){
                partitionsPermut[node->getPartitionId()] = (permutPartitionCounter++);
            }
            node->setPartitionId(partitionsPermut[node->getPartitionId()]);
        }
    }

    std::vector<int> getDistHistogram() const {
        std::vector<Node*> originalSources;
        for(auto& node : nodes){
            node->setPartitionId(-1);
            if(node->getPredecessors().size() == 0){
                originalSources.push_back(node);
            }
        }

        std::vector<int> maxDistFromTop(nodes.size(), -1);
        {
            std::vector<Node*> sources = originalSources;

            for(auto& node : sources){
                maxDistFromTop[node->getId()] = 0;
            }

            std::vector<int> counterRelease(nodes.size(), 0);
            while(sources.size()){
                Node* selectedNode = sources.back();
                sources.pop_back();

                // Add deps if released
                for(const auto& otherNode : selectedNode->getSuccessors()){
                    if(maxDistFromTop[otherNode->getId()] == -1){
                        maxDistFromTop[otherNode->getId()] = maxDistFromTop[selectedNode->getId()] + 1;
                    }
                    else{
                        maxDistFromTop[otherNode->getId()] = std::max(maxDistFromTop[selectedNode->getId()] + 1, maxDistFromTop[otherNode->getId()]);
                    }

                    counterRelease[otherNode->getId()] += 1;
                    assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                    if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                        sources.push_back(otherNode);
                    }
                }
            }
        }

        std::vector<int> distHist;
        for(const int dist : maxDistFromTop){
            if(int(distHist.size()) <= dist){
                distHist.resize(dist+1, 0);
            }
            distHist[dist] += 1;
        }

        return distHist;
    }
};

#endif
