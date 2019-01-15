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

#include "node.hpp"

class Graph{
    std::vector<std::unique_ptr<Node>> nodes;
public:
    explicit Graph(const int inNodes, const std::vector<std::pair<int,int>>& inDependencyList){
        nodes.resize(inNodes);
        for(int idxNode = 0 ; idxNode < int(nodes.size()) ; ++idxNode){
            nodes[idxNode].reset(new Node(idxNode));
        }

        for(const auto& dep : inDependencyList){
            assert(dep.first < dep.second);
            assert(dep.second < inNodes);
            nodes[dep.first]->addSuccessor(nodes[dep.second].get());
            nodes[dep.second]->addPredecessor(nodes[dep.first].get());
        }
    }


    Graph(const int inNodes,
          const std::vector<std::pair<int,int>>& inDependencyList,
          const std::vector<int>& inPartitionPerNode,
          const std::vector<int>& inCostPerNode) : Graph(inNodes, inDependencyList){
        assert(nodes.size() == inPartitionPerNode.size());
        assert(nodes.size() == inCostPerNode.size());

        for(int idxNode = 0 ; idxNode < int(nodes.size()) ; ++idxNode){
            nodes[idxNode]->setPartitionId(inPartitionPerNode[idxNode]);
            nodes[idxNode]->setCost(inCostPerNode[idxNode]);
        }
    }

    int getNbNodes() const{
        return int(nodes.size());
    }

    const Node* getNode(const int inIdxNode) const{
        assert(inIdxNode < int(nodes.size()));
        return nodes[inIdxNode].get();
    }

    void saveToDot(const std::string& inFilename){
        std::unordered_map<int, std::array<double, 3>> partitionColors;
        std::mt19937 rng(0);
        std::uniform_real_distribution<double> uni(0,1);

        std::ofstream dotFile(inFilename);
        dotFile << "digraph G {\n";

        for(const auto& node : nodes){
            for(const auto& otherNode : node->getSuccessors()){
                dotFile << node->getId() << " -> " << otherNode->getId() << ";\n";
            }

            if(partitionColors.find(node->getPartitionId()) == partitionColors.end()){
                partitionColors[node->getPartitionId()] = std::array<double, 3>{uni(rng), uni(rng), uni(rng)};
            }

            const auto& color = partitionColors[node->getPartitionId()];
            dotFile << node->getId() << " [label=\"" << node->getId() << "|" << node->getCost() << "\", style=filled,color=\"" << color[0] << " " << color[1] << " " << color[2] << "\"]\n";
        }


        dotFile << "}\n";
        dotFile.close();
    }

    void partitionRandom(const int maxSize){
        std::vector<Node*> sources;
        for(auto& node : nodes){
            if(node->getPredecessors().size() == 0){
                sources.push_back(node.get());
            }
        }

        std::vector<int> counterRelease(nodes.size(), 0);

        std::vector<Node*> orderedNodes;
        orderedNodes.reserve(nodes.size());

        std::mt19937 rng(0);
        while(sources.size()){
            std::uniform_int_distribution<int> uni(0,int(sources.size())-1);
            const int selectedNodePos = uni(rng);
            // Select a node
            Node* selectedNode = sources[selectedNodePos];
            orderedNodes.push_back(selectedNode);
            // Remove it from source
            sources[selectedNodePos] = sources[sources.size()-1];
            sources.pop_back();
            // Add deps if released
            for(const auto& otherNode : selectedNode->getSuccessors()){
                counterRelease[otherNode->getId()] += 1;
                assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                    sources.push_back(otherNode);
                }
            }
        }


        int currentPartitionId = 0;
        int currentPartitionSize = 0;
        for(auto& node : orderedNodes){
            node->setPartitionId(currentPartitionId);
            currentPartitionSize += 1;
            if(currentPartitionSize == maxSize){
                currentPartitionId += 1;
                currentPartitionSize = 0;
            }
        }
    }

    void partitionGreedy(const int maxSize, const bool useIntermediateStep = false){
        std::vector<int> minDistFromTop(nodes.size(), -1);
        std::vector<int> maxDistFromTop(nodes.size(), -1);

        std::vector<Node*> originalSources;
        for(auto& node : nodes){
            if(node->getPredecessors().size() == 0){
                originalSources.push_back(node.get());
                minDistFromTop[node->getId()] = 0;
                maxDistFromTop[node->getId()] = 0;
            }
        }
        {
            std::vector<Node*> sources = originalSources;
            std::vector<int> counterRelease(nodes.size(), 0);
            while(sources.size()){
                Node* selectedNode = sources.back();
                sources.pop_back();

                // Add deps if released
                for(const auto& otherNode : selectedNode->getSuccessors()){
                    if(minDistFromTop[otherNode->getId()] == -1){
                        minDistFromTop[otherNode->getId()] = minDistFromTop[selectedNode->getId()] + 1;
                    }
                    else{
                        minDistFromTop[otherNode->getId()] = std::min(minDistFromTop[selectedNode->getId()] + 1, minDistFromTop[otherNode->getId()]);
                    }

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

        {
            std::vector<Node*> sources = originalSources;
            std::vector<int> counterRelease(nodes.size(), 0);

            int currentPartitionId = 0;

            while(sources.size()){
                int idxSelectNode = 0;
                for(int idxNode = 1 ; idxNode < int(sources.size()) ; ++idxNode){
                    if(minDistFromTop[sources[idxNode]->getId()] < minDistFromTop[sources[idxSelectNode]->getId()]){
                        idxSelectNode = idxNode;
                    }
                }

                std::deque<Node*> partitionNodes;
                partitionNodes.push_front(sources[idxSelectNode]);
                std::swap(sources[idxSelectNode], sources[sources.size()-1]);
                sources.pop_back();

                int currentPartitionSize = 0;
                while(partitionNodes.size() && currentPartitionSize < maxSize){
                    Node* selectedNode = partitionNodes.front();
                    partitionNodes.pop_front();

                    selectedNode->setPartitionId(currentPartitionId);
                    currentPartitionSize += 1;

                    std::deque<Node*> nextNodes;
                    for(const auto& otherNode : selectedNode->getSuccessors()){
                        counterRelease[otherNode->getId()] += 1;
                        assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                        if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                            nextNodes.push_back(otherNode);
                        }
                    }

                    std::sort(nextNodes.begin(), nextNodes.end(), [&](const Node* n1, const Node* n2){
                        return minDistFromTop[n1->getId()] < minDistFromTop[n2->getId()]
                                || (minDistFromTop[n1->getId()] == minDistFromTop[n2->getId()]
                                    && maxDistFromTop[n1->getId()] <= maxDistFromTop[n2->getId()]);
                    });

                    if(useIntermediateStep){
                        std::deque<Node*> nextNextNodes;
                        while(nextNodes.size() && currentPartitionSize < maxSize){
                            Node* selectedNextNode = nextNodes.front();
                            nextNodes.pop_front();
                            selectedNextNode->setPartitionId(currentPartitionId);
                            currentPartitionSize += 1;

                            for(const auto& otherNode : selectedNextNode->getSuccessors()){
                                counterRelease[otherNode->getId()] += 1;
                                assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                                if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                                    nextNextNodes.push_back(otherNode);
                                }
                            }
                        }
                        partitionNodes.insert(nextNodes.end(), nextNextNodes.begin(), nextNextNodes.end());
                    }
                    partitionNodes.insert(partitionNodes.end(), nextNodes.begin(), nextNodes.end());
                }

                sources.insert(sources.end(), partitionNodes.begin(), partitionNodes.end());

                currentPartitionId += 1;
            }
        }
    }

    void partition(const int /*minSize*/, const int maxSize, const int /*degreeParallelism*/){
        std::vector<int> minDistFromTop(nodes.size(), -1);
        std::vector<int> maxDistFromTop(nodes.size(), -1);
        std::vector<int> maxDistFromRoot(nodes.size(), -1);

        std::vector<Node*> originalSources;
        std::vector<Node*> originalRoot;
        for(auto& node : nodes){
            node->setPartitionId(-1);
            if(node->getPredecessors().size() == 0){
                originalSources.push_back(node.get());
                minDistFromTop[node->getId()] = 0;
                maxDistFromTop[node->getId()] = 0;
            }
            if(node->getSuccessors().size() == 0){
                originalRoot.push_back(node.get());
                maxDistFromRoot[node->getId()] = 0;
            }
        }
        {
            std::vector<Node*> sources = originalSources;
            std::vector<int> counterRelease(nodes.size(), 0);
            while(sources.size()){
                Node* selectedNode = sources.back();
                sources.pop_back();

                // Add deps if released
                for(const auto& otherNode : selectedNode->getSuccessors()){
                    if(minDistFromTop[otherNode->getId()] == -1){
                        minDistFromTop[otherNode->getId()] = minDistFromTop[selectedNode->getId()] + 1;
                    }
                    else{
                        minDistFromTop[otherNode->getId()] = std::min(minDistFromTop[selectedNode->getId()] + 1, minDistFromTop[otherNode->getId()]);
                    }

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
        {
            std::vector<Node*> root = originalRoot;
            std::vector<int> counterRelease(nodes.size(), 0);
            while(root.size()){
                Node* selectedNode = root.back();
                root.pop_back();

                // Add deps if released
                for(const auto& otherNode : selectedNode->getPredecessors()){
                    if(maxDistFromRoot[otherNode->getId()] == -1){
                        maxDistFromRoot[otherNode->getId()] = maxDistFromRoot[selectedNode->getId()] + 1;
                    }
                    else{
                        maxDistFromRoot[otherNode->getId()] = std::max(maxDistFromRoot[selectedNode->getId()] + 1, maxDistFromRoot[otherNode->getId()]);
                    }

                    counterRelease[otherNode->getId()] += 1;
                    assert(counterRelease[otherNode->getId()] <= int(otherNode->getSuccessors().size()));
                    if(counterRelease[otherNode->getId()] == int(otherNode->getSuccessors().size())){
                        root.push_back(otherNode);
                    }
                }
            }
        }

        {
            std::set<Node*> sources(originalSources.begin(), originalSources.end());

            std::vector<int> counterRelease(nodes.size(), 0);

            int currentPartitionId = 0;

            while(sources.size()){
                int currentPartitionSize = 0;
                std::deque<Node*> partitionNodes;
                {
                    Node* startingNode = *sources.begin();
                    int nbReleases = 0;
                    for(Node* potentialStartingNode : sources){
                        int nbPotentialRelease = 0;
                        for(const auto& otherNode : potentialStartingNode->getSuccessors()){
                            assert(counterRelease[otherNode->getId()] < int(otherNode->getPredecessors().size()));
                            if(counterRelease[otherNode->getId()]+1 == int(otherNode->getPredecessors().size())){
                                nbPotentialRelease += 1;
                            }
                        }
                        if(nbReleases < nbPotentialRelease
                                || (nbReleases == nbPotentialRelease
                                      && maxDistFromRoot[startingNode->getId()] < maxDistFromRoot[potentialStartingNode->getId()])){
                            startingNode = potentialStartingNode;
                        }
                    }
                    sources.erase(startingNode);

                    startingNode->setPartitionId(currentPartitionId);
                    currentPartitionSize += 1;

                    partitionNodes.push_front(startingNode);
                }

                std::set<Node*> pathToPrivilegiate;
                while(partitionNodes.size()){
                    std::deque<Node*> notReadyNextNodes;
                    {
                        std::deque<Node*> nextNodes;
                        std::deque<Node*> nextNodesToPrivilegiate;
                        {
                            Node* selectedNode = partitionNodes.front();
                            partitionNodes.pop_front();

                            for(const auto& otherNode : selectedNode->getSuccessors()){
                                counterRelease[otherNode->getId()] += 1;
                                assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                                if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                                    if(pathToPrivilegiate.find(otherNode) == pathToPrivilegiate.end()){
                                        nextNodes.push_back(otherNode);
                                    }
                                    else{
                                        pathToPrivilegiate.erase(otherNode);
                                        nextNodesToPrivilegiate.push_back(otherNode);
                                    }
                                }
                                else{
                                    if(pathToPrivilegiate.find(otherNode) == pathToPrivilegiate.end()){
                                        notReadyNextNodes.push_back(otherNode);
                                    }
                                }
                            }
                        }

                        std::sort(nextNodes.begin(), nextNodes.end(), [&](const Node* n1, const Node* n2){
                            return minDistFromTop[n1->getId()] < minDistFromTop[n2->getId()]
                                    || (minDistFromTop[n1->getId()] == minDistFromTop[n2->getId()]
                                        && maxDistFromTop[n1->getId()] <= maxDistFromTop[n2->getId()]);
                        });


                        for(auto& toInsert : nextNodesToPrivilegiate){
                            if(currentPartitionSize < maxSize){
                                toInsert->setPartitionId(currentPartitionId);
                                currentPartitionSize += 1;
                                partitionNodes.push_front(toInsert);
                            }
                            else{
                                sources.insert(toInsert);
                            }
                        }
                        for(auto& toInsert : nextNodes){
                            if(currentPartitionSize < maxSize){
                                toInsert->setPartitionId(currentPartitionId);
                                currentPartitionSize += 1;
                                partitionNodes.push_back(toInsert);
                            }
                            else{
                                sources.insert(toInsert);
                            }
                        }
                    }

                    if(currentPartitionSize < maxSize && notReadyNextNodes.size() && pathToPrivilegiate.size() == 0){
                        std::set<Node*> sourcesToAdd;

                        for(auto& nodeNotReady: notReadyNextNodes){
                            pathToPrivilegiate.insert(nodeNotReady);
                        }

                        while(notReadyNextNodes.size()){
                            Node* backtrackNode = notReadyNextNodes.front();
                            notReadyNextNodes.pop_front();

                            for(const auto& backtrackPrevNode : backtrackNode->getPredecessors()){
                                if(counterRelease[backtrackPrevNode->getId()] == int(backtrackPrevNode->getPredecessors().size())){
                                    if(backtrackPrevNode->getPartitionId() == -1){
                                        assert(sources.find(backtrackPrevNode) != sources.end());
                                        sourcesToAdd.insert(backtrackPrevNode);
                                    }
                                }
                                else{
                                    assert(backtrackPrevNode->getPartitionId() == -1);
                                    pathToPrivilegiate.insert(backtrackPrevNode);
                                    notReadyNextNodes.push_back(backtrackPrevNode);
                                }
                            }
                        }

                        for(auto& nodeToAdd: sourcesToAdd){
                            assert(sources.find(nodeToAdd) != sources.end());
                            if(currentPartitionSize < maxSize){
                                nodeToAdd->setPartitionId(currentPartitionId);
                                currentPartitionSize += 1;

                                partitionNodes.push_front(nodeToAdd);
                                sources.erase(nodeToAdd);
                            }
                        }
                    }
                }

                currentPartitionId += 1;
            }
        }
    }


    std::pair<int,double> estimateDegreeOfParallelism() const{
        int maxSourcesSize = 0;
        int sumSourcesSize = 0;

        std::deque<Node*> sources;
        for(auto& node : nodes){
            if(node->getPredecessors().size() == 0){
                sources.push_back(node.get());
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

    Graph getPartitionGraph() const{
        std::vector<std::pair<int,int>> dependencyBetweenPartitions;
        std::vector<int> partitionCosts;

        for(const auto& node : nodes){
            for(const auto& otherNode : node->getSuccessors()){
                assert(node->getPartitionId() <= otherNode->getPartitionId());
                if(node->getPartitionId() != otherNode->getPartitionId()){
                    dependencyBetweenPartitions.emplace_back(std::pair<int,int>{node->getPartitionId(), otherNode->getPartitionId()});
                }
            }

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
};

#endif
