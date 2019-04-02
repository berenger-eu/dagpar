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
          const std::vector<int>& inCostPerNode) : Graph(inNodes, inDependencyList){
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

    void partitionRandom(const int maxSize){
        std::vector<Node*> sources;
        for(auto& node : nodes){
            if(node->getPredecessors().size() == 0){
                sources.push_back(node);
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

    void partitionGreedy(const int maxSize, const bool useAlternate = true, const bool useIntermediateStep = false){
        std::vector<int> minDistFromTop(nodes.size(), -1);
        std::vector<int> maxDistFromTop(nodes.size(), -1);
        std::vector<int> maxDistFromRoot(nodes.size(), -1);

        std::vector<Node*> originalSources;
        std::vector<Node*> originalRoot;
        for(auto& node : nodes){
            node->setPartitionId(-1);
            if(node->getPredecessors().size() == 0){
                originalSources.push_back(node);
                minDistFromTop[node->getId()] = 0;
                maxDistFromTop[node->getId()] = 0;
            }
            if(node->getSuccessors().size() == 0){
                originalRoot.push_back(node);
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
            std::vector<Node*> sources = originalSources;
            std::vector<int> counterRelease(nodes.size(), 0);

            int currentPartitionId = 0;

            while(sources.size()){
                int idxSelectNode = 0;
                int nbReleases = 0;
                for(int idxNode = 0 ; idxNode < int(sources.size()) ; ++idxNode){
                    if(useAlternate == false){
                        if(minDistFromTop[sources[idxNode]->getId()] < minDistFromTop[sources[idxSelectNode]->getId()]){
                            idxSelectNode = idxNode;
                        }
                    }
                    else{
                        int nbPotentialRelease = 0;
                        for(const auto& otherNode : sources[idxNode]->getSuccessors()){
                            assert(counterRelease[otherNode->getId()] < int(otherNode->getPredecessors().size()));
                            if(counterRelease[otherNode->getId()]+1 == int(otherNode->getPredecessors().size())){
                                nbPotentialRelease += 1;
                            }
                        }
                        if(nbReleases < nbPotentialRelease
                                || (nbReleases == nbPotentialRelease
                                      && maxDistFromRoot[sources[idxSelectNode]->getId()] < maxDistFromRoot[sources[idxNode]->getId()])){
                            idxSelectNode = idxNode;
                        }
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

    void partitionBacktrack(const int maxSize){
        std::vector<int> minDistFromTop(nodes.size(), -1);
        std::vector<int> maxDistFromTop(nodes.size(), -1);
        std::vector<int> maxDistFromRoot(nodes.size(), -1);

        std::vector<Node*> originalSources;
        std::vector<Node*> originalRoot;
        for(auto& node : nodes){
            node->setPartitionId(-1);
            if(node->getPredecessors().size() == 0){
                originalSources.push_back(node);
                minDistFromTop[node->getId()] = 0;
                maxDistFromTop[node->getId()] = 0;
            }
            if(node->getSuccessors().size() == 0){
                originalRoot.push_back(node);
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
                    Node* startingNode = (*sources.begin());
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
                    std::set<Node*> notReadyNextNodes;
                    {
                        std::deque<Node*> nextNodes;
                        std::deque<Node*> nextNodesToPrivilegiate;
                        {
                            Node* selectedNode = partitionNodes.front();
                            partitionNodes.pop_front();
                            assert(selectedNode->getPartitionId() == currentPartitionId);

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
                                        notReadyNextNodes.insert(otherNode);
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
                            Node* backtrackNode = (*notReadyNextNodes.begin());
                            notReadyNextNodes.erase(backtrackNode);

                            for(const auto& backtrackPrevNode : backtrackNode->getPredecessors()){
                                if(counterRelease[backtrackPrevNode->getId()] == int(backtrackPrevNode->getPredecessors().size())){
                                    if(backtrackPrevNode->getPartitionId() == -1){
                                        assert(sources.find(backtrackPrevNode) != sources.end());
                                        sourcesToAdd.insert(backtrackPrevNode);
                                    }
                                }
                                else if(notReadyNextNodes.find(backtrackPrevNode) == notReadyNextNodes.end()
                                        && pathToPrivilegiate.find(backtrackPrevNode) == pathToPrivilegiate.end()){
                                    assert(backtrackPrevNode->getPartitionId() == -1);
                                    pathToPrivilegiate.insert(backtrackPrevNode);
                                    notReadyNextNodes.insert(backtrackPrevNode);
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



    void partition(const int minSize, const int maxSize){
        assert(minSize <= maxSize);

        std::vector<Node*> originalSources;
        std::vector<Node*> originalRoot;
        for(auto& node : nodes){
            node->setPartitionId(-1);
            if(node->getPredecessors().size() == 0){
                originalSources.push_back(node);
            }
        }
        {
            std::vector<Node*> sources(originalSources);

            std::vector<int> counterRelease(nodes.size(), 0);

            int currentPartitionId = 0;

            while(sources.size()){
                int startingNodeIdx = -1;
                int nbReleases = 0;

                for(int potentialStartingNodeIdx = 0 ; potentialStartingNodeIdx < int(sources.size()) ; ++potentialStartingNodeIdx){
                    std::deque<Node*> partitionNodes;
                    {
                        Node* potentialStartingNode = sources[potentialStartingNodeIdx];
                        assert(potentialStartingNode->getPartitionId() == -1);
                        partitionNodes.push_front(potentialStartingNode);
                    }
                    std::unordered_map<int,int> counterReleaseOffset;

                    int currentPartitionSize = 0;
                    while(partitionNodes.size() && currentPartitionSize < maxSize){
                        Node* selectedNode = partitionNodes.front();
                        partitionNodes.pop_front();

                        currentPartitionSize += 1;

                        for(const auto& otherNode : selectedNode->getSuccessors()){
                            counterReleaseOffset[otherNode->getId()] += 1;
                            assert(counterRelease[otherNode->getId()] + counterReleaseOffset[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                            if(counterRelease[otherNode->getId()] + counterReleaseOffset[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                                partitionNodes.push_back(otherNode);
                            }
                        }
                    }

                    if(startingNodeIdx == -1 || nbReleases < currentPartitionSize){
                        startingNodeIdx = potentialStartingNodeIdx;
                        nbReleases = currentPartitionSize;
                        if(nbReleases == maxSize){
                            break;
                        }
                    }
                }

                if(minSize <= nbReleases || sources.size() == 1){
                    std::deque<Node*> partitionNodes;
                    {
                        Node* startingNode = sources[startingNodeIdx];
                        assert(startingNode->getPartitionId() == -1);
                        startingNode->setPartitionId(currentPartitionId);
                        partitionNodes.push_front(startingNode);
                        std::swap(sources[startingNodeIdx], sources[sources.size()-1]);
                        sources.pop_back();
                    }

                    int currentPartitionSize = 0;
                    while(partitionNodes.size() && currentPartitionSize < maxSize){
                        Node* selectedNode = partitionNodes.front();
                        partitionNodes.pop_front();

                        selectedNode->setPartitionId(currentPartitionId);
                        currentPartitionSize += 1;

                        for(const auto& otherNode : selectedNode->getSuccessors()){
                            counterRelease[otherNode->getId()] += 1;
                            assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                            if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                                assert(otherNode->getPartitionId() == -1);
                                partitionNodes.push_back(otherNode);
                            }
                        }
                    }

                    sources.insert(sources.end(), partitionNodes.begin(), partitionNodes.end());

                    currentPartitionId += 1;
                }
                else{
                    std::set<int> currentPartitionLocked;
                    std::deque<Node*> partitionNodes;
                    int currentPartitionSize = 0;
                    while(currentPartitionSize < minSize && int(sources.size())){
                        assert(partitionNodes.size() == 0);

                        std::vector<std::set<int>> releasableNodes(sources.size());
                        std::vector<std::set<int>> lockedNodes(sources.size());

                        for(int potentialStartingNodeIdx = 0 ; potentialStartingNodeIdx < int(sources.size()) ; ++potentialStartingNodeIdx){
                            std::deque<Node*> partitionNodesTest;
                            {
                                Node* potentialStartingNode = sources[potentialStartingNodeIdx];
                                partitionNodesTest.push_front(potentialStartingNode);
                            }
                            std::unordered_map<int,int> counterReleaseOffset;

                            int currentPartitionSizeTest = 0;
                            while(partitionNodesTest.size() && currentPartitionSizeTest < maxSize){
                                Node* selectedNode = partitionNodesTest.front();
                                partitionNodesTest.pop_front();

                                currentPartitionSizeTest += 1;

                                for(const auto& otherNode : selectedNode->getSuccessors()){
                                    counterReleaseOffset[otherNode->getId()] += 1;
                                    assert(counterRelease[otherNode->getId()] + counterReleaseOffset[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                                    if(counterRelease[otherNode->getId()] + counterReleaseOffset[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                                        partitionNodesTest.push_back(otherNode);
                                        releasableNodes[potentialStartingNodeIdx].insert(otherNode->getId());
                                    }
                                    else{
                                        lockedNodes[potentialStartingNodeIdx].insert(otherNode->getId());
                                    }
                                }
                            }
                        }

                        if(currentPartitionLocked.empty() && sources.size() >= 2){
                            int idxBest1 = 0;
                            int idxBest2 = 1;
                            double bestScore = 0;

                            for(int idxNode1 = 0 ; idxNode1 < int(sources.size())-1 ; ++idxNode1){
                                for(int idxNode2 = idxNode1+1 ; idxNode2 < int(sources.size()) ; ++idxNode2){
                                    std::set<int> locksIntersection;
                                    std::set_intersection(lockedNodes[idxNode1].begin(),lockedNodes[idxNode1].end(),
                                                          lockedNodes[idxNode2].begin(),lockedNodes[idxNode2].end(),
                                                      std::inserter(locksIntersection,locksIntersection.begin()));

                                    const int nbSameLocks = int(locksIntersection.size());
                                    const int score = nbSameLocks;// + int(releasableNodes[idxNode1].size()) + int(releasableNodes[idxNode2].size());

                                    if(bestScore < score){
                                        idxBest1 = idxNode1;
                                        idxBest2 = idxNode2;
                                    }
                                }
                            }

                            currentPartitionLocked = lockedNodes[idxBest1];
                            currentPartitionLocked.insert(lockedNodes[idxBest2].begin(), lockedNodes[idxBest2].end());

                            {
                                Node* startingNode1 = sources[idxBest1];
                                assert(startingNode1->getPartitionId() == -1);
                                startingNode1->setPartitionId(currentPartitionId);
                                partitionNodes.push_front(startingNode1);
                                std::swap(sources[idxBest1], sources[sources.size()-1]);
                                if(idxBest2 == int(sources.size()-1)){
                                    idxBest2 = idxBest1;
                                }
                                sources.pop_back();

                                Node* startingNode2 = sources[idxBest2];
                                assert(startingNode2->getPartitionId() == -1);
                                startingNode2->setPartitionId(currentPartitionId);
                                partitionNodes.push_front(startingNode2);
                                std::swap(sources[idxBest2], sources[sources.size()-1]);
                                sources.pop_back();
                            }
                        }
                        else{
                            int idxBest1 = 0;
                            double bestScore = 0;

                            for(int idxNode1 = 0 ; idxNode1 < int(sources.size()) ; ++idxNode1){
                                std::set<int> locksIntersection;
                                std::set_intersection(lockedNodes[idxNode1].begin(),lockedNodes[idxNode1].end(),
                                                      currentPartitionLocked.begin(),currentPartitionLocked.end(),
                                                  std::inserter(locksIntersection,locksIntersection.begin()));

                                const int nbSameLocks = int(locksIntersection.size());
                                const int score = nbSameLocks;// + int(releasableNodes[idxNode1].size()) + int(releasableNodes[idxNode2].size());

                                if(bestScore < score){
                                    idxBest1 = idxNode1;
                                }
                            }

                            currentPartitionLocked.insert(lockedNodes[idxBest1].begin(), lockedNodes[idxBest1].end());

                            Node* startingNode1 = sources[idxBest1];
                            assert(startingNode1->getPartitionId() == -1);
                            startingNode1->setPartitionId(currentPartitionId);
                            partitionNodes.push_front(startingNode1);
                            std::swap(sources[idxBest1], sources[sources.size()-1]);
                            sources.pop_back();
                        }

                        while(partitionNodes.size() && currentPartitionSize < maxSize){
                            Node* selectedNode = partitionNodes.front();
                            partitionNodes.pop_front();

                            selectedNode->setPartitionId(currentPartitionId);
                            currentPartitionSize += 1;

                            for(const auto& otherNode : selectedNode->getSuccessors()){
                                counterRelease[otherNode->getId()] += 1;
                                assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                                if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                                    assert(otherNode->getPartitionId() == -1);
                                    partitionNodes.push_back(otherNode);
                                }
                            }
                        }

                    }
                    assert(partitionNodes.empty() || currentPartitionSize == maxSize);
                    sources.insert(sources.end(), partitionNodes.begin(), partitionNodes.end());

                    currentPartitionId += 1;
                }
            }
        }
    }

    void partitionHorizontal(const int maxSize){
        std::vector<int> minDistFromTop(nodes.size(), -1);
        std::vector<int> maxDistFromTop(nodes.size(), -1);
        std::vector<int> minDistFromRoot(nodes.size(), -1);
        std::vector<int> maxDistFromRoot(nodes.size(), -1);

        std::vector<Node*> originalSources;
        std::vector<Node*> originalRoot;
        for(auto& node : nodes){
            node->setPartitionId(-1);
            if(node->getPredecessors().size() == 0){
                originalSources.push_back(node);
                minDistFromTop[node->getId()] = 0;
                maxDistFromTop[node->getId()] = 0;
            }
            if(node->getSuccessors().size() == 0){
                originalRoot.push_back(node);
                minDistFromRoot[node->getId()] = 0;
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
                    if(minDistFromRoot[otherNode->getId()] == -1){
                        minDistFromRoot[otherNode->getId()] = minDistFromRoot[selectedNode->getId()] + 1;
                    }
                    else{
                        minDistFromRoot[otherNode->getId()] = std::min(minDistFromRoot[selectedNode->getId()] + 1, minDistFromRoot[otherNode->getId()]);
                    }

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
            std::vector<Node*> sources(originalSources);

            std::vector<int> counterRelease(nodes.size(), 0);

            int currentPartitionId = 0;

            while(sources.size()){
                std::deque<Node*> releasedNodes;
                std::vector<std::set<int>> depsNodes(sources.size());
                std::set<int> partitionNodes;
                std::deque<Node*> sourcesSameWidth;
                {
                    int startingNodeIdx = -1;
                    for(int potentialStartingNodeIdx = 0 ; potentialStartingNodeIdx < int(sources.size()) ; ++potentialStartingNodeIdx){
                        Node* potentialStartingNode = sources[potentialStartingNodeIdx];
                        assert(potentialStartingNode->getPartitionId() == -1);

                        for(const auto& otherNode : potentialStartingNode->getSuccessors()){
                            depsNodes[potentialStartingNodeIdx].insert(otherNode->getId());
                        }

                        if(startingNodeIdx == -1
                                ||  minDistFromTop[potentialStartingNode->getId()] < minDistFromTop[sources[startingNodeIdx]->getId()]){
                            startingNodeIdx = potentialStartingNodeIdx;
                        }
                    }

                    for(const auto& otherNode : sources[startingNodeIdx]->getSuccessors()){
                        counterRelease[otherNode->getId()] += 1;
                        assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                        if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                            releasedNodes.push_back(otherNode);
                        }
                    }

                    const int idxSelectedDeep = minDistFromTop[sources[startingNodeIdx]->getId()];
                    partitionNodes = depsNodes[startingNodeIdx];

                    sources[startingNodeIdx]->setPartitionId(currentPartitionId);

                    depsNodes[startingNodeIdx] = std::move(depsNodes[depsNodes.size()-1]);
                    depsNodes.pop_back();
                    sources[startingNodeIdx] = std::move(sources[sources.size()-1]);
                    sources.pop_back();

                    int idxSources = 0;
                    while(idxSources < int(sources.size())){
                        if(minDistFromTop[sources[idxSources]->getId()] == idxSelectedDeep){
                            sourcesSameWidth.push_back(sources[idxSources]);
                            std::swap(sources[idxSources], sources[sources.size()-1]);
                            sources.pop_back();
                        }
                        else{
                            idxSources += 1;
                        }
                    }
                }

                int sizeCurrentPartition = 1;

                while(sizeCurrentPartition < maxSize && sourcesSameWidth.size()){
                    int bestIdx = -1;
                    int bestIntersectionSize = 0;

                    for(int potentialNeighborIdx = 0 ; potentialNeighborIdx < int(sourcesSameWidth.size()) ; ++potentialNeighborIdx){
                        std::set<int> locksIntersection;
                        std::set_intersection(partitionNodes.begin(),partitionNodes.end(),
                                              depsNodes[potentialNeighborIdx].begin(),depsNodes[potentialNeighborIdx].end(),
                                          std::inserter(locksIntersection,locksIntersection.begin()));

                        if(bestIdx == -1 ||  bestIntersectionSize < int(locksIntersection.size())){
                            bestIdx = potentialNeighborIdx;
                            bestIntersectionSize = int(locksIntersection.size());
                        }
                    }

                    for(const auto& otherNode : sourcesSameWidth[bestIdx]->getSuccessors()){
                        counterRelease[otherNode->getId()] += 1;
                        assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                        if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                            releasedNodes.push_back(otherNode);
                        }
                    }

                    partitionNodes.insert(depsNodes[bestIdx].begin(), depsNodes[bestIdx].end());

                    sourcesSameWidth[bestIdx]->setPartitionId(currentPartitionId);

                    depsNodes[bestIdx] = std::move(depsNodes[depsNodes.size()-1]);
                    depsNodes.pop_back();
                    sourcesSameWidth[bestIdx] = std::move(sourcesSameWidth[sourcesSameWidth.size()-1]);
                    sourcesSameWidth.pop_back();

                    sizeCurrentPartition += 1;
                }

                sources.insert(sources.end(), sourcesSameWidth.begin(), sourcesSameWidth.end());
                sources.insert(sources.end(), releasedNodes.begin(), releasedNodes.end());

                currentPartitionId += 1;
            }
        }
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
                    for(const auto& idxChild : proceedPartitionsInfo[testPartitionId].idsParentPartitionsSuccessors){
                        if(selectedNodeParentPartitionIds.find(idxChild) != selectedNodeParentPartitionIds.end()){
                            testPartitionIdIsLinkedToAParent = true;
                            break;
                        }
                        children.push_back(idxChild);
                    }

                    while(testPartitionIdIsLinkedToAParent == false && children.size()){
                        const int testPartitionIter = children.front();
                        children.pop_front();

                        for(const auto& idxChild : proceedPartitionsInfo[testPartitionIter].idsParentPartitionsSuccessors){
                            if(selectedNodeParentPartitionIds.find(idxChild) != selectedNodeParentPartitionIds.end()){
                                testPartitionIdIsLinkedToAParent = true;
                                break;
                            }
                            children.push_back(idxChild);
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
        std::vector<int> partitionCosts;

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

#ifdef USE_METIS
    void partitionMetis(const int partSize, const bool useKway = true){
        idx_t nvtxs = getNbNodes();
        idx_t ncon = 1;
        idx_t* vwgt = nullptr; // weights of vertices
        idx_t* vsize = nullptr; // size of vertices for computing total communication volum
        idx_t* adjwgt = nullptr; // weights of edges
        idx_t nparts = (getNbNodes() + partSize - 1)/partSize;
        real_t* tpwgts = nullptr; // desired weight
        real_t* ubvec = nullptr; // allowed load imbalance
        idx_t objval = 0; // result status
        std::unique_ptr<idx_t[]> part(new idx_t[nvtxs]());

        std::unique_ptr<idx_t[]> xadj(new idx_t[nvtxs+1]()); // Start indices in edgetab
        xadj[0] = 0;

        std::vector<idx_t> adjncy; // adjacency array for every vertex

        for(int idxNode = 0 ; idxNode < getNbNodes() ; ++idxNode){
            const auto& node = nodesNotTopological[idxNode];
            assert(node->getId() == idxNode);
            const auto& successors = node->getSuccessors();
            const int nbDeps = int(successors.size());
            xadj[idxNode+1] = xadj[idxNode] + nbDeps;

            assert(int(adjncy.size()) == xadj[idxNode]);
            int idxDep = 0;
            for(const auto& dep : successors){
                adjncy.push_back(dep->getId());
                idxDep += 1;
            };
            assert(int(adjncy.size()) == xadj[idxNode+1]);
        }

        idx_t options[METIS_NOPTIONS] = {0};
        METIS_SetDefaultOptions(options);
        //options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // METIS PTYPE RB
        //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // METIS_OBJTYPE_VOL
        options[METIS_OPTION_NUMBERING] = 0;
        options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;

        if(useKway){
            int result;
            if((result = METIS_PartGraphKway(&nvtxs, &ncon, xadj.get(), adjncy.data(),
                                    vwgt, vsize, adjwgt, &nparts, tpwgts,
                                    ubvec, options, &objval, part.get())) != METIS_OK){
                std::cout <<"[SCOTCH] Error " << result << std::endl;
                return;
            }
        }
        else{
            int result;
            if((result = METIS_PartGraphRecursive(&nvtxs, &ncon, xadj.get(), adjncy.data(),
                                                vwgt, vsize, adjwgt, &nparts, tpwgts,
                                                ubvec, options, &objval, part.get())) != METIS_OK){
                std::cout <<"[SCOTCH] Error " << result << std::endl;
                return;
            }
        }

        for(int idxNode = 0 ; idxNode < getNbNodes() ; ++idxNode){
            auto& node = nodesNotTopological[idxNode];
            node->setPartitionId(part[idxNode]);
        }

        // Order partitions
        std::deque<Node*> sources;
        for(auto& node : nodes){
            node->setPartitionId(node->getPartitionId() + getNbNodes());
            if(node->getPredecessors().size() == 0){
                sources.push_back(node);
            }
        }

        std::unordered_map<int,int> partitionsMapping;

        std::vector<int> counterRelease(nodes.size(), 0);

        while(sources.size()){
            Node* selectedNode = sources.front();
            sources.pop_front();

            if(partitionsMapping.find(selectedNode->getPartitionId()) == partitionsMapping.end()){
                partitionsMapping[selectedNode->getPartitionId()] = int(partitionsMapping.size());
            }
            selectedNode->setPartitionId(partitionsMapping[selectedNode->getPartitionId()]);

            for(const auto& otherNode : selectedNode->getSuccessors()){
                counterRelease[otherNode->getId()] += 1;
                assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                    sources.push_back(otherNode);
                }
            }
        }
    }
#endif

    // Derived from TEMPORAL  PARTITIONING  AND  SCHEDULING  DATA  FLOW  GRAPHS  FOR  RECONFIGURABLE
    // The fig 9 is not true and this is easy to see because the ready list
    // is a lifo
    void partitionTemporalPart(const int maxSize) const{
        std::vector<int> counterRelease(nodes.size(), 0);

        int current_size = 0;
        int i = 0;

        std::deque<Node*> sources;
        for(auto& node : nodes){
            if(node->getPredecessors().size() == 0){

                node->setPartitionId(i);
                current_size += 1;

                if(current_size == maxSize){
                    i += 1;
                    current_size = 0;
                }

                for(const auto& otherNode : node->getSuccessors()){
                    counterRelease[otherNode->getId()] += 1;
                    assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                    if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                        sources.push_back(otherNode);
                    }
                }
            }
        }

        while(sources.size()){
            Node* selectedNode = sources.back();
            sources.pop_back();

            selectedNode->setPartitionId(i);
            current_size += 1;
            if(current_size == maxSize){
                i += 1;
                current_size = 0;
            }

            for(const auto& otherNode : selectedNode->getSuccessors()){
                counterRelease[otherNode->getId()] += 1;
                assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                    sources.push_back(otherNode);
                }
            }
        }
    }

#ifdef USE_ACYCLIC
    // Derived from Acyclic Partitioning of Large Directed Acyclic Graphs
    void partitionAcyclic(const int maxSize){
        {// CompMatching
            std::vector<int> top(getNbNodes(), std::numeric_limits<int>::max());
            {
                std::deque<Node*> sources;
                for(auto& node : nodes){
                    if(node->getPredecessors().size() == 0){
                        top[node->getId()] = 0;
                        sources.push_back(node);
                    }
                }

                std::vector<int> counterRelease(nodes.size(), 0);

                while(sources.size()){
                    Node* selectedNode = sources.back();
                    sources.pop_back();

                    for(const auto& otherNode : selectedNode->getSuccessors()){
                        counterRelease[otherNode->getId()] += 1;
                        assert(counterRelease[otherNode->getId()] <= int(otherNode->getPredecessors().size()));
                        if(counterRelease[otherNode->getId()] == int(otherNode->getPredecessors().size())){
                            top[otherNode->getId()] = std::min(top[otherNode->getId()], top[selectedNode->getId()] + 1);
                            sources.push_back(otherNode);
                        }
                    }
                }
            }

            std::set<std::pair<int,int>> M;

            std::vector<bool> mark(getNbNodes(), false);

            for(int idxNode = 0 ; idxNode < getNbNodes() ; ++idxNode){
                Node* u = nodes[idxNode];

                if(mark[u->getId()]) continue;

                for(const auto& v : u->getPredecessors()){
                    if(mark[v->getId()]) continue;

                    if(top[u->getId()] != top[v->getId()]-1 && v->getPredecessors().size() != 1 && u->getSuccessors().size() != 1) continue;

                    M.insert(std::make_pair(v->getId(), u->getId()));

                    for(const auto& w : v->getSuccessors()){
                        if(top[v->getId()] == top[w->getId()]-1){
                            mark[v->getId()] = false;
                        }
                    }

                    mark[v->getId()] = true;
                    mark[u->getId()] = true;
                }

                for(const auto& v : u->getSuccessors()){
                    if(mark[v->getId()]) continue;

                    if(top[u->getId()] != top[v->getId()]-1 && v->getPredecessors().size() != 1 && u->getSuccessors().size() != 1) continue;

                    M.insert(std::make_pair(u->getId(), v->getId()));

                    for(const auto& w : u->getSuccessors()){
                        if(top[u->getId()] != top[w->getId()]-1){
                            mark[w->getId()] = false;
                        }
                    }

                    mark[v->getId()] = true;
                    mark[u->getId()] = true;
                }
            }
        }

        auto ComputeEdgeCut = [this](const std::vector<int>& part) -> int {
            int nbCuts = 0;
            for(auto& u : nodes){
                for(const auto& v : u->getSuccessors()){
                    if(part[v->getId()] != part[u->getId()]){
                        nbCuts += 1;
                    }
                }
            }
            return nbCuts;
        };

        auto CompGain = [this](const Node* u, const std::vector<int>& part, const int d) -> int {
            int gain = 0;

            for(const auto& v : u->getSuccessors()){
                const auto cuv = 1;
                const auto cvu = 1;
                if(part[v->getId()] == part[u->getId()]){
                    gain -= cuv;
                }
                else if(part[v->getId()] == d){
                    gain += cvu;
                }
            }

            for(const auto& v : u->getPredecessors()){
                const auto cuv = 1;
                const auto cvu = 1;
                if(part[v->getId()] == part[u->getId()]){
                    gain -= cuv;
                }
                else if(part[v->getId()] == d){
                    gain += cvu;
                }
            }


            return gain;
        };

        std::vector<int> part(getNbNodes(), -1);
        {// Greedy algorithm
            int k = (getNbNodes()+maxSize-1)/maxSize;
            const double lb = 1.1 * (getNbNodes() / k);
            std::vector<bool> free(getNbNodes(), true);

            for(int i = 0 ; i < k ; ++i){
                std::set<Node*> set;

                for(auto& u : nodes){
                    if(free[u->getId()]){
                        assert(part[u->getId()] == -1);
                        if(u->getPredecessors().size() == 0){
                            set.insert(u);
                        }
                        else{
                            bool ready = true;
                            for(auto& v : u->getPredecessors()){
                                if(free[v->getId()] == true){
                                    ready = false;
                                    break;
                                }
                            }
                            if(ready){
                                set.insert(u);
                                break;
                            }
                        }
                    }
                }

                auto cmp = [](const auto& n1, const auto& n2){
                    return n1.first > n2.first;
                };

                std::priority_queue<std::pair<int,Node*>,  std::vector<std::pair<int,Node*>>, decltype(cmp)> heap(cmp);

                for(auto& u : set){
                    const int gain = CompGain(u, part, i);
                    heap.push(std::make_pair(gain,u));
                }

                int size_Vi = 0;
                while(size_Vi < lb && heap.size()){
                    Node* u = heap.top().second;
                    heap.pop();

                    part[u->getId()] = i;
                    free[u->getId()] = false;

                    for(auto& v : u->getSuccessors()){
                        bool ready = true;

                        for(auto& w : v->getPredecessors()){
                            if(free[w->getId()] == true){
                                ready = false;
                                break;
                            }
                        }

                        if( ready ){
                            assert(free[v->getId()]);
                            const int gain = CompGain(v, part, i);
                            heap.push(std::make_pair(gain,v));
                        }
                    }

                    size_Vi += 1;
                }
            }
            for(auto& u : nodes){
                if(part[u->getId()] == -1){
                    part[u->getId()] = k;
                }
            }
        }

        auto UpdateHeap = [](auto& heap, auto&& newItem){
            std::vector<typename std::decay<decltype (heap.top())>::type> values;
            while(heap.size()){
                auto item = heap.top();
                heap.pop();
                if(item.first != newItem.first){
                    values.push_back(item);
                }
            }
            heap.push(newItem);
            for(auto& item : values){
                heap.push(item);
            }
        };

        auto Update = [&CompGain, &UpdateHeap](auto& _heap, std::vector<bool>& moveto, const std::vector<int>& _part,
                                  auto& _gain, auto& _u){
            const int _k = int(_part.size());
            int max = std::numeric_limits<int>::min();
            int min = std::numeric_limits<int>::max();

            if(_u->getPredecessors().size() == 0){
                max = std::max(_part[_u->getId()]-1, 1);
            }
            else{
                for(auto& _v : _u->getPredecessors()){
                    max = std::max(_part[_v->getId()], max);
                }
            }

            if(_u->getSuccessors().size() == 0){
                min = std::min(_part[_u->getId()]+1, _k);
            }
            else{
                for(auto& _v : _u->getPredecessors()){
                    min = std::min(_part[_v->getId()], min);
                }
            }

            if( max != _part[_u->getId()] && min == _part[_u->getId()]){
                moveto[_u->getId()] = max;
                _gain[_u->getId()][max] = CompGain(_u, _part, max);
                //_heap.push(std::make_pair(_u, _gain[_u->getId()][max]));
                UpdateHeap(_heap, std::make_pair(_u, _gain[_u->getId()][max]));
            }
            else if( min != _part[_u->getId()] && max == _part[_u->getId()]){
                moveto[_u->getId()] = min;
                _gain[_u->getId()][min] = CompGain(_u, _part, min);
                //_heap.push(std::make_pair(_u, _gain[_u->getId()][min]));
                UpdateHeap(_heap, std::make_pair(_u, _gain[_u->getId()][min]));
            }
            else if( min != _part[_u->getId()] && max != _part[_u->getId()]){
                const int gain1 = CompGain(_u, _part, min);
                const int gain2 = CompGain(_u, _part, max);
                if(gain1 > gain2){
                    moveto[_u->getId()] = min;
                    _gain[_u->getId()][min] = gain1;
                    //_heap.push(std::make_pair(_u, _gain[_u->getId()][min]));
                    UpdateHeap(_heap, std::make_pair(_u, _gain[_u->getId()][min]));
                }
                else{
                    moveto[_u->getId()] = max;
                    _gain[_u->getId()][max] = gain2;
                    //_heap.push(std::make_pair(_u, _gain[_u->getId()][max]));
                    UpdateHeap(_heap, std::make_pair(_u, _gain[_u->getId()][max]));
                }
            }
        };

        auto BuildQuotientGraph = [this](const std::vector<int>& _part) -> Graph {
            for(auto& _u : nodes){
                _u->setPartitionId(_part[_u->getId()]);
            }
            return getPartitionGraph();
        };

        auto CreateCycle = [this](auto& /*QG*/, auto& _u, int _k, const std::vector<int>& _part) -> bool {
            for(auto& _v : nodes){
                _v->setPartitionId(_part[_v->getId()]);
            }
            _u->setPartitionId(_k);

            return isPartitionGraphDag() == false;
        };

        auto UpdateQuotientGraph = [this](auto& _QG, auto& _u, int _k, const std::vector<int>& _part){
            for(auto& _v : nodes){
                _v->setPartitionId(_part[_v->getId()]);
            }
            _u->setPartitionId(_k);

            _QG = getPartitionGraph();
        };

        {// Acyclic kway refinement
            std::vector<std::vector<int>> gain(getNbNodes());

            std::vector<int> copy = part;
            std::vector<bool> moved(getNbNodes(), false);
            std::vector<int> moves;

            int ec = ComputeEdgeCut(part);
            int ecmin = ec;

            Graph QG = BuildQuotientGraph(part);

            std::vector<int> maxgain(getNbNodes(), 0);
            for(auto& u : nodes){
                gain[u->getId()].resize(part.size());
                for(int k = 0 ; k < int(part.size()); ++k){
                    gain[u->getId()][k] = CompGain(u, part, k);
                    maxgain[u->getId()] = std::max(maxgain[u->getId()], gain[u->getId()][k]);
                }
            }

            auto cmp = [](const auto& n1, const auto& n2){
                return n1.second > n2.second;
            };

            std::priority_queue<std::pair<Node*,int>,  std::vector<std::pair<Node*, int>>, decltype(cmp)> heap(cmp);
            for(auto& u : nodes){
                heap.push(std::make_pair(u, maxgain[u->getId()]));
            }

            int idx = 0;
            int ecidx = 0;

            while(!heap.empty()){
                Node* u = heap.top().first;
                heap.pop();

                std::vector<std::pair<int,int>> parts(part.size());
                for(int k = 0 ; k < int(part.size()); ++k){
                    parts[k].first = k;
                    parts[k].second = gain[u->getId()][k];
                }

                std::sort(parts.begin(), parts.end(), [](const auto& p1, const auto& p2){
                    return p1.second < p2.second;
                });

                int i = 0;
                int k = parts[0].first;

                while(i < int(part.size())  && CreateCycle(QG, u, k, part) && gain[u->getId()][k] >= gain[u->getId()][parts[2].first]){
                    i += 1;
                    k = parts[i].first;
                }

                moves.push_back(u->getId());
                idx += 1;
                assert(int(moves.size()) == idx);

                part[u->getId()] = k;
                UpdateQuotientGraph(QG, u, k, part);

                ec -= gain[u->getId()][k];
                if( ec < ecmin ){
                    ecmin = ec;
                    ecidx = idx;
                }

                for(const auto& v : u->getSuccessors()){
                    if(!moved[v->getId()]){
                        Update(heap, moved, part, gain, v);
                    }
                }

                for(const auto& v : u->getPredecessors()){
                    if(!moved[v->getId()]){
                        Update(heap, moved, part, gain, v);
                    }
                }
            }

            for(int i = idx- 1 ; i >= ecidx ; --i){
                part[moves[i]] = copy[moves[i]];
            }
        }

        {
            for(auto& u : nodes){
                u->setPartitionId(part[u->getId()]);
            }
            return;
        }
    }
#endif
};

#endif
