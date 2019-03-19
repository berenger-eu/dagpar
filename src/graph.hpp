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


    void partitionDiamond(const int maxSize){
        // This algorithm will work only if:
        // - there is one root
        // - each node has at most 3 pred/succ dependencies
        std::vector<Node*> originalSources;
        for(auto& node : nodes){
            node->setPartitionId(-1);
            if(node->getPredecessors().size() == 0){
                originalSources.push_back(node);
            }
            const bool depsAreValid = (int(node->getPredecessors().size() + node->getSuccessors().size()) <= 4)
                    && int(node->getPredecessors().size()) >= 0 && int(node->getSuccessors().size()) >= 0
                    && int(node->getPredecessors().size()) < 4 && int(node->getSuccessors().size()) < 4;
            if(!depsAreValid){
                std::cerr << "[GRAPH][ERROR] the number of dependencies are invlid for node " << node->getId() << "\n";
                std::cerr << "[GRAPH][ERROR] nb predecessors " << node->getPredecessors().size() << "\n";
                std::cerr << "[GRAPH][ERROR] nb successors " << node->getSuccessors().size() << "\n";
                std::cerr << "[GRAPH][ERROR] the code will continue...." << std::endl;
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

        int nbLevelsInHalfDiamond = 1;
        while(((nbLevelsInHalfDiamond+2)*(nbLevelsInHalfDiamond+1))/2 < maxSize/2){
            nbLevelsInHalfDiamond += 1;
        }

        int currentPartitionId = 0;

        for(int idxNode = 0 ; idxNode < int(maxDistFromTopWithNode.size()) ; ++idxNode){
            Node* selectedNode = maxDistFromTopWithNode[idxNode].second;

            if(selectedNode->getPartitionId() == -1){
                selectedNode->setPartitionId(currentPartitionId);

                const int selectedNodeDistanceFromTop = maxDistFromTop[selectedNode->getId()];

                std::vector<Node*> nodesAtPreviousLevel;
                nodesAtPreviousLevel.push_back(selectedNode);

                for(int idxIncreaseLevel = 1 ; idxIncreaseLevel <= nbLevelsInHalfDiamond ; ++idxIncreaseLevel){
                    std::vector<Node*> nodesAtCurrentLevel;

                    for(auto& parent : nodesAtPreviousLevel){
                        for(const auto& currentNode : parent->getSuccessors()){
                            if(currentNode->getPartitionId() == -1
                                    && maxDistFromTop[currentNode->getId()] == selectedNodeDistanceFromTop + idxIncreaseLevel){
                                currentNode->setPartitionId(currentPartitionId);
                                nodesAtCurrentLevel.push_back(currentNode);
                            }
                        }
                    }

                    nodesAtPreviousLevel = std::move(nodesAtCurrentLevel);
                }

                if(nodesAtPreviousLevel.size()){
                    std::vector<std::vector<std::set<Node*>>> allSubNodesPerLevel(nodesAtPreviousLevel.size());
                    for(int idxNodeLevel  = 0 ; idxNodeLevel < int(nodesAtPreviousLevel.size()) ; ++idxNodeLevel){
                        allSubNodesPerLevel[idxNodeLevel].resize(2 * nbLevelsInHalfDiamond + 1);

                        std::vector<Node*> subNodes;
                        subNodes.push_back(nodesAtPreviousLevel[idxNodeLevel]);

                        while(subNodes.size()){
                            Node* currentNode = subNodes.back();
                            subNodes.pop_back();

                            for(const auto& child : currentNode->getSuccessors()){
                                if(maxDistFromTop[child->getId()] < selectedNodeDistanceFromTop + 2 * nbLevelsInHalfDiamond + 1){
                                    const int levelPosition = maxDistFromTop[child->getId()] - selectedNodeDistanceFromTop - (nbLevelsInHalfDiamond + 1);
                                    assert(0 <= levelPosition && levelPosition < 2 * nbLevelsInHalfDiamond + 1);
                                    allSubNodesPerLevel[idxNodeLevel][levelPosition].insert(child);
                                    subNodes.push_back(child);
                                }
                            }
                        }
                    }

                    std::vector<Node*> completeIntersection;

                    for(int idxDecreaseLevel = 0 ; idxDecreaseLevel < nbLevelsInHalfDiamond ; ++idxDecreaseLevel){
                        std::set<Node*> currentLevelIntersection;

                        currentLevelIntersection = allSubNodesPerLevel[0][idxDecreaseLevel];

                        for(int idxNodeLevel  = 1 ; idxNodeLevel < int(nodesAtPreviousLevel.size()) ; ++idxNodeLevel){
                            std::set<Node*> intersect;
                            std::set_intersection(currentLevelIntersection.begin(),currentLevelIntersection.end(),
                                             allSubNodesPerLevel[idxNodeLevel][idxDecreaseLevel].begin(),allSubNodesPerLevel[idxNodeLevel][idxDecreaseLevel].end(),
                                              std::inserter(intersect,intersect.begin()));
                            currentLevelIntersection = std::move(intersect);
                        }

                        std::copy(currentLevelIntersection.begin(), currentLevelIntersection.end(), std::back_inserter(completeIntersection));
                        if(currentLevelIntersection.size() <= 1){
                            break;
                        }
                    }

                    for(auto& currentNode : completeIntersection){
                        if(currentNode->getPartitionId() == -1){
                            currentNode->setPartitionId(currentPartitionId);
                        }
                    }
                }

                currentPartitionId += 1;
            }
        }

        {
            std::set<std::pair<int,int>> depsBetweenPartitions;
            for(auto& node : nodes){
                for(const auto& otherNode : node->getSuccessors()){
                    if(node->getPartitionId() != otherNode->getPartitionId()){
                        depsBetweenPartitions.insert(std::pair<int,int>(node->getPartitionId(), otherNode->getPartitionId()));
                    }
                }
            }

            std::vector<std::pair<int,int>> depsBetweenPartitionsVector;
            depsBetweenPartitionsVector.reserve(depsBetweenPartitions.size());
            std::copy(depsBetweenPartitions.begin(), depsBetweenPartitions.end(), depsBetweenPartitionsVector.begin());

            Graph partionGraph(currentPartitionId, depsBetweenPartitionsVector);

            std::unordered_map<int,int> partitionTolopologicalMapping;

            for (int idxPartion = 0 ; idxPartion < partionGraph.getNbNodes() ; ++idxPartion) {
                partitionTolopologicalMapping[partionGraph.getNode(idxPartion)->getId()] = idxPartion;
            }

            for(auto& node : nodes){
                node->setPartitionId(partitionTolopologicalMapping[node->getPartitionId()]);
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
};

#endif
