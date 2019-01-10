#include <list>
#include <cassert>
#include <algorithm>

class Node{
    std::list<Node*> successors;
    std::list<Node*> predecessors;

    int id;
    int partitionId;

public:
    explicit Node(const int inId, const int inPartitionId = 0)
        : id(inId), partitionId(inPartitionId){
    }

    void setPartitionId(const int inPartitionId){
        partitionId = inPartitionId;
    }

    int getId() const{
        return  id;
    }

    int getPartitionId() const{
        return partitionId;
    }

    void addSuccessor(Node* next){
        assert(id < next->id);
        assert(std::find(successors.begin(), successors.end(), next) == successors.end());
        successors.push_back(next);
    }

    void addPredecessor(Node* prev){
        assert(prev->id < id);
        assert(std::find(predecessors.begin(), predecessors.end(), prev) == predecessors.end());
        predecessors.push_back(prev);
    }

    const std::list<Node*>& getSuccessors() const{
        return successors;
    }

    const std::list<Node*>& getPredecessors() const{
        return predecessors;
    }
};


#include <utility>
#include <vector>
#include <memory>

// For the export to dot file
#include <string>
#include <fstream>
#include <unordered_map>
#include <array>
// For colors
#include <random>

class Graph{
    std::vector<std::unique_ptr<Node>> nodes;
public:
    Graph(const std::vector<std::pair<int,int>>& inDependencyList){
        for(const auto& dep : inDependencyList){
            assert(dep.first < dep.second);
            if(nodes.size() <= dep.first){
                nodes.resize(dep.first+1);
            }
            if(!nodes[dep.first]){
                nodes[dep.first].reset(new Node(dep.first));
            }
            if(nodes.size() <= dep.second){
                nodes.resize(dep.second+1);
            }
            if(!nodes[dep.second]){
                nodes[dep.second].reset(new Node(dep.second));
            }

            nodes[dep.first]->addSuccessor(nodes[dep.second].get());
            nodes[dep.second]->addPredecessor(nodes[dep.first].get());
        }
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
            dotFile << node->getId() << " [style=filled,color=\"" << color[0] << " " << color[1] << " " << color[2] << "\"]\n";
        }


        dotFile << "}\n";
        dotFile.close();
    }

    void partition(const int minSize, const int maxSize, const int degreeParallelism){
        int currentPartitionId = 0;
        int currentPartitionSize = 0;
        for(auto& node : nodes){
            node->setPartitionId(currentPartitionId);
            currentPartitionSize += 1;
            if(currentPartitionSize == maxSize){
                currentPartitionId += 1;
                currentPartitionSize = 0;
            }
        }
    }


    Graph getPartitionGraph() const{
        std::vector<std::pair<int,int>> dependencyBetweenPartitions;

        for(const auto& node : nodes){
            for(const auto& otherNode : node->getSuccessors()){
                dependencyBetweenPartitions.emplace_back(std::pair<int,int>{node->getPartitionId(), otherNode->getPartitionId()});
            }
        }

        std::sort(dependencyBetweenPartitions.begin(), dependencyBetweenPartitions.end(),
                  [](const std::pair<int,int>& dep1, const std::pair<int,int>& dep2){
            return dep1.first < dep2.first || (dep1.first == dep2.first && dep1.second < dep2.second);
        });

        auto last = std::unique(dependencyBetweenPartitions.begin(), dependencyBetweenPartitions.end());
        dependencyBetweenPartitions.erase(last, dependencyBetweenPartitions.end());

        return Graph(dependencyBetweenPartitions);
    }
};



#include <iostream>

int main(){
    std::vector<std::pair<int,int>> someDeps{{0,1}, {0,2}, {1,2}};
    Graph aGraph(someDeps);
    aGraph.partition(1,1,3);
    aGraph.saveToDot("/tmp/agraph.dot");

    Graph depGraph = aGraph.getPartitionGraph();
    aGraph.saveToDot("/tmp/depgraph.dot");
    return 0;
}
