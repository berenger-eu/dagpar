#include <list>
#include <cassert>
#include <algorithm>

class Node{
    std::list<Node*> successors;
    std::list<Node*> predecessors;

    int id;
    int partitionId;
    int cost;

public:
    explicit Node(const int inId, const int inPartitionId = 0)
        : id(inId), partitionId(inPartitionId), cost(1){
    }

    void setPartitionId(const int inPartitionId){
        partitionId = inPartitionId;
    }

    void setCost(const int inCost){
        cost = inCost;
    }

    int getId() const{
        return  id;
    }

    int getPartitionId() const{
        return partitionId;
    }

    int getCost() const{
        return cost;
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


#include <cassert>
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
    explicit Graph(const std::vector<std::pair<int,int>>& inDependencyList){
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


    Graph(const std::vector<std::pair<int,int>>& inDependencyList,
          const std::vector<int>& inCostPerNode) : Graph(inDependencyList){
        assert(nodes.size() == inCostPerNode.size());

        for(int idxNode = 0 ; idxNode < nodes.size() ; ++idxNode){
            nodes[idxNode]->setCost(inCostPerNode[idxNode]);
        }
    }

    int getNbNodes() const{
        return nodes.size();
    }

    const Node* getNode(const int inIdxNode) const{
        assert(inIdxNode < nodes.size());
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
            dotFile << node->getId() << " [label=\"" << node->getCost() << "\", style=filled,color=\"" << color[0] << " " << color[1] << " " << color[2] << "\"]\n";
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
        std::vector<int> partitionCosts;

        for(const auto& node : nodes){
            for(const auto& otherNode : node->getSuccessors()){
                dependencyBetweenPartitions.emplace_back(std::pair<int,int>{node->getPartitionId(), otherNode->getPartitionId()});
            }

            if(partitionCosts.size() <= node->getPartitionId()){
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

        return Graph(dependencyBetweenPartitions, partitionCosts);
    }
};


#include <tuple>
#include <iostream>
#include <queue>

class Executor{
public:
    class Event{
        int workerId;
        int nodeId;
        int cost;
        int startingPoint;
    public:
        Event(const int inWorkerId, const int inNodeId, const int inCost, const int inStartingPoint)
            : workerId(inWorkerId), nodeId(inNodeId), cost(inCost), startingPoint(inStartingPoint){}
    };

    static std::tuple<int,std::vector<Event>> Execute(const Graph& inGraph, const int inNbWorkers){
        if(inGraph.getNbNodes() == 0){
            return std::make_tuple(0,std::vector<Event>());
        }

        struct Worker{
            int scheduledAvailable;
            int currentTaskId;
            int currentWorkerId;
            bool operator<(const Worker& other) const{
                return scheduledAvailable > other.scheduledAvailable;
            }
        };

        std::vector<int> idleWorkerCount;

        idleWorkerCount.resize(inNbWorkers);
        std::iota(idleWorkerCount.begin(), idleWorkerCount.end(), 0);

        std::vector<int> readyTasks;

        for(int idxNode = 0 ; idxNode < inGraph.getNbNodes() ; ++idxNode){
            const auto& node = inGraph.getNode(idxNode);
            if(node->getPredecessors().size() == 0){
                readyTasks.push_back(node->getId());
            }
        }

        std::cout << "[EXECUTION] Ready tasks at starting time => " << readyTasks.size() << std::endl;
        assert(readyTasks.size());

        std::vector<int> countPredecessorsOver(inGraph.getNbNodes(), 0);
        std::priority_queue<Worker> workers;
        int nbComputedTask = 0;
        int currentTime = 0;
        std::vector<Event> events;
        events.reserve(inGraph.getNbNodes());

        while(readyTasks.size() && idleWorkerCount.size()){
            const int readyTaskId = readyTasks.back();
            readyTasks.pop_back();

            const int workerId = idleWorkerCount.back();
            idleWorkerCount.pop_back();

            Worker wk{currentTime + inGraph.getNode(readyTaskId)->getCost(),
                     readyTaskId,
                     workerId};
            workers.push(wk);

            events.push_back(Event(workerId, readyTaskId, inGraph.getNode(readyTaskId)->getCost(), currentTime));

            nbComputedTask += 1;
        }

        assert(workers.size() != 0);

        while(nbComputedTask != inGraph.getNbNodes()){
            {
                assert(workers.size());
                Worker worker = workers.top();
                workers.pop();

                const int currentTaskId = worker.currentTaskId;
                assert(currentTime <= worker.scheduledAvailable);
                currentTime = worker.scheduledAvailable;

                // release dependencies
                for(const auto& successorNode : inGraph.getNode(currentTaskId)->getSuccessors()){
                    countPredecessorsOver[successorNode->getId()] += 1;
                    if(countPredecessorsOver[successorNode->getId()] == successorNode->getPredecessors().size()){
                        readyTasks.push_back(successorNode->getId());
                        // TODO inAllTasks[successorId].setReadyTime(currentTime);
                    }
                }

                // Make worker available again
                idleWorkerCount.push_back(worker.currentWorkerId);
            }

            while(readyTasks.size() && idleWorkerCount.size()){
                const int readyTaskId = readyTasks.back();
                readyTasks.pop_back();

                const int workerId = idleWorkerCount.back();
                idleWorkerCount.pop_back();

                Worker wk{currentTime + inGraph.getNode(readyTaskId)->getCost(),
                         readyTaskId,
                         workerId};
                workers.push(wk);

                events.push_back(Event(workerId, readyTaskId, inGraph.getNode(readyTaskId)->getCost(), currentTime));

                nbComputedTask += 1;
            }

            assert(workers.size() != 0);
        }

        assert(workers.size() != 0);
        while(workers.size() != 0){
            assert(workers.size());
            Worker worker = workers.top();
            workers.pop();

            const int currentTaskId = worker.currentTaskId;
            assert(currentTime <= worker.scheduledAvailable);
            currentTime = worker.scheduledAvailable;

            assert(inGraph.getNode(currentTaskId)->getSuccessors().size() == 0);
        }


        std::cout << "[EXECUTION] Total duration => " << currentTime << std::endl;

        return std::make_tuple(currentTime, events);
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

    int duration;
    std::vector<Executor::Event> events;
    std::tie(duration, events) = Executor::Execute(depGraph, 2);

    return 0;
}
