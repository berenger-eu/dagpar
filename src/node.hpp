// Please read the corresponding licence file
#ifndef NODE_HPP
#define NODE_HPP

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

#endif
