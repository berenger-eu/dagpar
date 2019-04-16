// Please read the corresponding licence file
#ifndef NODE_HPP
#define NODE_HPP

#include <vector>
#include <cassert>
#include <algorithm>

///////////////////////////////////////////////////////////
///
/// A node which has links to predecessors and successors
///
///////////////////////////////////////////////////////////
class Node{
    std::vector<Node*> successors;
    std::vector<Node*> predecessors;

    int id;
    int partitionId;
    double cost;

public:
    explicit Node(const int inId, const int inPartitionId = 0)
        : id(inId), partitionId(inPartitionId), cost(1){
    }

    void setPartitionId(const int inPartitionId){
        partitionId = inPartitionId;
    }

    void setCost(const double inCost){
        cost = inCost;
    }

    int getId() const{
        return  id;
    }

    int getPartitionId() const{
        return partitionId;
    }

    double getCost() const{
        return cost;
    }

    void addSuccessor(Node* next){
        //Is not true when index are not in order: assert(id < next->id);
        assert(next != this);
        assert(next->id != id);
        assert(std::find(successors.begin(), successors.end(), next) == successors.end());
        successors.push_back(next);
    }

    void addPredecessor(Node* prev){
        //Is not true when index are not in order: assert(prev->id < id);
        assert(prev != this);
        assert(prev->id != id);
        assert(std::find(predecessors.begin(), predecessors.end(), prev) == predecessors.end());
        predecessors.push_back(prev);
    }

    const std::vector<Node*>& getSuccessors() const{
        return successors;
    }

    const std::vector<Node*>& getPredecessors() const{
        return predecessors;
    }
};

#endif
