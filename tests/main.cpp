#include <list>
#include <cassert>
#include <algorithm>

class Node{
    std::list<Node*> successors;
    std::list<Node*> predecessors;

    int id;

public:
    explicit Node(const int inId) : id(inId){
    }

    int getId() const{
        return  id;
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

};

#include <iostream>

int main(){
    std::vector<std::pair<int,int>> someDeps{{0,1}, {0,2}, {1,2}};
    Graph aGraph(someDeps);

    return 0;
}
