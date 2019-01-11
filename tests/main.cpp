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
    static std::string ReplaceAllInString(std::string sentence, const std::string& itemToReplace, const std::string& inSubstitutionString){
        for(std::string::size_type idxFound = sentence.find(itemToReplace); idxFound != std::string::npos; idxFound = sentence.find(itemToReplace)){
            sentence.replace(idxFound, itemToReplace.size(), inSubstitutionString);
        }
        return sentence;
    }

public:
    class Event{
        int workerId;
        int nodeId;
        int cost;
        int startingPoint;
        int readyPoint;

    public:
        Event(){}

        void init(const int inNodeId, const int inCost){
            nodeId = inNodeId;
            cost = inCost;
        }

        void startTask(const int inWorkerId, const int inStartingPoint){
            workerId = inWorkerId;
            startingPoint = inStartingPoint;
        }

        void becomeReady(const int inReadyPoint){
            readyPoint = inReadyPoint;
        }

        // Use the same GET method names as the ones used in the heteroprio emulator

        int getThreadIdComputer() const{
            return workerId;
        }

        int getId() const{
            return nodeId;
        }

        int getDuration() const{
            return cost;
        }

        int getStartingTime() const{
            return startingPoint;
        }

        int getReadyTime() const{
            return readyPoint;
        }

        int getEndingTime() const{
            return getDuration() + getStartingTime();
        }

        int getCreationTime() const{
            return 0;
        }

        std::string getTaskName() const{
            return "Node " + std::to_string(nodeId);
        }
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
        std::vector<Event> events;
        events.resize(inGraph.getNbNodes());

        for(int idxNode = 0 ; idxNode < inGraph.getNbNodes() ; ++idxNode){
            const auto& node = inGraph.getNode(idxNode);
            assert(idxNode == node->getId());

            events[node->getId()].init(node->getId(), node->getCost());

            if(node->getPredecessors().size() == 0){
                readyTasks.push_back(node->getId());
                events[node->getId()].becomeReady(0);
            }
        }

        std::cout << "[EXECUTION] Ready tasks at starting time => " << readyTasks.size() << std::endl;
        assert(readyTasks.size());

        std::vector<int> countPredecessorsOver(inGraph.getNbNodes(), 0);
        std::priority_queue<Worker> workers;
        int nbComputedTask = 0;
        int currentTime = 0;

        while(readyTasks.size() && idleWorkerCount.size()){
            const int readyTaskId = readyTasks.back();
            readyTasks.pop_back();

            const int workerId = idleWorkerCount.back();
            idleWorkerCount.pop_back();

            Worker wk{currentTime + inGraph.getNode(readyTaskId)->getCost(),
                     readyTaskId,
                     workerId};
            workers.push(wk);

            events[readyTaskId].startTask(workerId, currentTime);

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
                        events[successorNode->getId()].becomeReady(currentTime);
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

                events[readyTaskId].startTask(workerId, currentTime);

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



    static void EventsToTrace(const std::string& outputFilename, const std::vector<Event>& tasksFinished, const int nbWorkers) {
        using time_type = int;
        const int startingTime = 0;

        std::ofstream svgfile(outputFilename);

        if(svgfile.is_open() == false){
            throw std::invalid_argument("Cannot open filename : " + outputFilename);
        }

        const long int vsizeperthread = std::max(100, std::min(200, 2000/nbWorkers));
        const long int vmargin = 100;
        const long int threadstrock = 5;
        const long int vmarginthreadstats = 10;
        const long int vdim = nbWorkers * (vsizeperthread+threadstrock) + 2*vmargin + vmarginthreadstats + 2*vsizeperthread;

        time_type duration = 0;
        for(const auto& atask : tasksFinished){
            duration = std::max(duration, atask.getEndingTime()) - startingTime;
        }

        const long int hdimtime = std::max(6000, int(log(duration)*200));
        const long int hmargin = 50;
        const long int hdim = hdimtime + 2*hmargin;

        svgfile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
        svgfile << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" << hdim << "\" height=\"" << vdim << "\">\n";
        svgfile << "  <title>Execution trace</title>\n";
        svgfile << "  <desc>\n";
        svgfile << "    Spetabaru traces for " << tasksFinished.size() << " tasks\n";
        svgfile << "  </desc>\n";
        svgfile << "\n";
        // Back
        svgfile << "  <rect width=\"" << hdim << "\" height=\"" << vdim << "\" x=\"0\" y=\"0\" fill=\"white\" />\n";


        svgfile << "  <line x1=\"" << hmargin << "\" y1=\"" << vmargin-10
                << "\" x2=\"" << hdimtime+hmargin << "\" y2=\"" << vmargin-10 << "\" stroke=\"" << "black" <<"\" stroke-width=\"" << "1" <<"\"/>\n";

        svgfile << "  <circle cx=\"" << hmargin << "\" cy=\"" << vmargin-10 << "\" r=\"" << "6" << "\" fill=\"" << "black" <<"\" />\n";
        svgfile << "  <circle cx=\"" << hdimtime+hmargin << "\" cy=\"" << vmargin-10 << "\" r=\"" << "6" << "\" fill=\"" << "black" << "\" />\n";

        for(time_type second = 1 ; second < duration ; second += 1.){
            svgfile << "  <line x1=\"" << static_cast<long  int>(double(hdimtime)*second/double(duration))+hmargin  << "\" y1=\"" << vmargin-14
                    << "\" x2=\"" << static_cast<long  int>(double(hdimtime)*second/double(duration))+hmargin << "\" y2=\"" << vmargin-6 << "\" stroke=\"" << "black" <<"\" stroke-width=\"" << "1" <<"\"/>\n";
        }
        for(time_type second10 = 10 ; second10 < duration ; second10 += 10.){
            svgfile << "  <line x1=\"" << static_cast<long  int>(double(hdimtime)*second10/double(duration))+hmargin  << "\" y1=\"" << vmargin-18
                    << "\" x2=\"" << static_cast<long  int>(double(hdimtime)*second10/double(duration))+hmargin << "\" y2=\"" << vmargin-2 << "\" stroke=\"" << "black" <<"\" stroke-width=\"" << "1" <<"\"/>\n";
        }

        {
            const std::string label = "Total time = " + std::to_string(duration) + " s";
            svgfile << "<text x=\"" << hmargin+hdimtime/2-int(label.size())*5 << "\" y=\"" << vmargin-20 <<
                   "\" font-size=\"30\" fill=\"black\">" << label << "</text>\n";
        }

        for(int idxThread = 0 ; idxThread < nbWorkers ; ++idxThread){
            svgfile << "  <rect width=\"" << hdimtime << "\" height=\"" << vsizeperthread
                    << "\" x=\"" << hmargin << "\" y=\"" << idxThread*(vsizeperthread+threadstrock) + vmargin << "\" style=\"fill:white;stroke:black;stroke-width:" << threadstrock << "\" />\n";

            const std::string label = "Worker " + std::to_string(idxThread);

            svgfile << "<text x=\"" << hmargin/2 << "\" y=\"" << idxThread*(vsizeperthread+threadstrock) + vmargin + label.size()*30 - vsizeperthread/2 <<
                       "\" font-size=\"30\" fill=\"black\" transform=\"rotate(-90, " << hmargin/2 << " "
                    << idxThread*(vsizeperthread+threadstrock) + vmargin + label.size()*30 - vsizeperthread/2 << ")\">" << label << "</text>\n";
        }

        std::unordered_map<std::string, std::string> colors;

        for(const auto& atask : tasksFinished){
            const long int idxThreadComputer = atask.getThreadIdComputer() + 1;
            assert(idxThreadComputer != -1);
            assert(0 < idxThreadComputer && idxThreadComputer <= nbWorkers);
            const long int ypos_start = (idxThreadComputer-1)*(vsizeperthread+threadstrock) + threadstrock/2 + vmargin;
            const long int ypos_end = ypos_start + vsizeperthread - threadstrock;
            const time_type taskStartTime = atask.getStartingTime() - startingTime;
            const time_type taskEndTime = atask.getEndingTime() - startingTime;
            const long int xpos_start = static_cast<long  int>(double(hdimtime)*taskStartTime/double(duration)) + hmargin;
            const long int xpos_end = std::max(xpos_start+1,static_cast<long  int>(double(hdimtime)*taskEndTime/double(duration)) + hmargin);
            const long int taskstrocke = 2;

            std::string strForColor = atask.getTaskName();

            const std::size_t ddpos = strForColor.find("--");
            if(ddpos != std::string::npos){
                strForColor = strForColor.substr(0, ddpos);
            }
            if(atask.getTaskName().length() && atask.getTaskName().at(atask.getTaskName().length()-1) == '\''){
                strForColor.append("'");
            }

            if(colors.find(strForColor) == colors.end()){
                size_t hashname =  std::hash<std::string>()(strForColor);

                std::mt19937 rng(hashname);
                std::uniform_real_distribution<double> uni(0,1);

                const int colorR = int(uni(rng)*200) + 50;
                const int colorG = int(uni(rng)*200) + 50;
                const int colorB = int(uni(rng)*200) + 50;
                colors[strForColor] = std::string("rgb(")
                        + std::to_string(colorR) + ","
                        + std::to_string(colorG) + ","
                        + std::to_string(colorB) + ")";
            }

            svgfile << "<g>\n";
            svgfile << "    <title id=\"" << atask.getId() << "\">" << ReplaceAllInString(ReplaceAllInString(ReplaceAllInString(atask.getTaskName(),"&", " REF "),"<","["),">","]")
                    << " -- Duration " << atask.getEndingTime() - atask.getStartingTime() << "s"
                    << " -- Enable = " << "TRUE" << "</title>\n";

            svgfile << "    <rect width=\"" << xpos_end-xpos_start << "\" height=\"" << ypos_end-ypos_start
                    << "\" x=\"" << xpos_start << "\" y=\"" << ypos_start
                    << "\" style=\"fill:" << colors[strForColor] << ";stroke:black" << ";stroke-width:" << taskstrocke << "\" />\n";

            svgfile << "</g>\n";
        }

        const long int offsetStat = nbWorkers*(vsizeperthread+threadstrock) + vmarginthreadstats;
        const char* statsNames[] = {"Submited", "Ready"};
        for(int idxStat = 0 ; idxStat < 2 ; ++idxStat){
            svgfile << "  <rect width=\"" << hdimtime << "\" height=\"" << vsizeperthread
                    << "\" x=\"" << hmargin << "\" y=\"" << offsetStat + idxStat*(vsizeperthread+threadstrock) + vmargin << "\" style=\"fill:white;stroke:black;stroke-width:" << threadstrock << "\" />\n";

            const std::string label = statsNames[idxStat];

            svgfile << "<text x=\"" << hmargin/2 << "\" y=\"" << offsetStat + idxStat*(vsizeperthread+threadstrock+50) + vmargin + label.size()*30 - vsizeperthread/2 <<
                       "\" font-size=\"30\" fill=\"black\" transform=\"rotate(-90, " << hmargin/2 << " "
                    << offsetStat + idxStat*(vsizeperthread+threadstrock+50) + vmargin + label.size()*30 - vsizeperthread/2 << ")\">" << label << "</text>\n";
        }

        std::vector<int> nbReady(hdimtime, 0);
        std::vector<int> nbSubmited(hdimtime, 0);

        for(const auto& atask : tasksFinished){
            const time_type taskSubmitedTime = atask.getCreationTime() - startingTime;
            const time_type taskReadyTime = atask.getReadyTime() - startingTime;
            const time_type taskStartTime = atask.getStartingTime() - startingTime;
            const long int xpos_submited = static_cast<long  int>(double(hdimtime)*taskSubmitedTime/double(duration));
            const long int xpos_ready = static_cast<long  int>(double(hdimtime)*taskReadyTime/double(duration));
            const long int xpos_start = static_cast<long  int>(double(hdimtime)*taskStartTime/double(duration));

            nbSubmited[xpos_submited] += 1;

            nbReady[xpos_ready] += 1;
            if(xpos_ready != xpos_start || xpos_ready == hdimtime-1){
                nbReady[xpos_start] -= 1;
            }
            else{
                nbReady[xpos_start+1] -= 1;
            }
        }

        int maxReady = 0;
        {
            int currentStat = 0;
            for(int idx = 0 ; idx < hdimtime ; ++idx){
                currentStat += nbReady[idx];
                maxReady = std::max(maxReady , currentStat);
            }
        }

        const std::reference_wrapper<const std::vector<int>> statVal[] = {nbSubmited, nbReady};
        const int maxStatVal[2] = {static_cast<int>(tasksFinished.size()), maxReady};

        for(int idxStat = 0 ; idxStat < 2 ; ++idxStat){
             svgfile << "<polyline points=\"";
             //20,20 40,25 60,40 80,120 120,140 200,180"
             const int maxStat = maxStatVal[idxStat];
             int currentStat = 0;
             for(int idx = 0 ; idx < hdimtime ; ++idx){
                 currentStat += statVal[idxStat].get()[idx];
                 const long int xpos = hmargin+idx;
                 const long int ypos = offsetStat + idxStat*(vsizeperthread+threadstrock) + vmargin
                                       + vsizeperthread
                                       - static_cast<long int>(double(vsizeperthread-threadstrock)*double(currentStat)/double(maxStat))
                                       - threadstrock/2;
                 svgfile << xpos << "," << ypos << " ";
             }
             svgfile << "\" style=\"fill:none;stroke:rgb(112,0,0);stroke-width:3\" />\n";
        }

        svgfile << "<text x=\"" << hmargin + hdimtime + 10
                << "\" y=\"" << offsetStat + 0*(vsizeperthread+threadstrock) + vmargin + 15 <<
                   "\" font-size=\"30\" fill=\"black\">" << tasksFinished.size() << "</text>\n";
        svgfile << "<text x=\"" << hmargin + hdimtime + 10
                << "\" y=\"" << offsetStat + 1*(vsizeperthread+threadstrock) + vmargin + 15 <<
                   "\" font-size=\"30\" fill=\"black\">" << maxReady << "</text>\n";

        svgfile << "</svg>\n";
    }


};



#include <iostream>

int main(){
    const int nbThreads = 2;
    std::vector<std::pair<int,int>> someDeps{{0,1}, {0,2}, {1,2}};
    Graph aGraph(someDeps);
    aGraph.partition(1,1,nbThreads);
    aGraph.saveToDot("/tmp/agraph.dot");

    Graph depGraph = aGraph.getPartitionGraph();
    aGraph.saveToDot("/tmp/depgraph.dot");

    int duration;
    std::vector<Executor::Event> events;
    std::tie(duration, events) = Executor::Execute(depGraph, nbThreads);
    Executor::EventsToTrace("/tmp/dep-graph-" + std::to_string(nbThreads) + "trace.svg", events, nbThreads);

    return 0;
}
