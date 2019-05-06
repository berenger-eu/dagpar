// Please read the corresponding licence file
#ifndef EXECUTOR_HPP
#define EXECUTOR_HPP

#include <tuple>
#include <iostream>
#include <queue>
#include <istream>
#include <random>

#include "graph.hpp"

///////////////////////////////////////////////////////////
///
/// To emulate the execution of a graph
///
///////////////////////////////////////////////////////////
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
        double cost;
        double startingPoint;
        double readyPoint;
        bool isWork;

    public:
        Event(){}

        void init(const int inNodeId, const double inCost, const bool inIsWork){
            nodeId = inNodeId;
            cost = inCost;
            isWork = inIsWork;
        }

        bool getIsWork() const{
            return isWork;
        }

        void startTask(const int inWorkerId, const double inStartingPoint){
            workerId = inWorkerId;
            startingPoint = inStartingPoint;
        }

        void becomeReady(const double inReadyPoint){
            readyPoint = inReadyPoint;
        }

        // Use the same GET method names as the ones used in the heteroprio emulator

        int getThreadIdComputer() const{
            return workerId;
        }

        int getId() const{
            return nodeId;
        }

        double getDuration() const{
            return cost;
        }

        double getStartingTime() const{
            return startingPoint;
        }

        double getReadyTime() const{
            return readyPoint;
        }

        double getEndingTime() const{
            return getDuration() + getStartingTime();
        }

        double getCreationTime() const{
            return 0;
        }

        std::string getTaskName() const{
            return "Node " + std::to_string(nodeId);
        }
    };

    static std::tuple<double,std::vector<Event>> Execute(const Graph& inGraph, const int inNbWorkers,
                                                      const double overheadPerTask = 0,
                                                      const double popOverhead = 0,
                                                      const double pushOverhead = 0){
        if(inGraph.getNbNodes() == 0){
            return std::make_tuple(0,std::vector<Event>());
        }

        struct Worker{
            double scheduledAvailable;
            int currentTaskId;
            int currentWorkerId;
            bool operator<(const Worker& other) const{
                return scheduledAvailable > other.scheduledAvailable;
            }
        };

        std::vector<int> idleWorkerCount;

        idleWorkerCount.resize(inNbWorkers);
        std::iota(idleWorkerCount.begin(), idleWorkerCount.end(), 0);

        std::vector<Event> events;
        events.resize(inGraph.getNbNodes());

        std::vector<Event> pushEvents;
        pushEvents.resize(inGraph.getNbNodes());

        std::vector<Event> popEvents;
        popEvents.resize(inGraph.getNbNodes());

        std::vector<int> readyTasks;

        double currentTime = 0;

        for(int idxNode = 0 ; idxNode < inGraph.getNbNodes() ; ++idxNode){
            const auto& node = inGraph.getNodeFromId(idxNode);
            assert(idxNode == node->getId());

            events[node->getId()].init(node->getId(), node->getCost() + overheadPerTask, true);
            pushEvents[node->getId()].init(node->getId(), pushOverhead, false);
            popEvents[node->getId()].init(node->getId(), popOverhead, false);

            if(node->getPredecessors().size() == 0){
                pushEvents[node->getId()].becomeReady(currentTime);
                pushEvents[node->getId()].startTask(0, currentTime);

                currentTime+= pushOverhead;
                events[node->getId()].becomeReady(currentTime);
                readyTasks.push_back(node->getId());
            }
        }

        std::cout << "[EXECUTION] Ready tasks at starting time => " << readyTasks.size() << std::endl;
        assert(readyTasks.size());

        std::vector<int> countPredecessorsOver(inGraph.getNbNodes(), 0);
        std::priority_queue<Worker> workers;
        int nbComputedTask = 0;

        while(readyTasks.size() && idleWorkerCount.size()){
            const int readyTaskId = readyTasks.back();
            readyTasks.pop_back();

            const int workerId = idleWorkerCount.back();
            idleWorkerCount.pop_back();

            popEvents[readyTaskId].becomeReady(currentTime);
            popEvents[readyTaskId].startTask(workerId, currentTime);

            currentTime += popOverhead;

            Worker wk{currentTime + inGraph.getNodeFromId(readyTaskId)->getCost() + overheadPerTask,
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
                currentTime = std::max(currentTime, worker.scheduledAvailable);

                // release dependencies
                for(const auto& successorNode : inGraph.getNodeFromId(currentTaskId)->getSuccessors()){
                    countPredecessorsOver[successorNode->getId()] += 1;
                    if(countPredecessorsOver[successorNode->getId()] == int(successorNode->getPredecessors().size())){                        
                        pushEvents[successorNode->getId()].becomeReady(currentTime);
                        pushEvents[successorNode->getId()].startTask(worker.currentWorkerId, currentTime);

                        currentTime+= pushOverhead;

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

                popEvents[readyTaskId].becomeReady(currentTime);
                popEvents[readyTaskId].startTask(workerId, currentTime);

                currentTime += popOverhead;

                Worker wk{currentTime + inGraph.getNodeFromId(readyTaskId)->getCost() + overheadPerTask,
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

            currentTime = std::max(currentTime, worker.scheduledAvailable);

            assert(inGraph.getNodeFromId(worker.currentTaskId)->getSuccessors().size() == 0);
        }


        std::cout << "[EXECUTION] Total duration => " << currentTime << std::endl;

        events.insert(events.end(), pushEvents.begin(), pushEvents.end());
        events.insert(events.end(), popEvents.begin(), popEvents.end());

        return std::make_tuple(currentTime, events);
    }



    static void EventsToTrace(const std::string& outputFilename, const std::vector<Event>& tasksFinished, const int nbWorkers) {
        using time_type = double;
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

        for(time_type second = 1 ; second < duration ; second += time_type(1)){
            svgfile << "  <line x1=\"" << static_cast<long  int>(double(hdimtime)*second/double(duration))+hmargin  << "\" y1=\"" << vmargin-14
                    << "\" x2=\"" << static_cast<long  int>(double(hdimtime)*second/double(duration))+hmargin << "\" y2=\"" << vmargin-6 << "\" stroke=\"" << "black" <<"\" stroke-width=\"" << "1" <<"\"/>\n";
        }
        for(time_type second10 = 10 ; second10 < duration ; second10 += time_type(10)){
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
            if(atask.getIsWork()){
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

#endif
