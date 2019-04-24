#include <iostream>
#include <functional>

#include "graph.hpp"
#include "executor.hpp"
#include "generator.hpp"
#include "utils.hpp"


int main(int argc, char** argv){
    std::vector<std::string> params;
    params.insert(params.end(), argv, argv+argc);

    const std::string helpContent = "[HELP] You must pass the graph generation method and a size.\n"
                                    "[HELP] $ ./main [generation method] [size]\n"
                                    "[HELP] Where generation method is among: tree, deptree, 2dgrid, doubletree\n";

    if(params.size() != 3){
        std::cout << "[ERROR] Invalid number of parameters.\n" << helpContent;
        return 1;
    }

    //////////////////////////////////////////////////////////////////////////

    std::pair<int, std::vector<std::pair<int,int>>> someDeps;

    const std::vector<std::string> methodNames{"tree", "deptree", "2dgrid", "doubletree"};
    const size_t choice = std::distance(methodNames.begin(), std::find(methodNames.begin(), methodNames.end(), params[1]));

    const int size = Utils::StrToNum<int>(params[2], -1);
    if(size == -1){
        std::cout << "[ERROR] cannot convert size for str: " << params[2] << std::endl;
        return 1;
    }

    switch(choice){
    case 0:
        std::cout << "[INFO] use tree\n";
        someDeps = Generator::GenerateDepTreeTasks(size);
        break;
    case 1:
        std::cout << "[INFO] use deptree\n";
        someDeps = Generator::GenerateDoubleDepTreeTasks(size);
        break;
    case 2:
        std::cout << "[INFO] use 2dgrid\n";
        someDeps = Generator::Generate2DGrid(size);
        break;
    case 3:
        std::cout << "[INFO] use doubletree\n";
        someDeps = Generator::GenerateDoubleDepTreeTasks(size);
        break;
    default:
        std::cout << "[ERROR] invalid method name " << params[1] << std::endl;
        return 2;
    }

    Graph aGraph(someDeps.first, someDeps.second);
    std::cout << "Number of nodes : " << aGraph.getNbNodes() << "\n";
    assert(aGraph.isDag());
    std::pair<int,double> degGraph = aGraph.estimateDegreeOfParallelism();
    std::cout << "Degree of parallelism one the " << methodNames[choice] + "-" + std::to_string(size) << " graph : " << degGraph.first << "  " << degGraph.second << "\n";
    aGraph.saveToDot("/tmp/agraph-" + methodNames[choice] + "-" + std::to_string(size) + ".dot");

    return 0;
}
