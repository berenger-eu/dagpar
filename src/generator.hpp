// Please read the corresponding licence file
#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <utility>
#include <set>

class Generator{
public:
    ///////////////////////////////////////////////////////////
    ///
    /// Methods to generate some graphs
    ///
    ///////////////////////////////////////////////////////////

    static std::pair<int, std::vector<std::pair<int,int>>> GenerateBinaryTreeTasks(const int inHeight){
        std::vector<std::pair<int,int>> someDeps;
        int cellsOffset = 0;
        for(int level = 1 ; level < inHeight ; ++level){
            const int nbCellsAtPreviousLevel = (1 << (level-1));
            const int nbCellsAtLevel = (1 << level);
            for(int idxCell = 0 ; idxCell < nbCellsAtLevel ; ++idxCell){
                someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell/2, cellsOffset + nbCellsAtPreviousLevel + idxCell});
            }
            cellsOffset += nbCellsAtPreviousLevel;
        }
        return std::pair<int, std::vector<std::pair<int,int>>>(someDeps.back().second+1,someDeps);
    }

    static std::pair<int, std::vector<std::pair<int,int>>> GenerateDepTreeTasks(const int inHeight){
        std::vector<std::pair<int,int>> someDeps;
        int cellsOffset = 0;
        for(int level = 1 ; level < inHeight ; ++level){
            const int nbCellsAtPreviousLevel = level;
            const int nbCellsAtLevel = (level+1);

            someDeps.push_back(std::pair<int,int>{cellsOffset, cellsOffset + nbCellsAtPreviousLevel});

            for(int idxCell = 1 ; idxCell < nbCellsAtLevel-1 ; ++idxCell){
                someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell-1, cellsOffset + nbCellsAtPreviousLevel + idxCell});
                someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell, cellsOffset + nbCellsAtPreviousLevel + idxCell});
            }

            someDeps.push_back(std::pair<int,int>{cellsOffset + nbCellsAtPreviousLevel-1, cellsOffset + nbCellsAtPreviousLevel + nbCellsAtLevel-1});

            cellsOffset += nbCellsAtPreviousLevel;
        }
        return std::pair<int, std::vector<std::pair<int,int>>>(someDeps.back().second+1,someDeps);
    }

    static std::pair<int, std::vector<std::pair<int,int>>> GenerateDoubleDepTreeTasks(const int inHeight){
        std::vector<std::pair<int,int>> someDeps;
        int cellsOffset = 0;
        for(int level = 1 ; level < (inHeight+1)/2 ; ++level){
            const int nbCellsAtPreviousLevel = level;
            const int nbCellsAtLevel = (level+1);

            someDeps.push_back(std::pair<int,int>{cellsOffset, cellsOffset + nbCellsAtPreviousLevel});

            for(int idxCell = 1 ; idxCell < nbCellsAtLevel-1 ; ++idxCell){
                someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell-1, cellsOffset + nbCellsAtPreviousLevel + idxCell});
                someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell, cellsOffset + nbCellsAtPreviousLevel + idxCell});
            }

            someDeps.push_back(std::pair<int,int>{cellsOffset + nbCellsAtPreviousLevel-1, cellsOffset + nbCellsAtPreviousLevel + nbCellsAtLevel-1});

            cellsOffset += nbCellsAtPreviousLevel;
        }

        if((inHeight & 1) == 0){
            const int nbCellsAtPreviousLevel = inHeight/2;
            const int nbCellsAtLevel = inHeight/2;

            for(int idxCell = 0 ; idxCell < nbCellsAtLevel ; ++idxCell){
                someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell, cellsOffset + nbCellsAtPreviousLevel + idxCell});
            }

            cellsOffset += nbCellsAtPreviousLevel;
        }

        for(int level = (inHeight+2)/2 ; level < inHeight ; ++level){
            const int nbCellsAtPreviousLevel = (inHeight+1) - level;
            const int nbCellsAtLevel = (inHeight+1) - level - 1;

            for(int idxCell = 0 ; idxCell < nbCellsAtLevel ; ++idxCell){
                someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell, cellsOffset + nbCellsAtPreviousLevel + idxCell});
                someDeps.push_back(std::pair<int,int>{cellsOffset + idxCell+1, cellsOffset + nbCellsAtPreviousLevel + idxCell});
            }

            cellsOffset += nbCellsAtPreviousLevel;
        }

        return std::pair<int, std::vector<std::pair<int,int>>>(someDeps.back().second+1,someDeps);
    }

    static std::pair<int, std::vector<std::pair<int,int>>> Generate2DGrid(const int inGridDim){
        std::vector<std::pair<int,int>> someDeps;
        for(int idxRow = 1 ; idxRow < inGridDim ; ++idxRow){
            for(int idxCol = 1 ; idxCol < inGridDim ; ++idxCol){
                someDeps.push_back(std::pair<int,int>{(idxCol*inGridDim)+idxRow-1, (idxCol*inGridDim)+idxRow});
                someDeps.push_back(std::pair<int,int>{((idxCol-1)*inGridDim)+idxRow, (idxCol*inGridDim)+idxRow});
            }
        }
        return std::pair<int, std::vector<std::pair<int,int>>>(inGridDim*inGridDim,someDeps);
    }

    static std::pair<int, std::vector<std::pair<int,int>>> LoadEdgeFile(const std::string& inFilename){
        std::ifstream edgeFile(inFilename);

        if(edgeFile.is_open() == false){
            std::cout << "[ERROR] Cannot load file " << inFilename << std::endl;
            std::cout << "[ERROR] Return empty graph" << std::endl;
            return std::pair<int, std::vector<std::pair<int,int>>>();
        }

        if(inFilename.length() >= 4 && inFilename.substr(inFilename.length() - 4) == ".dot"){
            std::pair<int, std::vector<std::pair<int,int>>> edges;
            int maxNodeId = -1;
            int minNodeId = std::numeric_limits<int>::max();

            std::string line;
            while (std::getline(edgeFile, line)){
                const auto arrowPos = line.find("->");
                if(arrowPos != std::string::npos){
                    std::istringstream iss(line);
                    char arrowPart1, arrowPart2;
                    int edgeSrc, edgeDst;
                    if (!(iss >> edgeSrc >> arrowPart1 >> arrowPart2 >> edgeDst)) {
                        std::cout << "[ERROR] Bad line " <<line << std::endl;
                        return std::pair<int, std::vector<std::pair<int,int>>>();
                    }

                    edges.second.emplace_back(edgeSrc, edgeDst);
                    maxNodeId = std::max(maxNodeId, std::max(edgeSrc, edgeDst));
                    minNodeId = std::min(minNodeId, std::min(edgeSrc, edgeDst));
                }
            }

            for(auto& pair : edges.second){
                pair.first -= minNodeId;
                pair.second -= minNodeId;
            }

            edges.first = maxNodeId+1-minNodeId;

            // Remove duplicate
            std::set<std::pair<int,int>> alreadyExist;
            int iter = 0;
            while(iter != int(edges.second.size())){
                if(alreadyExist.find(edges.second[iter]) == alreadyExist.end()){
                    alreadyExist.insert(edges.second[iter]);
                    iter += 1;
                }
                else{
                    std::swap(edges.second[iter], edges.second.back());
                    edges.second.pop_back();
                }
            }

            return edges;
        }
        else{
            std::pair<int, std::vector<std::pair<int,int>>> edges;
            int nbEdges;
            edgeFile >> nbEdges;

            std::cout << "[INFO] There are " << nbEdges << " edges" << std::endl;
            edges.second.reserve(nbEdges);

            int maxNodeId = -1;

            int edgeIdx, edgeSrc, edgeDst;
            while (edgeFile >> edgeIdx >> edgeSrc >> edgeDst){
                if(edgeIdx != int(edges.second.size())){
                    std::cout << "[ERROR] Bad idx, is " << edgeIdx << " should be " << edges.second.size() << std::endl;
                    std::cout << "[ERROR] Return empty graph" << std::endl;
                    return std::pair<int, std::vector<std::pair<int,int>>>();
                }
                if(nbEdges <= edgeIdx){
                    std::cout << "[ERROR] Bad idx, is " << edgeIdx << " which is largest than the expected number of edges" << std::endl;
                    std::cout << "[ERROR] Return empty graph" << std::endl;
                    return std::pair<int, std::vector<std::pair<int,int>>>();
                }

                edges.second.emplace_back(edgeSrc, edgeDst);
                maxNodeId = std::max(maxNodeId, std::max(edgeSrc, edgeDst));
            }

            edges.first = maxNodeId+1;

            return edges;
        }
    }


    static std::vector<double> GetCostIfExist(const std::string& inFilename){
        std::ifstream edgeFile(inFilename);

        if(edgeFile.is_open() == false){
            std::cout << "[ERROR] Cannot load file " << inFilename << std::endl;
            std::cout << "[ERROR] Return empty graph" << std::endl;
            return std::vector<double>();
        }

        if(inFilename.length() >= 4 && inFilename.substr(inFilename.length() - 4) == ".dot"){
            std::vector<double> costs;

            int minVecticeId = std::numeric_limits<int>::max();

            std::string line;
            while (std::getline(edgeFile, line)){
                const auto arrowPos = line.find("->");
                if(arrowPos == std::string::npos){
                    const auto sizePos = line.find("[size =\"");
                    if(sizePos != std::string::npos){
                        int verticeId;
                        if (!(std::istringstream(line) >> verticeId)) {
                            std::cout << "[ERROR] Bad line " << line << std::endl;
                            return std::vector<double>();
                        }

                        double cost;
                        if (!(std::istringstream(line.substr(arrowPos+7)) >> cost)) {
                            std::cout << "[ERROR] Bad line " <<line << std::endl;
                            return std::vector<double>();
                        }

                        if(int(costs.size()) <= verticeId){
                            costs.resize(verticeId+1, 0);
                        }
                        costs[verticeId] = cost;

                        minVecticeId = std::min(minVecticeId, verticeId);
                    }
                }
            }

            for(int idx = minVecticeId ; idx < int(costs.size()) ; ++idx){
                costs[idx-minVecticeId] = costs[idx];
            }
            costs.resize(std::max(0,int(costs.size()) - minVecticeId));

            return costs;
        }
        return std::vector<double>();
    }
};


#endif
