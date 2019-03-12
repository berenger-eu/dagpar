// Please read the corresponding licence file
#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <utility>

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
};


#endif
