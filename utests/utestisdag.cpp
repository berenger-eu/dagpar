#include "generator.hpp"
#include "graph.hpp"

#include "UTester.hpp"

#include <vector>

class TestIsDag : public UTester< TestIsDag > {
    using Parent = UTester< TestIsDag >;

    void CoreTest(){
        {
            std::vector<std::pair<int,int>> edges;
            edges.emplace_back(0, 1);
            edges.emplace_back(1, 0);

            Graph aGraph(2, edges);
            UASSERTETRUE(aGraph.isDag() == false);
        }
        {
            std::vector<std::pair<int,int>> edges;
            edges.emplace_back(0, 1);
            edges.emplace_back(1, 2);
            edges.emplace_back(2, 3);
            edges.emplace_back(3, 0);

            Graph aGraph(4, edges);
            UASSERTETRUE(aGraph.isDag() == false);
        }
        {
            std::vector<std::pair<int,int>> edges;
            edges.emplace_back(0, 1);
            edges.emplace_back(1, 2);
            edges.emplace_back(2, 3);
            edges.emplace_back(3, 2);

            Graph aGraph(4, edges);
            UASSERTETRUE(aGraph.isDag() == false);
        }
        {
            std::vector<std::pair<int,int>> edges;
            edges.emplace_back(0, 1);
            edges.emplace_back(1, 2);
            edges.emplace_back(2, 3);

            Graph aGraph(4, edges);
            UASSERTETRUE(aGraph.isDag() == true);
        }
    }

    void SetTests() {
        Parent::AddTest(&TestIsDag::CoreTest, "Test everything");
    }
};

// You must do this
TestClass(TestIsDag)

