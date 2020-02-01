This Software is a test project to study the acyclic partitioning of DAG (of tasks), see the related paper at https://peerj.com/articles/cs-247/

The official repository is https://gitlab.inria.fr/bramas/dagpar/

# Compilation

Simply go in the build dir and use cmake as usual:
```bash
mkdir build
cd build
cmake ..
make
```

If everything is fine, the tests should be correct:
```bash
make test
```

# Adding a test file

Any `.cpp` file added in `tests` will be automatically compiled after using `cmake ..` and `make` in the build dir.

# Where to make modifications

The `partition` method of the `Graph` class should do the partitioning.
Anyone interested to improve the method should modify it, or duplicate an existing method and look at the obtained results.

# Running the example

One can execute the `main` binary (from the `tests/main.cpp` file, compiled with `make` or `make main`).
This file generates a graph and tries to partition it with different strategies.
For each partition method a `.dot` file is generated and can be converted in to a pdf/png using the dot tool.
Also, each execution emulation will give an execution trace in `.svg`.

## Help

```bash
$ ./main --help
[HELP] Asked for help.
[HELP] You can only pass the graph generation or a filename.
[HELP] $ ./main [generation method]
[HELP] Where generation method is among: tree, deptree, 2dgrid, doubletree, filename
```

## Example

Each tested method will create `.dot`/`.svg` files and inform about an emulated execution time.

```bash
$ ./main tree
[INFO] use tree
nbThreads : 2 / partMinSize : 4 / partMaxSize : 4
Number of nodes : 36
Degree of parallelism one the original graph : 8  4.88889
Generate pdf of the original graph with: dot -Tpdf /tmp/agraph-original.dot -o /tmp/agraph-original.pdf
Generate pdf of the graph with: dot -Tpdf /tmp/agraph-rand.dot -o /tmp/agraph-rand.pdf
Degree of parallelism after random partitioning : 2  1.11111
Number of partitions : 9 -- avg part size : 4
[EXECUTION] Ready tasks at starting time => 1
[EXECUTION] Total duration => 32
...
```

Will create
```bash
$ ls /tmp/*.svg /tmp/*.dot
/tmp/agraph-backtrack.dot   /tmp/agraph-rand.dot                 /tmp/depgraph-diamond.dot
/tmp/agraph-diamond.dot     /tmp/dep-graph-2trace-backtrack.svg  /tmp/depgraph.dot
/tmp/agraph.dot             /tmp/dep-graph-2trace-greedy.svg     /tmp/depgraph-greedy.dot
/tmp/agraph-greedy.dot      /tmp/dep-graph-2trace-rand.svg       /tmp/depgraph-horizontal.dot
/tmp/agraph-horizontal.dot  /tmp/dep-graph-2trace.svg
/tmp/agraph-original.dot    /tmp/depgraph-backtrack.dot
```

Converting the `.dot` can be done one by one using the commands given in the output of the execution, or for all files at once with (be aware that if the graphes are huge this might take forever)
```bash
for fl in /tmp/*.dot ; do dot -Tpdf $fl -o $fl.pdf ; done
```

** Note that the files called `depgraph-...` represent the new clustered graphes. **

Currently, the size of the clustering and the number of threads in the emulated executions must be changed in the code:
```cpp
    const int nbThreads = 2;
    const int partMaxSize = 4;
```

## Chukrut example

The main binary can execute `.txt` and `.dot` files (it will simply looks at the dependency in the `.dot`).
Therefore a Chukrut data could be run with:
```bash
./main ../data/graph-140721493545840.dot
```

## Install Metis

```bash
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
tar -xvf metis-5.1.0.tar.gz
rm metis-5.1.0.tar.gz
cd metis-5.1.0
make config prefix=$(pwd)/install
make
make install
export METIS_DIR=$(pwd)/install

# Then
cmake .. -DMETIS_DIR=$METIS_DIR
```

## Install Scotch

```bash
wget https://gforge.inria.fr/frs/download.php/file/37622/scotch_6.0.6.tar.gz
tar -xvf scotch_6.0.6.tar.gz
rm scotch_6.0.6.tar.gz
cd scotch_6.0.6/src
cp Make.inc/Makefile.inc.x86-64_pc_linux2 ./Makefile.inc
make
make prefix=$(pwd)/install install
export SCOTCH_DIR=$(pwd)/install

# Then
cmake .. -DMETIS_DIR=$METIS_DIR -DSCOTCH_DIR=$SCOTCH_DIR -DUSE_SCOTCH_AS_METIS=ON
```


