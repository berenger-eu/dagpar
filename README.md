This Software is a test project to study the acyclic partitioning of DAG (of tasks).

# Compilation

Simply go in the build dir and use cmake as usual:
```
mkdir build
cd build
cmake ..
make
```

If everything is fine, the tests should be correct:
```
make test
```

# Adding a test file

Any `.cpp` file added in `tests` will be automatically compiled after using `cmake ..` and `make` in the build dir.

# Where to make modifications

The `partition` method of the `Graph` class should do the partitioning.
Anyone interested to improve the method should modify it or duplicate it and look at the obtained results.

# Running the example

One can execute the `main` binary (from the `tests/main.cpp` file, compiled with `make` or `make main`).
This file generates a graph and tries to partition it with different strategies.
For each partition method a `.dot` file is generated and can be converted in to a pdf/png using the dot tool.
Also, each execution emulation will give an execution trace in `.svg`.
