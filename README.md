This Software is a test project to study the acyclic partitioning of DAG (of tasks).

# Compilation

Simply go in the build dir and use cmake as usual:
```
mkdir build
cd build
cmake ..
make
```

# Adding a test file

Any `.cpp` file added in `tests` will be automatically compiled after using `cmake ..` and `make` in the build dir.

# Where to make modifications

The `partition` method of the `Graph` class should do the partitioning.
Anyone interested to improve the method should modify it or duplicate it and look at the obtained results.
