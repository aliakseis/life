# HashLife Game of Life (C++)

A high-performance C++ implementation of **Conway's Game of Life** based on the famous **HashLife** algorithm originally described by **Tomas G. Rokicki**.

This project is an adaptation of Rokicki's Java implementation, rewritten in C++ and extended to support advancing the simulation by an arbitrary number of generations in a single operation. It also includes an optional multithreaded execution mode for improved performance on multi-core systems.

## Features

* HashLife algorithm implementation
* Infinite universe simulation using quadtrees
* Structural sharing through node interning
* Memoized generation caching
* Arbitrary generation jumps
* Optional multithreaded execution
* Extremely large sparse universe support
* Canonical node deduplication using hash tables
* Efficient output generation

## Background

Traditional Game of Life implementations update every cell for every generation, resulting in a computational complexity proportional to the size of the simulated grid.

HashLife takes a completely different approach:

* The universe is represented as a recursive quadtree.
* Identical regions are represented by the same canonical node.
* Previously computed future states are cached.
* Large regions can be advanced by many generations recursively.
* Time complexity becomes dependent on pattern complexity rather than grid size.

This allows simulation of patterns that would be impractical using conventional cellular automata techniques.

## Architecture

### Node

The core of the implementation is the `Node` class.

Each node represents a square region of the universe and contains four child quadrants:

```text
+-----+-----+
| NW  | NE  |
+-----+-----+
| SW  | SE  |
+-----+-----+
```

Each node stores:

* Four child nodes
* Tree level
* Cached future states
* Alive-state information
* Hash-consing metadata

Nodes are immutable after creation, enabling safe sharing throughout the quadtree.

### Interning

To avoid duplication, nodes are canonicalized through a hash table.

When a node is created:

1. A hash is computed from its four children.
2. The hash table is searched.
3. Existing identical nodes are reused.
4. Otherwise a new canonical node is inserted.

This dramatically reduces memory consumption for repetitive structures.

### Quadtrees

The universe is stored as a recursive quadtree.

A level `n` node represents a square region of size:

```text
2^n × 2^n
```

As patterns grow, the universe automatically expands using:

```cpp
expandUniverse()
```

allowing effectively infinite simulation space.

## HashLife Algorithm

The heart of the implementation is:

```cpp
Node::nextGeneration()
```

HashLife recursively computes future generations by:

1. Splitting a region into overlapping subregions.
2. Computing future states for smaller regions.
3. Combining cached results.
4. Recursively advancing larger regions in time.

Computed results are stored inside:

```cpp
Node* result[2];
```

allowing future requests for the same state to be returned instantly.

### Memoization

Each node caches its future evolution.

When a state is requested:

```cpp
if (result[direct] != 0)
    return result[direct];
```

previously computed results are reused immediately.

This is the key optimization responsible for HashLife's extraordinary performance.

## Arbitrary Time Advancement

Unlike many HashLife examples that advance only by powers of two generations, this implementation supports advancing by any number of generations:

```cpp
universe.runSteps(1000);
```

Internally a supplemental bitmask controls partial recursive advancement:

```cpp
nextGeneration(stepSize - numSteps);
```

This allows efficient execution of non-power-of-two generation counts.

## Experimental Multithreaded Mode

When compiled with:

```cpp
#define MULTITHREADED
```

the simulation distributes recursive HashLife computations across worker threads.

### Parallel Execution

The helper class:

```cpp
NextGeneration_MT
```

creates asynchronous work items using:

```cpp
QueueUserWorkItem()
```

Each major sub-quadrant can evolve independently.

### Thread Safety

Thread-safe node interning is implemented using:

```cpp
InterlockedCompareExchange()
```

and

```cpp
std::atomic_bool
```

to prevent duplicate canonical node creation and protect cached results.

## Level-2 Optimization

The smallest active simulation level is handled by:

```cpp
slowSimulation()
```

Instead of recursively evaluating neighborhoods, the implementation uses:

```cpp
cachedOneGen[0x800]
```

which precomputes all possible:

```text
2^11 = 2048
```

neighbor configurations.

This converts the base case into extremely fast table lookups.

## Memory Optimization Tricks

The implementation contains several aggressive low-level optimizations.

### Encoded Leaf Nodes

Small leaf patterns are encoded directly inside pointer values:

```cpp
if (this < (void*)0x800)
```

This avoids allocating millions of tiny nodes.

### Structural Sharing

Identical regions are stored only once.

For highly repetitive patterns this can reduce memory requirements by orders of magnitude.

### Fixed-Size Hash Tables

The implementation uses preallocated hash storage:

```cpp
HASH_SIZE = 128 * 1024
```

in multithreaded mode and:

```cpp
HASH_SIZE = 64 * 1024
```

in single-threaded mode.

## Example Pattern

The sample program initializes a replicated glider-like pattern:

```cpp
SetBit(universe, 1, 0);
SetBit(universe, 2, 0);
SetBit(universe, 0, 1);
SetBit(universe, 1, 1);
SetBit(universe, 1, 2);
```

The pattern is duplicated across multiple distant locations and then evolved:

```cpp
universe.runSteps(1000);
```

## Output

After simulation the resulting universe is exported to:

```text
life.txt
```

containing a 1000×1000 snapshot centered around the origin.

Each cell is written as:

```text
0 = dead
1 = alive
```

## Performance

HashLife performance depends primarily on pattern complexity rather than universe size.

For highly regular or repetitive structures, the algorithm can effectively simulate:

* Millions of generations
* Billions of generations
* Extremely large sparse universes

while reusing previously computed evolution states.

The included benchmark measures average execution time over multiple runs and reports the result both to the console and to the generated output file.

## Original Work

Based on:

**Tomas G. Rokicki**
"An Algorithm for Compressing Space and Time"

Originally published in Dr. Dobb's Journal:

http://www.ddj.com/184406478

