# Space-Optimized Van Emde Boas Tree

Developed as part of the **CS201 – Data Structures** course, this project implements an **optimized and space-efficient Van Emde Boas (vEB) Tree** in **C++**.

The data structure supports the following operations:

- Insert  
- Delete  
- Membership Query  
- Successor  
- Predecessor  
- Minimum  
- Maximum  

All operations run in **O(log log U)** time, where **U** represents the universe size.

---

## Features

- Fully functional **Van Emde Boas Tree implementation**
- **Sparse cluster representation** using hash maps to reduce memory usage
- Multiple optimizations based on the universe size
- Interactive **menu-driven interface** for testing and experimentation

---

## Key Optimizations

### 1. Small Universe Mode (`U ≤ 64`)

For very small universes, the implementation switches to a **64-bit bitmask representation**.

- Operations are implemented using **bit manipulation**
- Uses built-in functions such as `__builtin_ctzll` and `__builtin_clzll`
- Extremely **fast and memory-efficient**

---

### 2. Medium Universe Mode (`64 < U ≤ 256`)

For moderately sized universes, the implementation uses a **compact byte array**.

- Efficient for small datasets
- Successor, predecessor, and min/max operations are handled via **simple linear scans**
- Avoids the overhead of recursive vEB structures

---

### 3. Large Universe Mode (Sparse vEB)

For large universes, a **sparse Van Emde Boas Tree** is used.

Key design decisions:

- Clusters are **created dynamically only when needed**
- A **custom chained hash map** stores cluster pointers
- A **summary tree** keeps track of non-empty clusters

This significantly **reduces memory usage compared to the classical vEB tree**.

---

## Data Structures Used

The implementation combines several techniques:

- Recursive **Van Emde Boas Tree structure**
- Custom **chained hash map** for sparse cluster storage
- **Bit manipulation** for fast operations in small universes

---

## Compilation and Execution

To compile and run the program:

```bash
g++ -std=gnu++17 src/veb_sparse.cpp -o veb
./veb
```