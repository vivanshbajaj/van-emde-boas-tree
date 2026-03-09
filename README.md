# Introduction
This repository contains my **CS201 – Data Structures** project:  
an **optimized, space-efficient implementation of the Van Emde Boas (vEB) Tree** in C++.

The implementation supports **Insert, Delete, Member, Successor, Predecessor, Minimum, and Maximum** operations in  
**O(log log U)** time, where `U` is the universe size.

---

## Features
- Fully functional **Van Emde Boas Tree**
- **Sparse representation** using hash maps for clusters
- Multiple optimizations:
  - **Bitmask mode** for small universes (`U ≤ 64`)
  - **Array mode** for medium universes (`64 < U ≤ 256`)
  - Recursive vEB structure for large universes
- Interactive **menu-driven interface** for testing

---

## Key Optimizations Explained

### 1️. Small Mode (`U ≤ 64`)
- Uses a **64-bit bitmask**
- All operations use bit tricks (`__builtin_ctzll`, `__builtin_clzll`)
- Extremely fast and memory-efficient

### 2️. Medium Mode (`64 < U ≤ 256`)
- Uses a **byte array**
- Simple linear scans for min/max/successor/predecessor

### 3️. Large Universe (Sparse vEB)
- Clusters are stored **only when needed**
- Uses a **custom hash map** to reduce memory usage
- Summary tree tracks non-empty clusters

---

## Data Structures Used

- Recursive **Van Emde Boas Tree**
- Custom **chained hash map** for sparse clusters
- Bit manipulation for fast operations

---

## How to Compile and Run

```bash
g++ -std=gnu++17 src/veb_sparse.cpp -o veb
./veb
