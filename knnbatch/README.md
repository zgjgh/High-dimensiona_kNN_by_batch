# batchkNN — Build & Run (Linux)

This guide explains how to build and run the project on Linux using CMake.

---

## 1) Build (out-of-source)

From the project root:

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
```

The executable will be created in the `build/` directory:
```
./batch_knn
```

---

## 2) Datasets

The repository includes two datasets under `datasets/`:
- `Normalized_WT.dat`
- `enron.txt`

These are sufficient to test the system. For larger datasets, you can refer to:  
https://github.com/DBAIWangGroup/nns_benchmark

If downloaded, place additional datasets under the `datasets/` directory and reference them with relative paths from the `build/` directory.

---

## 3) Usage

The program requires two parameters:

```
./batch_knn <dataset_path> <mode>
```

- `<dataset_path>`: Path to the dataset file. When running from inside `build/`, use relative paths such as:
  ```
  ../datasets/Normalized_WT.dat
  ../datasets/enron.txt
  ```
- `<mode>`: Update mode. Must be either `batch` or `single`.

### Examples

```bash
# Run on Normalized_WT.dat in single mode
./batch_knn ../datasets/Normalized_WT.dat single

# Run on enron.txt in batch mode
./batch_knn ../datasets/enron.txt batch
```

---

## 4) Requirements

- Linux
- CMake ≥ 3.12
- A C++14-compatible compiler (GCC or Clang)

On Ubuntu/Debian-based systems, if missing:
```bash
sudo apt-get update
sudo apt-get install -y cmake build-essential
```

---

## 5) Notes

- Always build inside a separate `build/` directory (out-of-source).
- Paths to datasets should be given relative to the `build/` directory.
- Use forward slashes (`/`) in paths.
- To clean and rebuild, remove the `build/` directory and repeat the steps above.
