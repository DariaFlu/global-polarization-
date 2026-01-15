# Global Polarization of Λ Hyperons at MPD (Afterburner)

Afterburner utilities for **Λ-hyperon global polarization** studies using **UrQMD**-based Monte-Carlo in **UniGen (UEvent/UParticle)** format.
The project builds a single shared library with ROOT dictionaries and provides a tiny runner (`burn`) to execute user ROOT macros that call the compiled functions.

---

## What’s inside

The repository contains C++ implementations (compiled into `libgp_macros.so`) and headers:

Main entry points (functions are declared in `include/*.hpp` and implemented in `src/*.cpp`):

1. **Λ enhancement / polarization vector**

* `void add_enhanced_lambda(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat)`

  * Produces an output UniGen file with an XYZ polarization vector for particles (Λ have non-zero vector).

2. **Λ decay simulation**

* `void simulate_lambda_decays(TString inputFile, TString outputFile, TString confInFile, Int_t enhanceStat)`

  * Does the same enhancement step + simulates Λ → p + π.
  * Stores polarization information into the proton.

3. **Global polarization analysis**

* `void calc_global_polarization(TString fileIn, Int_t NFiles, TString outFileName, Int_t enhancedFlag = -100)`

  * Measures global polarization vs kinematic bins (pT, y), fits angular distributions, etc.

Additional helpers:

* `read_unigen_root.*`, `getters.*`, `set_lambda_parameterization.*`, etc.

---

## Requirements

### Software

* **CMake ≥ 3.12**
* **C++17 compiler**
* **ROOT** (with components used by the project: Core, RIO, Tree, Hist, Physics, MathCore, GenVector, Gpad)
* **UniGenFormat** (MPDROOT / UniGen) providing:

  * header: `UParticle.h`
  * library: `libUniGenFormat.so`

### Input data

* UrQMD (or similar) data written in **UniGen ROOT** format (UEvent/UParticle).

---

## Build

### 1) Get the code

```bash
git clone https://github.com/DariaFlu/MPD_afterburner.git
cd MPD_afterburner
```

### 2) Prepare environment (recommended)

Load ROOT and UniGenFormat (via MPDROOT). Typical setup:

```bash
source /path/to/root/bin/thisroot.sh
export MPDROOT=/path/to/mpdroot/install
```

If UniGenFormat is not inside MPDROOT, you can provide explicit paths to CMake (see below).

### 3) Configure and compile

```bash
mkdir build && cd build
cmake ..
make
```

If everything is found, you will get:

* `build/libgp_macros.so` (library with dictionaries)
* `build/burn` (macro runner, if `BUILD_BURN=ON`, default)

### If CMake can’t find UniGenFormat

Either set:

```bash
export MPDROOT=/path/to/mpdroot/install
```

or configure explicitly:

```bash
cmake -S . -B build \
  -DUNIGEN_INCLUDE=/path/to/include \
  -DUNIGEN_LIBDIR=/path/to/lib
```

(where `UParticle.h` is inside `UNIGEN_INCLUDE`, and `libUniGenFormat.so` is inside `UNIGEN_LIBDIR`)

---

## How to use after compilation

### Workflow idea

You **do not run `src/*.cpp` directly**. They are compiled into `libgp_macros.so`.
Instead, you write a **ROOT macro** that includes the project headers and calls the compiled functions.

### 1) Create a macro in `macros/`

Create folder and a macro file:

```bash
mkdir -p macros
```

Example `macros/run_sim.C`:

```cpp
#include "simulate_lambda_decays.hpp"

void run_sim() {
  simulate_lambda_decays(
    "/path/to/input.mcini.root",
    "/path/to/output.root",
    "/path/to/qa_out.root",
    2
  );
}
```

### 2) Run your macro with `burn`

From repo root:

```bash
./build/burn macros/run_sim.C
```

Or from inside `build/`:

```bash
cd build
./burn ../macros/run_sim.C
```

---

## Notes

* The project generates ROOT dictionaries via `ROOT_GENERATE_DICTIONARY(...)` using:

  * `dict/GPDictTypes.hh`
  * `dict/GPDict_LinkDef.h`
* The dictionaries are compiled into `libgp_macros.so`, so ROOT can stream STL containers like:

  * `std::vector<UParticle>`
  * `std::vector<TVector3>`
  * `std::vector<ROOT::Math::XYZVector>`

---

## Quick troubleshooting

### “UniGenFormat not found”

* Ensure `MPDROOT` is set and points to an installation that contains:

  * `$MPDROOT/include/UParticle.h`
  * `$MPDROOT/lib/libUniGenFormat.so`

### Macro can’t include headers

* Keep `#include "simulate_lambda_decays.hpp"` style and run via `burn` from repo root (`./build/burn macros/...`) so includes resolve consistently.
* If you run from a different working directory, prefer absolute paths or run from repo root.

---
