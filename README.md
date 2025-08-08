# Differentiable ImpactX

## Prepare Developer Environment (Ubuntu 24.02 w/ Clang/LLVM 16)

### Once
```bash
# Ref: https://github.com/EnzymeAD/enzyme-dev-docker/blob/main/Dockerfile
sudo apt install llvm-dev clang-16 lld-16 zlib1g-dev libzstd-dev git automake autoconf cmake make lldb-16 ninja-build build-essential libtool llvm-16-dev libclang-16-dev libomp-16-dev libblas-dev libeigen3-dev libboost-dev python3
python3 -m pip install --break-system-package pipx pathlib2
python3 -m pipx install lit cmake

export CC="clang-16"
export CXX="clang++-16"

# Ref: https://enzyme.mit.edu/Installation/
git clone https://github.com/EnzymeAD/Enzyme ~/src/Enzyme
cd ~/src/Enzyme
cmake --fresh -G Ninja -S enzyme -B build -DLLVM_DIR=/usr/lib/llvm-16/lib/cmake/llvm -DLLVM_EXTERNAL_LIT=$(which lit)
cmake --build build -j 6

cd build
ninja check-enzyme

# linker and lld and lto enforcement in CMake is hard...
cd ..
mkdir mylld
cd mylld
ln -s /usr/lib/llvm-16/bin/lld-link lld
ln -s $(which ld.lld-16) ld.lld  # note: not sufficient yet... somehow hard-coded in compiler detection... use docker
                                 # manually linking /usr/bin/ld.lld-16 as /usr/bin/ld.lld works...
export PATH=$PWD:$PATH

# optional: create a venv for Python
rm -rf ~/src/venv-impactx-enzyme
python3 -m venv ~/src/venv-impactx-enzyme
source ~/src/venv-impactx-enzyme/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade scipy numpy packaging setuptools[core] wheel pytest pytest-benchmark matplotlib PyQt6
```


## Compile Options

From the [homepage](https://enzyme.mit.edu/getting_started/UsingEnzyme/):
> Enzyme supports differentiating C/C++ code through ClangEnzyme and LLDEnzyme as compiler and linker plugins, respectively.
> Clang gives our plugin more flexibility in adding and ordering optimization passes than LLD and therefore using ClangEnzyme could result in better performance than LLDEnzyme.
> However, ClangEnzyme can only differentiate one compilation unit at a time and will therefore fail if the function which you try to differentiate calls functions in other compilation units (generally other .c or .cpp files).
> In these cases we recommend the use of LLDEnzyme in combination with LTO.


## Compile (ClangEnzyme, one TU)

### Always

```bash
source ~/src/venv-impactx-enzyme/bin/activate
alias cmake=$HOME/.local/pipx/venvs/cmake/bin/cmake
export CC="clang-16"
export CXX="clang++-16"

# one TU: Clang Plugin
#   Extra Enzyme options, e.g., print https://enzyme.mit.edu/getting_started/UsingEnzyme/#semantic-options
#   optional add for verbose output: -mllvm -enzyme-print
export CXXFLAGS="-fplugin=$HOME/src/Enzyme/build/Enzyme/ClangEnzyme-16.so"
```

With the active developer env above, inside the ImpactX source dir:
```bash
cmake --fresh \
  -S .        \
  -B build    \
  -DImpactX_UNITY_BUILD=ON  \
  -DImpactX_MPI=OFF       \
  -DImpactX_COMPUTE=NOACC \
  -DImpactX_OPENPMD=OFF   \
  -DCMAKE_LINKER_TYPE=LLD \
  -DCMAKE_LINKER=/usr/lib/llvm-16/bin/lld-link
# optional:
#   -DImpactX_PYTHON=ON -DpyAMReX_IPO=OFF -DImpactX_PYTHON_IPO=OFF

cmake --build build -j 6

# optional:
cmake --build build -j 6 --target pip_install
```


## Compile (LLDEnzyme, multiple TUs)

TODO: redo this part using https://github.com/EnzymeAD/enzyme-dev-docker because ld, lld, ld.ldd with the non-system default is too tricky to get right in CMake (i.e. compiler detection).

### Always

```bash
export PATH=$HOME/.local/pipx/venvs/cmake/bin:$PATH
export CC="clang-16"
export CXX="clang++-16"

# many TU: LDD Plugin
#   https://github.com/EnzymeAD/Enzyme/blob/main/enzyme/Enzyme/CMakeLists.txt
#   https://enzyme.mit.edu/getting_started/UsingEnzyme/#semantic-options
export CXXFLAGS="-fuse-ld=/usr/lib/llvm-16/bin/lld-link -flto"  # -mllvm -enzyme-...
export LDFLAGS="-fuse-ld=/usr/lib/llvm-16/bin/lld-link -flto -Wl,-mllvm -Wl,-load=$HOME/src/Enzyme/build/Enzyme/LLDEnzyme-16.so -Wl,--load-pass-plugin=$HOME/src/Enzyme/build/Enzyme/LLDEnzyme-16.so"
```

With the active developer env above, inside the ImpactX source dir:
```bash
cmake --fresh \
  -S .        \
  -B build    \
  -DImpactX_UNITY_BUILD=ON  \
  -DImpactX_MPI=OFF       \
  -DImpactX_COMPUTE=NOACC \
  -DImpactX_OPENPMD=OFF   \
  -DCMAKE_LINKER_TYPE=LLD \
  -DCMAKE_LINKER=/usr/lib/llvm-16/bin/lld-link

# note: -DCMAKE_LINKER_TYPE=LLD appends general lld not necessarily the right version
#        bend ldd to ldd-16 and /usr/bin/ld.ldd to ld.ldd-16
# -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON

cmake --build build -j 6
```


## Test/Run

With the active developer env above:
```bash
./build/bin/impactx examples/fodo/input_fodo.in
```
