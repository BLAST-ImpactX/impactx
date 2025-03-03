# Differentiable ImpactX


## Prepare Developer Environment

### Once
```bash
conda env create -n clad -c conda-forge clad clangdev clang
```

### Always

```bash
conda activate clad

export CC="clang"
export CXX="clang++"
export CXXFLAGS="-I ${CONDA_PREFIX}/include -fplugin=${CONDA_PREFIX}/lib/clad.so"
```


## Compile

With the active developer env above:
```bash
cmake --fresh -S . -B build -DImpactX_UNITY_BUILD=ON -DImpactX_MPI=OFF -DImpactX_COMPUTE=NOACC -DImpactX_OPENPMD=OFF
cmake --build build -j 6
```


## Test/Run

With the active developer env above:
```bash
./build/bin/impactx examples/fodo/input_fodo.in
```
