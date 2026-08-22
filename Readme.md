# cpm-md

This repo combines Rust and Python via [PyO3](https://pyo3.rs/). It builds a native Rust extension (`pycpm2`) alongside a Python package, and also builds [HOOMD-blue](https://hoomd-blue.readthedocs.io/) as a dependency.

> **Note:** PyO3 does not yet support Python 3.14. Use **Python 3.12** — this has been tested and works.

## 1. Clone the repository

This repo uses git submodules, so clone with `--recursive`:

```bash
git clone --recursive <URL>
```

If you already cloned without `--recursive`, initialize the submodules afterward:

```bash
git submodule update --init --recursive
```

## 2. Install prerequisites

You'll need:

- **Python 3.12**
- **Rust / rustc**
- **HOOMD build dependencies**, at minimum:
  - `pybind11`
  - `cmake`
  - `eigen3`
  - `cereal`

On macOS, most of these can be installed via Homebrew:

```bash
brew install python@3.12 pybind11 cmake eigen@3
```

> HOOMD's own build process attempts to install some of these packages itself and may fail partway through (often with a `cmake` error about a missing package). If that happens, install the missing dependency yourself (e.g. via `brew`) and re-run `make`.

## 3. Build

From the repository root, start with:

```bash
make
```

### If the top-level `make` fails

Build HOOMD manually first:

```bash
cd lib/hoomd
source ../../.venv/bin/activate   # activate the root .venv before building
make install
```

Once HOOMD builds and installs successfully, go back to the repo root and run `make` again:

```bash
cd ../..
make
```

## Troubleshooting

- **PyO3 build error mentioning Python 3.14**: your `.venv` was created with the wrong Python version. Recreate it with Python 3.12 (see Makefile's `PY` variable) and rebuild.
- **HOOMD `cmake` errors about missing packages**: install the missing dependency manually (see step 2) and re-run `make`.
- **`numpy` `AttributeError` (e.g. `in1d`)**: likely a NumPy 2.x incompatibility with a dependency. Either pin `numpy<2` in `requirements.txt`, or update the offending code to use `numpy.isin` instead of the removed `numpy.in1d`.


Good luck!

