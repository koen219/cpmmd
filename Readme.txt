This repo runs on a combination of rust and python using pyo3. At moment of writing pyo3 is not supported for py3.14, I tested with 3.12 and that works.

Make sure the requirements for hoomd are installed. That includes at least:

- pybind11
- cmake
- eigen3
- cereal

For mac I could install them with 

```
    brew install pybind11 cmake eigen@3
```

Hoomd tries to install some of the packages itself and might fails in doing so. So if it fails with a cmake error not installing one of the required packages try installing them yourself and then doing make again. Good luck.
