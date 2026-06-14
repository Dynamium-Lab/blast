# Contributing

Contributions are welcome.

- Follow the coding style enforced by `.clang-format`. CI checks formatting on every
  pull request.
- Build and run the test suite before opening a PR (see the project `README` and
  `tests/CMakeLists.txt`).
- When adding public C++ API, document it with Doxygen comments (`/** @brief ... */`) so
  it renders in the [C++ API reference](../api/cpp/index.md).
- When adding Python bindings, include a docstring in the `.def(...)` call so it renders
  in the [Python API reference](../api/python/index.md).

## Building the documentation locally

```sh
python -m pip install -r docs/requirements.txt
python -m pip install .            # builds the nanobind extension for Python autodoc
doxygen docs/Doxyfile              # generates C++ XML (run from the repo root)
sphinx-build -b html docs docs/_build/html
```

Open `docs/_build/html/index.html` in a browser.
