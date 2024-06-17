<div align="center">

![](https://raw.githubusercontent.com/ViennaTools/ViennaCore/main/assets/logo.png)

<h1>ViennaCS</h1>

<!-- [![🧪 Tests](https://github.com/ViennaTools/ViennaCS/actions/workflows/build.yml/badge.svg)](https://github.com/ViennaTools/ViennaCS/actions/workflows/build.yml) -->

</div>


ViennaCS is a header-only C++ cell set library, which adds the possibility of using volumetric representations on top of existing level-set functionalities for surfaces. Combined with ray tracing techniques, this enables the simulation of particle scattering and ion implantation.

> [!NOTE]
> ViennaCS is under heavy development and improved daily. If you do have suggestions or find bugs, please let us know!

## Releases
Releases are tagged on the main branch and available in the [releases section](https://github.com/ViennaTools/ViennaCS/releases).

## Building

### Supported Operating Systems

* Linux (g++ / clang)

* macOS (XCode)

* Windows (Visual Studio)

### System Requirements

* C++17 Compiler with OpenMP support

### Dependencies (installed automatically)

* [ViennaCore](https://github.com/ViennaTools/viennacore)

* [ViennaLS](https://github.com/ViennaTools/ViennaLS)

* [ViennaRay](https://github.com/ViennaTools/ViennaRay)
  
* [pybind11](https://github.com/pybind/pybind11) (only for building Python libs)

The CMake configuration automatically checks if the dependencies are installed.
If the dependencies are not found on the system, they will be built from source. To use local installations of the dependencies, the `VIENNACS_LOOKUP_DIRS` variable can be set to the installation path of the dependencies.

## Installing

Since this is a header only project, it does not require any installation. However, we recommend the following procedure in order to set up all dependencies correctly:

```bash
git clone https://github.com/ViennaTools/ViennaCS.git
cd ViennaCS

cmake -B build && cmake --build build
cmake --install build --prefix "/path/to/your/custom/install/"
```

This will install the necessary headers and CMake files to the specified path. If `--prefix` is not specified, it will be installed to the standard path for your system, usually `/usr/local/` . 

## Building the Python package

The Python package can be built and installed using the `pip` command:

```bash
git clone https://github.com/ViennaTools/ViennaCS.git
cd ViennaCS

pip install .
```

> Some functionalities of the ViennaCS Python module only work in combination with the ViennaLS Python module. It is therefore recommended to additionally install the ViennaLS Python module on your system. Instructions to do so can be found in the [ViennaLS Git Repository](https://github.com/ViennaTools/viennals).

## Using the Python package

The 2D version of the library can be imported as follows:
```python
import viennacs2d as vcs
```

In order to switch to three dimensions, only the import needs to be changed:

```python
import viennacs3d as vcs
```

## Integration in CMake projects

We recommend using [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake) to consume this library.

* Installation with CPM
  ```cmake
  CPMAddPackage("gh:viennatools/viennacs@3.0.0")
  ```

* With a local installation
    > In case you have ViennaCS installed in a custom directory, make sure to properly specify the [`CMAKE_PREFIX_PATH`](https://cmake.org/cmake/help/latest/envvar/CMAKE_PREFIX_PATH.html#envvar:CMAKE_PREFIX_PATH).

    ```cmake
    list(APPEND CMAKE_PREFIX_PATH "/your/local/installation")

    find_package(ViennaCS)
    target_link_libraries(${PROJECT_NAME} PUBLIC ViennaTools::ViennaCS)
    ```

## Examples

### Building

The examples can be built using CMake:
> __Important__: Make sure all dependencies are installed and have been built previously

```bash
git clone https://github.com/ViennaTools/ViennaCS.git
cd ViennaCS

cmake -B build -DVIENNACS_BUILD_EXAMPLES=ON
cmake --build build
```

## Tests

ViennaCS uses CTest to run its tests.
In order to check whether ViennaCS runs without issues on your system, you can run:

```bash
git clone https://github.com/ViennaTools/ViennaCS.git
cd ViennaCS

cmake -B build -DVIENNACS_BUILD_TESTS=ON
cmake --build build
ctest -E "Benchmark|Performance" --test-dir build
```

## Contributing

If you want to contribute to ViennaCS, make sure to follow the [LLVM Coding guidelines](https://llvm.org/docs/CodingStandards.html).

Make sure to format all files before creating a pull request:
```bash
cmake -B build
cmake --build build --target format
```

## Authors

Current contributors: Tobias Reiter, Felix Strasser

Contact us via: viennatools@iue.tuwien.ac.at

ViennaCS was developed under the aegis of the 'Institute for Microelectronics' at the 'TU Wien'.
http://www.iue.tuwien.ac.at/

## License

See file [LICENSE](LICENSE) in the base directory.
