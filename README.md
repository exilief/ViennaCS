<div align="center">

![](https://raw.githubusercontent.com/ViennaTools/ViennaLS/master/assets/logo.png)

<h1>ViennaCS</h1>

[![🧪 Tests](https://github.com/ViennaTools/ViennaCS/actions/workflows/build.yml/badge.svg)](https://github.com/ViennaTools/ViennaCS/actions/workflows/build.yml)

</div>


ViennaCS is a header-only C++ cell set library, which adds the possibility of using volumetric representations on top of existing level-set functionalities for surfaces. Combined with ray tracing techniques, this enables the simulation of particle scattering and ion implantation.

> [!NOTE]
> ViennaCS is under heavy development and improved daily. If you do have suggestions or find bugs, please let us know!

## Releases
Releases are tagged on the main branch and available in the [releases section](https://github.com/ViennaTools/ViennaCS/releases).

## Building

### Supported Operating Systems

* Windows (Visual Studio)

* Linux (g++ / clang)

* macOS (XCode)

### System Requirements

* C++17 Compiler with OpenMP support

### Dependencies (installed automatically)

* [ViennaLS](https://github.com/ViennaTools/ViennaLS)

* [ViennaRay](https://github.com/ViennaTools/ViennaRay)

## Installing

Since this is a header only project, it does not require any installation. However, we recommend the following procedure in order to set up all dependencies correctly:

```bash
git clone https://github.com/ViennaTools/ViennaCS.git
cd ViennaCS

cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
cmake --build build
cmake --install build
```

This will install the necessary headers and CMake files to the specified path. If `-DCMAKE_INSTALL_PREFIX` is not specified, it will be installed to the standard path for your system, usually `/usr/local/`.

## Installing with dependencies already installed on the system

The CMake configuration automatically checks if the dependencies are installed. If CMake is unable to find them, the dependencies will be built from source.

## Running the Tests

ViennaCS uses CTest to run its tests.
In order to check whether ViennaCS runs without issues on your system, you can run:

```bash
git clone https://github.com/ViennaTools/ViennaCS.git
cd ViennaCS

cmake -B build -DVIENNACS_BUILD_TESTS=ON
cmake --build build
ctest -E "Benchmark|Performance" --test-dir build
```

## Building examples

The examples can be built using CMake:

```bash
cmake -B build -DVIENNACS_BUILD_EXAMPLES=ON
cmake --build build
```

## Integration in CMake projects

We recommend using [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake) to consume this library.

* Installation with CPM
  ```cmake
  CPMAddPackage("gh:viennatools/viennacs@1.0.0")
  ```

* With a local installation
    > In case you have ViennaCS installed in a custom directory, make sure to properly specify the `CMAKE_MODULE_PATH` or `PATHS` in your `find_package` call.

    ```cmake
    list(APPEND CMAKE_PREFIX_PATH "/your/local/installation")

    find_package(ViennaCS)
    target_link_libraries(${PROJECT_NAME} PUBLIC ViennaTools::ViennaCS)
    ```

## Contributing

If you want to contribute to ViennaCS, make sure to follow the [LLVM Coding guidelines](https://llvm.org/docs/CodingStandards.html). Before creating a pull request, make sure ALL files have been formatted by clang-format, which can be done using the `format-project.sh` script in the root directory.

## Authors

Current contributors: Tobias Reiter, Noah Karnel, Julius Piso

Contact us via: viennatools@iue.tuwien.ac.at

ViennaCS was developed under the aegis of the 'Institute for Microelectronics' at the 'TU Wien'.
http://www.iue.tuwien.ac.at/

## License

See file [LICENSE](LICENSE) in the base directory.
