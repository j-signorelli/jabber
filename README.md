# jabber
jabber is a planar wave acoustic forcing library and app-suite to support the modeling of freestream disturbances in high-speed wind tunnels for fluid dynamics simulations.

Please see the <a href="https://j-signorelli.github.io/jabber/index.html">full documentation</a> for more details.

## Building + Installing jabber
jabber uses CMake as a build system. For a default build of jabber:
```
mkdir build ; cd build
cmake ..
make -j 4
```
and to install at a specified `-DCMAKE_INSTALL_PREFIX`:
```
make install
```

## Compile-Time Options
Additional options that may be specified using `cmake -D<Option>=<Value> ..`:

| Option                    | Description                        | Default |
|---------------------------|------------------------------------|---------|
| `JABBER_ENABLE_OPENMP`    | Build jabber with OpenMP support   | OFF     |
| `JABBER_BUILD_APP_LIB`    | Build the jabber::app library. Required to build apps. | ON |
| `JABBER_ENABLE_MPI`       | Build (supported) apps with MPI. | OFF     |
| `JABBER_BUILD_TESTS`      | Build unit tests. Run by `./tests/unit-tests` in build directory. | OFF |
| `JABBER_BUILD_DOCS`       | Build documentation w/ Doxygen. In build directory: `make docs`. | OFF|
| `JABBER_BUILD_PSD`        | Build `jabber_psd` app. | ON |
| `JABBER_BUILD_PLOT`       | Build `jabber_plot` app. | ON |
| `JABBER_BUILD_PROFILE`    | Build `jabber_profile` app. | ON |
| `JABBER_BUILD_OUT`        | Build `jabber_out` app. | ON |
| `JABBER_BUILD_VIZ`        | Build `jabber_viz` app (requires MFEM). | OFF |
| `JABBER_BUILD_PARTICIPANT`| Build `jabber_participant` app (requires preCICE). | OFF |

## Using jabber
The jabber library can either be incorporated directly into an existing code, or the provided apps can be used on their own. To use jabber in a project with CMake, ensure that the `CMAKE_PREFIX_PATH` environment variable is set to the install directory of jabber and simply use `find_package`. An example is provided below:
```cmake
...
find_package(jabber 1.0.0 REQUIRED)
add_executable(example main.cpp)
target_link_libraries(example PRIVATE jabber)
...
```
where the header of `main.cpp` may look like:
```cpp
#include <jabber/jabber.hpp>
...
```
