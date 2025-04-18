*****************************
Building and installing DFTB+
*****************************

If you have problems with the build, you can find suggestions for some
frequently occurring scenarios in the `Troubleshooting <#troubleshooting>`_
section at the bottom.


Requirements
============

In order to compile DFTB+, you need the following software components:

* Fortran compiler supporting Fortran 2018 and OpenMP 4.0

* C compiler

* C++ compiler (when built with ELSI/PEXSI or ChIMES support)

* CMake (version 3.16 or newer)

* GNU make

* LAPACK/BLAS libraries (or compatible equivalents)

* Python (version >= 3.2) for the source preprocessor


Fortran compiler
----------------

The following Fortran compilers are known to build DFTB+ correctly:

* GNU >= 12.2

* Intel >= 2021.5

* NAG >= 7.2 (when built without OpenMP support)

Older versions of the compilers above are likely to fail due to missing Fortran
features and/or compiler bugs. Compilers by other vendors may work, but have not
been tested extensively (see also `Tested build environments
<#tested-build-environments>`_ and `Testing DFTB+ <#testing-dftb>`_).


Optional extra dependencies
---------------------------

Additionally there are optional requirements for some DFTB+ features:

* ScaLAPACK (version 2.0 or later) and a Fortran aware MPI framework, if you
  want to build the MPI-parallelised version of the code.

* In addition to ScaLAPACK, for MPI parallel builds it is recommended
  to use the `ELSI <https://wordpress.elsi-interchange.org/>`_ library
  for large scale systems (version 2.6.x of the library, with partial
  support of 2.5.0). If ELSI was compiled with PEXSI included, you
  will also need a C++ compiler.

* The ARPACK-ng library if using the excited state DFTB functionality. For
  MPI-parallel builds, the parallel version of ARPACK-ng (containing also
  PARPACK) is needed.

* The `MAGMA <http://icl.cs.utk.edu/magma/>`_ library for GPU
  accelerated computation (note that within ELSI, the ELPA library
  also supports distributed multiple GPUs if compiled with the correct
  options). The number of the available GPUs used by the MAGMA library
  is controlled at runtime by the `MAGMA_NUM_GPUS` shell variable
  (the usual default is 1).

* The `PLUMED2 <https://github.com/plumed/plumed2>`_ library for
  metadynamics simulations. If you build DFTB+ with MPI, the linked
  PLUMED library must also be MPI-aware (and must have been built with
  the same MPI-framework as DFTB+).


External library requirements
-----------------------------

* **Make sure that all external libraries are compiled with the same kind models
  for the numeric variables** (same integer size and floating point precision)
  as DFTB+. Also, they should preferably have been built with the same compiler
  and with similar compiler flags to DFTB+. (See the Troubleshooting section for
  further information.)

* External libraries in non-standard locations (as is typical on many
  HPC-systems using environment modules) can only be reliable found by CMake if
  their library path occurs in the ``CMAKE_PREFIX_PATH`` environment
  variable. **Make sure that your CMAKE_PREFIX_PATH environment variable
  contains all relevant library paths!**


Requirements for testing DFTB+
------------------------------

In order to execute the code tests and validate them against precalculated
results, you will additionally need:

* Python (version >= 3.2) with NumPy

* The Slater-Koster data used in the tests (see below)


Tested build environments
-------------------------

DFTB+ is regularly built and tested for both serial and MPI environments on the
following architectures:

+---------------+----------------------+-------------+------------------+-----+
| Architecture  | Compiler             | MPI         | Ext. libraries   |Notes|
+===============+======================+=============+==================+=====+
| x86_64 /      | GNU Fortran/C 12.2   | OpenMPI 4.1 | OpenBlas 0.3.21, |     |
| Linux         |                      |             | ScaLAPACK 2.2,   |     |
|               |                      |             | ELSI 2.9         |     |
+---------------+----------------------+-------------+------------------+-----+
| x86_64 /      | GNU Fortran/C 13.2   | OpenMPI 5.0 | OpenBlas 0.3.25, |     |
| Linux         |                      |             | ScaLAPACK 2.2,   |     |
|               |                      |             | ELSI 2.9         |     |
+---------------+----------------------+-------------+------------------+-----+
| x86_64 /      | GNU Fortran/C 13.2   | OpenMPI 5.0 | OpenBlas 0.3.29, |     |
| Linux         |                      |             | ScaLAPACK 2.2,   |     |
|               |                      |             | ELSI 2.11        |     |
+---------------+----------------------+-------------+------------------+-----+
| x86_64 /      | Intel Fortran/C      | IntelMPI    | MKL 2022.0,      |     |
| Linux         | 2021.5               | 2021.14     | ELSI 2.8         |     |
+---------------+----------------------+-------------+------------------+-----+
| x86_64 /      | Intel Fortran/C      | IntelMPI    | MKL 2024.2,      |     |
| Linux         | 2024.2               | 2021.14     | ELSI 2.9         |     |
+---------------+----------------------+-------------+------------------+-----+
| x86_64 /      | Intel Fortran/C      | IntelMPI    | MKL 2025.0,      |     |
| Linux         | 2025.0               | 2021.14     | ELSI 2.11        |     |
+---------------+----------------------+-------------+------------------+-----+
| x86_64 /      | NAG Fortran 7.2      | MPICH 4.2   | OpenBlas 0.3.26  | [1] |
| Linux         | GNU C 13.2           |             | ScaLAPACK 2.2    |     |
+---------------+----------------------+-------------+------------------+-----+

Notes:

[1] Only Debug build is tested regulary with OpenMP turned off and without ELSI.


Obtaining the source
====================

The source code of the last stable release can be downloaded from the `DFTB+
homepage <https://www.dftbplus.org/download/stable.html>`_.

Alternatively you can clone the `public git repository
<https://github.com/dftbplus/dftbplus>`_. The tagged revisions correspond to
stable releases, while the default branch contains the latest development
version. ::

  git clone https://github.com/dftbplus/dftbplus.git
  cd dftbplus

The project uses git-submodules for some external dependencies, which will be
automatically retrieved during configuration.


Optional extra components
-------------------------

Some optional software components are not distributed with the DFTB+ source code
and are also not retrieved automatically. If these are required, you can
download these components by using the `get_opt_externals` utility, e.g.::

  ./utils/get_opt_externals

This will download all license compatible optional external components. These
include the Slater-Koster (slako) data for testing the compiled code.

For more information see the detailed help for this tool by issuing
``./utils/get_opt_externals -h``.


Slater-Koster file locations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The DFTB+ code checks the shell variable `DFTBPLUS_PARAM_DIR` when
setting the path to check the Prefix keyword for finding data. If
unset, it assumes the local directory as the starting path.


Building
========

**Important note:** CMake caches its variables in the `CMakeCache.txt` file in
the build folder (e.g. ``_build/CMakeCache.txt``). Make sure to delete this file
before re-running CMake if you have changed variables in `config.cmake` or in
the toolchain files in the `sys/` folder. (Deleting the `CMakeCache.txt` file is
not necessary if you change a variable via the ``-D`` command line option.)

In order to build DFTB+ carry out the following steps:

* Inspect the `config.cmake` file and customise the global build parameters. (If
  you are unsure, leave the defaults as they are.)

* Invoke CMake to configure the build. Specify the installation destination
  (e.g. ``$HOME/opt/dftb+``) and pass an arbitrary folder (e.g. ``_build``) for
  the build and the directory containing the source files (e.g. ``.``) as
  arguments to CMake. Additionally define your Fortran and C compilers as
  environment variables, e.g. (in a BASH compatible shell)::

    FC=gfortran CC=gcc cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt/dftb+ -B _build .

  Based on the detected compilers, the build system will read further settings
  from a corresponding toolchain file in the `sys/` folder. Either from a
  compiler specific one (e.g. `gnu.cmake`, `intel.cmake`, etc.) or the generic
  one (`generic.cmake`) if the detected compiler combination does not correspond
  to any of the specific settings. The selected toolchain is indicated in the
  CMake output. (The toolchain file selection can be manually overridden by
  setting the ``TOOLCHAIN`` CMake variable.)

  You may adjust any CMake variable defined in `config.make` or in the
  toolchain files by either modifying the files directly or by setting
  (overriding) the variable via the ``-D`` command line option. For example, in
  order to use the MKL-library with the GNU-compiler, you would have to override
  the ``LAPACK_LIBRARY`` variable with the CMake command line argument ``-D``::

    -DLAPACK_LIBRARY="mkl_gf_lp64;mkl_gnu_thread;mkl_core"

  When needed, you can specify the complete path to a library or pass linker
  options as defined variables, e.g.::

    -DLAPACK_LIBRARY="/opt/openblas/libopenblas.a"
    -DLAPACK_LIBRARY="-Wl,--start-group -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group"

  By default CMake searches for the external libraries in the paths specified in
  the ``CMAKE_PREFIX_PATH`` environment variable. **Make sure that your
  CMAKE_PREFIX_PATH environment variable is set up correctly and contains
  all the relevant paths** when configuring the project, e.g. ::

    CMAKE_PREFIX_PATH=/opt/elsi:/opt/custom-openblas cmake [...] -B _build .

  Some of the external library finders also offer special ``_LIBRARY_DIR`` CMake
  variables for setting search paths, e.g. ::

    -DLAPACK_LIBRARY_DIR=/opt/custom-openblas

  Setting those variables is normally not necessary, provided the right search
  path is already present in the ``CMAKE_PREFIX_PATH`` environment variable.


* If the configuration was successful, start the build by ::

    cmake --build _build -- -j

  This will compile the code using several threads and showing only the most
  relevant information.

  If, for debugging purposes, you wish to see the exact compiling commands, you
  should execute a serial build with verbosity turned on instead::

    cmake --build _build -- VERBOSE=1

* Note: The code can be compiled with distributed memory parallelism (MPI), but
  for smaller shared memory machines, you may find that the performance is
  better when using OpenMP parallelism only and an optimised thread aware BLAS
  library is used.


CMake options for additional functionality
------------------------------------------

A subset of the cmake comand line options are listed below (for full details,
check the `config.cmake` file in the top repository directory).


Enabling optional code functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+-----------------+---------+------------------------------------------------+
|CMake option     |Default  |Notes                                           |
+=================+=========+================================================+
|-DWITH_ARPACK    |N        |Required for some types of excited state        |
|                 |         |calculations (the Stratmann solver does not     |
|                 |         |require this library).                          |
+-----------------+---------+------------------------------------------------+
|-DWITH_CHIMES    |N        |Required for the CHiMES repulsive potentials to |
|                 |         |be enabled.                                     |
+-----------------+---------+------------------------------------------------+
|-DWITH_MBD       |N        |Required for many-body dispersion to be enabled.|
+-----------------+---------+------------------------------------------------+
|-DWITH_PLUMED    |N        |Required for the PLUMED meta-dynamics library to|
|                 |         |be available for MD calculations.               |
+-----------------+---------+------------------------------------------------+
|-DWITH_POISSON   |N        |Required for Poisson solver (also enabled for   |
|                 |         |transport builds).                              |
+-----------------+---------+------------------------------------------------+
|-DWITH_SDFTD3    |N        |Required for 'simple' DFT-D3 dispersion to be   |
|                 |         |enabled.                                        |
+-----------------+---------+------------------------------------------------+
|-DWITH_SOCKETS   |N        |Required for i-PI socket interfaces to be       |
|                 |         |available.                                      |
+-----------------+---------+------------------------------------------------+
|-DWITH_TBLITE    |N        |Required for xTB hamiltonians to be enabled.    |
+-----------------+---------+------------------------------------------------+
|-DWITH_TRANSPORT |N        |Required for open-boundary calculations to be   |
|                 |         |enabled.                                        |
+-----------------+---------+------------------------------------------------+


Parallelism options
^^^^^^^^^^^^^^^^^^^

+-----------------+--------+--------------------------------------------------+
|CMake option     | Default| Notes                                            |
+=================+========+==================================================+
|-DWITH_OMP       | Y      |OpenMP parallelism enabled in the build.          |
+-----------------+--------+--------------------------------------------------+
|-DWITH_MPI       | N      |MPI parallelism enabled in the build.             |
+-----------------+--------+--------------------------------------------------+
|-DWITH_ELSI      | N      |Requires that MPI is also enabled.                |
+-----------------+--------+--------------------------------------------------+
|-DWITH_GPU       | N      |Depending on whether MPI+ELSI is enabled, this    |
|                 |        |will use GPU accelerated ELPA (if provided) or the|
|                 |        |MAGMA library (if provided) for GPU acceleration. |
+-----------------+--------+--------------------------------------------------+


Testing DFTB+
=============

* After successful compilation, change to the build folder and execute the code
  tests::

    pushd _build
    ctest
    popd

  You can also run the tests in parallel in order to speed this up.  If you use
  parallel testing, ensure that the number of OpenMP threads is reduced
  accordingly. As an example, assuming your workstation has 4 cores and you have
  set up the ``TEST_OMP_THREADS`` variable to ``2`` (in `config.cmake`), issue
  ::

    ctest -j2

  for an OpenMP compiled binary running two tests simultaneously, each using 2
  cores.

  If you want to test the MPI enabled binary with more than one MPI-process, you
  should set the ``TEST_MPI_PROCS`` variable accordingly.

  Testing with hybrid (MPI/OpenMP) parallelism can be specified by setting both,
  the ``TEST_MPI_PROCS`` and ``TEST_OMP_THREADS`` variables, e.g::

    set(TEST_MPI_PROCS "2" CACHE STRING "Nr. of processes used for testing")
    set(TEST_OMP_THREADS "2" CACHE STRING "Nr. of OMP-threads used for testing")

  Note that efficient production use of the code in this mode may require
  process affinity (settings will depend on your specific MPI implementation).

  The ``TEST_MPI_PROCS`` and ``TEST_OMP_THREADS`` cache variables can be updated
  or changed also after the compilation by invoking CMake with the appropriate
  ``-D`` options, e.g.::

    cmake -B _build -DTEST_MPI_PROCS=2 -DTEST_OMP_THREADS=2 .
    pushd _build; ctest; popd


Testing related CMake options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some of the command line options for CMake to modify testing behaviour before
building are below.

+------------------------+----------------------------------------------------+
|CMake option            |Notes                                               |
+========================+====================================================+
|-DTEST_MPI_PROCS        |Number of MPI processes used if testing an MPI build|
|                        |with ctest.                                         |
+------------------------+----------------------------------------------------+
|-DTEST_OMP_THREADS      |Number of openMP processes used if testing an openMP|
|                        |build with ctest.                                   |
+------------------------+----------------------------------------------------+
|-DTEST_RUNNER_TEMPLATE  |Modifies the DFTB+ invocation inside the test system|
|                        |(for example if a specific MPI loader is required to|
|                        |run calculations).                                  |
+------------------------+----------------------------------------------------+


Installing DFTB+
================

* The compiled executables, libraries, module files etc. can be copied into an
  installation directory by ::

    cmake --install _build

  where the destination directory can be configured by the variable
  ``CMAKE_INSTALL_PREFIX`` (in the `config.cmake` file). The default location is
  the `_install` subdirectory within the build directory.


Using DFTB+ as a library
========================

DFTB+ can be also be used as a library and linked into other simulation software
packages. In order to compile the library with its public API, make sure to set
the ``WITH_API`` option to ``TRUE`` in the CMake config file
`config.cmake`. When you install the program, it will also install the DFTB+
library, the C-include file and the Fortran module files, which are necessary
for linking DFTB+ with C and Fortran programs.


DFTB+ CMake library options
---------------------------

See `config.cmake` for details, but the main configuration options for building
are listed below.

+-------------------------+---------+-----------------------------------------+
|CMake option             |Default  |Notes                                    |
+=========================+=========+=========================================+
|-DBUILD_SHARED_LIBS      |N        |Build libdftbplus and other components as|
|                         |         |shared libraries.                        |
+-------------------------+---------+-----------------------------------------+
|-DENABLE_DYNAMIC_LOADING |N        |Use shared libraries externally (where   |
|                         |         |possible) for libdftbplus.               |
+-------------------------+---------+-----------------------------------------+
|-DINSTANCE_SAFE_BUILD    |N        |Compile libdftbplus as an instance safe  |
|                         |         |library (the build stops, if any         |
|                         |         |non-instance-safe components have been   |
|                         |         |selected)                                |
+-------------------------+---------+-----------------------------------------+
|-DWITH_API               |N        |Build the API bindings to use libdftbplus|
|                         |         |externally.                              |
+-------------------------+---------+-----------------------------------------+
|-DWITH_PYTHON            |N        |Build the Python3 bindings for           |
|                         |         |libdftbplus. Note that this should also  |
|                         |         |be built for shared libraries.           |
+-------------------------+---------+-----------------------------------------+


Linking the library in CMake based builds
-----------------------------------------

This is the preferred way of invoking the DFTB+ library into your project.  In
CMake based projects you can directly use the CMake export file of DFTB+, which
is installed in the `lib/cmake/dftbplus/` folder in the installation folder. It
exports the target ``DftbPlus::DftbPlus`` which you can use to obtain all
necessary compiler, include and linking options. Your projects `CMakeLists.txt`,
should like something like below::

  project(DftbPlusTest LANGUAGES Fortran C)
  find_package(DftbPlus REQUIRED)
  add_executable(testprogram testprogram.f90)
  target_link(testprogram DftbPlus::DftbPlus)

Note, that this will link all libraries in the correct order, which were
compiled during the DFTB+ build (e.g. libs-dftd3, libnegf, etc.). It will
additionally contain target dependencies for the external libraries needed to
create standalone applications with DFTB+ (e.g. ``LAPACK::LAPACK``,
``Scalapack::Scalapack``, ``Arpack::Arpack``, ``Plumed::Plumed``,
``Magma::Magma``, etc.). You can either use the CMake find-modules shipped with
the DFTB+ source to find those libraries (and to define the corresponding
targets) or create your own, provided they define the appropriate CMake targets.
The  arpack-ng and ELSI libraries offer CMake export files providing the
``ARPACK::ARPACK`` and ``elsi::elsi`` targets, respectively. Make sure, that
CMake can find the relevant export file if the DFTB+ library was compiled with
ELSI or ARPACK required (e.g., by setting up the environment variable
``CMAKE_PREFIX_PATH`` correctly). Note: you may need to install ELSI (not just
point the prefix path to its build system) to generate this file correctly.


Linking the library in non-CMake based builds
---------------------------------------------

Depending on the choice of external components and whether you want to link
DFTB+ to a C or a Fortran binary, you may need different compilation flags and
linker options. You can look up the necessary compiler flags and linker options
in the `dftbplus.pc` pkg-config file, which is usually installed into the
`lib/pkgconfig` folder in the installation directory. You can either inspect the
file directly, or use the ``pkg-config`` tool::

  export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:DFTBPLUS_INSTALL_FOLDER/lib/pkgconfig
  pkg-config --cflags dftbplus   # compilation flags (e.g. include options)
  pkg-config --libs dftbplus     # library linking options
  pkg-config --static --libs dftbplus   # library linking options for static linking

Note, that the flags and libraries shown are either for linking with Fortran or
with C, depending on the value of the configuration option
``PKGCONFIG_LANGUAGE``.

If you compile DFTB+ with ELSI, PLUMED or MAGMA-support, make sure that
pkg-config can also find the respective pkconfig files for these packages. If
you enable support for these components, their libraries are declared as
dependencies in the DFTB+ pkg-config file. Note, if you compile these libraries
themselves, you may have to follow their install processes to generate suitable
pkg-config files in their specified install location(s).

For external dependencies without pkg-config files (e.g. mbd, negf), the options
for linking those libraries can not be queried via pkg-config, and they must be
added manually.


Generating developer documentation
==================================

Developer documentation can be generated using the FORD source code
documentation generator by issuing ::

  cd doc/dftb+/ford && ford dftbplus-project-file.md

in the main source directory. The documentation will be created in the
`doc/dftb+/ford/doc` folder.


Developer build instructions
============================

You should avoid customizing the build by directly changing variables in the
CMake config files, as your changes may accidentally be checked in into the
repository. Instead, create a customized CMake config file, where you
pre-populate the appropriate cache variables. Then use the `-C` option to load
that file::

  FC=gfortran CC=gcc cmake -C custom.cmake -B _build .

The customized config file is read by CMake before the compiler detection
stage. If your config file contains toolchain dependent options, consider
defining the ``DFTBPPLUS_TOOLCHAIN`` environment variable and query it in your
config file.


Relevant CMake options
----------------------

+-------------------+--------------+------------------------------------------+
|CMake option       |Default       |Notes                                     |
+===================+==============+==========================================+
|-DLCOV_REPORT      |N             |Generate coverage reports if using        |
|                   |              |gfortran and lcov is available.           |
+-------------------+--------------+------------------------------------------+
|-CMAKE_BUILD_TYPE  |RelWithDebInfo|Can generate a binary with extra checking |
|                   |              |or profiling enabled.                     |
+-------------------+--------------+------------------------------------------+
|-DWITH_UNIT_TESTS  |N             |Build additional checking tests for       |
|                   |              |specific code routines.                   |
+-------------------+--------------+------------------------------------------+


Advanced build configuration (e.g. for packagers)
=================================================

Controlling the toolchain file selection
----------------------------------------

You can override the toolchain file, and select a different provided case,
passing the ``-DTOOLCHAIN`` option with the relevant name, e.g.::

  -DTOOLCHAIN=gnu

or by setting the toolchain name in the ``DFTBPLUS_TOOLCHAIN`` environment
variable. If you want to load an external toolchain file instead of one from the
source tree, you can specify the file path with the ``-DTOOLCHAIN_FILE`` option
::

  -DTOOLCHAIN_FILE=/some/path/myintel.cmake

or with the ``DFTBPLUS_TOOLCHAIN_FILE`` environment variable.

Similarly, you can also use an alternative build config file instead of
`config.cmake` in the source tree by specifying it with the
``-DBUILD_CONFIG_FILE`` option or by defining the ``DFTBPLUS_BUILD_CONFIG_FILE``
environment variable.


Preventing the download of external sources
-------------------------------------------

Depending on the value of the ``HYBRID_CONFIG_METHODS`` configuration variable,
some dependencies (e.g. mbd, negf, mpifx, scalapackfx) are automatically
downloaded during the configuration phase and built during the DFTB+ build
process. If you want to ensure that nothing gets downloaded during the build,
pass the variable definition ::

  -DHYBRID_CONFIG_METHODS="Find"

to CMake during the configuration. In this case, CMake will only try to find
those dependencies on the system (by searching in the standard system paths and
in the locations defined in the environment variable ``CMAKE_PREFIX_PATH``) and
stop if some components were not found.


Troubleshooting
===============

* **CMake finds the wrong compiler**

  CMake should be guided with the help of the environment variables ``FC``,
  ``CC`` (and eventually ``CXX``) to make sure it uses the right compilers,
  e.g. ::

    FC=gfortran CC=gcc cmake [...]


* **CMake fails to find a library / finds the wrong version of a library**

  In most cases this is due to a misconfigured ``CMAKE_PREFIX_PATH`` environment
  variable. It is essential, that ``CMAKE_PREFIX_PATH`` contains all paths
  (besides default system paths), which CMake should search when trying to find
  a library. Extend the library path if needed, e.g. ::

    CMAKE_PREFIX_PATH="/opt/somelib:${CMAKE_PREFIX_PATH}" cmake [...]


* **ScaLAPACK detection on Ubuntu 20.4 LTS fails**

  The OpenMPI version of ScaLAPACK on Ubuntu 20.4 LTS exports an incorrect CMake
  config file (as of October 2020), which refers to an non-existent
  library. Instead, set the library name with the ``SCALAPACK_LIBRARY`` variable
  explicitely, e.g. ::

    cmake -DSCALAPACK_LIBRARY=scalapack-openmpi [...]

  which should fix the problem.


* **My library settings in a "_LIBRARIES" variable are ignored**

  In order to be consistent with the naming scheme suggested by the CMake
  documentation, all library related cache variables have been changed to
  singular nouns, e.g. ::

    cmake -DSCALAPACK_LIBRARY=scalapack-openmpi [...]

  **instead** of the previous ::

    cmake -DSCALAPACK_LIBRARIES=scalapack-openmpi [...]


* **Fortran libraries compiled with the Intel compiler can not be linked**

  In order to enforce compliance with the Fortran 2003 standard (e.g. allowing
  the automatic allocation of arrays in expressions), DFTB+ passes the
  ``-standard-semantics`` option to the Intel compiler. All external modern
  Fortran dependencies (e.g. ELSI) must also be compiled by using the
  ``-standard-semantics`` or the ``-assume realloc_lhs`` option to ensure
  correct linking.
