*********************
Unit testing in DFTB+
*********************

The 24.1 release candidate, as of commit `62ad90f85` now uses the
`Fortuno <https://fortuno.readthedocs.io/en/latest/>`_ framework for
unit testing. The tests are located in a subdirectory of the `test/`
folder and included via cmake::

  test/src/dftbp/unit
  ├── test/src/dftbp/unit/CMakeLists.txt
  ├── test/src/dftbp/unit/common
  ├── test/src/dftbp/unit/dftb
  ├── test/src/dftbp/unit/include
  ├── test/src/dftbp/unit/io
  ├── test/src/dftbp/unit/math
  └── test/src/dftbp/unit/type

