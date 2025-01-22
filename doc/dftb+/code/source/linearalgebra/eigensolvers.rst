
.. _diagonalise:

Eigensolvers
============

At current, four eigensolvers are supplied, this will be increased in
the future.

#. LAPACK generalised eigensolver.

#. LAPACK divide and conquer generalised eigensolver.

#. LAPACK relatively robust eigensolver with wrapper code for the
   generalised case.

Eigensolvers / alternatives to be added are

#. FORTRAN 90 implementation of ewevge

#. Zoltan’s stored Cholesky optimisation (can this be done with current
   data structures?).

#. Dhillon’s new M\ :math:`^3`\ R relatively robust solver.

#. Inverse iteration from eigenvalues in :math:`O(N)` memory

#. divide and conquer.

#. Fermi expansion.

#. .....
