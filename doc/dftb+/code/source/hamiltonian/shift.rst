Potential
---------

This is more straightforward, since it becomes a contribution to the
sparse hamiltonian of the form

.. math::

   \begin{aligned}
     \Delta H^\sigma (m \in l \in a, m^\prime \in l \in a, \mathbf{R} ) &
     = & - ( U_l - J_l ) \rho^\sigma(m \in l \in a, m^\prime \in l \in a,
     \mathbf{R} )\\
     \Delta H^\sigma(m = m^\prime \in l \in a, \mathbf{0} ) &
     = & + \frac{( U_l - J_l )}{2}
   \end{aligned}

Where the second term shifts the trace of the density matrix for those
elements relevant to the potential.

Block populations
-----------------

The expressions above are written for the density matrix, instead we can
use a generalisation of Mulliken
populations (Han PRB 73 045110). The conventional
Mulliken populations given in Eqn. `[eqn:mulliken] <#eqn:mulliken>`__
are the *diagonal* of the generalised block populations:

.. math::

   \begin{aligned}
     q_{\mu\nu}^A &=& \frac{1}{2} \sum\limits_\tau \left( S_{\mu\tau} \rho_{\tau\nu} +
     \rho_{\mu\tau} S_{\tau\nu} \right); \mu,\nu\in A
   \end{aligned}

In the code, each on-site and diatomic block of the sparse matrices are
processed to produce this term. The resulting populations are symmetric
matrices, and in the case of spin, the analysis can be performed on each
channel of :math:`\rho` separately. For *imaginary* coefficients (see
section `8.4 <#sec:Pauli>`__) the skew-symmetric part of these matrices
is needed instead.

In the DFTB+U (and dual spin-orbit) extensions to the hamiltonian, these
block-like terms are used instead of the on-site part of the density
matrix in the energy expressions. To add terms to the hamiltonian using
block forms requires the reverse of the Mulliken process above to
calculate the contribution from the (block) on-site potential
(:math:`V`) contribution:

.. math::

   \begin{aligned}
     H_{\mu\nu} = \frac{1}{2} \left( S_{\mu\tau} V_{\tau \nu} +
     V_{\mu\upsilon} S_{\upsilon \nu} \right); \nu,\tau\in A;
     \mu,\upsilon \in B
   \end{aligned}
