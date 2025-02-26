.. _largeMat:

Large matrix data structure
===========================

The “large” matrices in the calculation are stored in a packed format,
these include the overlap, the non-SCC H\ :math:`^0` hamiltonian, the
SCC additions to the hamiltonian, the density matrix and energy-weighted
density matrix, etc. Hence these are all the matrices which take
:math:`O(N^2)` memory if stored in a square form. In the packed form,
for large calculations, they all take :math:`O(N)` memory in packed
format. It is only when these large matrices are unpacked, for example
for diagonalisation, that there is any :math:`O(N^2)` memory requirement
(this is in a section of code that can be easily replaced for example
with :math:`O(N)` methods or banded matrices to reduce total memory).

For these large matrices, only the non-zero elements of the lower
triangle matrices are stored, based on a cut-off criteria from the
Slater-Koster overlap data. Routines are supplied to convert between
this packed format and the conventional square matrices (see
``unpackHS`` and ``packHS`` interfaces in the code and
section `3.1 <#periodic>`__).

**Important :** Only the elements for which the overlap matrix would be
non-zero are stored – i.e., not the full density matrix or
energy-weighted density matrix. This is all that is required for total
energy and forces, but is not the *full matrix* – hence e.g. commutators
should not be constructed with the packed data, since full matrix-matrix
products require data which is not stored. To calculate the full density
matrix is a similar cost to diagonalising the matrices, being
:math:`O(N^3)`, with a smaller pre-factor than the eigensolver. But
calculating only the elements in the sparse matrix scales as
:math:`O(N^2)`.

| The packed format is as follows :
| For each atom ii starting at atom 1, there are nNeighbor(ii)
  neighbouring atoms. Each neighbour, jj, of atom ii is listed in
  the array ``iNeighbor(:,ii)``, where the 0 to nNeighbor(ii) elements
  of the first index hold the atom numbers. The form of the packed large
  array then runs as the following pseudocode :

.. container:: tt

   ::

      iCount = 1
      do ii = 1, nAtom
        do jj = 0, nNeighbor(ii)
          do kk = 1, ( mAngAtom(ii) + 1 )**2
            do ll = 1, ( mAngAtom( iNeighbor(jj,ii) ) + 1 )**2
              sparseMatrix(iCount) = someFunction
              iCount = iCount + 1
            end do
          endo do
        end do
      end do

Note:

#. The angular momentum for each atom listed in ``mAngAtom(:)`` starts
   with :math:`s = 0`.

#. In practice iCount is replaced by an index array (see code for
   examples) iPair that gives the offset in the large arrays where each
   diatomic/atomic block starts.

#. With periodic boundary conditions, the neighbour list contains extra
   atoms which are periodic images, hence ``mAngAtom(ii)`` should be
   indexed by ``mAngAtom(Img2CentCell(ii))`` instead, whenever there is
   a chance of a periodic case and ``ii`` :math:`>` ``nAtom``

#. The periodic image atoms listed as neighbours also obey the
   requirement that ``Img2CentCell(jj)`` :math:`\geqslant` ii.


To give two examples, one a molecule and one a periodic system, see
figures `2 <#molecule>`__ and `3 <#periodic>`__. In the case of the
molecule the index array ``Img2CentCell(:)`` has the following elements
:math:`(1,2,3,4,5)`\  , and in the periodic case the “bonding” is
:math:`(1,2,3,1,2,3,2,3)`. There are several points to note, in
particular about the periodic example:

#. These are *REAL* matrices. Only the square unpacked forms ever have
   complex elements (ignoring spin-orbit coupling).

#. The green numbered atoms are periodic images of the central cell, in
   surrounding cells at lattice vectors :math:`\mathbf{+R1}` and
   :math:`\mathbf{-R1}`.

#. There is no image of atom 1 in the cell at :math:`\mathbf{+R1}` as it
   would not interact with any atoms in the central cell. This is shown
   as an example that not all atoms in every periodic image cell bond to
   the central cell, only those where the overlap with their basis
   functions is non-zero.

#. The size of ``Img2CentCell(:)`` does not represent the number of
   atoms + image atoms, it may be larger due to optimisations in dynamic
   memory allocation. Use the variable ``nAllAtom`` if you want the
   actual number.

#. There is symmetry in the problem, if you need the upper triangle of a
   cell at :math:`\mathbf{+R_x}`, it is equal to the transpose of the
   lower triangle in the cell :math:`\mathbf{-R_x}`. This simplifies for
   the :math:`\mathbf{R_0}` cell (or molecules) to just the transpose
   inside that square matrix.


.. _neighbors:

Neighbour lists
===============

There are several types of neighbour map maintained internally in the
code.

For a basic non-SCC calculation, the code internally tracks all atoms
that are closer together than the cut-off given in the Slater-Koster
table. This data is stored in two arrays ``nNeighbor(:)`` and
``iNeighbor(:,:)``. nNeighbor holds the number of neighbouring atoms
around each of the real atoms in the calculation (there is no need to
have a neighbour map for interactions between periodic images of atoms).
``iNeighbor(:,ii)`` holds the ``nNeighbor(ii)`` atoms which surround
atom ii, where the 0\ :math:`^\mathrm{th}` neighbour is atom ii itself.

Additional neighbour maps are needed for SCC calculations.

For non-periodic calculations a neighbour list of atoms where the
short-range part of :math:`\gamma` is large enough to cause an effect is
maintained (this is the longest range neighbour list in the
calculation). ``iNeighbor(:,ii)`` also contains the atom numbers of
these “neighbours” of atom ii. The short-range part of :math:`\gamma` is
stored in a sparse format similar to the overlap (see ``shortGamma`` in
the code).

For periodic calculations, in addition to the neighbour list for the
short-range part of :math:`\gamma`, there are further lists for the
real-space part of the Ewald summation and the reciprocal
:math:`G`-vectors to include. See the code for details.

Spin
----

For spin-polarised calculations, since the hamiltonian is pairwise and
cut-off by the overlap matrix, the non-SCC neighbour information is all
that is needed to work out which interactions to include. The on-site
potential is calculated without needing any additional neighbour maps.

.. _periodic:

periodic
--------

Routines to convert between extended sparse :math:`n \times m`
matrices in real space and :math:`n \times n` matrices at arbitrary
:math:`\mathbf{k}` are available in `src/dftbp/dftb/periodic.F90`.

The aim of the code design is to build large matrices in real space
where there is a natural cut off as :math:`S \rightarrow 0`. Arrays of
all neighbouring atoms are supplied (iNeighbor and nNeighbor) to allow
construction of these arrays. Only the lower triangle is
constructed. See section `[periodcoords] <#periodcoords>`__ for a note
about coordinates. The :math:`k`-dependant folding is done using
appropriate phases (for molecules and at :math:`\mathbf{k} = (0,0,0)`
the phase is always :math:`1+0i`):

.. math:: N(\mathbf{k}) = \sum_k M \exp(i \mathbf{k}\cdot\mathbf{R})

\ and

.. math:: M(\mathbf{R}) = \sum_R N \exp(i \mathbf{k}\cdot\mathbf{R})

using two types of routines (pack and unpack). Since no symmetry
routines are implemented yet, all k-points unrelated by inversion must
be included to perform the transformations correctly.

Transformation between real-space and K-space quantities
--------------------------------------------------------

Transformation from real space to a given :math:`k`-point :

.. math::

   \label{eq:RToK}
     S_{\alpha\beta}(\ensuremath{\mathbf{k}}) = \sum_{\ensuremath{\mathbf{R}}} e^{i\ensuremath{\mathbf{k}}\ensuremath{\mathbf{R}}} S_{\alpha\beta}(\ensuremath{\mathbf{R}})

where

.. math:: S_{\alpha\beta}(\ensuremath{\mathbf{R}}) = \left< \phi_\alpha(0) | \phi_\beta(\ensuremath{\mathbf{R}}) \right>

Back transformation:

.. math::

   \label{eq:KToR}
     S_{\alpha\beta}(\ensuremath{\mathbf{R}}) = \sum_{l} \omega_l e^{-i\ensuremath{\mathbf{k}}_l \ensuremath{\mathbf{R}}}
     S_{\alpha\beta}(\ensuremath{\mathbf{k}}_l)

with :math:`\omega_l` being the :math:`k`-point weighting (note: without
symmetrising this sum, all :math:`k`-points must be included, apart from
inversionally related points).

Substituting `[eq:RToK] <#eq:RToK>`__ in `[eq:KToR] <#eq:KToR>`__ yields

.. math::

   \begin{aligned}
     S_{\alpha\beta}(\ensuremath{\mathbf{R}}) &=& \sum_{l} \omega_l e^{-i\ensuremath{\mathbf{k}}_l \ensuremath{\mathbf{R}}}
     \sum_{\ensuremath{\mathbf{R}}'} e^{i\ensuremath{\mathbf{k}}_l\ensuremath{\mathbf{R}}'} S_{\alpha\beta}(\ensuremath{\mathbf{R}}') \\
     &=&\sum_{\ensuremath{\mathbf{R}}'} \sum_{l} \omega_l e^{i\ensuremath{\mathbf{k}}_l(\ensuremath{\mathbf{R}}'-\ensuremath{\mathbf{R}})}
     S_{\alpha\beta}(\ensuremath{\mathbf{R}}') \\
     \label{eq:2}
     &=&\sum_{\ensuremath{\mathbf{R}}'} \delta_{\ensuremath{\mathbf{R}}\ensuremath{\mathbf{R}}'} S_{\alpha\beta}(\ensuremath{\mathbf{R}}') \\
     &=& S_{\alpha\beta}(\ensuremath{\mathbf{R}})
   \end{aligned}

Transition from `[eq:1] <#eq:1>`__ to `[eq:2] <#eq:2>`__ only holds, if
the K-point sampling is good enough:

.. math::

   \delta_{\ensuremath{\mathbf{R}}\ensuremath{\mathbf{R}}'} =  \frac{1}{V_{\text{BZ}}}\int_{\text{BZ}}
     e^{i\ensuremath{\mathbf{k}}(\ensuremath{\mathbf{R}}'-\ensuremath{\mathbf{R}})} d\ensuremath{\mathbf{k}}\approx \sum_{l} \omega_l e^{i\ensuremath{\mathbf{k}}_l(\ensuremath{\mathbf{R}}'-\ensuremath{\mathbf{R}})}

Mulliken analysis in real-space
-------------------------------

For periodic boundary conditions the basis functions are :

.. math::

   \label{eq:basisfuncs}
     \Phi_\alpha^\ensuremath{\mathbf{k}}= \frac{1}{\sqrt{N}} \sum_\ensuremath{\mathbf{R}}\phi_\alpha(\ensuremath{\mathbf{R}}) e^{i\ensuremath{\mathbf{k}}\ensuremath{\mathbf{R}}}

The wavefunction (molecular orbital) :math:`i` at K-point
:math:`\mathbf{k}` is then:

.. math::

   \label{eq:wavefuncdef}
     \psi_i^\ensuremath{\mathbf{k}}= \sum_\alpha c_{i\alpha}^\ensuremath{\mathbf{k}}\Phi_\alpha^\ensuremath{\mathbf{k}}

Charge on that orbital with the Mulliken partition is given by :

.. math::

   \label{eq:chargeperorb}
     q_i^\ensuremath{\mathbf{k}}= n_i^\ensuremath{\mathbf{k}}\left|\psi_i^\ensuremath{\mathbf{k}}\right|^2 = n_i^\ensuremath{\mathbf{k}}\sum_\alpha
     \sum_\beta c_{i\alpha}^{\ensuremath{\mathbf{k}}*} c_{i\beta}^\ensuremath{\mathbf{k}}\left<\Phi_\alpha^\ensuremath{\mathbf{k}}|
     \Phi_\beta^\ensuremath{\mathbf{k}}\right>

.. math::

   \begin{aligned}
     \left<\Phi_\alpha^\ensuremath{\mathbf{k}}| \Phi_\beta^\ensuremath{\mathbf{k}}\right> &=& \frac{1}{N}
     \sum_\ensuremath{\mathbf{R}}\sum_{\ensuremath{\mathbf{R}}'} e^{i\ensuremath{\mathbf{k}}(\ensuremath{\mathbf{R}}'-\ensuremath{\mathbf{R}})} \left<\phi_\alpha(\ensuremath{\mathbf{R}}) |
     \phi_\beta(\ensuremath{\mathbf{R}}') \right> \\ &=& \sum_\ensuremath{\mathbf{R}}e^{i\ensuremath{\mathbf{k}}\ensuremath{\mathbf{R}}}
     \left<\phi_\alpha(0) | \phi_\beta(\ensuremath{\mathbf{R}}) \right> \\ &=& \sum_\ensuremath{\mathbf{R}}
     e^{i\ensuremath{\mathbf{k}}\ensuremath{\mathbf{R}}} S_{\alpha\beta}(\ensuremath{\mathbf{R}})
   \end{aligned}

The total number of electrons (:math:`N_{\text e}`) in the system then
can be calculated as:

.. math::

   \begin{aligned}
     N_{\text e} &=& \sum_l \omega_l q_i^{\ensuremath{\mathbf{k}}_l} \\ &=& \sum_l \omega_l
     \sum_i n_i^{\ensuremath{\mathbf{k}}_l} \sum_\alpha \sum_\beta c_{i\alpha}^{\ensuremath{\mathbf{k}}_l*}
     c_{i\beta}^{\ensuremath{\mathbf{k}}_l}\sum_\ensuremath{\mathbf{R}}e^{i\ensuremath{\mathbf{k}}_l\ensuremath{\mathbf{R}}} S_{\alpha\beta}(\ensuremath{\mathbf{R}}) \\
     &=&\sum_\alpha \sum_l \omega_l \sum_\beta \sum_i n_i^{\ensuremath{\mathbf{k}}_l}
     c_{i\alpha}^{\ensuremath{\mathbf{k}}_l*} c_{i\beta}^{\ensuremath{\mathbf{k}}_l} \sum_\ensuremath{\mathbf{R}}e^{i\ensuremath{\mathbf{k}}_l\ensuremath{\mathbf{R}}}
     S_{\alpha\beta}(\ensuremath{\mathbf{R}})
   \end{aligned}

Introducing the k-dependent density matrix

.. math::

   P_{\alpha\beta}(\ensuremath{\mathbf{k}}) = \sum_i n_i^\ensuremath{\mathbf{k}}c_{i\alpha}^{\ensuremath{\mathbf{k}}*}
     c_{i\beta}^{\ensuremath{\mathbf{k}}}

we get for the charge on an arbitrary orbital :math:`\alpha`:

.. math::

   \begin{aligned}
     q_\alpha &=& \sum_l \omega_l \sum_\beta P_{\alpha\beta}^{\ensuremath{\mathbf{k}}_l}
     \sum_\ensuremath{\mathbf{R}}e^{i\ensuremath{\mathbf{k}}_l\ensuremath{\mathbf{R}}} S_{\alpha\beta}(\ensuremath{\mathbf{R}}) \\ &=& \sum_\ensuremath{\mathbf{R}}
     \sum_\beta S_{\alpha\beta}(\ensuremath{\mathbf{R}}) \sum_l \omega_l
     P_{\alpha\beta}^{\ensuremath{\mathbf{k}}_l} e^{i\ensuremath{\mathbf{k}}_l\ensuremath{\mathbf{R}}} \\
   \end{aligned}

Defining

.. math::

   \Tilde P_{\alpha\beta}(\ensuremath{\mathbf{R}}) = \sum_l \omega_l
     P_{\alpha\beta}^{\ensuremath{\mathbf{k}}_l} e^{i\ensuremath{\mathbf{k}}_l\ensuremath{\mathbf{R}}}

which is, because of the inversional symmetry of the BZ from
time-reversal symmetry, the same as

.. math::

   P_{\alpha\beta}(\ensuremath{\mathbf{R}}) = \sum_l \omega_l P_{\alpha\beta}^{\ensuremath{\mathbf{k}}_l}
     e^{-i\ensuremath{\mathbf{k}}_l\ensuremath{\mathbf{R}}}

(see equation `[eq:KToR] <#eq:KToR>`__), we get for the charge on
orbital :math:`\alpha`:

.. math::

   \begin{aligned}
     q_\alpha &=& \sum_\ensuremath{\mathbf{R}}\sum_\beta S_{\alpha\beta}(\ensuremath{\mathbf{R}}) \Tilde
     P_{\alpha\beta}(\ensuremath{\mathbf{R}}) \\ &=& \sum_{\beta'} S_{\alpha\beta'} \Tilde
     P_{\alpha\beta'} \label{eqn:mulliken}
   \end{aligned}

where the summation runs over the orbitals of the atoms in the central
cell and the image atoms. This practically translates to an element-wise
multiplication between the rectangular shaped real overlap and the
density matrix and a summation over the appropriate row.
