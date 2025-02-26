Non-collinear spin
==================

Notation
--------

For this section, I will denote spatial vectors as :math:`\vec{x}`,
:math:`2
\times 2` spin-matrices as :math:`\boldsymbol{x}`, and regular scalars
as :math:`x`. For spacial vectors, components will be labelled in Roman
letters (:math:`x`, :math:`y`, :math:`z`), but for spin-matrices they
will be in Greek (:math:`\alpha`, :math:`\beta`). For other general
matrix elements I’ve labelled them as :math:`i`, :math:`j`, :math:`k`,
etc.

I’ve superscripted the spin, and subscripted the spatial/general for
convenience, but this is *not* tensorial notation.

As always, atoms are labelled by Roman letters, a, b, etc.

Wavefunction
------------

Starting from a set of single-particle 2-component spinor wavefunctions,
the *ansatz* is made that they can be expanded in a basis as:
:math:`\psi_i(\vec{r}) = \sum_j \begin{pmatrix} c^\alpha_{ij} \\
c^\beta_{ij} \end{pmatrix} \phi_j(\vec{r})`, where :math:`\phi` are
spatial basis functions, and :math:`c` being expansion coefficients for
the basis function in up and down spin. The essential assumption is that
the spin is constant for each spatial basis function in the expansion,
i.e. :math:`c^\alpha_{ij}` and :math:`c^\beta_{ij}` are constants for
given values of :math:`i` and :math:`j`. Hence the spin is not
completely free to vary in space as with the modern DFT implementations.

The overlap matrix is then defined purely in spacial terms:
:math:`S_{ij} =
\langle \phi_i(\vec{r}) | \phi_j(\vec{r}) \rangle`, so it is only
defined for spin-diagonal parts.

From these functions, the single particle density matrix can be defined
as

.. math::

   \begin{aligned}
   \boldsymbol{\rho} &=& \sum_i n_i |\psi_i\rangle\langle\psi_i| \\
   \boldsymbol{\rho}^{\alpha\beta}_{ij} &=& \sum_k n_k c^\alpha_{ki} c^\beta_{kj}
   \end{aligned}

Hence this is a hermitian matrix. It is convenient to organise these
matrices using spin super-blocks in the hamiltonian
(Fig `5 <#ncspin>`__).

Pauli matrices
--------------

The charge, :math:`\rho(\vec{r})`, and magnetisation,
:math:`\vec{m}(\vec{r})`, parts of a 2-component spin density matrix (so
the ’Greek’ part) can be extracted by expanding in terms of Pauli
matrices Matrices (Landau 1986):

.. math::

   \begin{aligned}
   \boldsymbol{\rho} = \rho \mathbf{1} + \vec{m} \cdot \boldsymbol{\sigma}
   \end{aligned}

where the Pauli matrices are in the usual form

.. math::

   \begin{aligned}
   \mathbf{1} &=&  \begin{pmatrix}1&0\\0&1\end{pmatrix} \\
   \boldsymbol{\sigma}_x &=&  \begin{pmatrix}0&1\\1&0\end{pmatrix} \\
   \boldsymbol{\sigma}_y &=&  \begin{pmatrix}0&i\\-i&0\end{pmatrix} \\
   \boldsymbol{\sigma}_z &=&  \begin{pmatrix}1&0\\0&-1\end{pmatrix}
   \end{aligned}

with respect to spin-quantisation axes. Hence, since an arbitrary
hermitian matrix can be represented as a linear combination of these
matrices, we can use the well known results
that (Sticht JPCM 1 8155)

.. math::

   \begin{aligned}
   \rho &=& (\boldsymbol{\rho}^{\alpha\alpha} +
   \boldsymbol{\rho}^{\beta\beta}) \\ \tan \phi &=& -\frac{Im
   \boldsymbol{\rho}_{12}}{Re \boldsymbol{\rho}_{12}}\\ \tan \theta &=& 2
   [ (Im \boldsymbol{\rho}_{12})^2 + (Re \boldsymbol{\rho}_{12})^2 ]^{1/2} /
   (\boldsymbol{\rho}_{11} - \boldsymbol{\rho}_{22})
   \end{aligned}

Where the angles :math:`\theta` and :math:`\phi` are polar directors in
real space. Alternatively an arbitrary Pauli matrix can be broken down
into a linear combination of the 4 :math:`2 \times 2` matrices. If the
coefficients of this combination are purely real or imaginary (as is the
case here), the coefficients give the directions in real space.

.. math::

   \begin{aligned}
   \rho &=& (\boldsymbol{\rho}^{\alpha\alpha} +
   \boldsymbol{\rho}^{\beta\beta}) \\ m_x &=& Re(\boldsymbol{\rho}_{12} +
   \boldsymbol{\rho}_{21}) \\ m_y &=& Im(\boldsymbol{\rho}_{12} -
   \boldsymbol{\rho}_{21})\\ m_z &=& \boldsymbol{\rho}_{11} -
   \boldsymbol{\rho}_{22}\\
   \end{aligned}

Hence the Pauli-like Kohn-Sham hamiltonian looks
like (Nordstrom PRL 76 442):

.. math::

   \begin{aligned}
   H &=& \{ \nabla^2 + v(\vec{r}) \} \mathbf{1}
    + \{ \vec{b}(\vec{r}) + \xi \vec{l} \} \cdot \boldsymbol{\sigma}
   \end{aligned}

where :math:`\vec{b}` is the exchange correlation, and
:math:`\xi\cdot\vec{l}` a spin-orbit term. The part of the exchange
correlation that we can evaluate is locally collinear to the
magnetisation (Sticht JPCM 1 8155).

.. _`sec:Pauli`:

Symmetry of the Pauli-hamiltonian
---------------------------------

Using molecular boundary conditions for simplicity (i.e. no k-points and
hence phase factors). The usual Pauli matrices are:

.. math::

   \begin{aligned}
     \mathbf{1} &=& \begin{pmatrix}1&0\\0&1\end{pmatrix}
       \\ \boldsymbol{\sigma}_x &=& \begin{pmatrix}0&1\\1&0\end{pmatrix}
         \\ \boldsymbol{\sigma}_y
         &=& \begin{pmatrix}0&-i\\i&0\end{pmatrix}
           = \begin{pmatrix}0&i^{\dag}\\i&0\end{pmatrix}
             \\ \boldsymbol{\sigma}_z
             &=& \begin{pmatrix}1&0\\0&-1\end{pmatrix}
   \end{aligned}

The real part of the hamiltonian (:math:`2 \times 2` shown for
simplicity) is hermitian if its coefficients are symmetric:

.. math::

   \begin{aligned}
     \Re H &=& \begin{pmatrix}H_{aa}&H_{ab}\\H_{ab}&H_{bb}\end{pmatrix}
   \end{aligned}

Its products with the Pauli matrices are then:

.. math::

   \begin{aligned}
     \Re H \otimes \mathbf{1}
     &=& \begin{bmatrix} \begin{pmatrix}H_{aa}&H_{ab}\\ H_{ab}&H_{bb}\end{pmatrix}
       & 0 \\ 0
       & \begin{pmatrix}H_{aa}&H_{ab}\\H_{ab}&H_{bb}\end{pmatrix}
     \end{bmatrix}\\ \Re H
     \otimes \boldsymbol{\sigma}_x&=&
     \begin{bmatrix} 0
       & \begin{pmatrix}H_{aa}&H_{ab}\\ H_{ab}&H_{bb}\end{pmatrix}
       \\ \begin{pmatrix}H_{aa}&H_{ab}\\H_{ab}&H_{bb}\end{pmatrix} &
       0
     \end{bmatrix}\\ \Re H \otimes
     \boldsymbol{\sigma}_y&=& \begin{bmatrix} 0 &
       i^{\dag}\begin{pmatrix}H_{aa}&H_{ab}\\ H_{ab}&H_{bb}\end{pmatrix}
       \\ i\begin{pmatrix}H_{aa}&H_{ab}\\H_{ab}&H_{bb}\end{pmatrix}
       & 0
     \end{bmatrix}\\
     \Re H \otimes
     \boldsymbol{\sigma}_z&=& \begin{bmatrix} \begin{pmatrix}H_{aa}&H_{ab}\\
         H_{ab}&H_{bb}\end{pmatrix} & 0 \\ 0 &
       -\begin{pmatrix}H_{aa}&H_{ab}\\H_{ab}&H_{bb}\end{pmatrix}
     \end{bmatrix}
   \end{aligned}

All matrices are both overall hermitian and the coefficients are
block-wise hermitian (but it is only the symmetric part that matters for
molecular boundary conditions, in the case of k-points this should be
correctly treated as hermitian symmetry). Hence *symmetrising* the
blocks will leave the real coefficient part of the block unchanged,
while removing (as we’ll come to in a moment), the imaginary
coefficients.

Moving now to an imaginary hamiltonian (spin-orbit etc.), we need
products of Pauli matrices with an imaginary pre-factor:

.. math::

   \begin{aligned}
     i\mathbf{1}
     &=& \begin{pmatrix}i&0\\0&i\end{pmatrix}\\ i\boldsymbol{\sigma}_x
       &=& \begin{pmatrix}0&i\\i&0\end{pmatrix} \\ i\boldsymbol{\sigma}_y
         &=& \begin{pmatrix}0&1\\-1&0\end{pmatrix}
           \\ i\boldsymbol{\sigma}_z
           &=& \begin{pmatrix}i&0\\0&-i\end{pmatrix}
             = \begin{pmatrix}i&0\\0&i^{\dag}\end{pmatrix}
   \end{aligned}

The coefficients for the imaginary part of the hamiltonian must be skew
symmetric to be hermitian (i.e., :math:`H_{ab} = H_{ba}^{\dag}`,
:math:`\therefore \Im H_{ab} = - \Im H_{ba}`,
:math:`H_{aa} = - H_{aa} = 0`, :math:`H_{bb} = - H_{bb} = 0`)

.. math::

   \begin{aligned}
     \Im H &=& \begin{pmatrix}0&H_{ab}\\-H_{ab}&0\end{pmatrix}
   \end{aligned}

However, products with the Pauli matrices for an imaginary :math:`H` are
still overall hermitian:

.. math::

   \begin{aligned}
     \Im H \otimes i \mathbf{1} &=& \begin{bmatrix}
       i \begin{pmatrix}0&H_{ab}\\ -H_{ab}&0\end{pmatrix} & 0
       \\ 0 & i \begin{pmatrix}0&H_{ab}\\-H_{ab}&0\end{pmatrix}
     \end{bmatrix}\\ \Im H
     \otimes i \boldsymbol{\sigma}_x&=&
     \begin{bmatrix} 0
       & i \begin{pmatrix}0&H_{ab}\\ -H_{ab}&0\end{pmatrix}
       \\ i \begin{pmatrix}0&H_{ab}\\-H_{ab}&0\end{pmatrix} & 0
     \end{bmatrix}\\ \Im H \otimes i
     \boldsymbol{\sigma}_y &=& \begin{bmatrix} 0 &
       \begin{pmatrix}0&H_{ab}\\ -H_{ab}&0\end{pmatrix}
       \\ -\begin{pmatrix}0&H_{ab}\\-H_{ab}&0\end{pmatrix} & 0
     \end{bmatrix} \nonumber\\
     &=& \begin{bmatrix} 0 &
       \begin{pmatrix}0&H_{ab}\\ -H_{ab}&0\end{pmatrix}
       \\ \begin{pmatrix}0&-H_{ab}\\H_{ab}&0\end{pmatrix} & 0
     \end{bmatrix}\\
     \Im H \otimes i \boldsymbol{\sigma}_z&=& \begin{bmatrix}
       i \begin{pmatrix}0&H_{ab}\\ -H_{ab}&0\end{pmatrix} & 0
       \\ 0 &
       i^{\dag} \begin{pmatrix}0&H_{ab}\\-H_{ab}&0\end{pmatrix}
     \end{bmatrix}\\
   \end{aligned}

but the individual coefficient blocks are then skew symmetric (i.e., not
anti-hermitian).

Mulliken analysis
-----------------

The form for Mulliken analysis is a straightforward generalisation of
the spin-less Mulliken analysis :

.. math::

   \begin{aligned}
   \boldsymbol{q}^{\alpha\beta}_{i \in a} &=& \sum_j S_{ij}
   \boldsymbol{\rho}^{\alpha\beta}_{ij} \\
   \end{aligned}

Where the vectoral :math:`\boldsymbol{q}` can be split into the charge
and magnetisation parts as above, giving :math:`q` and :math:`\vec{m}`.

DFTB total energy
-----------------

Since we know that the exchange-correlation potential must be locally
parallel to the spin direction (or at least the part that can be
evaluated with existing functionals), this allows us to write the
non-collinear SDFTB energy (without spin-orbit or external fields) as
:

.. math::

   \begin{aligned}
   E &=& Tr(\rho H^0) + \sum_{ab}\sum_{i\in a,j \in b} \gamma_{ij} \Delta
   q_i \Delta q_j + \sum_a \sum_{l\in a, l^\prime\in a} W_{all^\prime}
   \vec{p}_{al} \cdot \vec{p}_{al^\prime} \nonumber \\ && +
   E_\mathrm{rep.}
   \end{aligned}

with

.. math::

   \begin{aligned}
   Tr(\rho H^0) &=& \sum_{\alpha\beta\;i\;j\;k} n_i
   c^\alpha_{ij} c^\beta_{ik} H^0_{jk} \delta_{\alpha\beta} \nonumber \\
   &=& \sum_{i\;j\;k} n_i (c^\alpha_{ij} c^\alpha_{ik} + c^\beta_{ij}
   c^\beta_{ik} ) H^0_{jk}\\ q_j &=& \sum_{\alpha \beta i k} n_i S_{jk}
   c^\alpha_{ij} c^\beta_{ik} \delta_{\alpha \beta} \nonumber \\ &=&
   \sum_{i k} n_i S_{jk}
   (c^\alpha_{ij} c^\alpha_{ik} + c^\beta_{ij} c^\beta_{ik})\\ p_{jx} &=&
   \sum_{i\;k} n_i Re (c^\alpha_{ij} c^\beta_{ik} + c^\beta_{ij}
   c^\alpha_{ik} ) S_{jk}\\ p_{jy} &=& \sum_{i\;k} n_i Im (c^\alpha_{ij}
   c^\beta_{ik} - c^\beta_{ij} c^\alpha_{ik} ) S_{jk}\\ p_{jz} &=&
   \sum_{i\;k} n_i (c^\alpha_{ij} c^\alpha_{ik} - c^\beta_{ij}
   c^\beta_{ik} ) S_{jk}%\\
   %\sum_a \sum_{l\in a, l^\prime\in a} W_{all^\prime}
   %\vec{p}_{al} \cdot \vec{p}_{al^\prime} &=& \nonumber
   \end{aligned}

DFTB\ :math:`^{\text{+}}` data structures
-----------------------------------------

Since we now find that we have to represent hermitian matrices in the
packed format (again still controlled by the overlap matrix cut-off), I
propose that the current sparse format is extended by increasing the
second index range for the spin to 4 and storing data as in
figure `6 <#ncsparse>`__.There are several points to note about the
figure:

#. The overall matrix must be hermitian, so the sparse matrices must be
   complex.

#. The blue elements are purely real in the absence of spin-orbit
   coupling, so there would be some saving in using a more complicated
   data structure with a separate real array for the on-diagonal and a
   complex array for the off diagonal.

#. The red/green elements are complex.

#. The same neighbour map can be used in all four blocks, as the
   functions have the same spacial overlap within their block (:math:`S`
   is spin independent, but doesn’t mix different spins together).

#. The green elements are redundant, as they map on-to existing red
   elements in the complex conjugate triangle, but these are only
   on-site and hence don’t waste much.

#. The lower triangle including the green elements should perhaps be
   stored as though it were in the upper block, i.e. the complex
   conjugate of the elements.
