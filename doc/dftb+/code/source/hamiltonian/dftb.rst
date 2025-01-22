.. _libSKrotations:

SK transforms
=============

SK integral tables are labelled X-Y as files, where interactions are
taken as :math:`l` on atom X interacts with :math:`l^\prime` on atom Y,
example : for O with no :math:`d` orbitals interacting with Si that has
:math:`d`.

==== ================ ======================================
\    :math:`sd\sigma` 
Si-O =0               no-:math:`d` orbital
O-Si :math:`\neq`\ 0  yes there is a :math:`d` orbital on Si
==== ================ ======================================

.. _SKtrans:

SK integral tables
------------------

new total list is (up to :math:`f` orbitals) :

.. math::

   ff^\sigma ff^\pi ff^\delta ff^\phi df^\sigma df^\pi df^\delta
   dd^\sigma dd^\pi dd^\delta pf^\sigma pf^\pi pd^\sigma pd^\pi pp^\sigma
   pp^\pi sf^\sigma sd^\sigma sp^\sigma ss^\sigma

old list is

.. math::

   dd^\sigma dd^\pi dd^\delta pd^\sigma pd^\pi pp^\sigma
   pp^\pi sd^\sigma sp^\sigma ss^\sigma

This loops as with pseudocode above in section `[SKloops] <#SKloops>`__,
and matches the old convention for SK files up to lmax is :math:`d`.

| Within code the table is stored as array:
| skover(row,column,species1,species2) and
  skham(row,column,species1,species2) – Note skover and skham have
  **20** columns for the full :math:`f` case, 9 for :math:`d`, 4 for
  :math:`p` and 1 for :math:`s`, padding everything with zeros when a
  given angular momentum is not present. An indexing array (see
  section\ `[SKindex] <#SKindex>`__) is supplied to translate between
  iSK2(:math:`l_1`,\ :math:`l_2`,\ :math:`|m|`) and elements of
  skham/skover.

symmetry
~~~~~~~~

For :math:`l_1 = l_2` the SK data is symmetric. For heteronuclear pairs,
the elements in the table between equal :math:`l` values are symmetric
between tables for species :math:`A` with :math:`B` and species
:math:`B` with :math:`A`. This leads to there being
:math:`((6 + 11 \times lmax + 6 \times lmax^2 +
lmax^3)/6)` unique elements for the homonuclear case, and
:math:`((6 + 13
\times lmax + 9 \times lmax^2 + 2 \times lmax^3)/6)` in the
heteronuclear case.

The resulting diatomic block also has some symmetry, the
:math:`l_1 = l_2` sub-blocks are symmetric. Additionally, there is
symmetry in the rotational transform for the :math:`l_1 > l_2` upper
triangle and the :math:`l_1 <
l_2` lower triangle, as the only difference in the transforms is the
:math:`-1^{(l_1 - l_2 + | l_1 - l_2 |)/2}` pre-factor.





.. _libperiodic:

lib_periodic
------------

Provides routines to convert between extended sparse :math:`n \times m`
matrices in real space and :math:`n \times n` matrices at arbitrary
:math:`\mathbf{k}`. The aim of the code design is to build all large
matrices in real space where there is a natural cut off as :math:`S
\rightarrow 0`. Arrays of all neighbouring atoms are supplied (iNeighbor
and nNeighbor) to allow construction of these arrays. Only the lower
triangle is constructed. See section `[periodcoords] <#periodcoords>`__
for a note about coordinates. The :math:`k`-dependant folding is done
using appropriate phases (for molecules and at
:math:`\mathbf{k} = (0,0,0)` the phase is always :math:`1+0i`):

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
