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
