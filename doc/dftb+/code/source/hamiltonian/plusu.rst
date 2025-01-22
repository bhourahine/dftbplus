DFTB+U
======

Lets first look at the LDA+U approximation to SIC, i.e. the Fully
localised limit (FLL), which is appropriate to lanthanides. Here the
energy correction is given as :math:`-\frac{U - J}{2} ( n_i^2 - n_i )`
for the states considered. Since :math:`0 \geq n_i \leq 1`, this
contribution is always positive. Starting from the expressions in
Petukhov (PRB-67-153106) for the FLL (eqn. 3 in the paper, with the
definitions for :math:`n_i` and :math:`N` substituted in)

.. math::

   \begin{aligned}
     \Delta E & = & -\frac{U-J}{2} \sum_\sigma \left[ Tr
       (\rho^\sigma\cdot\rho^\sigma - \rho^\sigma ) \right]\\
     \Delta V^{\sigma}_{m m^\prime} & = & - ( U - J )\left[ \rho_{m
         m^\prime}^\sigma - \frac{\delta_{m m^\prime}}{2} \right]
   \end{aligned}

Where :math:`\sigma` is a spin index and :math:`(U-J)` is the
re-normalised Coulomb term. The density matrices are restricted to only
the on-site blocks with the appropriate angular momenta, and since these
are hermitian matrices, the diagonal elements are always real.

These expression can be evaluated in the real-space formulation.

Energy
------

| This consists of two terms, which for atom :math:`a` are :
| 1\ :math:`^\mathrm{st}` term

  .. math::

     \begin{aligned}
       Tr^a( \rho^\sigma\cdot\rho^\sigma ) & = & \sum_{ij \in a} \sum_\mathbf{k}
       \omega_\mathbf{k} \rho^\sigma_{ij}(\mathbf{k})
       \rho^\sigma_{ji}(\mathbf{k})\\
       & = & \sum_{ij \in a} \sum_\mathbf{k} \omega_\mathbf{k}
       \rho^\sigma_{ij}(\mathbf{k}) \rho^\sigma_{ij}(\mathbf{k})^*\\
       & = & \sum_{ij \in a} \sum_\mathbf{k} \omega_\mathbf{k}
       \left(\sum_\mathbf{R}\rho^\sigma_{ij}(\mathbf{R})
       \exp(i\mathbf{k}\cdot\mathbf{R})\right)
       \left(\sum_\mathbf{R^\prime}\rho^\sigma_{ij}(\mathbf{R^\prime})
       \exp(-i\mathbf{k}\cdot\mathbf{R^\prime})\right)\\
       & = & \sum_{ij \in a} \sum_\mathbf{k} \omega_\mathbf{k}
       \left(\sum_\mathbf{RR^\prime} \rho^\sigma_{ij}(\mathbf{R})
       \rho^\sigma_{ij}(\mathbf{R^\prime})
       \exp(i\mathbf{k}\cdot(\mathbf{R}-\mathbf{R^\prime}))\right)\\
       & = & \sum_{ij \in a} \sum_\mathbf{RR^\prime} \rho^\sigma_{ij}(\mathbf{R})
       \rho^\sigma_{ij}(\mathbf{R^\prime}) \sum_\mathbf{k}
       \omega_\mathbf{k}
       \exp(i\mathbf{k}\cdot(\mathbf{R}-\mathbf{R^\prime}))
     \end{aligned}

  For a ’good’ set of :math:`k`-points, the final summation becomes
  :math:`\delta_{0,(\mathbf{R}-\mathbf{R^\prime})}`, hence only
  :math:`\mathbf{R} =
  \mathbf{R^\prime}` contributes and hence

  .. math::

     \begin{aligned}
       Tr^a( \rho^\sigma\cdot\rho^\sigma ) & = & \sum_{ij \in a} \sum_\mathbf{R}
       \rho^\sigma_{ij}(\mathbf{R})^2
     \end{aligned}

The second part of the expression is given by

.. math::

   \begin{aligned}
     Tr^a( \rho^\sigma ) & = & \sum_{i \in a} \sum_\mathbf{k}
     \omega_\mathbf{k} \rho^\sigma_{ii}(\mathbf{k})\\ & = & \sum_{i \in
     a} \sum_\mathbf{k} \omega_\mathbf{k} \sum_\mathbf{R}
     \rho^\sigma_{ii}(\mathbf{R}) \exp(i\mathbf{k}\cdot\mathbf{R})\\ & =
     & \sum_{i \in a} \sum_\mathbf{R} \rho^\sigma_{ii}(\mathbf{R})
     \sum_\mathbf{k} \omega_\mathbf{k} \exp(i\mathbf{k}\cdot\mathbf{R})
   \end{aligned}

and again, the final summation is :math:`\delta_{0,\mathbf{R}}`. Hence
The total energy expression is

.. math::

   \begin{aligned}
     Tr^a( \rho^\sigma ) & = & \sum_{i \in a} \sum_{\mathbf{R}=0}
     \rho^\sigma_{ii}(\mathbf{R})
   \end{aligned}

Giving a total energy contribution for the :math:`l`-th angular momentum
shell of all atoms :math:`a` of a given type :

.. math::

   \begin{aligned}
     \Delta E_l & = & -\frac{U_l-J_l}{2} \sum_a \sum_\sigma \left[ \left(
     \sum_\mathbf{R} \sum_{ij \in l \in a} \rho^\sigma_{ij}(\mathbf{R})^2
     \right) - \left( \sum_{i \in l \in a} \rho^\sigma_{ii}(\mathbf{0})
     \right)\right]
   \end{aligned}
