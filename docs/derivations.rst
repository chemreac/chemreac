Derivations
===========

Explicit formulae have been derived for the Smoluchowski equation
in isotropic media for different geometries and for its logarithmically
transformed version (concentration, time, space).


Law of mass action
------------------

The well established law of mass action give the rate of change of a
concentration  :math:`c_i` of species :math:`i`:

.. math ::

    \frac{\partial c_i}{\partial t} = \sum_l r_l S_{il}

:math:`t` is time, :math:`c_i` is the concentration of species :math:`i`, :math:`S_{il}` is the net stoichiometric coefficient of species :math:`i` in reaction :math:`l` and :math:`r_l` is the rate of reaction :math:`l`, which for a mass-action type of rate law is given by:

.. math ::

    r_l = \begin{cases} \kappa_l\prod_k c_k^{R_{kl}} &\mbox{if } \sum_k R_{kl} > 0 \\
    0 &\mbox{otherwise} \end{cases}

where :math:`\kappa_l` is the rate constant, :math:`R_{kl}` is the stoichiometric coefficient of species :math:`k` on the reactant side.

If we introduce the logarithmically transormed concentration
:math:`z`:

.. math ::

    z_i &= \log(c_i)

we have:

.. math ::

    \frac{\partial z_i}{\partial t} &= \frac{\frac{\partial c_i}{\partial t}}{c_i}

which can be expressed in :math:`z_i`:

.. math ::

    \frac{\partial z_i}{\partial t} &= e^{-z_i} \sum_l r_l S_{il}
    
where we may now express :math:`r_l` as:

.. math ::

    r_l = \begin{cases} \kappa_l e^{\sum_k R_{kl} z_k} &\mbox{if } \sum_k R_{kl} > 0 \\
        0 &\mbox{otherwise} \end{cases}
