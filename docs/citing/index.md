# Citing

``` eval_rst
.. meta::
   :description: How to cite acados in scientific publications, including BibTeX entries for the main paper and advanced features for nonlinear model predictive control (NMPC) like advanced-step real-time iteration, Gauss-Newton Runge-Kutta integrators, zero-order robust optimization, and multi-phase optimal control problems.
   :keywords: acados, acados citations, acados BibTeX, research papers, acados publications
```

If you are using `acados` in your scientific work, please cite the original journal publication.

```latex
@Article{Verschueren2021,
  Title                    = {acados -- a modular open-source framework for fast embedded optimal control},
  Author                   = {Robin Verschueren and Gianluca Frison and Dimitris Kouzoupis and Jonathan Frey and Niels van Duijkeren and Andrea Zanelli and Branimir Novoselnik and Thivaharan Albin and Rien Quirynen and Moritz Diehl},
  Journal                  = {Mathematical Programming Computation},
  Year                     = {2021},
}
```

We highly appreciate if you share your `acados` success stories with the world, and mention them on the [List of projects that feature `acados`](../list_of_projects/index.md).
Please submit a PR, make a forum post, or contact the developers if you want to showcase your project there!

## Publications on advanced `acados` features
If you are using some of the following advanced features of `acados`, please additionally cite the corresponding publications.

### Fast integrators with sensitivity propagation for use in CasADi
This paper demonstrates the efficiency of `acados` integrators, compared with integrators available in CasADi.
Additionally it describes the [`casados-integrator`](https://github.com/FreyJo/casados-integrators) a wrapper which allows using the `acados` integrators within a `CasADi` NLP solver, like IPOPT.

```latex
@InProceedings{Frey2023,
  Title                    = {Fast integrators with sensitivity propagation for use in {C}as{AD}i},
  Author                   = {Frey, Jonathan and De Schutter, Jochem and Diehl, Moritz},
  Booktitle                = ECC,
  Year                     = {2023},
}
```

### Multi-Phase Optimal Control Problems for Efficient Nonlinear Model Predictive Control with acados
Computationally efficient nonlinear model predictive control relies on elaborate discrete-time optimal control problem (OCP) formulations trading off accuracy with respect to the continuous-time problem and associated computational burden. Such formulations, however, are in general not easy to implement within specialized software frameworks tailored to numerical optimal control. This paper introduces a new multi-phase OCP interface for the open-source software acados allowing to conveniently formulate such problems and generate fast solvers that can be used for nonlinear model predictive control (NMPC). While multi-phase OCP (MOCP) formulations occur naturally in many applications, this paper focuses on MOCP formulations that can be used to efficiently approximate standard continuous-time OCPs in the context of NMPC.


This feature can be used via the `AcadosMultiphaseOcp` class.


```latex
@misc{Frey2024MultiPhase,
      title={Multi-Phase Optimal Control Problems for Efficient Nonlinear Model Predictive Control with acados},
      author={Jonathan Frey and Katrin Baumgärtner and Gianluca Frison and Moritz Diehl},
      year={2024},
      eprint={2408.07382},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2408.07382},
}
```

### Line-Search Funnel-Penalty Globalization
The line-search funnel–penalty method is a globalization strategy that combines 
a funnel mechanism with an ℓ₁ penalty function to determine whether a trial iterate 
is acceptable. If a trial step is rejected, the step size is reduced. The ℓ₁ penalty 
function guarantees that an acceptable step size can always be found.

This feature can be used by setting the option `globalization` to `FUNNEL_L1PEN_LINESEARCH`.

```latex
@PhdThesis{Kiessling2025,
  Title                    = {Algorithmic Advances in Feasible-Iterate and SQP Methods for Direct Optimal Control},
  Author                   = {David Kiessling},
  School                   = {KU Leuven},
  Year                     = {2025}
}
```

### SQP Solver with feasible QPs
This SQP solver guarantees a well-defined search direction at every iteration. 
When the QP is feasible, the solver produces the same iterates 
as the standard SQP solver in acados. If the QP is infeasible, the solver first 
attempts to solve the standard QP; upon failure, it switches to a Byrd–Omojokun 
strategy, which guarantees a well-defined search direction. The solver remains in Byrd–Omojokun 
mode until a heuristic indicates that feasible QPs are likely again. At this point 
the solver switches back to solving the standard QP.

This feature can be used by setting the option `nlp_solver_type` to `SQP_WITH_FEASIBLE_QP`.

```latex
@PhdThesis{Kiessling2025,
  Title                    = {Algorithmic Advances in Feasible-Iterate and SQP Methods for Direct Optimal Control},
  Author                   = {David Kiessling},
  School                   = {KU Leuven},
  Year                     = {2025}
}
```

### Advanced-Step Real-Time Iterations (AS-RTI)
Advanced-step real-time iterations provide an extension to the classic real-time iteration algorithm, which allows to performs additional multi-level iterations in the preparation phase, such as inexact or zero-order SQP iterations on a problem with a predicted state estimate.

This feature can be used by setting the options `as_rti_level` and `as_rti_iter`.

```latex
@Misc{Frey2024a,
  Title                    = {Advanced-Step Real-Time Iterations with Four Levels -- New Error Bounds and Fast Implementation in acados},
  Author                   = {Jonathan Frey and Armin Nurkanovic and Moritz Diehl},
  Year                     = {2024},
  Eprint                   = {2403.07101},
  Primaryclass             = {math.OC},
  Url                      = {https://arxiv.org/abs/2403.07101}
}
```

### Gauss-Newton Runge-Kutta (GNRK) integrators for efficient discretization of OCPs with long horizons and least-squares costs

The GNRK integration scheme can be used by setting the option `cost_discretization = 'INTEGRATOR'`.
This paper additionally demonstrates the effectiveness of using nonuniform discretization grids, and in particular in combination with GNRK.
```latex
@InProceedings{Frey2023b,
  Title                    = {{G}auss-{N}ewton {R}unge-{K}utta Integration for Efficient Discretization of Optimal Control Problems with Long Horizons and Least-Squares Costs},
  Author                   = {Jonathan Frey and Katrin Baumgärtner and Moritz Diehl},
  Booktitle                = {accepted for ECC 2024},
  Url                      = {https://arxiv.org/abs/2310.00618}
}
```

### Efficient Zero-Order Robust Optimization (zoRO) for Real-Time Model Predictive Control with acados
```latex
@InProceedings{Frey2024,
  Title                    = {Efficient Zero-Order Robust Optimization for Real-Time Model Predictive Control with acados},
  Author                   = {Jonathan Frey and Yunfan Gao and Florian Messerer and Amon Lahr and Melanie N Zeilinger and Moritz Diehl},
  Booktitle                = ECC,
  Year                     = {2024},
}
```

### Riccati-ZORO: An efficient algorithm for heuristic online optimization of internal feedback laws in robust and stochastic model predictive control
This feature can be used by setting `zoro_description.feedback_optimization_mode != "CONSTANT_FEEDBACK"`.
```latex
@Misc{Messerer2025,
  Title                    = {Riccati-{ZORO}: An efficient algorithm for heuristic online optimization of internal feedback laws in robust and stochastic model predictive control},

  Author                   = {Florian Messerer and Yunfan Gao and Jonathan Frey and Moritz Diehl},
  Year                     = {2025},

  Archiveprefix            = {arXiv},
  Eprint                   = {2511.10473},
  Primaryclass             = {math.OC},
  Url                      = {https://arxiv.org/abs/2511.10473}
}
```

### Structure exploiting implicit Runge-Kutta method: GNSF
The GNSF IRK integrator be used by setting the option `integrator_tpye = 'GNSF'`.
```latex
@InProceedings{Frey2019,
  Title                    = {Detecting and Exploiting {G}eneralized {N}onlinear {S}tatic {F}eedback
Structures in {DAE} Systems for {MPC}},
  Author                   = {Jonathan Frey and Rien Quirynen and Dimitris Kouzoupis and Gianluca Frison and
Jens Geisler and Axel Schild and Moritz Diehl},
  Booktitle                = ECC,
  Year                     = {2019},
}
```
