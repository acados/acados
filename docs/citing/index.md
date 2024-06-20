# Citing
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


### Advanced-Step Real-Time Iterations (AS-RTI)
Advanced-step real-time iterations provide an extension to the classic real-time iteration algorithm, which allows to performs additional multi-level iterations in the preparation phase, such as inexact or zero-order SQP iterations on a problem with a predicted state estimate.

This feature can be used by setting the options `as_rti_level` and `as_rti_level`.

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