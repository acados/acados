# Interfaces


## C API
The C API of acados is an efficient interface to the core functionalities of acados. 
It provides setters and getters that can be used to interact with the core of acados with 
negligible computation overhead. In order to learn about the acados C API, you 
can look at the examples in
[`acados/examples/c/`](https://github.com/acados/acados/tree/master/examples/c). 


## acados MATLAB (rapid prototyping)

This interface makes a broad set of `acados` functionalities available from Matlab or Octave 
for prototyping purpose. As of now, this closely tracks the latest developments in the core of acados, e.g.
exact Hessians, adjoint corrections, regularization, etc. However, for the time being, it will not be possible to 
generate a self-contained C library that can be deployed on an embedded system. For this purpose 
see the `acados emebedded` high-level interface below. 

Some examples for the use of this interface can be found in `<acados_root>/examples/acados_matlab_octave`

This interface uses the shared libraries created using the make command from the main acados folder

```bash
make shared_library
```

To run the examples, navigate into the selected folder, and there run the command
```
export ACADOS_INSTALL_DIR="<acados_root>"
source env.sh # Which can be found in the folder of one of the examples
```

If `ACADOS_INSTALL_DIR` is not specified, it will be assumed that the examples are run from the sub-folders of the current folder (i.e. acados main folder is 2 folders up from the current folder).

Afterwards, launch Matlab or Octave from the same shell.

If you want to run the examples in a different folder, please close the current shell and open a new one to repeat the procedure: this ensures the correct setting of the environment variables.


**Parameters for the Matlab interface**

The following table lists the main input parameters that you have to set for describing your problem and the methods to solve it. The table ends with the main output parameters.

<table>
<thead>
<tr class="header">
<th><strong><span dir="ltr">Input parameter</span></strong></th>
<th><strong><span dir="ltr">Description</span></strong></th>
<th><strong><span dir="ltr">Allowed values</span></strong></th>
<th><strong><span dir="ltr">Other info</span></strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><span dir="ltr">compile_mex</span></td>
<td><span dir="ltr"></span></td>
<td><p><span dir="ltr">true</span></p>
<p><span dir="ltr">false</span></p></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">codgen_model</span></td>
<td><span dir="ltr"></span></td>
<td><p><span dir="ltr">true</span></p>
<p><span dir="ltr">false</span></p></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">gnsf_detect_structure</span></td>
<td><p><span dir="ltr">It allows to detect linear structures in the dynamics to speed up the computations.</span></p>
<p><span dir="ltr"></span></p>
<p><span dir="ltr"><strong>gnsf</strong>: generalized nonlinear static feedback structure</span></p></td>
<td><p><span dir="ltr">true</span></p>
<p><span dir="ltr">false</span></p></td>
<td><span dir="ltr">References: <a href="https://cdn.syscop.de/publications/Frey2019.pdf"><span class="underline">Frey2019.pdf</span></a></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">param_scheme</span></td>
<td><span dir="ltr">It defines which scheme between single and multiple shooting schemes is chosen for solving the optimization problem.</span></td>
<td><p><span dir="ltr">single_shooting</span></p>
<p><span dir="ltr">multiple_shooting</span></p>
<p><span dir="ltr">multiple_shooting_unif_grid</span></p></td>
<td><p><span dir="ltr">References:</span></p>
<p><span dir="ltr"><a href="http://casadi.sourceforge.net/v3.1.0/users_guide/html/node8.html#SECTION00820000000000000000"><span class="underline">Casadi_user_guide</span></a></span></p>
<p><span dir="ltr"></span></p>
<p><span dir="ltr">Single shooting: sequential;</span></p>
<p><span dir="ltr">Multiple shooting: parallel.</span></p>
<p><span dir="ltr"></span></p>
<p><span dir="ltr">If you select ‘multiple_shooting’, then you can set the nodes with the following:</span></p>
<p><span dir="ltr">ocp_opts.set('param_scheme_shooting_nodes', shooting_nodes);</span></p></td>
</tr>
<tr class="odd">
<td><span dir="ltr">nlp_solver</span></td>
<td><p><span dir="ltr">It sets the solver for the NLP problem.</span></p>
<p><span dir="ltr"></span></p>
<p><span dir="ltr"><strong>sqp</strong>: sequential quadratic programming</span></p>
<p><span dir="ltr"><strong>rti</strong>: real time iteration</span></p></td>
<td><p><span dir="ltr">sqp</span></p>
<p><span dir="ltr">sqp_rti</span></p></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">nlp_solver_exact_hessian</span></td>
<td><span dir="ltr"></span></td>
<td><p><span dir="ltr">true</span></p>
<p><span dir="ltr">flase</span></p></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">regularize_method</span></td>
<td><span dir="ltr">Regularization method of the Hessian in the QP subproblems of SQP methods</span></td>
<td><p><span dir="ltr">no_regularize</span></p>
<p><span dir="ltr">project</span></p>
<p><span dir="ltr">project_reduc_hess</span></p>
<p><span dir="ltr">mirror</span></p>
<p><span dir="ltr">convexify</span></p></td>
<td><p><span dir="ltr">References: <a href="https://cdn.syscop.de/publications/Verschueren2017.pdf"><span class="underline">Verschueren2017.pdf</span></a></span></p>
<p><span dir="ltr"></span></p>
<p><span dir="ltr"><a href="https://www.math.uh.edu/~rohop/fall_06/Chapter4.pdf"><span class="underline">Slides_www.math.uh.edu</span></a></span></p>
<p><span dir="ltr"></span></p></td>
</tr>
<tr class="even">
<td><span dir="ltr">nlp_solver_max_iter</span></td>
<td><span dir="ltr">Maximum number of QP iterations allowed for the solution of the SQP</span></td>
<td><span dir="ltr">integer, i.e. 100</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">nlp_solver_tol_stat</span></td>
<td><span dir="ltr">Tolerance for the stationarity condition of the NLP problem</span></td>
<td><span dir="ltr">float, i.e. 1e-4</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">nlp_solver_tol_eq</span></td>
<td><span dir="ltr">Tolerance for the equality constraints of the NLP problem</span></td>
<td><span dir="ltr">float, i.e. 1e-4</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">nlp_solver_tol_ineq</span></td>
<td><span dir="ltr">Tolerance for the inequality constraints of the NLP problem</span></td>
<td><span dir="ltr">float, i.e. 1e-4</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">nlp_solver_tol_comp</span></td>
<td><span dir="ltr">Tolerance for the complementary conditions of the NLP problem</span></td>
<td><span dir="ltr">float, i.e. 1e-4</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">nlp_solver_ext_qp_res</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">qp_solver</span></td>
<td><p><span dir="ltr">It sets the solver for the QP subproblems of the NLP problem.</span></p>
<p><span dir="ltr"></span></p>
<p><span dir="ltr"><strong>hpipm</strong>: high-performance interior point method</span></p>
<p><span dir="ltr"></span></p></td>
<td><p><span dir="ltr">partial_condensing_hpipm</span></p>
<p><span dir="ltr">full_condensing_hpipm</span></p>
<p><span dir="ltr">full_condensing_qpoases</span></p></td>
<td><span dir="ltr"><a href="https://ch.mathworks.com/help/optim/ug/choosing-the-algorithm.html"><span class="underline">Mathworks_opt_algo</span></a>: how to select between (1) interior point, (2) sqp, (3) active set methods</span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">qp_solver_iter_max</span></td>
<td><span dir="ltr">Maximum number of iterations allowed for solving each QP subproblem.</span></td>
<td><span dir="ltr">integer, i.e. 100</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">qp_solver_cond_N</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">qp_solver_cond_ric_alg</span></td>
<td><span dir="ltr"></span></td>
<td><p><span dir="ltr">integer in {0,1}</span></p>
<p><span dir="ltr"></span></p>
<p><span dir="ltr">0: dont factorize hessian in the condensing;</span></p>
<p><span dir="ltr">1: factorize</span></p></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">qp_solver_ric_alg</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr">hpipm-specific</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">qp_solver_warm_start</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr">0: no warm start</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">sim_method</span></td>
<td><p><span dir="ltr">Integration method for the dynamics</span></p>
<p><span dir="ltr"></span></p>
<p><span dir="ltr"><strong>erk:</strong> explicit Runge Kutta</span></p>
<p><span dir="ltr"><strong>irk</strong>: implicit Runge Kutta</span></p></td>
<td><p><span dir="ltr">erk</span></p>
<p><span dir="ltr">irk</span></p>
<p><span dir="ltr">irk_gnsf</span></p></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">sim_method_num_stages</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr">integer in {1,2,4}</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">sim_method_num_steps</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr">integer</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">cost_type</span></td>
<td><span dir="ltr"></span></td>
<td><p><span dir="ltr">linear_ls</span></p>
<p><span dir="ltr">nonlinear_ls</span></p>
<p><span dir="ltr">ext_cost</span></p></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">dyn_type</span></td>
<td><span dir="ltr"></span></td>
<td><p><span dir="ltr">explicit</span></p>
<p><span dir="ltr">implicit</span></p>
<p><span dir="ltr">discrete</span></p></td>
<td><p><span dir="ltr">For ‘erk’:</span></p>
<p><span dir="ltr">ocp_model.set('dyn_type', 'explicit');</span></p>
<p><span dir="ltr">ocp_model.set('dyn_expr_f', model.expr_f_expl);</span></p>
<p><span dir="ltr"></span></p>
<p><span dir="ltr">For ‘irk’:</span></p>
<p><span dir="ltr">ocp_model.set('dyn_type', 'implicit');</span></p>
<p><span dir="ltr">ocp_model.set('dyn_expr_f', model.expr_f_impl);</span></p>
<p><span dir="ltr"></span></p>
<p><span dir="ltr">For ‘discrete’:</span></p>
<p><span dir="ltr">ocp_model.set('dyn_type', 'discrete');</span></p>
<p><span dir="ltr">ocp_model.set('dyn_expr_phi', model.expr_phi);</span></p></td>
</tr>
<tr class="odd">
<td><span dir="ltr">T</span></td>
<td><span dir="ltr">Prediction horizon [seconds]</span></td>
<td><span dir="ltr">double</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">nx</span></td>
<td><span dir="ltr">Dimension of the state vector (x)</span></td>
<td><span dir="ltr">integer</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">nu</span></td>
<td><span dir="ltr">Dimension of the input control vector (u)</span></td>
<td><span dir="ltr">integer</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">ny</span></td>
<td><span dir="ltr">Dimension of the output vector (y)</span></td>
<td><span dir="ltr">integer</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">ny_e</span></td>
<td><span dir="ltr">Dimension of the output vector y at the end time (y_e). It can be different from ny.</span></td>
<td><span dir="ltr">integer</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">nbx</span></td>
<td><span dir="ltr">Dimension of bounding constraints on x</span></td>
<td><span dir="ltr">integer</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">nbu</span></td>
<td><span dir="ltr">Dimension of bounding constraints on u</span></td>
<td><span dir="ltr">integer</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">ng</span></td>
<td><span dir="ltr">Dimension of general affine constraints</span></td>
<td><span dir="ltr">integer</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">ng_e</span></td>
<td><span dir="ltr">Dimension of general affine constraints at the end time</span></td>
<td><span dir="ltr">integer</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">nh</span></td>
<td><span dir="ltr">Dimension of nonlinear constraints</span></td>
<td><span dir="ltr">integer</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">nh_e</span></td>
<td><span dir="ltr">Dimension of nonlinear constraints at the end time</span></td>
<td><span dir="ltr">integer</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">W</span></td>
<td><span dir="ltr">Weight matrix in lagrange term</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">W_e</span></td>
<td><span dir="ltr">Weight matrix in mayer term</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">yr</span></td>
<td><span dir="ltr">Output reference in lagrange term</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">yr_e</span></td>
<td><span dir="ltr">Output reference in mayer term</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">x0</span></td>
<td><span dir="ltr">Initial condition for the state vector</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">Jbx</span></td>
<td><span dir="ltr">Matrix made of rows from the identity matrix. It encodes which of the state variables are bounded.</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">Jbu</span></td>
<td><span dir="ltr">Matrix made of rows from the identity matrix. It encodes which of the input control variables are bounded.</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">lbu</span></td>
<td><span dir="ltr">Lower bound for the control input vector</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">ubu</span></td>
<td><span dir="ltr">Upper bound for the control input vector</span></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
</tbody>
</table>

<span dir="ltr"></span>

<span dir="ltr"></span>

<table>
<thead>
<tr class="header">
<th><strong><span dir="ltr">Output variables</span></strong></th>
<th><strong><span dir="ltr">Description</span></strong></th>
<th><strong><span dir="ltr">Allowed values</span></strong></th>
<th><strong><span dir="ltr">Other info</span></strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><span dir="ltr">stat</span></td>
<td><p><span dir="ltr">Status of convergence.</span></p>
<p><span dir="ltr"></span></p>
<p><span dir="ltr">Composed by [iter, res_g, res_b, res_d, red_m, qp_stat, qp_iter, qp_res_g, qp_res_b, qp_res_d, qp_res_m]</span></p></td>
<td><span dir="ltr"></span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">iter</span></td>
<td><span dir="ltr">Index that counts the iterations of the NLP.</span></td>
<td><span dir="ltr">integer, from 0 up to a maximum of ‘nlp_solver_max_iter-1’</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">res_g</span></td>
<td><span dir="ltr">Residuals on stationarity</span></td>
<td><span dir="ltr">below threshold given in input</span></td>
<td><span dir="ltr"> </span></td>
</tr>
<tr class="even">
<td><span dir="ltr">res_b</span></td>
<td><span dir="ltr">Residuals on equality constraint feasibility</span></td>
<td><span dir="ltr">below threshold given in input</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">res_d</span></td>
<td><span dir="ltr">Residuals on inequality constraint feasibility</span></td>
<td><span dir="ltr">below threshold given in input</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="even">
<td><span dir="ltr">res_m</span></td>
<td><span dir="ltr">Residuals on complementary conditions</span></td>
<td><span dir="ltr">below threshold given in input</span></td>
<td><span dir="ltr"></span></td>
</tr>
<tr class="odd">
<td><span dir="ltr">qp_res_*</span></td>
<td><span dir="ltr">Same as before, but for the QP subproblems</span></td>
<td><span dir="ltr">below threshold given in input</span></td>
<td><span dir="ltr"></span></td>
</tr>
</tbody>
</table>

<span dir="ltr"></span>


## acados embedded - Python


`acados_template` is a Python package that can be used to specify optimal control problems from Python and to generate self-contained C code that uses the acados solvers to solve them.
In comparison with the MATLAB interface for rapid prototyping (see above), it supports less features, but it allows the user to generate a self-contained C library  
that can be easily deployed on an embedded system.

The framework is based on templated C files which are rendered from Python using the templating engine `Jinja2`.

The currently supported formulations reads as

```math
\begin{equation}
\begin{aligned}
&\underset{\begin{subarray}{c}
    x(\cdot),\,u(\cdot), \, z(\cdot)
\end{subarray}}{\min}	    &&\int_0^T l(x(\tau), u(\tau), z(\tau), p)\mathrm{d}\tau + m(x(T), z(T), p)\\ 
&\,\,\,\quad \text{s.t.}    &&x(0) - \bar{x}_0 = 0, &&\\
&&&F(x(t), \dot{x}(t), u(t), z(t), p) = 0, &&\quad t \in [0,\,T),\\
&&&\underline{h} \leq h(x(t), u(t), p) \leq \bar{h}, &&\quad t \in [0,\,T),\\
&&&\underline{x} \leq \Pi_{x}x(t) \leq \bar{x}, &&\quad t \in [0,\,T),\\
&&&\underline{u} \leq \Pi_{u}u(t) \leq \bar{u}, &&\quad t \in [0,\,T),\\
&&&\underline{c} \leq Cx(t) + Du(t)\leq \bar{c}, &&\quad t \in [0,\,T), \\
&&& && \\[-1em]
&&&\underline{h}^e \leq h^e(x(T), p) \leq \bar{h}^e, &&\\
&&&\underline{x}^e \leq \Pi_{x}^e x(T) \leq \bar{x}^{e}, &&\\
&&&\underline{c}^e \leq C^e x(T)\leq \bar{c}^e, &&\\
\end{aligned}
\end{equation}
```
```eval_rst
Where:

* :math:`l: \mathbb{R}^{n_x}\times\mathbb{R}^{n_u}\times\mathbb{R}^{n_z} \rightarrow \mathbb{R}` is the Lagrange objective term.
* :math:`m: \mathbb{R}^{n_x}\times\mathbb{R}^{n_z} \rightarrow \mathbb{R}` is the Mayer objective term.

* :math:`F: \mathbb{R}^{n_x}\times\mathbb{R}^{n_x}\times\mathbb{R}^{n_u}\times\mathbb{R}^{n_z}\times\mathbb{R}^{n_p} \rightarrow \mathbb{R}^{n_x+n_z}` describes the (potentially) fully implicit dynamics.

* :math:`h: \mathbb{R}^{n_x}\times\mathbb{R}^{n_u}\times\mathbb{R}^{n_z}\times\mathbb{R}^{n_p} \rightarrow \mathbb{R}^{n_h}` and :math:`h^e: \mathbb{R}^{n_x}\times\mathbb{R}^{n_z}\times\mathbb{R}^{n_p} \rightarrow \mathbb{R}^{n_{h_e}}` are general nonlinear functions.

* :math:`C,\,D,\,C^e,\,\Pi_x,\,\Pi_u,\,\Pi_x^e` are matrices of appropriate dimensions defining the polytopic and box constraints.

The `acados_template` interface makes some limiting assumptions on the problem formulation above, namely:
* :math:`l` must be in linear least-squares form :math:`l = \frac{1}{2}\| V_x x(t) + V_u u(t) + V_z z(t) - y_{\text{ref}}\|_W^2`
* :math:`m` must be in linear least-squares form :math:`m = \frac{1}{2}\| V^e_x x(t) - y_{\text{ref}}^e\|_{W^e}^2`
* Constraints cannot depend on algebraic variables (yet)
```


### Installation
1. Compile and install acados by running:
```bash
cd <acados_root>/build
cmake -DACADOS_WITH_QPOASES=ON ..
make install -j4
```

2. Install acados_template Python package by running
```
pip3 install <acados_root>/interfaces/acados_template
```

(Notice that, of course, you might need to use `pip` instead, if you run, for example, from within a Python virtual 
environment) You should now be able to import it as a Python module and use it as shown in the examples in `<acados_root>/examples/acados_template/python/<example_name>/generate_c_code.py`

For more information contact `@zanellia`

### optimal control description
``` eval_rst
.. automodule:: acados_template.casadi_functions
    :members:
    :private-members:
    :undoc-members:
```
``` eval_rst
.. automodule:: acados_template.acados_ocp_nlp
    :members:
    :private-members:
    :exclude-members: acados_ocp2json_layout, cast_ocp_nlp, dict2json_layout, dict2json_layout_rec, check_ra, json2dict_rec

```
## acados embedded - MATLAB
TODO!
