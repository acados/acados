# Getting Started

## problem formulation

```math
\begin{equation}
\begin{aligned}
&\underset{\begin{subarray}{c}
    x(\cdot),\,u(\cdot), \, z(\cdot)
\end{subarray}}{\min}	    &&\int_0^T l(x(\tau), u(\tau), z(\tau), p)\mathrm{d}\tau + m(x(T), z(T), p)\\ 
                            &\,\,\,\quad \text{s.t.}    &&x(0) - \bar{x}_0 = 0, &&\\
                            & 						    &&\underline{x}_0 \leq \Pi_{x_0}x(0) \leq \bar{x}_0, && \\
                            & 						    &&\underline{u}_0 \leq \Pi_{u_0}u(0) \leq \bar{u}_0, && \\
                            & 						    &&\underline{z}_0 \leq \Pi_{z_0}z(0) \leq \bar{z}_0, && \\
                            & 						    &&\underline{c}_0 \leq C_0x(0) + D_0u(0) + E_0z(0)\leq \bar{c}_0, && \\
                            &                           &&                                                   && \\[-1em]
                            & 						    &&F(x(t), \dot{x}(t), u(t), z(t), p) = 0, &&\quad t \in [0,\,T),\\
                            & 						    &&\underline{h} \leq g(h(x(t), u(t), z(t), p)) \leq \bar{h}, &&\quad t \in [0,\,T),\\
                            & 						    &&\underline{x} \leq \Pi_{x}x(t) \leq \bar{x}, &&\quad t \in (0,\,T),\\
                            & 						    &&\underline{u} \leq \Pi_{u}u(t) \leq \bar{u}, &&\quad t \in (0,\,T),\\
                            & 						    &&\underline{z} \leq \Pi_{z}z(t) \leq \bar{z}, && \quad t \in (0,\,T),\\
                            & 						    &&\underline{c} \leq Cx(t) + Du(t) + Ez(t)\leq \bar{c}, &&\quad t \in (0,\,T), \\
                            &                           &&                                                   && \\[-1em]
                            & 						    &&F_T(x(T), z(T), p) = 0, &&\\
                            & 						    &&\underline{h}_T \leq g_T(h_T(x(T), z(T), p)) \leq \bar{h}_T, &&\\
                            & 						    &&\underline{x}_T \leq \Pi_{x_T}x(T) \leq \bar{u}_{T}, &&\\
                            & 						    &&\underline{z}_T \leq \Pi_{z_T}z(T) \leq \bar{z}_T, && \\
                            & 						    &&\underline{c}_T \leq C_Tx(T) + E_Tz(T)\leq \bar{c}_T, &&\\
\end{aligned}
\end{equation}
```
