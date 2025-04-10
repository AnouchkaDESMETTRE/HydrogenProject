This README documents our review of the article "Absorption and desorption of hydrogen in long metal hydride tank equipped with phase change material jacket" (Rabienataj Darzi et al., International Journal of Hydrogen Energy, April 2016) and outlines our project plan to model similar aspects of hydrogen storage.

#  Review of "Absorption and desorption of hydrogen in long metal hydride tank equipped with phase change material jacket" 

The article titled "**Absorption and desorption of hydrogen in long metal hydride tank equipped with phase change material jacket**" was published in the International Journal of Hydrogen Energy in **April 2016** and was authored by **A. Ali Rabienataj Darzi, H. Hassanzadeh Afrouzi, A. Moshfegh, and M. Farhadi**. We are reviewing this study to understand its findings for application in our own research or project.

This study numerically investigated the **hydrogen absorption and desorption processes** in a long tubular **LaNi5 metal hydride tank (MHT)** integrated with a **Rubitherm phase change material (PCM) jacket** for PEM fuel cell hydrogen supply. The authors analyzed the influence of different **H2 supply pressures (10, 15, and 20 bar)**, different **discharge pressures (1.5, 1.75, and 2 bar)**, and **metal hydride bed porosities (0.4, 0.5, and 0.6)** on transient and local temperature distributions across the H2-MHT system and PCM jacket. The time-dependent changes of the **hydrogen to metal (H/M) ratio** and **PCM melt fraction** were also investigated until equilibrium was reached.

Key findings of the study include:

*   It was found that **system temperature, PCM melt fraction, and H/M ratio reach steady state at different rates**. Systems with **higher supply pressure in absorption, lower discharge pressure in desorption, and higher bed porosity** approach steady state faster.
*   The **MHT charges with hydrogen much faster under high supply pressures** and **discharges much faster under lower discharge pressures**.
*   Inserting **metal foam in the PCM jacket enhances thermal conductivity** and **significantly reduces the charging and discharging time**.
*   **PCM melt fraction and H/M ratio reached steady state at different rates during absorption**. The H/M ratio approaches unity (full saturation) quicker at high supply pressures, indicating a shorter charging time.
*   The **full discharge times** during desorption occurred at 300, 420, and 600 minutes for discharge pressures of 1.5, 1.75, and 2 bars, respectively. Systems with lower discharge pressure approach steady state faster.
*   Increasing the **bed porosity led to faster hydrogen charging** but a **decrease in H/M ratio during desorption**. Higher bed porosity leads to a lower PCM melt fraction but faster equilibrium. Lower porosity means a higher amount of hydride mass requiring more heat for release during desorption.
*   Inserting metal foam in the PCM jacket enhances the melting and solidification rate of the PCM, desirably reducing the charging and discharging times.

The study concluded that the operating parameters and bed properties significantly affect system performance, and that integrating metal foam into the PCM is an effective strategy for improving heat transfer and reducing cycle times.

---
# Project Plan

Our work will be structured in three main phases, each focusing on a critical aspect of hydrogen storage modeling, drawing inspiration from the reviewed study and considering its conditions and results.

**1/ Convection-Diffusion-Reaction Modeling for Hydrogen Gas and Metal-Oxide Solid Density**

* Physics of the Problem:
    * We will model the evolution of hydrogen gas concentration (ρ_g) within the porous bed and the concentration of hydrogen absorbed in the metal hydride tank (ρ_s). This initial model will focus on transport and reaction in a transient regime. We will start with a simplified approach, not fully coupling the temperature.
    * This model aims to capture the fundamental behavior of hydrogen absorption and desorption, including diffusion, convection, and reaction kinetics.

* Partial Differential Equations (PDEs) to Solve:
    * **For the density of the hydrogen gas (ρ_g):**
        ```
        ε ∂ρ_g/∂t = D ∂²ρ_g/∂z² + u_z ∂ρ_g/∂z + ṁ(ρ_s, t)
        ```
        given ρ_g(t=0) = ρ_g,0 and boundary conditions.
        Where:
        * `ε` is the porosity of the bed.
        * `t` is time.
        * `D` is the diffusion coefficient of hydrogen in the porous bed.
        * `u_z` is the velocity of the gas in the z-direction.
        * `ṁ(ρ_s, t)` is a source/sink term representing the absorption/desorption reaction.
    * **For the density of hydrogen in the solid metal-oxide (ρ_s):**
        ```
        (1 - ε) ∂ρ_s/∂t = ṁ(ρ_s, t)
        ```
        given ρ_s(t=0) = ρ_s,0 and boundary conditions.
        Where:
        * `ρ_s` represents the density of hydrogen in the solid.
        * `ṁ(ρ_s, t)` is the reaction term, coupling the two equations. A simplified reaction rate law (e.g., first-order) will be used initially.
          
* Boundary and Initial Conditions (Absorption case):
    * Boundary Conditions:
        * ρ_s requires no boundary conditions because the PDE for ρ_s is first-order in time and zero-order in space.
        * At the tank inlet (x=0):
            * **Initial Injection Phase (0 <= t <= t_inject):** During the initial injection period, we may use :
                * Prescribed gas density: ρ_g(0, t) = ρ_in(t).
                * Prescribed flux: -D \* ∇ρ_g(0, t) ⋅ n = q_in(t) (representing the injection rate).
            * **Post-Injection Phase (t > t_inject):** After the injection stops, we will apply a homogeneous Neumann condition:
                * Zero flux: -D \* ∇ρ_g(0, t) ⋅ n = 0, The zero-flux condition at the inlet (z=0) for t > t_inject reflects the scenario where the hydrogen injection is stopped after an initial period, ensuring no further mass transfer across the boundary.
       * At the tank outlet (x=L): Zero flux condition for gas density: -D \* ∇ρ_g(L, t) ⋅ n = 0
    * Initial Conditions:
        * Initial gas density: ρ_g(x, 0) = ρ_g0(x).
        * Initial solid density: ρ_s(x, 0) = ρ_s0(x) (zero for absorption, maximum for desorption, based on Darzi et al.).

* Implementation in Julia:
    * DifferentialEquations.jl: For time integration of the system of ordinary differential equations (ODEs) resulting from spatial discretization.
    * Ferrite.jl: For spatial discretization using the Finite Element Method (FEM).
        * Mesh Generation: We will start with 1D mesh generation using `generate_grid` in Ferrite.jl. For 2D and axisymmetric 3D (which we plan to explore later), we will use Gmsh to generate unstructured meshes and import them into Ferrite.jl.

* Time-Stepping:
    * We will use time-stepping schemes from DifferentialEquations.jl. Implicit methods are likely to be preferred for stability, especially for stiff systems arising from the reaction term.

* Expected Types of Results:
    * Spatiotemporal profiles of hydrogen gas density (ρ_g) along the tank.
    * Temporal evolution of the average hydrogen density in the solid (ρ_s).
    * Qualitative comparison of the absorption/desorption rate with Darzi et al.'s results based on parameters like temperature or pressure.

* Ferrite.jl and DifferentialEquations.jl Tutorials:
    * Ferrite.jl:
        * [Ferrite.jl Documentation and Transient Heat Equation](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/transient_heat_equation/)
    * DifferentialEquations.jl:
        * [DifferentialEquations.jl Documentation](https://docs.sciml.ai/DiffEqDocs/stable/)


**2/ Laminar Flow Modeling of Hydrogen Gas Through the Reactor**

* Physics of the Problem:
    * We will model the flow of hydrogen gas through the porous bed. Initially, we will consider laminar flow. We aim to determine the gas velocity field (`v_g`) for use in the convection-diffusion-reaction model.
    * Darzi et al. use a more comprehensive approach, including momentum conservation, which we will incorporate.
   
    * **Continuity Equation (Conservation of Mass):**
        * From Darzi et al. (Equation 2):
            ```
            ε ∂ρ_g/∂t + ∇ ⋅ (ρ_g v_g) = ṁ
            ```
            where:
            * `ε` is the porosity.
            * `ρ_g` is the gas density.
            * `t` is time.
            * `u_r`, `u_z` are velocity components.
            * `ṁ` is the mass reaction rate.
        
    * **Momentum Equation (Conservation of Momentum):**
        * From Darzi et al. (Equations 6 and 7):
            ```
            ρ_g ∂v_g/∂t + ρ_g (v_g ⋅ ∇) v_g = -∇p + μ ∇²v_g - S
            ```
            where:
            * `p_g` is the gas pressure.
            * `μ` is the viscosity.
            * `S` represents pressure losses due to the porous medium.
     
    * **Energy Equation (Conservation of Energy):**
        * From Darzi et al. (Equation 10):
            ```
            (ρC_p)_eff ∂T/∂t + (ε ρ_g C_{p,g} v_g ⋅ ∇) T = ∇ ⋅ (k_eff ∇T) - ṁ ((1-ε) ΔH + T (C_{p,g} - C_s))
            ```
            where:
            * `T` is the temperature.
            * `(ρC_p)_eff` is the effective volumetric heat capacity.
            * `C_{p,g}` is the heat capacity of hydrogen gas.
            * `ΔH` is the enthalpy of absorption/desorption.
            * `C_s` is the heat capacity of the solid.
            * `k_eff` is the effective thermal conductivity.

* Boundary and Initial Conditions (Adapted from Darzi et al.):
    * ...

* Implementation in Julia:
    * ...

* Time-Stepping:
    * ...

* Expected Types of Results:
    * Pressure, velocity, and temperature fields.
    * Flow rate as a function of pressure drop.
    * Impact of porosity and permeability on flow and heat transfer.
    * Temperature distribution within the reactor.

* Ferrite.jl and DifferentialEquations.jl Tutorials:
    * ...

**3/  Combination of the Two Models (Coupling)**

* Coupling Approach:
    * We will explore coupling the convection-diffusion-reaction model with the flow and energy models.
    * **One-Way Coupling (Initial):**
        1.  Solve the flow model (including energy) to obtain the velocity field (v_g) and temperature field (T).
        2.  Use v_g and T in the convection and reaction terms of the mass and energy conservation equations in the convection-diffusion-reaction model.
        3.  The mass reaction rate (ṁ) depends on ρ_g, ρ_s, and T.
    * **Two-Way Coupling (Advanced):**
        1.  Iteratively solve all models, exchanging information (e.g., pressure and temperature changes from reaction affecting flow and heat transfer).
        2.  Include density and energy changes from the reaction (ṁ and heat release/absorption) in the flow and energy model's mass and energy conservation equations.

* Coupled Iterations or Time Steps:
    * Sequential or iterative coupling within/between time steps.

* Integration of Thermal Effects (Now Integrated):
    * The energy equation is now a fundamental part of the flow model and will be coupled with the convection-diffusion-reaction model.

**Expected Results of the Coupled Model:**

* Spatiotemporal evolution of hydrogen gas (ρ_g) and solid (ρ_s) densities, temperature (T), and velocity (v_g), considering the flow and heat transfer, aiming to qualitatively reproduce trends observed by Darzi et al. regarding the impact of pressure, porosity, and temperature.
* Predictions of storage capacity, cycle times, and temperature distributions.
    
