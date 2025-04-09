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
    * We will model the evolution of hydrogen gas concentration (ρ_g) within the porous bed and the concentration of hydrogen absorbed in the solid metal-oxide (ρ_s). This initial model will focus on transport and reaction in a transient regime. We will start with a simplified approach, not fully coupling the temperature.
    * This model aims to capture the fundamental behavior of hydrogen absorption and desorption, including diffusion, convection, and reaction kinetics.

* Partial Differential Equations (PDEs) to Solve:
    * **For the density of the hydrogen gas (ρ_g):**
        $$
        ∂(ε * ρ_g)/∂t + ∇ ⋅ (-D_eff * ∇ρ_g - v_g * ρ_g) = S
        $$
        where:
        * `ε` is the porosity of the bed.
        * `t` is time.
        * `D_eff` is the effective diffusion coefficient of hydrogen in the porous bed.
        * `v_g` is the velocity of the gas (initially simplified or from a separate flow model).
        * `S` is a source/sink term representing the absorption/desorption reaction.
    * **For the density of hydrogen in the solid metal-oxide (ρ_s):**
        ```
        ∂((1-ε) * ρ_s)/∂t = -S
        ```
        where:
        * `ρ_s` represents the density of hydrogen in the solid.
        * `-S` is the reaction term, coupling the two equations. A simplified reaction rate law (e.g., first-order) will be used initially.

* Boundary and Initial Conditions (Adapted from Darzi et al.):
    * Boundary Conditions:
        * At the tank inlet (x=0): We need to specify boundary conditions for ρ_g. Since the PDE for ρ_g is second-order, we require two boundary conditions. For a 1D domain (0 <= x <= L), we will start with either:
            * Prescribed gas density: ρ_g(0, t) = ρ_in(t) (This could be related to Pin from Darzi et al. using an equation of state).
            * Prescribed flux: -D_eff \* ∇ρ_g(0, t) ⋅ n = q_in(t) (Where n is the inward normal vector).
            * A combination of both (Robin boundary condition).
            * * **Note:** Pin and Pout are not directly used here, but ρ_in and q_in will be informed by them.
        * At the tank outlet (x=L): Zero flux condition for gas density: -D_eff \* ∇ρ_g(L, t) ⋅ n = 0.
        * ρ_s requires no boundary conditions because the PDE for ρ_s is first-order in time and zero-order in space. This means its evolution is determined locally by the reaction term.
    * Initial Conditions:
        * Initial gas density: ρ_g(x, 0) = ρ_g0(x).
        * Initial solid density: ρ_s(x, 0) = ρ_s0(x) (zero for absorption, maximum for desorption, based on Darzi et al.).

* Implementation in Julia:
    * DifferentialEquations.jl: For time integration of the system of ordinary differential equations (ODEs) resulting from spatial discretization.
    * Ferrite.jl: For spatial discretization using the Finite Element Method (FEM).
        * Mesh Generation: We will start with 1D mesh generation using `generate_grid` in Ferrite.jl. For 2D and axisymmetric 3D (which we plan to explore later), we will use Gmsh or similar meshing tools to generate unstructured meshes and import them into Ferrite.jl.

* Time-Stepping:
    * We will use time-stepping schemes from DifferentialEquations.jl. Implicit methods (e.g., BDF) are likely to be preferred for stability, especially for stiff systems arising from the reaction term.

* Expected Types of Results:
    * Spatiotemporal profiles of hydrogen gas density (ρ_g) along the tank.
    * Temporal evolution of the average hydrogen density in the solid (ρ_s).
    * Qualitative comparison of the absorption/desorption rate with Darzi et al.'s results based on parameters like porosity (ε) (which will be our initial study parameter).

* Ferrite.jl and DifferentialEquations.jl Tutorials:
    * Ferrite.jl:
        * [Ferrite.jl Documentation](https://ferrite-fem.github.io/Ferrite.jl/stable/)
        * [Ferrite.jl Examples](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/)
    * DifferentialEquations.jl:
        * [DifferentialEquations.jl Documentation](https://docs.sciml.ai/DiffEqDocs/stable/)
        * [DifferentialEquations.jl Tutorials](https://docs.sciml.ai/DiffEqTutorials.jl/stable/)

**2/ Laminar Flow Modeling of Hydrogen Gas Through the Reactor**

* Physics of the Problem:
    * We will model the laminar flow of hydrogen gas through the porous bed. Initially, we will use Darcy's law to describe the flow. The goal is to determine the gas velocity field (`v_g`) for use in the convection-diffusion-reaction model.

* Equations to Solve:
    * **Darcy's Law for flow in porous media:**
        ```
        v_g = - (K / μ) * ∇p
        ```
        where:
        * `v_g` is the gas velocity vector.
        * `K` is the permeability of the porous bed (a function of porosity).
        * `μ` is the dynamic viscosity of the hydrogen gas.
        * `∇p` is the pressure gradient.
    * **Continuity equation for the gas (initially steady-state):**
        ```
        ∇ ⋅ (ρ_g * v_g) = 0
        ```
        where:
        * `ρ_g` is the density of the gas.

* Boundary and Initial Conditions (Adapted from Darzi et al.):
    * Boundary Conditions:
        * At the tank inlet (x=0): Prescribed hydrogen gas pressure: p(0) = p_in (corresponding to Pin).
        * At the tank outlet (x=L): Prescribed hydrogen gas pressure: p(L) = p_out (corresponding to Pout).
        * If we extend to 2D/3D, we'll add no-slip boundary conditions on the walls.
    * Initial Conditions: Uniform initial pressure field: p(x, 0) = p_0 (for transient flow calculations, which may be a later extension).

* Implementation in Julia:
    * DifferentialEquations.jl: May be used if we extend to time-dependent flow.
    * Ferrite.jl: For spatial discretization of the pressure field using FEM.
        * Mesh Generation: Similar to Part 1, we'll use `generate_grid` for 1D and Gmsh for 2D/3D, imported into Ferrite.jl.

* Time-Stepping:
    * For steady-state flow, we'll solve the algebraic system of equations. For transient flow (future), we'll use time-stepping methods from DifferentialEquations.jl.

* Expected Types of Results:
    * Pressure profile (p) and gas velocity field (v_g) along the tank for different inlet/outlet pressures (p_in, p_out) related to Darzi et al.'s study.
    * Analysis of the impact of porosity (ε) on the flow rate via its effect on permeability (K), qualitatively comparing trends with Darzi et al.'s observations.

* Ferrite.jl and DifferentialEquations.jl Tutorials:
    * (Same as in Part 1)

**3/ Combination of the Two Models (Coupling)**

* Coupling Approach:
    * We will couple the convection-diffusion-reaction model with the flow model. Initially, we will explore a weak coupling approach:
        1.  Solve the flow model to obtain the gas velocity field (v_g).
        2.  Use this v_g in the convection term of the gas density equation in the convection-diffusion-reaction model.
        3.  The reaction term (S) will depend on both ρ_g and ρ_s.
    * We will investigate how pressure changes due to absorption/desorption affect the flow.

* Coupled Iterations or Time Steps:
    * The coupling will be performed either iteratively within each time step or sequentially between time steps. We will evaluate the stability and accuracy of both approaches.

* Integration of Thermal Effects (Later Stage):
    * We will integrate an energy equation to account for heat release/absorption during reactions and heat transfer. This will involve coupling temperature to the flow and reaction models.

* MOF Adsorption Parameters (Future Extension):
    * If we extend the model to MOFs, we will incorporate adsorption models (e.g., Langmuir) to replace the reaction kinetics used for metal hydrides.

**Expected Results of the Coupled Model:**

* Spatiotemporal evolution of hydrogen gas (ρ_g) and solid (ρ_s) densities, considering the flow, aiming to qualitatively reproduce trends observed by Darzi et al. regarding the impact of pressure and porosity.
* Predictions of storage capacity and cycle times.
* (In future stages) Impact of thermal effects and MOF adsorption parameters.

This adapted README incorporates the specific information from Darzi et al.'s article and addresses the tutor's feedback by:

* Clarifying the use of Pin and Pout and the boundary conditions for ρ_g.
* Explaining why ρ_s needs no boundary conditions.
* Detailing mesh generation plans for 2D and 3D.
* Providing links to relevant tutorials.
