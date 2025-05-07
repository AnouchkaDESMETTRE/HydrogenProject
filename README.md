This README documents our review of the article "Absorption and desorption of hydrogen in long metal hydride tank equipped with phase change material jacket" (Rabienataj Darzi et al., International Journal of Hydrogen Energy, April 2016) and outlines our project plan to model similar aspects of hydrogen storage.

#  Review of "Absorption and desorption of hydrogen in long metal hydride tank equipped with phase change material jacket" 

The article titled "**Absorption and desorption of hydrogen in long metal hydride tank equipped with phase change material jacket**" was published in the International Journal of Hydrogen Energy in **April 2016** and was authored by **A. Ali Rabienataj Darzi, H. Hassanzadeh Afrouzi, A. Moshfegh, and M. Farhadi**. We are reviewing this study to understand its findings for application in our own research or project.

This study numerically investigated the **hydrogen absorption and desorption processes** in a long tubular **LaNi5 metal hydride tank (MHT)** integrated with a **Rubitherm phase change material (PCM) jacket** for PEM fuel cell hydrogen supply. The authors analyzed the influence of different **H2 supply pressures (10, 15, and 20 bar)**, different **discharge pressures (1.5, 1.75, and 2 bar)**, and **metal hydride bed porosities (0.4, 0.5, and 0.6)** on transient and local temperature distributions across the H2-MHT system and PCM jacket. The time-dependent changes of the **hydrogen to metal (H/M) ratio** and **PCM melt fraction** were also investigated until equilibrium was reached.

Key findings of the Darzi's Thesis include:

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

* Coupled Densities Equations:
    * **For the density of the hydrogen gas (ρ_g):**
        ```
        ε ∂ρ_g/∂t = D (∂²ρ_g/∂x² + ∂²ρ_g/∂z²) + u_x ∂ρ_g/∂x + u_z ∂ρ_g/∂z + ṁ(ρ_s, t)
        ```
  * **For the density of hydrogen in the solid metal-oxide `(ρ_s)`:**
        ```
        (1 - ε) ∂ρ_s/∂t = ṁ(ρ_s, t)
        ```
        Where:
        * `ρ_g` is the density of the gas (hydrogen). Its value varies as a function of time (`t`) and space (`x`, `z`) as the equation describes its variation. The gas is assumed to behave as an ideal gas.
        * `ρ_s` represents the density of hydrogen in the solid,  it also varies as a function of time (`t`) and space (`x`, `z`).
        * `ε` is the porosity of the bed. This structural property significantly impacts absorption and desorption times, as well as heat transfer. The values used in this study include `0.4`, `0.5`, and `0.6`.
        * `t` is time.
        * `u_x` and `u_z` are the velocities of hydrogen in the longitudinal and radial directions, respectively. These velocities are variables that are obtained by solving the momentum equations.
        * `ṁ(ρ_s, t)` is the mass reaction rate of hydrogen. It represents the mass of hydrogen absorbed or desorbed per unit volume and time. This rate is calculated using separate kinetic equations: 
        
    * **For Absorption:**
        ```
        ṁ = C_a exp(-E_a / (RT)) ln(p_g / p_{eq,a}) (ρ_{sat} - ρ_s)
        ```

    * **For Desorption:**
        ```
        ṁ = C_d exp(-E_d / (RT)) (p_g - p_{eq,d}) / p_{eq,d} (ρ_s - ρ_{emp})
        ```
        
        Where:
      
        * `C_a` is the absorption rate coefficient (Value: 59.187 s⁻¹).
      
        * `C_d` is the desorption rate coefficient (Value: 9.57 s⁻¹).
      
        * `E_a` is the activation energy for absorption (Value: 21179.6 J/mol).
      
        * `E_d` is the activation energy for desorption (Value: 16473 J/mol).
      
        * `R` is the universal gas constant (Value: 8.314 J/mole K).
      
        * `T` is Temperature (Variable obtained by solving the heat transfer equation - Equation 10 in Darzie's Thesis).
   
        * `p_g` is Absolute pressure (Variable, see momentum equations).
      
        * `p_{eq,a}` is the equilibrium pressure for absorption and `p_{eq,d}`for desorption, calculated using the Van't Hoff relation: `p_{eq} = p_{ref} exp(A - B/T)` with `A = 10.7`, `B = 3704.6`, and `p_{ref} = 1 MPa` for `p_{eq,a}` and `A = 10.57`, `B = 3704.6`, and `p_{ref} = 1 MPa` for `p_{eq,d}`.
            
        * `ρ_{sat}` is the saturated bed density (Value: 7259 kg/m³).
      
        * `ρ_{emp}` (r_emp) is the density of the metal hydride without hydrogen (empty bed density) (Value: 7164 kg/m³).
      
        * `ρ_s` is the solid bed density (Variable, its evolution is given by an additional equation).
      
But a simplified reaction rate law (e.g., first-order) will be used initially.

      
* Initial and Boundary Conditions:
    * Initial Conditions :
         * Initial Solid Density `ρ_s(t=0)`:
             * Absorption: Initial solid density is equal to the density of the metal hydride bed without hydrogen `ρ_{emp}`
             * Desorption: Initial solid density is equal to the saturated density of the metal hydride bed  `ρ_{sat}`.
         * Initial Temperature `T(t=0)`:
             * Absorption: Initial temperature of the entire system (PCM and MHT tank) is set to 301.15 K.
             * Desorption: Initial temperature of the entire system is set to 305.15 K.
         * Initial Hydrogen Pressure `P(t=0)`: The initial hydrogen pressure is assumed to be equal to the equilibrium pressure at the initial temperature of the system. This equilibrium pressure is calculated using the Van't Hoff equation:
    ```
    Ln(P_{eq}/P_{ref}) = A - B/T
    ```
Where:
* `P_{eq}`: Equilibrium pressure of the MHT.
* `P_{ref}`: Reference pressure (1 MPa).
* `T`: Equilibrium temperature between hydrogen and the MHT.
* `A, B`: Constants specific to the absorption: A = 10.7 et B = 3704.6 - desorption: A = 10.57 et B = 3704.6 process.

For the initial condition of the simulation, the initial hydrogen pressure is assumed to be the equilibrium pressure at the initial system temperature `T (t=0)`. We use the Van't Hoff equation with `T = T_{0}` to calculate this initial pressure `P_{initial} = P_{eq}(T_{initial})`.

Knowing the initial pressure and initial temperature, the initial density of the gaseous hydrogen `ρ_g` can then be determined using the ideal gas law :  `ρ_{g, initial} = P_{initial}/(R_{specific}xT_{initial}`, with R the specific constante of hydrogen gaz : R=4124$ J/mole K.

   * Boundary Conditions:
    
While explicit boundary conditions for the gas phase density `ρ_g` and the solid phase density `ρ_s` are not directly listed, their boundary behavior is inherently linked to those imposed on pressure `P`, temperature `T`, and velocity components `u_x, u_z` used in solving the conservation of mass and momentum equations; however, initial boundary conditions will still be imposed to solve these initial density equations.

Initially, we will assume that `ρ_s` is spatially independent. This is a simplification, as `ρ_s` actually depends on the mass reaction rate `ṁ`, which in turn is a function of temperature and pressure – both of which vary spatially. However, our initial focus will be solely on the first two density equations, employing a simplified form of the mass reaction rate. Consequently, we choose to impose a non-homogeneous Dirichlet boundary condition on `ρ_g` at the left inlet boundary for hydrogen. This absorption case is the first scenario we are testing.

* Implementation in Julia:
    * **Spatial Discretization (Finite Element Method - FEM)**

The spatial discretization of the domain is achieved using the Ferrite.jl library and the Finite Element Method (FEM). The main steps involved are:

1.  *Mesh Generation:* We generate a simple rectangular mesh for testing purposes using the `generate_grid` function in Ferrite.jl. We create a structured quadrilateral mesh over an arbitrary rectangular domain, which serves as the basis for the finite element analysis. For more complex 2D or axisymmetric 3D geometries (intended for later exploration), external tools like Gmsh are used to generate unstructured meshes, which are then imported into Ferrite.jl.

2.  *Definition of Finite Element Spaces:* For each variable in the problem (`ρ_s` represented by `:u1` and `ρ_g` represented by `:u2`), a finite element space is defined using Lagrange elements (`Lagrange{RefQuadrilateral, degree}()`). The `degree` parameter specifies the polynomial order of the basis functions used for interpolation within each element.

3.  *Quadrature Rules:* For the numerical integration of the weak forms of the equations over each element, quadrature rules (`QuadratureRule{RefQuadrilateral}(2*degree+1)`) are defined. The order of the quadrature rule is chosen based on the degree of the polynomials to be integrated to ensure exact or sufficiently accurate integration.

4.  *CellValues:* The `CellValues` objects (`cellvalues_u1` and `cellvalues_u2`) are created for each finite element space and quadrature rule. These objects efficiently compute the values of the basis functions and their gradients at the quadrature points within each cell of the mesh.

5.  *DofHandler (Degrees of Freedom Management):* The `DofHandler` is responsible for associating the degrees of freedom (DOFs) with the nodes of the mesh for each variable. It assigns a unique global index to each DOF, enabling the construction of global solution vectors.

6.  *Assembly of Matrices:* The functions `assemble_mass_matrix` and `assemble_stiffness_matrix` iterate over each cell of the mesh and contribute to the construction of the global mass matrix (`M`) and stiffness matrix (`K`). They utilize the values of the basis functions and their gradients as well as the determinant of the Jacobian (`getdetJdV`), to perform numerical integration over each element.

7.  *Boundary Conditions:* Dirichlet boundary conditions are applied using a `ConstraintHandler` (`ch`). Facet boundaries where the conditions are imposed are defined, and the boundary values are specified (here for `u2` on the "left" boundary).

    * **Numerical Time Integration**
    
The numerical time integration of the system of ordinary differential equations (ODEs) resulting from the spatial discretization is performed using the DifferentialEquations.jl library. The main steps are:

1.  *Definition of the Right-Hand Side (RHS) Function:* The function `rhs!(du, u, p, t)` defines the system of ODEs. It takes as input the time derivatives vector (`du`), the current solution vector (`u`), the problem parameters (`p`), and the current time (`t`). Within this function, boundary conditions are applied, and the time derivatives are computed based on the current solution and the spatial matrices (`K`). The coupling between the variables `u1` and `u2` is also implemented directly within this function with `ṁ`.

2.  *Definition of Initial Conditions:* A vector `uinit` is created and initialized with the initial conditions for all variables across the entire mesh using the `setup_initial_conditions!` function.

3.  *Definition of the Time Span:* The time interval of the simulation (`tspan`) and the final time (`Tend`) are specified.

4.  *Definition of the ODE Problem:* an `ODEProblem` object is created by providing the RHS function (`rhs`), the initial conditions (`uinit`), the time span (`tspan`), and the problem parameters (`p`).

5.  *Choice of Time Stepper:* A time integration solver is chosen to integrate the system of ODEs. In the code, `Rodas5P` is used, which is an implicit solver suitable for stiff problems.

6.  *Solving the Problem:* The `solve` function (or initializing an integrator with `init` and iterating via `intervals`) is used to perform the time integration and obtain the solution `sol` (or the integrator `integrator`).

* Expected Types of Results:
The solution is then processed and visualized, by writing the results to VTK files for visualization with Paraview.

    * Spatiotemporal profiles of hydrogen gas density (ρ_g) along the tank.
    * Temporal evolution of the average hydrogen density in the solid (ρ_s).
    * Qualitative comparison of the absorption/desorption rate with Darzi et al.'s results based on parameters like temperature or pressure.

* Ferrite.jl and DifferentialEquations.jl Tutorials:
    * Ferrite.jl:
        * [Ferrite.jl Documentation and Transient Heat Equation](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/transient_heat_equation/)
    * DifferentialEquations.jl:
        * [DifferentialEquations.jl Mathematical Specification of an ODE Problem](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction)


**2/ Laminar Flow Modeling of Hydrogen Gas Through the Reactor**

* Physics of the Problem:
    * We will model the flow of hydrogen gas through the porous bed. Initially, we will consider laminar flow. We aim to determine the gas velocity field (`v_g`) for use in the convection-diffusion-reaction model.
    * Darzi et al. use a more comprehensive approach, including momentum conservation, which we will incorporate.
   
   * **Continuity Equation (Conservation of Mass):**
       ```
        div(u) = 0
       ```
       where `u(z)` is the unknown velocity. In 1D, this implies a spatially constant velocity.

   * **Momentum Equation (Conservation of Momentum):**
       ```
       u(z) ∂u/∂z = 1/Re ∂^2u/∂z^2 + f(z)
       ```
       where:
       * `u(z)` is the unknown velocity.
       * `f(z)` is the pressure gradient.
       * `Re > 0` is the Reynolds number.
     
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
    
