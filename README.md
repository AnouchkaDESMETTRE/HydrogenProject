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

* Physics of the Problem: ...

* Partial Differential Equations (PDEs) to Solve: ...

* Boundary and Initial Conditions (Adapted from Darzi et al.):
    * Boundary Conditions:
        * At the tank inlet (x=0): We will simulate a gas concentration condition corresponding to the supply pressure (`Pin`) for absorption and the outlet pressure (`Pout`) for desorption, using a simplified relationship (e.g., ideal gas law or an approximation of the sorption isotherm).
        * At the tank outlet (x=L): Zero flux condition for gas concentration.
        * Initial conditions for the solid (zero concentration for absorption, maximum for desorption) based on Darzi et al.'s indications.
    * Initial Conditions: ....

* Implementation in Julia:
   *  DifferentialEquations.jl
   * Spatial Discretization by Finite Element Method (FEM)
   
   * Time-Stepping: ...

* Expected Types of Results:
    * Spatiotemporal profiles of hydrogen gas concentration along the tank.
    * Temporal evolution of the average hydrogen concentration absorbed in the solid.
    * Qualitative comparison of the absorption/desorption rate with Darzi et al.'s results based on parameters like porosity (which will be our initial study parameter).

**2/ Laminar Flow Modeling of Hydrogen Gas Through the Reactor**

* Physics of the Problem:...

* Equations to Solve: ...

* Boundary and Initial Conditions (Adapted from Darzi et al.):
    * Boundary Conditions:
        * At the tank inlet (x=0): Prescribed hydrogen gas pressure (`Pin` or `Pout`).
        * At the tank outlet (x=L): Prescribed hydrogen gas pressure (corresponding to `Pout` or `Pin`).
    * Initial Conditions: Uniform initial pressure field.
      
* Implementation in Julia:
   *  DifferentialEquations.jl
   * Spatial Discretization by Finite Element Method (FEM)
   
   * Time-Stepping: ...

* Expected Types of Results:
    * Pressure profile and gas velocity field along the tank for different inlet/outlet pressures (related to the values studied by Darzi et al.).
    * Analysis of the impact of porosity (via permeability) on the flow rate, qualitatively comparing trends with Darzi et al.'s observations.

**3/ Combination of the Two Models (Coupling)**

* Coupling Approach: 

* Coupled Iterations or Time Steps: 

* Integration of Thermal Effects (Later Stage):

* MOF Adsorption Parameters (Future Extension): 

**Expected Results of the Coupled Model:**

* Spatiotemporal evolution of hydrogen gas and solid concentration/density considering the flow, aiming to qualitatively reproduce the trends observed by Darzi et al. regarding the impact of pressure and porosity.

This adapted README integrates the specific information from Darzi et al.'s article regarding boundary conditions, initial conditions, and key results, directly linking them to the setup and validation of our own model.
