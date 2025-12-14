## 🧰 How to Use This Template    

Click the green **"Use this template"** button at the top of the page, then choose **"Create a new repository"**.   

This will create your own copy of this project, which you can modify freely — no need to fork!   

 
<p align="center">
  <img src="./images/SHG-banner.png" alt="SHG Logo">
</p>


<h1 align="center">SHG-PW-G-Simulation</h1>

<div align="center">

| **Term** | **Definition** |
|----------|----------------|
| **SHG** | Second Harmonic Generation |
| **PW** | Pulsed Wave |
| **G** | Gaussian |
</div>

&nbsp;

<div align="center">

Article title:       
**The efficiency changes of pulsed Gaussian second harmonic generation in KTP crystal: investigating the influences of pulse energy, laser spot size, cooling temperature by emphasizing the interaction length scale**
</div>

&nbsp;

---

***Table of Contents***

<div>
  &nbsp;&nbsp;&nbsp;&nbsp;<a href="#1-about-this-repository"><i><b>1. About this repository</b></i></a>
</div>
&nbsp;

<div>
  &nbsp;&nbsp;&nbsp;&nbsp;<a href="#2-getting-started"><i><b>2. Getting Started</b></i></a>
</div>
&nbsp;

<div>
  &nbsp;&nbsp;&nbsp;&nbsp;<a href="#3-how-to-cite-us"><i><b>3. How to Cite Us</b></i></a>
</div>
&nbsp;


<div>
  &nbsp;&nbsp;&nbsp;&nbsp;<a href="#4-contact-information"><i><b>4. Contact Information</b></i></a>
</div>
&nbsp;

---    

## 1. About this repository

This repository contains the **Computational Toolkit for Pulsed Gaussian Second Harmonic Generation in KTP Crystal**, an open-source Fortran implementation developed to investigate the efficiency changes of type-II second harmonic generation (SHG) under repetitively pulsed Gaussian laser pumping. This toolkit addresses the research problem described in the article: **"The efficiency changes of pulsed Gaussian second harmonic generation in KTP crystal: investigating the influences of pulse energy, laser spot size, cooling temperature by emphasizing the interaction length scale"**

This toolkit implements a coupled field-phase-heat equation model that simultaneously solves the thermal effects in type II second harmonic generation of Gaussian pulsed waves in a KTP crystal. The model incorporates thermally induced phase mismatch (TIPM) through the thermal dispersion relations of the ordinary and extraordinary refractive indices, providing a comprehensive framework for analyzing how pulse energy, laser spot size, and cooling temperature influence SHG efficiency with particular emphasis on the interaction length scale.

The toolkit provides:
- **Coupled equation solver** for simultaneous solution of SHG field, phase, and heat equations using the Finite Difference Method (FDM)
- **Pulsed Gaussian beam simulation** with realistic pulse characteristics (pulse duration, repetition frequency, pulse energy)
- **Thermal effects modeling** including spatiotemporal temperature distribution and thermally induced phase mismatch
- **Interaction length analysis** with temperature-dependent variations along radial and longitudinal directions
- **Time evolution analysis** from transient to steady-state conditions across multiple pulses
- **Gaussian beam propagation** with absorption coefficients for fundamental and second-harmonic waves
- **KTP crystal properties** with temperature-dependent material parameters and realistic cooling mechanisms (conduction, convection, radiation)
- **Home-computer compatible** numerical procedures that enable efficient computation on standard personal computers

The implementation has been validated through comprehensive numerical studies, demonstrating how SHG efficiency varies with pulse energy, reaching local maxima near 70% efficiency. The model successfully captures the cyclical behavior where higher pulse energy initially enhances SHG efficiency, but the associated temperature rise induces phase mismatch, reducing efficiency and driving partial reconversion of the second-harmonic wave to the fundamental wave. This toolkit provides a complete computational framework for analyzing thermal effects and optimizing SHG efficiency in pulsed Gaussian laser systems.



```
Folder PATH listing
+---citation                       <-- Contains research paper PDFs
│       1_Heat-Equation_Continuou… <-- Analytical heat equation paper
│       2_Heat-Equation_Continuou… <-- Heat equation paper
│       3_Heat-Equation_Pulsed-Wa… <-- Pulsed wave heat equation
│       4_Phase-Mismatch_Pulsed-W… <-- Phase mismatch paper
│       5_Ideal_Continuous-Wave_G… <-- Ideal continuous wave paper
│       6_Ideal_Pulsed-Wave_Besse… <-- Ideal pulsed wave paper
│       7_Coupled_Continuous-Wave… <-- Coupled continuous wave paper
│       README.md                  <-- Citation guidelines
│
+---images                         <-- Contains project images
│       SHG-banner.png             <-- Project banner image
│
+---results                        <-- Contains simulation output files
│       E045_f4000_Np1_tp50_Elec1… <-- Electric field 12r data
│       E045_f4000_Np1_tp50_Elec1… <-- Electric field 12t data
│       E045_f4000_Np1_tp50_Elec1… <-- Electric field 12z data
│       E045_f4000_Np1_tp50_Elec2… <-- Electric field 22r data
│       E045_f4000_Np1_tp50_Elec2… <-- Electric field 22t data
│       E045_f4000_Np1_tp50_Elec2… <-- Electric field 22z data
│       E045_f4000_Np1_tp50_Elec3… <-- Electric field 32r data
│       E045_f4000_Np1_tp50_Elec3… <-- Electric field 32t data
│       E045_f4000_Np1_tp50_Elec3… <-- Electric field 32z data
│       E045_f4000_Np1_tp50_ibest… <-- Best current data
│       E045_f4000_Np1_tp50_Phase… <-- Phase minimum data
│       E045_f4000_Np1_tp50_Pr.pl… <-- Power radial component
│       E045_f4000_Np1_tp50_Psi2p… <-- Psi2 picks data
│       E045_f4000_Np1_tp50_Psi3p… <-- Psi3 picks data
│       E045_f4000_Np1_tp50_Pt.pl… <-- Power theta component
│       E045_f4000_Np1_tp50_Pz.pl… <-- Power z component
│       E045_f4000_Np1_tp50_Tempm… <-- Maximum temperature data
│       E045_f4000_Np1_tp50_Tr.pl… <-- Temperature radial component
│       E045_f4000_Np1_tp50_Tt.pl… <-- Temperature theta component
│       E045_f4000_Np1_tp50_Tz.pl… <-- Temperature z component
│
+---src                            <-- Contains source code
│       Code_SHG-PW-G-Simulation…  <-- Main Fortran simulation code
│
        LICENSE                    <-- License information
        README.md                  <-- Project documentation
```

## 2. Getting Started

### 2.1. Prerequisites
- **Intel Fortran Compiler** (ifort) or compatible Fortran compiler
- **Text Editor** (VS Code, Cursor, or any Fortran-capable editor)
- **PDF Reader** (for accessing research papers and documentation)
- **Git** (for cloning the repository)
- **Make** (for building the project, optional but recommended)

### 2.2. Quick Start

1. **Clone the repository**
   ```bash
   git clone https://github.com/Second-Harmonic-Generation/SHG-PW-G-Simulation.git
   cd SHG-PW-G-Simulation
   ```

2. **Explore the Research Papers**
   - Review the `citation/` folder for the main research article and supporting references
   - Read the `README.md` files in each subdirectory for detailed explanations
   - Consult the article for detailed theoretical background and validation results

3. **Compile and Run the Code**
   ```bash
   cd src/
   ifort -o shg_simulation Code_SHG-PW-G-Simulation.f90
   ./shg_simulation
   ```

4. **Analyze Results**
   - Check the `results/` folder for generated plot data files (.plt format)
   - Use your preferred plotting software to visualize the results
   - Analyze interaction length variations, temperature distributions, and efficiency changes
   - Compare with the theoretical predictions and findings reported in the research article

5. **Development Workflow**
   - Edit the Fortran source code in `src/Code_SHG-PW-G-Simulation.f90`
   - Modify parameters such as pulse energy, spot size, and cooling temperature as needed for your specific analysis
   - Recompile and run to generate new results
   - Document your findings and modifications


## 3. How to Cite Us
Please refer to the [**citation**](./citation/) folder for accurate citations. It contains essential guidelines for accurate referencing, ensuring accurate acknowledgement of our work.


  
## 4. Contact Information

For questions not addressed in the resources above, please connect with [Mostafa Rezaee](https://www.linkedin.com/in/mostafa-rezaee/) on LinkedIn for personalized assistance.
