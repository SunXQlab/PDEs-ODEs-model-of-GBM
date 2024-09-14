# MMOTMIUIRIG
Multiscale Modeling of Tumor-Macrophage Interactions Underlying Immunotherapy Resistance<br> 
We designed a hybrid analytical-numerical approach for simulating the coupled PDE-ODE model. The **MMOTMIUIRIG** model consists of four main functions:<br> 
1.**'sum1.m'**: This function simulates tumor growth, drug treatment, and the development of drug resistance by solving a mathematical model of tumors.<br> 
2. **'ParaEstimate_1_ODE.m'**: This function estimates parameters using an improved genetic algorithm.<br> 
3. **'LHS_PRCC2.m'**: This function performs sensitivity analysis of parameters.<br> 
4.  **'Stochastic_Simulation_2.m'**: This function runs a random simulation for an in-silico cohort.<br> 
The input and output of each function are described below.<br> 


# 1.sum1

By solving a mathematical model of tumors, we can simulate the processes of tumor growth, drug treatment, and the development of drug resistance.<br> 
**select_CSF1R_I, select_EGFR_I, select_IGF1R_I**: These parameters represent the delivery methods for CSF1R inhibitor, EGFR inhibitor, and IGF1R inhibitor, respectively. The values 0, 1, 2, 3, 4, and 5 correspond to the following boundary dosing methods:<br> 

-   **0**: No boundary dosing<br> 
-   **1**: Constant boundary dosing<br> 
-   **2**: Sinusoidal boundary dosing<br> 
-   **3**: Periodic binary boundary dosing<br> 
-   **4**: Binary boundary dosing<br> 
-   **5**: Adaptive boundary dosing<br> 

This function is used to solve a dynamic model and simulate the processes of tumor growth, drug treatment, and development of drug resistance.<br> 
## Input

**Parameters_nondimensional_5_gai_jin_fangcheng_5**: The parameter values of the model<br>
**pdex4pde**: The format of partial differential equations (PDEs)<br> 
**odex4ode**: The format of ordinary differential equations (ODEs)<br> 
**pdex4ic**: The initial conditions for partial differential equations  (PDEs) and ordinary differential equations (ODEs) refer to the initial values needed to solve these equations<br> 
**pdex4bc**:The boundary conditions for partial differential equations  (PDEs)<br> 

## Output

Model Solving is:<br> 
**sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options)** or<br> 
**sol = pdepe3(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options)**<br> 
sol: The values of variables changing over time in space in partial differential equations (PDEs) and ordinary differential equations (ODEs).<br> 
**u1** = sol(:,:,1): Tumor cell density<br> 
**u2** = sol(:,:,2): M1 macrophage density<br> 
**u3** = sol(:,:,3): M2 macrophage density<br> 
**u4** = sol(:,:,4): CSF1 Concentration<br> 
**u5** = sol(:,:,5): EGF Concentration<br> 
**u6** = sol(:,:,6): IGF1 Concentration<br> 
**u7** = sol(:,:,7): EGFR Concentration<br> 
**u8** = sol(:,:,8): ERK Concentration<br> 
**u9** = sol(:,:,9): IGF1R Concentration<br> 
**u10** = sol(:,:,10): AKT Concentration<br> 
**C_T_mean**=sum(u1(:,:),2)/100: Average cancer cell density<br> 
**M1_mean**=sum(u2(:,:),2)/100: Average M1 macrophage density<br> 
**M2_mean**=sum(u3(:,:),2)/100: Average M2 macrophage density<br> 
Figures illustrating the density or concentration of 13 variables changing over time.

# 2.ParaEstimate_1_ODE

This function estimates parameters using an improved genetic algorithm.<br> 

## Input
**group** = 10: Using a genetic algorithm to estimate 10 sets of parameters<br> 
**NVAR**  = 22: Number of variables, according to the original setting of Para0<br> 
**shu_min** = 0.000000001: First-time minimum search<br>  
**shu_max** = 50: First-time maximum search<br> 
**FieldD** = [rep(PRECI,[1, NVAR]); [shu_min*ones(1,NVAR);shu_max*ones(1,NVAR)] rep([1; 0; 1 ;1], [1, NVAR])]: Decoding function<br> 
**Objfun1**: Fitness Evaluation Function. Evaluating fitness by calculating a fitness score for each individual (chromosome).<br> 
**Conditions_ODES_2**: Solvingg ODE (Ordinary Differential Equation)<br> 
**Conditions**: 11 conditions<br> 
**ExpData**: The experimental data corresponding to the 11 conditions<br> 
**intraODEsys2**: The format of ordinary differential equations<br> 

## Output
**chrom_value** = bs2rv(Chrom,FieldD):  The estimated values of the parameters<br> 
Four comparison figures of experimental data and model data.<br> 

# 3.LHS_PRCC2

 This function performs sensitivity analysis of parameters.

## Input
**Parameters_nondimensional_6_gai_jin_fangcheng_6**<br> 
**select_CSF1R_I**: 0 corresponds to no boundary dosing, and 1 corresponds to constant boundary dosing with CSF1R inhibitors.<br> 
**lb**: lower bounds for parameters<br> 
**ub**: upper bounds for parameters<br> 

## Output
**C3**: Individual's survival time<br> 
**CC:** For each parameter, the Pearson correlation coefficient between the 1000 sets of sample values and simulated survival times<br> 
**Pvalue**: Pearson correlation coefficient p-values<br> 
Figures describing Pearson correlation coefficients and Pearson correlation coefficient p-values<br> 

# 4.Stochastic_Simulation_2

This is a program that performs stochastic simulation of two key parameters in the model.<br> 

## Input
**Parameters_nondimensional_6_gai_jin_fangcheng_6**
**select_CSF1R_I**: 0 corresponds to no boundary dosing, and 1 corresponds to constant boundary dosing with CSF1R inhibitors.<br> 
**p**: Population size<br> 
**mu**: Means of normal distributions or bimodal distributions<br> 
**sigma**: Variances  of normal distributions or bimodal distributions<br> 
**rrr**: Weights or proportions of each peak<br> 

## Output
**C3**: Individual's survival time
**C5**: The number of surviving individuals over a period of time
Figures of stochastic simulation based on an in-silico cohort and comparison to survival data.

