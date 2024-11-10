# Metabolic-Thermodynamic modelling of anaerobic oxidation of butane (AOB) and reverse AOB (rAOB)

In scenario 1-3, the electron carrier XH<sub>2</sub>/X between butane-oxidizing archaea and sulfate reducing bacteria was considered to have a redox potential of -0.22 V, equivalent to the average redox potential of sulfate reduction. This setup reflects the proposed nanowire-based extracellular electron transfer mechanism. In scenario 4-5, the electron carrier XH<sub>2</sub>/X was considered to close the redox potential of H<sub>2</sub>/H<sup>+</sup>.


## Scenario 1 (modelling of AOB without energy investment at endergonic steps)

### 1. Model setup
* Oxidation of F<sub>420</sub>H<sub>2</sub> and Fdred are considered as the only sites for energy conservation during the AOB.
* No energy is invested at the step of butyl-CoM/butyryl-CoA conversion and electron transport between X and the MQ pool. 

### 2. Run the model
```
python 13_AOB_M1_SRB.py  
```

### 3. Model results 
* The log10 of the limit concentration of metabolites under scenario 1
	*13_AOB_M1_SRB_no_Einvest__res_log10C_nonfeasible.txt*

* The energy conservation and investment predicted in AOB pathway under scenario 1
	*13_AOB_M1_SRB_no_Einvest__res_dG_nonfeasible.txt*

### 4. Visualization of the results
![Figure 1. Limited concentration of AOB metabolite under scenario 1](https://github.com/SongCanChen11/BackFluxDuringAOB/blob/main/MetabolicModel/13_AOB_github.png)
*Figure 1. The limited concentration of AOB metabolites predicted under scenario 1. Check out our paper for more details (https://www.nature.com/articles/s41467-024-53932-9).*
<br />
<br />
<br />


## Scenario 2 (modelling of AOB with energy investment at endergonic steps)

### 1. Model setup
* Oxidation of F<sub>420</sub>H<sub>2</sub> and Fdred are considered as candidate sites for energy conservation during the AOB.
* In addition, energy is invested at the step of butyl-CoM/butyryl-CoA conversion and electron transport between X and the MQ pool in the form of proton motive force. 
* Energy dissipation is enabled at exergonic reactions in the AOB pathway, such as the step of methylene-tetrahydromethanopterin oxidation.

### 2. Run the model
```
python 17_AOB_M1_SRB.py  
```

### 3. Model results 
* The log10 of the limit concentration of metabolites
	*17_AOB_M1_SRB__res_log10C.txt*

* The energy conservation and investment of the AOB
	*17_AOB_M1_SRB__res_dG.txt*


### 4. Visualization of the results
![Figure 2. Limited concentration of AOB metabolite under scenario 2](https://github.com/SongCanChen11/BackFluxDuringAOB/blob/main/MetabolicModel/17_AOB_github.png)
*Figure 2. The limited concentration of AOB metabolites under scenario 2. Check out our paper for more details (https://www.nature.com/articles/s41467-024-53932-9).*
<br />
<br />
<br />


## Scenario 3 (modelling of rAOB by reversing the sites for energy investment and energy harvest of AOB in Scenario 2)

### 1. Model setup
* Reduction of F<sub>420</sub> and Fd<sub>ox</sub> are considered as sites that need energy investment during the rAOB.
* Energy is harvest at the step of butyl-CoM/butyryl-CoA conversion and electron transport between X and the MQ pool in the form of proton motive force. 
* Energy dissipation is enabled at exergonic reactions in the rAOB pathway, such as the step of acetyl-CoA synthesis.

### 2. Run the model
```
python 18_rAOB_M1_SRB.py  
```

### 3. Model results 
* The log10 of the limit concentration of metabolites
	*18_rAOB_M1_SRB__res_log10C.txt*

* The energy conservation and investment of the AOB
	*18_rAOB_M1_SRB__res_dG.txt*

### 4. Visualization of the results
![Figure 3. Limited concentration of rAOB metabolite under scenario 3](https://github.com/SongCanChen11/BackFluxDuringAOB/blob/main/MetabolicModel/18_rAOB_github.png)
*Figure 3. The limited concentration of rAOB metabolites under scenario 3. Check out our paper for more details (https://www.nature.com/articles/s41467-024-53932-9).*

<br />
<br />
<br />
	

## Scenario 4 (modelling of AOB by assuming the redox potential of XH<sub>2</sub>/X is close to H2/H+)

### 1. Model setup
* Oxidation of F<sub>420</sub>H<sub>2</sub> and Fdred are considered as candidate sites for energy conservation during the AOB.
* In addition, energy is invested at the step of butyl-CoM/butyryl-CoA conversion and electron transport between X and the MQ pool in the form of proton motive force. 


### 2. Run the model
```
python 15_AOB_M1_H2.py
```

### 3. Model results 
* The log10 of the limit concentration of metabolites
	*15_AOB_M1_H2__res_log10C.txt*

* The energy conservation and investment of the AOB
	*15_AOB_M1_H2__res_dG.txt*
  
### 4. Visualization of the results
![Figure 4. Limited concentration of AOB metabolite under scenario 4](https://github.com/SongCanChen11/BackFluxDuringAOB/blob/main/MetabolicModel/15_AOB_log10C_github.png)
*Figure 4. The limited concentration of AOB metabolites under scenario 4. Check out our paper for more details (https://www.nature.com/articles/s41467-024-53932-9).*


<br />
<br />
<br />

## Scenario 5 (modelling of rAOB by assuming the redox potential of XH<sub>2</sub>/X is close to H<sub>2</sub>/H<sup>+</sup>)

### 1. Model setup
* Reduction of F<sub>420</sub> and Fd<sub>ox</sub> are considered as sites that need energy investment during the rAOB.
* Energy is harvest at the step of butyl-CoM/butyryl-CoA conversion and electron transport between X and the MQ pool in the form of proton motive force. 


### 2. Run the model
```
python 16_rAOB_M1_SRB.py  
```

### 3. Model results 
* The log10 of the limit concentration of metabolites
	*16_rAOB_M1_H2__res_log10C.txt*

* The energy conservation and investment of the AOB
	*16_rAOB_M1_H2__res_dG.txt*

### 4. Visualization of the results
![Figure 5. Limited concentration of rAOB metabolite under scenario 5](https://github.com/SongCanChen11/BackFluxDuringAOB/blob/main/MetabolicModel/16_rAOB_log10C_github.png)
*Figure 5. The limited concentration of rAOB metabolites under scenario 5. Check out our paper for more details (https://www.nature.com/articles/s41467-024-53932-9).*

<br />
<br />
<br />

