Codes for optimizing oral contraceptive timing

This model was developed using the Gavina et al. model as the baseline, as described in the manuscript: Optimizing oral contraceptive timing: Morning intake reduces doses and enhances efficacy

To use the model: 
1. Navigate to the folder constant and run the file main_model.m to perform simulations
2. Navigate to the folder src/julia_codes, then run the file run_optimization_clock.sh to perform optimization for dosing times from 1:00 to 24:00.
3. Navigate to the folder constant/Cosine fitting, then run the file cosine_fit.m to fit cosine curve to normalized data
4. Navigate to the folder constant/Parameter estimation drug PK, then run the file PK_est.m to estimate PK parameters
5. Navigate to the folder constant/Parameter estimation circadian, then run the file circ_par_est to estimate circadian parameters of model.
   


--------------------------------------------

Code of Gavina et al. model of hormonal contraception is updated to investigate the role of dosing time in contraceptive efficacy.

The Gavina et al. model (https://github.com/3r3nd/menstrual-cycle-project) was originally developed by modifying the Margolskee et al. model (Margolskee, A. and Selgrade, J.F., 2011. Dynamics and bifurcation of a model for hormonal control of the menstrual cycle with inhibin delay. Mathematical Biosciences, 234(2), pp.95-107) to depict exogenous progesterone’s contraceptive effect.
