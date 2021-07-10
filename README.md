# drugkinetics-resistance
## Description
This repository contains code that was used to generate the results in 'The risk of drug resistance during long-acting antimicrobial therapy' by Anjalika Nande and Alison L Hill.
## Contents
* sinusoidal_model_preexistence_rescue.nb : This Mathematica notebook contains code that was used to calculate the risks of resistance due to pre-existing and rescue mutants under the sinusoidal model of drug efficacy.
* pharma_model_preexistence_rescue.nb : This Mathematica notebook contains code that was used to calculate the risks of resistance due to pre-existing and rescue mutants under the pharmacologically-inspired model of drug efficacy.
* pharma_latency.nb : The code used to calculate the risks of resistance due to latency reactivation events is in this Mathematica notebook. This uses the pharmacologically-inspired model of drug efficacy.
* pharma_adherence.nb : This Mathematica notebook contains the code that was used to calculate the risks of resistance due to pre-existing and rescue mutants in presence of imperfect adherence under the pharmacologically-inspired model of drug efficacy. 
* pharma_model_twostep_mutant.nb :  This risk of resistance due to a two-step rescue mutant under the pharmacologically-inspired model of drug efficacy was calculated using the code in this notebook.
* stochastic_preexistence.py : This file contains code in Python 3.0 that was used to stochastically simulate the full model in order to get the establishment probability of one pre-existing resistant mutant.
* stochastic_rescue.py : This file contains code in Python 3.0 that was used to stochastically simulate the full model in order to obtain the risks of resistance associated with rescue mutants.
* stochastic_latency.py :  This file contains code in Python 3.0 that was used to stochastically simulate the full model in order to obtain the risks of resistance associated with events of latency reactivation. 
