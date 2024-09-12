# A-STEP Framework: Quantifying Heterogeneous Treatment Effects

This repository contains the code for the A-STEP framework, which quantifies heterogeneous treatment effects as published in [XXXXX]. The repository is organized into three primary folders: `Modeling`, `PSM`, and `SHAP`. Each folder contains scripts designed to facilitate model training, reproduction of results, and various analyses presented in the publication.

## Repository Structure

### 1. Modeling
This folder contains the core implementation of the A-STEP framework. Below is a brief overview of the key scripts:

- **`1. Train_from_scratch.R`**  
  Script for users who wish to train the model from scratch. This will involve training the model using the original datasets and parameters.
  
- **`2. Reload_trained_models.R`**  
  Script for replicating the results published in the associated paper. This script utilizes pre-trained hyperparameters stored in the file `trained_models.RData`. You can download the pre-trained model from Zenodo at [https://zenodo.org/records/13736278](https://zenodo.org/records/13736278).
  
- **`3. Interaction_plots.R`**  
  Script used to generate interaction plots as shown in the published paper.
  
- **`4. Plot_PFS.R`**  
  Script used to generate Kaplan-Meier (KM) survival curves, as presented in the paper.
  
- **`PrepareData_SHAP.R`**  
  Script used to preprocess the data for SHAP (SHapley Additive exPlanations) analysis.

### 2. PSM
This folder provides an example of how users can conduct Propensity Score Matching (PSM) to create 1:1 matched pairs, as performed in the published study.

### 3. SHAP
This folder contains the scripts used for conducting SHAP analysis, which was performed to interpret the model results in the paper.

## Data Access
The data utilized in this study can be provided upon reasonable request. Please visit the following Zenodo link to request access: [Zenodo: https://zenodo.org/records/13368111](https://zenodo.org/records/13368111).


---

## Reference
To reference the A-STEP framework, please cite:

```bibtex
@article{AStep,
  title={A-STEP: An Attention-based Scoring for Treatment Effect Prediction in Immunotherapy-Treated Advanced-Stage NSCLC Patients},
  author={},
  journal={},
  year={Year},
  volume={Volume},
  pages={Pages},
  doi={DOI}
}
```
---

