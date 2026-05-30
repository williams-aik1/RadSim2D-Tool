# RadSim2D Tool: Monte Carlo Simulation for CT-Based Radiation Dose Planning Using Patient-Specific Tissue Modeling

![Python](https://img.shields.io/badge/Python-3.10-blue)
![MonteCarlo](https://img.shields.io/badge/Monte%20Carlo-Simulation-orange)
![MedicalPhysics](https://img.shields.io/badge/Medical-Physics-green)
![Streamlit](https://img.shields.io/badge/Streamlit-App-red)

> CT-based physics simulation framework for deterministic and Monte Carlo radiation dose modeling in heterogeneous biological tissue.

---

## Overview

RadSim2D Tool simulates radiation dose deposition in human tissue using real CT (DICOM) images. It compares a deterministic beam model with a Monte Carlo photon transport simulation to demonstrate how stochastic interactions and tissue heterogeneity influence dose distribution.

The system is implemented as an interactive Streamlit web application, enabling real-time CT upload, tissue segmentation, dose visualization, and statistical analysis.

---

## Objectives

- Simulate radiation transport using patient-specific CT images  
- Compare deterministic vs Monte Carlo photon transport models  
- Analyze dose deposition across heterogeneous tissue types  
- Evaluate tumor targeting efficiency and underdose probability  

---

## Features

- 📂 Upload CT DICOM (.dcm) images  
- 🧬 Automatic tissue segmentation (Air / Soft Tissue / Bone)  
- ⚛️ Deterministic vs Monte Carlo dose simulation  
- 🎯 Tumor region detection and beam targeting  
- 🔥 Real-time dose heatmap visualization  
- 📊 Statistical dose analysis dashboard  
- 📉 Underdose probability computation  
- 🌐 Interactive Streamlit interface  

---

## Outputs

### CT Image Visualization

<p align="center">
  <img src="https://github.com/user-attachments/assets/c58b31e8-d1b0-4717-b6dc-81d0d633e773" width="300"/>
  <img src="https://github.com/user-attachments/assets/5f80ab86-71d6-4aa8-96d5-ccbe189f6505" width="400"/>
</p>

---

### Tissue Segmentation Map

<p align="center">
  <img src="https://github.com/user-attachments/assets/2e8cbf3a-d247-45c0-849f-e1455582f900" width="400"/>
</p>

---

### Radiation Dose Heatmap

<p align="center">
  <img src="https://github.com/user-attachments/assets/b95736ca-50ab-4dc8-9a55-a1a7d54a0cc4" width="400"/>
</p>

---

### Tumor Analysis

<p align="center">
  <img src="https://github.com/user-attachments/assets/4a681290-4c46-40d7-8b9b-8cacbf0c203f" width="400"/>
</p>

---

### Dose Distribution Histogram

<p align="center">
  <img src="https://github.com/user-attachments/assets/de68e5c8-2681-4179-b40e-e8f3b95d8d9f" width="450"/>
</p>

---

### Random vs Targeted Beam Comparison

<p align="center">
  <img src="https://github.com/user-attachments/assets/fef18c13-3b01-404a-8ac6-cbe67784d453" width="500"/>
</p>

---

## Technical Approach

### 1. Hounsfield Unit Conversion

CT pixel values are converted into physical density:
