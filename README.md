# RadSim2D Tool: Monte Carlo Simulation for CT-Based Radiation Dose Planning Using Patient-Specific Tissue Modeling

![Python](https://img.shields.io/badge/Python-3.10-blue)
![MonteCarlo](https://img.shields.io/badge/Monte%20Carlo-Simulation-orange)
![MedicalPhysics](https://img.shields.io/badge/Medical-Physics-green)
![Streamlit](https://img.shields.io/badge/Streamlit-App-red)

> CT-based physics simulation framework for deterministic and Monte Carlo radiation dose modeling in heterogeneous biological tissue.

---

## Overview

**RadSim2D Tool** is an interactive CT-based radiation transport simulator designed for educational and research exploration of dose deposition in heterogeneous human tissue.

The system compares:

- Deterministic exponential attenuation beam modeling  
- Monte Carlo photon transport simulation  

using real CT (DICOM) images to demonstrate how tissue heterogeneity and stochastic scattering influence radiation dose distribution.

The framework is implemented as a **Streamlit web application**, enabling real-time visualization, parameter tuning, and comparative analysis.

---

## Objectives

- Simulate radiation transport using patient-specific CT anatomy  
- Compare deterministic vs Monte Carlo dose models  
- Evaluate the effect of tissue heterogeneity on dose deposition  
- Quantify tumor targeting efficiency and dose variability  
- Provide an interactive educational tool for medical physics  

---

## Features

- 📂 Upload CT DICOM (.dcm) images  
- 🧬 Automatic tissue segmentation (Air / Soft Tissue / Bone)  
- ⚛️ Deterministic vs Monte Carlo photon transport simulation  
- 🎯 Tumor region detection and beam targeting  
- 🔥 Real-time dose heatmap visualization  
- 📊 Statistical dose analysis dashboard  
- 📉 Underdose probability estimation  
- 🧪 Interactive parameter control (photons, beam spread, thresholds)  

---

## System Pipeline

RadSim2D follows a modular computational workflow:

1. **Data Acquisition** → Load CT DICOM slice  
2. **HU Conversion** → Convert pixel intensities to Hounsfield Units  
3. **Tissue Segmentation** → Classify anatomical regions  
4. **Physics Mapping** → Assign attenuation coefficients  
5. **Dose Simulation** → Deterministic + Monte Carlo models  
6. **Statistical Analysis** → Tumor dose metrics & uncertainty  

---

## Technical Methodology

### 1. Hounsfield Unit Conversion

CT intensities are converted using scanner calibration:

\[
HU = pixel \cdot slope + intercept
\]

This standardizes CT values into physical density representation.

---

### 2. Tissue Segmentation

| Tissue Type   | HU Range        |
|--------------|-----------------|
| Air          | HU < -500       |
| Soft Tissue  | -500 ≤ HU < 300 |
| Bone         | HU ≥ 300        |

Attenuation coefficients:

| Tissue       | μ (attenuation) |
|--------------|-----------------|
| Air          | 0.02            |
| Soft Tissue  | 0.20            |
| Bone         | 0.50            |

---

### 3. Deterministic Beam Model

Radiation attenuation follows Beer–Lambert law:

\[
I(x) = I_0 e^{-\mu x}
\]

Dose approximation:

\[
D(x) = I_0 (1 - e^{-\mu x})
\]

Characteristics:
- Straight-line propagation  
- No scattering  
- Computationally efficient  
- Smooth dose gradients  

---

### 4. Monte Carlo Photon Transport

Photon propagation is modeled stochastically:

Interaction probability:

\[
P_{interaction} = 1 - e^{-\mu}
\]

Key features:
- Random interaction sampling  
- Energy deposition at collision points  
- Lateral scattering of photons  
- Energy decay per interaction step  

This produces realistic heterogeneous dose fields consistent with stochastic radiation transport behavior.

---

### 5. Tumor Modeling

- Tumor region defined within soft tissue  
- Gaussian beam initialization around tumor axis  
- Targeted photon sampling improves dose deposition  
- Region-based statistical evaluation performed  

Key metrics:
- Mean tumor dose  
- Standard deviation  
- Underdose probability  
- Coverage efficiency  

---

## Key Results

- Targeted beam improves tumor dose by ~40× compared to random beam  
- Monte Carlo simulation produces heterogeneous, realistic dose distributions  
- Soft tissue receives highest energy deposition  
- Bone strongly attenuates photon transport  
- Underdose probability ≈ 0.4 under single-beam conditions  

---

## Visualization Outputs

RadSim2D generates:

- CT slice visualization  
- Hounsfield Unit map  
- Tissue segmentation map  
- Deterministic dose heatmap  
- Monte Carlo dose distribution  
- Tumor overlay visualization  
- Dose histogram analysis  
- Random vs targeted beam comparison  

---

## Tech Stack

- Python  
- NumPy  
- Matplotlib  
- Pydicom  
- Pandas  
- Streamlit  

---

## Running the Application

### Install dependencies

```bash
pip install -r requirements.txt
