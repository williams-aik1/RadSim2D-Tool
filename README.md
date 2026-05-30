# RadSim2D Tool: Monte Carlo Simulation for CT-Based Radiation Dose Planning Using Patient-Specific Tissue Modeling

![Python](https://img.shields.io/badge/Python-3.10-blue)
![MonteCarlo](https://img.shields.io/badge/Monte%20Carlo-Simulation-orange)
![MedicalPhysics](https://img.shields.io/badge/Medical-Physics-green)
![Streamlit](https://img.shields.io/badge/Streamlit-App-red)

> CT-based physics simulation framework for deterministic and Monte Carlo radiation dose modeling in heterogeneous biological tissue.

---

## Overview

RadSim2D Tool simulates radiation dose deposition in human tissue using real CT (DICOM) images. It compares deterministic and Monte Carlo photon transport models to study how stochastic interactions and tissue heterogeneity affect dose distribution.

The system is implemented as an interactive Streamlit web application enabling:
- CT upload
- Tissue segmentation
- Dose simulation
- Tumor targeting analysis
- Statistical evaluation

---

## Objectives

- Simulate radiation transport using patient-specific CT images  
- Compare deterministic vs Monte Carlo photon transport  
- Analyze dose deposition across heterogeneous tissues  
- Evaluate tumor targeting efficiency and underdose probability  

---

## Features

- Upload CT DICOM (.dcm) images  
- Automatic tissue segmentation (Air / Soft Tissue / Bone)  
- Deterministic vs Monte Carlo dose simulation  
- Tumor region detection and beam targeting  
- Real-time dose heatmap visualization  
- Statistical dose analysis dashboard  
- Underdose probability computation  
- Streamlit interactive interface  

---

## Theory & Methods

### 1. Hounsfield Unit (HU) Conversion

CT pixel values are converted using:

```text
HU = pixel_value × slope + intercept
````

This converts scanner-dependent values into physical density representation.

---

### 2. Tissue Segmentation

Tissues are classified using HU thresholds:

| Tissue Type | HU Range        |
| ----------- | --------------- |
| Air         | HU < -500       |
| Soft Tissue | -500 ≤ HU < 300 |
| Bone        | HU ≥ 300        |

---

### 3. Radiation Attenuation Model

Each tissue has an attenuation coefficient:

| Tissue      | μ    |
| ----------- | ---- |
| Air         | 0.02 |
| Soft Tissue | 0.20 |
| Bone        | 0.50 |

Radiation follows exponential attenuation:

```text
I = I₀ · exp(-μx)
```

---

### 4. Deterministic Model

* Straight beam propagation
* Energy deposition via exponential attenuation
* Minimal scattering approximation

---

### 5. Monte Carlo Model

Photon transport is simulated stochastically:

* Interaction probability:

```text
P(interaction) = 1 - exp(-μ)
```

* Energy deposition occurs at interaction sites
* Random lateral scattering per step
* Energy decay per interaction event

---

### 6. Tumor Targeting

* Tumor defined in soft tissue region
* Gaussian beam centered on tumor
* Increased interaction probability inside tumor
* Comparison with random beam

---

## Results & Insights

* Targeted beam increases tumor dose significantly (~10–40× depending on parameters)
* Monte Carlo model produces realistic heterogeneous dose fields
* Soft tissue absorbs most energy deposition
* Bone shows minimal dose due to high attenuation
* Tumor underdose probability ≈ 0.3–0.5 depending on beam spread

---

## Tech Stack

* Python
* NumPy
* Matplotlib
* Pydicom
* Pandas
* Streamlit

---

## Running the App

### Install dependencies

```bash
pip install -r requirements.txt
```

### Run locally

```bash
streamlit run app.py
```

### Run in Colab

```bash
!streamlit run app.py & npx localtunnel --port 8501
```

---

## Project Structure

```text
RadSim2D/
│── app.py
│── core/
│   │── deterministic.py
│   │── monte_carlo.py
│── requirements.txt
│── README.md
│── sample_data/
```

---

## Future Improvements

* 3D CT volume simulation
* Dose Volume Histogram (DVH)
* Beam angle optimization
* GPU-accelerated Monte Carlo engine
* Energy-dependent photon physics
* Clinical radiotherapy integration

---

## Applications

* Medical physics education
* Monte Carlo research simulation
* Radiation therapy planning concepts
* Biomedical engineering prototyping

---

## Disclaimer

This project is for **educational and research purposes only**.
It is NOT intended for clinical diagnosis or treatment planning.

All outputs are normalized simulation units.

---

## Author

**RadSim2D Tool developed by:**
Williams Stonard Kaphika

---

## Contact

📧 Email: [kaphika.ws@gmail.com](mailto:kaphika.ws@gmail.com)
🔗 LinkedIn: [https://www.linkedin.com/in/williamskaphika](https://www.linkedin.com/in/williamskaphika)
