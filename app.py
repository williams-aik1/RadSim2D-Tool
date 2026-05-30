import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pydicom

# ==========================================================
# CONFIG
# ==========================================================

st.set_page_config(
    page_title="RadSim2D Tool",
    layout="wide"
)

AUTHOR = "Williams Stonard Kaphika"

# ==========================================================
# HEADER
# ==========================================================

st.title("RadSim2D Tool")
st.subheader("CT-Based Monte Carlo Radiation Dose Simulation Framework")

st.markdown(
    f"""
**Developer:** {AUTHOR}

This framework demonstrates deterministic and Monte Carlo
radiation transport modelling using CT-derived tissue
information.
"""
)

st.warning(
    """
DISCLAIMER

This software is intended solely for education, research,
and demonstration purposes.

It is NOT clinically validated and must NOT be used for:

* patient treatment planning
* clinical decision making
* medical diagnosis
* radiation prescription

All dose values are relative simulation units.
"""
)

# ==========================================================
# SIDEBAR
# ==========================================================

st.sidebar.header("Simulation Controls")

n_photons = st.sidebar.slider(
    "Number of Photons",
    1000,
    50000,
    15000,
    step=1000
)

beam_sigma = st.sidebar.slider(
    "Beam Spread (σ)",
    1,
    15,
    3
)

dose_threshold = st.sidebar.slider(
    "Underdose Threshold",
    0.05,
    1.0,
    0.30
)

uploaded_file = st.file_uploader(
    "Upload CT DICOM",
    type=["dcm"]
)

# ==========================================================
# MONTE CARLO
# ==========================================================

@st.cache_data(show_spinner=True)
def monte_carlo(
    mu_map,
    tumor_mask,
    center_x,
    n_photons,
    sigma
):
    rows, cols = mu_map.shape
    dose = np.zeros((rows, cols))

    for _ in range(n_photons):

        x = int(np.random.normal(center_x, sigma))
        x = np.clip(x, 0, cols - 1)

        y = 0
        energy = 1.0

        while energy > 0.001 and y < rows:

            mu = mu_map[y, x]
            p_interact = 1 - np.exp(-mu)

            if np.random.random() < p_interact:

                deposit = energy * (
                    0.20 if tumor_mask[y, x] else 0.10
                )

                dose[y, x] += deposit
                energy *= 0.90

            x += np.random.choice(
                [-1, 0, 1],
                p=[0.1, 0.8, 0.1]
            )

            x = np.clip(x, 0, cols - 1)
            y += 1

    return dose


# ==========================================================
# DETERMINISTIC
# ==========================================================

def deterministic(mu_map):

    rows, cols = mu_map.shape
    dose = np.zeros_like(mu_map)

    I0 = 100

    for c in range(cols):

        intensity = I0

        for r in range(rows):

            mu = mu_map[r, c]

            attenuation = np.exp(-mu)

            deposited = intensity * (1 - attenuation)

            dose[r, c] += deposited

            intensity *= attenuation

            if intensity < 1e-6:
                break

    return dose


# ==========================================================
# MAIN
# ==========================================================

if uploaded_file is None:
    st.info("Upload a CT DICOM file to begin.")
    st.stop()

# ==========================================================
# LOAD DICOM
# ==========================================================

ds = pydicom.dcmread(uploaded_file)

image = ds.pixel_array.astype(float)
image = image[::2, ::2]

slope = getattr(ds, "RescaleSlope", 1.0)
intercept = getattr(ds, "RescaleIntercept", 0.0)

hu = image * slope + intercept

# ==========================================================
# SEGMENTATION
# ==========================================================

tissue = np.zeros_like(hu)

tissue[hu < -500] = 0
tissue[(hu >= -500) & (hu < 300)] = 1
tissue[hu >= 300] = 2

mu_map = np.zeros_like(hu)

mu_map[tissue == 0] = 0.02
mu_map[tissue == 1] = 0.20
mu_map[tissue == 2] = 0.50

rows, cols = tissue.shape

# ==========================================================
# TUMOR REGION
# ==========================================================

soft = np.argwhere(tissue == 1)

if len(soft) == 0:
    st.error("No soft tissue detected.")
    st.stop()

r_center, c_center = soft[len(soft) // 2]

tumor_mask = np.zeros_like(tissue, dtype=bool)

r1 = max(0, r_center - 10)
r2 = min(rows, r_center + 10)

c1 = max(0, c_center - 10)
c2 = min(cols, c_center + 10)

tumor_mask[r1:r2, c1:c2] = True

# ==========================================================
# RUN SIMULATIONS
# ==========================================================

dose_det = deterministic(mu_map)

dose_mc = monte_carlo(
    mu_map,
    tumor_mask,
    c_center,
    n_photons,
    beam_sigma
)

dose_norm = dose_mc / (dose_mc.max() + 1e-8)

tumor_dose = dose_norm[tumor_mask]

# ==========================================================
# METRICS
# ==========================================================

mean_tumor = np.mean(tumor_dose)
std_tumor = np.std(tumor_dose)

underdose_prob = np.mean(
    tumor_dose < dose_threshold
)

coverage_95 = np.mean(
    tumor_dose >= 0.95
)

# ==========================================================
# DASHBOARD
# ==========================================================

m1, m2, m3, m4 = st.columns(4)

m1.metric(
    "Mean Tumor Dose",
    f"{mean_tumor:.3f}"
)

m2.metric(
    "Dose Std",
    f"{std_tumor:.3f}"
)

m3.metric(
    "Underdose Probability",
    f"{underdose_prob:.3f}"
)

m4.metric(
    "Coverage ≥95%",
    f"{coverage_95:.3f}"
)

# ==========================================================
# TABS
# ==========================================================

tab1, tab2, tab3, tab4 = st.tabs(
    [
        "CT",
        "Tissue",
        "Deterministic",
        "Monte Carlo"
    ]
)

with tab1:

    fig, ax = plt.subplots()

    ax.imshow(hu, cmap="gray")
    ax.set_title("CT Slice (HU)")
    ax.axis("off")

    st.pyplot(fig)

with tab2:

    fig, ax = plt.subplots()

    ax.imshow(tissue, cmap="viridis")

    ax.contour(
        tumor_mask,
        levels=[0.5]
    )

    ax.set_title("Tissue Segmentation")
    ax.axis("off")

    st.pyplot(fig)

with tab3:

    fig, ax = plt.subplots()

    ax.imshow(dose_det, cmap="hot")

    ax.set_title("Deterministic Dose")
    ax.axis("off")

    st.pyplot(fig)

with tab4:

    col1, col2 = st.columns(2)

    with col1:

        fig, ax = plt.subplots()

        ax.imshow(dose_mc, cmap="hot")

        ax.contour(
            tumor_mask,
            levels=[0.5]
        )

        ax.set_title("Monte Carlo Dose")
        ax.axis("off")

        st.pyplot(fig)

    with col2:

        fig, ax = plt.subplots()

        ax.hist(
            tumor_dose,
            bins=25
        )

        ax.set_title(
            "Tumor Dose Histogram"
        )

        ax.set_xlabel("Normalized Dose")
        ax.set_ylabel("Frequency")

        st.pyplot(fig)

# ==========================================================
# FOOTER
# ==========================================================

st.markdown("---")

st.markdown(
    f"""
**RadSim2D Tool**

Author: {AUTHOR}

Educational and research software only.  
Not for clinical use.
"""
)
