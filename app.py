import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pydicom

# =========================================================
# CONFIG
# =========================================================
st.set_page_config(
    page_title="RadSim2D Tool",
    layout="wide"
)

AUTHOR_NAME = "Williams Stonard Kaphika"
PROJECT_NAME = "RadSim2D Tool"

# =========================================================
# HEADER
# =========================================================
st.title(f"{PROJECT_NAME}: CT-Based Dose Simulation")
st.caption(f"Developed by {AUTHOR_NAME}")

st.info(
    "Educational research tool only — NOT for clinical use."
)

# =========================================================
# SIDEBAR
# =========================================================
st.sidebar.header("Controls")

n_photons = st.sidebar.slider(
    "Number of Photons",
    1000, 20000, 5000, step=1000
)

beam_sigma = st.sidebar.slider(
    "Beam Spread (Sigma)",
    1, 10, 3
)

# =========================================================
# FILE UPLOAD
# =========================================================
uploaded_file = st.file_uploader("Upload CT DICOM", type=["dcm"])

# =========================================================
# MONTE CARLO MODEL
# =========================================================
@st.cache_data(show_spinner=True)
def monte_carlo(mu_map, r_center, c_center, n_photons, beam_sigma, tumor_mask):

    rows, cols = mu_map.shape
    dose_map = np.zeros((rows, cols))

    for _ in range(n_photons):

        x = int(np.random.normal(c_center, beam_sigma))
        x = np.clip(x, 0, cols - 1)

        y = 0
        energy = 1.0

        while energy > 0.001 and y < rows:

            mu = mu_map[y, x]
            p_interact = 1 - np.exp(-mu)

            if np.random.random() < p_interact:
                dose_map[y, x] += energy * (0.2 if tumor_mask[y, x] else 0.1)
                energy *= 0.9

            x += np.random.choice([-1, 0, 1], p=[0.1, 0.8, 0.1])
            x = np.clip(x, 0, cols - 1)

            y += 1

    return dose_map


# =========================================================
# DETERMINISTIC MODEL
# =========================================================
def deterministic(mu_map):

    rows, cols = mu_map.shape
    dose = np.zeros_like(mu_map)

    I0 = 100

    for col in range(cols):

        beam_col = col
        intensity = I0

        for row in range(rows):

            mu = mu_map[row, beam_col]
            att = np.exp(-mu)

            deposited = intensity * (1 - att)
            dose[row, beam_col] += deposited

            intensity *= att

            beam_col = np.clip(
                beam_col + np.random.choice([-1, 0, 1]),
                0, cols - 1
            )

            if intensity < 1e-6:
                break

    return dose


# =========================================================
# MAIN APP
# =========================================================
if uploaded_file is not None:

    try:
        ds = pydicom.dcmread(uploaded_file)
        image = ds.pixel_array.astype(float)
    except Exception as e:
        st.error(f"Invalid DICOM file: {e}")
        st.stop()

    # Downsample for speed
    image = image[::2, ::2]

    slope = getattr(ds, "RescaleSlope", 1.0)
    intercept = getattr(ds, "RescaleIntercept", 0.0)

    hu = image * slope + intercept

    # Tissue segmentation
    tissue = np.zeros_like(hu)
    tissue[hu < -500] = 0
    tissue[(hu >= -500) & (hu < 300)] = 1
    tissue[hu >= 300] = 2

    # Attenuation map
    mu_map = np.zeros_like(tissue, dtype=float)
    mu_map[tissue == 0] = 0.02
    mu_map[tissue == 1] = 0.20
    mu_map[tissue == 2] = 0.50

    rows, cols = mu_map.shape

    # Tumor selection (safe)
    soft = np.argwhere(tissue == 1)

    if len(soft) == 0:
        st.error("No soft tissue detected.")
        st.stop()

    r_center, c_center = soft[len(soft) // 2]

    tumor_rows = range(max(0, r_center - 10), min(rows, r_center + 10))
    tumor_cols = range(max(0, c_center - 10), min(cols, c_center + 10))

    tumor_mask = np.zeros_like(tissue, dtype=bool)
    tumor_mask[np.ix_(tumor_rows, tumor_cols)] = True

    # Run models
    dose_det = deterministic(mu_map)
    dose_mc = monte_carlo(
        mu_map,
        r_center,
        c_center,
        n_photons,
        beam_sigma,
        tumor_mask
    )

    # Normalize
    dose_norm = dose_mc / (np.max(dose_mc) + 1e-8)
    tumor = dose_norm[np.ix_(tumor_rows, tumor_cols)].ravel()

    # Tabs
    tab1, tab2, tab3 = st.tabs(["CT", "Deterministic", "Monte Carlo"])

    # =====================================================
    # CT VIEW
    # =====================================================
    with tab1:
        col1, col2 = st.columns(2)

        with col1:
            fig, ax = plt.subplots()
            ax.imshow(hu, cmap="gray")
            ax.set_title("CT (HU)")
            ax.axis("off")
            st.pyplot(fig)

        with col2:
            fig, ax = plt.subplots()
            ax.imshow(tissue, cmap="viridis")
            ax.set_title("Tissue Map")
            ax.axis("off")
            st.pyplot(fig)

    # =====================================================
    # DETERMINISTIC
    # =====================================================
    with tab2:
        col1, col2 = st.columns(2)

        with col1:
            fig, ax = plt.subplots()
            ax.imshow(dose_det, cmap="hot")
            ax.set_title("Deterministic Dose")
            ax.axis("off")
            st.pyplot(fig)

        with col2:
            st.metric("Mean Dose", float(np.mean(dose_det)))
            st.metric("Max Dose", float(np.max(dose_det)))
            st.metric("Min Dose", float(np.min(dose_det)))

    # =====================================================
    # MONTE CARLO
    # =====================================================
    with tab3:
        col1, col2 = st.columns(2)

        with col1:
            fig, ax = plt.subplots()
            ax.imshow(dose_mc, cmap="hot")
            ax.set_title("Monte Carlo Dose")
            ax.axis("off")
            st.pyplot(fig)

        with col2:
            st.metric("Mean Tumor Dose", float(np.mean(tumor)))
            st.metric("Underdose Probability", float(np.mean(tumor < 0.3)))

            fig, ax = plt.subplots()
            ax.hist(tumor, bins=20)
            ax.set_title("Tumor Dose Distribution")
            st.pyplot(fig)

else:
    st.info("Upload a CT DICOM file to begin simulation.")
