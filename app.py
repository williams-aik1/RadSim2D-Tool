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
st.title(f"{PROJECT_NAME}")
st.caption("CT-Based Deterministic vs Monte Carlo Energy Deposition Simulator")

st.info(
    "Educational research simulator only — NOT for clinical use."
)

# =========================================================
# SIDEBAR
# =========================================================
st.sidebar.header("Simulation Controls")

n_photons = st.sidebar.slider("Number of Photons", 1000, 20000, 5000, step=1000)
beam_sigma = st.sidebar.slider("Beam Spread (Sigma)", 1, 10, 3)

st.sidebar.markdown("---")

# =========================================================
# WATERMARK (BOTTOM LEFT UNDER SIDEBAR)
# =========================================================
st.sidebar.markdown(
    f"""
    <div style="
        position: fixed;
        bottom: 10px;
        left: 20px;
        font-size: 12px;
        opacity: 0.6;
    ">
        Developed by <b>{AUTHOR_NAME}</b>
    </div>
    """,
    unsafe_allow_html=True
)

# =========================================================
# FILE UPLOAD
# =========================================================
uploaded_file = st.file_uploader("Upload CT DICOM (.dcm)", type=["dcm"])

# =========================================================
# PHYSICS MODELS
# =========================================================

def deterministic(mu_map):
    rows, cols = mu_map.shape
    dose = np.zeros_like(mu_map)

    for x in range(cols):
        intensity = 1.0

        for y in range(rows):
            mu = mu_map[y, x]
            intensity *= np.exp(-mu)
            dose[y, x] = 1 - intensity

    return dose


def monte_carlo(mu_map, tumor_mask, n_photons, sigma, seed=42):

    rng = np.random.default_rng(seed)

    rows, cols = mu_map.shape
    dose = np.zeros((rows, cols))
    variance = np.zeros((rows, cols))

    for _ in range(n_photons):

        x = int(np.clip(rng.normal(cols // 2, sigma), 0, cols - 1))
        y = 0
        energy = 1.0

        while y < rows - 1 and energy > 1e-3:

            mu = mu_map[y, x]

            step = -np.log(rng.random() + 1e-12) / (mu + 1e-12)
            y_new = min(rows - 1, int(y + step))

            deposited = energy * (0.2 if tumor_mask[y_new, x] else 0.1)

            dose[y_new, x] += deposited
            variance[y_new, x] += deposited ** 2

            energy *= (1 - deposited)

            x += rng.choice([-1, 0, 1], p=[0.15, 0.7, 0.15])
            x = np.clip(x, 0, cols - 1)

            y = y_new

    uncertainty = np.sqrt(variance / max(n_photons, 1))

    return dose, uncertainty


# =========================================================
# MAIN PIPELINE
# =========================================================
if uploaded_file is not None:

    ds = pydicom.dcmread(uploaded_file)
    image = ds.pixel_array.astype(float)

    image = image[::2, ::2]

    slope = getattr(ds, "RescaleSlope", 1.0)
    intercept = getattr(ds, "RescaleIntercept", 0.0)

    hu = image * slope + intercept

    # -----------------------------
    # TISSUE MODEL
    # -----------------------------
    tissue = np.zeros_like(hu)
    tissue[hu < -500] = 0
    tissue[(hu >= -500) & (hu < 300)] = 1
    tissue[hu >= 300] = 2

    mu_map = np.zeros_like(tissue, dtype=float)
    mu_map[tissue == 0] = 0.02
    mu_map[tissue == 1] = 0.20
    mu_map[tissue == 2] = 0.50

    rows, cols = mu_map.shape

    # -----------------------------
    # ROI (tumor surrogate)
    # -----------------------------
    soft = np.argwhere(tissue == 1)

    if len(soft) == 0:
        st.error("No soft tissue detected.")
        st.stop()

    r_center, c_center = soft[len(soft) // 2]

    tumor_mask = np.zeros_like(tissue, dtype=bool)
    tumor_mask[r_center-10:r_center+10, c_center-10:c_center+10] = True

    # -----------------------------
    # RUN SIMULATIONS
    # -----------------------------
    dose_det = deterministic(mu_map)
    dose_mc, uncertainty = monte_carlo(
        mu_map,
        tumor_mask,
        n_photons,
        beam_sigma
    )

    dose_norm = dose_mc / (np.max(dose_mc) + 1e-8)

    tumor_vals = dose_norm[tumor_mask]

    # =====================================================
    # TABS
    # =====================================================
    tab1, tab2, tab3, tab4 = st.tabs(
        ["CT + Tissue", "Deterministic", "Monte Carlo", "Analytics"]
    )

    # -----------------------------
    # TAB 1
    # -----------------------------
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
            ax.set_title("Tissue Segmentation")
            ax.axis("off")
            st.pyplot(fig)

    # -----------------------------
    # TAB 2
    # -----------------------------
    with tab2:
        col1, col2 = st.columns(2)

        with col1:
            fig, ax = plt.subplots()
            ax.imshow(dose_det, cmap="hot")
            ax.set_title("Deterministic Energy Deposition")
            ax.axis("off")
            st.pyplot(fig)

        with col2:
            st.metric("Mean", float(np.mean(dose_det)))
            st.metric("Max", float(np.max(dose_det)))
            st.metric("Min", float(np.min(dose_det)))

    # -----------------------------
    # TAB 3
    # -----------------------------
    with tab3:
        col1, col2 = st.columns(2)

        with col1:
            fig, ax = plt.subplots()
            ax.imshow(dose_norm, cmap="hot")
            ax.set_title("Monte Carlo Energy Deposition")
            ax.axis("off")
            st.pyplot(fig)

        with col2:
            st.metric("Mean Tumor Dose", float(np.mean(tumor_vals)))
            st.metric("Underdose Probability", float(np.mean(tumor_vals < 0.3)))

            fig, ax = plt.subplots()
            ax.hist(tumor_vals, bins=20)
            ax.set_title("Tumor Dose Distribution")
            st.pyplot(fig)

    # -----------------------------
    # TAB 4 (RESEARCH OUTPUT)
    # -----------------------------
    with tab4:
        st.subheader("Statistical Summary")

        st.write({
            "Tumor Mean Dose": float(np.mean(tumor_vals)),
            "Tumor Std": float(np.std(tumor_vals)),
            "Underdose Probability": float(np.mean(tumor_vals < 0.3)),
        })

        st.subheader("Uncertainty Map")
        fig, ax = plt.subplots()
        ax.imshow(uncertainty, cmap="magma")
        ax.set_title("Monte Carlo Uncertainty")
        ax.axis("off")
        st.pyplot(fig)

else:
    st.info("Upload a CT DICOM file to begin simulation.")
