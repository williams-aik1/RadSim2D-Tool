import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pydicom

# -----------------------------
# Page Config
# -----------------------------
st.set_page_config(
    layout="wide",
    page_title="RadSim2D Tool"
)

st.title("RadSim2D Tool: Interactive Dose Planning")

# -----------------------------
# Sidebar Controls
# -----------------------------
st.sidebar.header("Simulation Controls")

n_photons = st.sidebar.slider("Number of Photons", 1000, 20000, 5000, step=1000)
beam_sigma = st.sidebar.slider("Beam Focus (Sigma)", 1, 10, 3)

# -----------------------------
# Upload CT
# -----------------------------
uploaded_file = st.file_uploader("Upload CT DICOM (.dcm)", type=["dcm"])


# =========================================================
# Monte Carlo Simulation (cached for performance)
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

                if tumor_mask[y, x]:
                    deposited = energy * 0.2
                else:
                    deposited = energy * 0.1

                dose_map[y, x] += deposited
                energy *= 0.9

            x += np.random.choice([-1, 0, 1], p=[0.1, 0.8, 0.1])
            x = np.clip(x, 0, cols - 1)

            y += 1

    return dose_map


# =========================================================
# Run App
# =========================================================
if uploaded_file is not None:

    ds = pydicom.dcmread(uploaded_file)
    image = ds.pixel_array.astype(float)

    # downsample for speed
    image = image[::2, ::2]

    # HU conversion (safe fallback)
    slope = getattr(ds, "RescaleSlope", 1.0)
    intercept = getattr(ds, "RescaleIntercept", 0.0)
    hu = image * slope + intercept

    # -----------------------------
    # Tissue Segmentation
    # -----------------------------
    tissue = np.zeros_like(hu)

    tissue[hu < -500] = 0
    tissue[(hu >= -500) & (hu < 300)] = 1
    tissue[hu >= 300] = 2

    # -----------------------------
    # Attenuation Map
    # -----------------------------
    mu_map = np.zeros_like(tissue, dtype=float)
    mu_map[tissue == 0] = 0.02
    mu_map[tissue == 1] = 0.20
    mu_map[tissue == 2] = 0.50

    rows, cols = mu_map.shape

    # -----------------------------
    # Tumor region (SAFE selection)
    # -----------------------------
    soft = np.argwhere(tissue == 1)

    if len(soft) == 0:
        st.error("No soft tissue detected in image.")
        st.stop()

    r_center, c_center = soft[len(soft)//2]

    tumor_rows = range(max(0, r_center-10), min(rows, r_center+10))
    tumor_cols = range(max(0, c_center-10), min(cols, c_center+10))

    tumor_mask = np.zeros_like(tissue, dtype=bool)
    tumor_mask[np.ix_(tumor_rows, tumor_cols)] = True

    # =========================================================
    # Deterministic model
    # =========================================================
    def deterministic(mu_map):
        I0 = 100
        dose = np.zeros_like(mu_map)

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


    dose_det = deterministic(mu_map)
    dose_mc = monte_carlo(mu_map, r_center, c_center, n_photons, beam_sigma, tumor_mask)

    # =========================================================
    # UI Layout
    # =========================================================
    tab1, tab2, tab3 = st.tabs(["CT + Tissue", "Deterministic", "Monte Carlo"])

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
            ax.set_title("Tissue Map")
            ax.axis("off")

            for r in tumor_rows:
                for c in tumor_cols:
                    ax.plot(c, r, 'r.', markersize=1)

            st.pyplot(fig)

    # -----------------------------
    # TAB 2
    # -----------------------------
    with tab2:
        col1, col2 = st.columns(2)

        with col1:
            fig, ax = plt.subplots()
            ax.imshow(dose_det, cmap="hot")
            ax.set_title("Deterministic Dose")
            ax.axis("off")
            st.pyplot(fig)

        with col2:
            st.write("Mean:", np.mean(dose_det))
            st.write("Max:", np.max(dose_det))
            st.write("Min:", np.min(dose_det))

    # -----------------------------
    # TAB 3
    # -----------------------------
    with tab3:
        col1, col2 = st.columns(2)

        with col1:
            fig, ax = plt.subplots()
            ax.imshow(dose_mc, cmap="hot")
            ax.set_title("Monte Carlo Dose")
            ax.axis("off")

            for r in tumor_rows:
                for c in tumor_cols:
                    ax.plot(c, r, 'b.', markersize=1)

            st.pyplot(fig)

        with col2:
            dose_norm = dose_mc / (np.max(dose_mc) + 1e-8)
            tumor = dose_norm[np.ix_(tumor_rows, tumor_cols)].ravel()

            st.write("Mean Tumor Dose:", np.mean(tumor))
            st.write("Underdose Probability:", np.mean(tumor < 0.3))

            fig, ax = plt.subplots()
            ax.hist(tumor, bins=20)
            ax.set_title("Tumor Dose Distribution")
            st.pyplot(fig)

else:
    st.info("Upload a CT DICOM file to begin simulation.")
