# Elite-Grade `app.py` for RadSim2D

```python
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pydicom
from matplotlib.patches import Rectangle

# =========================================================
# PAGE CONFIG
# =========================================================
st.set_page_config(
    page_title="RadSim2D Tool",
    page_icon="☢️",
    layout="wide"
)

# =========================================================
# HEADER
# =========================================================
st.title("☢️ RadSim2D Tool")
st.subheader(
    "CT-Based Monte Carlo Radiation Dose Simulation Framework"
)

st.markdown(
    """
Developed by **Williams Stonard Kaphika**  
Independent Researcher in Computational Medical Physics
"""
)

# =========================================================
# DISCLAIMER
# =========================================================
st.warning(
    """
    DISCLAIMER:

    This software is intended strictly for educational, research,
    and prototyping purposes only.

    RadSim2D is NOT a clinically validated treatment planning system
    and must NOT be used for medical diagnosis, radiation treatment
    planning, or clinical decision-making.

    All dose values are relative simulation metrics and do not
    represent calibrated clinical dosimetry.
    """
)

# =========================================================
# SIDEBAR
# =========================================================
st.sidebar.header("Simulation Controls")

n_photons = st.sidebar.slider(
    "Number of Photons",
    1000,
    50000,
    15000,
    step=1000
)

beam_sigma = st.sidebar.slider(
    "Beam Spread (Sigma)",
    1,
    15,
    3
)

beam_energy = st.sidebar.selectbox(
    "Photon Energy",
    [
        "Low (50 keV)",
        "Medium (100 keV)",
        "High (150 keV)"
    ]
)

beam_mode = st.sidebar.selectbox(
    "Beam Strategy",
    [
        "Single Beam",
        "Multi-Angle"
    ]
)

show_contours = st.sidebar.checkbox(
    "Show Isodose Contours",
    value=True
)

# =========================================================
# ENERGY FACTORS
# =========================================================
energy_scale = {
    "Low (50 keV)": 0.8,
    "Medium (100 keV)": 1.0,
    "High (150 keV)": 1.2
}

energy_factor = energy_scale[beam_energy]

# =========================================================
# FILE UPLOAD
# =========================================================
uploaded_file = st.file_uploader(
    "Upload CT DICOM (.dcm)",
    type=["dcm"]
)

# =========================================================
# MONTE CARLO ENGINE
# =========================================================
@st.cache_data(show_spinner=True)
def monte_carlo(
    mu_map,
    tumor_mask,
    n_photons,
    beam_sigma,
    energy_factor,
    multi_angle=False
):

    rows, cols = mu_map.shape

    dose_map = np.zeros((rows, cols))
    variance_map = np.zeros((rows, cols))

    beam_angles = [0]

    if multi_angle:
        beam_angles = [0, 30, -30]

    for angle in beam_angles:

        theta0 = np.deg2rad(angle)

        for _ in range(n_photons // len(beam_angles)):

            x = cols // 2 + int(np.random.normal(0, beam_sigma))
            x = np.clip(x, 0, cols - 1)

            y = 0
            theta = theta0
            energy = 1.0 * energy_factor

            while energy > 0.001:

                xi = np.random.random()

                mu = mu_map[int(y), int(x)] * energy_factor

                step = -np.log(max(xi, 1e-8)) / max(mu, 1e-8)

                dx = np.sin(theta) * step
                dy = np.cos(theta) * step

                x += dx
                y += dy

                x_int = int(np.clip(x, 0, cols - 1))
                y_int = int(np.clip(y, 0, rows - 1))

                if y_int >= rows - 1 or x_int <= 0 or x_int >= cols - 1:
                    break

                p_interact = 1 - np.exp(-mu)

                if np.random.random() < p_interact:

                    if tumor_mask[y_int, x_int]:
                        deposited = energy * 0.25
                    else:
                        deposited = energy * 0.08

                    dose_map[y_int, x_int] += deposited
                    variance_map[y_int, x_int] += deposited ** 2

                    scatter = np.random.normal(0, 0.08)
                    theta += scatter

                    energy *= np.random.uniform(0.88, 0.97)

    uncertainty = np.sqrt(
        variance_map / max(n_photons, 1)
    )

    return dose_map, uncertainty

# =========================================================
# DETERMINISTIC MODEL
# =========================================================
def deterministic(mu_map, energy_factor):

    rows, cols = mu_map.shape

    dose = np.zeros_like(mu_map)

    I0 = 100 * energy_factor

    for col in range(cols):

        intensity = I0
        beam_col = col

        for row in range(rows):

            mu = mu_map[row, beam_col] * energy_factor

            attenuation = np.exp(-mu)

            deposited = intensity * (1 - attenuation)

            dose[row, beam_col] += deposited

            intensity *= attenuation

            beam_col = np.clip(
                beam_col + np.random.choice([-1, 0, 1]),
                0,
                cols - 1
            )

            if intensity < 1e-6:
                break

    return dose

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

    # =====================================================
    # TISSUE SEGMENTATION
    # =====================================================
    tissue = np.zeros_like(hu)

    tissue[hu < -500] = 0
    tissue[(hu >= -500) & (hu < 300)] = 1
    tissue[hu >= 300] = 2

    # =====================================================
    # ATTENUATION MAP
    # =====================================================
    mu_map = np.zeros_like(tissue, dtype=float)

    mu_map[tissue == 0] = 0.02
    mu_map[tissue == 1] = 0.20
    mu_map[tissue == 2] = 0.50

    rows, cols = mu_map.shape

    # =====================================================
    # TUMOR REGION
    # =====================================================
    soft = np.argwhere(tissue == 1)

    if len(soft) == 0:
        st.error("No soft tissue detected.")
        st.stop()

    r_center, c_center = soft[len(soft)//2]

    tumor_size = 20

    tumor_rows = range(
        max(0, r_center - tumor_size//2),
        min(rows, r_center + tumor_size//2)
    )

    tumor_cols = range(
        max(0, c_center - tumor_size//2),
        min(cols, c_center + tumor_size//2)
    )

    tumor_mask = np.zeros_like(tissue, dtype=bool)

    tumor_mask[np.ix_(tumor_rows, tumor_cols)] = True

    # =====================================================
    # RUN SIMULATIONS
    # =====================================================
    dose_det = deterministic(mu_map, energy_factor)

    dose_mc, uncertainty = monte_carlo(
        mu_map,
        tumor_mask,
        n_photons,
        beam_sigma,
        energy_factor,
        multi_angle=(beam_mode == "Multi-Angle")
    )

    dose_norm = dose_mc / (np.max(dose_mc) + 1e-8)

    tumor_dose = dose_norm[np.ix_(tumor_rows, tumor_cols)].ravel()

    # =====================================================
    # METRICS
    # =====================================================
    mean_dose = np.mean(tumor_dose)
    max_dose = np.max(tumor_dose)
    min_dose = np.min(tumor_dose)
    underdose_probability = np.mean(tumor_dose < 0.3)
    coverage = np.mean(tumor_dose > 0.5)

    # =====================================================
    # TABS
    # =====================================================
    tab1, tab2, tab3, tab4 = st.tabs([
        "CT + Tissue",
        "Deterministic",
        "Monte Carlo",
        "Analytics"
    ])

    # =====================================================
    # TAB 1
    # =====================================================
    with tab1:

        col1, col2 = st.columns(2)

        with col1:

            fig, ax = plt.subplots(figsize=(6,6))
            ax.imshow(hu, cmap="gray")
            ax.set_title("CT Image (HU)")
            ax.axis("off")

            st.pyplot(fig)

        with col2:

            fig, ax = plt.subplots(figsize=(6,6))
            ax.imshow(tissue, cmap="viridis")
            ax.set_title("Tissue Segmentation")
            ax.axis("off")

            rect = Rectangle(
                (min(tumor_cols), min(tumor_rows)),
                tumor_size,
                tumor_size,
                linewidth=2,
                edgecolor='red',
                facecolor='none'
            )

            ax.add_patch(rect)

            st.pyplot(fig)

    # =====================================================
    # TAB 2
    # =====================================================
    with tab2:

        col1, col2 = st.columns([2,1])

        with col1:

            fig, ax = plt.subplots(figsize=(7,7))
            img = ax.imshow(dose_det, cmap="hot")
            ax.set_title("Deterministic Dose Distribution")
            ax.axis("off")

            plt.colorbar(img, ax=ax)

            st.pyplot(fig)

        with col2:

            st.metric("Mean Dose", f"{np.mean(dose_det):.4f}")
            st.metric("Maximum Dose", f"{np.max(dose_det):.4f}")
            st.metric("Minimum Dose", f"{np.min(dose_det):.4f}")

    # =====================================================
    # TAB 3
    # =====================================================
    with tab3:

        col1, col2 = st.columns([2,1])

        with col1:

            fig, ax = plt.subplots(figsize=(7,7))

            img = ax.imshow(dose_norm, cmap="hot")

            if show_contours:
                ax.contour(
                    dose_norm,
                    levels=[0.2, 0.5, 0.8],
                    colors='cyan',
                    linewidths=1
                )

            rect = Rectangle(
                (min(tumor_cols), min(tumor_rows)),
                tumor_size,
                tumor_size,
                linewidth=2,
                edgecolor='blue',
                facecolor='none'
            )

            ax.add_patch(rect)

            ax.set_title("Monte Carlo Dose Distribution")
            ax.axis("off")

            plt.colorbar(img, ax=ax)

            st.pyplot(fig)

        with col2:

            st.metric("Mean Tumor Dose", f"{mean_dose:.3f}")
            st.metric("Maximum Tumor Dose", f"{max_dose:.3f}")
            st.metric("Minimum Tumor Dose", f"{min_dose:.3f}")
            st.metric("Underdose Probability", f"{underdose_probability:.3f}")
            st.metric("Coverage (>0.5)", f"{coverage:.3f}")

    # =====================================================
    # TAB 4
    # =====================================================
    with tab4:

        col1, col2 = st.columns(2)

        with col1:

            fig, ax = plt.subplots(figsize=(6,4))
            ax.hist(tumor_dose, bins=25)
            ax.set_title("Tumor Dose Histogram")
            ax.set_xlabel("Normalized Dose")
            ax.set_ylabel("Voxel Count")

            st.pyplot(fig)

        with col2:

            sorted_dose = np.sort(tumor_dose)
            dvh = 1 - np.arange(len(sorted_dose))/len(sorted_dose)

            fig, ax = plt.subplots(figsize=(6,4))
            ax.plot(sorted_dose, dvh)
            ax.set_title("Dose Volume Histogram (DVH)")
            ax.set_xlabel("Normalized Dose")
            ax.set_ylabel("Volume Fraction")
            ax.grid(True)

            st.pyplot(fig)

        st.subheader("Monte Carlo Uncertainty Map")

        fig, ax = plt.subplots(figsize=(8,6))
        img = ax.imshow(uncertainty, cmap="magma")
        ax.set_title("Statistical Uncertainty")
        ax.axis("off")

        plt.colorbar(img, ax=ax)

        st.pyplot(fig)

else:

    st.info("Upload a DICOM CT file to begin simulation.")

# =========================================================
# FOOTER
# =========================================================
st.markdown("---")

st.caption(
    "RadSim2D Tool © 2026 | Williams Stonard Kaphika"
)

```
