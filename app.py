import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pydicom

from core import deterministic, monte_carlo
from utils import to_hu, segment_tissue, tumor_region

# =========================================================
# CONFIG
# =========================================================
st.set_page_config(page_title="RadSim2D", layout="wide")

AUTHOR = "Williams Stonard Kaphika"

st.title("RadSim2D Tool")
st.caption("CT-based Monte Carlo Radiation Simulation")

st.info("Educational research tool — not for clinical use.")

# =========================================================
# SIDEBAR
# =========================================================
st.sidebar.header("Controls")

n_photons = st.sidebar.slider("Photons", 1000, 20000, 5000)
sigma = st.sidebar.slider("Beam Spread", 1, 10, 3)

st.sidebar.markdown("---")

# ✅ REQUIRED: bottom-left watermark
st.sidebar.markdown(
    f"""
    <div style="position:fixed; bottom:10px; left:20px; opacity:0.6; font-size:12px;">
    Developed by <b>{AUTHOR}</b>
    </div>
    """,
    unsafe_allow_html=True
)

# =========================================================
# FILE UPLOAD
# =========================================================
file = st.file_uploader("Upload CT DICOM", type=["dcm"])

if file:

    ds = pydicom.dcmread(file)
    img = ds.pixel_array.astype(float)

    img = img[::2, ::2]

    hu = to_hu(img, ds)

    tissue, mu_map = segment_tissue(hu)

    tumor_mask = tumor_region(tissue)

    if tumor_mask is None:
        st.error("No soft tissue found.")
        st.stop()

    # =====================================================
    # RUN MODELS
    # =====================================================
    dose_det = deterministic(mu_map)
    dose_mc, unc = monte_carlo(mu_map, tumor_mask, n_photons, sigma)

    dose_norm = dose_mc / (np.max(dose_mc) + 1e-8)

    tumor_vals = dose_norm[tumor_mask]

    # =====================================================
    # TABS
    # =====================================================
    t1, t2, t3, t4 = st.tabs(["CT", "Deterministic", "Monte Carlo", "Analytics"])

    with t1:
        col1, col2 = st.columns(2)

        with col1:
            plt.imshow(hu, cmap="gray")
            plt.title("CT (HU)")
            plt.axis("off")
            st.pyplot(plt.gcf())

        with col2:
            plt.imshow(tissue, cmap="viridis")
            plt.title("Tissue")
            plt.axis("off")
            st.pyplot(plt.gcf())

    with t2:
        plt.imshow(dose_det, cmap="hot")
        plt.title("Deterministic Dose")
        plt.axis("off")
        st.pyplot(plt.gcf())

    with t3:
        plt.imshow(dose_norm, cmap="hot")
        plt.title("Monte Carlo Dose")
        plt.axis("off")
        st.pyplot(plt.gcf())

    with t4:
        st.metric("Mean Tumor Dose", float(np.mean(tumor_vals)))
        st.metric("Underdose Probability", float(np.mean(tumor_vals < 0.3)))

        plt.imshow(unc, cmap="magma")
        plt.title("Uncertainty Map")
        plt.axis("off")
        st.pyplot(plt.gcf())

else:
    st.info("Upload a DICOM file to begin.")
