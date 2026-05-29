import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pydicom
from matplotlib.patches import Rectangle

st.set_page_config(page_title="RadSim2D Tool", layout="wide")

st.title("☢️ RadSim2D Tool")

st.sidebar.header("Controls")

n_photons = st.sidebar.slider("Photons", 1000, 20000, 5000, step=1000)
beam_sigma = st.sidebar.slider("Beam Sigma", 1, 10, 3)
show_contours = st.sidebar.checkbox("Contours", True)

uploaded_file = st.file_uploader("Upload DICOM", type=["dcm"])


@st.cache_data
def monte_carlo(mu_map, tumor_mask, n_photons, beam_sigma):
    rows, cols = mu_map.shape
    dose = np.zeros_like(mu_map)

    for _ in range(n_photons):
        x = np.random.normal(cols // 2, beam_sigma)
        x = int(np.clip(x, 0, cols - 1))

        for y in range(rows):
            mu = mu_map[y, x]
            if np.random.random() < (1 - np.exp(-mu)):
                dose[y, x] += 1 if tumor_mask[y, x] else 0.5

            x += np.random.choice([-1, 0, 1])
            x = np.clip(x, 0, cols - 1)

    return dose


if uploaded_file:

    ds = pydicom.dcmread(uploaded_file)
    img = ds.pixel_array.astype(float)

    img = img[::2, ::2]

    slope = getattr(ds, "RescaleSlope", 1.0)
    intercept = getattr(ds, "RescaleIntercept", 0.0)

    hu = img * slope + intercept

    tissue = np.zeros_like(hu)
    tissue[(hu >= -500) & (hu < 300)] = 1
    tissue[hu >= 300] = 2

    mu_map = np.zeros_like(hu, dtype=float)
    mu_map[tissue == 0] = 0.02
    mu_map[tissue == 1] = 0.2
    mu_map[tissue == 2] = 0.5

    soft = np.argwhere(tissue == 1)

    if len(soft) == 0:
        st.error("No soft tissue found")
        st.stop()

    r0, c0 = soft[len(soft)//2]

    tumor_mask = np.zeros_like(tissue, dtype=bool)
    tumor_mask[r0-10:r0+10, c0-10:c0+10] = True

    dose = monte_carlo(mu_map, tumor_mask, n_photons, beam_sigma)

    tab1, tab2 = st.tabs(["CT", "Dose"])

    with tab1:
        fig, ax = plt.subplots()
        ax.imshow(hu, cmap="gray")
        ax.set_title("CT")
        ax.axis("off")
        st.pyplot(fig)

    with tab2:
        fig, ax = plt.subplots()
        ax.imshow(dose, cmap="hot")

        ax.add_patch(Rectangle((c0-10, r0-10), 20, 20,
                               edgecolor="red", fill=False))

        st.pyplot(fig)

else:
    st.info("Upload DICOM to start")
