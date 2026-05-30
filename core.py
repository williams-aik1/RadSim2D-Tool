import numpy as np

def deterministic(mu_map):
    rows, cols = mu_map.shape
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
            beam_col = np.clip(beam_col + np.random.choice([-1, 0, 1]), 0, cols - 1)

            if intensity < 1e-6:
                break

    return dose


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
                deposited = energy * (0.2 if tumor_mask[y, x] else 0.1)
                dose_map[y, x] += deposited
                energy *= 0.9

            x += np.random.choice([-1, 0, 1], p=[0.1, 0.8, 0.1])
            x = np.clip(x, 0, cols - 1)
            y += 1

    return dose_map
