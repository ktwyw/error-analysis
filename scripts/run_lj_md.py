"""Molecular dynamics of a Lennard-Jones liquid in pure NumPy.

Generates the data analysed in notebooks/part5_examples/23_md_lennard_jones.ipynb:

    data/lj_md_timeseries.csv   time series of U/N, E/N, T and P
    data/lj_md_msd_blocks.csv   mean-square displacement in consecutive blocks
    data/lj_md_parameters.json  run parameters

Model: truncated-and-shifted Lennard-Jones potential (cutoff 2.5 sigma),
N = 256 particles, reduced density 0.8, temperature close to 1.0, cubic
periodic box, velocity-Verlet integration in the NVE ensemble after an
equilibration with velocity rescaling.  All quantities are in reduced
Lennard-Jones units.

Run from the repository root:  python scripts/run_lj_md.py
(one to two minutes on a laptop).
"""

import json
import pathlib
import time

import numpy as np

# ------------------------------------------------------------------ parameters
N_CELLS = 4            # fcc cells per edge -> N = 4 * N_CELLS**3 = 256
RHO = 0.8
T_TARGET = 1.0
R_CUT = 2.5
DT = 0.005
N_EQUIL = 5000         # velocity rescaling every 10 steps
N_RELAX = 2000         # NVE, discarded
N_PROD = 40000         # NVE, analysed
SAMPLE_EVERY = 1       # thermodynamic sampling interval (steps)
N_MSD_BLOCKS = 10
MSD_EVERY = 10         # position sampling interval for the MSD (steps)
MSD_MAX_LAG = 200      # frames -> 10 time units
SEED = 20261112

rng = np.random.default_rng(SEED)
N = 4 * N_CELLS**3
L = (N / RHO) ** (1.0 / 3.0)
V = L**3
IU, JU = np.triu_indices(N, k=1)
RC2 = R_CUT**2
U_SHIFT = 4.0 * (R_CUT**-12 - R_CUT**-6)


def fcc_lattice():
    a = L / N_CELLS
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    cells = np.array([[i, j, k] for i in range(N_CELLS)
                      for j in range(N_CELLS) for k in range(N_CELLS)])
    return ((cells[:, None, :] + basis[None, :, :]).reshape(-1, 3)) * a


def forces(r):
    """Forces, potential energy and virial for positions r (unwrapped)."""
    d = r[IU] - r[JU]
    d -= L * np.round(d / L)
    r2 = np.einsum("ij,ij->i", d, d)
    m = r2 < RC2
    d, r2 = d[m], r2[m]
    inv2 = 1.0 / r2
    inv6 = inv2**3
    fr = 24.0 * inv6 * (2.0 * inv6 - 1.0) * inv2      # f / r
    fvec = fr[:, None] * d
    F = np.empty((N, 3))
    iu, ju = IU[m], JU[m]
    for k in range(3):
        F[:, k] = (np.bincount(iu, weights=fvec[:, k], minlength=N)
                   - np.bincount(ju, weights=fvec[:, k], minlength=N))
    U = np.sum(4.0 * inv6 * (inv6 - 1.0) - U_SHIFT)
    W = np.sum(fr * r2)                               # sum of f_ij * r_ij
    return F, U, W


def temperature(v):
    return np.sum(v**2) / (3 * N - 3)


def main():
    t0 = time.time()
    r = fcc_lattice()
    v = rng.normal(size=(N, 3))
    v -= v.mean(axis=0)
    v *= np.sqrt(T_TARGET / temperature(v))
    F, U, W = forces(r)

    def step(r, v, F):
        v = v + 0.5 * DT * F
        r = r + DT * v
        F, U, W = forces(r)
        v = v + 0.5 * DT * F
        return r, v, F, U, W

    for i in range(N_EQUIL):
        r, v, F, U, W = step(r, v, F)
        if i % 10 == 0:
            v *= np.sqrt(T_TARGET / temperature(v))
    for i in range(N_RELAX):
        r, v, F, U, W = step(r, v, F)

    n_samples = N_PROD // SAMPLE_EVERY
    series = np.empty((n_samples, 5))
    block_steps = N_PROD // N_MSD_BLOCKS
    frames_per_block = block_steps // MSD_EVERY
    frames = np.empty((frames_per_block, N, 3))
    msd = np.empty((MSD_MAX_LAG + 1, N_MSD_BLOCKS))

    for i in range(N_PROD):
        if i % MSD_EVERY == 0:
            frames[(i % block_steps) // MSD_EVERY] = r
        r, v, F, U, W = step(r, v, F)
        if (i + 1) % SAMPLE_EVERY == 0:
            T = temperature(v)
            K = 0.5 * np.sum(v**2)
            P = RHO * T + W / (3.0 * V)
            series[(i + 1) // SAMPLE_EVERY - 1] = (
                (i + 1) * DT, U / N, (U + K) / N, T, P)
        if (i + 1) % block_steps == 0:
            b = (i + 1) // block_steps - 1
            msd[0, b] = 0.0
            for lag in range(1, MSD_MAX_LAG + 1):
                dr = frames[lag:] - frames[:-lag]
                msd[lag, b] = np.mean(np.sum(dr**2, axis=2))

    out = pathlib.Path(__file__).resolve().parents[1] / "data"
    out.mkdir(exist_ok=True)
    np.savetxt(out / "lj_md_timeseries.csv", series, delimiter=",",
               fmt=["%.3f", "%.5f", "%.6f", "%.5f", "%.5f"],
               header="time,U_per_N,E_per_N,T,P", comments="")
    lag_time = np.arange(MSD_MAX_LAG + 1) * MSD_EVERY * DT
    np.savetxt(out / "lj_md_msd_blocks.csv", np.column_stack([lag_time, msd]),
               delimiter=",", fmt="%.6f",
               header="lag_time," + ",".join(
                   f"msd_block_{b}" for b in range(N_MSD_BLOCKS)), comments="")
    params = dict(
        model="Lennard-Jones, truncated and shifted", N=N, density=RHO,
        box_length=L, r_cut=R_CUT, dt=DT, n_equil=N_EQUIL, n_relax=N_RELAX,
        n_prod=N_PROD, sample_every=SAMPLE_EVERY, n_msd_blocks=N_MSD_BLOCKS,
        msd_block_length=block_steps * DT, msd_every=MSD_EVERY, seed=SEED,
        units="reduced Lennard-Jones units")
    (out / "lj_md_parameters.json").write_text(json.dumps(params, indent=2) + "\n")
    print(f"done in {time.time() - t0:.0f} s;  <T> = {series[:, 3].mean():.4f}, "
          f"<P> = {series[:, 4].mean():.4f}, "
          f"energy drift = {series[-1, 2] - series[0, 2]:+.2e} per particle")


if __name__ == "__main__":
    main()
