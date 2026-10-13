"""
Transmittance and Reflectance Simulation of a Polystyrene Sphere Monolayer
-------------------------------------------------------------------------
This script uses Meep (Finite-Difference Time-Domain solver) to compute the
transmittance (T), reflectance (R), and absorbance (A = 1 - T - R) spectra of
a monolayer of hexagonal close-packed (HCP) polystyrene spheres.

Physical & Geometric Parameters:
- Sphere diameter D = 500 nm (0.5 um)
- Sphere center-to-center pitch a = D = 0.5 um
- Lattice symmetry: 2D triangular / hexagonal lattice
- Unit cell: Rectangular supercell representation with
    s_x = a = 0.5 um
    s_y = sqrt(3) * a = 0.866025 um
- Material: Polystyrene (refractive index n = 1.59)
- Background medium: Air (n = 1.0)
- Spectral range: Visible to near-IR (wavelengths 0.35 um to 1.0 um)
"""

import math
import numpy as np

try:
    import meep as mp
except ImportError:
    mp = None


def run_simulation(
    diameter=0.5,        # Sphere diameter in um (500 nm)
    index_sphere=1.59,   # Refractive index of polystyrene
    resolution=50,       # Grid resolution (pixels per um)
    wvl_min=0.35,        # Minimum wavelength in um
    wvl_max=1.00,        # Maximum wavelength in um
    nfreq=100,           # Number of frequency/wavelength points
    dpml=1.0,            # PML boundary layer thickness in um
    dair=1.0,            # Air padding above and below monolayer in um
    output_prefix="hcp_monolayer"
):
    """
    Runs normalization (empty background) and monolayer simulations,
    calculates T(lambda), R(lambda), A(lambda), saves data to CSV and generates a plot.
    """
    if mp is None:
        raise RuntimeError("Meep Python module (pymeep) is required to run this simulation.")

    a = diameter  # Close-packed pitch = sphere diameter
    sx = a
    sy = math.sqrt(3) * a
    sz = 2 * dpml + 2 * dair + diameter

    cell_size = mp.Vector3(sx, sy, sz)

    # Frequency range
    fcen = 0.5 * (1.0 / wvl_min + 1.0 / wvl_max)
    df = (1.0 / wvl_min - 1.0 / wvl_max)

    # Monolayer material
    mat_sphere = mp.Medium(index=index_sphere)

    # PML boundary conditions along z
    boundary_layers = [mp.PML(thickness=dpml, direction=mp.Z)]

    # Source setup: Gaussian plane wave located near bottom PML
    src_z = -0.5 * sz + dpml + 0.2 * dair
    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=df),
            component=mp.Ex,
            center=mp.Vector3(0, 0, src_z),
            size=mp.Vector3(sx, sy, 0),
        )
    ]

    # Monitor locations
    refl_z = -0.5 * sz + dpml + 0.6 * dair
    tran_z = 0.5 * sz - dpml - 0.6 * dair

    # 1. Normalization Run (Empty Cell)
    print("--- Running Normalization Simulation (Empty Cell) ---")
    sim_norm = mp.Simulation(
        cell_size=cell_size,
        resolution=resolution,
        boundary_layers=boundary_layers,
        sources=sources,
        k_point=mp.Vector3(0, 0, 0),
    )

    refl_region = mp.FluxRegion(center=mp.Vector3(0, 0, refl_z), size=mp.Vector3(sx, sy, 0))
    tran_region = mp.FluxRegion(center=mp.Vector3(0, 0, tran_z), size=mp.Vector3(sx, sy, 0))

    refl_flux_norm = sim_norm.add_flux(fcen, df, nfreq, refl_region)
    tran_flux_norm = sim_norm.add_flux(fcen, df, nfreq, tran_region)

    sim_norm.run(
        until_after_sources=mp.stop_when_fields_decayed(20, mp.Ex, mp.Vector3(0, 0, tran_z), 1e-6)
    )

    inc_flux_data = mp.get_fluxes(tran_flux_norm)
    refl_flux_data_saved = sim_norm.get_flux_data(refl_flux_norm)

    # 2. Monolayer Geometry Setup
    # For a rectangular cell (sx = a, sy = sqrt(3)*a), placing 1 sphere at center
    # and 4 spheres at four corners forms an exact HCP monolayer.
    radius = diameter / 2.0
    spheres = [
        mp.Sphere(radius=radius, center=mp.Vector3(0, 0, 0), material=mat_sphere),
        mp.Sphere(radius=radius, center=mp.Vector3(0.5 * sx, 0.5 * sy, 0), material=mat_sphere),
        mp.Sphere(radius=radius, center=mp.Vector3(-0.5 * sx, 0.5 * sy, 0), material=mat_sphere),
        mp.Sphere(radius=radius, center=mp.Vector3(0.5 * sx, -0.5 * sy, 0), material=mat_sphere),
        mp.Sphere(radius=radius, center=mp.Vector3(-0.5 * sx, -0.5 * sy, 0), material=mat_sphere),
    ]

    print("\n--- Running Monolayer Simulation ---")
    sim = mp.Simulation(
        cell_size=cell_size,
        resolution=resolution,
        boundary_layers=boundary_layers,
        geometry=spheres,
        sources=sources,
        k_point=mp.Vector3(0, 0, 0),
    )

    refl_flux = sim.add_flux(fcen, df, nfreq, refl_region)
    tran_flux = sim.add_flux(fcen, df, nfreq, tran_region)

    # Subtract incident flux data to get reflected flux
    sim.load_minus_flux_data(refl_flux, refl_flux_data_saved)

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(20, mp.Ex, mp.Vector3(0, 0, tran_z), 1e-6)
    )

    res_tran = mp.get_fluxes(tran_flux)
    res_refl = mp.get_fluxes(refl_flux)
    freqs = mp.get_flux_freqs(tran_flux)

    # 3. Calculate Transmittance, Reflectance, and Absorbance
    wvls = []
    T = []
    R = []
    A = []

    for i in range(nfreq):
        wvl = 1.0 / freqs[i]
        tran_val = res_tran[i] / inc_flux_data[i]
        refl_val = -res_refl[i] / inc_flux_data[i]
        abs_val = 1.0 - tran_val - refl_val

        wvls.append(wvl)
        T.append(tran_val)
        R.append(refl_val)
        A.append(abs_val)

    # Convert wavelengths to nanometers for presentation
    wvls_nm = [w * 1000.0 for w in wvls]

    # Save to CSV
    csv_file = f"{output_prefix}_spectra.csv"
    with open(csv_file, "w") as f:
        f.write("wavelength_nm,transmittance,reflectance,absorbance\n")
        for i in range(nfreq):
            f.write(f"{wvls_nm[i]:.2f},{T[i]:.6f},{R[i]:.6f},{A[i]:.6f}\n")
    print(f"\nSpectra data saved to {csv_file}")

    # Plot results using matplotlib if available
    try:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 5))
        plt.plot(wvls_nm, T, "b-", label="Transmittance (T)", linewidth=2)
        plt.plot(wvls_nm, R, "r-", label="Reflectance (R)", linewidth=2)
        plt.plot(wvls_nm, A, "g--", label="Absorbance / Loss (1 - T - R)", linewidth=1.5)
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Fraction")
        plt.title("Spectra of HCP Polystyrene Sphere Monolayer (D = 500 nm, n = 1.59)")
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.legend(loc="best")
        plt.tight_layout()
        plot_file = f"{output_prefix}_spectra.png"
        plt.savefig(plot_file, dpi=300)
        print(f"Spectra plot saved to {plot_file}")
    except ImportError:
        print("Matplotlib not installed; skipping plot generation.")


if __name__ == "__main__":
    if mp is not None:
        run_simulation()
    else:
        print("Meep (pymeep) is not installed in the current environment.")
        print("Please install pymeep (e.g. via conda: 'conda install -c conda-forge pymeep') to execute this simulation.")
