# Simulation Guide: Monolayer of Hexagonal Close-Packed (HCP) Polystyrene Spheres

This guide explains how to use Meep (FDTD electromagnetic solver) to simulate the transmittance ($T$), reflectance ($R$), and absorbance ($A$) spectra of a hexagonal close-packed (HCP) monolayer of polystyrene spheres (diameter $D = 500\text{ nm}$, refractive index $n = 1.59$).

---

## 1. Physical Model & Geometry

* **Sphere Diameter ($D$):** $500\text{ nm} = 0.5\ \mu\text{m}$
* **Lattice Pitch ($a$):** $a = D = 0.5\ \mu\text{m}$ (close-packed monolayer)
* **Unit Cell Representation:**
  A 2D triangular/hexagonal lattice with pitch $a$ can be represented using an equivalent rectangular supercell of dimensions:
  $$s_x = a = 0.5\ \mu\text{m}$$
  $$s_y = \sqrt{3} \, a \approx 0.866025\ \mu\text{m}$$
* **Sphere Positions in Rectangular Supercell:**
  - One central sphere at $(0, 0, 0)$
  - Four corner spheres at $(\pm 0.5 s_x, \pm 0.5 s_y, 0)$
* **Boundary Conditions:**
  - Periodic boundary conditions along $x$ and $y$.
  - Perfectly Matched Layers (PML) along $z$ (top and bottom).
* **Source:**
  - Normally incident Gaussian plane wave pulse (linear polarization $E_x$).
* **Spectral Range:**
  - Wavelength range: $350\text{ nm} \le \lambda \le 1000\text{ nm}$.

---

## 2. File Structure

* `hcp_polystyrene_monolayer.py`: Python script to execute the normalization and monolayer simulations.
* `hcp_monolayer_spectra.csv`: Exported transmittance, reflectance, and absorbance data.
* `hcp_monolayer_spectra.png`: Plot of $T(\lambda)$, $R(\lambda)$, and $A(\lambda)$.

---

## 3. How to Run the Simulation

### Prerequisites
Install Meep with Python bindings (`pymeep`). The recommended method is via Conda:
```bash
conda create -n meep -c conda-forge pymeep numpy matplotlib
conda activate meep
```

### Running the Script
Run the Python script directly:
```bash
python3 hcp_polystyrene_monolayer.py
```

To run in parallel using MPI:
```bash
mpirun -np 4 python3 hcp_polystyrene_monolayer.py
```

---

## 4. Customizing Parameters

You can modify function parameters in `hcp_polystyrene_monolayer.py`:

```python
run_simulation(
    diameter=0.5,        # Sphere diameter in um (e.g. 0.5 for 500 nm)
    index_sphere=1.59,   # Refractive index of polystyrene
    resolution=50,       # Grid resolution (pixels / um)
    wvl_min=0.35,        # Minimum wavelength in um
    wvl_max=1.00,        # Maximum wavelength in um
    nfreq=100,           # Number of wavelength sample points
    dpml=1.0,            # PML thickness in um
    dair=1.0,            # Air padding thickness in um
    output_prefix="hcp_monolayer"
)
```

---

## 5. Output Data Format

The CSV output (`hcp_monolayer_spectra.csv`) contains:
- `wavelength_nm`: Wavelength in nanometers.
- `transmittance`: Fraction of incident power transmitted ($T$).
- `reflectance`: Fraction of incident power reflected ($R$).
- `absorbance`: Estimated loss fraction ($A = 1 - T - R$).
