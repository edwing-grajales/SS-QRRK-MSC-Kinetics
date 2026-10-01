# SS-QRRK/MSC Pressure-Dependent Rate Constants

Python scripts that compute pressure-dependent rate constants with the system-specific quantum RRK theory combined with the modified strong collision model (SS-QRRK/MSC). Each run reports results for three definitions of the collision efficiency: the original SS-QRRK/MSC, Gilbert, and Dean.

| Script | Reaction type | Results |
|---|---|---|
| `Full_SS-QRRK_Unimolecular_v4.2.py` | Unimolecular (thermally activated): AB → products | k(T, P) |
| `Full_SS-QRRK_Chemical_Activation_v2B3.py` | Chemically activated: A + B ⇌ AB\* → P, with AB\* + M → AB | k<sub>stab</sub>(T, P), k<sub>p</sub>(T, P), k<sub>rev</sub>(T, P) |

## Requirements

Python 3.10 or later, plus:

```bash
pip install numpy scipy pandas openpyxl xlsxwriter lmfit matplotlib
```

Tested with Python 3.10 (NumPy 1.26, pandas 2.2) and Python 3.11 (NumPy 2.4, pandas 3.0).

## Quick start

1. Put the script and the input workbook in the same folder.
2. Set the workbook name near the top of the script:
   ```python
   file = 'my_input.xlsx'
   ```
   The defaults are the example inputs from the cited papers: `jpcatest_wd.xlsx` (unimolecular) and `H-Assisted_SM.xlsx` (chemical activation). In the unimolecular script, `rate` sets the name of the rate-constant column (default `'kf'`).
3. Run the script:
   ```bash
   python Full_SS-QRRK_Unimolecular_v4.2.py
   python Full_SS-QRRK_Chemical_Activation_v2B3.py
   ```
   Running from Spyder or another IDE also works.

The input workbook is only read. All results are written to new files.

## Input workbook

The workbook has three sheets. Empty cells are skipped, so columns can have different lengths.

**`Rates`**: high-pressure-limit (HPL) rate constants

| Column | Unimolecular | Chemical activation |
|---|---|---|
| `T` | Temperature [K] | Temperature [K] |
| `kf` (name set by `rate`) | AB → products [s⁻¹] | – |
| `k1` | – | A + B → AB (bimolecular; any units, which the outputs keep) |
| `km1` | – | AB → A + B [s⁻¹] |
| `k2` | – | AB → P [s⁻¹] |

**`Energy`**: form of the modified-Arrhenius fit, `0` = endothermic, `1` = exothermic. The unimolecular script reads column `energy`; the chemical-activation script reads `energym1` (for `km1`) and `energy2` (for `k2`).

**`SS-QRRK`**: properties of the adduct AB and of the bath gas

| Column | Content |
|---|---|
| `Fr` | Vibrational frequencies of AB [cm⁻¹], real modes only |
| `E_Z` | Zero-point energy of AB [kcal mol⁻¹] |
| `alfa_d` | Energy transferred per deactivating collision at 300 K, α [cm⁻¹] |
| `Td` | Temperature exponent: α(T) = `alfa_d` · (T/300)<sup>`Td`</sup> |
| `PaA` | Lennard-Jones data of AB, one value per row in this order: σ [Å], ε/k<sub>B</sub> [K], molar mass [g mol⁻¹] |
| `PaB` | Bath gas, in this order: σ [Å], ε/k<sub>B</sub> [K], molar mass [g mol⁻¹], T<sub>c</sub> [K], P<sub>c</sub> [bar] |
| `P` | Pressures [atm] |

The scripts always use 200 energy levels. A `MaxI` column in the sheet is ignored.

## Output

**Unimolecular**
- `P-Dependent_<rate>.xlsx`, with sheets `SS-QRRK`, `Gilbert` and `Dean`. Columns: `T (K)`, one column per pressure, and `HPL`.
- `Fit_<rate>.xlsx`: T, A<sub>∞</sub> and E<sub>a</sub> from the fit.
- `<rate>.png`: the fit and its residuals.

**Chemical activation**
- `P-Dependent <input file>.xlsx`: one `k_stab` and one `k_p` sheet for each efficiency model (`SS-QRRK_MSC`, `_G` for Gilbert, `_D` for Dean).
  - Each `k_stab` sheet also has a `Bc` column with the collision efficiency.
  - The `k_m1` sheets (Gilbert and Dean only) give k<sub>rev</sub> = k1 − k<sub>stab</sub> − k<sub>p</sub>.
- `Fit_km1.xlsx` and `Fit_k2.xlsx`: fit results.
- `ln(kₘ₁).png` and `ln(k₂).png`: fits and residuals.

## Method in brief

1. **Arrhenius fit.** Each HPL rate constant is fitted to a modified-Arrhenius form, giving A<sub>∞</sub>(T) and the local activation energy E<sub>a</sub>(T). E<sub>a</sub>(T) is used as the threshold energy E<sub>0</sub>.
2. **F<sub>E</sub>.** F<sub>E</sub> is obtained by numerical integration over the Whitten–Rabinovitch density of states.
3. **Collision efficiency.** The three definitions are:
   - SS-QRRK/MSC: β<sub>c</sub> = [α / (α + F<sub>E</sub>k<sub>B</sub>T)]²
   - Gilbert: β<sub>c</sub>/Δ, with Δ = 1 − exp(−E<sub>0</sub>/F<sub>E</sub>k<sub>B</sub>T)·[1 + E<sub>0</sub>/(α + F<sub>E</sub>k<sub>B</sub>T)]
   - Dean: β<sub>c</sub>/Δ, with Δ = Δ<sub>1</sub> − [F<sub>E</sub>k<sub>B</sub>T/(α + F<sub>E</sub>k<sub>B</sub>T)]·Δ<sub>2</sub>
4. **Collision rate and bath gas.** The Lennard-Jones collision rate uses Troe's Ω<sup>(2,2)\*</sup> fit with a hard-sphere diameter of 2<sup>1/6</sup>σ. The bath-gas concentration comes from the Redlich–Kwong equation of state.
5. **Rate constants.** The QRRK sums over energy levels give k(T, P), or k<sub>stab</sub> and k<sub>p</sub> for chemical activation.

The references give full details.

## What changed in v4.2 and v2B3

The algorithms are the same as in the published versions (v4 and v2B1). There are two kinds of changes.

- **Compatibility.** The scripts run with current pandas (≥ 2.0), NumPy (≥ 1.23, including 2.x) and lmfit (≥ 1.2). The input workbook is no longer modified; earlier versions could overwrite it or empty it.
- **Robust integration of F<sub>E</sub> and the Δ factors.** Earlier versions integrated with `scipy.integrate.quad` from 0 to 20 000 kcal mol⁻¹. At low temperature the integrand is only a few k<sub>B</sub>T wide, so `quad` missed it, returned F<sub>E</sub> ≪ 1, and β<sub>c</sub> came out as 1.
  - The new routine integrates in steps of k<sub>B</sub>T/4 where the integrand is non-negligible.
  - It works in logarithms, which also prevents overflow for adducts with 73 or more vibrational modes.
  - Results change only where the old integration failed: at low temperature, up to about 500 K in the test systems.

Set `ROBUST_INTEGRALS = False` at the top of a script to reproduce the previous version exactly.

## Notes and limitations

- **Chemical activation needs a bound adduct.** E<sub>0</sub> from `km1` must be positive and well above k<sub>B</sub>T at every temperature. As a consistency check, E<sub>a</sub>(k1) − E<sub>a</sub>(km1) should be close to minus the well depth of AB.
- **The QRRK sums use 200 energy levels in double precision,** as in the published versions. For large adducts at high temperature the sums can be truncated. In the unimolecular output, check that the `HPL` column reproduces the input rate constants.
- **The physical constants** are those of the published versions.

## How to cite

If you use these scripts, please cite:

- Unimolecular: E. Grajales-González, M. Monge-Palacios, S. M. Sarathy, "Collision Efficiency Parameter Influence on Pressure-Dependent Rate Constant Calculations Using the SS-QRRK Theory," *J. Phys. Chem. A* **2020**, 124, 6277–6286. https://doi.org/10.1021/acs.jpca.0c02943
- Chemical activation: E. Grajales-González, M. Monge-Palacios, S. M. Sarathy, "A theoretical study of the Ḣ- and HOȮ-assisted propen-2-ol tautomerizations: Reactive systems to evaluate collision efficiency definitions on chemically activated reactions using SS-QRRK theory," *Combust. Flame* **2021**, 225, 485–498. https://doi.org/10.1016/j.combustflame.2020.11.015

The QRRK/MSC formulation for chemical activation follows A. Y. Chang, J. W. Bozzelli, A. M. Dean, "Kinetic Analysis of Complex Chemical Activation and Unimolecular Dissociation Reactions using QRRK Theory and the Modified Strong Collision Approximation," *Z. Phys. Chem.* **2000**, 214, 1533–1568. https://doi.org/10.1524/zpch.2000.214.11.1533

## License

Released under the [MIT License](LICENSE).
