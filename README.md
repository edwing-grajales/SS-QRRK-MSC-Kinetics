# SS-QRRK/MSC Pressure-Dependent Rate Constants

These Python scripts calculate pressure-dependent rate coefficients using the System-Specific Quantum Rice–Ramsperger–Kassel (SS-QRRK) theory in combination with the modified strong collision (MSC) model. The strength of this methodology lies in its ability to incorporate tunneling, multi-structure effects, and torsional anharmonicity through the high-pressure limit data. Each run provides results for three definitions of collision efficiency within the MSC model: the original SS-QRRK as described by Bao et al. [1,2], Gilbert et al. [3], and Dean et al. [4]. The last two reduce the deviations of the original definition at elevated temperatures for systems larger than C<sub>2</sub> hydrocarbons, which are particularly relevant to high-temperature hydrocarbon oxidation.

| Script | Reaction type | Results |
|---|---|---|
| `ss_qrrk_unimolecular.py` | Unimolecular (thermally activated): AB → products | *k*(*T*, *P*) |
| `ss_qrrk_chemical_activation.py` | Chemically activated: A + B ⇌ AB\* → P, with AB\* + M → AB | *k*<sub>stab</sub>(*T*, *P*), *k*<sub>p</sub>(*T*, *P*), *k*<sub>rev</sub>(*T*, *P*) |

## Requirements

Python 3.10 or later, plus:

```bash
pip install numpy scipy pandas openpyxl xlsxwriter lmfit matplotlib
```

Tested with Python 3.10 (NumPy 1.26, pandas 2.2) and Python 3.11 (NumPy 2.4, pandas 3.0).

## Quick start

1. Copy the Python script and the workbook with the inputs into the same folder.
2. Specify the workbook name in the script. Each script includes the line:
   ```python
   file = 'my_input.xlsx'
   ```
   The default names are the input examples of the cited papers: `jpcatest_wd.xlsx` (unimolecular) and `H-Assisted_SM.xlsx` (chemical activation).
3. The beginning of each script also loads the workbook, its sheets, and the variables shown in the examples. Modify them as needed.
4. Run the script:
   ```bash
   python ss_qrrk_unimolecular.py
   python ss_qrrk_chemical_activation.py
   ```
   Running from Spyder or another IDE also works.
5. The input workbook is only read. All results are written to new files.

## Input workbook

The workbook has three sheets. Empty cells are skipped, so columns can have different lengths.

- **`Rates`**: high-pressure-limit (HPL) rate constants

| Column | Unimolecular | Chemical activation |
|---|---|---|
| `T` | Temperature [K] | Temperature [K] |
| `kf`&nbsp;(name&nbsp;set&nbsp;by&nbsp;`rate`) | AB&nbsp;→&nbsp;P&nbsp;[s⁻¹] | – |
| `k1` | – | A + B → AB (bimolecular; any units, which the outputs keep) |
| `km1` | – | AB → A + B [s⁻¹] |
| `k2` | – | AB → P [s⁻¹] |

- **`Energy`**: form of the modified-Arrhenius fit, `0` = endothermic, `1` = exothermic. The unimolecular script reads column `energy`; the chemical-activation script reads `energym1` (for `km1`) and `energy2` (for `k2`).

- **`SS-QRRK`**: properties of the adduct AB and of the bath gas

| Column | Content |
|---|---|
| `Fr` | Vibrational frequencies of AB [cm⁻¹], real modes only |
| `E_Z` | Zero-point energy of AB [kcal mol⁻¹] |
| `alfa_d` | Energy transferred per deactivating collision at 300 K, *α* [cm⁻¹] |
| `Td` | Temperature exponent: *α*(*T*) = `alfa_d` · (*T*/300)<sup>`Td`</sup> |
| `PaA` | Lennard-Jones data of AB, one value per row in this order: *σ* [Å], *ε*/*k*<sub>B</sub> [K], molar mass [g mol⁻¹] |
| `PaB` | Bath gas, in this order: *σ* [Å], *ε*/*k*<sub>B</sub> [K], molar mass [g mol⁻¹], *T*<sub>c</sub> [K], *P*<sub>c</sub> [bar] |
| `P` | Pressures [atm] |

The scripts always use 200 energy levels. A `MaxI` column in the sheet is ignored.

## Output

**Unimolecular**
- `P-Dependent_<rate>.xlsx`, with sheets `SS-QRRK`, `Gilbert` and `Dean`. Columns: `T (K)`, one column per pressure, and `HPL`.
- `Fit_<rate>.xlsx`: *T*, *A*<sub>∞</sub> and *E*<sub>a</sub> from the fit.
- `<rate>.png`: the fit and its residuals.

**Chemical activation**
- `P-Dependent <input file>.xlsx`: one `k_stab` and one `k_p` sheet for each efficiency model (`SS-QRRK_MSC`, `_G` for Gilbert, `_D` for Dean).
    - Each `k_stab` sheet also has a `Bc` column with the collision efficiency.
    - The `k_m1` sheets (Gilbert and Dean only) give *k*<sub>rev</sub> = *k*<sub>1</sub> − *k*<sub>stab</sub> − *k*<sub>p</sub>.
- `Fit_km1.xlsx` and `Fit_k2.xlsx`: fit results.
- `ln(kₘ₁).png` and `ln(k₂).png`: fits and residuals.

## Method in brief

The calculation goes from the high-pressure-limit input to the pressure-dependent rate constants in five steps:

1. **Arrhenius fit.** Each HPL rate constant is first fitted to a modified Arrhenius expression. At each temperature, the fit gives the pre-exponential factor *A*<sub>∞</sub>(*T*) and the local activation energy *E*<sub>a</sub>(*T*), which the SS-QRRK method uses as the threshold energy *E*<sub>0</sub> [1,2].
2. **Energy dependence of the density of states.** With *E*<sub>0</sub> determined, the factor *F*<sub>*E*</sub> (the normalized number of states above the threshold) is calculated by numerically integrating the Whitten–Rabinovitch density of states of species AB.
3. **Collision efficiency.** *F*<sub>*E*</sub> and the energy transferred per deactivating collision, *α*(*T*), define the collision efficiency of the MSC model. The original expression [1,2] is *β*<sub>c</sub> = [*α* / (*α* + *F*<sub>*E*</sub>*k*<sub>B</sub>*T*)]². This expression underestimates *β*<sub>c</sub> when *F*<sub>*E*</sub> is large (high temperatures, large molecules), so the two alternatives divide it by a correction factor *Δ*:
    - Gilbert et al. [3]: *Δ* = 1 − exp(−*E*<sub>0</sub>/*F*<sub>*E*</sub>*k*<sub>B</sub>*T*)·[1 + *E*<sub>0</sub>/(*α* + *F*<sub>*E*</sub>*k*<sub>B</sub>*T*)]
    - Dean et al. [4]: *Δ* = *Δ*<sub>1</sub> − [*F*<sub>*E*</sub>*k*<sub>B</sub>*T*/(*α* + *F*<sub>*E*</sub>*k*<sub>B</sub>*T*)]·*Δ*<sub>2</sub>, where *Δ*<sub>1</sub> and *Δ*<sub>2</sub> are integrals of the Boltzmann-weighted density of states below *E*<sub>0</sub>.
4. **Collision rate and bath gas.** The deactivation rate is *β*<sub>c</sub>*k*<sub>LJ</sub>[M]. The Lennard-Jones collision rate constant *k*<sub>LJ</sub> uses Troe's fit to *Ω*<sup>(2,2)\*</sup> with a hard-sphere diameter of 2<sup>1/6</sup>*σ*. The bath-gas concentration [M] comes from the Redlich–Kwong equation of state.
5. **Rate constants.** Finally, the QRRK sums over the energy levels of AB combine these quantities into *k*(*T*, *P*) for unimolecular reactions, or into *k*<sub>stab</sub>, *k*<sub>p</sub> and *k*<sub>rev</sub> for chemically activated reactions.

Check the references for complete details.

## Notes and limitations

- **Chemical activation needs a chemically bound adduct.** The method may not apply to weakly bound intermediate complexes. In terms of the inputs, *E*<sub>0</sub> from `km1` must be positive and well above *k*<sub>B</sub>*T* at every temperature. As a consistency check, *E*<sub>a</sub>(*k*<sub>1</sub>) − *E*<sub>a</sub>(*k*<sub>−1</sub>) should be close to minus the well depth of AB.
- **The QRRK sums use 200 energy levels in double precision,** as in the published versions. For large adducts at high temperature, the sums can be truncated. In the unimolecular output, check that the `HPL` column reproduces the input rate constants.
- **The physical constants** are the same as in the published versions.

## How to cite

If you use these scripts, please cite:

- Unimolecular: E. Grajales-González, M. Monge-Palacios, S. M. Sarathy, "Collision Efficiency Parameter Influence on Pressure-Dependent Rate Constant Calculations Using the SS-QRRK Theory," *J. Phys. Chem. A* **2020**, 124, 6277–6286. https://doi.org/10.1021/acs.jpca.0c02943
- Chemical activation: E. Grajales-González, M. Monge-Palacios, S. M. Sarathy, "A theoretical study of the Ḣ- and HOȮ-assisted propen-2-ol tautomerizations: Reactive systems to evaluate collision efficiency definitions on chemically activated reactions using SS-QRRK theory," *Combust. Flame* **2021**, 225, 485–498. https://doi.org/10.1016/j.combustflame.2020.11.015

## References

1. J. L. Bao, X. Zhang, D. G. Truhlar, "Predicting pressure-dependent unimolecular rate constants using variational transition state theory with multidimensional tunneling combined with system-specific quantum RRK theory: a definitive test for fluoroform dissociation," *Phys. Chem. Chem. Phys.* **2016**, 18, 16659–16670. https://doi.org/10.1039/C6CP02765B
2. J. L. Bao, J. Zheng, D. G. Truhlar, "Kinetics of Hydrogen Radical Reactions with Toluene Including Chemical Activation Theory Employing System-Specific Quantum RRK Theory Calibrated by Variational Transition State Theory," *J. Am. Chem. Soc.* **2016**, 138, 2690–2704. https://doi.org/10.1021/jacs.5b11938
3. R. G. Gilbert, K. Luther, J. Troe, "Theory of Thermal Unimolecular Reactions in the Fall-off Range. II. Weak Collision Rate Constants," *Ber. Bunsenges. Phys. Chem.* **1983**, 87, 169–177. https://doi.org/10.1002/bbpc.19830870218
4. A. Y. Chang, J. W. Bozzelli, A. M. Dean, "Kinetic Analysis of Complex Chemical Activation and Unimolecular Dissociation Reactions using QRRK Theory and the Modified Strong Collision Approximation," *Z. Phys. Chem.* **2000**, 214, 1533–1568. https://doi.org/10.1524/zpch.2000.214.11.1533

## License

Released under the [MIT License](LICENSE).

