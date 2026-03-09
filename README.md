# Hill48 Anisotropic Plasticity UMAT with Voce Hardening

An Abaqus UMAT subroutine implementing **Hill48 anisotropic plasticity** with **Voce isotropic hardening**, using the **Generalized Radial-Return (GRR) algorithm** in the material eigenspace.

## Method

The implementation follows the GRR approach of Versino & Bennett (2018), which:
1. Decomposes the Hill48 anisotropy tensor into its eigenspace via spectral decomposition
2. Transforms the trial deviatoric stress into this eigenspace
3. Performs a scalar radial-return in eigenspace (Newton-Raphson)
4. Maps the corrected stress back to physical space

This avoids the 6x6 matrix inversion per iteration required by standard implicit Hill48 return-mapping algorithms.

## Hardening Law

**Voce isotropic hardening:**

$$\sigma_y(\alpha) = \sigma_{Y0} + Q_\infty \left(1 - e^{-b\alpha}\right)$$

where $\alpha$ is the equivalent plastic strain, $\sigma_{Y0}$ is the initial yield stress, $Q_\infty$ is the saturation stress increment, and $b$ controls the hardening rate.

## Material Properties

| PROPS Index | Parameter | Description |
|:-----------:|-----------|-------------|
| 1 | E | Young's modulus |
| 2 | nu | Poisson's ratio |
| 3 | Sigma_Y0 | Initial yield stress |
| 4 | Q_inf | Voce saturation stress increment |
| 5 | b | Voce hardening rate |
| 6 | F | Hill48 parameter |
| 7 | G | Hill48 parameter |
| 8 | H | Hill48 parameter |
| 9 | L | Hill48 parameter |
| 10 | M | Hill48 parameter |
| 11 | N | Hill48 parameter |

## State Variables

| STATEV Index | Variable |
|:------------:|----------|
| 1 | Equivalent plastic strain (alpha) |
| 2 | Current yield stress |
| 3-8 | Plastic strain tensor (11, 22, 33, 12, 13, 23) |
| 9 | Convergence flag (0 = converged, 1 = failed) |

## Usage in Abaqus

In your Abaqus input file:
```
*MATERIAL, NAME=HILL48_VOCE
*DEPVAR
9
*USER MATERIAL, CONSTANTS=11
** E,    nu,   SY0,  Q_inf, b,     F,    G,    H
 70000., 0.33, 100., 50.,   10.,   0.5,  0.5,  0.5,
** L,    M,    N
 1.5,   1.5,  1.5
```

Compile and link:
```bash
abaqus job=myJob user=Hill48.f
```

## Notes

- 3D only (NTENS=6). For plane stress/strain, a wrapper or modification is needed.
- The tangent `DDSDDE` currently returns the elastic stiffness (or a scaled approximation), not the fully consistent algorithmic tangent. This is sufficient for explicit or modified Newton solvers but may slow convergence with full Newton-Raphson at the global level.
- LAPACK's `DSYEV` routine is required for eigendecomposition. Link against Intel MKL or equivalent.

## Reference

Versino, D. and Bennett, K.C. (2018). *Generalized radial-return mapping algorithm for anisotropic von Mises plasticity framed in material eigenspace.* Int. J. Numer. Meth. Engng, 116: 202-222. [DOI: 10.1002/nme.5921](https://doi.org/10.1002/nme.5921)

## Author

Mohammad Hasaninia
Computational Advanced Manufacturing and Materials Laboratory (CAMML)
Department of Mechanical Engineering, University of Wyoming
