# Hill48 Anisotropic Plasticity UMAT with Voce Hardening

An Abaqus UMAT subroutine implementing **Hill48 anisotropic plasticity** with **Voce isotropic hardening**, using the **Generalized Radial-Return (GRR) algorithm** in the material eigenspace.

## Method

The implementation follows the GRR approach of Versino & Bennett (2018), which:
1. Decomposes the Hill48 anisotropy tensor into its eigenspace via spectral decomposition
2. Transforms the trial deviatoric stress into this eigenspace
3. Performs a scalar radial-return in eigenspace (Newton-Raphson)
4. Maps the corrected stress back to physical space

This avoids the 6x6 matrix inversion per iteration required by standard implicit Hill48 return-mapping algorithms.

The Hill48 equivalent stress is defined by

$$
\sigma_{\mathrm{eq}}^2 =
F(\sigma_{22}-\sigma_{33})^2 +
G(\sigma_{33}-\sigma_{11})^2 +
H(\sigma_{11}-\sigma_{22})^2 +
2L\sigma_{23}^2 + 2M\sigma_{13}^2 + 2N\sigma_{12}^2.
$$

Yield occurs when the equivalent stress reaches the current Voce yield stress.
Associated flow uses the Hill48 yield surface to determine the plastic strain increment.

## Shear convention

Abaqus supplies physical shear stress components and engineering shear strain increments in the order
`11, 22, 33, 12, 13, 23`. The return map multiplies the trial shear stresses by √2 to work in
Mandel coordinates. It divides the corrected shear stresses by √2 on return and stores plastic
shear strains in engineering form. With all Hill coefficients set to the isotropic values in the
example below, the yield function reduces to von Mises.

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

Save `Hill48.f90` alongside your Abaqus input deck. The file contains the UMAT, return mapping, and tangent calculation. It requires Abaqus/Standard, a supported Fortran compiler, and LAPACK. Submit the job with:
```bash
abaqus job=myJob user=Hill48.f90
```

## Notes

- Full three-dimensional stress only (`NDI=3`, `NSHR=3`, `NTENS=6`). Plane stress and plane strain require a separate implementation.
- The UMAT returns the elastic tangent before yield and the algorithmic GRR tangent for plastic increments.
- If the eigensolve or plastic iteration fails, the UMAT requests a smaller increment through `PNEWDT`. Abaqus/Standard needs automatic incrementation to retry.
- The stored plastic strain components have no explicit `DROT` transport in this implementation. Check material-axis and state-output conventions before using it for large-rotation histories.
- LAPACK's `DSYEV` routine is required for eigendecomposition. Link against Intel MKL or equivalent.

The standalone wrapper changes have been reviewed in source. They have not been run in Abaqus.

## Reference

Versino, D. and Bennett, K.C. (2018). *Generalized radial-return mapping algorithm for anisotropic von Mises plasticity framed in material eigenspace.* Int. J. Numer. Meth. Engng, 116: 202-222. [DOI: 10.1002/nme.5921](https://doi.org/10.1002/nme.5921)

## Author

Mohammad Hasaninia
Computational Advanced Manufacturing and Materials Laboratory (CAMML)
Department of Mechanical Engineering, University of Wyoming
