# Principal Angles Between Arithmetic Weighted-Composition Ranges in H²

This repository contains the source and computational verification materials for the preprint:

**Shad Claiborne**, *Principal Angles Between Arithmetic Weighted-Composition Ranges in H²*.

ORCID: **0009-0005-9443-3679**

## Main result

Let

[
W_n f(z)=(1+z+cdots+z^{n-1})f(z^n),qquad V_n=W_n/sqrt n,
]

and let (mathcal H_n=operatorname{Ran}V_n) in the Hardy space (H^2).
For coprime integers (2le m<n), the paper reduces the compression geometry of
(mathcal H_m) and (mathcal H_n) to an arithmetic weighted path Laplacian.
It derives universal spectral invariants, a sharp extremal theorem, and a closed Hahn-polynomial spectrum when (nequivpm1pmod m).

## Status

Preprint v1.0 prepared for Zenodo publication.

This paper does **not** claim a proof of the Riemann Hypothesis. The results in this repository are unconditional operator-theoretic statements.

## DOI

Zenodo DOI: [10.5281/zenodo.22885535](https://doi.org/10.5281/zenodo.22885535)

## Author

Shad Claiborne — Independent Researcher  
ORCID: 0009-0005-9443-3679

## Files

- `paper/main.tex` — LaTeX manuscript source
- `verification/verify_spectrum.py` — independent finite-dimensional numerical checks
- `notes/computational-verification.md` — verification scope and interpretation

## Reproducibility

The numerical checks are supplementary and are not part of the mathematical proofs.

## Licensing

- Manuscript text and documentation: **CC BY 4.0**
- Verification code: **MIT License**
