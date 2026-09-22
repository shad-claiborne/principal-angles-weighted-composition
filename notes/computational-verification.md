# Computational verification

The script `verification/verify_spectrum.py` performs finite-dimensional checks of the formulas proved in the manuscript.

It constructs the block-indicator bases directly, forms the cross-Gram matrix independently, and compares its compression Gramian with the arithmetic weighted path Laplacian formula. It also checks the trace and product identities, the Hahn spectrum in the classes `n ≡ ±1 (mod m)`, and the strict extremal inequality outside those classes over a user-selectable finite range.

These calculations are sanity checks only. They are not used as proofs of the theorems.
