#!/usr/bin/env python3
"""Numerical checks for the weighted-composition range geometry preprint.

This code verifies finite-dimensional matrix identities. It is not a proof.
"""

from __future__ import annotations

import argparse
import math
import numpy as np


def cross_gram(m: int, n: int) -> np.ndarray:
    """Cross-Gram matrix for normalized m-block and n-block indicators in C^(mn)."""
    N = m * n
    U = np.zeros((N, m), dtype=float)
    V = np.zeros((N, n), dtype=float)
    for j in range(m):
        U[j * n : (j + 1) * n, j] = 1.0 / math.sqrt(n)
    for i in range(n):
        V[i * m : (i + 1) * m, i] = 1.0 / math.sqrt(m)
    return V.T @ U


def arithmetic_weights(m: int, n: int) -> np.ndarray:
    if math.gcd(m, n) != 1:
        raise ValueError("m and n must be coprime")
    return np.array([(j * n % m) * (m - (j * n % m)) for j in range(1, m)], dtype=float)


def weighted_path_laplacian(weights: np.ndarray) -> np.ndarray:
    m = len(weights) + 1
    L = np.zeros((m, m), dtype=float)
    for j, w in enumerate(weights):
        L[j, j] += w
        L[j + 1, j + 1] += w
        L[j, j + 1] -= w
        L[j + 1, j] -= w
    return L


def check_pair(m: int, n: int, tol: float = 1e-10) -> dict[str, float]:
    C = cross_gram(m, n)
    G = C.T @ C
    w = arithmetic_weights(m, n)
    L = weighted_path_laplacian(w)
    formula = np.eye(m) - L / (m * n)

    gram_error = float(np.max(np.abs(G - formula)))
    mus = np.linalg.eigvalsh(L)

    expected_trace = m * (m * m - 1) / 3
    trace_error = abs(float(np.trace(L)) - expected_trace)

    nonzero = mus[1:]
    expected_pseudodet = m * math.factorial(m - 1) ** 2
    pseudodet_error = abs(float(np.prod(nonzero)) - expected_pseudodet) / expected_pseudodet

    extremal = (n % m) in (1, m - 1)
    mu1 = float(mus[1])
    if extremal:
        hahn = np.array([k * (k + 1) for k in range(m)], dtype=float)
        hahn_error = float(np.max(np.abs(mus - hahn)))
        extremal_error = max(0.0, abs(mu1 - 2.0) - tol)
    else:
        hahn_error = 0.0
        extremal_error = max(0.0, mu1 - 2.0 - tol)

    return {
        "gram_error": gram_error,
        "trace_error": trace_error,
        "pseudodet_rel_error": pseudodet_error,
        "hahn_error": hahn_error,
        "mu1": mu1,
        "extremal_violation": extremal_error,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-m", type=int, default=20)
    parser.add_argument("--max-n", type=int, default=80)
    parser.add_argument("--tol", type=float, default=1e-9)
    args = parser.parse_args()

    worst = {"gram_error": 0.0, "trace_error": 0.0, "pseudodet_rel_error": 0.0, "hahn_error": 0.0}
    strict_gap = float("inf")
    pairs = 0

    for m in range(2, args.max_m + 1):
        for n in range(m + 1, args.max_n + 1):
            if math.gcd(m, n) != 1:
                continue
            pairs += 1
            out = check_pair(m, n, args.tol)
            for key in worst:
                worst[key] = max(worst[key], out[key])
            if n % m not in (1, m - 1):
                strict_gap = min(strict_gap, 2.0 - out["mu1"])
            if out["gram_error"] > args.tol or out["hahn_error"] > args.tol or out["extremal_violation"] > 0:
                raise AssertionError((m, n, out))

    print(f"checked coprime pairs: {pairs}")
    print(f"max Gramian error: {worst['gram_error']:.3e}")
    print(f"max trace error: {worst['trace_error']:.3e}")
    print(f"max pseudodeterminant relative error: {worst['pseudodet_rel_error']:.3e}")
    print(f"max Hahn-spectrum error: {worst['hahn_error']:.3e}")
    if strict_gap < float('inf'):
        print(f"smallest observed strict gap 2-mu1 outside ±1 classes: {strict_gap:.6g}")


if __name__ == "__main__":
    main()
