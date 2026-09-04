# 🧩 Lattice-Based ECDSA Private Key Recovery Toolkit

> **A Unified Framework for ECDSA Nonce Weakness Exploitation on secp256k1**

[![SageMath](https://img.shields.io/badge/SageMath-9.0%2B-006400?logo=sagemath&logoColor=white)](https://www.sagemath.org/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Curve](https://img.shields.io/badge/curve-secp256k1-orange)](https://en.bitcoin.it/wiki/Secp256k1)

---

## 📋 Table of Contents

1. [Abstract](#-abstract)
2. [Mathematical Preliminaries](#-mathematical-preliminaries)
3. [Modules Overview](#-modules-overview)
   - [minerva.sage — Centered MSB Bias Attack](#1-minervasage--centered-msb-bias-attack-cve-2019-15809)
   - [polynonce.sage — Polynomial Recurrence Attack](#2-polynoncesage--polynomial-recurrence-attack)
   - [bias_scan.sage — Statistical Detector](#3-bias_scansage--statistical-detector-with-null-hypothesis-testing)
   - [hnp_lattice.sage — Targeted HNP Executor](#4-hnp_latticesage--targeted-hnp-executor)
4. [Integrated Workflow](#-integrated-workflow-detection--execution)
5. [Practical Considerations](#-practical-considerations)
6. [Conclusion](#-conclusion)

---

## 🎯 Abstract

This work presents a comprehensive, modular cryptanalytic toolkit designed to recover **private signing keys** from ECDSA signatures generated with flawed nonces (*k*) on the **secp256k1** elliptic curve — the curve that underpins the Bitcoin protocol.

The framework addresses **three distinct classes** of nonce vulnerabilities:

| Vulnerability Class | Attack Type | Module |
|---|---|---|
| Short bit-length bias | Minerva (CVE-2019-15809) | `minerva.sage` |
| Polynomial recurrence | Kudelski 2023 | `polynonce.sage` |
| Arbitrary static MSB/LSB bias | HNP Lattice | `bias_scan.sage` → `hnp_lattice.sage` |

The system is architecturally divided into a **statistical detection layer** and an **algebraic execution layer**, ensuring computational efficiency when processing massive signature datasets (e.g., full Bitcoin blockchain dumps). All modules are implemented in **SageMath**, leveraging its native support for finite fields, lattice reduction (LLL/BKZ), and polynomial algebra.
---

## 📐 Mathematical Preliminaries

Let **Fₙ** denote the prime field defined by the order of the secp256k1 group:

```
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
```

For each valid ECDSA signature **(rᵢ, sᵢ)** on a message hash **zᵢ**, the signing equation is:

```
sᵢ · kᵢ ≡ zᵢ + rᵢ · d  (mod n)
```

where **kᵢ** is the ephemeral nonce and **d** is the private key. Rearranging yields the affine representation:

```
kᵢ ≡ aᵢ + bᵢ · d  (mod n),   aᵢ = zᵢ · sᵢ⁻¹,   bᵢ = rᵢ · sᵢ⁻¹
```

All four tools exploit weaknesses in the generation of **kᵢ** by reducing the problem to either a **Hidden Number Problem (HNP)** or a **polynomial root-finding problem** over **Fₙ**.

---

## 🧠 Modules Overview

### 1. `minerva.sage` — Centered MSB Bias Attack (CVE-2019-15809)

**Vulnerability class:** The nonce is generated with a reduced bit-length **L** (e.g., 128 or 192 bits) due to timing side-channels or poor RNG implementation.

**Mathematical formulation:** Instead of the standard interval `0 ≤ k < 2ᴸ`, this attack applies a *centering transformation*:

```
kᵢ = 2ˡ⁻¹ + eᵢ,   |eᵢ| < 2ˡ⁻¹ = B
```

Centering halves the error bound, producing significantly shorter lattice vectors than the conventional MSB approach.

**Lattice construction (HNP):** Given **m** signatures, a lattice **Λ** of dimension `m + 2` is constructed with a block-diagonal matrix encoding the relations `eᵢ = aᵢ + bᵢ · d`. The lattice is reduced using **BKZ** (block size β, default 25) to find a short vector whose last coordinate is exactly **±B**, yielding a candidate **d**. The sliding window mechanism (`--max-m 48`) targets localised biases within long signature sequences.

---

### 2. `polynonce.sage` — Polynomial Recurrence Attack

**Vulnerability class:** Consecutive nonces satisfy a deterministic polynomial recurrence of degree **L**:

```
kᵢ₊₁ = c₀ + c₁·kᵢ + c₂·kᵢ² + ⋯ + cₗ·kᵢˡ  (mod n)
```

This pattern is typical of broken PRNGs or incorrectly seeded LCGs.

**Algebraic reduction:** For a window of **L+3** consecutive signatures, each **kᵢ** is expressed as a linear polynomial in **X = d**:

```
kᵢ(X) = aᵢ + bᵢ·X ∈ Fₙ[X]
```

We construct an **(L+2) × (L+2)** matrix **M**, where row **i** is:

```
[1, kᵢ(X), kᵢ²(X), …, kᵢˡ(X), kᵢ₊₁(X)]
```

If the recurrence holds, the determinant **det(M)** must vanish identically in **Fₙ[X]**. The script computes the roots of this polynomial; each root is tested as a candidate private key by verifying `d·G = Qₚᵤᵦ`. The algorithm supports degrees **L = 1** (LCG), **L = 2**, and **L = 3** by default.

---

### 3. `bias_scan.sage` — Statistical Detector with Null Hypothesis Testing

**Vulnerability class:** Arbitrary MSB or LSB bias (nonce has **b** fixed leading/trailing zero bits). Unlike the previous two modules, this script performs **detection without attempting full key recovery**.

**Gaussian Heuristic (GH) baseline:** For a lattice **Λ** of dimension **dim** and determinant **det(Λ)**, the expected length of the shortest vector in a random lattice is:

```
GH = √(dim / (2πe)) · det(Λ)^(1/dim)
```

For each candidate bias **b** and type (MSB/LSB), we construct the HNP lattice and compute the ratio **ρ = λ₁(Λ) / GH**, where **λ₁** is the actual shortest vector length found via LLL. A ratio **ρ < 1** suggests lattice structure.

**Null model (permutation test):** To eliminate false positives caused by lattice scaling artefacts, we generate a *null distribution* by permuting the **zᵢ** values (message hashes). This destroys the algebraic relation `kᵢ = aᵢ + bᵢ·d` while preserving the marginal distributions of **rᵢ, sᵢ, zᵢ**. For each permutation, we recompute **ρ**. The **Z-score** is defined as:

```
Z = (μₙᵤₗₗ − ρᵣₑₐₗ) / σₙᵤₗₗ
```

A signal is classified as **genuine** only when `Z ≥ 3.5` (configurable) and `ρᵣₑₐₗ < 0.8`, ensuring a false-positive rate below **0.02%**.

---

### 4. `hnp_lattice.sage` — Targeted HNP Executor

**Role:** This module is the algebraic execution engine that directly recovers the private key from a window flagged as genuine by `bias_scan.sage`.

#### 1:1 Parameter Compatibility

The lattice construction, bias levels, minimal signature requirement function:

```
m_need(b) = max(4, ⌊(256 / b) · (4 / 3)⌋)
```

and window stepping logic (`step = ⌊m · 0.5⌋`) are implemented **identically** to the detector. This guarantees that any window reported by the detector can be re-targeted exactly using the `--bias`, `--type`, and `--start` parameters.

#### Lattice Construction (HNP Reduction)

We define the bound `B = 2²⁵⁶⁻ᵇ` and a scaling factor `C = max(1, ⌈n/B⌉)`. To solve for **d**, we construct a lattice **Λ** of dimension `dim = m + 2` spanned by the rows of the following matrix **M**:

```
M = | C·n·Iₘ    0        0   |
    | C·b       1        0   |
    | -C·a       0        W  |
```

where:
- **Iₘ** is the `m × m` identity matrix
- **a = (a₁, …, aₘ)** and **b = (b₁, …, bₘ)**
- **W = B** (the final diagonal coordinate)

For the correct private key **d**, there exists a short vector in **Λ** whose last coordinate is exactly **±W**, and whose penultimate coordinate encodes the candidate **d**:

```
v = ((C·n)·εᵢ, 1, W)
```

#### Lattice Reduction & Candidate Extraction

The lattice is reduced using a cascade strategy:

| Condition | Algorithm | Default |
|-----------|-----------|---------|
| `m + 2 ≤ β + 20` | **BKZ**-β | β = 30 |
| Larger dimensions | **LLL** | — |

After reduction, each row vector is evaluated. A row is a valid HNP solution iff its last coordinate **|v_(dim−1)| = W**. For such a row, we compute a candidate private key:

```
d′ = ± v_(dim−2) (mod n)
```

**Cryptographic verification:** `d′·G ≡ Qₚᵤᵦ` confirms successful recovery. This design ensures **zero false positives** — recovery is binary and exact.

---

## 🔄 Integrated Workflow: Detection → Execution

The framework implements a **two-stage pipeline** optimised for massive datasets (e.g., billions of signatures):

```
┌─────────────────────────────────────────────────────────────┐
│                    STAGE 1: DETECTION                       │
│                    bias_scan.sage                           │
│                                                             │
│  Scans all candidate biases, windows, and key files         │
│  Computes GH ratio & null-model statistics (lightweight)    │
│  Outputs JSON report of genuine bias flags                  │
│  Format: (key, bias, type, window_start)                    │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│                    STAGE 2: EXECUTION                        │
│                    hnp_lattice.sage                          │
│                                                             │
│  Takes coordinates from Stage 1                             │
│  Performs BKZ/LLL reduction exclusively on flagged window   │
│  Recovers literal private key (algebraic, costly)           │
└─────────────────────────────────────────────────────────────┘
```

### Resilience & Resumption

All modules support the `--resume` flag, persisting processed filenames (`.progress`) and recovered keys (`.hits.jsonl`). This allows safe interruption and resumption across arbitrary computational environments — a critical feature for long-running cryptanalysis campaigns.
---

## ⚙️ Practical Considerations

### Environment

| Requirement | Details |
|-------------|---------|
| **Runtime** | SageMath 9.0+ |
| **Operations** | Finite-field arithmetic, BKZ/LLL routines, polynomial root finding over Fₙ |
| **Curve** | secp256k1 (hard-coded prime **p**, group order **n**, generator **G**) |
| **Precision** | Native arbitrary-precision arithmetic (no overflow) |

### Input Format

JSON dictionaries mapping public keys (as `"x_y"` strings) to arrays of triples `[r_hex, s_hex, z_hex]`.

### Performance

| Operation | Throughput |
|-----------|------------|
| Detection (`bias_scan`) | ~10–20 MB of signature data / second |
| Execution (`hnp_lattice`) | < 2 seconds per window (dimension ≤ 50) |
| Complexity | O(dim³) for basis reduction |

### Verification

The script does **not** rely on statistical confidence intervals; recovery is **binary** — either the private key satisfies `dG = Q` or it is rejected. This guarantees **zero false positives**.

---

## 🏁 Conclusion

This toolkit provides a rigorous, academically grounded solution for exploiting three major classes of ECDSA nonce weaknesses. By decoupling statistical detection from algebraic key extraction, and by enforcing strict false-positive controls via null-hypothesis testing, the framework delivers both **scalability** and **precision**.

It is suitable for:
- 🔍 Auditing Bitcoin wallet implementations
- 📊 Analysing historical blockchain data for compromised addresses
- 📖 Serving as a reference implementation for future lattice-based cryptanalysis on secp256k1

---

<div align="center">
  <sub>Built with ❤️ for the cryptographic research community</sub>
  <br>
  <sub>Disclaimer: This software is intended for academic research and authorized security audits only.</sub>
</div>
