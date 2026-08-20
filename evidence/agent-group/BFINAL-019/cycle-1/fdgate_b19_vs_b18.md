# BFINAL-019 FD gate vs BFINAL-018 (bit-level column comparison)

FD columns (FD_J, FD_gDP, FD_gV) and ADJ_gDP / ADJ_gV columns: string-identical in all 13 rows.
ADJ_J changed by exactly the two new assembly terms. Per-direction (h=0.001):

| Dir | FD_J (b19=b18) | ADJ_J b18 | ADJ_J b19 | delta (new terms) | signJ b18 | signJ b19 | signDP | signV |
|---|---|---|---|---|---|---|---|---|
| D1 | +0.042856 | -0.037489 | -0.034754 | +0.002735 | 0 | 0 | 1 | 1 |
| D2 | -0.003469 | -0.024719 | -0.023814 | +0.000905 | 1 | 1 | 1 | 1 |
| D3 | -0.000283 | -0.012713 | -0.013641 | -0.000929 | 1 | 1 | 1 | 1 |

gDP end-to-end ratios (ADJ/FD, b19): D1 2.130  D2 2.113  D3 2.169 (identical to B18; reproduced bit-identically)

Consistency of the new-terms net effect vs the B18 offline estimate (frozen-lambda state, b18_decomp.json):
- b19 in-pass, same-lambda net:  D1 +0.002735  D2 +0.000905  D3 -0.000929
- b18 lab Prow+Gx (other state): D1 +0.004714  D2 +0.005444  D3 -0.003566
  (lab Prow: +0.006834/+0.003007/-0.011188, Gx: -0.002119/+0.002438/+0.007621)
  -> ALL THREE net-effect signs agree across the two states (D1+/D2+/D3-); magnitudes differ by
     ~2-4x, consistent with the documented state sensitivity of the near-cancelling b_TC functional.
