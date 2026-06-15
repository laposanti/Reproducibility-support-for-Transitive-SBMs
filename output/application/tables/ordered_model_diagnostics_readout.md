# Ordered-model diagnostics audit

Run: `application_run_20260414_104327`
EBF prior simulations per dataset: `20000`

## Definitions used here

- Partition backward mass is `sum(A_ij)` over cross-block directed edges with `z_i > z_j`, divided by total cross-block mass for the percentage.
- Empirical WST/SST conformity is a same-K posterior block-level conditional triple rate computed from observed block flow fractions.
- DC-SBM model-implied conformity uses the same conditional triple score on each draw's model-implied directional probability matrix `rho`, after ordering blocks by mean success probability.
- The EBF uses the exact ordered region: WST means all ordered upper-triangle entries are at least 0.5; SST means WST plus the monotone triple inequalities.
- Predictive improvement is judged from LOO. The hierarchy diagnostics say whether the unordered posterior lies in an ordered region; they do not by themselves prove better prediction.

## Comparison with current files

Metrics checked: `15`; mismatching rows: `15`.

Mismatching metrics: bf_sst_0, bf_wst_0, p_post_wst, p_prior_sst, p_prior_wst.

The partition-level backward mass, cross-block mass, empirical block conformity, and DC-SBM model-implied conformity all agree with the current overview. The remaining differences are the corrected EBF quantities.

## Dataset-level readout

| Dataset | DCSBM WST post. | DCSBM SST post. | EBF WST | EBF SST | Ordered vs DC-SBM LOO | Reading |
|---|---:|---:|---:|---:|---|---|
| Bighorn sheep | 27.6% | 0.0% | 4.70 | NA | SST vs DC-SBM: Delta ELPD=43.0, abs t=2.41 | partial WST mass under the DC-SBM posterior and no exact SST mass; ordered predicts clearly better than DC-SBM. |
| Spotted hyenas | 8.2% | 0.0% | 9.16 | NA | WST vs DC-SBM: Delta ELPD=397.8, abs t=4.69 | weak WST mass under the DC-SBM posterior and no exact SST mass; ordered predicts clearly better than DC-SBM. |
| Mountain goats | 100.0% | 3.2% | 4.02 | 3.61 | SST vs DC-SBM: Delta ELPD=37.8, abs t=4.57 | substantial WST mass under the DC-SBM posterior and very little exact SST mass; ordered predicts clearly better than DC-SBM. |
| Stat. journals | 0.0% | 0.0% | NA | NA | WST vs DC-SBM: Delta ELPD=-156.1, abs t=1.47 | no WST mass under the DC-SBM posterior and no exact SST mass; ordered and DC-SBM are predictively close. |
| Japanese macaques | 56.1% | 0.0% | >11220.4 | NA | SST vs DC-SBM: Delta ELPD=-83.4, abs t=2.97 | substantial WST mass under the DC-SBM posterior and no exact SST mass; DC-SBM predicts clearly better than ordered. |
| High school | 0.0% | 0.0% | NA | NA | SST vs DC-SBM: Delta ELPD=150.4, abs t=5.14 | no WST mass under the DC-SBM posterior and no exact SST mass; ordered predicts clearly better than DC-SBM. |

## Overall conclusion

The corrected diagnostics separate two claims that should stay separate. First, a DCSBM posterior can contain ordered structure, measured by posterior mass in the WST/SST regions and the EBF. Second, an ordered model can predict better than the unordered DCSBM, measured by LOO. The first claim is structural; the second is predictive.

The strongest predictive cases for ordered models are the sheep, goats, hyenas, and high school. The citations and macaques do not support a simple 'ordered beats unordered' conclusion: citations are best left to the DC-SBM predictively, and macaques retain a strong directional structure but the DC-SBM predicts better. Exact SST support inside the DCSBM posterior is generally scarce, so Toeplitz SST wins should be read mainly as useful regularisation or shape restriction unless the LOO comparison also clearly separates it from DC-SBM.
