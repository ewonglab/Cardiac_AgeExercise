# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Shared density-controlled co-binding functions (feature matrix, per-TF logistic regression, BH-FDR) imported by the snATAC XGBoost / co-binding scripts.
# Steps: 1) Build per-site co-binder features and merged-footprint density; 2) fit density-controlled logit models; 3) apply Benjamini-Hochberg FDR.
# Provenance: verbatim copy of TOBIAS_Cobinding/density_core.py; only this header was added.
# =======================================
#!/usr/bin/env python3
"""
Shared core for density-controlled cobinding regression.

Extracted verbatim from density_controlled_cobinding.py so multiple analyses
(age/exercise contrasts, exercise-rescue, …) use exactly the same, validated
footprint-region density method. See that script's header for the full rationale.

Key point: density is measured PER FOOTPRINT (merged physical regions), not per
motif name — TOBIAS calls ~7 paralogue names per footprint, so a name count
over-states crowding and is collinear with paralogue-family signals.
"""
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm

# AP-1 genes are excluded from the co-binder set in every analysis.
ALL_AP1 = ['Fos','FosB','c-Fos','Jun','JunB','JunD','c-Jun',
           'FOSL1','FOSL2','JUND','JUNB','Fosb','Fosl1','Fosl2']


def bh_fdr(pvalues):
    n     = len(pvalues)
    order = np.argsort(pvalues)
    ranks = np.argsort(order) + 1
    q     = np.minimum(1.0, pvalues * n / ranks)
    q_sorted = q[order]
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q[order] = q_sorted
    return q


def build_feature_matrix(full_df, labeled_df, cell_type, cobinder_col):
    """
    Returns:
      feat_df    : binary (n_sites × n_TFs) — TF co-bound OUTSIDE the Fos motif (outcome).
      labels     : site class (1 = positive, 0 = control).
      dens_total : (n_sites,) number of distinct MERGED physical footprint regions
                   formed by all co-binder footprints outside the Fos motif.
      excl_mat   : (n_sites × n_TFs) int — merged regions occupied EXCLUSIVELY by that TF,
                   so density_excl(X) = dens_total - excl_mat[:,X] (non-circular exclusion).
    """
    cell_data = full_df[
        (full_df['CellType'] == cell_type) & (full_df[cobinder_col] == True)
    ].copy()
    cell_data = cell_data[~cell_data['name'].isin(ALL_AP1)]
    all_tfs   = sorted(cell_data['name'].unique())
    tf_to_idx = {tf: i for i, tf in enumerate(all_tfs)}
    peak_tf_coords = {
        peak_id: {tf: list(zip(g['start'].values, g['end'].values))
                  for tf, g in peak_data.groupby('name')}
        for peak_id, peak_data in cell_data.groupby('peak_coord')
    }
    n_sites  = len(labeled_df)
    n_tfs    = len(all_tfs)
    features = np.zeros((n_sites, n_tfs), dtype=np.int8)
    excl_mat = np.zeros((n_sites, n_tfs), dtype=np.int16)
    dens_tot = np.zeros(n_sites, dtype=np.int32)
    labels   = np.zeros(n_sites, dtype=np.int32)
    site_idx = 0
    for peak_coord, peak_rows in labeled_df.groupby('peak_coord'):
        peak_tfs = peak_tf_coords.get(peak_coord, {})
        for _, fos_row in peak_rows.iterrows():
            fos_start, fos_end = fos_row['start'], fos_row['end']
            fps = []   # (start, end, tf_idx) for footprints outside the Fos motif
            for tf, coords in peak_tfs.items():
                j = tf_to_idx[tf]; present = False
                for ts, te in coords:
                    if te <= fos_start or ts >= fos_end:
                        fps.append((ts, te, j)); present = True
                if present:
                    features[site_idx, j] = 1
            if fps:
                fps.sort()
                cur_e   = fps[0][1]
                cur_tfs = {fps[0][2]}
                n_reg   = 0
                for ts, te, j in fps[1:]:
                    if ts < cur_e:
                        cur_e = max(cur_e, te); cur_tfs.add(j)
                    else:
                        n_reg += 1
                        if len(cur_tfs) == 1:
                            excl_mat[site_idx, next(iter(cur_tfs))] += 1
                        cur_e = te; cur_tfs = {j}
                n_reg += 1
                if len(cur_tfs) == 1:
                    excl_mat[site_idx, next(iter(cur_tfs))] += 1
                dens_tot[site_idx] = n_reg
            labels[site_idx] = fos_row['label']
            site_idx += 1
    feat_df = pd.DataFrame(features[:site_idx], columns=all_tfs)
    return feat_df, labels[:site_idx], dens_tot[:site_idx], excl_mat[:site_idx]


def fit_one_tf_statsmodels(y, labels, density_excl):
    """logit(TF_present) ~ site_class + density_excl via MLE. Returns (log_OR, SE, p)."""
    X = sm.add_constant(
        np.column_stack([labels.astype(float), density_excl]),
        has_constant='add')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            res = sm.Logit(y, X).fit(disp=0, maxiter=200)
            return res.params[1], res.bse[1], res.pvalues[1]
        except Exception:
            return np.nan, np.nan, np.nan


def compute_density_controlled(feat_df, labels, dens_total, excl_mat):
    """Per TF: Logit with footprint-region density_excl covariate."""
    dens_total = dens_total.astype(float)
    rows = []
    for jcol, tf in enumerate(feat_df.columns):
        y = feat_df[tf].values.astype(int)
        if y.sum() < 5 or (1 - y).sum() < 5:
            rows.append({'TF':tf, 'cond_log_OR':np.nan, 'cond_se':np.nan,
                         'cond_OR':np.nan, 'cond_p':np.nan}); continue
        density_excl = dens_total - excl_mat[:, jcol]
        log_or, se, p = fit_one_tf_statsmodels(y, labels, density_excl)
        rows.append({'TF':tf, 'cond_log_OR':log_or, 'cond_se':se,
                     'cond_OR': np.exp(log_or) if np.isfinite(log_or) else np.nan,
                     'cond_p':p})
    result = pd.DataFrame(rows)
    valid  = result['cond_p'].notna()
    result.loc[valid, 'cond_q'] = bh_fdr(result.loc[valid, 'cond_p'].values)
    result['cond_q']           = result['cond_q'].fillna(np.nan)
    result['cond_significant'] = result['cond_q'].fillna(1) < 0.05
    return result
