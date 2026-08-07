# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Identify fixed-window TF co-binding features of dynamic AP-1 sites.
# Steps: 1) Build per-TF proportion features; 2) fit cross-validation models; 3) export marker and SHAP summaries.
# =======================================
# Source reformatted from snATAC_AgeEx/xgb_window50_tf_contrasts.py.
# Set the variables in config/.env before running this script.
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.paths import publication_path

#!/usr/bin/env python3
"""
TF-identity ("proportional") co-binder model, restricted to a fixed +/-50bp
window around the Fos motif (same window fix validated in
xgb_density_window50bp.py). This redoes the per-TF marker analysis from
xgb_fosonly_contrasts.py (which searches the whole ATAC peak) to check
whether the marker set changes once the peak-width confound is removed by
construction -- in particular whether any TFs emerge as significant POSITIVE
markers of dynamic sites (AUC_proportion > 0.5, qval < 0.05), which the
whole-peak analysis never found (only depletion/stable markers).

Outputs (prefix xgb_window50_tf_):
  xgb_window50_tf_summary_6cell.csv
  xgb_window50_tf_summary_roc_interp.csv / _base_fpr.npy
  xgb_window50_tf_ac_{contrast}_pertf.csv
  xgb_window50_tf_ac_{contrast}_shap.csv
  xgb_window50_tf_ac_all_pertf.csv
  xgb_window50_tf_dynamic_markers.csv   -- significant positive (dynamic) markers, if any
"""
import os, warnings, time
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import xgboost as xgb
import shap
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_curve, roc_auc_score, f1_score
from density_core import ALL_AP1, bh_fdr
from tf_families import assign_family
warnings.filterwarnings('ignore')
os.chdir(publication_path("TOBIAS_COBINDING_DIR"))

FOS_GENE = 'Fos'
FLANK    = 50
RELIABLE = ['Fibroblast','Cardiomyocyte','Coronary_EC','Endocardial_EC',
            'AllMacrophage','Pericytes']
CONTRASTS = {
    'age_gain':      {'cobinder_col':'in_MAS',
                      'pos_rule': lambda df: df['in_MAS'] & ~df['in_YS'],
                      'label':'Age-up'},
    'age_down':      {'cobinder_col':'in_YS',
                      'pos_rule': lambda df: df['in_YS']  & ~df['in_MAS'],
                      'label':'Age-down'},
    'exercise_down': {'cobinder_col':'in_MAS',
                      'pos_rule': lambda df: df['in_MAS'] & ~df['in_MAE'],
                      'label':'Exercise-down'},
    'exercise_up':   {'cobinder_col':'in_MAE',
                      'pos_rule': lambda df: df['in_MAE'] & ~df['in_MAS'],
                      'label':'Exercise-up'},
}
STABLE_RULE = lambda df: df['in_MAS'] & df['in_YS'] & df['in_MAE']

base_fpr = np.linspace(0, 1, 300)
cv5      = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)


def create_labels(fos_ct, pos_rule):
    stable  = fos_ct[STABLE_RULE(fos_ct)].copy(); stable['label']  = 0
    dynamic = fos_ct[pos_rule(fos_ct)].copy();    dynamic['label'] = 1
    return pd.concat([dynamic, stable], ignore_index=True)


def build_prop_features_window(full_df, labeled_df, cell_type, cobinder_col,
                                flank=FLANK):
    """Per-TF co-binder proportions (count / total footprints), restricted to
    a fixed +/-flank bp window around the Fos motif. AP-1 excluded."""
    cell_data = full_df[(full_df['CellType'] == cell_type) &
                        (full_df[cobinder_col] == True)].copy()
    cell_data = cell_data[~cell_data['name'].isin(ALL_AP1)]
    all_tfs   = sorted(cell_data['name'].unique())
    tf_to_idx = {tf: i for i, tf in enumerate(all_tfs)}
    peak_tf_coords = {
        pk: {tf: list(zip(g['start'].values, g['end'].values))
             for tf, g in pdata.groupby('name')}
        for pk, pdata in cell_data.groupby('peak_coord')}

    n_sites = len(labeled_df); n_tfs = len(all_tfs)
    counts = np.zeros((n_sites, n_tfs), dtype=np.int16)
    n_fp   = np.zeros(n_sites, dtype=np.int32)
    labels = np.zeros(n_sites, dtype=np.int32)
    si = 0
    for peak_coord, rows in labeled_df.groupby('peak_coord'):
        ptfs = peak_tf_coords.get(peak_coord, {})
        for _, fr in rows.iterrows():
            fs, fe = fr['start'], fr['end']
            ws, we = fs - flank, fe + flank
            fps = []
            for tf, coords in ptfs.items():
                j = tf_to_idx[tf]
                for ts, te in coords:
                    # must not overlap the Fos footprint, and must overlap
                    # the fixed +/-flank window
                    if (te <= fs or ts >= fe) and (ts < we and te > ws):
                        counts[si, j] += 1
                        fps.append((ts, te))
            if fps:
                fps.sort(); cur_e = fps[0][1]; nreg = 0
                for ts, te in fps[1:]:
                    if ts < cur_e: cur_e = max(cur_e, te)
                    else: nreg += 1; cur_e = te
                n_fp[si] = nreg + 1
            labels[si] = fr['label']; si += 1

    X      = pd.DataFrame(counts[:si], columns=all_tfs)
    X_prop = X.div(np.maximum(n_fp[:si], 1).astype(float), axis=0)
    return X_prop, n_fp[:si], labels[:si]


def make_model(y):
    npos = int(y.sum()); nneg = int(len(y) - npos)
    return xgb.XGBClassifier(
        n_estimators=350, max_depth=5, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8, min_child_weight=3,
        eval_metric='auc', tree_method='hist', n_jobs=-1,
        scale_pos_weight=nneg / max(npos, 1), random_state=42)


print("Loading data..."); t0 = time.time()
df  = pd.read_csv('Tobias_df_Classified_Full.csv.gz', index_col=0)
fos = df[df['name'] == FOS_GENE].copy()
print(f"  {len(df):,} total rows | {len(fos):,} Fos-only rows [{time.time()-t0:.0f}s]\n")

summary_rows = []
interp_rows  = []
all_pertf    = {c: [] for c in CONTRASTS}
all_shap     = {c: [] for c in CONTRASTS}

for contrast, cfg in CONTRASTS.items():
    print(f"=== {cfg['label']} ===")
    cobinder_col = cfg['cobinder_col']

    for ct in RELIABLE:
        fos_ct = fos[fos['CellType'] == ct]
        lab    = create_labels(fos_ct, cfg['pos_rule'])
        if lab['label'].nunique() < 2:
            print(f"  {ct:18s}  skipped (only one class)")
            continue
        t0 = time.time()

        X_prop, n_fp, y = build_prop_features_window(df, lab, ct, cobinder_col)
        dy_mask = y == 1; st_mask = y == 0
        n_dy = int(dy_mask.sum()); n_st = int(st_mask.sum())

        oof         = cross_val_predict(make_model(y), X_prop, y,
                                        cv=cv5, method='predict_proba')[:, 1]
        fpr_, tpr_, th_ = roc_curve(y, oof)
        auc_prop    = roc_auc_score(y, oof)
        ji          = np.argmax(tpr_ - fpr_)
        f1_prop     = f1_score(y, (oof >= th_[ji]).astype(int))

        summary_rows.append({'contrast': contrast, 'model': 'prop_window50', 'CellType': ct,
                             'AUC': auc_prop, 'F1': f1_prop,
                             'n_stable': n_st, 'n_dynamic': n_dy})
        interp_tpr = np.interp(base_fpr, fpr_, tpr_)
        interp_rows.append({'contrast': contrast, 'CellType': ct, 'AUC': auc_prop,
                            **{f'tpr_{i:03d}': interp_tpr[i]
                               for i in range(len(base_fpr))}})

        # ── Per-TF Mann-Whitney AUC (dynamic vs stable proportions) ──────────
        pvals, aucs_tf, med_dy, med_st = [], [], [], []
        for tf in X_prop.columns:
            U, p = mannwhitneyu(X_prop.loc[dy_mask, tf].values,
                                X_prop.loc[st_mask, tf].values,
                                alternative='two-sided')
            pvals.append(p)
            aucs_tf.append(U / (n_dy * n_st))
            med_dy.append(float(np.median(X_prop.loc[dy_mask, tf].values)))
            med_st.append(float(np.median(X_prop.loc[st_mask, tf].values)))
        pvals = np.array(pvals); qvals = bh_fdr(pvals)
        for i, tf in enumerate(X_prop.columns):
            all_pertf[contrast].append({
                'contrast': contrast, 'CellType': ct, 'TF': tf,
                'family': assign_family(tf),
                'AUC_proportion': aucs_tf[i],
                'median_prop_dynamic': med_dy[i],
                'median_prop_stable':  med_st[i],
                'pval': pvals[i], 'qval': qvals[i],
                'significant': qvals[i] < 0.05})

        # ── SHAP (full-data fit) ─────────────────────────────────────────────
        mdl  = make_model(y).fit(X_prop, y)
        sv   = shap.TreeExplainer(mdl).shap_values(X_prop)
        mean_abs = np.abs(sv).mean(axis=0)
        for j, f in enumerate(X_prop.columns):
            if X_prop[f].std() > 0 and np.std(sv[:, j]) > 0:
                d = float(np.sign(np.corrcoef(X_prop[f].values, sv[:, j])[0, 1]))
            else:
                d = 0.0
            all_shap[contrast].append({
                'contrast': contrast, 'CellType': ct,
                'feature': f, 'family': assign_family(f),
                'mean_abs_shap': mean_abs[j],
                'direction_toward': 'dynamic' if d > 0 else 'stable'})

        n_dynamic_sig = sum(1 for i in range(len(X_prop.columns))
                            if qvals[i] < 0.05 and aucs_tf[i] > 0.5)
        print(f"  {ct:18s}  n={n_dy}/{n_st}  AUC={auc_prop:.3f}  F1={f1_prop:.3f}  "
              f"n_TFs={len(X_prop.columns)}  sig_dynamic_markers={n_dynamic_sig}"
              f"  [{time.time()-t0:.0f}s]")
    print()

# ── Save ──────────────────────────────────────────────────────────────────────
summary = pd.DataFrame(summary_rows)
summary.to_csv('xgb_window50_tf_summary_6cell.csv', index=False)
pd.DataFrame(interp_rows).to_csv('xgb_window50_tf_summary_roc_interp.csv', index=False)
np.save('xgb_window50_tf_summary_base_fpr.npy', base_fpr)

for contrast in CONTRASTS:
    pd.DataFrame(all_pertf[contrast]).sort_values(
        ['CellType','AUC_proportion'], ascending=[True,False]).to_csv(
        f'xgb_window50_tf_ac_{contrast}_pertf.csv', index=False)
    pd.DataFrame(all_shap[contrast]).sort_values(
        ['CellType','mean_abs_shap'], ascending=[True,False]).to_csv(
        f'xgb_window50_tf_ac_{contrast}_shap.csv', index=False)

pt_all = pd.concat([pd.DataFrame(all_pertf[c]) for c in CONTRASTS], ignore_index=True)
pt_all.to_csv('xgb_window50_tf_ac_all_pertf.csv', index=False)

print("Saved all xgb_window50_tf_* outputs\n")
print("Mean AUC by contrast (proportional, +/-50bp window):")
print(summary.groupby('contrast')['AUC'].mean().round(3).to_string())

# ── Positive (dynamic) marker check ─────────────────────────────────────────
dyn_hits = pt_all[(pt_all['AUC_proportion'] > 0.5) & (pt_all['qval'] < 0.05)]
dyn_hits.to_csv('xgb_window50_tf_dynamic_markers.csv', index=False)

print(f"\nSignificant (qval<0.05) POSITIVE dynamic markers "
      f"(AUC_proportion>0.5), any contrast/cell type: {len(dyn_hits)} hits")
if len(dyn_hits) > 0:
    agg = (dyn_hits.groupby('TF')
           .agg(n_contrasts=('contrast', 'nunique'),
                n_cells=('CellType', 'count'),
                mean_auc=('AUC_proportion', 'mean'),
                family=('family', 'first'))
           .sort_values('mean_auc', ascending=False))
    print(agg.to_string())
else:
    print("None found -- consistent with whole-peak analysis "
          "(xgb_fosonly_ac_all_pertf.csv): no positive dynamic-site markers.")

stab_hits = pt_all[(pt_all['AUC_proportion'] < 0.5) & (pt_all['qval'] < 0.05)]
print(f"\nFor comparison, significant NEGATIVE (stable) markers: {len(stab_hits)} hits")

print("\n=== Done ===")
