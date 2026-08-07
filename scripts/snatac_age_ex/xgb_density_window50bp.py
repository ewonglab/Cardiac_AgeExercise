# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Train the reported fixed ±50 bp co-binding-density XGBoost classifier.
# Steps: 1) Build fixed-window features; 2) fit stratified cross-validation models; 3) save ROC, SHAP and confound-check outputs.
# =======================================
# Source reformatted from snATAC_AgeEx/xgb_density_window50bp.py.
# Set the variables in config/.env before running this script.
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.paths import publication_path

#!/usr/bin/env python3
"""
Fixed +/-50bp window density model (Fos-only, fixed).

Addresses reviewer concern: the original build_density_features() counts
co-binding TF footprints ANYWHERE within the same ATAC peak, so the search
area for co-binders is literally the peak width (median 747bp, ~76x the
11bp Fos motif width) -- a direct confound. This script instead restricts
the co-binder search to a fixed +/-50bp window around the Fos motif
(111bp total, constant regardless of peak width), removing the confound
by construction rather than by post-hoc statistical matching.

Both the old whole-peak features and the new window-restricted features are
computed on the same sites so a before/after confound check (correlation of
features vs peak_width) can be reported directly.

Outputs (prefix xgb_density_window50_):
  xgb_density_window50_roc.csv
  xgb_density_window50_base_fpr.npy
  xgb_density_window50_shap.csv
  xgb_density_window50_confound_check.csv
"""
import os, warnings, time
import numpy as np
import pandas as pd
import xgboost as xgb
import shap
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve
from density_core import ALL_AP1
warnings.filterwarnings('ignore')
os.chdir(publication_path("TOBIAS_COBINDING_DIR"))

FOS_GENE = 'Fos'
FLANK    = 50
RELIABLE = ['Fibroblast','Cardiomyocyte','Coronary_EC','Endocardial_EC',
            'AllMacrophage','Pericytes']
CONTRASTS = {
    'age_gain':      {'cobinder_col':'in_MAS',
                      'pos_rule': lambda df: df['in_MAS'] & ~df['in_YS']},
    'age_down':      {'cobinder_col':'in_YS',
                      'pos_rule': lambda df: df['in_YS']  & ~df['in_MAS']},
    'exercise_down': {'cobinder_col':'in_MAS',
                      'pos_rule': lambda df: df['in_MAS'] & ~df['in_MAE']},
    'exercise_up':   {'cobinder_col':'in_MAE',
                      'pos_rule': lambda df: df['in_MAE'] & ~df['in_MAS']},
}
STABLE_RULE = lambda df: df['in_MAS'] & df['in_YS'] & df['in_MAE']
MODEL_FEAT_NAMES = ['n_footprints', 'n_bound_motifs']
cv5 = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
BASE_FPR = np.linspace(0, 1, 300)


def _merge_and_count(fps):
    """Given list of (start,end) footprint tuples, return n_footprints
    (merged non-overlapping regions)."""
    if not fps:
        return 0
    fps = sorted(fps); cur_e = fps[0][1]; nreg = 0
    for ts, te in fps[1:]:
        if ts < cur_e: cur_e = max(cur_e, te)
        else: nreg += 1; cur_e = te
    return nreg + 1


def build_density_features_full(full_df, labeled_df, cell_type, cobinder_col):
    """Original whole-peak search (search area = entire ATAC peak)."""
    cell_data = full_df[(full_df['CellType']==cell_type) &
                        (full_df[cobinder_col]==True)].copy()
    cell_data = cell_data[~cell_data['name'].isin(ALL_AP1)]
    peak_tf_coords = {
        pk: {tf: list(zip(g['start'].values, g['end'].values))
             for tf, g in pdata.groupby('name')}
        for pk, pdata in cell_data.groupby('peak_coord')}
    rows_out = []
    for peak_coord, rows in labeled_df.groupby('peak_coord'):
        ptfs = peak_tf_coords.get(peak_coord, {})
        for _, fr in rows.iterrows():
            fs, fe = fr['start'], fr['end']
            fps = []; n_motifs = 0
            for tf, coords in ptfs.items():
                for ts, te in coords:
                    if te <= fs or ts >= fe:
                        fps.append((ts, te)); n_motifs += 1
            rows_out.append({'n_footprints': _merge_and_count(fps),
                             'n_bound_motifs': n_motifs,
                             'peak_width': fr['peak_end']-fr['peak_start'],
                             'label': fr['label']})
    d = pd.DataFrame(rows_out)
    return d[['n_footprints','n_bound_motifs','peak_width']], d['label'].values


def build_density_features_window(full_df, labeled_df, cell_type, cobinder_col,
                                   flank=FLANK):
    """Fixed +/-flank bp window around the Fos motif (search area constant,
    independent of peak width)."""
    cell_data = full_df[(full_df['CellType']==cell_type) &
                        (full_df[cobinder_col]==True)].copy()
    cell_data = cell_data[~cell_data['name'].isin(ALL_AP1)]
    peak_tf_coords = {
        pk: {tf: list(zip(g['start'].values, g['end'].values))
             for tf, g in pdata.groupby('name')}
        for pk, pdata in cell_data.groupby('peak_coord')}
    rows_out = []
    for peak_coord, rows in labeled_df.groupby('peak_coord'):
        ptfs = peak_tf_coords.get(peak_coord, {})
        for _, fr in rows.iterrows():
            fs, fe = fr['start'], fr['end']
            ws, we = fs - flank, fe + flank
            fps = []; n_motifs = 0
            for tf, coords in ptfs.items():
                for ts, te in coords:
                    # exclude footprints overlapping the Fos motif itself,
                    # and require overlap with the fixed +/-flank window
                    if (te <= fs or ts >= fe) and (ts < we and te > ws):
                        fps.append((ts, te)); n_motifs += 1
            rows_out.append({'n_footprints': _merge_and_count(fps),
                             'n_bound_motifs': n_motifs,
                             'peak_width': fr['peak_end']-fr['peak_start'],
                             'label': fr['label']})
    d = pd.DataFrame(rows_out)
    return d[['n_footprints','n_bound_motifs','peak_width']], d['label'].values


def make_model(y):
    npos = int(y.sum()); nneg = int(len(y)-npos)
    return xgb.XGBClassifier(
        n_estimators=350, max_depth=5, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8, min_child_weight=3,
        eval_metric='auc', tree_method='hist', n_jobs=-1,
        scale_pos_weight=nneg/max(npos,1), random_state=42)


def safe_corr(a, b):
    if np.std(a) == 0 or np.std(b) == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


print("Loading data..."); t0 = time.time()
df  = pd.read_csv('Tobias_df_Classified_Full.csv.gz', index_col=0)
fos = df[df['name']==FOS_GENE].copy()
print(f"  {len(df):,} rows [{time.time()-t0:.0f}s]\n")

roc_rows  = []
shap_rows = []
confound_rows = []

for contrast, cfg in CONTRASTS.items():
    print(f"=== {contrast} ===")
    for ct in RELIABLE:
        t0 = time.time()
        fos_ct   = fos[fos['CellType']==ct]
        dyn_rows = fos_ct[cfg['pos_rule'](fos_ct)].copy(); dyn_rows['label']=1
        stb_rows = fos_ct[STABLE_RULE(fos_ct)].copy();    stb_rows['label']=0
        lab = pd.concat([dyn_rows, stb_rows], ignore_index=True)
        if lab['label'].nunique()<2: continue

        X_full,   y = build_density_features_full(df, lab, ct, cfg['cobinder_col'])
        X_window, _ = build_density_features_window(df, lab, ct, cfg['cobinder_col'])

        # ── Confound check: correlation of features vs peak_width, before/after ──
        confound_rows.append({
            'contrast': contrast, 'CellType': ct, 'n': len(y),
            'corr_nfootprints_peakwidth_full':   safe_corr(X_full['n_footprints'].values,   X_full['peak_width'].values),
            'corr_nfootprints_peakwidth_window': safe_corr(X_window['n_footprints'].values, X_window['peak_width'].values),
            'corr_nmotifs_peakwidth_full':       safe_corr(X_full['n_bound_motifs'].values,   X_full['peak_width'].values),
            'corr_nmotifs_peakwidth_window':     safe_corr(X_window['n_bound_motifs'].values, X_window['peak_width'].values),
        })

        # ── Train only on the window-restricted features ──
        X = X_window[MODEL_FEAT_NAMES]
        oof = cross_val_predict(make_model(y), X, y, cv=cv5,
                                method='predict_proba')[:,1]
        auc = roc_auc_score(y, oof)
        fpr, tpr, _ = roc_curve(y, oof)
        tpr_interp = np.interp(BASE_FPR, fpr, tpr)

        row = {'contrast': contrast, 'CellType': ct, 'AUC': auc}
        for i, v in enumerate(tpr_interp):
            row[f'tpr_{i:03d}'] = v
        roc_rows.append(row)

        # SHAP on full-data fit
        mdl = make_model(y).fit(X, y)
        sv  = shap.TreeExplainer(mdl).shap_values(X)
        mean_abs = np.abs(sv).mean(axis=0)
        for j, feat in enumerate(MODEL_FEAT_NAMES):
            if X[feat].std()>0 and np.std(sv[:,j])>0:
                d = float(np.sign(np.corrcoef(X[feat].values, sv[:,j])[0,1]))
            else:
                d = 0.0
            shap_rows.append({
                'contrast': contrast, 'CellType': ct,
                'feature': feat, 'mean_abs_shap': mean_abs[j],
                'direction': 'dynamic' if d>0 else 'stable',
                'signed_shap': mean_abs[j]*(1 if d>0 else -1), 'AUC': auc})

        cr = confound_rows[-1]
        print(f"  {ct:18s}  AUC={auc:.3f}  "
              f"corr(n_fp,pw) full={cr['corr_nfootprints_peakwidth_full']:.3f}->window={cr['corr_nfootprints_peakwidth_window']:.3f}  "
              f"corr(n_motifs,pw) full={cr['corr_nmotifs_peakwidth_full']:.3f}->window={cr['corr_nmotifs_peakwidth_window']:.3f}"
              f"  [{time.time()-t0:.0f}s]")
    print()

roc_df = pd.DataFrame(roc_rows)
roc_df.to_csv('xgb_density_window50_roc.csv', index=False)
np.save('xgb_density_window50_base_fpr.npy', BASE_FPR)

shap_df = pd.DataFrame(shap_rows)
shap_df.to_csv('xgb_density_window50_shap.csv', index=False)

confound_df = pd.DataFrame(confound_rows)
confound_df.to_csv('xgb_density_window50_confound_check.csv', index=False)

print("Summary AUC:")
print(roc_df.groupby('contrast')['AUC'].mean().round(3).to_string())

print("\nMean |SHAP| by feature (averaged over contrasts and cell types):")
feat_summary = shap_df.groupby('feature')[['mean_abs_shap','signed_shap']].mean().round(4)
print(feat_summary.to_string())

print("\nConfound check -- mean correlation with peak_width, full-peak vs +/-50bp window:")
print(confound_df[['corr_nfootprints_peakwidth_full','corr_nfootprints_peakwidth_window',
                   'corr_nmotifs_peakwidth_full','corr_nmotifs_peakwidth_window']]
      .mean().round(3).to_string())

print("\nSaved: xgb_density_window50_roc.csv, xgb_density_window50_base_fpr.npy, "
      "xgb_density_window50_shap.csv, xgb_density_window50_confound_check.csv")
print("=== Done ===")
