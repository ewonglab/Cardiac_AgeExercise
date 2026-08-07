# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Plot the reported ROC curves for fixed-window co-binding density models.
# Steps: 1) Load cross-validation ROC tables; 2) aggregate curves by contrast; 3) write the final classifier figure.
# =======================================
# Source reformatted from snATAC_AgeEx/xgb_density_window50_contrasts_roc.py.
# Set the variables in config/.env before running this script.
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.paths import publication_path

#!/usr/bin/env python3
"""
+/-50bp window density model — 4 ROC curves, one per contrast (dynamic set),
mean +/- SD across 6 cell types. Style matches xgb_results_figure.py.
"""
import os
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
os.chdir(publication_path("TOBIAS_COBINDING_DIR"))

BASE_FPR = np.linspace(0, 1, 300)
tpr_cols = [f'tpr_{i:03d}' for i in range(300)]

CONTRAST_ORDER = ['age_gain','age_down','exercise_down','exercise_up']
CONTRAST_LBLS  = {'age_gain':'Age gain','age_down':'Age loss',
                  'exercise_down':'Exercise loss','exercise_up':'Exercise gain'}
CONTRAST_COLORS = {'age_gain':'#C62828','age_down':'#1565C0',
                   'exercise_down':'#2E7D32','exercise_up':'#E65100'}

plt.rcParams.update({'pdf.fonttype':42, 'ps.fonttype':42,
                     'font.family':'sans-serif', 'font.size':6,
                     'axes.labelsize':6, 'axes.titlesize':7,
                     'xtick.labelsize':5.5, 'ytick.labelsize':5.5})

df = pd.read_csv('xgb_density_window50_roc.csv')

fig, ax = plt.subplots(1, 1, figsize=(2.6, 2.5))
fig.subplots_adjust(left=0.20, right=0.95, top=0.85, bottom=0.18)

for contrast in CONTRAST_ORDER:
    grp = df[df['contrast'] == contrast]
    if grp.empty:
        continue
    tprs = grp[tpr_cols].values
    mean_tpr = tprs.mean(0); std_tpr = tprs.std(0)
    mean_auc = grp['AUC'].mean()
    col = CONTRAST_COLORS[contrast]
    lbl = f"{CONTRAST_LBLS[contrast]} (AUC={mean_auc:.2f})"
    ax.plot(BASE_FPR, mean_tpr, color=col, lw=1.2, label=lbl)
    ax.fill_between(BASE_FPR,
                    np.clip(mean_tpr - std_tpr, 0, 1),
                    np.clip(mean_tpr + std_tpr, 0, 1),
                    alpha=0.10, color=col)

ax.plot([0,1],[0,1], 'k--', lw=0.6, alpha=0.4)
ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.set_xlabel('False positive rate', fontsize=6)
ax.set_ylabel('True positive rate', fontsize=6)
ax.set_title('+/-50bp window model\n(mean over 6 cell types)',
             fontsize=7, fontweight='bold', loc='left', pad=4)
ax.legend(fontsize=4.5, loc='lower right', framealpha=0.8,
          handlelength=1.0, handletextpad=0.4)
ax.spines[['top','right']].set_visible(False)

plt.savefig('xgb_density_window50_contrasts_roc.pdf', bbox_inches='tight')
plt.savefig('xgb_density_window50_contrasts_roc.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved xgb_density_window50_contrasts_roc.pdf/.png")
