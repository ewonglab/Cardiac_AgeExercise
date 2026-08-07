# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Score and remove candidate doublets in AAV dnAP-1 snRNA-seq libraries.
# Steps: 1) Read per-sample matrices; 2) calculate Scrublet scores; 3) write doublet calls and diagnostic figures.
# =======================================
# Source reformatted from dnAP-1_snRNA/Scrublet.ipynb.
# Set the variables in config/.env before running this script.
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.paths import publication_path

import os
import scanpy as sc
import scrublet as scr
import pandas as pd
import concurrent.futures
import matplotlib
matplotlib.use('Agg')  # Forces matplotlib to not use any Xwindows backend
import matplotlib.pyplot as plt

# Define the input structure
base_dir = publication_path("DNAP1_SNRNA_ROOT", "data/CellRanger")
pools = ["Pool_MS", "Pool_SW", "Pool_VW"]

def process_sample_scrublet(args):
    """
    Reads a 10x matrix, runs Scrublet, and returns a DataFrame of doublet scores.
    Accepts a tuple of (pool_name, sample_path).
    """
    # Unpack the variables from the tuple
    pool, sample_path = args 
    
    sample_name = os.path.basename(sample_path)
    matrix_dir = os.path.join(sample_path, "sample_filtered_feature_bc_matrix")
    
    try:
        # Load the Gene Expression matrix
        adata = sc.read_10x_mtx(matrix_dir, var_names='gene_symbols', cache=False)
        
        # Initialize and run Scrublet
        expected_rate = (adata.n_obs / 1000) * 0.008
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_rate)
        
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        
        # Plot the histogram
        fig, axs = scrub.plot_histogram()
        
        # --- NEW CODE: Annotate the threshold on the plots ---
        if scrub.threshold_ is not None:
            for ax in axs:
                # Place the text anchored to the threshold line on the x-axis, 
                # and dynamically near the top of the y-axis
                ax.annotate(f'Thresh: {scrub.threshold_:.3f}', 
                            xy=(scrub.threshold_, 1), 
                            xytext=(5, -15), # Shifts text slightly right and down
                            xycoords=('data', 'axes fraction'),
                            textcoords='offset points',
                            color='red', 
                            fontweight='bold',
                            fontsize=10)
        # -----------------------------------------------------

        plot_filename = f"{pool}_{sample_name}_scrublet_histogram.pdf"
        fig.savefig(plot_filename, bbox_inches='tight')
        plt.close(fig)  # Crucial: closes the figure to free up memory
        
        # Create a clean DataFrame to export
        results_df = pd.DataFrame({
            'Barcode': adata.obs_names,
            'Scrublet_Score': doublet_scores,
            'Scrublet_Predicted_Doublet': predicted_doublets,
            'Sample_ID': sample_name,
            'Pool_ID': pool  # Added to ensure the CSV is trackable
        })
        
        return results_df
        
    except Exception as e:
        print(f"Failed processing {pool} - {sample_name}: {e}")
        return None

def main():
    # 1. Gather all demultiplexed sample directories as tuples: (pool, path)
    sample_tasks = []
    for pool in pools:
        per_sample_dir = os.path.join(base_dir, pool, "outs", "per_sample_outs")
        if os.path.exists(per_sample_dir):
            for sample in os.listdir(per_sample_dir):
                full_path = os.path.join(per_sample_dir, sample)
                if os.path.isdir(full_path):
                    sample_tasks.append((pool, full_path))

    # 2. Run the samples in parallel
    all_results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        # Map the function across all sample tasks (our list of tuples)
        for result_df in executor.map(process_sample_scrublet, sample_tasks):
            if result_df is not None:
                all_results.append(result_df)

    # 3. Concatenate and save
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        out_file = "all_samples_scrublet_scores.csv"
        final_df.to_csv(out_file, index=False)
        print(f"Successfully processed {len(all_results)} samples. Results saved to {out_file}.")
    else:
        print("No results were generated.")

if __name__ == '__main__':
    main()

os.getcwd()
os.chdir(publication_path("DNAP1_SNRNA_ROOT", "Code/"))
