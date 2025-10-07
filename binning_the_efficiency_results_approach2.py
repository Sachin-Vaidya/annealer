import sys
import os
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba, to_hex
from scipy.optimize import curve_fit
from PIL import Image

def fit_and_plot_logistic(filename_prefix, data_file, x_col, y_col, y_col_list_logistic_fit, y_col_binned_logistic_fit, y_col_list_mean, y_col_list_no_logistic_Valid_CorrectValid, output_png, color, x_label='X', y_label='Y', reuse_bin_edges=None, title="Binned Data and Logistic Fit", label="Binned means", ymin=None, ymax=None, multiply_cols=None):
    ymin, ymax, axis_scale = axes_features(filename_prefix, y_col)
    
    df = pd.read_csv(data_file, sep='\t', header=None)
    
    x = df[x_col].values
    if multiply_cols is not None:
        # Calculate product of specified columns for y
        y = df[multiply_cols[0]].values * df[multiply_cols[1]].values
    else:
        y = df[y_col].values
    
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    y = y[sort_idx]
    
    # Define logistic function
    def logistic(x, L, k, x0):
        return L / (1 + np.exp(-k * (x - x0)))
    
    # Initialize variables for binning and fitting
    best_sse = np.inf
    best_num_bins = None
    best_popt = None
    best_bin_mids = None
    best_bin_means = None
    best_bin_stderr = None
    best_bin_edges = None
    best_bin_edges_left = None
    best_bin_edges_right = None
    
    
    if reuse_bin_edges is None:
        # Iterate over possible number of bins to find best fit
        possible_bins = range(6, 31)
        for num_bins in possible_bins:
            try:
                # Quantile binning (equal count per bin, handles skew)
                bin_labels, bin_edges = pd.qcut(x, num_bins, labels=False, retbins=True, duplicates='drop')
                actual_bins = len(bin_edges) - 1
                if actual_bins < 4:  # Need enough points for meaningful fit
                    continue
        
                bin_mids = []
                bin_means = []
                bin_stderr = []
                for i in range(actual_bins):
                    mask = (bin_labels == i)
                    if mask.sum() > 0:
                        bin_mids.append(np.mean(x[mask]))
                        bin_means.append(np.mean(y[mask]))
                        bin_stderr.append(np.std(y[mask]) / np.sqrt(mask.sum()))  # Standard error
                bin_mids = np.array(bin_mids)
                bin_means = np.array(bin_means)
                bin_stderr = np.array(bin_stderr)
                bin_edges_left = bin_edges[:-1]
                bin_edges_right = bin_edges[1:]
                
                if y_col not in y_col_list_logistic_fit or multiply_cols is not None: #y_col_list_logistic_fit = [6, 7, 9] for SA and DA
                    try:
                        # Initial guess for parameters
                        p0 = [0.5, 0.1, np.median(bin_mids)]
                    
                        # Fit
                        bounds = ([0, 0, min(x)], [1, 50, max(x)])  # L>0, k>0, x0 within bin_mids
                        popt, _ = curve_fit(logistic, bin_mids, bin_means, p0=p0, bounds=bounds, maxfev=10000, ftol=1e-6, xtol=1e-6, gtol=1e-6)
                    
                        # Compute SSE
                        fitted = logistic(bin_mids, *popt)
                        sse = np.sum((bin_means - fitted)**2)
                    except:
                        popt = None
                        sse = np.inf
                '''else:
                    # For y_col=5 or 6, compute SSE based on mean line
                    mean_y = np.mean(y)
                    sse = np.sum((bin_means - mean_y)**2)
                    popt = None'''
                
                if sse < best_sse:
                    best_sse = sse
                    best_num_bins = num_bins
                    best_popt = popt
                    best_bin_mids = bin_mids
                    best_bin_means = bin_means
                    best_bin_stderr = bin_stderr
                    best_bin_edges_left = bin_edges_left
                    best_bin_edges_right = bin_edges_right
                    best_bin_edges = bin_edges
            except:
                pass  # Skip if fit fails
    else:
        # Use provided bin edges
        bin_edges = reuse_bin_edges
        bin_labels = pd.cut(x, bins=bin_edges, labels=False, include_lowest=True).astype(float)
        actual_bins = len(bin_edges) - 1
        bin_mids = []
        bin_means = []
        bin_stderr = []
        for i in range(actual_bins):
            mask = (bin_labels == i)
            if mask.sum() > 0:
                bin_mids.append(np.mean(x[mask]))
                bin_means.append(np.mean(y[mask]))
                if multiply_cols is not None:
                    # Calculate error for product: sqrt((col6 * err_col5)^2 + (col5 * err_col6)^2)
                    col5_vals = df[multiply_cols[0]].values[sort_idx][mask]
                    col6_vals = df[multiply_cols[1]].values[sort_idx][mask]
                    col5_mean = np.mean(col5_vals)
                    col6_mean = np.mean(col6_vals)
                    col5_stderr = np.std(col5_vals) / np.sqrt(mask.sum())
                    col6_stderr = np.std(col6_vals) / np.sqrt(mask.sum())
                    stderr = np.sqrt((col6_mean * col5_stderr)**2 + (col5_mean * col6_stderr)**2)
                    bin_stderr.append(stderr)
                else:
                    bin_stderr.append(np.std(y[mask]) / np.sqrt(mask.sum()))  # Standard error
                #bin_stderr.append(np.std(y[mask]) / np.sqrt(mask.sum()))  # Standard error
        bin_mids = np.array(bin_mids)
        bin_means = np.array(bin_means)
        bin_stderr = np.array(bin_stderr)
        bin_edges_left = bin_edges[:-1]
        bin_edges_right = bin_edges[1:]
        best_num_bins = actual_bins
        best_bin_mids = bin_mids
        best_bin_means = bin_means
        best_bin_stderr = bin_stderr
        best_bin_edges_left = bin_edges_left
        best_bin_edges_right = bin_edges_right
        best_bin_edges = bin_edges
        # Perform logistic fit only if y_col is not 5 or 6, or if multiply_cols is specified
        if y_col==y_col_binned_logistic_fit or multiply_cols is not None:  #y_col_binned_logistic_fit = 6 for SA correct|valid
            try:
                p0 = [0.5, 0.1, np.median(bin_mids)]
                bounds = ([0, 0, min(x)], [1, 50, max(x)])
                best_popt, _ = curve_fit(logistic, bin_mids, bin_means, p0=p0, bounds=bounds, maxfev=10000, ftol=1e-6, xtol=1e-6, gtol=1e-6)
                fitted = logistic(bin_mids, *best_popt)
                best_sse = np.sum((bin_means - fitted)**2)
            except:
                best_popt = None
                best_sse = np.inf
        elif y_col in y_col_list_mean: #y_col_list_mean [7, 9] for SA and DA anneal time
            mean_y = np.mean(y)
            best_sse = np.sum((bin_means - mean_y)**2)
            best_popt = None

    if best_popt is not None:
        # Monte Carlo error estimation for L and x0
        N = 10000
        L_samples = []
        k_samples = []
        x0_samples = []
        for _ in range(N):
            y_sample = np.random.normal(best_bin_means, best_bin_stderr)
            try:
                popt_sample, _ = curve_fit(logistic, best_bin_mids, y_sample, p0=best_popt, bounds=bounds, maxfev=10000, ftol=1e-6, xtol=1e-6, gtol=1e-6)
                L_samples.append(popt_sample[0])
                k_samples.append(popt_sample[1])
                x0_samples.append(popt_sample[2])
            except:
                pass
        std_L = np.std(L_samples) if L_samples else np.nan
        std_k = np.std(k_samples) if k_samples else np.nan
        std_x0 = np.std(x0_samples) if x0_samples else np.nan
    else:
        y_mean = np.mean(y)
    
    print(f"Best number of bins: {best_num_bins}")
    if best_popt is not None:
        print(f"Best SSE: {best_sse}")
        print(f"Fitted parameters: L={best_popt[0]:.16f} ± {std_L:.16f}, k={best_popt[1]:.16f} ± {std_k:.16f}, x0={best_popt[2]:.16f} ± {std_x0:.16f}")
    else:
        print(f"Best SSE (mean line): {best_sse}")
    
    # Print raw data points for the best binning configuration
    if best_num_bins is not None:  # Ensure a valid binning was found
        print(f"\nRaw data points for best binning ({best_num_bins} bins):")
        # Recompute bin labels for the best number of bins
        bin_labels, bin_edges = pd.qcut(x, best_num_bins, labels=False, retbins=True, duplicates='drop') if reuse_bin_edges is None else (pd.cut(x, bins=bin_edges, labels=False, include_lowest=True).astype(float), bin_edges)
        actual_bins = len(bin_edges) - 1
        for i in range(actual_bins):
            mask = (bin_labels == i)
            bin_x_values = x[mask]
            bin_y_values = y[mask]
            bin_mean = np.mean(bin_y_values) if len(bin_y_values) > 0 else np.nan
            if len(bin_x_values) == len(bin_y_values):
                print(f"Bin {i+1} (edges: [{bin_edges[i]:.16f}, {bin_edges[i+1]:.16f}]): {np.mean(bin_x_values):.16f}, {bin_mean:.16f}")
            else:
                print(f"  Error: Mismatched lengths (x: {len(bin_x_values)}, y: {len(bin_y_values)})")
        print(f"\nBin Centers and Mean Values for Best Binning ({best_num_bins} bins):")
        for i, (center, mean) in enumerate(zip(best_bin_mids, best_bin_means)):
            print(f"Bin {i+1}: Center = {center:.16f}, Mean y = {mean:.16f}")
    else:
        print("No valid binning configuration found.")
    
        # Create fainter color by reducing alpha
    base_color = to_rgba(color)
    fainter_color = (base_color[0], base_color[1], base_color[2], 0.6)  # Set alpha to 0.6 for fainter effect
    fainter_color_hex = to_hex(fainter_color, keep_alpha=True)
    
    # Plot
    plt.figure(figsize=(10, 10))
    # Plot binned means with error bars
    plt.errorbar(
        best_bin_mids, 
        best_bin_means, 
        yerr=best_bin_stderr, 
        xerr=[best_bin_mids - best_bin_edges_left, best_bin_edges_right - best_bin_mids],
        fmt='^', 
        color=color, 
        ecolor=color, 
        capsize=3, 
        label=label
    )
    # Plot fitted curve or mean line, but skip for y_col=5 or 6 unless multiply_cols is specified
    if y_col not in y_col_list_no_logistic_Valid_CorrectValid or multiply_cols is not None:  #y_col_list_no_logistic_Valid_CorrectValid = [5, 6] meaning don't fit logistic for SA valid and correct|valid
        x_plot = np.linspace(min(x), max(x), 10000)
        if 'logistic' in title.lower() and best_popt is not None:
            # Construct the label for the fitted curve
            fit_label = f"{best_popt[0]:.4f}/(1+exp(-{best_popt[1]:.4f}(x - {best_popt[2]:.4f})))"
            plt.plot(x_plot, logistic(x_plot, *best_popt), color=fainter_color_hex, linewidth=2, linestyle='--', label=None)
        else:
            mean_y = np.mean(y)  # Use raw data mean
            plt.axhline(y=mean_y, color=fainter_color_hex, linewidth=2, linestyle='--', label=f'Binned Mean y = {mean_y:.4f}')
    
    plt.xlabel(x_label, loc='right', fontsize=20)
    plt.ylabel(y_label, loc='top', fontsize=20)
    title = plt.title(title)
    title.set_visible(False)
    plt.legend(fontsize=16, frameon=False)
    plt.xscale('log')
    plt.yscale(axis_scale)
    plt.xlim(0.5, x.max())
    plt.ylim(ymin, ymax)
    if best_popt is not None:
        x0 = best_popt[2]
        L = best_popt[0]
        plt.vlines(x0, 0, 0.6 * L, color=color, linestyle='dotted', alpha=0.7)
        plt.plot(x0, ymin, marker='v', color=color, markersize=8, fillstyle='none')
    plt.grid(False)
    #plt.xticks(fontsize=20)
    #plt.yticks(fontsize=20)
    plt.tick_params(axis='x', which='both', labelsize=20)
    plt.tick_params(axis='y', which='both', labelsize=20)
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    if best_popt is not None:
        return best_bin_edges, best_popt[0], std_L, best_popt[1], std_k, best_popt[2], std_x0, np.nan, best_bin_mids, best_bin_means, best_bin_stderr, best_bin_edges_left, best_bin_edges_right, color, label, fainter_color_hex
    else:
        return best_bin_edges, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, y_mean, best_bin_mids, best_bin_means, best_bin_stderr, best_bin_edges_left, best_bin_edges_right, color, label, fainter_color_hex

def plot_combined_correct_dunn_index(filename_prefix, data_file, x_col, y_col,
    bin_mids_SA, bin_means_SA, bin_stderr_SA, bin_edges_left_SA, bin_edges_right_SA, color_SA, label_SA, fainter_color_SA, popt_SA,
    x_col_SA=None,  # New: override x_col for SA data
    bin_mids_DA=None, bin_means_DA=None, bin_stderr_DA=None, bin_edges_left_DA=None, bin_edges_right_DA=None, color_DA=None, label_DA=None, fainter_color_DA=None, popt_DA=None,
    x_col_DA=None,  # New: override x_col for DA data
    bin_mids_QA=None, bin_means_QA=None, bin_stderr_QA=None, bin_edges_left_QA=None, bin_edges_right_QA=None, 
    color_QA=None, label_QA=None, fainter_color_QA=None, popt_QA=None,
    x_col_QA=None,  # New: override x_col for QA data
    output_png='Convergence_efficiency_correct_dunn_index.png', x_label='Dunn Index', y_label='Convergence Efficiency', 
    title='Binned Data and Logistic Fit'):
    
    ymin, ymax, axis_scale = axes_features(filename_prefix, y_col)
    
    # Load df once for x data
    df = pd.read_csv(data_file, sep='\t', header=None)
    
    # Compute per-dataset x_max (fallback to shared x_col if not overridden)
    x_max_SA = df[x_col_SA if x_col_SA is not None else x_col].max()
    x_max_DA = df[x_col_DA].max() if (bin_mids_DA is not None and x_col_DA is not None) else None
    x_max_QA = df[x_col_QA].max() if (bin_mids_QA is not None and x_col_QA is not None) else None
    
    # Overall x_max for shared axis/curves
    all_x_maxes = [x_max_SA]
    if x_max_DA is not None:
        all_x_maxes.append(x_max_DA)
    if x_max_QA is not None:
        all_x_maxes.append(x_max_QA)
    overall_x_max = max(all_x_maxes)
    
    # For curves: use overall min/max (fallback to SA if single dataset)
    x_data_min = min(df[x_col_SA if x_col_SA is not None else x_col].min(), 
                     df[x_col_QA].min() if x_max_QA is not None else float('inf'))
    x_data_min = min(x_data_min, df[x_col_DA].min() if x_max_DA is not None else float('inf'))
    x_plot = np.linspace(x_data_min, overall_x_max, 10000)
    
    plt.figure(figsize=(10, 10))
    
    # Plot SA binned data (always present)
    plt.errorbar(
        bin_mids_SA,
        bin_means_SA,
        yerr=bin_stderr_SA,
        xerr=[bin_mids_SA - bin_edges_left_SA, bin_edges_right_SA - bin_mids_SA],
        fmt='^',
        color=color_SA,
        ecolor=color_SA,
        capsize=3,
        label=label_SA
    )
    
    # Plot DA binned data (conditional)
    if bin_mids_DA is not None:
        plt.errorbar(
            bin_mids_DA,
            bin_means_DA,
            yerr=bin_stderr_DA,
            xerr=[bin_mids_DA - bin_edges_left_DA, bin_edges_right_DA - bin_mids_DA],
            fmt='o',
            color=color_DA,
            ecolor=color_DA,
            capsize=3,
            label=label_DA
        )
    
    # Plot QA binned data (conditional)
    if bin_mids_QA is not None:
        plt.errorbar(
            bin_mids_QA,
            bin_means_QA,
            yerr=bin_stderr_QA,
            xerr=[bin_mids_QA - bin_edges_left_QA, bin_edges_right_QA - bin_mids_QA],
            fmt='v',
            color=color_QA,
            ecolor=color_QA,
            capsize=3,
            label=label_QA
        )
    
    # Plot fitted curves
    def logistic(x, L, k, x0):
        return L / (1 + np.exp(-k * (x - x0)))
    
    if popt_SA is not None:
        plt.plot(x_plot, logistic(x_plot, *popt_SA), color=fainter_color_SA, linewidth=2, linestyle='--', label=None)
    if popt_DA is not None:
        plt.plot(x_plot, logistic(x_plot, *popt_DA), color=fainter_color_DA, linewidth=2, linestyle='--', label=None)
    if popt_QA is not None:
        plt.plot(x_plot, logistic(x_plot, *popt_QA), color=fainter_color_QA, linewidth=2, linestyle='--', label=None)
        
    if popt_SA is not None:
        plt.vlines(popt_SA[2], 0, 0.6 * popt_SA[0], color=color_SA, linestyle='dotted', alpha=0.7)
        plt.plot(popt_SA[2], ymin, marker='v', color=color_SA, markersize=10, fillstyle='none')
    if popt_DA is not None:
        plt.vlines(popt_DA[2], 0, 0.6 * popt_DA[0], color=color_DA, linestyle='dotted', alpha=0.7)
        plt.plot(popt_DA[2], ymin, marker='v', color=color_DA, markersize=10, fillstyle='none')
    if popt_QA is not None:
        plt.vlines(popt_QA[2], 0, 0.6 * popt_QA[0], color=color_QA, linestyle='dotted', alpha=0.7)
        plt.plot(popt_QA[2], ymin, marker='v', color=color_QA, markersize=10, fillstyle='none')
    
    plt.xlabel(x_label, loc='right', fontsize=20)
    plt.ylabel(y_label, loc='top', fontsize=20)
    title_obj = plt.title(title)
    title_obj.set_visible(False)
    handles, labels = plt.gca().get_legend_handles_labels()
    new_labels = [filename_prefix.replace('_', ', ').replace('V', ' V').replace('T', ' T')] + labels
    new_handles = [plt.Line2D([0], [0], color='none', linewidth=0)] + handles
    plt.legend(new_handles, new_labels, fontsize=16, frameon=False)
    plt.xscale('log')
    plt.yscale(axis_scale)
    plt.xlim(0.5, overall_x_max)  # Updated to use overall
    plt.ylim(ymin, ymax)
    plt.grid(False)
    #plt.xticks(fontsize=20)
    #plt.yticks(fontsize=20)
    plt.tick_params(axis='x', which='both', labelsize=20)
    plt.tick_params(axis='y', which='both', labelsize=20)
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    
def axes_features(filename_prefix, y_col):
    if filename_prefix in ["2Vertices_10Tracks", "2Vertices_16Tracks", "2Vertices_28Tracks", "3Vertices_15Tracks"]:
        if y_col == 1 or y_col == 5 or y_col == 8:
            return 0, 1.4, 'linear'
        elif y_col == 2 or y_col == 6:
            return 0, 1.9, 'linear'
        else:
            return None, None, 'linear'
    elif filename_prefix in ["4Vertices_16Tracks", "4Vertices_20Tracks", "5Vertices_15Tracks"]:
        if y_col == 1 or y_col == 5 or y_col == 8:
            return 1e-4, 10, 'log'
        elif y_col == 2 or y_col == 6:
            return 1e-4, 10, 'log'
        else:
            return None, None, 'linear'
        

def main():
    if len(sys.argv) != 10:
        print("Usage: python script.py <THREADS> <STAGES> <SAMPLES_PER_THREAD> <SWEEPS> <DA_SWEEPS> <input_filename> <output_filenames_prefix> <name_or_extension>")
        sys.exit(1)
    
    THREADS = int(sys.argv[1])
    STAGES = int(sys.argv[2])
    SAMPLES_PER_THREAD = int(sys.argv[3])
    SWEEPS = int(sys.argv[4])
    DA_SWEEPS = int(sys.argv[5])
    if THREADS <= 0 or STAGES <= 0 or SAMPLES_PER_THREAD <= 0 or SWEEPS <= 0 or DA_SWEEPS <= 0:
        print("Error: All numeric arguments must be positive integers.")
        sys.exit(1)
    
    filename_base = sys.argv[6]
    filename_prefix = sys.argv[7]
    extension = sys.argv[8]
    prefix_part = sys.argv[9]
    prefix = prefix_part + "/" + filename_prefix
    
    dir_name_str = os.path.join(prefix, f"{THREADS}threads_{STAGES}stages_{SAMPLES_PER_THREAD}SamplesPerThread_{SWEEPS}sweeps_{DA_SWEEPS}DAsweeps")
    
    
    for filename_png in os.listdir(dir_name_str):
        if filename_png.endswith('.png'):
            file_path_png = os.path.join(dir_name_str, filename_png)
            os.remove(file_path_png)
            print(f"{file_path_png} has been deleted.")
            
    
    os.makedirs(dir_name_str, exist_ok=True)
    # Delete all .png files in the output directory if they exist
    for file in os.listdir(dir_name_str):
        if file.endswith('.png'):
            os.remove(os.path.join(dir_name_str, file))
    
    print("Directories created or already existed: " + dir_name_str)
    
    original1 = "ConvergenceEfficiency_and_TimePerAnneal.txt"
    renamed_original1 = os.path.join(dir_name_str, filename_prefix + "_" + original1)
    
    print(renamed_original1)
    
    if os.path.exists(renamed_original1):
        print("renamed_original1 file exists.")
    else:
        print("Neither original1 nor renamed_original1 file exists.")
        sys.exit(1)
        
        
    original2 = "fit_parameters_and_errors.txt"
    original2_new = "fit_parameters_and_errors_for_slope.txt"
    
    flag_for_slope_eval = False
    
    if os.path.exists(original2):
        print("original2 file exists.")
        
        if flag_for_slope_eval:
            os.rename(original2, original2_new)
        else:
            original2_new = original2
            
        original2 = original2_new
        
        renamed_original2 = os.path.join(prefix_part, original2)
        
        # Make sure target directory exists
        os.makedirs(prefix_part, exist_ok=True)
        
        shutil.move(original2, renamed_original2)
    else:
        print("Original2 file doesn't exist.")
        
        if flag_for_slope_eval:
            original2 = original2_new
        else:
            original2 = original2
        
        renamed_original2 = os.path.join(prefix_part, original2)
        
        # Make sure target directory exists
        os.makedirs(prefix_part, exist_ok=True)
        
    print(renamed_original2)
    
    
    if not flag_for_slope_eval:
        base_dir = filename_base.split('/')[0]
        original3 = f"{base_dir}/QA_valid_correctvalid_efficiency.txt"
        renamed_original3 = original3
        
        print(renamed_original3)
        
        if os.path.exists(renamed_original3):
            print("renamed_original3 file exists.")
        else:
            print("Neither original3 nor renamed_original3 file exists.")
            sys.exit(1)
            
        if filename_prefix in ["2Vertices_10Tracks", "2Vertices_16Tracks", "3Vertices_15Tracks"]:
            avg_time_QA = 45
        else:
            avg_time_QA = 80
            
            
    # Load the data for check
    df1 = pd.read_csv(renamed_original1, sep='\t', header=None)
    df3 = pd.read_csv(renamed_original3, sep='\t', header=None)
    
    # Check for x values match in SA_DA and QA
    dunn1 = np.sort(df1[4].values)
    dunn3 = np.sort(df3[0].values)
    
    diffs = np.abs(dunn1 - dunn3)
    max_diff = np.max(diffs)
    print(f"Maximum absolute difference in sorted Dunn indices: {max_diff}")
    
    filename_initials = filename_prefix+'_SAsweeps_'+str(SWEEPS)+'_DAsweeps_'+str(DA_SWEEPS)

    # Original plots
    bin_edges_5, L_5, std_L_5, k_5, std_k_5, x0_5, std_x0_5, avg_time_5, bin_mids_5, bin_means_5, bin_stderr_5, bin_edges_left_5, bin_edges_right_5, color_5, label_5, fainter_color_5 = fit_and_plot_logistic(
        filename_prefix=filename_prefix,
        data_file=renamed_original1,
        x_col=4,
        y_col=5,
        y_col_list_logistic_fit=[6, 7, 9], 
        y_col_binned_logistic_fit=6, 
        y_col_list_mean=[7, 9], 
        y_col_list_no_logistic_Valid_CorrectValid=[5, 6], 
        output_png=os.path.join(dir_name_str,filename_initials+'_Convergence_efficiency_SA_valid_dunn_index.png'),
        color='red',
        x_label='Dunn Index',
        y_label='SA Valid Solutions Efficiency',
        title='Binned Data and Logistic Fit',
        label='CPU Simulated Annealing'
    )

    bin_edges_6, L_6, std_L_6, k_6, std_k_6, x0_6, std_x0_6, avg_time_6, bin_mids_6, bin_means_6, bin_stderr_6, bin_edges_left_6, bin_edges_right_6, color_6, label_6, fainter_color_6 = fit_and_plot_logistic(
        filename_prefix=filename_prefix,
        data_file=renamed_original1,
        x_col=4,
        y_col=6,
        y_col_list_logistic_fit=[6, 7, 9], 
        y_col_binned_logistic_fit=6, 
        y_col_list_mean=[7, 9], 
        y_col_list_no_logistic_Valid_CorrectValid=[5, 6], 
        output_png=os.path.join(dir_name_str,filename_initials+'_Convergence_efficiency_SA_correctvalid_dunn_index.png'),
        color='red',
        x_label='Dunn Index',
        y_label='SA Efficiency (Correct | Valid)',
        title='Binned Data and Logistic Fit',
        reuse_bin_edges=bin_edges_5,
        label='CPU Simulated Annealing'
    )
    
    # New plot for product of columns 5 and 6 (SA)
    bin_edges_56, L_56, std_L_56, k_56, std_k_56, x0_56, std_x0_56, avg_time_56, bin_mids_56, bin_means_56, bin_stderr_56, bin_edges_left_56, bin_edges_right_56, color_56, label_56, fainter_color_56 = fit_and_plot_logistic(
        filename_prefix=filename_prefix,
        data_file=renamed_original1,
        x_col=4,
        y_col=5,  # y_col is set to 5 but will be overridden by multiply_cols
        y_col_list_logistic_fit=[6, 7, 9], 
        y_col_binned_logistic_fit=6, 
        y_col_list_mean=[7, 9], 
        y_col_list_no_logistic_Valid_CorrectValid=[5, 6], 
        output_png=os.path.join(dir_name_str,filename_initials+'_Convergence_efficiency_SA_correct_dunn_index.png'),
        color='red',
        x_label='Dunn Index',
        y_label='SA Convergence Efficiency',
        title='Binned Data and Logistic Fit',
        reuse_bin_edges=bin_edges_5,
        label='CPU Simulated Annealing',
        multiply_cols=[5, 6]
    )
    
    with open(renamed_original2, 'a') as file:
        file.write(f"{filename_prefix}\tSA\t{THREADS}\t{STAGES}\t{SAMPLES_PER_THREAD}\t{SWEEPS}\t{L_56:.16f}\t{std_L_56:.16f}\t{k_56:.16f}\t{std_k_56:.16f}\t{x0_56:.16f}\t{std_x0_56:.16f}\t")
    
    bin_edges_7, L_7, std_L_7, k_7, std_k_7, x0_7, std_x0_7, avg_time_7, bin_mids_7, bin_means_7, bin_stderr_7, bin_edges_left_7, bin_edges_right_7, color_7, label_7, fainter_color_7 = fit_and_plot_logistic(
        filename_prefix=filename_prefix,
        data_file=renamed_original1,
        x_col=4,
        y_col=7,
        y_col_list_logistic_fit=[6, 7, 9], 
        y_col_binned_logistic_fit=6, 
        y_col_list_mean=[7, 9], 
        y_col_list_no_logistic_Valid_CorrectValid=[5, 6], 
        output_png=os.path.join(dir_name_str,filename_initials+'_Convergence_efficiency_SA_anneal_time.png'),
        color='green',
        x_label='Dunn Index',
        y_label='SA Anneal Time',
        title='Binned Data and Mean Line',
        reuse_bin_edges=bin_edges_5,
        label='Binned means'
    )
    
    with open(renamed_original2, 'a') as file:
        file.write(f"{avg_time_7:.16f}\n")
    
    bin_edges_8, L_8, std_L_8, k_8, std_k_8, x0_8, std_x0_8, avg_time_8, bin_mids_8, bin_means_8, bin_stderr_8, bin_edges_left_8, bin_edges_right_8, color_8, label_8, fainter_color_8 = fit_and_plot_logistic(
        filename_prefix=filename_prefix,
        data_file=renamed_original1,
        x_col=4,
        y_col=8,
        y_col_list_logistic_fit=[6, 7, 9], 
        y_col_binned_logistic_fit=6, 
        y_col_list_mean=[7, 9], 
        y_col_list_no_logistic_Valid_CorrectValid=[5, 6], 
        output_png=os.path.join(dir_name_str,filename_initials+'_Convergence_efficiency_DA_correct_dunn_index.png'),
        color='#b48a00',
        x_label='Dunn Index',
        y_label='DA Convergence Efficiency',
        title='Binned Data and Logistic Fit',
        label='CPU Deterministic Annealing'
    )
    
    with open(renamed_original2, 'a') as file:
        file.write(f"{filename_prefix}\tDA\t{THREADS}\t{STAGES}\t{SAMPLES_PER_THREAD}\t{DA_SWEEPS}\t{L_8:.16f}\t{std_L_8:.16f}\t{k_8:.16f}\t{std_k_8:.16f}\t{x0_8:.16f}\t{std_x0_8:.16f}\t")
    
    bin_edges_9, L_9, std_L_9, k_9, std_k_9, x0_9, std_x0_9, avg_time_9, bin_mids_9, bin_means_9, bin_stderr_9, bin_edges_left_9, bin_edges_right_9, color_9, label_9, fainter_color_9 = fit_and_plot_logistic(
        filename_prefix=filename_prefix,
        data_file=renamed_original1,
        x_col=4,
        y_col=9,
        y_col_list_logistic_fit=[6, 7, 9], 
        y_col_binned_logistic_fit=6, 
        y_col_list_mean=[7, 9], 
        y_col_list_no_logistic_Valid_CorrectValid=[5, 6], 
        output_png=os.path.join(dir_name_str,filename_initials+'_Convergence_efficiency_DA_anneal_time.png'),
        color='purple',
        x_label='Dunn Index',
        y_label='DA Anneal Time',
        title='Binned Data and Mean Line',
        reuse_bin_edges=bin_edges_8,
        label='Binned means'
    )
    
    with open(renamed_original2, 'a') as file:
        file.write(f"{avg_time_9:.16f}\n")
        
    
    if not flag_for_slope_eval:
        bin_edges_1, L_1, std_L_1, k_1, std_k_1, x0_1, std_x0_1, avg_time_1, bin_mids_1, bin_means_1, bin_stderr_1, bin_edges_left_1, bin_edges_right_1, color_1, label_1, fainter_color_1 = fit_and_plot_logistic(
            filename_prefix=filename_prefix,
            data_file=renamed_original3,
            x_col=0,
            y_col=1,
            y_col_list_logistic_fit=[2], 
            y_col_binned_logistic_fit=2, 
            y_col_list_mean=[], 
            y_col_list_no_logistic_Valid_CorrectValid=[1, 2], 
            output_png=os.path.join(dir_name_str,filename_initials+'_Convergence_efficiency_QA_valid_dunn_index.png'),
            color='black',
            x_label='Dunn Index',
            y_label='QA Valid Solutions Efficiency',
            title='Binned Data and Logistic Fit',
            label='QPU Quantum Annealing'
        )
        
        bin_edges_2, L_2, std_L_2, k_2, std_k_2, x0_2, std_x0_2, avg_time_2, bin_mids_2, bin_means_2, bin_stderr_2, bin_edges_left_2, bin_edges_right_2, color_2, label_2, fainter_color_2 = fit_and_plot_logistic(
            filename_prefix=filename_prefix,
            data_file=renamed_original3,
            x_col=0,
            y_col=2,
            y_col_list_logistic_fit=[2], 
            y_col_binned_logistic_fit=2, 
            y_col_list_mean=[], 
            y_col_list_no_logistic_Valid_CorrectValid=[1, 2], 
            output_png=os.path.join(dir_name_str,filename_initials+'_Convergence_efficiency_QA_correctvalid_dunn_index.png'),
            color='black',
            x_label='Dunn Index',
            y_label='QA Efficiency (Correct | Valid)',
            title='Binned Data and Logistic Fit',
            reuse_bin_edges=bin_edges_1,
            label='QPU Quantum Annealing'
        )
            
        # New plot for product of columns 5 and 6 (QA)
        bin_edges_12, L_12, std_L_12, k_12, std_k_12, x0_12, std_x0_12, avg_time_12, bin_mids_12, bin_means_12, bin_stderr_12, bin_edges_left_12, bin_edges_right_12, color_12, label_12, fainter_color_12 = fit_and_plot_logistic(
            filename_prefix=filename_prefix,
            data_file=renamed_original3,
            x_col=0,
            y_col=1,  
            y_col_list_logistic_fit=[2], 
            y_col_binned_logistic_fit=2, 
            y_col_list_mean=[], 
            y_col_list_no_logistic_Valid_CorrectValid=[1, 2], 
            output_png=os.path.join(dir_name_str,filename_initials+'_Convergence_efficiency_QA_correct_dunn_index.png'),
            color='black',
            x_label='Dunn Index',
            y_label='QA Convergence Efficiency',
            title='Binned Data and Logistic Fit',
            reuse_bin_edges=bin_edges_1,
            label='QPU Quantum Annealing',
            multiply_cols=[1, 2]
        )
            
        with open(renamed_original2, 'a') as file:
            file.write(f"{filename_prefix}\tQA\t{THREADS}\t{STAGES}\t{SAMPLES_PER_THREAD}\t{SWEEPS}\t{L_12:.16f}\t{std_L_12:.16f}\t{k_12:.16f}\t{std_k_12:.16f}\t{x0_12:.16f}\t{std_x0_12:.16f}\t{avg_time_QA:.16f}\n")
        
        # Call combined plot function
        # First combined plot: SA valid (col 5) + QA valid (col 1)
        plot_combined_correct_dunn_index(
            filename_prefix=filename_prefix,
            data_file=renamed_original1,
            x_col=4,
            y_col=5,
            bin_mids_SA=bin_mids_5, bin_means_SA=bin_means_5, bin_stderr_SA=bin_stderr_5,
            bin_edges_left_SA=bin_edges_left_5, bin_edges_right_SA=bin_edges_right_5,
            color_SA=color_5, label_SA=label_5, fainter_color_SA=fainter_color_5,
            x_col_SA=4,
            popt_SA=None,  # Override to skip curve
            # DA skipped (defaults to None)
            bin_mids_QA=bin_mids_1, bin_means_QA=bin_means_1, bin_stderr_QA=bin_stderr_1,
            bin_edges_left_QA=bin_edges_left_1, bin_edges_right_QA=bin_edges_right_1,
            color_QA=color_1, label_QA=label_1, fainter_color_QA=fainter_color_1,
            x_col_QA=0,
            popt_QA=None,  # Override to skip curve
            output_png=os.path.join(dir_name_str, filename_initials+'_Convergence_efficiency_valid_dunn_index.png'),
            x_label='Dunn Index', y_label='Valid solutions Efficiency',
            title='Binned Data and Logistic Fit'
        )
        
        # Second combined plot: SA correct|valid (col 6) + QA correct|valid (col 2)
        plot_combined_correct_dunn_index(
            filename_prefix=filename_prefix,
            data_file=renamed_original1,
            x_col=4,
            y_col=6,
            bin_mids_SA=bin_mids_6, bin_means_SA=bin_means_6, bin_stderr_SA=bin_stderr_6,
            bin_edges_left_SA=bin_edges_left_6, bin_edges_right_SA=bin_edges_right_6,
            color_SA=color_6, label_SA=label_6, fainter_color_SA=fainter_color_6,
            x_col_SA=4,
            popt_SA=None,  # Already None, but explicit for clarity
            # DA skipped (defaults to None)
            bin_mids_QA=bin_mids_2, bin_means_QA=bin_means_2, bin_stderr_QA=bin_stderr_2,
            bin_edges_left_QA=bin_edges_left_2, bin_edges_right_QA=bin_edges_right_2,
            color_QA=color_2, label_QA=label_2, fainter_color_QA=fainter_color_2,
            x_col_QA=0,
            popt_QA=None,  # Already None, but explicit for clarity
            output_png=os.path.join(dir_name_str, filename_initials+'_Convergence_efficiency_correctvalid_dunn_index.png'),
            x_label='Dunn Index', y_label='Efficiency (Correct|Valid)',
            title='Binned Data and Logistic Fit'
        )
        
        plot_combined_correct_dunn_index(
            filename_prefix=filename_prefix,
            data_file=renamed_original1,
            x_col=4,
            y_col=5,
            bin_mids_SA=bin_mids_56, bin_means_SA=bin_means_56, bin_stderr_SA=bin_stderr_56,
            bin_edges_left_SA=bin_edges_left_56, bin_edges_right_SA=bin_edges_right_56,
            color_SA=color_56, label_SA=label_56, fainter_color_SA=fainter_color_56,
            x_col_SA=4,
            popt_SA=[L_56, k_56, x0_56],  # Already None, but explicit for clarity
            bin_mids_DA=bin_mids_8, bin_means_DA=bin_means_8, bin_stderr_DA=bin_stderr_8,
            bin_edges_left_DA=bin_edges_left_8, bin_edges_right_DA=bin_edges_right_8,
            color_DA=color_8, label_DA=label_8, fainter_color_DA=fainter_color_8,
            x_col_DA=4,
            popt_DA=[L_8, k_8, x0_8],  # Already None, but explicit for clarity
            bin_mids_QA=bin_mids_12, bin_means_QA=bin_means_12, bin_stderr_QA=bin_stderr_12,
            bin_edges_left_QA=bin_edges_left_12, bin_edges_right_QA=bin_edges_right_12,
            color_QA=color_12, label_QA=label_12, fainter_color_QA=fainter_color_12,
            x_col_QA=0,
            popt_QA=[L_12, k_12, x0_12],  # Already None, but explicit for clarity
            output_png=os.path.join(dir_name_str, filename_initials+'_Convergence_efficiency_correct_dunn_index.png'),
            x_label='Dunn Index',
            y_label='Convergence Efficiency',
            title='Binned Data and Logistic Fit'
        )

    else:
        # Call combined plot function
        plot_combined_correct_dunn_index(
            filename_prefix=filename_prefix,
            data_file=renamed_original1,
            x_col=4,
            y_col=5,
            bin_mids_SA=bin_mids_56, bin_means_SA=bin_means_56, bin_stderr_SA=bin_stderr_56,
            bin_edges_left_SA=bin_edges_left_56, bin_edges_right_SA=bin_edges_right_56,
            color_SA=color_56, label_SA=label_56, fainter_color_SA=fainter_color_56,
            x_col_SA=4,
            popt_SA=[L_56, k_56, x0_56],  # Already None, but explicit for clarity
            bin_mids_DA=bin_mids_8, bin_means_DA=bin_means_8, bin_stderr_DA=bin_stderr_8,
            bin_edges_left_DA=bin_edges_left_8, bin_edges_right_DA=bin_edges_right_8,
            color_DA=color_8, label_DA=label_8, fainter_color_DA=fainter_color_8,
            x_col_DA=4,
            popt_DA=[L_8, k_8, x0_8],  # Already None, but explicit for clarity
            # QA dropped
            output_png=os.path.join(dir_name_str, filename_initials+'_Convergence_efficiency_correct_dunn_index.png'),
            x_label='Dunn Index',
            y_label='Convergence Efficiency',
            title='Binned Data and Logistic Fit'
        )
    
        
    for filename_png in os.listdir(dir_name_str):
        if filename_png.endswith('.png') and filename_png not in [filename_initials + "_Convergence_efficiency_valid_dunn_index.png",
        filename_initials + "_Convergence_efficiency_correct_dunn_index.png", filename_initials + "_Convergence_efficiency_correctvalid_dunn_index.png",
        filename_initials + "_Convergence_efficiency_SA_anneal_time.png", filename_initials + "_Convergence_efficiency_DA_anneal_time.png"]:
            file_path_png = os.path.join(dir_name_str, filename_png)
            os.remove(file_path_png)
            print(f"{file_path_png} has been deleted.")

if __name__ == "__main__":
    main()