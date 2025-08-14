import sys
import os
import shutil
import math
import numpy as np
from scipy.optimize import curve_fit
import subprocess

def logistic(x, L, k, x0):
    return L / (1 + np.exp(-k * (x - x0)))

def sorted_x_y(x, y):
    sorted_indices = sorted(range(len(x)), key=lambda ii: x[ii])
    x_sorted = [x[ii] for ii in sorted_indices]
    y_sorted = [y[ii] for ii in sorted_indices]
    return x_sorted, y_sorted

def compute_mse(x, y, L, k, x0):
    mse = 0.0
    for i in range(len(x)):
        expected = logistic(x[i], L, k, x0)
        mse += (y[i] - expected) ** 2
    return mse / len(x) if len(x) > 0 else 0.0

def fit_logistic(x_sorted, y_sorted, p0=[1.0, 0.1, 30.0]):
    try:
        params, _ = curve_fit(logistic, x_sorted, y_sorted, p0=p0, bounds=([0.5, 0.0, -np.inf], [1.5, 1.0, np.inf]), maxfev=1000)
    except RuntimeError:
        print("Curve fit failed, using default parameters")
        params = p0
    L, k, x0 = params
    if L < 0.5 or L > 1.5 or k <= 0:
        print(f"Warning: Invalid parameters detected (L={L}, k={k}, x0={x0}), using default [1, 0.1, 30]")
        L, k, x0 = 1.0, 0.1, 30.0
    return L, k, x0

def main():
    if len(sys.argv) != 9:
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
    prefix = sys.argv[7]
    extension = sys.argv[8]

    dir_name_str = prefix + "/" + str(THREADS) + "threads_" + \
                   str(STAGES) + "stages_" + \
                   str(SAMPLES_PER_THREAD) + "SamplesPerThread_" + \
                   str(SWEEPS) + "sweeps_" + \
                   str(DA_SWEEPS) + "DAsweeps"
    os.makedirs(dir_name_str, exist_ok=True)
    print("Directories created or already existed: " + dir_name_str)

    original1 = "ConvergenceEfficiency_and_TimePerAnneal.txt"
    original2 = "ConvergenceEfficiency_and_DunnIndex_Binning.txt"
    renamed_original1 = dir_name_str + "/" + prefix + "_" + original1
    renamed_original2 = dir_name_str + "/" + prefix + "_SA_" + original2
    renamed_original3 = dir_name_str + "/" + prefix + "_DA_" + original2

    if os.path.exists(renamed_original1):
        print("renamed_original1 file exists.")
    else:
        print("Neither original1 nor renamed_original1 file exists.")
        sys.exit(1)

    if os.path.exists(original2):
        if '.' in original2:
            temp_file = original2.rsplit('.', 1)[0] + "_tmp.txt"
        else:
            temp_file = original2 + "_tmp.txt"
        print(f"original2 file exists. Copying to {temp_file} and renaming to {renamed_original2}")
        shutil.copyfile(original2, temp_file)
        os.rename(temp_file, renamed_original2)

    if os.path.exists(original2):
        print(f"original2 file exists. Renaming to {renamed_original3}")
        os.rename(original2, renamed_original3)

    # Read data
    data = []
    with open(renamed_original1, 'r') as file:
        lines = file.readlines()
    for line in lines:
        try:
            row = [float(v) for v in line.split()]
            if len(row) != 10 or any(math.isnan(v) or math.isinf(v) for v in row):
                print("Error: Invalid data format in line: " + line.strip())
                continue
            data.append(row)
        except ValueError:
            print("Error: Invalid data format in line: " + line.strip())
            continue
    if not data:
        print("Error: No valid data read from " + renamed_original1)
        sys.exit(1)

    NUM_COLUMNS = 5
    bin_means = [[] for _ in range(NUM_COLUMNS)]
    bin_stderr = [[] for _ in range(NUM_COLUMNS)]
    saved_bins_x = [[] for _ in range(NUM_COLUMNS)]
    saved_bins_y = [[] for _ in range(NUM_COLUMNS)]
    saved_bin_edges = [[] for _ in range(NUM_COLUMNS)]

    min_bins = 5
    
    L_param = [0,0]
    k_param = [0,0]
    x0_param = [0,0]
    
    for ind in range(NUM_COLUMNS):
        x_data = [row[4] for row in data]
        y_data = [row[5 + ind] for row in data]

        bins_x = []
        bins_y = []
        bin_edges = []

        if ind == 0 or ind == 3:
            ind_new = ind // 3
            
            x_sorted, y_sorted = sorted_x_y(x_data, y_data)
            
            L, k, x0 = fit_logistic(x_sorted, y_sorted)
            L_param[ind_new], k_param[ind_new], x0_param[ind_new] = L, k, x0
            
            # Compute final MSE
            mse = compute_mse(x_sorted, y_sorted, L, k, x0)
            print(f"Final MSE: {mse:.6f}")
            
            # Adaptive binning
            num_bins_target = 8
            min_points_per_bin = 5
            bin_edges = [x_sorted[0] - 0.01]
            current_bin_x = []
            current_bin_y = []
            for i in range(len(x_sorted)):
                current_bin_x.append(x_sorted[i])
                current_bin_y.append(y_sorted[i])
                if len(current_bin_x) >= min_points_per_bin:
                    mean_y = np.mean(current_bin_y)
                    mid_x = (current_bin_x[0] + current_bin_x[-1]) / 2.0
                    expected_y = logistic(mid_x, L, k, x0)
                    mse_current = (mean_y - expected_y) ** 2
                    print(f"Bin attempt at index {i}: MSE = {mse_current:.6f}")
                    if mse_current < 0.0005 or i == len(x_sorted) - 1 or len(bins_x) == num_bins_target - 1:
                        bins_x.append(list(current_bin_x))
                        bins_y.append(list(current_bin_y))
                        bin_edges.append(current_bin_x[-1] + 0.01)
                        current_bin_x = []
                        current_bin_y = []
            if current_bin_x:
                bins_x.append(list(current_bin_x))
                bins_y.append(list(current_bin_y))
                bin_edges.append(current_bin_x[-1] + 0.01)

            # Ensure at least min_bins by splitting largest
            while len(bins_x) < min_bins:
                sizes = [len(b) for b in bins_x]
                if not sizes:
                    break
                largest_idx = np.argmax(sizes)
                max_size = sizes[largest_idx]
                if max_size < 2 * min_points_per_bin:
                    break
                split_index = max_size // 2
                if split_index < min_points_per_bin or (max_size - split_index) < min_points_per_bin:
                    break
                left_x = bins_x[largest_idx][:split_index]
                left_y = bins_y[largest_idx][:split_index]
                right_x = bins_x[largest_idx][split_index:]
                right_y = bins_y[largest_idx][split_index:]
                new_edge = (left_x[-1] + right_x[0]) / 2.0
                bin_edges.insert(largest_idx + 1, new_edge)
                bins_x[largest_idx] = left_x
                bins_y[largest_idx] = left_y
                bins_x.insert(largest_idx + 1, right_x)
                bins_y.insert(largest_idx + 1, right_y)

            if len(bins_x) < min_bins:
                print(f"Warning: Only {len(bins_x)} bins created for ind={ind}.")

            # Debug bin statistics
            print(f"Bin Edges for ind={ind}: " + " ".join([f"{e:.6f}" for e in bin_edges]))
            for i in range(len(bins_y)):
                mean = np.mean(bins_y[i])
                variance = np.var(bins_y[i], ddof=1) if len(bins_y[i]) > 1 else 0.0
                stderr = np.sqrt(variance) / np.sqrt(len(bins_y[i])) if len(bins_y[i]) > 1 else 0.0
                print(f"Bin {bin_edges[i]:.6f}-{bin_edges[i+1]:.6f}: Mean = {mean:.6f}, Std Error = {stderr:.6f}, Points = {len(bins_y[i])}")

            saved_bins_x[ind] = bins_x
            saved_bins_y[ind] = bins_y
            saved_bin_edges[ind] = bin_edges
        else:
            base_ind = 0 if ind < 3 else 3
            bin_edges = saved_bin_edges[base_ind]
            bins_y = []
            for bb in range(len(bin_edges) - 1):
                lower = bin_edges[bb]
                upper = bin_edges[bb + 1]
                bin_y_temp = [y_data[ii] for ii in range(len(x_data)) if lower <= x_data[ii] < upper]
                bins_y.append(bin_y_temp)
            saved_bins_y[ind] = bins_y
            saved_bin_edges[ind] = bin_edges

        # Compute bin means and stderr
        for i in range(len(saved_bins_y[ind])):
            if not saved_bins_y[ind][i]:
                continue
            mean = np.mean(saved_bins_y[ind][i])
            bin_means[ind].append(mean)
            variance = np.var(saved_bins_y[ind][i], ddof=1) if len(saved_bins_y[ind][i]) > 1 else 0.0
            stderr = np.sqrt(variance) / np.sqrt(len(saved_bins_y[ind][i])) if len(saved_bins_y[ind][i]) > 1 else 0.0
            bin_stderr[ind].append(stderr)

    # Write to outFile2
    with open(renamed_original2, 'a') as outFile2:
        for j in range(len(saved_bin_edges[0]) - 1):
            center = (saved_bin_edges[0][j] + saved_bin_edges[0][j + 1]) / 2.0
            lower = saved_bin_edges[0][j]
            upper = saved_bin_edges[0][j + 1]
            outFile2.write(f"{THREADS}\t{STAGES}\t{SAMPLES_PER_THREAD}\t{SWEEPS}\t{DA_SWEEPS}\t{center:.6f}\t{lower:.6f}\t{upper:.6f}")
            for k in range(3):
                mean = bin_means[k][j] if j < len(bin_means[k]) else 0.0
                err = bin_stderr[k][j] if j < len(bin_stderr[k]) else 0.0
                outFile2.write(f"\t{mean:.6f}\t{err:.6f}")
            outFile2.write("\n")

    # Write to outFile3
    with open(renamed_original3, 'a') as outFile3:
        for j in range(len(saved_bin_edges[3]) - 1):
            center = (saved_bin_edges[3][j] + saved_bin_edges[3][j + 1]) / 2.0
            lower = saved_bin_edges[3][j]
            upper = saved_bin_edges[3][j + 1]
            outFile3.write(f"{THREADS}\t{STAGES}\t{SAMPLES_PER_THREAD}\t{SWEEPS}\t{DA_SWEEPS}\t{center:.6f}\t{lower:.6f}\t{upper:.6f}")
            for k in range(3, 5):
                mean = bin_means[k][j] if j < len(bin_means[k]) else 0.0
                err = bin_stderr[k][j] if j < len(bin_stderr[k]) else 0.0
                outFile3.write(f"\t{mean:.6f}\t{err:.6f}")
            outFile3.write("\n")
            
    print(f"Start fit the valid*correct|valid")
    L, k, x0 = fit_logistic([row[4] for row in data], [row[6] for row in data])
    L2, k2, x02 = fit_logistic([row[4] for row in data], [row[5]*row[6] for row in data])
    print(f"{L}, {k}, {x0}")
    print(f"{L2}, {k2}, {x02}")
    print(f"End fit the valid*correct|valid")

    # Generate gnuplot scripts and run them
    gnuplot_script_da1 = os.path.join(dir_name_str, "plot_da_k3.gp")
    gnuplot_script_da2 = os.path.join(dir_name_str, "plot_da_k4.gp")
    gnuplot_script_sa1 = os.path.join(dir_name_str, "plot_sa_k0_k1.gp")
    gnuplot_script_sa2 = os.path.join(dir_name_str, "plot_sa_k2_k0.gp")
    gnuplot_script_sa_da_correct = os.path.join(dir_name_str, "plot_sa_da_correct.gp")
    
    x_data = [row[4] for row in data]
    max_x_data = round(max(x_data))
    
    plot_range = f"set xrange [0.5:{max_x_data}]\nset yrange [-0.01:1.4]"
    plot_range2 = f"set xrange [0.5:{max_x_data}]"
    
    logscale_x = True
    logscale_y = False
    
    logscale_cmds = ""
    if logscale_x:
        logscale_cmds += "set logscale x\n"
    if logscale_y:
        logscale_cmds += "set logscale y\n"

    # DA Plot 1: Column 6 (x), Column 9 (y, k=3) with errorbars
    with open(gnuplot_script_da1, 'w') as gp:
        gp.write(f"""
set terminal pngcairo size 800,800 enhanced font 'Arial,12'
set output '{dir_name_str}/{prefix}_DA_efficiency_vs_dunn_index.png'
set title 'DA efficiency vs dunn index'
set xlabel 'Dunn index'
set ylabel 'DA correct efficiency'
{logscale_cmds}{plot_range}
set style data errorbars
set bars small

# Define the logistic functions
f_da(x) = {L_param[1]}/(1 + exp(-{k_param[1]} * (x - {x0_param[1]})))

plot '{renamed_original3}' using 6:9:7:8:($9-$10):($9+$10) with xyerrorbars title 'DA correct efficiency' pt 7 ps 1 lc rgb 'blue', \\
     f_da(x) title sprintf('DA correct fit: %.3f / (1 + exp(-%.3f * (x - %.3f)))', {L_param[1]}, {k_param[1]}, {x0_param[1]}) lc rgb 'web-blue'
""")

    # DA Plot 2: Column 6 (x), Column 11 (y, k=4) with errorbars
    with open(gnuplot_script_da2, 'w') as gp:
        gp.write(f"""
set terminal pngcairo size 800,800 enhanced font 'Arial,12'
set output '{dir_name_str}/{prefix}_DA_avg_anneal_time_vs_dunn_index.png'
set title 'DA anneal time vs dunn index'
set xlabel 'Dunn index'
set ylabel 'DA average anneal time'
{logscale_cmds}{plot_range2}
set style data errorbars
set bars small
plot '{renamed_original3}' using 6:11:7:8:($11-$12):($11+$12) with xyerrorbars title 'DA average anneal time' pt 7 ps 1
""")

    # SA Plot 1: Column 6 (x), Columns 8 (k=0) and 10 (k=1) on y, same plot
    with open(gnuplot_script_sa1, 'w') as gp:
        gp.write(f"""
set terminal pngcairo size 800,800 enhanced font 'Arial,12'
set output '{dir_name_str}/{prefix}_SA_efficiencies_vs_dunn_index.png'
set title 'SA efficiencies vs dunn index'
set xlabel 'Dunn index'
set ylabel 'SA valid & correct|valid efficiencies'
{logscale_cmds}{plot_range}
set style data errorbars
set bars small

# Define the logistic functions
f_sa_valid(x) = {L_param[0]}/(1+exp(-{k_param[0]}*(x-{x0_param[0]})))
f_sa_correctvalid(x) = {L}/(1 + exp(-{k} * (x - {x0})))

plot '{renamed_original2}' using 6:9:7:8:($9-$10):($9+$10) with xyerrorbars title 'SA valid efficiency' pt 7 ps 1 lc rgb 'blue', \\
     f_sa_valid(x) w l title sprintf('SA valid fit: %.3f / (1 + exp(-%.3f * (x - %.3f)))', {L_param[0]}, {k_param[0]}, {x0_param[0]}) lc rgb 'web-blue', \\
     '{renamed_original2}' using 6:11:7:8:($11-$12):($11+$12) with xyerrorbars title 'SA correct|valid efficiency' pt 9 ps 1 lc rgb 'red', \\
     f_sa_correctvalid(x) w l title sprintf('SA correct|valid fit: %.3f / (1 + exp(-%.3f * (x - %.3f)))', {L}, {k}, {x0}) lc rgb 'dark-salmon'
""")

    # SA Plot 2: Column 6 (x), Columns 12 (k=2) and 8 (k=0) on y, same plot
    with open(gnuplot_script_sa2, 'w') as gp:
        gp.write(f"""
set terminal pngcairo size 800,800 enhanced font 'Arial,12'
set output '{dir_name_str}/{prefix}_SA_avg_anneal_time_vs_dunn_index.png'
set title 'SA anneal time vs dunn index'
set xlabel 'Dunn index'
set ylabel 'SA average anneal time'
{logscale_cmds}{plot_range2}
set style data errorbars
set bars small
plot '{renamed_original2}' using 6:13:7:8:($13-$14):($13+$14) with xyerrorbars title 'SA average anneal time' pt 7 ps 1
""")

    # SA and DA correct efficiency
    with open(gnuplot_script_sa_da_correct, 'w') as gp:
        gp.write(f"""
set terminal pngcairo size 800,800 enhanced font 'Arial,12'
set output '{dir_name_str}/{prefix}_convergence_efficiencies_vs_dunn_index.png'
set title 'Convergence efficiencies vs dunn index'
set xlabel 'Dunn index'
set ylabel 'Convergence efficiency'
{logscale_cmds}{plot_range}
set style data errorbars
set bars small

# Define the logistic functions
f_sa(x) = {L2}/(1 + exp(-{k2} * (x - {x02})))
f_da(x) = {L_param[1]}/(1 + exp(-{k_param[1]} * (x - {x0_param[1]})))

plot '{renamed_original2}' using 6:($9*$11):7:8:(($9*$11)-((($10*$11)**2)+(($9*$12)**2))**0.5):(($9*$11)+((($10*$11)**2)+(($9*$12)**2))**0.5) with xyerrorbars title 'SA correct efficiency' pt 7 ps 1 lc rgb 'blue', \\
     f_sa(x) title sprintf('SA correct fit: %.3f / (1 + exp(-%.3f * (x - %.3f)))', {L2}, {k2}, {x02}) lc rgb 'web-blue', \\
     '{renamed_original3}' using 6:9:7:8:($9-$10):($9+$10) with xyerrorbars title 'DA correct efficiency' pt 7 ps 1 lc rgb 'red', \\
     f_da(x) title sprintf('DA correct fit: %.3f / (1 + exp(-%.3f * (x - %.3f)))', {L_param[1]}, {k_param[1]}, {x0_param[1]}) lc rgb 'dark-salmon'
""")

    # Run gnuplot commands
    try:
        for script in [gnuplot_script_da1, gnuplot_script_da2, gnuplot_script_sa1, gnuplot_script_sa2, gnuplot_script_sa_da_correct]:
            subprocess.run(['gnuplot', script], check=True)
            print(f"Generated plot: {script}")
    except subprocess.CalledProcessError as e:
        print(f"Error running gnuplot for {script}: {e}")
    except FileNotFoundError:
        print("Error: gnuplot not found. Please install it in your conda environment with 'conda install gnuplot'.")

    os.remove(gnuplot_script_da1)
    os.remove(gnuplot_script_da2)
    os.remove(gnuplot_script_sa1)
    os.remove(gnuplot_script_sa2)
    os.remove(gnuplot_script_sa_da_correct)

if __name__ == "__main__":
    main()