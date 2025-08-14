#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <filesystem>
#include <ceres/ceres.h>
//#include <glog/logging.h>
//#include <Eigen/Core>

// Logistic function
double logistic(double x, double L, double k, double x0) {
    return L / (1.0 + std::exp(-k * (x - x0)));
}

// Compute MSE for logistic fit
double compute_mse(const std::vector<double>& x, const std::vector<double>& y, double L, double k, double x0) {
    double mse = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double expected = logistic(x[i], L, k, x0);
        mse += (y[i] - expected) * (y[i] - expected);
    }
    return mse / x.size();
}

// Cost functor for Ceres
struct LogisticCost {
    LogisticCost(double x, double y) : x_(x), y_(y) {}
    template <typename T>
    bool operator()(const T* const L, const T* const k, const T* const x0, T* residual) const {
        T predicted = *L / (T(1.0) + exp(-(*k) * (T(x_) - *x0)));
        residual[0] = T(y_) - predicted;
        return true;
    }
private:
    const double x_, y_;
};

// Ceres-based optimization for logistic fit
void ceres_curve_fit(std::vector<double>& params, const std::vector<double>& x, const std::vector<double>& y, double L_min, double L_max) {
    if (x.empty() || y.empty()) {
        std::cerr << "Error: Empty x or y data passed to ceres_curve_fit." << std::endl;
        return;
    }

    // Initialize parameters
    params = {1.0, 0.1, 30.0}; // Match Python's p0
    double L = params[0], k = params[1], x0 = params[2];

    // Setup Ceres problem
    ceres::Problem problem;
    for (size_t i = 0; i < x.size(); ++i) {
        problem.AddResidualBlock(
            new ceres::AutoDiffCostFunction<LogisticCost, 1, 1, 1, 1>(new LogisticCost(x[i], y[i])),
            nullptr, &L, &k, &x0);
    }

    // Set bounds to match Python
    problem.SetParameterLowerBound(&L, 0, L_min); // L >= 0.5
    problem.SetParameterUpperBound(&L, 0, L_max); // L <= 1.5
    problem.SetParameterLowerBound(&k, 0, 0.0);   // k >= 0
    problem.SetParameterUpperBound(&k, 0, 1.0);   // k <= 1
    // x0 is unbounded, no bounds set

    // Configure solver
    ceres::Solver::Options options;
    options.minimizer_type = ceres::TRUST_REGION;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 1000;
    options.function_tolerance = 1e-10;
    options.parameter_tolerance = 1e-8;
    ceres::Solver::Summary summary;

    // Run solver
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << std::endl;

    // Validate parameters
    if (L < 0.5 || L > 1.5 || k <= 0) {
        std::cerr << "Warning: Invalid parameters detected (L=" << L << ", k=" << k << ", x0=" << x0 << "), using default [1, 0.1, 30]" << std::endl;
        L = 1.0;
        k = 0.1;
        x0 = 30.0;
    }

    // Update output parameters
    params = {L, k, x0};

    // Compute final MSE for debugging
    double mse = compute_mse(x, y, L, k, x0);
    std::cout << "Final MSE: " << std::setprecision(6) << mse << std::endl;
}

// Compute MSE for bin edges
double compute_bin_mse(const std::vector<double>& bin_edges, const std::vector<double>& x, const std::vector<double>& y, double L, double k, double x0) {
    double mse = 0.0;
    size_t start = 0;
    size_t bin_count = 0;
    for (size_t i = 1; i < bin_edges.size(); ++i) {
        double edge = bin_edges[i];
        size_t idx = std::lower_bound(x.begin(), x.end(), edge) - x.begin();
        if (idx > start && idx <= x.size()) {
            std::vector<double> bin_y(y.begin() + start, y.begin() + idx);
            if (!bin_y.empty()) {
                double mean_y = 0.0;
                for (double val : bin_y) mean_y += val;
                mean_y /= bin_y.size();
                double mid_x = (x[start] + x[idx - 1]) / 2.0;
                double expected_y = logistic(mid_x, L, k, x0);
                mse += (mean_y - expected_y) * (mean_y - expected_y);
                ++bin_count;
            }
            start = idx;
        }
    }
    return bin_count > 0 ? mse / bin_count : mse;
}

// Function to check if file exists
bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

// Function to create directories recursively
bool createDirectoriesRecursively(const std::string& path) {
    try {
        std::filesystem::create_directories(path);
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error creating directories: " << e.what() << std::endl;
        return false;
    }
}

int main(int argc, char* argv[]) {
    // Initialize Google logging for Ceres
    //google::InitGoogleLogging(argv[0]);

    if (argc != 9) {
        std::cerr << "Usage: " << argv[0] << " <THREADS> <STAGES> <SAMPLES_PER_THREAD> <SWEEPS> <DA_SWEEPS> <input_filename> <output_filenames_prefix> <name_or_extension>" << std::endl;
        return 1;
    }

    int THREADS = std::atoi(argv[1]);
    int STAGES = std::atoi(argv[2]);
    int SAMPLES_PER_THREAD = std::atoi(argv[3]);
    int SWEEPS = std::atoi(argv[4]);
    int DA_SWEEPS = std::atoi(argv[5]);
    if (THREADS <= 0 || STAGES <= 0 || SAMPLES_PER_THREAD <= 0 || SWEEPS <= 0 || DA_SWEEPS <= 0) {
        std::cerr << "Error: All numeric arguments must be positive integers." << std::endl;
        return 1;
    }

    std::string filename_base = argv[6];
    std::string prefix = argv[7];
    std::string extension = argv[8];

    std::string dir_name_str = prefix + "/" + std::to_string(THREADS) + "threads_" + 
                               std::to_string(STAGES) + "stages_" + 
                               std::to_string(SAMPLES_PER_THREAD) + "SamplesPerThread_" + 
                               std::to_string(SWEEPS) + "sweeps_" + 
                               std::to_string(DA_SWEEPS) + "DAsweeps";
    if (!createDirectoriesRecursively(dir_name_str)) {
        std::cerr << "Failed to create directories: " << dir_name_str << std::endl;
        return 1;
    }
    std::cout << "Directories created or already existed: " << dir_name_str << std::endl;

    std::string original1 = "ConvergenceEfficiency_and_TimePerAnneal.txt";
    std::string original2 = "ConvergenceEfficiency_and_DunnIndex_Binning.txt";
    std::string renamed_original1 = dir_name_str + "/" + prefix + "_" + original1;
    std::string renamed_original2 = dir_name_str + "/" + prefix + "_SA_" + original2;
    std::string renamed_original3 = dir_name_str + "/" + prefix + "_DA_" + original2;

    if (file_exists(renamed_original1)) {
        std::cout << "renamed_original1 file exists." << std::endl;
    } else {
        std::cerr << "Neither original1 nor renamed_original1 file exists." << std::endl;
        return 1;
    }

    if (file_exists(original2)) {
        std::string temp_file;
        if (original2.size() >= 4 && original2.substr(original2.size() - 4) == ".txt") {
            temp_file = original2.substr(0, original2.size() - 4) + "_tmp.txt";
        } else {
            temp_file = original2 + "_tmp.txt";
        }
        std::cout << "original2 file exists. Copying to " << temp_file << " and renaming to " << renamed_original2 << std::endl;
        try {
            if (!std::filesystem::copy_file(original2, temp_file, std::filesystem::copy_options::overwrite_existing)) {
                std::cerr << "Error copying original2 to " << temp_file << std::endl;
                return 1;
            }
            if (std::rename(temp_file.c_str(), renamed_original2.c_str()) != 0) {
                std::perror("Error renaming temporary file to renamed_original2");
                std::filesystem::remove(temp_file);
                return 1;
            }
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Filesystem error: " << e.what() << std::endl;
            if (file_exists(temp_file)) {
                std::filesystem::remove(temp_file);
            }
            return 1;
        }
    } else if (!file_exists(renamed_original2)) {
        std::cerr << "Neither original2 nor renamed_original2 file exists." << std::endl;
        return 1;
    }

    if (file_exists(original2)) {
        std::cout << "original2 file exists. Renaming to " << renamed_original3 << std::endl;
        if (std::rename((original2).c_str(), renamed_original3.c_str()) != 0) {
            std::perror("Error renaming original2 file");
            return 1;
        }
    } else if (!file_exists(renamed_original3)) {
        std::cerr << "Neither original2 nor renamed_original3 file exists." << std::endl;
        return 1;
    }

    std::ofstream outFile2(renamed_original2, std::ios::app);
    if (!outFile2) {
        std::cerr << "Failed to open " << renamed_original2 << " for appending." << std::endl;
        return 1;
    }
    outFile2 << std::fixed << std::setprecision(6);

    std::ofstream outFile3(renamed_original3, std::ios::app);
    if (!outFile3) {
        std::cerr << "Failed to open " << renamed_original3 << " for appending." << std::endl;
        return 1;
    }
    outFile3 << std::fixed << std::setprecision(6);

    std::ifstream file(renamed_original1);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << renamed_original1 << std::endl;
        return 1;
    }

    std::vector<std::vector<double>> data;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> row(10);
        bool valid = true;
        for (int i = 0; i < 10; ++i) {
            if (!(iss >> row[i]) || std::isnan(row[i]) || std::isinf(row[i])) {
                std::cerr << "Error: Invalid data format in line: " << line << std::endl;
                valid = false;
                break;
            }
        }
        if (valid) {
            data.push_back(row);
        }
    }
    file.close();
    if (data.empty()) {
        std::cerr << "Error: No valid data read from " << renamed_original1 << std::endl;
        return 1;
    }

    const int NUM_COLUMNS = 5;
    std::vector<std::vector<double>> bin_means(NUM_COLUMNS), bin_stderr(NUM_COLUMNS);
    std::vector<std::vector<std::vector<double>>> saved_bins_x(NUM_COLUMNS), saved_bins_y(NUM_COLUMNS);
    std::vector<std::vector<double>> saved_bin_edges(NUM_COLUMNS);

    const int min_bins = 5;
    for (int ind = 0; ind < NUM_COLUMNS; ++ind) {
        std::vector<double> x_data, y_data;
        for (const auto& row : data) {
            x_data.push_back(row[4]);
            y_data.push_back(row[5 + ind]);
        }

        std::vector<std::vector<double>> bins_x, bins_y;
        std::vector<double> bin_edges;

        if (ind == 0 || ind == 3) {
            // Logistic fit for ind == 0 (columns 4, 8) or ind == 3
            std::vector<double> params = {1.0, 0.1, 30.0};
            ceres_curve_fit(params, x_data, y_data, 0.5, 1.5);
            double L = params[0], k = params[1], x0 = params[2];
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "Processing Column Pair (4, " << 5 + ind << "):\n";
            std::cout << "Estimated parameters: L=" << L << ", k=" << k << ", x0=" << x0 << std::endl;

            if (ind == 0) {
                // Adaptive binning for ind == 0 (columns 4, 8) matching Python
                if (x_data.empty()) {
                    std::cerr << "Error: No data for ind == 0" << std::endl;
                    return 1;
                }

                // Sort data by x_data
                std::vector<size_t> indices(x_data.size());
                for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;
                std::sort(indices.begin(), indices.end(), [&x_data](size_t a, size_t b) { return x_data[a] < x_data[b]; });
                std::vector<double> x_sorted(x_data.size()), y_sorted(y_data.size());
                for (size_t i = 0; i < indices.size(); ++i) {
                    x_sorted[i] = x_data[indices[i]];
                    y_sorted[i] = y_data[indices[i]];
                }

                // Debug: Output first few sorted points
                std::cout << "First 5 sorted x_data: ";
                for (size_t i = 0; i < std::min<size_t>(5, x_sorted.size()); ++i) {
                    std::cout << std::fixed << std::setprecision(6) << x_sorted[i] << " ";
                }
                std::cout << "\nFirst 5 sorted y_data: ";
                for (size_t i = 0; i < std::min<size_t>(5, y_sorted.size()); ++i) {
                    std::cout << std::fixed << std::setprecision(6) << y_sorted[i] << " ";
                }
                std::cout << std::endl;

                const int num_bins = 8;
                const int min_points_per_bin = 16;
                bin_edges = {x_sorted[0] - 0.01};
                std::vector<double> current_bin_x, current_bin_y;

                for (size_t i = 0; i < x_sorted.size(); ++i) {
                    current_bin_x.push_back(x_sorted[i]);
                    current_bin_y.push_back(y_sorted[i]);
                    if (current_bin_x.size() >= static_cast<size_t>(min_points_per_bin)) {
                        double mean_y = 0.0;
                        for (double val : current_bin_y) mean_y += val;
                        mean_y /= current_bin_y.size();
                        double mid_x = (current_bin_x[0] + current_bin_x.back()) / 2.0;
                        double expected_y = logistic(mid_x, L, k, x0);
                        double mse = (mean_y - expected_y) * (mean_y - expected_y);
                        std::cout << "Bin attempt at index " << i << ": MSE = " << std::setprecision(6) << mse << std::endl;
                        if (mse < 0.0005 || i == x_sorted.size() - 1 || bins_x.size() == static_cast<size_t>(num_bins - 1)) {
                            bins_x.push_back(current_bin_x);
                            bins_y.push_back(current_bin_y);
                            bin_edges.push_back(current_bin_x.back() + 0.01);
                            current_bin_x.clear();
                            current_bin_y.clear();
                        }
                    }
                }
                if (!current_bin_x.empty() && current_bin_x.size() >= 8) {
                    bins_x.push_back(current_bin_x);
                    bins_y.push_back(current_bin_y);
                    bin_edges.push_back(current_bin_x.back() + 0.01);
                }

                // Ensure at least min_bins bins
                while (bins_x.size() < static_cast<size_t>(min_bins)) {
                    size_t largest_bin_idx = 0;
                    size_t max_size = 0;
                    for (size_t i = 0; i < bins_x.size(); ++i) {
                        if (bins_x[i].size() > max_size) {
                            max_size = bins_x[i].size();
                            largest_bin_idx = i;
                        }
                    }
                    if (max_size < 2 * min_points_per_bin) break;
                    size_t split_index = max_size / 2;
                    if (split_index < min_points_per_bin || (max_size - split_index) < min_points_per_bin) break;

                    std::vector<double> left_x(bins_x[largest_bin_idx].begin(), bins_x[largest_bin_idx].begin() + split_index);
                    std::vector<double> left_y(bins_y[largest_bin_idx].begin(), bins_y[largest_bin_idx].begin() + split_index);
                    std::vector<double> right_x(bins_x[largest_bin_idx].begin() + split_index, bins_x[largest_bin_idx].end());
                    std::vector<double> right_y(bins_y[largest_bin_idx].begin() + split_index, bins_y[largest_bin_idx].end());
                    double new_edge = (left_x.back() + right_x.front()) / 2.0;
                    bin_edges.insert(bin_edges.begin() + largest_bin_idx + 1, new_edge);
                    bins_x[largest_bin_idx] = left_x;
                    bins_y[largest_bin_idx] = left_y;
                    bins_x.insert(bins_x.begin() + largest_bin_idx + 1, right_x);
                    bins_y.insert(bins_y.begin() + largest_bin_idx + 1, right_y);
                }

                // Warn if requirements not met
                if (bins_x.size() < static_cast<size_t>(min_bins)) {
                    std::cerr << "Warning: Only " << bins_x.size() << " bins created for ind=" << ind << ".\n";
                }

                // Debug: Print bin statistics
                std::cout << "Bin Edges for ind=" << ind << ": ";
                for (double edge : bin_edges) std::cout << std::fixed << std::setprecision(6) << edge << " ";
                std::cout << std::endl;
                for (size_t i = 0; i < bins_y.size(); ++i) {
                    if (bins_y[i].size() < min_points_per_bin) {
                        std::cerr << "Warning: Bin " << i << " for ind=" << ind << " has only " << bins_y[i].size() << " points" << std::endl;
                    }
                    double mean = 0.0;
                    for (double val : bins_y[i]) mean += val;
                    mean /= bins_y[i].size();
                    double variance = 0.0;
                    for (double val : bins_y[i]) variance += (val - mean) * (val - mean);
                    double stderr = bins_y[i].size() > 1 ? std::sqrt(variance / (bins_y[i].size() - 1)) / std::sqrt(bins_y[i].size()) : 0.0;
                    std::cout << "Bin " << std::fixed << std::setprecision(6) << bin_edges[i] << "-" << bin_edges[i + 1]
                              << ": Mean = " << std::setprecision(6) << mean << ", Std Error = " << stderr
                              << ", Points = " << bins_y[i].size() << std::endl;
                }

                // Save bin data for external plotting
                std::ofstream plot_data("bin_data_ind0.txt");
                for (size_t i = 0; i < bins_y.size(); ++i) {
                    double mean = 0.0;
                    for (double val : bins_y[i]) mean += val;
                    mean /= bins_y[i].size();
                    double variance = 0.0;
                    for (double val : bins_y[i]) variance += (val - mean) * (val - mean);
                    double stderr = bins_y[i].size() > 1 ? std::sqrt(variance / (bins_y[i].size() - 1)) / std::sqrt(bins_y[i].size()) : 0.0;
                    double center = (bin_edges[i] + bin_edges[i + 1]) / 2.0;
                    plot_data << std::fixed << std::setprecision(6) << center << "\t" << mean << "\t" << stderr << "\n";
                }
                plot_data.close();
            } else {
                // Original binning for ind == 3
                const int num_bins = 8;
                const int min_points_per_bin = 15;
                bin_edges = {0.0};
                std::vector<double> current_bin_x, current_bin_y;
                size_t start_idx = 0;

                for (size_t i = 0; i < x_data.size(); ++i) {
                    current_bin_x.push_back(x_data[i]);
                    current_bin_y.push_back(y_data[i]);
                    bool force_bin = (bins_x.size() < min_bins && i < x_data.size() - 1 &&
                                      current_bin_x.size() >= static_cast<size_t>(min_points_per_bin));
                    if (current_bin_x.size() >= static_cast<size_t>(min_points_per_bin)) {
                        double mean_y = 0.0;
                        for (double val : current_bin_y) mean_y += val;
                        mean_y /= current_bin_y.size();
                        double mid_x = (current_bin_x[0] + current_bin_x.back()) / 2.0;
                        double expected_y = logistic(mid_x, L, k, x0);
                        double mse = (mean_y - expected_y) * (mean_y - expected_y);
                        if (mse < 0.0005 || i == x_data.size() - 1 || bins_x.size() == static_cast<size_t>(num_bins - 1) || force_bin) {
                            bins_x.push_back(current_bin_x);
                            bins_y.push_back(current_bin_y);
                            bin_edges.push_back(static_cast<double>(i + 1));
                            current_bin_x.clear();
                            current_bin_y.clear();
                            start_idx = i + 1;
                        }
                    }
                }
                if (!current_bin_x.empty()) {
                    bins_x.push_back(current_bin_x);
                    bins_y.push_back(current_bin_y);
                    bin_edges.push_back(static_cast<double>(x_data.size()));
                }

                while (bins_x.size() < static_cast<size_t>(min_bins) && !x_data.empty()) {
                    size_t largest_bin_idx = 0;
                    size_t max_size = bins_x[0].size();
                    for (size_t i = 1; i < bins_x.size(); ++i) {
                        if (bins_x[i].size() > max_size) {
                            max_size = bins_x[i].size();
                            largest_bin_idx = i;
                        }
                    }
                    if (max_size < static_cast<size_t>(min_points_per_bin)) break;
                    size_t split_point = bins_x[largest_bin_idx].size() / 2;
                    if (split_point < static_cast<size_t>(min_points_per_bin)) break;

                    std::vector<double> new_bin_x(bins_x[largest_bin_idx].begin() + split_point, bins_x[largest_bin_idx].end());
                    std::vector<double> new_bin_y(bins_y[largest_bin_idx].begin() + split_point, bins_y[largest_bin_idx].end());
                    bins_x[largest_bin_idx].resize(split_point);
                    bins_y[largest_bin_idx].resize(split_point);
                    bins_x.insert(bins_x.begin() + largest_bin_idx + 1, new_bin_x);
                    bins_y.insert(bins_y.begin() + largest_bin_idx + 1, new_bin_y);
                    bin_edges.insert(bin_edges.begin() + largest_bin_idx + 1, bin_edges[largest_bin_idx] + split_point);
                }
            }

            saved_bins_x[ind] = bins_x;
            saved_bins_y[ind] = bins_y;
            saved_bin_edges[ind] = bin_edges;
        } else if (ind == 1 || ind == 2) {
            if (saved_bins_x[0].empty()) {
                std::cerr << "Error: Bins for ind=0 not initialized for ind=" << ind << std::endl;
                return 1;
            }
            bins_x = saved_bins_x[0];
            bins_y.clear();
            bin_edges = saved_bin_edges[0];
            size_t start = 0;
            for (size_t i = 1; i < bin_edges.size(); ++i) {
                size_t idx = static_cast<size_t>(bin_edges[i]);
                if (idx <= x_data.size()) {
                    std::vector<double> bin_y(y_data.begin() + start, y_data.begin() + idx);
                    bins_y.push_back(bin_y);
                    start = idx;
                }
            }
            saved_bins_x[ind] = bins_x;
            saved_bins_y[ind] = bins_y;
            saved_bin_edges[ind] = bin_edges;
        } else if (ind == 4) {
            if (saved_bins_x[3].empty()) {
                std::cerr << "Error: Bins for ind=3 not initialized for ind=4" << std::endl;
                return 1;
            }
            bins_x = saved_bins_x[3];
            bins_y.clear();
            bin_edges = saved_bin_edges[3];
            size_t start = 0;
            for (size_t i = 1; i < bin_edges.size(); ++i) {
                size_t idx = static_cast<size_t>(bin_edges[i]);
                if (idx <= x_data.size()) {
                    std::vector<double> bin_y(y_data.begin() + start, y_data.begin() + idx);
                    bins_y.push_back(bin_y);
                    start = idx;
                }
            }
            saved_bins_x[ind] = bins_x;
            saved_bins_y[ind] = bins_y;
            saved_bin_edges[ind] = bin_edges;
        }

        for (size_t i = 0; i < bins_y.size(); ++i) {
            if (bins_y[i].empty()) continue;
            double mean = 0.0;
            for (double val : bins_y[i]) mean += val;
            mean /= bins_y[i].size();
            bin_means[ind].push_back(mean);

            double variance = 0.0;
            for (double val : bins_y[i]) variance += (val - mean) * (val - mean);
            double stderr = bins_y[i].size() > 1 ? std::sqrt(variance / (bins_y[i].size() - 1)) / std::sqrt(bins_y[i].size()) : 0.0;
            bin_stderr[ind].push_back(stderr);
        }
    }

    // Write to outFile2 with bounds checking
    for (size_t j = 0; j < saved_bins_x[0].size() && j + 1 < saved_bin_edges[0].size(); ++j) {
        outFile2 << THREADS << "\t" << STAGES << "\t" << SAMPLES_PER_THREAD << "\t" 
                 << SWEEPS << "\t" << DA_SWEEPS << "\t" 
                 << (saved_bin_edges[0][j + 1] + saved_bin_edges[0][j])/2 << "\t"
                 << saved_bin_edges[0][j] << "\t"
                 << saved_bin_edges[0][j + 1];
        for (int k = 0; k < 3; ++k) {
            outFile2 << "\t" << (j < bin_means[k].size() ? bin_means[k][j] : 0.0) 
                     << "\t" << (j < bin_stderr[k].size() ? bin_stderr[k][j] : 0.0);
        }
        outFile2 << std::endl;
    }
    outFile2.close();

    // Write to outFile3 with bounds checking
    for (size_t j = 0; j < saved_bins_x[3].size() && j + 1 < saved_bin_edges[3].size(); ++j) {
        outFile3 << THREADS << "\t" << STAGES << "\t" << SAMPLES_PER_THREAD << "\t" 
                 << SWEEPS << "\t" << DA_SWEEPS << "\t" 
                 << (saved_bin_edges[3][j + 1] + saved_bin_edges[3][j])/2 << "\t"
                 << saved_bin_edges[3][j] << "\t"
                 << saved_bin_edges[3][j + 1];
        for (int k = 3; k < 5; ++k) {
            outFile3 << "\t" << (j < bin_means[k].size() ? bin_means[k][j] : 0.0) 
                     << "\t" << (j < bin_stderr[k].size() ? bin_stderr[k][j] : 0.0);
        }
        outFile3 << std::endl;
    }
    outFile3.close();

    return 0;
}