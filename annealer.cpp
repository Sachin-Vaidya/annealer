#include <iomanip>
#include <iostream>
#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <thread>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/histogram/utility/clopper_pearson_interval.hpp>
#include <boost/histogram/fwd.hpp>
#include <boost/math/distributions/beta.hpp>
#include <set>

#include "annealer.hh"
#include "vertexing.hh"
#include "detanneal.hh"

ftype Temp_0;
ftype Temp_f;

int THREADS;
int STAGES;
int SAMPLES_PER_THREAD;
int SWEEPS;
int DA_SWEEPS;

ftype Temperature_SA, Temperature_DA;
ftype SA_timer_start;
ftype SA_timer_stop;

double hotFlipProbability = 0.5;
double coldFlipProbability = 0.01;

std::vector<std::vector<int>> CorrectGivenValidSolutionOrNot, ValidSolutionOrNot;
int CorrectGivenValid_SA_Solutions, Valid_SA_Solutions;
int Correct_DA_Solutions, Valid_DA_Solutions;
std::vector<ftype> SATimePerThreadsAnneal;

std::set<std::set<int>> valueSetsInMap_truth;
std::vector<int> int_truth;
vector<vector<int>> sorted_truth_clusters;

std::vector<std::vector<int>> convertAndSortClusters(const std::set<std::set<int>>& cluster_sets) {
    std::vector<std::vector<int>> clusters_vec;
    for (const auto& cluster : cluster_sets) {
        std::vector<int> sorted_cluster(cluster.begin(), cluster.end());
        std::sort(sorted_cluster.begin(), sorted_cluster.end());
        clusters_vec.push_back(std::move(sorted_cluster));
    }
    std::sort(clusters_vec.begin(), clusters_vec.end());
    return clusters_vec;
}

// Function to check if a file exists
bool file_exists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

ostream& operator << (ostream& os, const solution_t& x) {
    for (auto xi : x) os << xi << ' ';
    return os;
}

void reassignUniqueIDs(std::vector<int>& v) {
    std::unordered_map<int, int> valueToID;
    int nextID = 0;

    for (int& x : v) {
        auto it = valueToID.find(x);
        if (it == valueToID.end()) {
            valueToID[x] = nextID;
            x = nextID;
            ++nextID;
        } else {
            x = it->second;
        }
    }
}

// ostream& operator << (ostream& os, const qubo_t& Q) {
//     for (const auto& entry : Q) {
//         os << '[' << entry.first.first << ' ' << entry.first.second << "] : " << entry.second << endl;
//     }
//     return os;
// }

struct settings {
    int max_iter;
    ftype T_0;
    ftype T_f;
    scheduler_t temp_scheduler;
    unsigned seed;
    problem_context context;
    bool dolog = true;
    // bool OTF = true;
    int da_sweeps;
};

struct result {
    solution_t solution;
    ftype energy;
};

inline int qubo_size(const qubo_t& Q) {
    int n = 0;
    for (const auto& entry : Q) {
        n = max(n, entry.first.first);
    }
    return n + 1; // 0 indexed
}

QUBO::QUBO(event_t& event) {
    this->OTF = true;
    this->event = &event;
    this->max_D = get_max_D(event);
}

QUBO::QUBO(qubo_t& Q) {
    this->OTF = false;

    n = qubo_size(Q);

    offset.resize(n + 1, 0); // n+1 to store end too.

    for (const auto& [idx, val] : Q) {
        offset[idx.first+1]++;
        if (idx.first == idx.second) continue;
        offset[idx.second+1]++;
    }

    // currently, offset[i] = size of row i-1. prefix sum to get the start of row i

    for (int i = 1; i < offset.size(); i++) {
        offset[i] += offset[i - 1];
    }

    affectedby_flat.resize(offset[n]); // total size of the flat list

    vector<int> insert_pos = offset;

    for (const auto& [idx, val] : Q) {
        int i = idx.first;
        int j = idx.second;
        affectedby_flat[insert_pos[i]++] = {j, val};
        if (i == j) continue;
        affectedby_flat[insert_pos[j]++] = {i, val};
    }

    //std::cout << "conversion done\n";
    //std::cout << "size of flat list = " << affectedby_flat.size() << "\n";
}

const std::vector<std::pair<int, double>>& QUBO::getAffectedByFlat() const {
    return affectedby_flat;
}

const std::vector<int>& QUBO::getOffset() const {
    return offset;
}

ftype QUBO::evaluate(const solution_t& x) const {
    if (this->OTF) {
        return evaluate_full_OTF(x, *this->event, this->max_D);
    }

    ftype value = 0.0;
    for (int i = 0; i < x.size(); i++) {
        if (x[i]) {
            for (int k = offset[i]; k < offset[i + 1]; k++) {
                const auto& [j, Q_ij] = affectedby_flat[k];
                if (x[j] && j >= i) value += Q_ij; // only count each pair once
            }
        }
    }
    return value;
}

ftype QUBO::evaluateDiff(const solution_t& x, int flip_idx) const {
    if (this->OTF) {
        return evaluate_diff_on_the_fly(x, *this->event, flip_idx, this->max_D);
    }

    ftype diff = 0.0; // first, find what would be the value if this bit was on
    // #pragma clang loop vectorize_width(2)
    // #pragma clang loop interleave_count(2)
    for (int k = offset[flip_idx]; k < offset[flip_idx + 1]; k++) {
        const auto& [j, Q_ij] = affectedby_flat[k];
        
        //diff += Q_ij*x[j];
        diff += (j == flip_idx) ? Q_ij : Q_ij*x[j];
    }
    //return diff * (1 - 2 * x[flip_idx]);
    return x[flip_idx] ? -diff : diff; // if on, turn off. if off, turn on.
}

// problem specific anneal


//=================================== SA without and with sweeps from Sachin's code =========================================//

//only one bit flip
/*
result sim_anneal(const QUBO& Q, const settings s, const solution_t init_guess = {}, int num_stage=1) { // intentionally get copy of settings
    mt19937 gen(s.seed);
    uniform_real_distribution<> dis(0.0, 1.0);

    int nT = s.context.event.nT;
    int nV = s.context.event.nV;

    uniform_int_distribution<> t_dis(0, nT - 1);
    uniform_int_distribution<> v_dis(0, nV - 1);

    solution_t x(nT * nV, 0);

    auto bit_idx = [nT](int t, int v) {
        return t + nT * v;
    };

    // helper.
    vector<int> track_to_vertex(nT, -1);

    // for each track, assign random vertex
    for (int t = 0; t < nT; t++) {
        int v = v_dis(gen);
        track_to_vertex[t] = v;
        x[bit_idx(t, v)] = 1;
    }

    if (!init_guess.empty()) {
        x = init_guess;

        for (int t = 0; t < nT; t++) {
            for (int v = 0; v < nV; v++) {
                if (x[bit_idx(t, v)]) {
                    track_to_vertex[t] = v;
                    // break;
                }
            }
        }

        //std::cout << "Using init guess\n";
    }

    ftype f_x = Q.evaluate(x);

    solution_t best_x = x;
    ftype best_f_x = f_x;

    ftype T = s.T_0;
    for (int iter = 0; iter < s.max_iter; iter++) {
        if (iter % 1000 == 0 && s.dolog) {
            //std::cout << "Iter: " << iter << " Energy: " << f_x << " T: " << T << '\n';
        }

        int t = t_dis(gen);
        int v = v_dis(gen);
        int bit = bit_idx(t, v);

        solution_t x_prime = x;
        ftype delta = Q.evaluateDiff(x_prime, bit);
        x_prime[bit] = 1 - x[bit]; // toggle bit

        ftype f_x_prime = f_x + delta;

        if (f_x_prime < f_x || (T != 0 && dis(gen) < exp((f_x - f_x_prime) / T))) {
            x = x_prime;
            f_x = f_x_prime;

            // Optionally update track_to_vertex[t], but careful:
            // Since multiple bits per track can be set, maybe just set to -1
            //track_to_vertex[t] = -1; // or recompute after loop
        }

        if (f_x < best_f_x) {
            best_x = x;
            best_f_x = f_x;
        }

        T = s.temp_scheduler(s.T_0, s.T_f, T, iter, s.max_iter, num_stage);
        Temperature_SA = T;

        //std::cout << "\nTemperature_SA calculation verification: " << num_stage << ", " << T << ", " << Temperature_SA << endl;
    }

    // re-evaluate best solution to eliminate float point errors
    best_f_x = Q.evaluate(best_x);

    return {best_x, best_f_x};
}
*/

//Sweep through entire bit string and flip them
/*
result sim_anneal(const QUBO& Q, const settings s, const solution_t init_guess = {}, int num_stage = 1) {
    mt19937 gen(s.seed);
    uniform_real_distribution<> dis(0.0, 1.0);

    int nT = s.context.event.nT;
    int nV = s.context.event.nV;

    // track-major index function
    auto bit_idx = [nV](int t, int v) {
        return nV * t + v;
    };

    solution_t x(nT * nV, 0);

    // Initialize random or with init_guess
    if (!init_guess.empty()) {
        x = init_guess;
    } else {
        for (int t = 0; t < nT; ++t) {
            int v = gen() % nV;
            x[bit_idx(t, v)] = 1;
        }
    }

    ftype f_x = Q.evaluate(x);
    solution_t best_x = x;
    ftype best_f_x = f_x;

    ftype T = s.T_0;

    for (int iter = 0; iter < s.max_iter; ++iter) {
        for (int t = 0; t < nT; ++t) {
            for (int v = 0; v < nV; ++v) {
                int bit = bit_idx(t, v);

                // Try flipping this bit in-place
                ftype delta = Q.evaluateDiff(x, bit);
                ftype f_x_prime = f_x + delta;

                bool accept = (f_x_prime < f_x) || (T > 0 && dis(gen) < exp((f_x - f_x_prime) / T));

                if (accept) {
                    x[bit] = 1 - x[bit]; // Flip bit
                    f_x = f_x_prime;

                    if (f_x < best_f_x) {
                        best_x = x;
                        best_f_x = f_x;
                    }
                }
                // No need to revert the bit — flip only if move is accepted
            }
        }

        T = s.temp_scheduler(s.T_0, s.T_f, T, iter, s.max_iter, num_stage);
        Temperature_SA = T;
    }

    best_f_x = Q.evaluate(best_x); // re-evaluate in case of rounding
    return {best_x, best_f_x};
}
*/

//=================================== ====================================== =========================================//


/*
//=================================== SA without and with sweeps from Ishan's code =========================================//

//SA without sweep

result sim_anneal(const QUBO& Q, const settings s, const solution_t init_guess = {}, int num_stage = 1) { // intentionally get copy of settings
    // int n = qubo_size(Q);

    mt19937 gen(s.seed);
    uniform_real_distribution<> dis(0.0, 1.0);
    uniform_int_distribution<> flip_dis(0, s.context.event.nT*s.context.event.nV - 1);

    solution_t x(s.context.event.nT*s.context.event.nV, 0);
    for (auto &xi : x) {
        xi = dis(gen) < 0.5 ? 0 : 1;
    }

    if (!init_guess.empty()) x = init_guess;

    ftype f_x = Q.evaluate(x);

    solution_t best_x = x;
    ftype best_f_x = f_x;

    ftype T = s.T_0;
    for (int iter = 0; iter < s.max_iter; iter++) {
        if (iter % 10000 == 0 && s.dolog) {
            //std::cout << "Iter: " << iter << " Energy: " << best_f_x << " T: " << T << '\n';
        }

        solution_t x_prime = x;
        int flip_idx = flip_dis(gen);

        ftype f_x_prime = Q.evaluateDiff(x, flip_idx) + f_x;

        x_prime[flip_idx] = !x_prime[flip_idx];

        if (f_x_prime < f_x || dis(gen) < exp((f_x - f_x_prime) / T)) {
            x = x_prime;
            f_x = f_x_prime;
        }

        if (f_x < best_f_x) {
            best_x = x;
            best_f_x = f_x;
        }

        // Lower temperature after full sweep
        T = s.temp_scheduler(s.T_0, s.T_f, T, iter, s.max_iter, num_stage);
        Temperature_SA = T;

        //std::cout << "\nTemperature_SA calculation verification: " << num_stage << ", " << T << ", " << Temperature_SA << endl;
    }

    best_f_x = Q.evaluate(best_x); // cleanup
    return {best_x, best_f_x};
}
*/

//SA with sweep

result sim_anneal(const QUBO& Q, const settings s, const solution_t init_guess = {}, int num_stage = 1) {
    //mt19937 gen(42);
    mt19937 gen(s.seed);
    uniform_real_distribution<> dis(0.0, 1.0);
    
    solution_t x(s.context.event.nT*s.context.event.nV, 0);
    if (!init_guess.empty()) {
        x = init_guess;
    }
    else {
        for (auto &xi : x) {
            xi = ftype(dis(gen)/1.0) < 0.5 ? 0 : 1;
        }
    }

    ftype f_x = Q.evaluate(x);
    solution_t best_x = x;
    ftype best_f_x = f_x;

    ftype T = s.T_0;

    //std::cout << "Sweeps: " << s.max_iter << ", bitstring size: " << static_cast<int>(x.size()) << endl;

    for (int iter = 0; iter < s.max_iter; iter++) {

        //std::cout << "Inverse Temperature: " << 1./T << endl;
        for (int i = 0; i < static_cast<int>(x.size()); i++) {
            
            ftype f_x_diff = Q.evaluateDiff(x, i);
            
            ftype  f_x_prime = f_x + f_x_diff;
            //std::cout << f_x << ", " << f_x_diff << ", " << f_x_prime << std::endl;
            
            if (f_x_diff < 0 || ftype(dis(gen)/1.0) < exp(-f_x_diff / T)) {
                x[i] = !x[i];  // Flip bit
                f_x = f_x_prime;
            }
            
            /*// Create a copy of x
            solution_t x_prime = x;
            // Flip the i-th bit in the copy
            x_prime[i] = !x[i];
            
            ftype f_x_prime = Q.evaluate(x_prime);
            
            if (f_x_prime < f_x || ftype(dis(gen)/1.0) < exp(-(f_x_prime - f_x) / T)) {
                x[i] = !x[i];  // Flip bit
                f_x = f_x_prime;
            }*/

            /*if (f_x < best_f_x) {
                best_x = x;
                best_f_x = f_x;
            }*/
        }

        T = s.temp_scheduler(s.T_0, s.T_f, T, iter, s.max_iter, num_stage);
        Temperature_SA = T;

        //std::cout << "\nTemperature_SA calculation verification: " << num_stage << ", " << T << ", " << Temperature_SA << endl;
        //std::cout << "Inverse Temperature: " << 1./T << endl;
    }
    
    best_x = x;
    best_f_x = Q.evaluate(best_x); // cleanup
    
    return {best_x, best_f_x};
}

//SA with sweep with better optimization
/*result sim_anneal(const QUBO& Q, const settings s, const solution_t init_guess = {}, int num_stage = 1) {
    mt19937 gen(s.seed);
    uniform_real_distribution<> dis(0.0, 1.0);

    solution_t x;
    if (!init_guess.empty()) {
        x = init_guess;
    } else {
        x.resize(s.context.event.nT * s.context.event.nV);
        for (auto &xi : x)
            xi = dis(gen) < 0.5 ? 1 : 0;
    }

    ftype f_x = Q.evaluate(x);
    solution_t best_x = x;
    ftype best_f_x = f_x;

    std::vector<ftype> T_schedule(s.max_iter);
    for (int i = 0; i < s.max_iter; ++i)
        T_schedule[i] = s.temp_scheduler(s.T_0, s.T_f, s.T_0, i, s.max_iter, num_stage);

    for (int iter = 0; iter < s.max_iter; ++iter) {
        ftype T = T_schedule[iter];
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            ftype f_x_prime = Q.evaluateDiff(x, i) + f_x;

            if (f_x_prime < f_x || dis(gen) < exp((f_x - f_x_prime) / T)) {
                x[i] = !x[i];
                f_x = f_x_prime;

                if (f_x < best_f_x) {
                    best_x = x;
                    best_f_x = f_x;
                }
            }
        }
        Temperature_SA = T;
    }

    return {best_x, best_f_x};
}*/

//=================================== ====================================== =========================================//


void assert_lower_triangular(const qubo_t& Q) {
    for (const auto& entry : Q) {
        if (entry.first.first < entry.first.second) {
            cerr << "Error: QUBO is not in lower triangular form." << endl;
            exit(1);
        }
    }
}

void assert_upper_triangular(const qubo_t& Q) {
    for (const auto& entry : Q) {
        if (entry.first.first > entry.first.second) {
            cerr << "Error: QUBO is not in upper triangular form." << endl;
            exit(1);
        }
    }
}

bool is_valid_solution(const auto& vert, int nT, int nV, bool flag_cluster_major, vector<int>& int_vert_assignment) {
    int_vert_assignment.assign(nT, -1);

    for (int l = 0; l < nT; ++l) {
        int sum = 0, selected_p = -1;
        for (int p = 0; p < nV; ++p) {
            int val = flag_cluster_major ? vert[p * nT + l] : vert[p + l * nV];
            if (val == 1) {
                selected_p = p;
                sum += val;
            }
        }
        if (sum != 1) {
            return false;  // Invalid if not exactly one selected
        } else {
            int_vert_assignment[l] = selected_p;
        }
    }
    return true;
}

bool is_correct_solution(const vector<int>& int_vert_assignment, const vector<vector<int>>& sorted_truth_clusters, int nV) {
    vector<vector<int>> clusters(nV);
    for (int t = 0; t < int_vert_assignment.size(); ++t) {
        clusters[int_vert_assignment[t]].push_back(t);
    }

    // Remove empty clusters, sort internally and externally
    vector<vector<int>> filtered_clusters;
    for (auto& c : clusters) {
        if (!c.empty()) {
            sort(c.begin(), c.end());
            filtered_clusters.push_back(std::move(c));
        }
    }
    sort(filtered_clusters.begin(), filtered_clusters.end());

    return filtered_clusters == sorted_truth_clusters;
}

vector<result> multithreaded_sim_anneal(const QUBO& Q, const settings s, int num_threads, int samples_per_thread, const vector<solution_t> init_guess, string filename, int num_stage) {
    vector<result> results(num_threads * samples_per_thread);
    vector<solution_t> local_init_guess;

    // Handle init_guess fallback logic
    if (init_guess.size() != num_threads) {
        if (!init_guess.empty()) {
            local_init_guess = vector<solution_t>(num_threads, init_guess[0]);
        } else {
            local_init_guess = vector<solution_t>(num_threads, solution_t{});
        }
    }

    omp_set_num_threads(num_threads);

    event_t event = loadTracks(filename);
    ftype ground = ground_state(Q, s.context.event);

    CorrectGivenValid_SA_Solutions = 0;
    Valid_SA_Solutions = 0;

    const int nT = s.context.event.nT;
    const int nV = s.context.event.nV;

    CorrectGivenValidSolutionOrNot.assign(num_threads, vector<int>(samples_per_thread, 0));
    ValidSolutionOrNot.assign(num_threads, vector<int>(samples_per_thread, 0));

    for (int j = 0; j < samples_per_thread; j++) {
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int index = tid * samples_per_thread + j;

            settings s_copy = s;
            s_copy.seed += tid * samples_per_thread + j;

            result& r = results[index];
            r = sim_anneal(Q, s_copy, local_init_guess[tid], num_stage);

            const auto& vert = r.solution;
            vector<int> int_vert_assignment(nT, -1);
            bool valid = true;

            // Assignment check
            for (int l = 0; l < nT && valid; ++l) {
                int sum = 0, selected_p = -1;
                for (int p = 0; p < nV; ++p) {
                    int val = flag_cluster_major ? vert[p * nT + l] : vert[p + l * nV];
                    if (val == 1)
                    {   
                        sum += val;
                        selected_p = p;
                    }
                }
                if (sum != 1) {
                    valid = false;
                } else {
                    int_vert_assignment[l] = selected_p;
                }
            }

            if (valid) {
                ValidSolutionOrNot[tid][j] = 1;

                // Collect clusters via vector instead of map
                vector<vector<int>> clusters(nV);
                for (int t = 0; t < nT; ++t) {
                    clusters[int_vert_assignment[t]].push_back(t);
                }

                // Filter empty, sort inner and outer vectors
                vector<vector<int>> filtered_clusters;
                for (auto& c : clusters) {
                    if (!c.empty()) {
                        std::sort(c.begin(), c.end());
                        filtered_clusters.push_back(std::move(c));
                    }
                }
                std::sort(filtered_clusters.begin(), filtered_clusters.end());

                bool is_correct = (filtered_clusters == sorted_truth_clusters);
                CorrectGivenValidSolutionOrNot[tid][j] = is_correct ? 1 : 0;
            } else {
                ValidSolutionOrNot[tid][j] = 0;
                CorrectGivenValidSolutionOrNot[tid][j] = 0;
            }
        }
    }

    for (int i = 0; i < num_threads; i++) {
        for (int j = 0; j < samples_per_thread; j++) {
            CorrectGivenValid_SA_Solutions += CorrectGivenValidSolutionOrNot[i][j];
            Valid_SA_Solutions += ValidSolutionOrNot[i][j];
        }
    }

    std::sort(results.begin(), results.end(), [](const result& a, const result& b) {
        return a.energy < b.energy;
    });

    return results;
}

/*// returns sorted results of length num_threads * samples_per_thread. first elem is best.
vector<result> multithreaded_sim_anneal(const QUBO& Q, const settings s, int num_threads, int samples_per_thread, const vector<solution_t> init_guess, string filename, int num_stage) {
    vector<result> results(num_threads * samples_per_thread);
    vector<solution_t> local_init_guess;

    // Handle init_guess fallback logic
    if (init_guess.size() != num_threads) {
        if (!init_guess.empty()) {
            local_init_guess = vector<solution_t>(num_threads, init_guess[0]);
        } else {
            local_init_guess = vector<solution_t>(num_threads, solution_t{});
        }
    }

    omp_set_num_threads(num_threads);

    event_t event = loadTracks(filename);
    ftype ground = ground_state(Q, s.context.event);

    CorrectGivenValid_SA_Solutions = 0;
    Valid_SA_Solutions = 0;

    const int nT = s.context.event.nT;
    const int nV = s.context.event.nV;

    CorrectGivenValidSolutionOrNot.assign(num_threads, vector<int>(samples_per_thread, 0));
    ValidSolutionOrNot.assign(num_threads, vector<int>(samples_per_thread, 0));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        for (int j = 0; j < samples_per_thread; j++) 
        {
            int index = tid * samples_per_thread + j;
            //int index = tid * samples_per_thread;
            
            settings s_copy = s;
            s_copy.seed += index;

            result& r = results[index];
            r = sim_anneal(Q, s_copy, local_init_guess[tid], num_stage);

            const auto& vert = r.solution;
            vector<int> int_vert_assignment(nT, -1);
            bool valid = is_valid_solution(vert, nT, nV, flag_cluster_major, int_vert_assignment);
            
            if (valid) {
                ValidSolutionOrNot[tid][j] = 1;
                bool correct = is_correct_solution(int_vert_assignment, sorted_truth_clusters, nV);
                CorrectGivenValidSolutionOrNot[tid][j] = correct ? 1 : 0;
            } else {
                ValidSolutionOrNot[tid][j] = 0;
                CorrectGivenValidSolutionOrNot[tid][j] = 0;
            }
        }
    }

    for (int i = 0; i < num_threads; i++) {
        for (int j = 0; j < samples_per_thread; j++) {
            CorrectGivenValid_SA_Solutions += CorrectGivenValidSolutionOrNot[i][j];
            Valid_SA_Solutions += ValidSolutionOrNot[i][j];
        }
    }

    std::sort(results.begin(), results.end(), [](const result& a, const result& b) {
        return a.energy < b.energy;
    });

    return results;
}*/

/*vector<result> multithreaded_sim_anneal(const QUBO& Q, const settings s, int num_threads, int samples_per_thread, const vector<solution_t> init_guess, string filename, int num_stage) {
    vector<result> results(num_threads * samples_per_thread);
    vector<solution_t> local_init_guess;

    // Handle init_guess fallback logic
    if (init_guess.size() != num_threads) {
        if (!init_guess.empty()) {
            local_init_guess = vector<solution_t>(num_threads, init_guess[0]);
        } else {
            local_init_guess = vector<solution_t>(num_threads, solution_t{});
        }
    }

    omp_set_num_threads(num_threads);

    //event_t event = loadTracks(filename);
    //ftype ground = ground_state(Q, s.context.event);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        for (int j = 0; j < samples_per_thread; j++) 
        {
            int index = tid * samples_per_thread + j;
            
            settings s_copy = s;
            s_copy.seed += index;

            result& r = results[index];
            r = sim_anneal(Q, s_copy, local_init_guess[tid], num_stage);
        }
    }
    
    std::sort(results.begin(), results.end(), [](const result& a, const result& b) {
        return a.energy < b.energy;
    });

    // Write results to JSON file when num_stage == STAGES - 1
    if (num_stage == STAGES - 1) {
        // Replace "/serializedEvents.json" with "/results_flattened_prob.json"
        string output_filename = filename;
        size_t pos = output_filename.find("/serializedEvents.json");
        if (pos != string::npos) {
            output_filename.replace(pos, string("/serializedEvents.json").length(), "/results_flattened_prob.json");
        } else {
            output_filename += "/results_flattened_prob.json"; // Fallback if extension not found
        }

        std::ofstream file(output_filename);
        file << "[\n";
        for (size_t i = 0; i < results.size(); ++i) {
            const auto& r = results[i];
            std::ostringstream oss;
            oss << "[" << r.energy << ",[";
            for (size_t j = 0; j < r.solution.size(); ++j) {
                oss << static_cast<int>(r.solution[j]);
                if (j < r.solution.size() - 1) oss << ",";
                
                //std::cout << static_cast<int>(r.solution[j]) << ",";
            }
            //std::cout << std::endl;
            
            oss << "]]";
            file << "    " << oss.str();
            if (i < results.size() - 1) file << ",";
            file << "\n";
        }
        file << "]";
        file.close();
    }

    return results;
}

std::pair<int, int> compute_solution_efficiencies(const string& filename, const settings& s, bool flag_cluster_major) {
    // Construct output filename
    string input_filename = filename;
    size_t pos = input_filename.find("/serializedEvents.json");
    if (pos != string::npos) {
        input_filename.replace(pos, string("/serializedEvents.json").length(), "/results_flattened_prob.json");
    } else {
        input_filename += "/results_flattened_prob.json"; // Fallback
    }

    // Read JSON file
    std::ifstream file(input_filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << input_filename << std::endl;
        return {0.0, 0.0};
    }

    // Read file content into a string
    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();
    std::string content = buffer.str();

    // Basic JSON parsing
    int Valid_SA_Solutions = 0;
    int CorrectGivenValid_SA_Solutions = 0;
    int total_solutions = 0;

    const int nT = s.context.event.nT;
    const int nV = s.context.event.nV;

    // Remove outer brackets and split entries
    if (content.length() < 2 || content.front() != '[' || content.back() != ']') {
        std::cerr << "Error: Invalid JSON format in " << input_filename << std::endl;
        return {0.0, 0.0};
    }

    // Remove [ and ]
    content = content.substr(1, content.length() - 2);

    // Split by "],\n" to get individual [energy,[solution]] entries
    std::vector<std::string> entries;
    size_t start = 0;
    while (true) {
        size_t end = content.find("],\n", start);
        if (end == std::string::npos) {
            entries.push_back(content.substr(start));
            break;
        }
        entries.push_back(content.substr(start, end - start + 1)); // Include ]
        start = end + 3; // Skip ],\n
    }

    total_solutions = entries.size();
    
    //std::cout << "total_solutions size = " << total_solutions << std::endl;

    for (const auto& entry : entries) {
        // Skip whitespace and [
        size_t pos = entry.find('[');
        if (pos == std::string::npos) continue;
        pos++;

        // Extract energy (before first comma)
        size_t comma_pos = entry.find(',', pos);
        if (comma_pos == std::string::npos) continue;
        std::string energy_str = entry.substr(pos, comma_pos - pos);
        ftype energy;
        try {
            energy = std::stod(energy_str);
        } catch (...) {
            std::cerr << "Error: Invalid energy format in entry: " << entry << std::endl;
            continue;
        }

        // Extract solution (between [ and ])
        pos = entry.find('[', comma_pos);
        if (pos == std::string::npos) continue;
        pos++;
        size_t end_pos = entry.find(']', pos);
        if (end_pos == std::string::npos) continue;
        std::string solution_str = entry.substr(pos, end_pos - pos);

        // Parse solution as vector<int>
        std::vector<int> solution;
        std::stringstream ss(solution_str);
        std::string num;
        while (std::getline(ss, num, ',')) {
            try {
                solution.push_back(std::stoi(num));
            } catch (...) {
                std::cerr << "Error: Invalid solution number in entry: " << std::endl;
                solution.clear();
                break;
            }
        }

        if (solution.empty()) continue;

        // Validity check
        std::vector<int> int_vert_assignment(nT, -1);
        bool valid = is_valid_solution(solution, nT, nV, flag_cluster_major, int_vert_assignment);
        if (valid) {
            Valid_SA_Solutions++;
            // Correctness check
            bool correct = is_correct_solution(int_vert_assignment, sorted_truth_clusters, nV);
            if (correct) {
                CorrectGivenValid_SA_Solutions++;
            }
        }
    }
    
    //std::cout << Valid_SA_Solutions << ", " << CorrectGivenValid_SA_Solutions << std::endl;

    return {Valid_SA_Solutions, CorrectGivenValid_SA_Solutions};
}*/

// replaces bottom half with top half
vector<solution_t> best_effort_unique(const vector<result>& results, int n) {
    if (results.empty()) return {};
    // takes sorted list and replaces bottom half with top half

    vector<solution_t> newresults(n);
    for (int i = 0; i < n; i++) {
        newresults[i] = results[i % (results.size() / 2)].solution;
    }

    return newresults;
}

//----------num_stages is actually num_stages-----------//
vector<result> stage_rejoin_sa(const QUBO& Q, const settings s, int num_threads, int num_stages, int samples_per_thread = 1, string filename = "") {
    settings modified_settings = s;
    modified_settings.max_iter /= num_stages;
    vector<result> results;
    
    for (int i = 0; i < num_stages; i++) {
        results = multithreaded_sim_anneal(Q, modified_settings, num_threads, samples_per_thread, best_effort_unique(results, num_threads), filename, i);
        //std::cout << "Branch " << i << " best energy: " << results[0].energy << '\n';
        //std::cout << "Branch " << i << " worst energy: " << results.back().energy << '\n';
        //std::cout << "Temperature_SA at the end of stage " << i << ": " << Temperature_SA <<endl;
    }
    //std::cout << "Max Iter for stage: " << modified_settings.max_iter <<endl;
    return results;
}

ftype linear_scheduler(ftype T_0, ftype T_f, ftype T, int iter, int max_iter, int num_stage) {
    ftype num_stage_fraction = ftype(num_stage * 1. / STAGES);
    ftype iter_fraction = ftype((iter + 1) * 1. / SWEEPS);
    //return pow(1 / T_0 - (1 / T_0 - 1 / T_f) * (num_stage * (max_iter - 1) + iter) / (max_iter - 1) / STAGES,-1.);
    return pow(1. / T_0 + (1. / T_f - 1. / T_0) * (num_stage_fraction + iter_fraction),-1.);
}

scheduler_t make_geometric_scheduler(ftype alpha) {
    return [alpha](ftype T_0, ftype T_f, ftype T, int iter, int max_iter, int num_stage) {
        return T * alpha;
    };
}

// void trial(solution_t x, const QUBO& Q) {
//     //std::cout << "For solution: " << x << endl;
//     //std::cout << "Energy: " << Q.evaluate(x) << endl << endl;
// }

void present_results(const vector<result>& results, bool show_sols = true, int precision = 5) {
    map<long, map<solution_t, int>> counts; // energy -> solution -> count
    auto d_to_l = [precision](ftype d) {
        return static_cast<long>(round(d * pow(10, precision)));
    };
    auto l_to_d = [precision](long l) {
        return static_cast<ftype>(l) / pow(10, precision);
    };

    for (const auto& r : results) {
        counts[d_to_l(r.energy)][r.solution]++;
    }

    //std::cout << '\n';

    //std::cout << fixed << setprecision(precision);

    //std::cout << "Best energy: " << results[0].energy << '\n';
    //std::cout << "Worst energy: " << results.back().energy << '\n';
    //std::cout << "Best solution: " << results[0].solution << '\n';

    //std::cout << '\n';

    /*for (const auto& [energy, sols] : counts) {
        int total = 0;
        for (const auto& [sol, count] : sols) {
            total += count;
        }
        //std::cout << "Energy: " << l_to_d(energy) << " (" << total << "x)" << '\n';
        if (show_sols) {
            for (const auto& [sol, count] : sols) {
                //std::cout << "\t Sol: " << sol << " (" << count << "x)" << '\n';
            }
        }
    }*/
}

// qubo_t randgen_qubo(int n) {
//     qubo_t Q;
//     random_device rd;
//     unsigned seed = rd();
//     mt19937 gen(seed);
//     uniform_real_distribution<> dis(-10.0, 10.0);
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j <= i; j++) {
//             // Q[{i, j}] = dis(gen);
//             Q.push_back({{i, j}, dis(gen)});
//         }
//     }
//     return Q;
// }

void logfile_append(const string& filename, const string& problem, ftype ARI, ftype energy_diff, ftype mse, ftype ground) {
    ofstream logfile;
    logfile.open(filename, ios::app);
    logfile << problem << ',' << ARI << ',' << energy_diff << ',' << mse << ',' << ground << '\n';
    logfile.close();
}

struct clustering_result {
    string filename;
    vector<int> assignment;
    vector<int> truth;
    vector<ftype> vertices;
    vector<ftype> truth_vertices;
    ftype energy;
    ftype ground;
    ftype valid_eff;
    ftype correct_eff;
    ftype seconds;
    /*vector<int> da_assignment;
    vector<ftype> da_vertices;
    ftype da_energy;*/
};

void json_init(ofstream& jsonfile) {
    jsonfile << "{\n";
    jsonfile << "\"results\": [\n";
}

void json_append(ofstream& jsonfile, const clustering_result& result) {
    jsonfile << "{\n";
    jsonfile << "\"filename\": \"" << result.filename << "\",\n";

    auto write_ints = [&jsonfile](string name, const vector<int>& vec) {
        jsonfile << "\"" << name << "\": [";
        for (int i = 0; i < vec.size(); i++) {
            jsonfile << vec[i];
            if (i != vec.size() - 1) jsonfile << ", ";
        }
        jsonfile << "],\n";
    };

    auto write_floats = [&jsonfile](string name, const vector<ftype>& vec) {
        jsonfile << "\"" << name << "\": [";
        for (int i = 0; i < vec.size(); i++) {
            jsonfile << vec[i];
            if (i != vec.size() - 1) jsonfile << ", ";
        }
        jsonfile << "],\n";
    };

    write_ints("assignment", result.assignment);
    write_ints("truth", result.truth);
    write_floats("vertices", result.vertices);
    write_floats("truth_vertices", result.truth_vertices);

    jsonfile << "\"energy\": " << result.energy << ",\n";
    jsonfile << "\"ground\": " << result.ground << ",\n";
    jsonfile << "\"seconds\": " << result.seconds << "\n";
    jsonfile << "}";
}

void json_close(ofstream& jsonfile) {
    jsonfile << "]\n";
    jsonfile << "}\n";
}

void setTemperatureRange(double hotFlipProbability, double coldFlipProbability, ftype* T_0, ftype* T_f, const QUBO& Q) {
    ftype minCoupling = std::numeric_limits<ftype>::max();
    int minCouplingIndex = -1;  // Index of minCoupling

    ftype maxCouplings = 0.0;
    int maxCouplingIndex = -1;  // Index of maxCouplings

    // Iterate over QUBO terms to find minCoupling
    const auto& affected = Q.getAffectedByFlat();
    for (size_t i = 0; i < affected.size(); ++i) {
        ftype coupling = std::abs(affected[i].second);
        if (coupling > 0.1 && coupling < minCoupling) {
            minCoupling = coupling;
            minCouplingIndex = static_cast<int>(i);
        }
    }

    // Compute maxCouplings by summing absolute couplings for each bit
    std::vector<ftype> sumCouplings(Q.getOffset().size() - 1, 0.0);
    for (int i = 0; i < Q.getOffset().size() - 1; ++i) {
        for (int k = Q.getOffset()[i]; k < Q.getOffset()[i + 1]; ++k) {
            sumCouplings[i] += affected[k].second;
        }
        if (sumCouplings[i] > maxCouplings) {
            maxCouplings = sumCouplings[i];
            maxCouplingIndex = i;
        }
    }

    // Set T_0 and T_f as inverse of betas, ensuring positive and finite values
    ftype beta_0 = -std::log(hotFlipProbability);
    ftype beta_f = -std::log(coldFlipProbability);
    *T_0 = maxCouplings / beta_0;
    *T_f = minCoupling / beta_f;

    // Ensure T_0 and T_f are positive and finite
    if (*T_0 <= 0.0 || !std::isfinite(*T_0)) {
        *T_0 = 1.0; // Fallback to a reasonable positive value
    }
    if (*T_f <= 0.0 || !std::isfinite(*T_f)) {
        *T_f = 0.1; // Fallback to a reasonable positive value
    }

    // Print debug info
    //std::cout << "Min coupling: " << minCoupling << " at affectedByFlat index " << minCouplingIndex << std::endl;
    //std::cout << "Max coupling sum: " << maxCouplings << " at variable index " << maxCouplingIndex << std::endl;
}

pair<clustering_result, clustering_result> run_vertexing(string filename) {
    bool OTF = false; // on the fly evaluation vs store qubo terms. memory and speed tradeoff
    // empirically OTF actually ends up being faster since we don't spend time constructing the qubo (!!)

    event_t event = loadTracks(filename);

    //std::cout << "Loaded " << event.nT << " tracks\n";
    //std::cout << "Loaded " << event.nV << " vertices\n";

    auto timer_start = chrono::high_resolution_clock::now();

    QUBO Q = OTF ? QUBO(event) : event_to_qubo(event);

    //std::cout << Q << std::endl;

    random_device rd;

    settings s = {
        .max_iter = SWEEPS,
        // .T_0 = 0.26*2,///10000000,
        .T_0 = 10,
        // .T_0 = 400,
        .T_f = 0.1,
        // .temp_scheduler = make_geometric_scheduler(0.999999),
        .temp_scheduler = linear_scheduler,
        .seed = rd(),
        .context = {.event=event, .max_D = get_max_D(event)},
        // .seed = 0,
        .da_sweeps = DA_SWEEPS,
    };
    
    //setTemperatureRange(hotFlipProbability, coldFlipProbability, &s.T_0, &s.T_f, Q);

    Temp_0 = s.T_0;
    Temp_f = s.T_f;

    ftype seconds;

    vector<result> results;
    result best;

    auto time_just_before_SA = chrono::high_resolution_clock::now();

    results = stage_rejoin_sa(Q, s, THREADS, STAGES, SAMPLES_PER_THREAD, filename); // adjust threads, stages, samples per thread
    //results = multithreaded_sim_anneal(Q, s, THREADS, SAMPLES_PER_THREAD, {}, filename, 0);

    auto time_just_after_SA = chrono::high_resolution_clock::now();
    
    //auto [valid_eff, correct_eff] = compute_solution_efficiencies(filename, s, flag_cluster_major);
    auto [valid_eff, correct_eff] = std::make_pair(Valid_SA_Solutions, CorrectGivenValid_SA_Solutions);

    best = results[0];

    //std::cout << "\nBranch rejoin (approach B) results: " << '\n';

    present_results(results, false);

    // s.max_iter *= 2;
    // s.temp_scheduler = linear_scheduler;
    // s.seed = rd();
    // s.dolog = false;
    // results = multithreaded_sim_anneal(Q, s, 8, 1); // threads, samples per thread

    // results = multithreaded_sim_anneal(Q, s, 8, 4); // threads, samples per thread
    // best = results[0];

    // if (results[0].energy < best.energy) {
    //     //std::cout << "choosing multithreaded result\n";
    //     best = results[0];
    // }

    // //std::cout << "\nMultithreaded (approach C) results: " << '\n';

    // present_results(results, false);

    // present_results(results);

    vector<int> assignment = interpret(best.solution, event.nT, event.nV);

    //std::cout << "Assignment: \n";

    // for (int i = 0; i < assignment.size(); i++) {
    //     //std::cout << "Track " << i << " -> Vertex " << assignment[i] << '\n';
    //     //std::cout << "track position: " << event.trackData[i].first << " vertex position: " << event.trackData[assignment[i]].first << '\n';
    // }

    // map<int, vector<int>> vertex_to_tracks;
    // for (int i = 0; i < assignment.size(); i++) {
    //     vertex_to_tracks[assignment[i]].push_back(i);
    // }

    // for (const auto& [vertex, tracks] : vertex_to_tracks) {
    //     //std::cout << "Vertex " << vertex << " tracks (" << tracks.size() << "): \n";
    //     for (int track : tracks) {
    //         //std::cout << track << " position: " << event.trackData[track].first
    //         // << " error: " << event.trackData[track].second
    //         << '\n';
    //     }
    //     //std::cout << '\n';
    // }

    ftype ari = print_score(assignment, event);

    ftype ground = ground_state(Q, event);

    //std::cout << "Ground state: " << ground << '\n';
    //std::cout << "Best energy: " << best.energy << '\n';
    //std::cout << "ratio: " << best.energy / ground << '\n';
    //std::cout << "diff: " << best.energy - ground << '\n';

    vector<ftype> vertices = assignment_to_vertices(assignment, event);
    ftype mse = vertex_mse(vertices, event);
    //std::cout << "MSE: " << mse << '\n';

    vector<int> truth;
    for (int i = 0; i < event.nT; i++) {
        truth.push_back(i / round(event.nT/event.nV));
    }

    //std::cout << "event.nT = " << event.nT << ", event.nV = " << event.nV <<endl;

    //ftype seconds = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - timer_start).count() / 1000.0;

    seconds = std::chrono::duration_cast<std::chrono::microseconds>(time_just_after_SA - time_just_before_SA).count();

    clustering_result sa_result = {filename, assignment, truth, vertices, event.vertices, best.energy, ground, valid_eff, correct_eff, seconds};

    // json_append(sa_stream, sa_result);

    // logfile_append("sa.csv", filename, ari, best.energy - ground, mse, ground);

    // return 0; // skip da.

    auto timer_mid = chrono::high_resolution_clock::now();
    //std::cout << "SA time elapsed: " << seconds << " seconds\n";

    //std::cout << "running da\n";

    auto time_just_before_DA = chrono::high_resolution_clock::now();

    pair<vector<int>, vector<ftype>> da_result = runDA(event,DA_SWEEPS);
    
    auto time_just_after_DA = chrono::high_resolution_clock::now();

    assignment = da_result.first;
    reassignUniqueIDs(assignment);
    
    // Validate assignment size and indices
    bool valid_assignment = assignment.size() == s.context.event.nT;
    if (valid_assignment) {
        for (int i = 0; i < s.context.event.nT; ++i) {
            if (assignment[i] < 0 || assignment[i] >= s.context.event.nV) {
                valid_assignment = false;
                break;
            }
        }
    }
    
    bool flag_correct_DA = false; // Default to false
    if (valid_assignment) {
        int max_val = *std::max_element(assignment.begin(), assignment.end());
        
        if (s.context.event.nV == max_val + 1) {
            vector<vector<int>> num_tracks_assigned_to_assigned_vertex_DA(s.context.event.nV);
            for (int i_internal = 0; i_internal < s.context.event.nT; i_internal++) {
                num_tracks_assigned_to_assigned_vertex_DA[assignment[i_internal]].push_back(i_internal);
            }
            
            vector<vector<int>> filtered_clusters_DA;
            for (auto& c : num_tracks_assigned_to_assigned_vertex_DA) {
                if (!c.empty()) {
                    sort(c.begin(), c.end());
                    filtered_clusters_DA.push_back(std::move(c));
                }
            }
            sort(filtered_clusters_DA.begin(), filtered_clusters_DA.end());
            
            flag_correct_DA = (filtered_clusters_DA == sorted_truth_clusters);
        }
    }
    
    Correct_DA_Solutions = flag_correct_DA;

    ari = print_score(assignment, event);

    ftype energy = energy_from_assignment(assignment, Q, event.nT, event.nV);

    //std::cout << "Ground state: " << ground << '\n';
    //std::cout << "DA energy: " << energy << '\n';
    //std::cout << "ratio: " << energy / ground << '\n';
    //std::cout << "diff: " << energy - ground << '\n';

    // vector<ftype> vertices = assignment_to_vertices(assignment, event);
    vertices = da_result.second;

    ftype da_mse = vertex_mse(vertices, event);
    //std::cout << "DA MSE: " << da_mse << '\n';

    auto timer_end = chrono::high_resolution_clock::now();

    seconds = std::chrono::duration_cast<std::chrono::microseconds>(time_just_after_DA - time_just_before_DA).count();
    
    std::tie(valid_eff, correct_eff) = std::make_pair(Correct_DA_Solutions, Correct_DA_Solutions);

    clustering_result da_clustering_result = {filename, assignment, truth, vertices, event.vertices, energy, ground, valid_eff, correct_eff, seconds};
    // json_append(da_stream, da_clustering_result);

    // logfile_append("da.csv", filename, ari, energy - ground, da_mse, ground);

    //std::cout << "SA time elapsed: " << chrono::duration_cast<chrono::seconds>(timer_mid - timer_start).count() << " seconds\n";
    //std::cout << "DA time elapsed: " << chrono::duration_cast<chrono::seconds>(timer_end - timer_mid).count() << " seconds\n";
    //std::cout << "Total time elapsed: " << chrono::duration_cast<chrono::seconds>(timer_end - timer_start).count() << " seconds\n";

    return {sa_result, da_clustering_result};
}

bool directoryExists(const std::string& path) {
    struct stat info;
    return stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode);
}

bool createDirectoriesRecursively(const std::string& path) {
    if (directoryExists(path)) {
        return true; // Already exists
    }

    size_t pos = 0;
    do {
        pos = path.find_first_of('/', pos + 1);
        std::string subdir = path.substr(0, pos == std::string::npos ? path.length() : pos);

        if (!directoryExists(subdir)) {
            if (mkdir(subdir.c_str(), 0777) != 0) {
                if (errno == EEXIST) {
                    // Directory was created by another process/thread after the check, safe to continue
                    continue;
                } else {
                    perror(("mkdir failed for " + subdir).c_str());
                    return false;
                }
            }
        }
    } while (pos != std::string::npos);

    return true;
}

ftype get_t_crit(ftype alpha, size_t df) {
    boost::math::students_t dist(df);
    return boost::math::quantile(dist,1.0 - alpha/2.0);
}

std::pair<ftype, std::pair<ftype, ftype>> calculate_se_mean_bin(const std::vector<ftype>& y) {
    size_t n = y.size();

    if (n == 0) {
        return std::make_pair(0.0, std::make_pair(0.0, 0.0));
    }

    // Handle the n=1 case specifically
    else if (n == 1) {
        // The t-critical value is undefined for 0 degrees of freedom (n-1).
        // A confidence interval cannot be formed. Return 0.0 for the t-value and se_mean.
        // Or, you could return std::numeric_limits<ftype>::infinity() for the t-value
        // to signify an infinitely wide interval.
        ftype y_mean = y[0];
        return std::make_pair(0.0, std::make_pair(y_mean, 0.0));
    }

    else {
        ftype y_mean = 0.0;
        for (const ftype& yi : y) {
            y_mean += yi;
        }
        y_mean /= n;

        ftype sum_sq_diff = 0.0;
        for (const ftype& yi : y) {
            sum_sq_diff += (yi - y_mean) * (yi - y_mean);
        }

        // n > 1 is guaranteed by the check above
        ftype s_y = std::sqrt(sum_sq_diff / (n - 1));
        ftype se_mean = s_y / std::sqrt(n);

        // This call is now safe because n > 1
        ftype t_crit = get_t_crit(0.05, n - 1);

        return std::make_pair(t_crit, std::make_pair(y_mean, se_mean));
    }
}


/*// Clopper Pearson
std::tuple<ftype, ftype, ftype> calculate_se_mean_bin(const std::vector<ftype>& y) {
    const ftype confidence_level = 0.95;

    ftype alpha_CI = 1 - confidence_level;
    const ftype z_alpha = 1.96;

    ftype lower_bound = 0.0;
    ftype upper_bound = 1.0;

    size_t n = y.size();

    if (n == 0) {
        return std::make_tuple(0.0, 0.0, 0.0);
    }

    // Handle the n=1 case specifically
    else if (n == 1) {
        // The t-critical value is undefined for 0 degrees of freedom (n-1).
        // A confidence interval cannot be formed. Return 0.0 for the t-value and se_mean.
        // Or, you could return std::numeric_limits<ftype>::infinity() for the t-value
        // to signify an infinitely wide interval.
        ftype y_mean = y[0];
        return std::make_tuple(y_mean, y_mean, y_mean);
    }

    else {
        ftype y_mean = 0.0;
        for (const ftype& yi : y) {
            y_mean += yi;
        }

        y_mean /= n;

        ftype sd = z_alpha * sqrt(y_mean*(1-y_mean)/n);

        lower_bound = y_mean - sd;

        upper_bound = y_mean + sd;

        //lower_bound = (y_mean != 0.0) ? boost::math::quantile(boost::math::beta_distribution<>(y_mean * SAMPLES_PER_THREAD * THREADS, SAMPLES_PER_THREAD * THREADS - y_mean * SAMPLES_PER_THREAD * THREADS + 1), alpha_CI / 2.0) : 0.0;

        //upper_bound = (y_mean < 1) ? boost::math::quantile(boost::math::beta_distribution<>(y_mean * SAMPLES_PER_THREAD * THREADS + 1, SAMPLES_PER_THREAD * THREADS - y_mean * SAMPLES_PER_THREAD * THREADS), 1.0 - (alpha_CI / 2.0)) : 1.0;

        return std::make_tuple(y_mean, upper_bound, lower_bound);
    }
}*/


int main(int argc, char* argv[]) {
    if (argc != 9) { // argv[0] is program name, argv[1..5] are THREADS, STAGES, SAMPLES_PER_THREAD, SWEEPS, DA_SWEEPS
        std::cerr << "Usage: ./annealer <THREADS> <STAGES> <SAMPLES_PER_THREAD> <input_filename> <output_filenames_prefix> <name_or_extension>" << std::endl;
        return 1;
    }

    // Convert command-line arguments to integers
    THREADS = std::atoi(argv[1]);
    STAGES = std::atoi(argv[2]);
    SAMPLES_PER_THREAD = std::atoi(argv[3]);
    SWEEPS = std::atoi(argv[4]);
    DA_SWEEPS = std::atoi(argv[5]);

    // Print received values (for debugging or logging)
    //std::cout << "Received: THREADS=" << THREADS << ", STAGES=" << STAGES << ", SAMPLES_PER_THREAD=" << SAMPLES_PER_THREAD << ", SWEEPS=" << SWEEPS << std::endl;

    ftype SA_TIME_AVG = 0.;

    ftype SA_Convergence_Efficiency = 0.;

    int total_CorrectGivenValid_SA_Solutions = 0;

    int total_Valid_SA_Solutions = 0;

    ftype DA_TIME_AVG = 0.;

    ftype DA_Convergence_Efficiency = 0.;

    int total_Correct_DA_Solutions = 0;

    int total_Valid_DA_Solutions = 0;

    string filename_base = argv[6];

    string extension = argv[8];

    std::string dir_name_str = std::string(argv[7]) + "/" + std::to_string(THREADS) + "threads_" + std::to_string(STAGES) + "stages_" + std::to_string(SAMPLES_PER_THREAD) + "SamplesPerThread_" + std::to_string(SWEEPS) + "sweeps_" + std::to_string(DA_SWEEPS) + "DAsweeps";
    if (createDirectoriesRecursively(dir_name_str)) {
        //std::cout << "Directories created or already existed.\n";
    } else {
        std::cerr << "Failed to create directories.\n";
    }

    std::string original1 = "ConvergenceEfficiency_and_TimePerAnneal.txt";
    //std::string renamed_original1 = std::string(argv[7])+"/"+std::string(argv[7])+"_"+original1;
    std::string renamed_original1 = dir_name_str+"/"+std::string(argv[7])+"_"+original1;

    // Check if original1 file exists
    if (file_exists(original1)) {
        //std::cout << "original1 file exists. Renaming to renamed_original1.\n";
        if (std::rename(original1.c_str(), renamed_original1.c_str()) != 0) {
            std::perror("Error renaming file");
            return 1;
        }
        else {
            //std::cout << "File moved and renamed successfully" << endl;
        }
    } else if (file_exists(renamed_original1)) {
        //std::cout << "original1 file not found. Using existing renamed_original1.\n";
    } else {
        std::cerr << "Neither original1 nor renamed_original1 file exists. Nothing to append to.\n";
        return 1;
    }

    // Open the renamed_original1 file in append mode
    std::ofstream outFile1(renamed_original1, std::ios::app);
    if (!outFile1) {
        std::cerr << "Failed to open file for appending.\n";
        return 1;
    }

    /*std::string original2 = "ConvergenceEfficiency_and_DunnIndex_Binning.txt";
    //std::string renamed_original2 = std::string(argv[7])+"/"+std::string(argv[7])+"_"+original2;
    std::string renamed_original2 = dir_name_str+"/"+std::string(argv[7])+"_"+original2;

    // Check if original2 file exists
    if (file_exists(original2)) {
        //std::cout << "original2 file exists. Renaming to renamed_original2.\n";
        if (std::rename(original2.c_str(), renamed_original2.c_str()) != 0) {
            std::perror("Error renaming file");
            return 1;
        }
        else {
            //std::cout << "File moved and renamed successfully" << endl;
        }
    } else if (file_exists(renamed_original2)) {
        //std::cout << "original2 file not found. Using existing renamed_original2.\n";
    } else {
        std::cerr << "Neither original2 nor renamed_original2 file exists. Nothing to append to.\n";
        return 1;
    }

    // Open the renamed_original2 file in append mode
    std::ofstream outFile2(renamed_original2, std::ios::app);
    if (!outFile2) {
        std::cerr << "Failed to open file for appending.\n";
        return 1;
    }

    std::vector<ftype> dunn_index_bins_edges = {0, 2, 4, 6, 10, 20, 40, 60};
    //std::vector<ftype> dunn_index_bins_edges = {0, 2, 4, 6, 10, 15};
    //std::vector<ftype> dunn_index_bins_edges = {0, 2, 4, 6, 10, 15, 25};
    std::vector<ftype> dunn_index_bins_centers;
    for (int i=0;i<dunn_index_bins_edges.size()-1;i++) {
        dunn_index_bins_centers.push_back((dunn_index_bins_edges[i+1]+dunn_index_bins_edges[i])/2);

    }
    dunn_index_bins_centers.push_back(80);
    //dunn_index_bins_centers.push_back(30);

    std::vector<ftype> dunn_sample_counts = std::vector<ftype>(dunn_index_bins_centers.size(),0);

    std::vector<ftype> t_dist_Valid_SA_val = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    std::vector<ftype> t_dist_CorrectGivenValid_SA_val = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    std::vector<ftype> t_dist_Correct_DA_val = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    
    std::vector<ftype> t_dist_SA_TIME = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    std::vector<ftype> t_dist_DA_TIME = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    
    std::vector<ftype> SAConvEff_Valid_Solutions_dunn_upper = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    std::vector<ftype> SAConvEff_Valid_Solutions_dunn_lower = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    std::vector<ftype> SAConvEff_CorrectGivenValid_Solutions_dunn_upper = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    std::vector<ftype> SAConvEff_CorrectGivenValid_Solutions_dunn_lower = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    std::vector<ftype> DAConvEff_Correct_Solutions_dunn_upper = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    std::vector<ftype> DAConvEff_Correct_Solutions_dunn_lower = std::vector<ftype>(dunn_index_bins_centers.size(),0);

    std::vector<ftype> SAConvEff_Valid_Solutions_dunn = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    std::vector<ftype> SAConvEff_Valid_Solutions_dunn_sample_sd = std::vector<ftype>(dunn_index_bins_centers.size(),0);

    std::vector<ftype> SAConvEff_CorrectGivenValid_Solutions_dunn = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    std::vector<ftype> SAConvEff_CorrectGivenValid_Solutions_dunn_sample_sd = std::vector<ftype>(dunn_index_bins_centers.size(),0);

    std::vector<std::vector<ftype>> SAConvEff_Valid_Solutions_dunn_binned(dunn_index_bins_centers.size());

    std::vector<std::vector<ftype>> SAConvEff_CorrectGivenValid_Solutions_dunn_binned(dunn_index_bins_centers.size());
    
    std::vector<ftype> SA_TIME_AVG_dunn = std::vector<ftype>(dunn_index_bins_centers.size(),0);    
    std::vector<ftype> SA_TIME_AVG_dunn_sample_sd = std::vector<ftype>(dunn_index_bins_centers.size(),0); 
    std::vector<std::vector<ftype>> SA_TIME_AVG_binned(dunn_index_bins_centers.size());

    std::vector<ftype> DAConvEff_Correct_Solutions_dunn = std::vector<ftype>(dunn_index_bins_centers.size(),0);
    std::vector<ftype> DAConvEff_Correct_Solutions_dunn_sample_sd = std::vector<ftype>(dunn_index_bins_centers.size(),0);

    std::vector<std::vector<ftype>> DAConvEff_Correct_Solutions_dunn_binned(dunn_index_bins_centers.size());
    
    std::vector<ftype> DA_TIME_AVG_dunn = std::vector<ftype>(dunn_index_bins_centers.size(),0); 
    std::vector<ftype> DA_TIME_AVG_dunn_sample_sd = std::vector<ftype>(dunn_index_bins_centers.size(),0); 
    std::vector<std::vector<ftype>> DA_TIME_AVG_binned(dunn_index_bins_centers.size());*/

    int num_files = 100;

    //std::cout << "Running " << num_files << " files\n";

    string SA_out_filename = dir_name_str + "/sa_w_omp_dunn.json";
    string DA_out_filename = dir_name_str + "/da_dunn.json";

    ofstream sa_stream(SA_out_filename);
    ofstream da_stream(DA_out_filename);

    json_init(sa_stream);
    json_init(da_stream);

    for (int i = 1; i <= num_files; i++) {
        //std::cout << "Running file " << i << '\n';

        string filename = filename_base+to_string(i)+extension;

        event_t event = loadTracks(filename);
        
        vector<uint8_t> truth;
        for (int i_int = 0; i_int < event.nT; i_int++) {
            truth.push_back(i_int / round(event.nT/event.nV));
        }

        int_truth.clear();
        for (uint8_t val1 : truth) {
            int_truth.push_back(static_cast<int>(val1));
        }

        std::map<int,std::set<int>> num_tracks_assigned_to_truth_vertex;
        for (int i_internal=0; i_internal<event.nT; i_internal++) {
            num_tracks_assigned_to_truth_vertex[int_truth[i_internal]].insert(i_internal);
        }

        valueSetsInMap_truth.clear();
        for (const auto& pairs : num_tracks_assigned_to_truth_vertex) {
            valueSetsInMap_truth.insert(pairs.second);
        }

        sorted_truth_clusters.clear();
        for (const auto& cluster : valueSetsInMap_truth) {
            vector<int> tmp(cluster.begin(), cluster.end());
            std::sort(tmp.begin(), tmp.end());
            sorted_truth_clusters.push_back(std::move(tmp));
        }
        std::sort(sorted_truth_clusters.begin(), sorted_truth_clusters.end());

        pair<clustering_result, clustering_result> p = run_vertexing(filename);
        //if (i!=1) sa_stream << ",\n"; json_append(sa_stream, p.first);
        //if (i!=1) da_stream << ",\n"; json_append(da_stream, p.second);

        //3Vertices_5Tracks is from older folder which actually contained total number of tracks instead of tracks/vertex

        SA_TIME_AVG = p.first.seconds;
        DA_TIME_AVG = p.second.seconds;

        sa_stream << ",\n"; json_append(sa_stream, p.first);
        da_stream << ",\n"; json_append(da_stream, p.second);

        total_CorrectGivenValid_SA_Solutions = p.first.correct_eff;
        total_Valid_SA_Solutions = p.first.valid_eff;

        sa_stream << "\n time for sa for each problem: " << p.first.seconds << endl;
        sa_stream << "Temperature_SA: " << Temperature_SA << endl;
        sa_stream << "Number of Correct Solutions given valid: " << total_CorrectGivenValid_SA_Solutions <<endl;
        sa_stream << "Number of valid Solutions: " << total_Valid_SA_Solutions <<endl;

        //std::cout << "Done with file " << i << '\n';

        //std::cout << "SA TIME AVG = " << SA_TIME_AVG << ", total SAs = " << SAMPLES_PER_THREAD*STAGES*num_files <<endl;

        ftype SAConvEff_Valid_Solutions = ftype(total_Valid_SA_Solutions/ftype(SAMPLES_PER_THREAD*THREADS));
        ftype SAConvEff_CorrectGivenValid_Solutions = (total_Valid_SA_Solutions != 0) ? ftype(total_CorrectGivenValid_SA_Solutions/ftype(total_Valid_SA_Solutions)) : 0;

        SA_TIME_AVG = ftype(SA_TIME_AVG/(ftype(THREADS*SAMPLES_PER_THREAD)));

        sa_stream << "\n total correct solutions: " << total_CorrectGivenValid_SA_Solutions << ", SA Convergence Efficiency correct solutions: " << SAConvEff_CorrectGivenValid_Solutions << ", SA Convergence Efficiency valid solutions: " << SAConvEff_Valid_Solutions;
        //sa_stream << "\n Average SA time: " << SAavgTimePerAnneal;

        total_Correct_DA_Solutions = p.second.correct_eff;

        da_stream << "\n time for da for each problem: " << p.second.seconds << endl;
        da_stream << "Temperature_DA: " << Temperature_DA << endl;
        da_stream << "Number of Correct Solutions given valid: " << total_Correct_DA_Solutions <<endl;

        //std::cout << "Done with file " << i << '\n';

        //std::cout << "DA TIME AVG = " << DA_TIME_AVG << ", total DAs = " << SAMPLES_PER_THREAD*STAGES*num_files <<endl;

        ftype DAConvEff_Correct_Solutions = total_Correct_DA_Solutions;

        da_stream << "\n total correct solutions: " << total_Correct_DA_Solutions << ", DA Convergence Efficiency correct solutions: " << DAConvEff_Correct_Solutions;
        //da_stream << "\n Average DA time: " << DAavgTimePerAnneal;

        /*ftype SAavgTimePerAnneal = 0.0;
        for (int j = 0; j < SAMPLES_PER_THREAD; j++) {
            SAavgTimePerAnneal += SATimePerThreadsAnneal[j];
        }
        SAavgTimePerAnneal = ftype(SAavgTimePerAnneal/(ftype(THREADS*SAMPLES_PER_THREAD*num_files)));*/

        //Remove this below if stages != 1
        /*ftype SAavgTimePerAnneal_sigma = 0.;
        for (int j = 0; j < SAMPLES_PER_THREAD; j++) {
            //std::cout << "printint times: " <<SATimePerThreadsAnneal[j]/THREADS <<", "<<SAavgTimePerAnneal<<endl;
            SAavgTimePerAnneal_sigma += pow(SATimePerThreadsAnneal[j]/THREADS - SAavgTimePerAnneal,2.);
        }
        SAavgTimePerAnneal_sigma = sqrt(ftype(SAavgTimePerAnneal_sigma*THREADS/ftype(num_files*SAMPLES_PER_THREAD*THREADS - 1)));*/
        //Remove this above if stages != 1

        json_close(sa_stream);
        json_close(da_stream);

        ftype dunn_idx = calculate_dunn_index(p.second.truth, event, event.nV, event.vertices);
        
        std::cout << dunn_idx << "\t" << Temp_0 << "\t" << Temp_f << std::endl;

        /*bool Condition;

        for (int j=0 ; j<dunn_index_bins_edges.size() ; j++) {
            Condition = false;
            if (j==dunn_index_bins_edges.size()-1) {
                if (dunn_idx >= dunn_index_bins_edges[j]) {
                    Condition = true;
                }
            }
            else {
                if (dunn_idx >= dunn_index_bins_edges[j] && dunn_idx < dunn_index_bins_edges[j+1]) {
                    Condition = true;
                }
            }

            if (Condition) {
                SAConvEff_Valid_Solutions_dunn_binned[j].push_back(SAConvEff_Valid_Solutions);

                SAConvEff_CorrectGivenValid_Solutions_dunn_binned[j].push_back(SAConvEff_CorrectGivenValid_Solutions);

                DAConvEff_Correct_Solutions_dunn_binned[j].push_back(DAConvEff_Correct_Solutions);
                
                SA_TIME_AVG_binned[j].push_back(SA_TIME_AVG);
                
                DA_TIME_AVG_binned[j].push_back(DA_TIME_AVG);
            }
        }*/

        // Append data
        outFile1 << THREADS <<"\t"<< STAGES << "\t" << SWEEPS << "\t" << DA_SWEEPS << "\t" << dunn_idx << "\t" << SAConvEff_Valid_Solutions << "\t" << SAConvEff_CorrectGivenValid_Solutions << "\t" << SA_TIME_AVG  << "\t" << DAConvEff_Correct_Solutions << "\t" << DA_TIME_AVG << endl;
    }
    outFile1.close();

    /*for (int j = 0; j < dunn_index_bins_centers.size(); j++) {
        // Call once for SA Valid Solutions
        auto sa_valid_result = calculate_se_mean_bin(SAConvEff_Valid_Solutions_dunn_binned[j]);
        t_dist_Valid_SA_val[j] = sa_valid_result.first;
        SAConvEff_Valid_Solutions_dunn[j] = sa_valid_result.second.first;
        SAConvEff_Valid_Solutions_dunn_sample_sd[j] = sa_valid_result.second.second;

        // Call once for SA Correct Given Valid Solutions
        auto sa_correct_result = calculate_se_mean_bin(SAConvEff_CorrectGivenValid_Solutions_dunn_binned[j]);
        t_dist_CorrectGivenValid_SA_val[j] = sa_correct_result.first;
        SAConvEff_CorrectGivenValid_Solutions_dunn[j] = sa_correct_result.second.first;
        SAConvEff_CorrectGivenValid_Solutions_dunn_sample_sd[j] = sa_correct_result.second.second;
        
        // Call once for DA Correct Solutions
        auto da_correct_result = calculate_se_mean_bin(DAConvEff_Correct_Solutions_dunn_binned[j]);
        t_dist_Correct_DA_val[j] = da_correct_result.first;
        DAConvEff_Correct_Solutions_dunn[j] = da_correct_result.second.first;
        DAConvEff_Correct_Solutions_dunn_sample_sd[j] = da_correct_result.second.second;
        
        auto sa_time_avg_binned = calculate_se_mean_bin(SA_TIME_AVG_binned[j]);
        t_dist_SA_TIME[j] = sa_time_avg_binned.first;
        SA_TIME_AVG_dunn[j] = sa_time_avg_binned.second.first;
        SA_TIME_AVG_dunn_sample_sd[j] = sa_time_avg_binned.second.second;
        
        auto da_time_avg_binned = calculate_se_mean_bin(DA_TIME_AVG_binned[j]);
        t_dist_DA_TIME[j] = da_time_avg_binned.first;
        DA_TIME_AVG_dunn[j] = da_time_avg_binned.second.first;
        DA_TIME_AVG_dunn_sample_sd[j] = da_time_avg_binned.second.second;
        
        SAConvEff_Valid_Solutions_dunn_upper[j] = SAConvEff_Valid_Solutions_dunn[j] + SAConvEff_Valid_Solutions_dunn_sample_sd[j];
        SAConvEff_Valid_Solutions_dunn_lower[j] = SAConvEff_Valid_Solutions_dunn[j] - SAConvEff_Valid_Solutions_dunn_sample_sd[j];
        SAConvEff_CorrectGivenValid_Solutions_dunn_upper[j] = SAConvEff_CorrectGivenValid_Solutions_dunn[j] + SAConvEff_CorrectGivenValid_Solutions_dunn_sample_sd[j];
        SAConvEff_CorrectGivenValid_Solutions_dunn_lower[j] = SAConvEff_CorrectGivenValid_Solutions_dunn[j] - SAConvEff_CorrectGivenValid_Solutions_dunn_sample_sd[j];
        DAConvEff_Correct_Solutions_dunn_upper[j] = DAConvEff_Correct_Solutions_dunn[j] + DAConvEff_Correct_Solutions_dunn_sample_sd[j];
        DAConvEff_Correct_Solutions_dunn_lower[j] = DAConvEff_Correct_Solutions_dunn[j] - DAConvEff_Correct_Solutions_dunn_sample_sd[j];
        
        //std::cout << "calculate_se_mean_bin " << SAConvEff_Valid_Solutions_dunn_binned[j][0] << "\t" << SAConvEff_Valid_Solutions_dunn[j] << endl;
    }

    for (int j=0;j<dunn_index_bins_centers.size();j++) {
        if (j<dunn_index_bins_centers.size()-1) {
            outFile2 << THREADS <<"\t"<< STAGES << "\t" << SWEEPS << "\t" << DA_SWEEPS << "\t" << t_dist_Valid_SA_val[j] << "\t" << dunn_index_bins_centers[j] << "\t" << dunn_index_bins_edges[j] << "\t" << dunn_index_bins_edges[j+1] << "\t" << SAConvEff_Valid_Solutions_dunn[j] << "\t" << SAConvEff_Valid_Solutions_dunn_sample_sd[j] << "\t" << SAConvEff_CorrectGivenValid_Solutions_dunn[j] << "\t" << SAConvEff_CorrectGivenValid_Solutions_dunn_sample_sd[j] << "\t" << SA_TIME_AVG_dunn[j] << "\t" << SA_TIME_AVG_dunn_sample_sd[j] << "\t" << DAConvEff_Correct_Solutions_dunn[j] << "\t" << DAConvEff_Correct_Solutions_dunn_sample_sd[j] << "\t" << DA_TIME_AVG_dunn[j] << "\t" << DA_TIME_AVG_dunn_sample_sd[j] << endl;
        }
        else {
            outFile2 << THREADS <<"\t"<< STAGES << "\t" << SWEEPS << "\t" << DA_SWEEPS << "\t" << t_dist_Valid_SA_val[j] << "\t" << dunn_index_bins_centers[j] << "\t" << dunn_index_bins_edges[j] << "\t" << 100 << "\t" << SAConvEff_Valid_Solutions_dunn[j] << "\t" << SAConvEff_Valid_Solutions_dunn_sample_sd[j] << "\t" << SAConvEff_CorrectGivenValid_Solutions_dunn[j] << "\t" << SAConvEff_CorrectGivenValid_Solutions_dunn_sample_sd[j] << "\t" << SA_TIME_AVG_dunn[j] << "\t" << SA_TIME_AVG_dunn_sample_sd[j] << "\t" << DAConvEff_Correct_Solutions_dunn[j] << "\t" << DAConvEff_Correct_Solutions_dunn_sample_sd[j] << "\t" << DA_TIME_AVG_dunn[j] << "\t" << DA_TIME_AVG_dunn_sample_sd[j] << endl;
        }
    }*/

    /*for (int j = 0; j < dunn_index_bins_centers.size(); j++) {
        // Call once for SA Valid Solutions
        tie(SAConvEff_Valid_Solutions_dunn[j], SAConvEff_Valid_Solutions_dunn_upper[j], SAConvEff_Valid_Solutions_dunn_lower[j]) = calculate_se_mean_bin(SAConvEff_Valid_Solutions_dunn_binned[j]);

        // Call once for SA Correct Given Valid Solutions
        tie(SAConvEff_CorrectGivenValid_Solutions_dunn[j], SAConvEff_CorrectGivenValid_Solutions_dunn_upper[j], SAConvEff_CorrectGivenValid_Solutions_dunn_lower[j]) = calculate_se_mean_bin(SAConvEff_CorrectGivenValid_Solutions_dunn_binned[j]);

        // Call once for DA Correct Solutions
        tie(DAConvEff_Correct_Solutions_dunn[j], DAConvEff_Correct_Solutions_dunn_upper[j], DAConvEff_Correct_Solutions_dunn_lower[j]) = calculate_se_mean_bin(DAConvEff_Correct_Solutions_dunn_binned[j]);
    }


    for (int j=0;j<dunn_index_bins_centers.size();j++) {
        if (j<dunn_index_bins_centers.size()-1) {
            outFile2 << THREADS <<"\t"<< STAGES << "\t" << SWEEPS << "\t" << t_dist_Valid_SA_val[j] << "\t" << dunn_index_bins_centers[j] << "\t" << dunn_index_bins_edges[j] << "\t" << dunn_index_bins_edges[j+1] << "\t" << SAConvEff_Valid_Solutions_dunn[j] << "\t" << SAConvEff_Valid_Solutions_dunn_upper[j] << "\t" << SAConvEff_Valid_Solutions_dunn_lower[j] << "\t" << SAConvEff_CorrectGivenValid_Solutions_dunn[j] << "\t" << SAConvEff_CorrectGivenValid_Solutions_dunn_upper[j] << "\t" << SAConvEff_CorrectGivenValid_Solutions_dunn_lower[j] << "\t" << DAConvEff_Correct_Solutions_dunn[j] << "\t" << DAConvEff_Correct_Solutions_dunn_upper[j] << "\t" << DAConvEff_Correct_Solutions_dunn_lower[j] << endl;
        }
        else {
            outFile2 << THREADS <<"\t"<< STAGES << "\t" << SWEEPS << "\t" << t_dist_Valid_SA_val[j] << "\t" << dunn_index_bins_centers[j] << "\t" << dunn_index_bins_edges[j] << "\t" << 100 << "\t" << SAConvEff_Valid_Solutions_dunn[j] << "\t" << SAConvEff_Valid_Solutions_dunn_upper[j] << "\t" << SAConvEff_Valid_Solutions_dunn_lower[j] << "\t" << SAConvEff_CorrectGivenValid_Solutions_dunn[j] << "\t" << SAConvEff_CorrectGivenValid_Solutions_dunn_upper[j] << "\t" << SAConvEff_CorrectGivenValid_Solutions_dunn_lower[j] << "\t" << DAConvEff_Correct_Solutions_dunn[j] << "\t" << DAConvEff_Correct_Solutions_dunn_upper[j] << "\t" << DAConvEff_Correct_Solutions_dunn_lower[j] << endl;
        }
    }*/

    //outFile2.close();

    return 0;
}
