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

#include "annealer.hh"
#include "vertexing.hh"
#include "detanneal.hh"

int THREADS;
int STAGES;
int SAMPLES_PER_THREAD;
int SWEEPS;

double Temperature;
double SA_timer_start;
double SA_timer_stop;

std::vector<std::vector<int>> CorrectSolutionOrNot;
int Correct_Solutions;
std::vector<ftype> SATimePerThreadsAnneal;

// Function to check if a file exists
bool file_exists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

ostream& operator << (ostream& os, const solution_t& x) {
    for (auto xi : x) os << xi << ' ';
    return os;
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

    //cout << "conversion done\n";
    //cout << "size of flat list = " << affectedby_flat.size() << "\n";
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
        if (x[j] || j == flip_idx) diff += Q_ij;
        // diff += (x[j] || j == flip_idx) * Q_ij;
    }
    return diff * (1 - 2 * x[flip_idx]);
    //return x[flip_idx] ? -diff : diff; // if on, turn off. if off, turn on.
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

        //cout << "Using init guess\n";
    }

    ftype f_x = Q.evaluate(x);

    solution_t best_x = x;
    ftype best_f_x = f_x;

    ftype T = s.T_0;
    for (int iter = 0; iter < s.max_iter; iter++) {
        if (iter % 1000 == 0 && s.dolog) {
            //cout << "Iter: " << iter << " Energy: " << f_x << " T: " << T << '\n';
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
        Temperature = T;
        
        cout << "\nTemperature calculation verification: " << num_stage << ", " << T << ", " << Temperature << endl;
    }

    // re-evaluate best solution to eliminate float point errors
    best_f_x = Q.evaluate(best_x);

    return {best_x, best_f_x};
}


//Sweep through entire bit string and flip them

result sim_anneal(const QUBO& Q, const settings s, const solution_t init_guess = {}, int num_stage = 1) {
    mt19937 gen(s.seed);
    uniform_real_distribution<> dis(0.0, 1.0);

    int nT = s.context.event.nT;
    int nV = s.context.event.nV;

    auto bit_idx = [nT](int t, int v) {
        return t + nT * v;
    };

    solution_t x(nT * nV, 0);
    vector<int> track_to_vertex(nT, -1);

    // Random initial guess if not provided
    for (int t = 0; t < nT; t++) {
        int v = gen() % nV;
        track_to_vertex[t] = v;
        x[bit_idx(t, v)] = 1;
    }

    if (!init_guess.empty()) {
        x = init_guess;
        for (int t = 0; t < nT; t++) {
            for (int v = 0; v < nV; v++) {
                if (x[bit_idx(t, v)]) {
                    track_to_vertex[t] = v;
                }
            }
        }
    }

    ftype f_x = Q.evaluate(x);
    solution_t best_x = x;
    ftype best_f_x = f_x;

    ftype T = s.T_0;

    for (int iter = 0; iter < s.max_iter; iter++) {
        if (iter % 1000 == 0 && s.dolog) {
            // cout << "Iter: " << iter << " Energy: " << f_x << " T: " << T << '\n';
        }

        // SWEEP through all bits
        for (int t = 0; t < nT; ++t) {
            for (int v = 0; v < nV; ++v) {
                int bit = bit_idx(t, v);
                
                solution_t x_prime = x;
                ftype delta = Q.evaluateDiff(x_prime, bit);
                x_prime[bit] = 1 - x[bit]; // toggle bit
                ftype f_x_prime = f_x + delta;

                // Metropolis criterion
                if (f_x_prime < f_x || (T != 0 && dis(gen) < exp((f_x - f_x_prime) / T))) {
                    x = x_prime;
                    f_x = f_x_prime;

                    // Optionally update track_to_vertex, but unclear validity for multiple 1s
                    //track_to_vertex[t] = -1;
                }

                if (f_x < best_f_x) {
                    best_x = x;
                    best_f_x = f_x;
                }
            }
        }

        // Lower temperature after full sweep
        T = s.temp_scheduler(s.T_0, s.T_f, T, iter, s.max_iter, num_stage);
        Temperature = T;

        cout << "\nTemperature calculation verification: " << num_stage << ", " << T << ", " << Temperature << endl;
    }

    best_f_x = Q.evaluate(best_x); // cleanup
    return {best_x, best_f_x};
}

//=================================== ====================================== =========================================//
*/

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

    double f_x = Q.evaluate(x);

    solution_t best_x = x;
    double best_f_x = f_x;

    double T = s.T_0;
    for (int iter = 0; iter < s.max_iter; iter++) {
        if (iter % 10000 == 0 && s.dolog) {
            cout << "Iter: " << iter << " Energy: " << best_f_x << " T: " << T << '\n';
        }

        solution_t x_prime = x;
        int flip_idx = flip_dis(gen);

        double f_x_prime = Q.evaluateDiff(x, flip_idx) + f_x;

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
        Temperature = T;

        cout << "\nTemperature calculation verification: " << num_stage << ", " << T << ", " << Temperature << endl;
    }

    best_f_x = Q.evaluate(best_x); // cleanup
    return {best_x, best_f_x};
}
*/

//SA with sweep

result sim_anneal(const QUBO& Q, const settings s, const solution_t init_guess = {}, int num_stage = 1) {
    mt19937 gen(s.seed);
    uniform_real_distribution<> dis(0.0, 1.0);

    solution_t x(s.context.event.nT * s.context.event.nV, 0);
    for (auto &xi : x) {
        xi = dis(gen) < 0.5 ? 0 : 1;
    }

    if (!init_guess.empty()) x = init_guess;

    double f_x = Q.evaluate(x);
    solution_t best_x = x;
    double best_f_x = f_x;

    double T = s.T_0;

    for (int iter = 0; iter < s.max_iter; iter++) {
        if (iter % 1000 == 0 && s.dolog) {
            cout << "Iter: " << iter << " Energy: " << best_f_x << " T: " << T << '\n';
        }

        for (int i = 0; i < static_cast<int>(x.size()); i++) {
            double f_x_prime = Q.evaluateDiff(x, i) + f_x;

            if (f_x_prime < f_x || dis(gen) < exp((f_x - f_x_prime) / T)) {
                x[i] = !x[i];  // Flip bit
                f_x = f_x_prime;
            }

            if (f_x < best_f_x) {
                best_x = x;
                best_f_x = f_x;
            }
        }

        T = s.temp_scheduler(s.T_0, s.T_f, T, iter, s.max_iter, num_stage);
        Temperature = T;

        cout << "\nTemperature calculation verification: " << num_stage << ", " << T << ", " << Temperature << endl;
    }

    best_f_x = Q.evaluate(best_x); // cleanup
    return {best_x, best_f_x};
}

//=================================== ====================================== =========================================//

/*
//==================================================== SA* ============================================================//
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

        cout << "Using init guess\n";
    }

    ftype f_x = Q.evaluate(x);

    solution_t best_x = x;
    ftype best_f_x = f_x;

    ftype T = s.T_0;
    for (int iter = 0; iter < s.max_iter; iter++) {
        if (iter % 1000 == 0 && s.dolog) {
            cout << "Iter: " << iter << " Energy: " << f_x << " T: " << T << '\n';
        }

        int t = t_dis(gen); // pick random track to change vertex of
        
        int old_v = track_to_vertex[t];
        int new_v = v_dis(gen); // pick random new vertex for track

        // if (old_v == new_v) continue;

        int old_bit = bit_idx(t, old_v);
        int new_bit = bit_idx(t, new_v);

        solution_t x_prime = x;

        ftype delta = Q.evaluateDiff(x_prime, old_bit);

        x_prime[old_bit] = 0;

        delta += Q.evaluateDiff(x_prime, new_bit);

        x_prime[new_bit] = 1;

        ftype f_x_prime = f_x + delta;

        if (f_x_prime < f_x || (T != 0 && dis(gen) < exp((f_x - f_x_prime) / T))) {
            x = x_prime;
            f_x = f_x_prime;
            track_to_vertex[t] = new_v;
        }

        if (f_x < best_f_x) {
            best_x = x;
            best_f_x = f_x;
        }

        T = s.temp_scheduler(s.T_0, s.T_f, T, iter, s.max_iter, num_stage);
        Temperature = T;
        
        cout << "\nTemperature calculation verification: " << num_stage << ", " << T << ", " << Temperature << endl;
    }

    // re-evaluate best solution to eliminate float point errors
    best_f_x = Q.evaluate(best_x);

    return {best_x, best_f_x};
}
//==================================================== === ============================================================//
*/

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

// returns sorted results of length num_threads * samples_per_thread. first elem is best.
vector<result> multithreaded_sim_anneal(const QUBO& Q, const settings s, int num_threads, int samples_per_thread = 1, const vector<solution_t> init_guess = {}, string filename = "", int num_stage=1) {
    vector<thread> threads;
    vector<result> results(num_threads * samples_per_thread);
    vector<solution_t> local_init_guess = init_guess;

    if (init_guess.size() != num_threads) {
        if (!init_guess.empty()) {
            local_init_guess = vector<solution_t>(num_threads, init_guess[0]);
        } else {
            local_init_guess = vector<solution_t>(num_threads, solution_t{});
        }
    }

    /*
    // threads: manually
    for (int i = 0; i < num_threads; i++) {
        threads.emplace_back([&results, i, &Q, s, samples_per_thread, local_init_guess](){
            for (int j = 0; j < samples_per_thread; j++) {
                settings s_copy = s;
                s_copy.seed += i * samples_per_thread + j; // different seed for each thread
                results[i * samples_per_thread + j] = sim_anneal(Q, s_copy, local_init_guess[i]);
            }
        });
    }
    */
    
    // Using OpenMP
    omp_set_num_threads(num_threads);
    
    event_t event = loadTracks(filename);
    ftype ground = ground_state(Q, event);
    
    Correct_Solutions = 0;
    
    SATimePerThreadsAnneal = std::vector<ftype>(samples_per_thread, 0.0);
    CorrectSolutionOrNot = std::vector<std::vector<int>>(num_threads, std::vector<int>(samples_per_thread, 0));
    
    for (int j = 0; j < samples_per_thread; j++) 
    {
        SATimePerThreadsAnneal[j] = omp_get_wtime();
        #pragma omp parallel
        {
            int tid = omp_get_thread_num(); // actual thread ID
            settings s_copy = s;
            s_copy.seed += tid * samples_per_thread + j; // different seed for each thread
            
            #ifdef __linux__
            int cpu = sched_getcpu();
            printf("Thread %d handling, sample %d running on CPU %d\n", tid, j, cpu);
            #endif

            
            
            
            results[tid * samples_per_thread + j] = sim_anneal(Q, s_copy, local_init_guess[tid], num_stage);
            
            
            
            
            
            //what's given below
            //results[i * samples_per_thread + j] = sim_anneal(Q, s_copy, local_init_guess[i], num_stage);
            
            
            
            
            if (abs(ground - results[tid * samples_per_thread + j].energy) < 1.e-9) {
                CorrectSolutionOrNot[tid][j] = 1;
            }
        }
        SATimePerThreadsAnneal[j] = omp_get_wtime() - SATimePerThreadsAnneal[j];
    }
    
    for (int i = 0; i < num_threads; i++) {
        for (int j = 0; j < samples_per_thread; j++) {
            Correct_Solutions += CorrectSolutionOrNot[i][j];
        }
    }
    
    /*Correct_Solutions = 0;
    for (int i = 0; i < num_threads; i++) {
        for (int j = 0; j < samples_per_thread; j++) {
            cout << "Correct Solutions: " << i <<", " << j << ", " << ground << ", " << results[i * samples_per_thread + j].energy <<endl;
            if (abs(ground - results[i * samples_per_thread + j].energy) < 1.e-9) {
                Correct_Solutions += 1;    
            }
        }
    }*/
    
    //for (auto& t : threads) t.join();

    sort(results.begin(), results.end(), [](const result& a, const result& b) {
        return a.energy < b.energy;
    });

    return results;
}

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

//----------num_branches is actually num_stages-----------//
vector<result> branch_rejoin_sa(const QUBO& Q, const settings s, int num_threads, int num_branches, int samples_per_thread = 1, string filename = "") {
    settings modified_settings = s;
    modified_settings.max_iter /= num_branches;
    vector<result> results;
    for (int i = 0; i < num_branches; i++) {
        results = multithreaded_sim_anneal(Q, modified_settings, num_threads, samples_per_thread, best_effort_unique(results, num_threads), filename, i);
        //cout << "Branch " << i << " best energy: " << results[0].energy << '\n';
        //cout << "Branch " << i << " worst energy: " << results.back().energy << '\n';
        cout << "Temperature at the end of stage " << i << ": " << Temperature <<endl;
    }
    cout << "Max Iter for stage: " << modified_settings.max_iter <<endl;
    return results;
}

ftype linear_scheduler(ftype T_0, ftype T_f, ftype T, int iter, int max_iter, int num_stage) {
    //return T_0 - double(T_0 / double((max_iter - 1)*STAGES)) * (num_stage*(max_iter - 1)+iter);
    return pow(1 / T_0 - (1 / T_0 - 1 / T_f) * (num_stage * (max_iter - 1) + iter) / (max_iter - 1) / STAGES,-1.);
    //cout << "checking variable values: " << T_0 << "\t" << max_iter << "\t" << num_stage << "\t" << iter << endl;
}

scheduler_t make_geometric_scheduler(ftype alpha) {
    return [alpha](ftype T_0, ftype T_f, ftype T, int iter, int max_iter, int num_stage) {
        return T * alpha;
    };
}


// void trial(solution_t x, const QUBO& Q) {
//     cout << "For solution: " << x << endl;
//     cout << "Energy: " << Q.evaluate(x) << endl << endl;
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

    //cout << '\n';

    //cout << fixed << setprecision(precision);

    //cout << "Best energy: " << results[0].energy << '\n';
    //cout << "Worst energy: " << results.back().energy << '\n';
    // cout << "Best solution: " << results[0].solution << '\n';

    //cout << '\n';

    for (const auto& [energy, sols] : counts) {
        int total = 0;
        for (const auto& [sol, count] : sols) {
            total += count;
        }
        //cout << "Energy: " << l_to_d(energy) << " (" << total << "x)" << '\n';
        if (show_sols) {
            for (const auto& [sol, count] : sols) {
                //cout << "\t Sol: " << sol << " (" << count << "x)" << '\n';
            }
        }
    }
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
    ftype seconds;
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

pair<clustering_result, clustering_result> run_vertexing(string filename) {
    bool OTF = false; // on the fly evaluation vs store qubo terms. memory and speed tradeoff
    // empirically OTF actually ends up being faster since we don't spend time constructing the qubo (!!)

    event_t event = loadTracks(filename);

    //cout << "Loaded " << event.nT << " tracks\n";
    //cout << "Loaded " << event.nV << " vertices\n";

    auto timer_start = chrono::high_resolution_clock::now();

    QUBO Q = OTF ? QUBO(event) : event_to_qubo(event);

    // cout << Q;

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
    };

    vector<result> results;
    result best;
    
    auto time_just_before_SA = chrono::high_resolution_clock::now();
    
    //results = branch_rejoin_sa(Q, s, thread::hardware_concurrency(), 4, 1); // threads, branches, samples per thread
    results = branch_rejoin_sa(Q, s, THREADS, STAGES, SAMPLES_PER_THREAD, filename); // adjust threads, branches, samples per thread
    
    auto time_just_after_SA = chrono::high_resolution_clock::now();
    
    best = results[0];

    //cout << "\nBranch rejoin (approach B) results: " << '\n';

    present_results(results, false);

    // s.max_iter *= 2;
    // s.temp_scheduler = linear_scheduler;
    // s.seed = rd();
    // s.dolog = false;
    // results = multithreaded_sim_anneal(Q, s, 8, 1); // threads, samples per thread



    // results = multithreaded_sim_anneal(Q, s, 8, 4); // threads, samples per thread
    // best = results[0];

    // if (results[0].energy < best.energy) {
    //     cout << "choosing multithreaded result\n";
    //     best = results[0];
    // }

    // cout << "\nMultithreaded (approach C) results: " << '\n';

    // present_results(results, false);



    // present_results(results);

    vector<int> assignment = interpret(best.solution, event.nT, event.nV);

    //cout << "Assignment: \n";

    // for (int i = 0; i < assignment.size(); i++) {
    //     cout << "Track " << i << " -> Vertex " << assignment[i] << '\n';
    //     cout << "track position: " << event.trackData[i].first << " vertex position: " << event.trackData[assignment[i]].first << '\n';
    // }

    // map<int, vector<int>> vertex_to_tracks;
    // for (int i = 0; i < assignment.size(); i++) {
    //     vertex_to_tracks[assignment[i]].push_back(i);
    // }

    // for (const auto& [vertex, tracks] : vertex_to_tracks) {
    //     cout << "Vertex " << vertex << " tracks (" << tracks.size() << "): \n";
    //     for (int track : tracks) {
    //         cout << track << " position: " << event.trackData[track].first 
    //         // << " error: " << event.trackData[track].second
    //         << '\n';
    //     }
    //     cout << '\n';
    // }

    ftype ari = print_score(assignment, event);

    ftype ground = ground_state(Q, event);

    //cout << "Ground state: " << ground << '\n';
    //cout << "Best energy: " << best.energy << '\n';
    //cout << "ratio: " << best.energy / ground << '\n';
    //cout << "diff: " << best.energy - ground << '\n';


    vector<ftype> vertices = assignment_to_vertices(assignment, event);
    ftype mse = vertex_mse(vertices, event);
    //cout << "MSE: " << mse << '\n';

    vector<int> truth;
    for (int i = 0; i < event.nT; i++) {
        truth.push_back(i / round(event.nT/event.nV));
    }

    //ftype seconds = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - timer_start).count() / 1000.0;
    ftype seconds = chrono::duration_cast<chrono::milliseconds>(time_just_after_SA - time_just_before_SA).count() / 1000.0;
    
    clustering_result sa_result = {filename, assignment, truth, vertices, event.vertices, best.energy, ground, seconds};
    // json_append(sa_stream, sa_result);

    // logfile_append("sa.csv", filename, ari, best.energy - ground, mse, ground);

    // return 0; // skip da.

    auto timer_mid = chrono::high_resolution_clock::now();
    //cout << "SA time elapsed: " << seconds << " seconds\n";

    cout << "running da\n";

    pair<vector<int>, vector<ftype>> da_result = runDA(event);

    vector<int> da_assignment = da_result.first;

    ari = print_score(da_assignment, event);

    ftype da_energy = energy_from_assignment(da_assignment, Q, event.nT, event.nV);

    //cout << "Ground state: " << ground << '\n';
    //cout << "DA energy: " << da_energy << '\n';
    //cout << "ratio: " << da_energy / ground << '\n';
    //cout << "diff: " << da_energy - ground << '\n';

    // vector<ftype> da_vertices = assignment_to_vertices(da_assignment, event);
    vector<ftype> da_vertices = da_result.second;
    
    ftype da_mse = vertex_mse(da_vertices, event);
    //cout << "DA MSE: " << da_mse << '\n';

    auto timer_end = chrono::high_resolution_clock::now();

    seconds = chrono::duration_cast<chrono::milliseconds>(timer_end - timer_mid).count() / 1000.0;

    clustering_result da_clustering_result = {filename, da_assignment, truth, da_vertices, event.vertices, da_energy, ground, seconds};
    // json_append(da_stream, da_clustering_result);

    // logfile_append("da.csv", filename, ari, da_energy - ground, da_mse, ground);

    // cout << "SA time elapsed: " << chrono::duration_cast<chrono::seconds>(timer_mid - timer_start).count() << " seconds\n";
    // cout << "DA time elapsed: " << chrono::duration_cast<chrono::seconds>(timer_end - timer_mid).count() << " seconds\n";
    //cout << "Total time elapsed: " << chrono::duration_cast<chrono::seconds>(timer_end - timer_start).count() << " seconds\n";

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

int main(int argc, char* argv[]) {
    if (argc != 8) { // argv[0] is program name, argv[1..4] are THREADS, STAGES, SAMPLES_PER_THREAD, SWEEPS
        std::cerr << "Usage: ./annealer <THREADS> <STAGES> <SAMPLES_PER_THREAD> <input_filename> <output_filenames_prefix> <name_or_extension>" << std::endl;
        return 1;
    }

    // Convert command-line arguments to integers
    THREADS = std::atoi(argv[1]);
    STAGES = std::atoi(argv[2]);
    SAMPLES_PER_THREAD = std::atoi(argv[3]);
    SWEEPS = std::atoi(argv[4]);

    // Print received values (for debugging or logging)
    std::cout << "Received: THREADS=" << THREADS << ", STAGES=" << STAGES 
              << ", SAMPLES_PER_THREAD=" << SAMPLES_PER_THREAD << ", SWEEPS=" << SWEEPS << std::endl;
    
    double SA_TIME_AVG = 0.;
    
    double SA_Convergence_Efficiency = 0.;
    
    int total_Correct_Solutions = 0;
    
    string filename_base = argv[5];
    
    string extension = argv[7];
    
    std::string dir_name_str = std::string(argv[6]) + "/" + std::to_string(THREADS) + "threads_" + std::to_string(STAGES) + "stages_" + std::to_string(SAMPLES_PER_THREAD) + "SamplesPerThread_" + std::to_string(SWEEPS) + "sweeps";
    if (createDirectoriesRecursively(dir_name_str)) {
        std::cout << "Directories created or already existed.\n";
    } else {
        std::cerr << "Failed to create directories.\n";
    }
    
    int num_files = 1;

    cout << "Running " << num_files << " files\n";
    
    string SA_out_filename = dir_name_str + "/sa_w_omp_min_dunn.json";
    string DA_out_filename = dir_name_str + "/da_min_dunn.json";
    
    ofstream sa_stream(SA_out_filename);
    ofstream da_stream(DA_out_filename);

    json_init(sa_stream);
    json_init(da_stream);

    for (int i = 1; i <= num_files; i++) {
        cout << "Running file " << i << '\n';

        string filename = filename_base+to_string(i)+extension;
        
        pair<clustering_result, clustering_result> p = run_vertexing(filename);
        //if (i!=1) sa_stream << ",\n"; json_append(sa_stream, p.first);
        //if (i!=1) da_stream << ",\n"; json_append(da_stream, p.second);
        
        //3Vertices_5Tracks is from older folder which actually contained total number of tracks instead of tracks/vertex
        
        SA_TIME_AVG += p.first.seconds;
        
        sa_stream << ",\n"; json_append(sa_stream, p.first);
        da_stream << ",\n"; json_append(da_stream, p.second);
        
        total_Correct_Solutions += Correct_Solutions;
        
        sa_stream << "\n time for sa for each problem: " << p.first.seconds << endl;
        sa_stream << "Temperature: " << Temperature << endl;
        sa_stream << "Number of Correct Solutions: " << Correct_Solutions <<endl;
        
        cout << "Done with file " << i << '\n';
    }
    
    //cout << "SA TIME AVG = " << SA_TIME_AVG << ", total SAs = " << SAMPLES_PER_THREAD*STAGES*num_files <<endl;
    
    ftype SAConvEff = double(total_Correct_Solutions/double(num_files*SAMPLES_PER_THREAD*THREADS));
    
    ftype SAavgTimePerAnneal = 0.0;
    for (int j = 0; j < SAMPLES_PER_THREAD; j++) {
        SAavgTimePerAnneal += SATimePerThreadsAnneal[j];
    }
    SAavgTimePerAnneal = double(SAavgTimePerAnneal/(double(THREADS*SAMPLES_PER_THREAD*num_files)));
    
    SA_TIME_AVG = double(SA_TIME_AVG/(double(THREADS*SAMPLES_PER_THREAD*num_files)));
    
    //Remove this below if stages != 1
    ftype SAavgTimePerAnneal_sigma = 0.;
    for (int j = 0; j < SAMPLES_PER_THREAD; j++) {
        cout << "printint times: " <<SATimePerThreadsAnneal[j]/THREADS <<", "<<SAavgTimePerAnneal<<endl;
        SAavgTimePerAnneal_sigma += pow(SATimePerThreadsAnneal[j]/THREADS - SAavgTimePerAnneal,2.);
    }
    SAavgTimePerAnneal_sigma = sqrt(double(SAavgTimePerAnneal_sigma*THREADS/double(num_files*SAMPLES_PER_THREAD*THREADS - 1)));
    //Remove this above if stages != 1
    
    
    
    
    
    
    sa_stream << "\n total correct solutions: " << total_Correct_Solutions << ", SA Convergence Efficiency: " << SAConvEff;
    sa_stream << "\n Average SA time: " << SAavgTimePerAnneal;
    
    json_close(sa_stream);
    json_close(da_stream);
    
    std::string original = "SA_ConvergenceEfficiency_and_TimePerAnneal.txt";
    std::string renamed_original = std::string(argv[6])+"/"+std::string(argv[6])+"SA_ConvergenceEfficiency_and_TimePerAnneal.txt";

    // Check if original file exists
    if (file_exists(original)) {
        std::cout << "Original file exists. Renaming to renamed_original.\n";
        if (std::rename(original.c_str(), renamed_original.c_str()) != 0) {
            std::perror("Error renaming file");
            return 1;
        }
        else {
            std::cout << "File moved and renamed successfully" << endl;
        }
    } else if (file_exists(renamed_original)) {
        std::cout << "Original file not found. Using existing renamed_original.\n";
    } else {
        std::cerr << "Neither original nor renamed_original file exists. Nothing to append to.\n";
        return 1;
    }

    // Open the renamed_original file in append mode
    std::ofstream outFile(renamed_original, std::ios::app);
    if (!outFile) {
        std::cerr << "Failed to open file for appending.\n";
        return 1;
    }

    // Append data
    outFile << THREADS <<"\t"<< STAGES << "\t" << SWEEPS << "\t" << SAConvEff << "\t" << SA_TIME_AVG << "\t" << SAavgTimePerAnneal << "\t" << SAavgTimePerAnneal_sigma << endl;
    outFile.close();
    
    return 0;
}
