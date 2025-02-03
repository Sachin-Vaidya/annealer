#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <thread>

#include "annealer.hh"
#include "vertexing.hh"
#include "detanneal.hh"

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
    scheduler_t temp_scheduler;
    unsigned seed;
    problem_context context;
    bool dolog = true;
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

QUBO::QUBO(qubo_t& Q) {
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

    cout << "conversion done\n";
    cout << "size of flat list = " << affectedby_flat.size() << "\n";
}

ftype QUBO::evaluate(const solution_t& x) const {
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
    ftype diff = 0.0; // first, find what would be the value if this bit was on
// #pragma clang loop vectorize_width(2)
// #pragma clang loop interleave_count(2)
    for (int k = offset[flip_idx]; k < offset[flip_idx + 1]; k++) {
        const auto& [j, Q_ij] = affectedby_flat[k];
        if (x[j] || j == flip_idx) diff += Q_ij;
        // diff += (x[j] || j == flip_idx) * Q_ij;
    }
    return x[flip_idx] ? -diff : diff; // if on, turn off. if off, turn on.
}

// problem specific anneal
result sim_anneal(const QUBO& Q, const settings s, const solution_t init_guess = {}) { // intentionally get copy of settings
    // int n = qubo_size(Q);

    mt19937 gen(s.seed);
    uniform_real_distribution<> dis(0.0, 1.0);

    uniform_int_distribution<> t_dis(0, s.context.nT - 1);
    uniform_int_distribution<> v_dis(0, s.context.nV - 1);

    solution_t x(s.context.nT * s.context.nV, 0);

    auto bit_idx = [s](int t, int v) {
        return t + s.context.nT * v;
    };

    // helper.
    vector<int> track_to_vertex(s.context.nT, -1);

    // for each track, assign random vertex
    for (int t = 0; t < s.context.nT; t++) {
        int v = v_dis(gen);
        track_to_vertex[t] = v;
        x[bit_idx(t, v)] = 1;
    }

    if (!init_guess.empty()) {
        x = init_guess;
        cout << "Using init guess\n";
    }

    ftype f_x = Q.evaluate(x);

    solution_t best_x = x;
    ftype best_f_x = f_x;

    ftype T = s.T_0;
    for (int iter = 0; iter < s.max_iter; iter++) {
        if (iter % 10000 == 0 && s.dolog) {
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

        T = s.temp_scheduler(s.T_0, T, iter, s.max_iter);
    }

    best_f_x = Q.evaluate(best_x); // re-evaluate best solution to remove float point errors

    return {best_x, best_f_x};
}

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
vector<result> multithreaded_sim_anneal(const QUBO& Q, const settings s, int num_threads, int samples_per_thread = 1, const vector<solution_t> init_guess = {}) {
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

    for (int i = 0; i < num_threads; i++) {
        threads.emplace_back([&results, i, &Q, s, samples_per_thread, local_init_guess](){
            for (int j = 0; j < samples_per_thread; j++) {
                settings s_copy = s;
                s_copy.seed += i * samples_per_thread + j; // different seed for each thread
                results[i * samples_per_thread + j] = sim_anneal(Q, s_copy, local_init_guess[i]);
            }
        });
    }

    for (auto& t : threads) t.join();

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


vector<result> branch_rejoin_sa(const QUBO& Q, const settings s, int num_threads, int num_branches, int samples_per_thread = 1) {
    settings modified_settings = s;
    modified_settings.max_iter /= num_branches;
    vector<result> results;
    for (int i = 0; i < num_branches; i++) {
        results = multithreaded_sim_anneal(Q, modified_settings, num_threads, samples_per_thread, best_effort_unique(results, num_threads));
        cout << "Branch " << i << " best energy: " << results[0].energy << '\n';
        cout << "Branch " << i << " worst energy: " << results.back().energy << '\n';
    }
    return results;
}

ftype linear_scheduler(ftype T_0, ftype T, int iter, int max_iter) {
    return T_0 - (T_0 / max_iter) * iter;
}

scheduler_t make_geometric_scheduler(ftype alpha) {
    return [alpha](ftype T_0, ftype T, int iter, int max_iter) {
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

    cout << '\n';

    cout << fixed << setprecision(precision);

    cout << "Best energy: " << results[0].energy << '\n';
    cout << "Worst energy: " << results.back().energy << '\n';
    // cout << "Best solution: " << results[0].solution << '\n';

    cout << '\n';

    for (const auto& [energy, sols] : counts) {
        int total = 0;
        for (const auto& [sol, count] : sols) {
            total += count;
        }
        cout << "Energy: " << l_to_d(energy) << " (" << total << "x)" << '\n';
        if (show_sols) {
            for (const auto& [sol, count] : sols) {
                cout << "\t Sol: " << sol << " (" << count << "x)" << '\n';
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

int run_vertexing(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage: ./annealer <filename>\n";
        return 1;
    }

    string filename = argv[1];
    event_t event = loadTracks(filename);

    cout << "Loaded " << event.nT << " tracks\n";
    cout << "Loaded " << event.nV << " vertices\n";

    QUBO Q = event_to_qubo(event);

    // cout << Q;

    random_device rd;

    settings s = {.max_iter = 800000,
    // .T_0 = 0.26*2,///10000000,
    .T_0 = 0,
    // .T_0 = 400,
    // .temp_scheduler = make_geometric_scheduler(0.999999),
    .context = {.nT = event.nT, .nV = event.nV},
    .temp_scheduler = linear_scheduler,
    .seed = rd()
    // .seed = 0,
    };

    vector<result> results;
    result best;

    // results = branch_rejoin_sa(Q, s, 8, 4, 1); // threads, branches, samples per thread

    // best = results[0];

    // cout << "\nBranch rejoin (approach B) results: " << '\n';

    // present_results(results, false);

    // s.max_iter *= 2;
    // s.temp_scheduler = linear_scheduler;
    // s.seed = rd();
    // s.dolog = false;
    // results = multithreaded_sim_anneal(Q, s, 8, 1); // threads, samples per thread
    
    results = multithreaded_sim_anneal(Q, s, 8, 1); // threads, samples per thread
    best = results[0];

    if (results[0].energy < best.energy) {
        cout << "choosing multithreaded result\n";
        best = results[0];
    }

    cout << "\nMultithreaded (approach C) results: " << '\n';

    present_results(results, false);

    // present_results(results);

    vector<int> assignment = interpret(best.solution, event.nT, event.nV);

    cout << "Assignment: \n";

    // for (int i = 0; i < assignment.size(); i++) {
    //     cout << "Track " << i << " -> Vertex " << assignment[i] << '\n';
    //     cout << "track position: " << event.trackData[i].first << " vertex position: " << event.trackData[assignment[i]].first << '\n';
    // }

    map<int, vector<int>> vertex_to_tracks;
    for (int i = 0; i < assignment.size(); i++) {
        vertex_to_tracks[assignment[i]].push_back(i);
    }

    for (const auto& [vertex, tracks] : vertex_to_tracks) {
        cout << "Vertex " << vertex << " tracks (" << tracks.size() << "): \n";
        for (int track : tracks) {
            cout << track << " position: " << event.trackData[track].first 
            // << " error: " << event.trackData[track].second
            << '\n';
        }
        cout << '\n';
    }

    print_score(assignment, event);

    ftype ground = ground_state(Q, event);

    cout << "Ground state: " << ground << '\n';
    cout << "Best energy: " << best.energy << '\n';
    cout << "ratio: " << best.energy / ground << '\n';
    cout << "diff: " << best.energy - ground << '\n';

    cout << "running da\n";

    vector<int> da_assignment = runDA(event);

    print_score(da_assignment, event);

    ftype da_energy = energy_from_assignment(da_assignment, Q, event.nT, event.nV);

    cout << "Ground state: " << ground << '\n';
    cout << "DA energy: " << da_energy << '\n';
    cout << "ratio: " << da_energy / ground << '\n';
    cout << "diff: " << da_energy - ground << '\n';

    return 0;
}

int main(int argc, char* argv[]) {
    run_vertexing(argc, argv);
}
