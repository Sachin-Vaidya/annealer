#include "vertexing.hh"

ostream& operator << (ostream& os, const solution_t& x) {
    for (auto xi : x) os << xi << ' ';
    return os;
}

ostream& operator << (ostream& os, const qubo_t& Q) {
    for (const auto& entry : Q) {
        os << '[' << entry.first.first << ' ' << entry.first.second << "] : " << entry.second << endl;
    }
    return os;
}

int qubo_size(const qubo_t& Q) {
    int n = 0;
    for (const auto& entry : Q) {
        n = max(n, entry.first.first);
    }
    return n + 1; // 0 indexed
}

struct QUBO {
    int n;
    qubo_t Q;
    // todo: more memory efficient data structure?
    vector<vector<pair<int, double>>> affectedby; // list of js that are affected by flipping a certain bit

    QUBO(qubo_t Q) : Q(Q) {
        n = qubo_size(Q);
        affectedby.resize(n, {});
        // for (const auto& entry : Q) {
        for (const auto& [idx, val] : Q) {
            affectedby[idx.first].emplace_back(idx.second, val);
            if (idx.first == idx.second) continue;
            affectedby[idx.second].emplace_back(idx.first, val);
        }
    }

    double evaluate(const solution_t& x) const {
        double value = 0.0;
        for (const auto& [idx, val] : Q)
            if (x[idx.first] && x[idx.second]) value += val;
        return value;
    }

    double evaluateDiff(const solution_t& x, int flip_idx) const {
        double diff = 0.0; // first, find what would be the value if this bit was on
        for (const auto& [j, Q_ij] : affectedby[flip_idx]) if (x[j] || j == flip_idx) diff += Q_ij;
        return x[flip_idx] ? -diff : diff; // if on, turn off. if off, turn on.
    }

    friend ostream& operator << (ostream& os, const QUBO& Q) {
        os << "QUBO of size " << Q.n << ":\n";
        for (const auto& entry : Q.Q) {
            os << '[' << entry.first.first << ' ' << entry.first.second << "] : " << entry.second << endl;
        }
        return os;
    }
};

struct settings {
    int max_iter;
    double T_0;
    scheduler_t temp_scheduler;
    unsigned seed;
};

struct result {
    solution_t solution;
    double energy;
};

result sim_anneal(const QUBO& Q, const settings s, const solution_t init_guess = {}) { // intentionally get copy of settings
    // int n = qubo_size(Q);

    mt19937 gen(s.seed);
    uniform_real_distribution<> dis(0.0, 1.0);
    uniform_int_distribution<> flip_dis(0, Q.n - 1);

    solution_t x(Q.n, 0);
    for (auto xi : x) {
        xi = dis(gen) < 0.5 ? 0 : 1;
    }

    if (!init_guess.empty()) x = init_guess;

    double f_x = Q.evaluate(x);

    solution_t best_x = x;
    double best_f_x = f_x;

    double T = s.T_0;
    for (int iter = 0; iter < s.max_iter; iter++) {
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

        T = s.temp_scheduler(s.T_0, T, iter, s.max_iter);
    }

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

// todo: we shouldn't need to do this...
// hack for now.
void upper_to_lower(qubo_t& Q) {
    assert_upper_triangular(Q);
    for (auto& entry : Q) {
        if (entry.first.first < entry.first.second) {
            swap(entry.first.first, entry.first.second);
        }
    }
}

// returns sorted results of length num_threads * samples_per_thread. first elem is best.
vector<result> multithreaded_sim_anneal(const QUBO& Q, const settings s, int num_threads, int samples_per_thread = 1, const vector<solution_t> init_guess = {}) {
    assert_lower_triangular(Q.Q);
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

// takes in a sorted list of results and returns at least n solutions with as many unique solutions as possible that all have the best energy
vector<solution_t> best_effort_unique(const vector<result>& results, int n) {
    if (results.empty()) return {};

    vector<solution_t> result(n);
    map<solution_t, bool> best_unique_solutions;

    for (int i = 0; i < results.size(); i++) {
        if (results[i].energy == results[0].energy) {
            best_unique_solutions[results[i].solution] = true;
        } else {
            break;
        }
    }

    int num_best_unique = best_unique_solutions.size();

    int i = 0;
    for (const auto& [sol, _] : best_unique_solutions) {
        if (i == n) break;
        result[i++] = sol;
    }
    for (int j = 0; i < n; i++, j++) { // copy the start to fill the rest
        result[i] = result[j];
    }

    // debug

    cout << "Num best unique: " << num_best_unique << '\n';
    // for (const auto& sol : result) {
    //     cout << sol << '\n';
    // }
    // cout << '\n';

    return result;
}


vector<result> branch_rejoin_sa(const QUBO& Q, const settings s, int num_threads, int num_branches, int samples_per_thread = 1) {
    auto modified_settings = s;
    modified_settings.max_iter /= num_branches;
    vector<result> results;
    for (int i = 0; i < num_branches; i++) {
        results = multithreaded_sim_anneal(Q, modified_settings, num_threads, samples_per_thread, best_effort_unique(results, num_threads));
        cout << "Branch " << i << " best energy: " << results[0].energy << '\n';
        cout << "Branch " << i << " worst energy: " << results.back().energy << '\n';
    }
    return results;
}

double linear_scheduler(double T_0, double T, int iter, int max_iter) {
    return T_0 - (T_0 / max_iter) * iter;
}

scheduler_t make_geometric_scheduler(double alpha) {
    return [alpha](double T_0, double T, int iter, int max_iter) {
        return T * alpha;
    };
}


void trial(solution_t x, const QUBO& Q) {
    cout << "For solution: " << x << endl;
    cout << "Energy: " << Q.evaluate(x) << endl << endl;
}

void present_results(const vector<result>& results, bool show_sols = true, int precision = 5) {
    map<long, map<solution_t, int>> counts; // energy -> solution -> count
    auto d_to_l = [precision](double d) {
        return static_cast<long>(round(d * pow(10, precision)));
    };
    auto l_to_d = [precision](long l) {
        return static_cast<double>(l) / pow(10, precision);
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

qubo_t condense(vector<vector<double>> sparse) {
    qubo_t dense;

    if (sparse.size() != sparse[0].size()) {
        cerr << "Error: matrix is not square." << endl;
        exit(1);
    }

    for (int i = 0; i < sparse.size(); i++) {
        for (int j = 0; j <= i; j++) { // only add lower triangular
            if (sparse[i][j] != 0) { // only add non-zero entries
                // dense[{i, j}] = sparse[i][j];
                dense.push_back({{i, j}, sparse[i][j]});
            }
        }
    }

    return dense;
}

vector<vector<double>> sparsen(const qubo_t& dense) {
    int n = qubo_size(dense);
    vector<vector<double>> sparse(n, vector<double>(n, 0.0));
    for (const auto& entry : dense) {
        sparse[entry.first.first][entry.first.second] = entry.second;
    }
    return sparse;
}

QUBO randgen_qubo(int n) {
    qubo_t Q;
    random_device rd;
    unsigned seed = rd();
    mt19937 gen(seed);
    uniform_real_distribution<> dis(-10.0, 10.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            // Q[{i, j}] = dis(gen);
            Q.push_back({{i, j}, dis(gen)});
        }
    }
    return Q;
}

int run_vertexing(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage: ./annealer <filename>\n";
        return 1;
    }

    string filename = argv[1];
    event_t event = loadTracks(filename);

    cout << "Loaded " << event.nT << " tracks\n";
    cout << "Loaded " << event.nV << " vertices\n";
    // for (int i = 0; i < trackData.trackData.size(); i++) {
    //     cout << trackData.trackData[i].first << " " << trackData.trackData[i].second << endl;
    // }

    qubo_t qubo = event_to_qubo(event);

    cout << "QUBO size: " << qubo_size(qubo) << '\n';

    upper_to_lower(qubo);

    cout << "converted to lower triangular\n";

    auto Q = QUBO(qubo);

    // cout << Q;

    random_device rd;
    // unsigned seed = rd();

    settings s = {.max_iter = 40000, .T_0 = 100.0,
    .temp_scheduler = make_geometric_scheduler(0.999),
    // .temp_scheduler = linear_scheduler,
    .seed = rd()};

    vector<result> results = branch_rejoin_sa(Q, s, 8, 4, 4); // threads, branches, samples per thread

    result best = results[0];

    cout << "\nBranch rejoin (approach B) results: " << '\n';

    present_results(results, false);

    // present_results(results);

    // s.max_iter *= 2;

    // s.temp_scheduler = linear_scheduler;

    s.seed = rd();

    results = multithreaded_sim_anneal(Q, s, 8, 4); // threads, samples per thread

    // best = results[0];

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
        cout << "Vertex " << vertex << " tracks: \n";
        for (int track : tracks) {
            cout << track << " position: " << event.trackData[track].first 
            // << " error: " << event.trackData[track].second
            << '\n';
        }
        cout << '\n';
    }

    return 0;
}

int main(int argc, char* argv[]) {
    run_vertexing(argc, argv);
}

// int main() {
//     // cout << fixed << setprecision(numeric_limits<double>::max_digits10);
//     // cout << fixed << setprecision(5);
//     // qubo_t Q = condense({
//     //     { -2,  1,  1,  0,  0 },
//     //     {  1, -2,  0,  1,  0 },
//     //     {  1,  0, -3,  1,  1 },
//     //     {  0,  1,  1, -3,  1 },
//     //     {  0,  0,  1,  1, -2 }
//     // });

//     // qubo_t Q = condense({
//     //     {-17, 10, 10, 10, 0, 20},
//     //     {10, -18, 10, 10, 10, 20},
//     //     {10, 10, -29, 10, 20, 20},
//     //     {10, 10, 10, -19, 10, 10},
//     //     {0, 10, 20, 10, -17, 10},
//     //     {20, 20, 20, 10, 10, -28}
//     // });

//     auto t = "qubo.txt";

//     // t = "2v.json";

//     auto Q = QUBO(parse_qubo(read_file(t)));

//     // auto Q = randgen_qubo(6000);

//     // cout << Q;

//     // auto Q = randgen_qubo(10);

//     random_device rd;
//     // unsigned seed = rd();

//     settings s = {.max_iter = 40000, .T_0 = 100.0,
//     .temp_scheduler = make_geometric_scheduler(0.999), 
//     // .temp_scheduler = linear_scheduler,
//     .seed = rd()};

//     vector<result> results = branch_rejoin_sa(Q, s, 8, 4, 4); // threads, branches, samples per thread

//     result best = results[0];

//     cout << "\nBranch rejoin (approach B) results: " << '\n';
//     present_results(results, false);
//     // present_results(results);

//     // s.max_iter *= 2;
//     // s.temp_scheduler = linear_scheduler;

//     s.seed = rd();

//     results = multithreaded_sim_anneal(Q, s, 8, 4); // threads, samples per thread

//     best = results[0];

//     cout << "\nMultithreaded (approach C) results: " << '\n';
//     present_results(results, false);
//     // present_results(results);

//     vector<int> assignment = interpret(best.solution, Q.n/2, 2); // nT, nV

//     cout << "Assignment: ";
//     for (int i = 0; i < assignment.size(); i++) {
//         cout << "Track " << i << " -> Vertex " << assignment[i] << '\n';
//     }

//     return 0;
// }
