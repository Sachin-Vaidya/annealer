#include <iostream>
#include <ostream>
#include <vector>
#include <random>
#include <thread>
#include <algorithm>
#include <functional>
#include <map>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include "ThreadPool.h"

using namespace std;

using solution_t = vector<bool>;
// using qubo_t = map<pair<int, int>, double>;
using qubo_t = vector<pair<pair<int, int>, double>>;
using scheduler_t = function<double(double T_0, double T, int iter, int max_iter)>;
using event_t = struct {
    int nV;
    int nT;
    vector<pair<double, double>> trackData; // z position, error
};

inline int qubo_size(const qubo_t& Q) {
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

    ThreadPool& threadPool;

    // QUBO(qubo_t Q) : Q(Q) {
    //     n = qubo_size(Q);
    //     affectedby.resize(n, {});
    //     // for (const auto& entry : Q) {
    //     for (const auto& [idx, val] : Q) {
    //         affectedby[idx.first].emplace_back(idx.second, val);
    //         if (idx.first == idx.second) continue;
    //         affectedby[idx.second].emplace_back(idx.first, val);
    //     }
    // }
    QUBO(qubo_t Q, ThreadPool& tp) : Q(Q), threadPool(tp) {
        n = qubo_size(Q);
        affectedby.resize(n, {});
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


    // double evaluateDiff(const solution_t& x, int flip_idx) const {
    //     const auto& affected = affectedby[flip_idx];
    //     // TODO: make sure same as num threads in pool.
    //     int num_threads = POOL_SIZE; 
    //     int chunk_size = (int) ((affected.size() + num_threads - 1) / num_threads);

    //     if (chunk_size == 0) {
    //         // just do sequential
    //         double diff = 0.0;
    //         for (auto & kv : affected) {
    //             const int j = kv.first;
    //             const double Q_ij = kv.second;
    //             if (x[j] || j == flip_idx) diff += Q_ij;
    //         }
    //         return x[flip_idx] ? -diff : diff;
    //     }

    //     vector<future<double>> futures;
    //     futures.reserve(num_threads);

    //     for (int t = 0; t < num_threads; t++) {
    //         int start = t * chunk_size;
    //         int end = min((t + 1) * chunk_size, (int) affected.size());
    //         if (start >= (int) affected.size()) break;

    //         futures.push_back(threadPool.enqueue([&affected, &x, flip_idx, start, end]() {
    //             double local_sum = 0.0;
    //             for (int i = start; i < end; i++) {
    //                 const auto& [j, Q_ij] = affected[i];
    //                 if (x[j] || j == flip_idx) local_sum += Q_ij;
    //             }
    //             return local_sum;
    //         }));
    //     }

    //     double total_diff = 0.0;
    //     for (auto &f : futures) {
    //         total_diff += f.get();
    //     }

    //     return x[flip_idx] ? -total_diff : total_diff;
    // }


    friend ostream& operator << (ostream& os, const QUBO& Q) {
        os << "QUBO of size " << Q.n << ":\n";
        for (const auto& entry : Q.Q) {
            os << '[' << entry.first.first << ' ' << entry.first.second << "] : " << entry.second << endl;
        }
        return os;
    }
};

vector<int> interpret(const solution_t &solution, const int nT, const int nV);

int run_vertexing(int argc, char* argv[]);
event_t loadTracks(string filename);
qubo_t event_to_qubo(const event_t &event);
double adjustedRandIndex(const vector<int> &a, const vector<int> &b);
void print_score(const vector<int> &assignment, const event_t &event);

double ground_state(const QUBO &qubo, const event_t &event);