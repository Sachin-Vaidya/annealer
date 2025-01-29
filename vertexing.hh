#pragma once

#include <vector>
#include <functional>
#include <vector>
#include <map>
#include <iostream>

using namespace std;

// todo: types and qubo struct should be in a separate header file

using ftype = float;
using solution_t = vector<uint8_t>;
using qubo_t = map<pair<int, int>, ftype>;
using scheduler_t = function<ftype(ftype T_0, ftype T, int iter, int max_iter)>;

struct event_t {
    int nV;
    int nT;
    vector<pair<ftype, ftype>> trackData; // z position, error
};

struct problem_context {
    int nT;
    int nV;
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

    // todo: more memory efficient data structure?
    // vector<vector<pair<int, ftype>>> affectedby; // list of js that are affected by flipping a certain bit

    vector<pair<int, ftype>> affectedby_flat; // list of js that are affected by flipping a certain bit
    vector<int> offset; // offset of each row in the flat list

    QUBO(qubo_t& Q) {
        n = qubo_size(Q);
        // affectedby.resize(n, {});
        // for (const auto& [idx, val] : Q) {
        //     affectedby[idx.first].emplace_back(idx.second, val);
        //     if (idx.first == idx.second) continue;
        //     affectedby[idx.second].emplace_back(idx.first, val);
        // }

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

    // void affectedby_stats_print() {
    //     cout << "affectedby stats: \n";
    //     cout << "n = " << n << " \n";
    //     cout << "affectedby size = " << affectedby.size() << " \n";
    //     for (int i = 0; i < n; i++) {
    //         int count = affectedby[i].size();
    //         cout << count << " ";
    //     }
    // }

    ftype evaluate(const solution_t& x) const {
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

    ftype evaluateDiff(const solution_t& x, int flip_idx) const {
        ftype diff = 0.0; // first, find what would be the value if this bit was on
// #pragma clang loop vectorize_width(2)
// #pragma clang loop interleave_count(2)
        // for (const auto& [j, Q_ij] : affectedby[flip_idx])
        //     if (x[j] || j == flip_idx) diff += Q_ij;


        for (int k = offset[flip_idx]; k < offset[flip_idx + 1]; k++) {
            const auto& [j, Q_ij] = affectedby_flat[k];
            if (x[j] || j == flip_idx) diff += Q_ij;
            // diff += (x[j] || j == flip_idx) * Q_ij;
        }
        return x[flip_idx] ? -diff : diff; // if on, turn off. if off, turn on.
    }
};

vector<int> interpret(const solution_t &solution, const int nT, const int nV);

int run_vertexing(int argc, char* argv[]);
event_t loadTracks(string filename);
QUBO event_to_qubo(const event_t &event);
ftype adjustedRandIndex(const vector<int> &a, const vector<int> &b);
void print_score(const vector<int> &assignment, const event_t &event);

ftype ground_state(const QUBO &qubo, const event_t &event);
ftype energy_from_assignment(const vector<int> &assignment, const QUBO &qubo, const int nT, const int nV);