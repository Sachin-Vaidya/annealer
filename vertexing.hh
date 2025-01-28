#pragma once

#include <vector>
// #include <algorithm>
#include <functional>
#include <vector>
#include <map>

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
    vector<vector<pair<int, ftype>>> affectedby; // list of js that are affected by flipping a certain bit

    QUBO(qubo_t& Q) {
        n = qubo_size(Q);
        affectedby.resize(n, {});
        for (const auto& [idx, val] : Q) {
            affectedby[idx.first].emplace_back(idx.second, val);
            if (idx.first == idx.second) continue;
            affectedby[idx.second].emplace_back(idx.first, val);
        }
    }

    ftype evaluate(const solution_t& x) const {
        ftype value = 0.0;
        for (int i = 0; i < x.size(); i++)
            if (x[i])
                for (const auto& [j, Q_ij] : affectedby[i])
                    if (x[j] && j >= i) value += Q_ij; // only count each pair once
        return value;
    }

    ftype evaluateDiff(const solution_t& x, int flip_idx) const {
        ftype diff = 0.0; // first, find what would be the value if this bit was on
// #pragma clang loop vectorize_width(2)
// #pragma clang loop interleave_count(2)
        for (const auto& [j, Q_ij] : affectedby[flip_idx])
            if (x[j] || j == flip_idx) diff += Q_ij;
            // diff += (x[j] || j == flip_idx) * Q_ij;
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