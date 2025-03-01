#pragma once

#include <vector>
#include "annealer.hh"

using namespace std;

struct event_t {
    int nV;
    int nT;
    const vector<pair<ftype, ftype>> trackData; // z position, error
    const vector<ftype> vertices; // actual vertex z position
};

struct problem_context {
    const event_t &event;
    ftype max_D;
};

// todo: types and qubo struct should be in a separate header file

vector<int> interpret(const solution_t &solution, const int nT, const int nV);

int run_vertexing(int argc, char* argv[]);
event_t loadTracks(string filename);
QUBO event_to_qubo(const event_t &event);
ftype adjustedRandIndex(const vector<int> &a, const vector<int> &b);
ftype vertex_mse(const vector<ftype> &vertices, const event_t &event);
vector<ftype> assignment_to_vertices(const vector<int> &assignment, const event_t &event);
ftype print_score(const vector<int> &assignment, const event_t &event);

ftype ground_state(const QUBO &qubo, const event_t &event);
ftype energy_from_assignment(const vector<int> &assignment, const QUBO &qubo, const int nT, const int nV);

ftype get_max_D(const event_t &event);
ftype evaluate_diff_on_the_fly(const solution_t &x, const event_t &event, int flip_idx, ftype max_D);