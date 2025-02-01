#pragma once

#include <vector>
#include "annealer.hh"

using namespace std;

struct event_t {
    int nV;
    int nT;
    vector<pair<ftype, ftype>> trackData; // z position, error
};

// todo: types and qubo struct should be in a separate header file

vector<int> interpret(const solution_t &solution, const int nT, const int nV);

int run_vertexing(int argc, char* argv[]);
event_t loadTracks(string filename);
QUBO event_to_qubo(const event_t &event);
ftype adjustedRandIndex(const vector<int> &a, const vector<int> &b);
void print_score(const vector<int> &assignment, const event_t &event);

ftype ground_state(const QUBO &qubo, const event_t &event);
ftype energy_from_assignment(const vector<int> &assignment, const QUBO &qubo, const int nT, const int nV);