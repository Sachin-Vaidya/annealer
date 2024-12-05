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

vector<int> interpret(const solution_t &solution, const int nT, const int nV);

int run_vertexing(int argc, char* argv[]);
event_t loadTracks(string filename);
qubo_t event_to_qubo(const event_t &event);