#pragma once

#include <vector>
#include <functional>
#include <vector>
#include <map>

using namespace std;

using ftype = float;
using solution_t = vector<uint8_t>;
using qubo_t = map<pair<int, int>, ftype>;
using scheduler_t = function<ftype(ftype T_0, ftype T, int iter, int max_iter)>;

struct problem_context {
    int nT;
    int nV;
};

class QUBO {
    private:
        int n;

        vector<pair<int, ftype>> affectedby_flat; // list of js that are affected by flipping a certain bit
        vector<int> offset; // offset of each row in the flat list
    public:
        QUBO(qubo_t& Q);
        ftype evaluate(const solution_t& x) const;
        ftype evaluateDiff(const solution_t& x, int flip_idx) const;
};