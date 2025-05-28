#pragma once
#include <vector>
#include <functional>
#include <vector>
#include <map>

using namespace std;

using ftype = double; // use double for actual math stuff
using solution_t = vector<uint8_t>;
using qubo_t = map<pair<int, int>, ftype>;
using scheduler_t = function<ftype(ftype T_0, ftype T_f, ftype T, int iter, int max_iter)>;

struct event_t;

class QUBO {
    private:
        int n;

        bool OTF;
        event_t* event; // optional pointer to event
        ftype max_D;

        // use float for params to save memory
        vector<pair<int, float>> affectedby_flat; // list of js that are affected by flipping a certain bit
        vector<int> offset; // offset of each row in the flat list
    public:
        QUBO(event_t& event); // OTF mode.
        QUBO(qubo_t& Q);
        ftype evaluate(const solution_t& x) const;
        ftype evaluateDiff(const solution_t& x, int flip_idx) const;
};
