#include "vertexing.hh"

vector<int> interpret(const solution_t &solution, const int nT, const int nV) {
    // returns a vector of integers, where the i-th element is the vertex that
    // the i-th track is assigned to.

    // solution must be of size nT * nV
    assert(solution.size() == nT * nV);

    vector<int> assignment(nT, -1);

    for (int i = 0; i < solution.size(); i++) {
        if (solution[i]) {
            if (assignment[i % nT] != -1) { // track already assigned
                cout << "track " << i % nT << " already assigned!\n";
                // return vector<int>(); // invalid solution
            }
            int track = i % nT;
            int vertex = i / nT;
            assignment[track] = vertex;
        }
    }

    if (find(assignment.begin(), assignment.end(), -1) != assignment.end()) {
        cout << "not all tracks assigned!\n";
        // return vector<int>(); // invalid solution
        return assignment;
    }

    return assignment;
}

ftype energy_from_assignment(const vector<int> &assignment, const QUBO &qubo, const int nT, const int nV) {
    solution_t solution(qubo.n, 0);

    for (int i = 0; i < assignment.size(); i++) {
        int vertex = assignment[i];
        solution[i + nT * vertex] = 1;
    }

    return qubo.evaluate(solution);
}

void print_score(const vector<int> &assignment, const event_t &event) {
    // print the ARI score of the assignment
    vector<int> true_labels(event.nT);
    for (int i = 0; i < event.nT; i++) {
        true_labels[i] = i/30; // 30 tracks per vertex. TODO: magic number bad
    }

    // debug print true labels and assignment
    cout << "true labels: ";
    for (int i = 0; i < true_labels.size(); i++) {
        cout << true_labels[i] << " ";
    }
    cout << endl;

    cout << "assignment: ";
    for (int i = 0; i < assignment.size(); i++) {
        cout << assignment[i] << " ";
    }
    cout << endl;

    cout << "ARI: " << adjustedRandIndex(true_labels, assignment) << endl;
}

ftype ground_state(const QUBO &qubo, const event_t &event) {
    solution_t solution(qubo.n, 0);

    auto idx = [event](int track, int vertex) {
        return track + event.nT * vertex;
    };

    for (int i = 0; i < event.nT; i++) {
        int vertex = i / 30; // 30 tracks per vertex. TODO: magic number bad
        solution[idx(i, vertex)] = 1;
    }

    return qubo.evaluate(solution);
}

// https://arxiv.org/pdf/1903.08879
qubo_t event_to_qubo(const event_t &event) {
    int nV = event.nV;
    int nT = event.nT;
    const auto& trackData = event.trackData;

    const ftype scale = 1.5;

    map<pair<int, int>, ftype> qubo_map;

    auto D = [](pair<ftype, ftype> i, pair<ftype, ftype> j) {
        return abs(i.first - j.first) / sqrt(i.second * i.second + j.second * j.second);
    };

    auto g = [scale](ftype x, ftype m = 5) {
        // return 1.0 - exp(-m * x);
        // return x + log(1.0 + x);
        // return x;

        return scale * (1.0 - exp(-m * x));
        // cout << x << '\n';
        // return 1.0 - exp(-x);
        // cout << 1.0 - exp(-x) << '\n';
        // return 1.0 - exp(-x);
    };

    auto idx = [nT](int track, int vertex) {
        return track + nT * vertex;
    };

    // ftype lambda = 1.0;
    ftype lambda = 2;
    // ftype lambda = 1.2;

    ftype max_D = 0.0;

    // for (int k = 0; k < nV; ++k) {
    //     for (int i = 0; i < nT; ++i) {
    //         for (int j = i + 1; j < nT; ++j) {
    //             ftype D_ij = D(trackData[i], trackData[j]);
    //             max_D = max(max_D, D_ij);
    //             qubo_map[{idx(i, k), idx(j, k)}] += g(D_ij);
    //         }
    //     }
    // }

    for (int k = 0; k < nV; ++k) {
        for (int i = 0; i < nT; ++i) {
            for (int j = i + 1; j < nT; ++j) {
                ftype D_ij = D(trackData[i], trackData[j]);
                max_D = max(max_D, D_ij);
            }
        }
    }

    // int num_nonzero = 0;

    for (int k = 0; k < nV; ++k) {
        for (int i = 0; i < nT; ++i) {
            for (int j = i + 1; j < nT; ++j) {
                ftype D_ij = D(trackData[i], trackData[j]);
                // qubo_map[{idx(i, k), idx(j, k)}] += g(D_ij);
                // qubo_map[{idx(i, k), idx(j, k)}] += g(D_ij/max_D);
                qubo_map[{idx(i, k), idx(j, k)}] += g(D_ij/max_D);
            }
        }
    }

    // lambda *= max_D;

    // penalty
    for (int i = 0; i < nT; ++i) {
        for (int k = 0; k < nV; ++k) {
            qubo_map[{idx(i, k), idx(i, k)}] -= lambda;
            for (int k2 = k + 1; k2 < nV; ++k2) {
                qubo_map[{idx(i, k), idx(i, k2)}] += 2 * lambda;
            }
        }
    }

    qubo_t qubo(qubo_map.begin(), qubo_map.end());

    cout << "qubo num terms: " << qubo.size() << '\n';
    cout << "max possible: " << nT * nV * nT * nV << '\n';

    return qubo;
}


// source: andrew wildridge. lightly modified
event_t loadTracks(string filename) {
    ifstream trackFile(filename.c_str());
    string line;
    getline(trackFile, line);

    vector<pair<ftype, ftype>> trackData;
    vector<ftype> vertices;

    int nVertices = 0;

    size_t startPos = 0;
    vector<size_t> trackEndPositions;
    trackEndPositions.push_back(2);
    string vertexDelimiter = "]], [";
    while (true) {
        startPos = line.find(vertexDelimiter, startPos);
        if (startPos == string::npos) {
            break;
        }
        startPos += vertexDelimiter.length();
        trackEndPositions.push_back(startPos);
        nVertices += 1;
    }
    nVertices += 1; // for the last vertex that isn't counted

    string data;
    string trackArrayBeginStr = "[[";
    string trackArrayEndStr = "]]";
    string trackDelimiter = "], [";
    string trackErrorDelimiter = ", ";

    for (int i = 0; i < nVertices; ++i) {
        if (i != nVertices - 1) {
            data = line.substr(trackEndPositions[i], line.find(vertexDelimiter, trackEndPositions[i]) + 1 - trackEndPositions[i]);
        } else {
            data = line.substr(trackEndPositions[i], line.length() - trackEndPositions[i] - 1);
        }

        string vertex = data.substr(0, data.find(trackArrayBeginStr) - 2);
        string trackDataStr = data.substr(data.find(trackArrayBeginStr), data.find(trackArrayEndStr) + trackArrayEndStr.length() - data.find(trackArrayBeginStr));
        vertices.push_back(stod(vertex));

        size_t currentPos = 0;
        while (currentPos != string::npos) {
            currentPos = trackDataStr.find(trackArrayBeginStr) + trackArrayBeginStr.length();

            string track_x = trackDataStr.substr(currentPos, trackDataStr.find(trackErrorDelimiter) - currentPos);
            ftype x = stod(track_x);
            string track_errorx = trackDataStr.substr(trackDataStr.find(trackErrorDelimiter) + trackErrorDelimiter.length(),
                                                      trackDataStr.find(trackDelimiter, currentPos) - (trackDataStr.find(trackErrorDelimiter) + trackErrorDelimiter.length()));
            ftype errorx = stod(track_errorx);

            trackData.emplace_back(x, errorx);

            if (trackDataStr.find(trackDelimiter) == string::npos) {
                break;
            }

            trackDataStr = trackDataStr.substr(trackDataStr.find(trackDelimiter), trackDataStr.length() - trackDataStr.find(trackDelimiter));
            trackDataStr.replace(0, trackDelimiter.length(), trackArrayBeginStr);
        }
    }

    return event_t{nVertices, (int) trackData.size(), trackData};
}

// thank you O1 for the ARI implementation

inline ftype comb2(unsigned long x) {
    return (x < 2) ? 0.0 : (static_cast<ftype>(x) * (x - 1)) / 2.0;
}

ftype adjustedRandIndex(const vector<int>& labels_true, const vector<int>& labels_pred) {
    // 1. Basic sanity checks
    if (labels_true.size() != labels_pred.size() || labels_true.empty()) {
        return 0.0;
    }
    size_t N = labels_true.size();

    // 2. Map each unique label to an index
    unordered_map<int, int> true_label_to_index;
    unordered_map<int, int> pred_label_to_index;
    int next_true_index = 0;
    int next_pred_index = 0;
    
    for (auto lbl : labels_true) {
        if (true_label_to_index.find(lbl) == true_label_to_index.end()) {
            true_label_to_index[lbl] = next_true_index++;
        }
    }
    for (auto lbl : labels_pred) {
        if (pred_label_to_index.find(lbl) == pred_label_to_index.end()) {
            pred_label_to_index[lbl] = next_pred_index++;
        }
    }

    int n_true_clusters = next_true_index;
    int n_pred_clusters = next_pred_index;

    // 3. Build the contingency table
    vector<vector<unsigned long>> contingency(n_true_clusters, vector<unsigned long>(n_pred_clusters, 0ULL));
    
    for (size_t i = 0; i < N; ++i) {
        int true_idx = true_label_to_index[labels_true[i]];
        int pred_idx = pred_label_to_index[labels_pred[i]];
        contingency[true_idx][pred_idx]++;
    }

    // 4. Compute row sums, column sums
    vector<unsigned long> row_sums(n_true_clusters, 0ULL);
    vector<unsigned long> col_sums(n_pred_clusters, 0ULL);

    for (int i = 0; i < n_true_clusters; ++i) {
        for (int j = 0; j < n_pred_clusters; ++j) {
            row_sums[i] += contingency[i][j];
            col_sums[j] += contingency[i][j];
        }
    }

    // 5. Compute the components of ARI
    ftype indexVal = 0.0L;
    for (int i = 0; i < n_true_clusters; ++i) {
        for (int j = 0; j < n_pred_clusters; ++j) {
            indexVal += comb2(contingency[i][j]);
        }
    }

    ftype sum_comb2_row = 0.0L;
    for (int i = 0; i < n_true_clusters; ++i) {
        sum_comb2_row += comb2(row_sums[i]);
    }

    ftype sum_comb2_col = 0.0L;
    for (int j = 0; j < n_pred_clusters; ++j) {
        sum_comb2_col += comb2(col_sums[j]);
    }

    ftype comb2N = comb2(N);
    if (comb2N == 0.0L) {
        return 0.0;
    }

    ftype expectedIndex = (sum_comb2_row * sum_comb2_col) / comb2N;
    ftype maxIndex = 0.5L * (sum_comb2_row + sum_comb2_col);

    ftype denominator = maxIndex - expectedIndex;
    if (denominator < 1e-15L && denominator > -1e-15L) {
        // Degenerate case
        return 0.0;
    }

    ftype ARI = (indexVal - expectedIndex) / denominator;
    return static_cast<ftype>(ARI);
}