#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <assert.h>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <cmath>
#include <limits>

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


// returns a vector of the z positions of the vertices. same vertex indices as that of the assignment
// we simply average the z positions of the tracks assigned to each vertex
/*
vector<ftype> assignment_to_vertices(const vector<int> &assignment, const event_t &event) {
    vector<ftype> vertices(event.nV, 0);
    vector<int> counts(event.nV, 0);

    assert(assignment.size() == event.nT);

    for (int i = 0; i < event.nT; i++) {
        int vertex = assignment[i]; // index of vertex track i is assigned to
        vertices[vertex] += event.trackData[i].first; // z position of this track
        counts[vertex]++;
    }

    for (int i = 0; i < vertices.size(); i++) {
        if (counts[i] == 0) {
            vertices[i] = -1;
            continue;
        }
        vertices[i] /= counts[i];
    }

    // remove all -1s
    vertices.erase(remove(vertices.begin(), vertices.end(), -1), vertices.end());

    return vertices;
}
*/

vector<ftype> assignment_to_vertices(const vector<int> &assignment, const event_t &event) {
    vector<ftype> vertices(event.nV, 0);
    vector<int> counts(event.nV, 0);

    assert(assignment.size() == event.nT);

    for (int i = 0; i < event.nT; i++) {
        int vertex = assignment[i]; // index of vertex track i is assigned to
        
        // Check for invalid assignment (e.g., -1 or out of range)
        if (vertex < 0 || vertex >= event.nV) {
            // skip or handle as unassigned track
            continue;
        }
        
        vertices[vertex] += event.trackData[i].first; // z position of this track
        counts[vertex]++;
    }

    for (int i = 0; i < vertices.size(); i++) {
        if (counts[i] == 0) {
            vertices[i] = -1;
            continue;
        }
        vertices[i] /= counts[i];
    }

    // remove all -1s
    vertices.erase(remove(vertices.begin(), vertices.end(), -1), vertices.end());

    return vertices;
}


// both inputs should be sorted.
ftype dp_matching(const vector<ftype> &A, const vector<ftype> &B) {
    if (B.size() < A.size()) { // swap A and B if B is smaller
        return dp_matching(B, A);
    }

    // rest assumes A is smaller or equal to B

    // dp[i][j] is the minimum mse between the first i matches when considering the first j elems of B
    vector<vector<ftype>> dp(A.size() + 1, vector<ftype>(B.size() + 1, 0));

    for (int i = 0; i <= A.size(); i++) {
        for (int j = i; j <= B.size(); j++) {
            if (i == 0) {
                dp[i][j] = 0;
            } else if (j == 0) {
                dp[i][j] = 0;
            } else {
                ftype e = A[i - 1] - B[j - 1];
                if (i == j) { // we _have_ to match these
                    dp[i][j] = dp[i - 1][j - 1] + e * e;
                    continue;
                }
                dp[i][j] = min(dp[i - 1][j - 1] + e * e, // pick this matching
                               dp[i][j - 1]);            // or skip this matching
            }
        }
    }

    return sqrt(dp[A.size()][B.size()] / A.size()); // this is actually rmse
}


// this is a bit dubious because how do we know which vertex is supposed to be which?
// currently we use a dp matching algo to get the min mse.
// one could still argue that this metric may be biased towards solutions with missing elems.
ftype vertex_mse(const vector<ftype> &vertices, const event_t &event) {
    vector<ftype> vertices_copy = vertices;
    vector<ftype> event_vertices = event.vertices;

    sort(vertices_copy.begin(), vertices_copy.end());
    sort(event_vertices.begin(), event_vertices.end());

    // cout both
    cout << "vertices: ";
    for (int i = 0; i < vertices.size(); i++) {
        cout << vertices_copy[i] << " ";
    }
    cout << "\nevent vertices: ";
    for (int i = 0; i < event_vertices.size(); i++) {
        cout << event_vertices[i] << " ";
    }
    cout << "\n";

    return dp_matching(vertices_copy, event_vertices);

    // ftype mse = 0;
    // for (int i = 0; i < vertices_copy.size(); i++) {
    //     ftype e = vertices_copy[i] - event_vertices[i];
    //     mse += e * e;
    // }

    // return sqrt(mse / vertices_copy.size()); // this is actually rmse
}

ftype energy_from_assignment(const vector<int> &assignment, const QUBO &qubo, const int nT, const int nV) {
    solution_t solution(nT*nV, 0);

    for (int i = 0; i < assignment.size(); i++) {
        int vertex = assignment[i];
        // Critical: Ensure 'vertex' is within expected bounds
        if (vertex < 0 || vertex >= nV) {
            // Handle error: print a message, throw an exception, or adjust vertex
            // For debugging, a simple print is good:
            std::cerr << "Error: 'vertex' value " << vertex << " is out of bounds for nV=" << nV << std::endl;
            // You might even want to skip this iteration or exit.
            continue; // Skip current iteration if vertex is invalid
        }
        /*
        //cluster-major
        solution[nT * vertex + i] = 1; // Or the corrected index formula
        */
        //track-major
        solution[vertex + nV * i] = 1;
    }

    return qubo.evaluate(solution);
}

ftype print_score(const vector<int> &assignment, const event_t &event) {
    // print the ARI score of the assignment
    vector<int> true_labels(event.nT);
    for (int i = 0; i < event.nT; i++) {
        true_labels[i] = i / round(event.nT/event.nV); // event.nT tracks per vertex. TODO: magic number bad
    }

    // // debug print true labels and assignment
    // cout << "true labels: ";
    // for (int i = 0; i < true_labels.size(); i++) {
    //     cout << true_labels[i] << " ";
    // }
    // cout << endl;

    // cout << "assignment: ";
    // for (int i = 0; i < assignment.size(); i++) {
    //     cout << assignment[i] << " ";
    // }
    // cout << endl;
    
    cout<<"Number of tracks per vertex X Number of vertices = " << event.nT << ", Number of vertices = "<<event.nV<<endl;

    ftype ari = adjustedRandIndex(true_labels, assignment);

    cout << "ARI: " << ari << endl;
    
    return ari;
}

// https://arxiv.org/pdf/1903.08879
QUBO event_to_qubo(const event_t &event) {
    int nV = event.nV;
    int nT = event.nT;
    const auto& trackData = event.trackData;

    const ftype scale = 1.5;
    //const ftype scale = 1.;

    qubo_t qubo_map;

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

    /*
    auto idx = [nT](int track, int vertex) {
        //cluster-major
        return track + nT * vertex;
    };
    */
    
    auto idx = [nV](int track, int vertex) {
        //track-major
        return nV * track + vertex;
    };

    ftype lambda = 1.2;

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
                qubo_map[{idx(i, k), idx(j, k)}] += g(D_ij / max_D);
                //qubo_map[{idx(j, k), idx(i, k)}] += g(D_ij);
            }
        }
    }

    //lambda *= max_D/scale;
    //lambda /= max_D;
    //lambda /= log(max_D);
    //lambda *= exp(-1.0);
    //lambda *= scale/10;
    
    std::cout << "max D = " << max_D << endl;

    // penalty
    for (int i = 0; i < nT; ++i) {
        for (int k = 0; k < nV; ++k) {
            qubo_map[{idx(i, k), idx(i, k)}] -= lambda;
            for (int k2 = k + 1; k2 < nV; ++k2) {
                qubo_map[{idx(i, k), idx(i, k2)}] += 2 * lambda;
            }
        }
    }
    
    cout << "qubo num terms: " << qubo_map.size() << '\n';
    cout << "max possible: " << nT * nV << "^2\n";
    
    /*std::cout << "start" << endl;
    for (const auto& entry : qubo_map) {
        auto indices = entry.first;
        ftype value = entry.second;
        std::cout << "(" << indices.first << "," << indices.second << ") : " << value << '\n';
    }
    std::cout << "end" << endl;*/

    return QUBO(qubo_map);
}

ftype get_max_D(const event_t &event) {
    auto D = [](pair<ftype, ftype> i, pair<ftype, ftype> j) {
        return abs(i.first - j.first) / sqrt(i.second * i.second + j.second * j.second);
    };

    ftype max_D = 0.0;
    for (int i = 0; i < event.nT; i++) {
        for (int j = i + 1; j < event.nT; j++) {
            max_D = max(max_D, D(event.trackData[i], event.trackData[j]));
        }
    }

    cout << "max_D: " << max_D << '\n' << endl;

    return max_D;
}

ftype evaluate_full_OTF(const solution_t &x, const event_t &event, ftype max_D) {
    ftype value = 0.0;
    solution_t curr_x(x.size(), 0);
    for (int i = 0; i < x.size(); i++) {
        if (x[i]) {
            ftype delta = evaluate_diff_on_the_fly(curr_x, event, i, max_D);
            value += delta;
            curr_x[i] = 1;
        }
    }
    return value;
}

ftype evaluate_diff_on_the_fly(const solution_t &x, const event_t &event, int flip_idx, ftype max_D) {
    ftype lambda = 1.2;
    
    int nT = event.nT, nV = event.nV;
    int track = flip_idx % nT;
    int vertex = flip_idx / nT;

    auto D = [](pair<ftype, ftype> i, pair<ftype, ftype> j) {
        return abs(i.first - j.first) / sqrt(i.second * i.second + j.second * j.second);
    };

    const ftype scale = 1.5;
    //const ftype scale = 1.;

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

    /*
    auto idx = [nT](int track, int vertex) {
        //cluster-major
        return track + nT * vertex;
    };
    */
    
    auto idx = [nV](int track, int vertex) {
        //track-major
        return nV * track + vertex;
    };

    ftype diff = 0.0;
    for (int j = 0; j < nT; j++) {
        if (j == track) continue;
        int other_idx = idx(j, vertex);
        ftype term = g(D(event.trackData[track], event.trackData[j]) / max_D);
        //ftype term = g(D(event.trackData[track], event.trackData[j]));
        diff += term * x[other_idx];
    }
    for (int v = 0; v < nV; v++) {
        int current_idx = idx(track, v);
        if (v == vertex) {
            diff -= lambda;
        } else {
            diff += 2 * lambda * x[current_idx];
        }
    }
    return x[flip_idx] ? -diff : diff;
}


ftype ground_state(const QUBO &qubo, const event_t &event) {
    solution_t solution(event.nT*event.nV, 0);

    auto idx = [event](int track, int vertex) {
        /*
        //cluster-major
        return track + event.nT * vertex;
        */
        //track-major
        return event.nV * track + vertex;
    };

    for (int i = 0; i < event.nT; i++) {
        int vertex = i / round(event.nT/event.nV); //TODO: magic number bad
        solution[idx(i, vertex)] = 1;
    }
    return qubo.evaluate(solution);
}

/////////////////////// ==================== Dunn Index start ===================== /////////////////////////////

ftype track_distance(const pair<ftype, ftype>& i, const pair<ftype, ftype>& j) {
    if (i.second == 0 && j.second == 0) {
        return (i.first != j.first) ? std::numeric_limits<ftype>::infinity() : 0.0;
    }
    ftype denominator = std::sqrt(i.second * i.second + j.second * j.second);
    if (denominator == 0) {
        return std::numeric_limits<ftype>::infinity(); // Avoid division by zero
    }
    return std::abs(i.first - j.first);// / denominator;
}

ftype get_vertex_distance(ftype centroid1, ftype centroid2) {
    return std::abs(centroid1 - centroid2);
}

// Computes the maximum intra-cluster distance (diameter) for a given cluster
// 'cluster_track_indices' contains the global indices of tracks belonging to this cluster
ftype get_max_intra_cluster_D(const vector<int>& cluster_track_indices, const vector<pair<ftype, ftype>>& all_track_data) {
    if (cluster_track_indices.size() < 2) {
        return 0.0; // Diameter is 0 for clusters with less than 2 tracks
    }

    ftype max_d = 0.0;
    for (size_t i = 0; i < cluster_track_indices.size(); ++i) {
        for (size_t j = i + 1; j < cluster_track_indices.size(); ++j) {
            const auto& track1 = all_track_data[cluster_track_indices[i]];
            const auto& track2 = all_track_data[cluster_track_indices[j]];
            max_d = std::max(max_d, track_distance(track1, track2));
        }
    }
    return max_d;
}

// Computes the minimum inter-cluster distance between two clusters
// 'cluster1_track_indices' and 'cluster2_track_indices' contain global indices for their respective tracks
/*ftype get_min_inter_cluster_D(const vector<int>& cluster1_track_indices,const vector<int>& cluster2_track_indices,const vector<pair<ftype, ftype>>& all_track_data) {
    ftype min_d = std::numeric_limits<ftype>::infinity();
    for (int track_idx1 : cluster1_track_indices) {
        for (int track_idx2 : cluster2_track_indices) {
            const auto& track1 = all_track_data[track_idx1];
            const auto& track2 = all_track_data[track_idx2];
            min_d = std::min(min_d, track_distance(track1, track2));
        }
    }
    return min_d;
}*/


// --- Main Dunn Index Calculation Function ---
// --- Main Dunn Index Calculation Function (Modified to accept pre-calculated centroids) ---
// cluster_centroids: A vector where each element is the centroid of a cluster.
// The index of the centroid in this vector should correspond to the cluster_id.
ftype calculate_dunn_index(const std::vector<int>& assignment, const event_t& event, int nV, const std::vector<ftype>& cluster_centroids) {
    // 1. Organize tracks by cluster
    std::vector<std::vector<int>> clusters_track_indices(nV);
    for (int i = 0; i < event.nT; ++i) {
        int cluster_id = assignment[i];
        if (cluster_id >= 0 && cluster_id < nV) { // Ensure valid cluster ID
            clusters_track_indices[cluster_id].push_back(i);
        }
    }

    // Identify valid clusters (non-empty ones) and their corresponding centroid IDs
    std::vector<int> valid_cluster_ids;
    for (int i = 0; i < nV; ++i) {
        if (!clusters_track_indices[i].empty()) {
            valid_cluster_ids.push_back(i);
        }
    }

    // Check if there are at least 2 non-empty clusters
    int num_valid_clusters = valid_cluster_ids.size();
    if (num_valid_clusters < 2) {
        // Dunn Index is undefined for less than 2 clusters
        return std::numeric_limits<ftype>::infinity();
    }

    // 2. Find the maximum intra-cluster distance among all valid clusters
    ftype max_intra_d = 0.0; // This will store the max diameter found across all clusters
    for (int cluster_id : valid_cluster_ids) {
        max_intra_d = std::max(max_intra_d, get_max_intra_cluster_D(clusters_track_indices[cluster_id], event.trackData));
    }

    if (max_intra_d == 0) {
        // If the largest intra-cluster distance is 0, it means all tracks within
        // all clusters are identical. This implies perfectly compact clusters,
        // making the Dunn index infinite, indicating perfect clustering.
        return std::numeric_limits<ftype>::infinity();
    }

    // 3. Find the minimum inter-cluster distance (vertex distance) among all pairs of valid clusters
    ftype min_inter_d = std::numeric_limits<ftype>::infinity(); // This will store the min separation
    for (int k1_idx = 0; k1_idx < num_valid_clusters; ++k1_idx) {
        int cluster_id1 = valid_cluster_ids[k1_idx];
        for (int k2_idx = k1_idx + 1; k2_idx < num_valid_clusters; ++k2_idx) {
            int cluster_id2 = valid_cluster_ids[k2_idx];
            
            // Ensure centroids exist for these IDs before accessing
            if (cluster_id1 < cluster_centroids.size() && cluster_id2 < cluster_centroids.size()) {
                 min_inter_d = std::min(min_inter_d, get_vertex_distance(cluster_centroids[cluster_id1], cluster_centroids[cluster_id2]));
            } else {
                // Handle error: centroid for a valid cluster ID was not provided
                // This might indicate an issue with how cluster_centroids is populated.
                // For robustness, you might want to log this or return an error value.
                // For now, it will simply skip this pair if a centroid is missing.
            }
        }
    }
    
    // If no valid inter-cluster distances were found (e.g., due to centroid mismatch)
    // and min_inter_d remains infinity, it's an issue.
    // If it's infinity because there were < 2 valid clusters, that's handled above.

    // 4. Calculate the Dunn Index
    return min_inter_d / max_intra_d;
}

/////////////////////// ==================== Dunn Index end ===================== /////////////////////////////


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

    return event_t{nVertices, (int) trackData.size(), trackData, vertices};
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
