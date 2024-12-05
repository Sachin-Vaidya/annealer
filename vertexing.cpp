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
                return vector<int>(); // invalid solution
            }
            int track = i % nT;
            int vertex = i / nT;
            assignment[track] = vertex;
        }
    }

    if (find(assignment.begin(), assignment.end(), -1) != assignment.end()) {
        cout << "not all tracks assigned!\n";
        return vector<int>(); // invalid solution
        return assignment;
    }

    return assignment;
}


// https://arxiv.org/pdf/1903.08879
qubo_t event_to_qubo(const event_t &event) {
    int nV = event.nV;
    int nT = event.nT;
    const auto& trackData = event.trackData;

    map<pair<int, int>, double> qubo_map;

    auto D = [](pair<double, double> i, pair<double, double> j) {
        return abs(i.first - j.first) / sqrt(i.second * i.second + j.second * j.second);
    };

    auto g = [](double x, double m = 5) {
        // return 1.0 - exp(-m * x);
        return x;
    };

    auto idx = [nT](int track, int vertex) {
        return track + nT * vertex;
    };

    double lambda = 1.2;
    // double lambda = 1.2;

    double max_D = 0.0;

    for (int k = 0; k < nV; ++k) {
        for (int i = 0; i < nT; ++i) {
            for (int j = i + 1; j < nT; ++j) {
                double D_ij = D(trackData[i], trackData[j]);
                max_D = max(max_D, D_ij);
                qubo_map[{idx(i, k), idx(j, k)}] += g(D_ij);
            }
        }
    }

    lambda *= max_D;

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
    return qubo;
}


// source: andrew wildridge. lightly modified
event_t loadTracks(string filename) {
    ifstream trackFile(filename.c_str());
    string line;
    getline(trackFile, line);

    vector<pair<double, double>> trackData;
    vector<double> vertices;

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
            double x = stod(track_x);
            string track_errorx = trackDataStr.substr(trackDataStr.find(trackErrorDelimiter) + trackErrorDelimiter.length(),
                                                      trackDataStr.find(trackDelimiter, currentPos) - (trackDataStr.find(trackErrorDelimiter) + trackErrorDelimiter.length()));
            double errorx = stod(track_errorx);

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