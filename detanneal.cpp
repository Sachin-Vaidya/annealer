#include <cstdio>
#include <math.h>
#include <utility>
#include <vector>
#include <stdio.h>
#include <ctime>

#include "vertexing.hh"
#include "detanneal.hh"

using namespace std;

ftype getDistortion(ftype& x, ftype& errorX, ftype& y) {
    ftype sigma = (errorX != 0) ? (x - y) / errorX : 0; // Added zero-check
    return sigma * sigma;
}

ftype getClusterProbability(int& j, const int& N, vector<vector<ftype>>& associationMatrix) {
    ftype probability = 0.0;
    for (int i = 0; i < N; ++i) {
        probability += associationMatrix[j][i];
    }
    probability /= N;
    return probability;
}

vector<ftype> getCriticalTemperatures(vector<ftype>& X, vector<ftype>& errorX, vector<ftype>& Y, vector<vector<ftype>>& associationMatrix, const int& N) {
    vector<ftype> criticalTemperatures = {};
    for (int j = 0; j < Y.size(); ++j) {
        ftype prob_y = getClusterProbability(j, N, associationMatrix);
        ftype criticalTemperature = 0;
        if (prob_y == 0) { // Added zero-check
            criticalTemperatures.push_back(0);
            continue;
        }
        for (int i = 0; i < N; ++i) {
            criticalTemperature += associationMatrix[j][i] * getDistortion(X[i], errorX[i], Y[j]);
        }
        criticalTemperature /= (N * prob_y); // Changed to Code 1’s formulation for consistency
        criticalTemperatures.push_back(criticalTemperature);
    }
    return criticalTemperatures;
}

bool merge(vector<ftype>& Y, vector<ftype>& X, vector<ftype>& errorX, vector<ftype>& clusterProbabilities, vector<vector<ftype>>& associationMatrix, ftype& T, const int& N) {
    for (int j = 0; (j + 1) < Y.size(); ++j) {
        for (int k = j + 1; k < Y.size(); ++k) {
            if (fabs(Y[j] - Y[k]) < 0.002) {
                ftype totalProbability = clusterProbabilities[j] + clusterProbabilities[k];
                ftype newCentroid = (totalProbability > 0)
                    ? (clusterProbabilities[j] * Y[j] + clusterProbabilities[k] * Y[k]) / totalProbability
                    : 0.5 * (Y[j] + Y[k]);
                ftype newCriticalTemp = 0.0;
                for (int i = 0; i < N; ++i) {
                    newCriticalTemp += (associationMatrix[j][i] + associationMatrix[k][i]) * getDistortion(X[i], errorX[i], newCentroid);
                }
                newCriticalTemp = (totalProbability > 0) ? newCriticalTemp / (N * totalProbability) : 0; // Added zero-check and changed to Code 1’s formulation
                if (newCriticalTemp > T) {
                    continue;
                }
                clusterProbabilities[j] = totalProbability;
                Y[j] = newCentroid;
                for (int i = 0; i < N; ++i) {
                    associationMatrix[j][i] += associationMatrix[k][i];
                }
                Y.erase(Y.begin() + k);
                clusterProbabilities.erase(clusterProbabilities.begin() + k);
                associationMatrix.erase(associationMatrix.begin() + k);
                return true;
            }
        }
    }
    return false;
}

ftype getPartitionFunction(ftype& x, ftype& errorx, vector<ftype>& Y, ftype& T, vector<ftype>& clusterProbabilities, vector<ftype>& partitionComponents) {
    ftype partition_func = 0.0;
    for (int i = 0; i < Y.size(); ++i) {
        ftype partitionComponent = clusterProbabilities[i] * exp(-1 * getDistortion(x, errorx, Y[i]) / T);
        partitionComponents.push_back(partitionComponent);
        partition_func += partitionComponent;
    }
    return partition_func;
}

void updateAssociationProbabilityMatrix(vector<vector<ftype>>& associationMatrix, vector<ftype>& X, vector<ftype>& errorX, vector<ftype>& Y, ftype& T, const int& N, vector<ftype>& clusterProbabilities, vector<ftype>& partitionComponents) {
    for (int i = 0; i < N; ++i) {
        ftype partitionFunc = getPartitionFunction(X[i], errorX[i], Y, T, clusterProbabilities, partitionComponents);
        for (int j = 0; j < Y.size(); ++j) {
            associationMatrix[j][i] = (partitionFunc > 0) ? partitionComponents[j] / partitionFunc : 0; // Added zero-check
        }
        partitionComponents.clear();
    }
}

ftype getCentroid(vector<ftype>& X, vector<ftype>& jthClusterAssociationProbs, ftype& jthClusterProbability, const int& N) {
    ftype centroid = 0.0;
    for (int i = 0; i < N; ++i) {
        centroid += X[i] * jthClusterAssociationProbs[i];
    }
    return (jthClusterProbability > 0) ? centroid / (N * jthClusterProbability) : 0; // Added zero-check
}

vector<ftype> update(vector<ftype>& clusterCentroids, vector<vector<ftype>>& associationMatrix, vector<ftype>& clusterProbabilities, vector<ftype>& X, vector<ftype>& errorX, ftype& T, const int& N, vector<ftype>& partitionComponents) {
    updateAssociationProbabilityMatrix(associationMatrix, X, errorX, clusterCentroids, T, N, clusterProbabilities, partitionComponents);
    vector<ftype> deltas = {};
    for (int j = 0; j < clusterProbabilities.size(); ++j) {
        clusterProbabilities[j] = getClusterProbability(j, N, associationMatrix);
        ftype newCentroid = getCentroid(X, associationMatrix[j], clusterProbabilities[j], N);
        ftype deltaCentroid = clusterCentroids[j] - newCentroid;
        deltas.push_back(deltaCentroid * deltaCentroid);
        clusterCentroids[j] = newCentroid;
    }
    return deltas;
}

bool split(vector<ftype>& X, vector<ftype>& errorX, vector<ftype>& clusterCentroids, vector<vector<ftype>>& associationMatrix, vector<ftype>& clusterProbabilities, ftype& T, const int& N, ftype delta = 1e-3) {
    vector<ftype> criticalTemps = getCriticalTemperatures(X, errorX, clusterCentroids, associationMatrix, N);
    bool split = false;
    for (int k = 0; k < criticalTemps.size(); ++k) {
        if (T <= criticalTemps[k]) {
            split = true;
            ftype leftClusterProb = 0.0;
            ftype rightClusterProb = 0.0;
            ftype leftTotalWeight = 0.0;
            ftype rightTotalWeight = 0.0;
            ftype leftCentroid = 0.0;
            ftype rightCentroid = 0.0;
            for (int i = 0; i < N; ++i) {
                ftype probabilty = associationMatrix[k][i];
                ftype errorx = errorX[i];
                ftype x = X[i];
                ftype weight = (errorx != 0) ? probabilty / (errorx * errorx) : 0; // Added zero-check
                if (x < clusterCentroids[k]) {
                    leftClusterProb += probabilty;
                    leftTotalWeight += weight;
                    leftCentroid += weight * x;
                } else {
                    rightClusterProb += probabilty;
                    rightTotalWeight += weight;
                    rightCentroid += weight * x;
                }
            }
            leftCentroid = (leftTotalWeight > 0) ? leftCentroid / leftTotalWeight : (clusterCentroids[k] - delta);
            rightCentroid = (rightTotalWeight > 0) ? rightCentroid / rightTotalWeight : (clusterCentroids[k] + delta);
            if (rightCentroid - leftCentroid > delta) {
                ftype totalProb = leftClusterProb + rightClusterProb;
                ftype leftProbRatio = (totalProb > 0) ? leftClusterProb / totalProb : 0.5; // Added zero-check
                ftype rightProbRatio = (totalProb > 0) ? rightClusterProb / totalProb : 0.5; // Added zero-check
                clusterProbabilities.push_back(leftProbRatio * clusterProbabilities[k]);
                vector<ftype> newAssociationProbs = {};
                for (int i = 0; i < N; ++i) {
                    newAssociationProbs.push_back(leftProbRatio * associationMatrix[k][i]);
                    associationMatrix[k][i] = rightProbRatio * associationMatrix[k][i];
                }
                associationMatrix.push_back(newAssociationProbs);
                clusterCentroids.push_back(leftCentroid);
                clusterProbabilities[k] = rightProbRatio * clusterProbabilities[k];
                clusterCentroids[k] = rightCentroid;
            }
        }
    }
    return split;
}

pair<vector<int>, vector<ftype>> runDA(event_t event, int nSWEEPS) {
    vector<pair<ftype, ftype>> data = event.trackData;
    vector<ftype> X(data.size());
    vector<ftype> errorX(data.size());
    for (int i = 0; i < data.size(); ++i) {
        X[i] = data[i].first;
        errorX[i] = data[i].second;
    }
    std::clock_t time_start = std::clock();
    ftype Tmin = 4;
    ftype betaMax = 1.0 / Tmin;
    ftype Tstop = 1.0;
    ftype coolingFactor = 0.6;
    ftype nSweeps = nSWEEPS;
    bool useLinearCooling = true;
    ftype delta = 3.3e-5;
    int maxIterations = 350;
    ftype convergenceCriteria = 1e-9;
    int Kmax = event.nV; // Changed to use event.nV instead of hardcoded 2
    int N = X.size();
    vector<ftype> clusterProbabilities = { 1.0 };
    vector<ftype> partitionComponents = {};
    vector<ftype> clusterCentroids = { 0.0 };
    for (int j = 0; j < clusterCentroids.size(); ++j) {
        for (int i = 0; i < N; ++i) {
            clusterCentroids[j] += X[i];
        }
        clusterCentroids[j] /= N;
    }
    vector<vector<ftype>> associationMatrix = { vector<ftype>(N, 1.0) };
    ftype T = getCriticalTemperatures(X, errorX, clusterCentroids, associationMatrix, N)[0];
    ftype beta = 1.0 / T;
    ftype deltaBeta = (betaMax - beta) / nSweeps;
    while (T > Tmin) {
        for (int n = 0; n < maxIterations; ++n) {
            vector<ftype> deltas = update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
            ftype sum = 0.0;
            for (int j = 0; j < deltas.size(); ++j) {
                sum += deltas[j];
            }
            if (sum < convergenceCriteria) {
                break;
            }
        }
        while (merge(clusterCentroids, X, errorX, clusterProbabilities, associationMatrix, T, N)) {
            update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
        }
        split(X, errorX, clusterCentroids, associationMatrix, clusterProbabilities, T, N, delta);
        if (useLinearCooling) {
            beta += deltaBeta;
            T = 1.0 / beta;
        } else {
            T *= coolingFactor;
        }
    }
    update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
    while (merge(clusterCentroids, X, errorX, clusterProbabilities, associationMatrix, T, N)) {
        update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
    }
    unsigned int ntry = 0;
    while (split(X, errorX, clusterCentroids, associationMatrix, clusterProbabilities, T, N, delta) && ntry++ < 10 && clusterCentroids.size() < Kmax) { // Added Kmax check
        for (int n = 0; n < maxIterations; ++n) {
            vector<ftype> deltas = update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
            ftype sum = 0.0;
            for (int j = 0; j < deltas.size(); ++j) {
                sum += deltas[j];
            }
            if (sum < convergenceCriteria) {
                break;
            }
        }
        while (merge(clusterCentroids, X, errorX, clusterProbabilities, associationMatrix, T, N)) {
            update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
        }
    }
    T = Tstop;
    for (int n = 0; n < maxIterations; ++n) {
        vector<ftype> deltas = update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
        ftype sum = 0.0;
        for (int j = 0; j < deltas.size(); ++j) {
            sum += deltas[j];
        }
        if (sum < convergenceCriteria) {
            break;
        }
    }
    while (merge(clusterCentroids, X, errorX, clusterProbabilities, associationMatrix, T, N)) {
        update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
    }
    vector<ftype> deltas = update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
    std::clock_t time_stop = std::clock();
    printf("runDA completed in %.3f ms\n", (time_stop - time_start) * 1000.0 / CLOCKS_PER_SEC);
    vector<int> clusterAssignment(N, -1);
    for (int i = 0; i < N; ++i) {
        ftype maxProb = 0.0;
        int maxCluster = -1;
        for (int j = 0; j < clusterProbabilities.size(); ++j) {
            if (associationMatrix[j][i] > maxProb) {
                maxProb = associationMatrix[j][i];
                maxCluster = j;
            }
        }
        clusterAssignment[i] = maxCluster;
    }
    return {clusterAssignment, clusterCentroids};
}
