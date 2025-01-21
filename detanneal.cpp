// code here lightly edited from a private deterministic annealer implementation from andrew wildridge.

#include <cstdio>
#include <math.h>
#include <utility>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <algorithm>
#include <ctime>

#include "vertexing.hh"
#include "detanneal.hh"

using namespace std;

ftype getDistortion(ftype& x, ftype& errorX, ftype& y) {
	ftype sigma = (x - y) / errorX;
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
		for (int i = 0; i < N; ++i) {
			criticalTemperature += associationMatrix[j][i] * getDistortion(X[i], errorX[i], Y[j]) / prob_y;
		}
		criticalTemperature /= N;
		criticalTemperatures.push_back(criticalTemperature);
	}

	return criticalTemperatures;
}

bool merge(vector<ftype>& Y, vector<ftype>& X, vector<ftype>& errorX, vector<ftype>& clusterProbabilities, vector<vector<ftype>>& associationMatrix, ftype& T, const int& N) {
	for (int j = 0; (j + 1) < Y.size(); ++j) {
		for (int k = j + 1; k < Y.size(); ++k) {
			if (fabs(Y[j] - Y[k]) < 0.002) { // Need to merge
				// First check to make sure that the merged critical temperature isn't higher than current temperature
				ftype totalProbability = clusterProbabilities[j] + clusterProbabilities[k];

				ftype newCentroid = 0.0;
				if (totalProbability > 0) {
					newCentroid = (clusterProbabilities[j] * Y[j] + clusterProbabilities[k] * Y[k]) / totalProbability;
				}
				else {
					newCentroid = 0.5 * (Y[j] + Y[k]);
				}

				ftype newCriticalTemp = 0.0;
				for (int i = 0; i < N; ++i) {
					newCriticalTemp += (associationMatrix[j][i] + associationMatrix[k][i]) * getDistortion(X[i], errorX[i], newCentroid) / totalProbability;
				}
				newCriticalTemp /= N;

				//printf("Attempting to merge clusters %i and %i, new critical temperature is %f. \n", j, k, newCriticalTemp);

				// If merged cluster critical temperature is greater than the current temperature, skip this merge
				if (newCriticalTemp > T) {
					continue;
				}

				clusterProbabilities[j] = totalProbability;
				Y[j] = newCentroid;
				for (int i = 0; i < N; ++i) {
					associationMatrix[j][i] += associationMatrix[k][i];
				}
				//printf("Merged clusters %i and %i \n", j, k);
				// Now delete the merged cluster
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
			associationMatrix[j][i] = partitionComponents[j] / partitionFunc;
		}
		partitionComponents.clear();
	}
}

ftype getCentroid(vector<ftype>& X, vector<ftype>& jthClusterAssociationProbs, ftype& jthClusterProbability, const int& N) {
	ftype centroid = 0.0;
	for (int i = 0; i < N; ++i) {
		centroid += X[i] * jthClusterAssociationProbs[i];
	}
	centroid /= (N * jthClusterProbability);

	return centroid;
}

/**
 * Updates the association matrix, cluster centroids, and cluster probabilities while returning the squared difference in the cluster centroids
 * @param clusterCentroids The vector containing all of the cluster centroids
 * @param associationMatrix The matrix telling what track is associated with what cluster
 * @param clusterProbabilities The vector containing all of the probabilities that we would split into each cluster
 * @param X The data that is being clustered
 * @param errorX The error on the data being clustered
 * @param T The current temperature
 *
 * @return The vector containing the squared difference in cluster centroids between old and new centroids
 */
vector<ftype> update(vector<ftype>& clusterCentroids, vector<vector<ftype>>& associationMatrix, vector<ftype>& clusterProbabilities, vector<ftype>& X, vector<ftype>& errorX, ftype& T, const int& N, vector<ftype>& partitionComponents) {
	// Compute association probabilities first
	updateAssociationProbabilityMatrix(associationMatrix, X, errorX, clusterCentroids, T, N, clusterProbabilities, partitionComponents);

	// Now recompute cluster probabilities and cluster centroids
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
	//stable_sort(criticalTemps.begin(), criticalTemps.end(), std::greater<ftype>() );
	/*printf("The current cluster critical temperatures are: \n");
	for (int i = 0; i < criticalTemps.size(); ++i) {
		printf("\t Temperature %i: %f \n", i, criticalTemps[i]);
	}
	printf("\n");*/

	bool split = false;

	for (int k = 0; k < criticalTemps.size(); ++k) {
		if (T <= criticalTemps[k]) { // need to split that cluster
			split = true;
			//printf("Splitting the %ith cluster. \n", k);
			ftype leftClusterProb = 0.0; // new cluster formed from tracks whose z < old centroid
			ftype rightClusterProb = 0.0; // new cluster formed from tracks whose z > old centroid

			ftype leftTotalWeight = 0.0;
			ftype rightTotalWeight = 0.0;
			ftype leftCentroid = 0.0;
			ftype rightCentroid = 0.0;

			for (int i = 0; i < N; ++i) {
				ftype probabilty = associationMatrix[k][i];
				ftype errorx = errorX[i];
				ftype x = X[i];
				ftype weight = probabilty / (errorx * errorx);

				if (x < clusterCentroids[k]) {
					leftClusterProb += probabilty;
					leftTotalWeight += weight;
					leftCentroid += weight * x;
				}
				else {
					rightClusterProb += probabilty;
					rightTotalWeight += weight;
					rightCentroid += weight * x;
				}
			}

			if (leftTotalWeight > 0) {
				leftCentroid = leftCentroid / leftTotalWeight;
			}
			else {
				leftCentroid = clusterCentroids[k] - delta;
			}

			if (rightTotalWeight > 0) {
				rightCentroid = rightCentroid / rightTotalWeight;
			}
			else {
				rightCentroid = clusterCentroids[k] + delta;
			}

			// TODO: Reduce split size if there is not enough room
			// reduce split size if there is not enough room
			/*if( ( ik   > 0       ) && ( y[ik-1].z>=z1 ) ){ z1=0.5*(y[ik].z+y[ik-1].z); }
			if( ( ik+1 < y.size()) && ( y[ik+1].z<=z2 ) ){ z2=0.5*(y[ik].z+y[ik+1].z); }*/

			if (rightCentroid - leftCentroid > delta) {
				clusterProbabilities.push_back(leftClusterProb * clusterProbabilities[k] / (leftClusterProb + rightClusterProb));

				vector<ftype> newAssociationProbs = {};
				for (int i = 0; i < N; ++i) {
					newAssociationProbs.push_back(leftClusterProb * associationMatrix[k][i] / (leftClusterProb + rightClusterProb));
				}
				associationMatrix.push_back(newAssociationProbs);
				clusterCentroids.push_back(leftCentroid);

				clusterProbabilities[k] = rightClusterProb * clusterProbabilities[k] / (leftClusterProb + rightClusterProb);
				clusterCentroids[k] = rightCentroid;
				for (int i = 0; i < N; ++i) {
					associationMatrix[k][i] = rightClusterProb * associationMatrix[k][i] / (leftClusterProb + rightClusterProb);
				}
			}

			// old splitting code
			/*
			clusterCentroids.push_back(clusterCentroids[k] + delta);
			vector<ftype> newAssociationProbs = {};
			for (int i = 0; i < N; ++i) {
				newAssociationProbs.push_back(associationMatrix[k][i] / 2.0);
			}
			associationMatrix.push_back(newAssociationProbs);
			clusterProbabilities.push_back(clusterProbabilities[k] / 2.0);
			clusterProbabilities[k] /= 2.0;
			for (int i = 0; i < N; ++i) {
				associationMatrix[k][i] /= 2.0;
			}*/
		}
	}
	return split;
}


vector<int> runDA(event_t event) {
    vector<pair<ftype, ftype>> data = event.trackData;

    vector<ftype> X(data.size());
    vector<ftype> errorX(data.size());

    for (int i = 0; i < data.size(); ++i) {
        X[i] = data[i].first;
        errorX[i] = data[i].second;
    }

	// vector<ftype> X = data.trackData
	// vector<ftype> errorX = data.second;

	std::clock_t time_start = std::clock();
	// Step 1: Set Limts
	ftype Tmin = 4;
	ftype betaMax = 1.0 / Tmin;
	ftype Tstop = 1.0;
	ftype coolingFactor = 0.6;
	ftype nSweeps = 20000;
	bool useLinearCooling = true; // uses linear cooling in beta (1 / T is linearly increased)
	ftype delta = 3.3e-5;
	int maxIterations = 350;
	ftype convergenceCriteria = 1e-9; // squared difference between new centroid and old centroid
	int Kmax = 2; //TODO: Set this to nVertices
	int N = X.size();

	// Step 2: Initialize
	vector<ftype> clusterProbabilities = { 1.0 };
	vector<ftype> partitionComponents = {};

	vector<ftype> clusterCentroids = { 0.0 };
	for (int j = 0; j < clusterCentroids.size(); ++j) {
		for (int i = 0; i < N; ++i) {
			clusterCentroids[j] += X[i];
		}
		clusterCentroids[j] /= N;
	}

	vector<vector<ftype>> associationMatrix = { {} };
	for (int i = 0; i < N; ++i) {
		associationMatrix[0].push_back(1.0);
	}

	// Set initial temperature to first critical temperature
	ftype T = getCriticalTemperatures(X, errorX, clusterCentroids, associationMatrix, N)[0];
	ftype beta = 1.0 / T;

	ftype deltaBeta = (betaMax - beta) / nSweeps;
	//printf("The first critical temperature is at %f \n", T);

	// Annealing Loop
	while (T > Tmin) {
		//printf("Current temperature is %f. \n", T);

		// Get the state into equilibrium
		for (int n = 0; n < maxIterations; ++n) {
			vector<ftype> deltas = update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);

			// Check for convergence
			ftype sum = 0.0;
			for (int j = 0; j < deltas.size(); ++j) {
				sum += deltas[j];
			}

			if (sum < convergenceCriteria) {
				break;
			}
		}

		// Check for merging
		while (merge(clusterCentroids, X, errorX, clusterProbabilities, associationMatrix, T, N)) {
			update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
		}

		/*
				printf("The current centroid locations are: \n");
				for (int i = 0; i < clusterCentroids.size(); ++i) {
					printf("\t Centroid %i: %f \n", i, clusterCentroids[i]);
				}
				printf("\n");

				printf("The current centroid probabilities are: \n");
				for (int i = 0; i < clusterProbabilities.size(); ++i) {
					printf("\t Probability %i: %f \n", i, clusterProbabilities[i]);
				}
				printf("\n");*/

				// Step 7 split clusters that have a critical temperature above the new temperature
		split(X, errorX, clusterCentroids, associationMatrix, clusterProbabilities, T, N, delta);

		// Step 6 cool the temperature
		if (useLinearCooling) {
			beta += deltaBeta;
			T = 1.0 / beta;
		}
		else {
			T *= coolingFactor;
		}

		/*if (T < 10) {
			printf("The current association probabilities: \n");
			for (int j = 0; j < associationMatrix.size(); ++j) {
				for (int i = 0; i < associationMatrix[j].size(); ++i) {
					printf("The probability that the %ith track is associated with the %ith cluster is %f. \n", i, j, associationMatrix[j][i]);
				}
			}
		}*/
	}

	// Do your final splitting before you anneal down to the final/stopping temperature
	// There is no splitting after this

	update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents); // update first
	// Check for merging
	while (merge(clusterCentroids, X, errorX, clusterProbabilities, associationMatrix, T, N)) {
		update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
	}
	unsigned int ntry = 0;
	while (split(X, errorX, clusterCentroids, associationMatrix, clusterProbabilities, T, N, delta) && ntry++ < 10) {
		// Get the state into equilibrium
		for (int n = 0; n < maxIterations; ++n) {
			vector<ftype> deltas = update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);

			// Check for convergence
			ftype sum = 0.0;
			for (int j = 0; j < deltas.size(); ++j) {
				sum += deltas[j];
			}

			if (sum < convergenceCriteria) {
				break;
			}
		}

		merge(clusterCentroids, X, errorX, clusterProbabilities, associationMatrix, T, N);
		update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents); // update first
	}

	T = Tstop;
	// Step 5
	// Get the state into equilibrium for the final temperature
	for (int n = 0; n < maxIterations; ++n) {
		vector<ftype> deltas = update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);

		// Check for convergence
		ftype sum = 0.0;
		for (int j = 0; j < deltas.size(); ++j) {
			sum += deltas[j];
		}

		if (sum < convergenceCriteria) {
			break;
		}
	}

	// Check for merging at the end
	while (merge(clusterCentroids, X, errorX, clusterProbabilities, associationMatrix, T, N)) {
		vector<ftype> deltas = update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);
	}

	// do a final update on paramaters since some might have merged
	vector<ftype> deltas = update(clusterCentroids, associationMatrix, clusterProbabilities, X, errorX, T, N, partitionComponents);

	std::clock_t time_stop = std::clock();

	printf("The current centroid locations are: \n");
	for (int i = 0; i < clusterCentroids.size(); ++i) {
		printf("\t Centroid %i: %f \n", i, clusterCentroids[i]);
	}
	printf("\n");

	printf("The current centroid probabilities are: \n");
	for (int i = 0; i < clusterProbabilities.size(); ++i) {
		printf("\t Probability %i: %f \n", i, clusterProbabilities[i]);
	}
	printf("\n");


    // // print association matrix to 1 decimal place
    // for (int j = 0; j < associationMatrix.size(); ++j) {
    //     for (int i = 0; i < associationMatrix[j].size(); ++i) {
    //         printf("%.2f ", associationMatrix[j][i]);
    //     }
    //     printf("\n");
    // }

    // // print size of association matrix 
    // cout << associationMatrix.size() << " " << associationMatrix[0].size() << endl;


	// //printf("%i", clusterProbabilities.size());
	// std::ofstream responseFile("serializableResponse.json");
	// ftype energy = 0.0;
	// responseFile << "[";
	// responseFile << "[" << energy << ", [";
	// vector<string> bitstring;
	// for (int i = 0; i < N; ++i) {
	// 	for (int j = 0; j < clusterProbabilities.size(); ++j) {
	// 		responseFile << "\"" << ((int)round(associationMatrix[j][i])) << "\"";
	// 		if (i == N - 1 && j == clusterProbabilities.size() - 1) {
	// 			responseFile << "]";
	// 		}
	// 		else {
	// 			responseFile << ", ";
	// 		}
	// 	}
	// }
	// responseFile << "], " << (time_stop - time_start) << "]";

    // assoc matrix contains mapping from cluster number -> probability of track being in that cluster
    // we will process it to create a vector of length track count, where each element is the cluster number

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

    // for (int i = 0; i < clusterAssignment.size(); ++i) {
    //     printf("%i, ", clusterAssignment[i]);
    // }

    return clusterAssignment;
}