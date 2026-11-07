#pragma once
#include <random>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <fstream>

using namespace std;
using Matrix = vector<vector<int>>;

extern mt19937 rng;

vector<int> genWord(int K);
vector<int> randomPositions(int N, int count);
vector<int> encode(const vector<int>& word, const Matrix& G);
vector<double> modulate(const vector<int>& code);
double sigmaFromEbN0(double ebn0dB, int K, int N);
vector<double> awgn(const vector<double>& signal, double sigma);
vector<double> gaussianNoise(const vector<double>& signal, double sigma, int mode = 1, int errors = 1, vector<int> errorPositions = vector<int>());
void getDataset(const Matrix& G, double sigma1Train, double sigma2Train, double sigma1Test, double sigma2Test,
    int datasetSize = 36000, int mode = 1, int gaussMode = 1, int errors = 1,
    const string& trainFile = "trainset.txt", const string& testFile = "testset.txt",
    const vector<int>& errorPositionsTrain = vector<int>(), const vector<int>& errorPositionsTest = vector<int>());