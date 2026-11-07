#include "dataset.h"

mt19937 rng(12345);
using namespace std;

vector<int> genWord(int K)
{
    uniform_int_distribution<int> bit(0, 1);
    vector<int> word(K);
    for (int i = 0; i < K; i++)
        word[i] = bit(rng);
    return word;
}

vector<int> randomPositions(int N, int count)
{
    vector<int> pos(N);
    iota(pos.begin(), pos.end(), 0);
    shuffle(pos.begin(), pos.end(), rng);
    pos.resize(count);
    return pos;
}

vector<int> encode(const vector<int>& word, const Matrix& G)
{
    vector<int> code(G[0].size(), 0);
    for (size_t i = 0; i < G.size(); i++)
    {
        if (word[i])
        {
            for (size_t j = 0; j < code.size(); j++)
                code[j] ^= G[i][j];
        }
    }
    return code;
}

vector<double> modulate(const vector<int>& code)
{
    vector<double> signal(code.size());
    for (size_t i = 0; i < code.size(); i++)
        signal[i] = (code[i] == 0) ? 1.0 : -1.0;
    return signal;
}

double sigmaFromEbN0(double ebn0dB, int K, int N)
{
    double rate = (double)K / N;
    return sqrt(1.0 / (2.0 * rate * pow(10.0, ebn0dB / 10.0)));
}

vector<double> awgn(const vector<double>& signal, double sigma)
{
    normal_distribution<double> normal(0.0, sigma);
    vector<double> noisy(signal.size());
    for (size_t i = 0; i < signal.size(); i++)
        noisy[i] = signal[i] + normal(rng);
    return noisy;
}

vector<double> gaussianNoise(const vector<double>& signal, double sigma, int mode, int errors, vector<int> errorPositions) // mode 1 - normal noise, mode 2 - noisy error bit
{
    normal_distribution<double> normal0(0.0, sigma);
    normal_distribution<double> normal2(2.0, sigma);

    if (errorPositions.empty())
        errorPositions = randomPositions(signal.size(), errors);
    errors = errorPositions.size();

    vector<bool> err_pos(signal.size(), false);
    for (int pos : errorPositions)
        err_pos[pos] = true;

    vector<double> noisy(signal.size());
    int error = 0;

    do
    {
        error = 0;
        for (size_t i = 0; i < signal.size(); i++)
        {
            double bit = signal[i];
            double noise = normal0(rng);
            double noisyBit = (mode == 1) ? bit + noise : bit;

            if (err_pos[i])
            {
                noise = normal2(rng);
                noisyBit = bit + noise;
                if (bit * noisyBit > 0)
                    noisyBit = bit - noise;
            }
            if (bit * noisyBit < 0)
                error++;
            noisy[i] = noisyBit;
        }
    } while (error != errors);

    return noisy;
}

void writeSet(const Matrix& G, double sigma1, double sigma2, int size, int mode, int gaussMode, int errors, const string& filename, const vector<int>& errorPositions)
{
    ofstream out(filename);
    out.precision(8);

    int K = G.size();
    uniform_real_distribution<double> uniformSigma(sigma1, sigma2);
    double sigma = (mode == 3) ? uniformSigma(rng) : sigma1;

    for (int n = 0; n < size; n++)
    {
        vector<int> word = genWord(K);
        vector<double> signal = gaussianNoise(modulate(encode(word, G)), sigma, gaussMode, errors, errorPositions);

        out << sigma << ' ';
        for (int bit : word)
            out << bit;
        out << ' ';
        for (double value : signal)
            out << value << ' ';
        out << '\n';

        if (mode == 2)
        {
            sigma += 0.01;
            if (sigma > sigma2)
                sigma = sigma1;
        }
        if (mode == 3)
            sigma = uniformSigma(rng);
    }
}

void getDataset(const Matrix& G, double sigma1Train, double sigma2Train, double sigma1Test, double sigma2Test,
    int datasetSize, int mode, int gaussMode, int errors,
    const string& trainFile, const string& testFile,
    const vector<int>& errorPositionsTrain, const vector<int>& errorPositionsTest)
{
    int trainSize = (int)(datasetSize * 0.7);
    int testSize = (int)(datasetSize * 0.3);

    writeSet(G, sigma1Train, sigma2Train, trainSize, mode, gaussMode, errors, trainFile, errorPositionsTrain);
    writeSet(G, sigma1Test, sigma2Test, testSize, mode, gaussMode, errors, testFile, errorPositionsTest);
}