#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include<aff3ct.hpp>
#include <random>

mt19937 rng(random_device{}());
using namespace std;
using Sparse = aff3ct::tools::Sparse_matrix;
using Handler = aff3ct::tools::LDPC_matrix_handler;
using NMS = aff3ct::tools::Update_rule_NMS<float>;
using DecoderNMS = aff3ct::module::Decoder_LDPC_BP_flooding<int, float, NMS>;

using Matrix = vector<vector<int>>;

Matrix parseQC(const std::string& filename)
{
    std::ifstream file(filename);

    int z, n, m;
    file >> z >> n >> m;

    Matrix H(m * z, std::vector<int>(n * z, 0));

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int shift;
            file >> shift;
            if (shift != -1)
            {
                for (int k = 0; k < z; k++)
                {
                    int col = (k + shift) % z;
                    H[i * z + k][j * z + col] = 1;
                }
            }
        }
    }
    return H;
}

Matrix parseAlist(const string& filename)
{
    ifstream file(filename);

    int n, m, max_col_deg, max_row_deg;
    file >> n >> m;
    file >> max_col_deg >> max_row_deg;

    vector<int> col_weights(n);
    for (int i = 0; i < n; ++i)
        file >> col_weights[i];

    vector<int> row_weights(m);
    for (int i = 0; i < m; ++i)
        file >> row_weights[i];

    Matrix H(m, vector<int>(n, 0));

    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < max_col_deg; k++)
        {
            int x;
            file >> x;
            if (k < col_weights[j])
                H[x - 1][j] = 1;
        }
    }
    return H;
}

Sparse toSparse(const Matrix& H)
{
    Sparse S(H.size(), H[0].size());
    for (size_t i = 0; i < H.size(); i++)
    {
        for (size_t j = 0; j < H[i].size(); j++)
        {
            if (H[i][j])
                S.add_connection(i, j);
        }
    }
    return S;
}

Matrix toDense(const Sparse& S)
{
    Matrix M(S.get_n_rows(), vector<int>(S.get_n_cols(), 0));
    for (size_t i = 0; i < S.get_n_rows(); i++)
    {
        for (auto j : S.get_row_to_cols()[i])
            M[i][j] = 1;
    }
    return M;
}

Matrix get_G(const Matrix& H, vector<uint32_t>& infoPos)
{
    Sparse Hs = toSparse(H).turn(aff3ct::tools::Matrix::Way::HORIZONTAL);
    Sparse G = Handler::transform_H_to_G_identity(Hs, infoPos);
    return toDense(G.turn(aff3ct::tools::Matrix::Way::HORIZONTAL));
}

vector<int> decodeNMS(DecoderNMS& decoder, const vector<float>& code, int K)
{
    vector<int> bits(K);
    decoder.decode_siho(code.data(), bits.data());
    return bits;
}

void saveMatrix(const Matrix& H, const string& filename)
{
    ofstream out(filename);
    for (const auto& row : H)
    {
        for (size_t j = 0; j < row.size(); j++)
            out << row[j] << (j + 1 == row.size() ? "" : " ");
        out << "\n";
    }
}

vector<int> genWord(int K)
{
    uniform_int_distribution<int> bit(0, 1);
    vector<int> word(K);
    for (int i = 0; i < K; i++)
        word[i] = bit(rng);
    return word;
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

vector<double> gaussianNoise(const vector<double>& signal, double sigma, int mode = 1, int errors = 1, const vector<int>& errorPositions = vector<int>()) // mode 1 - normal noise, mode 2 - noisy error bit
{
    normal_distribution<double> normal0(0.0, sigma);
    normal_distribution<double> normal2(2.0, sigma);
    uniform_real_distribution<double> uniform01(0.0, 1.0);

    vector<double> noisy(signal.size());
    int error;

    while (error != errors)
    {
        error = 0;

        if (!errorPositions.empty())
        {
            uniform_int_distribution<size_t> pick(0, errorPositions.size() - 1);
            int pos = errorPositions[pick(rng)];

            for (size_t i = 0; i < signal.size(); i++)
            {
                double bit = signal[i];
                double noise = normal0(rng);
                double noisyBit = (mode == 1) ? bit + noise : bit;

                if ((int)i == pos)
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
        }
        else
        {
            int flag = errors;
            double p = 3.0 * errors / signal.size();
            int state = 1;

            for (size_t i = 0; i < signal.size(); i++)
            {
                double bit = signal[i];
                double randNum = uniform01(rng);
                if (state == 1 && randNum < p)
                    state = 2;

                double noise = normal0(rng);
                double noisyBit = (mode == 1) ? bit + noise : bit;

                if (state == 2)
                {
                    noise = normal2(rng);
                    noisyBit = bit + noise;
                    if (bit * noisyBit > 0)
                        noisyBit = bit - noise;
                    flag--;
                    state = (flag != 0) ? 1 : 3;
                }
                if (bit * noisyBit < 0)
                    error++;
                noisy[i] = noisyBit;
            }
        }
    }

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
    int datasetSize = 36000, int mode = 1, int gaussMode = 1, int errors = 1,
    const string& trainFile = "trainset.txt", const string& testFile = "testset.txt",
    const vector<int>& errorPositionsTrain = vector<int>(), const vector<int>& errorPositionsTest = vector<int>()) 
{
    int trainSize = (int)(datasetSize * 0.7);
    int testSize = (int)(datasetSize * 0.3);

    writeSet(G, sigma1Train, sigma2Train, trainSize, mode, gaussMode, errors, trainFile, errorPositionsTrain);
    writeSet(G, sigma1Test, sigma2Test, testSize, mode, gaussMode, errors, testFile, errorPositionsTest);
}

int main()
{
    int iterations = 50;
    float normFactor = 0.75f;
    
    Matrix H = parseAlist("matrix.alist");
    saveMatrix(H, "matrix.txt");

    vector<uint32_t> infoPos;
    Matrix G = get_G(H, infoPos);

    int K = G.size();
    int N = H[0].size();

    Sparse Hv = toSparse(H).turn(aff3ct::tools::Matrix::Way::VERTICAL);
    DecoderNMS decoder(K, N, iterations, Hv, infoPos, NMS(normFactor));

    return 0;
}