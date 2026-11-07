#include "parser.h"
#include "dataset.h"
#include "experiments.h"

vector<float> toLLR(const vector<double>& signal, double sigma)
{
    vector<float> llr(signal.size());
    for (size_t i = 0; i < signal.size(); i++)
        llr[i] = (float)(2.0 * signal[i] / (sigma * sigma));
    return llr;
}

int checkNoiseless(DecoderNMS& decoder, const Matrix& G, int words)
{
    int K = G.size();
    int bad = 0;

    for (int n = 0; n < words; n++)
    {
        vector<int> word = genWord(K);
        vector<float> llr = toLLR(modulate(encode(word, G)), 1.0);
        if (decodeNMS(decoder, llr, K) != word)
            bad++;
    }
    return bad;
}

void runCurve(DecoderNMS& decoder, const Matrix& G, int N, double ebn0Start, double ebn0End, double ebn0Step, int minFrameErrors, int maxFrames, const string& filename)
{
    int K = G.size();
    int points = (int)round((ebn0End - ebn0Start) / ebn0Step) + 1;

    ofstream out(filename);
    out.precision(8);

    for (int p = 0; p < points; p++)
    {
        double ebn0 = ebn0Start + p * ebn0Step;
        double sigma = sigmaFromEbN0(ebn0, K, N);
        long long frames = 0;
        long long frameErrors = 0;
        long long bitErrors = 0;

        while (frameErrors < minFrameErrors && frames < maxFrames)
        {
            vector<int> word = genWord(K);
            vector<double> noisy = awgn(modulate(encode(word, G)), sigma);
            vector<int> decoded = decodeNMS(decoder, toLLR(noisy, sigma), K);

            int errs = 0;
            for (int i = 0; i < K; i++)
            {
                if (decoded[i] != word[i])
                    errs++;
            }
            bitErrors += errs;
            if (errs > 0)
                frameErrors++;
            frames++;
        }
        out << ebn0 << ' ' << sigma << ' ' << (double)bitErrors / (frames * K) << ' ' << (double)frameErrors / frames << ' ' << frames << '\n';
    }
}

vector<int> checksOf(const Matrix& H, int v)
{
    vector<int> checks;
    for (size_t r = 0; r < H.size(); r++)
    {
        if (H[r][v])
            checks.push_back(r);
    }
    return checks;
}

vector<int> varsOf(const Matrix& H, int c)
{
    vector<int> vars;
    for (size_t v = 0; v < H[c].size(); v++)
    {
        if (H[c][v])
            vars.push_back(v);
    }
    return vars;
}

int pairClass(const Matrix& H, int i, int j)
{
    vector<int> checksI = checksOf(H, i);
    vector<int> checksJ = checksOf(H, j);

    vector<int> common;
    for (int a : checksI)
    {
        if (H[a][j])
            common.push_back(a);
    }

    if (common.size() >= 2)
        return SHORT_CYCLE;

    if (common.size() == 1)
    {
        for (int a : checksI)
        {
            if (a == common[0])
                continue;
            for (int v : varsOf(H, a))
            {
                if (v == i || v == j)
                    continue;
                for (int b : checksJ)
                {
                    if (b != common[0] && H[b][v])
                        return SHORT_CYCLE;
                }
            }
        }
        return SHARED_CHECK;
    }

    for (int a : checksI)
    {
        for (int v : varsOf(H, a))
        {
            if (v == i)
                continue;
            for (int b : checksJ)
            {
                if (H[b][v])
                    return SHORT_PATH;
            }
        }
    }
    return DISTANT;
}

