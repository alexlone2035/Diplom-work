#include "parser.h"
#include "dataset.h"
#include "experiments.h"

int main()
{
    int iterations = 50;
    float normFactor = 0.75f;

    Matrix H = parseAlist("matrix.alist");
    saveMatrix(H, "matrix.txt");

    vector<uint32_t> infoPos;
    Matrix G = get_G(H, infoPos);
    saveMatrix(G, "generator.txt");

    int K = G.size();
    int N = H[0].size();
    int M = H.size();

    ofstream info("code_info.txt");
    info << "N = " << N << "\nM = " << M << "\nK = " << K << "\nrank = " << N - K << "\n";
    info.close();

    Sparse Hv = toSparse(H).turn(aff3ct::tools::Matrix::Way::VERTICAL);
    DecoderNMS decoder(K, N, iterations, Hv, infoPos, NMS(normFactor));

    int words = 1000;
    int bad = checkNoiseless(decoder, G, words);
    ofstream check("nms_check.txt");
    check << bad << " of " << words << " noiseless words decoded incorrectly\n";
    check.close();

    runCurve(decoder, G, N, 0.0, 6.0, 0.5, 100, 100000, "nms_curve.txt");

    getDataset(G, 0.1, 0.5, 0.3, 0.3, 1000, 3, 1, 1);

    return 0;
}