#include "parser.h"

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
