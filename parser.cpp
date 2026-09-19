#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

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

int main()
{
    Matrix H_alist = parseQC("matrix.qc");
    saveMatrix(H_alist, "matrix.txt");

    return 0;
}