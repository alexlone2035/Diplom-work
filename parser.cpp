#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

using Matrix = vector<vector<int>>;

Matrix parseAlist(const string& filename) 
{
    ifstream file(filename);

    int N, M, max_col_deg, max_row_deg;
    file >> N >> M;
    file >> max_col_deg >> max_row_deg;

    vector<int> col_weights(N);
    for (int i = 0; i < N; ++i) 
        file >> col_weights[i];

    vector<int> row_weights(M);
    for (int i = 0; i < M; ++i) 
        file >> row_weights[i];

    Matrix H(M, vector<int>(N, 0));

    for (int j = 0; j < N; j++) 
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
    Matrix H_alist = parseAlist("matrix.alist");
    saveMatrix(H_alist, "matrix.txt");

    return 0;
}