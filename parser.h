#pragma once
#include <vector>
#include <string>
#include<aff3ct.hpp>
#include <random>
#include <fstream>

using namespace std;
using Sparse = aff3ct::tools::Sparse_matrix;
using Handler = aff3ct::tools::LDPC_matrix_handler;
using NMS = aff3ct::tools::Update_rule_NMS<float>;
using DecoderNMS = aff3ct::module::Decoder_LDPC_BP_flooding<int, float, NMS>;

using Matrix = vector<vector<int>>;

Matrix parseQC(const std::string& filename);
Matrix parseAlist(const string& filename);
Sparse toSparse(const Matrix& H);
Matrix toDense(const Sparse& S);
Matrix get_G(const Matrix& H, vector<uint32_t>& infoPos);
vector<int> decodeNMS(DecoderNMS& decoder, const vector<float>& code, int K);
void saveMatrix(const Matrix& H, const string& filename);