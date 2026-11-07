#pragma once

const int SHARED_CHECK = 0;
const int SHORT_PATH = 1;
const int DISTANT = 2;
const int SHORT_CYCLE = 3;

vector<float> toLLR(const vector<double>& signal, double sigma);
int checkNoiseless(DecoderNMS& decoder, const Matrix& G, int words);
void runCurve(DecoderNMS& decoder, const Matrix& G, int N, double ebn0Start, double ebn0End, double ebn0Step, int minFrameErrors, int maxFrames, const string& filename);
vector<int> checksOf(const Matrix& H, int v);
vector<int> varsOf(const Matrix& H, int c);
int pairClass(const Matrix& H, int i, int j);