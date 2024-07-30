#include <stdio.h>
#include<vector>
#include<string.h>
#include<iostream>
#include<algorithm>
#include<fstream>
#include <cmath>
#include <sys/resource.h>
#include <chrono>
#define SIZE_T int

long getMemoryUsage() 
{
  struct rusage usage;
  if(0 == getrusage(RUSAGE_SELF, &usage))
    return usage.ru_maxrss; // bytes
  else
    return 0;
}

using namespace std;

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

SIZE_T getMem(){
    #ifdef _WIN32
        PROCESS_MEMORY_COUNTERS pmc;
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
            return(pmc.WorkingSetSize); //Value in Bytes!
        else
            return 0;
    #endif
    #ifdef linux
        FILE* file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != NULL){
            if (strncmp(line, "VmSize:", 7) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    #endif
    #if __APPLE__
        return getMemoryUsage() / 1000;
    #endif
    
}

typedef vector<vector<float> > matrix;

class HMM {
    public:
        int nS;
        matrix A;
        matrix B;
        vector<float> pi;
        int max_mem_used;

        HMM() {

        }

        HMM(matrix A, matrix B, vector<float> pi) {
            this->A = A;
            this->B = B;
            this->pi = pi;
            this->nS = A.size();
            this->max_mem_used = 0;
        }

        ~HMM(){
            for (int i = 0; i < A.size(); ++i) {
                A[i].clear();
                A[i].shrink_to_fit();
            }
            for (int i = 0; i < B.size(); ++i) {
                B[i].clear();
                B[i].shrink_to_fit();
            }
            A.clear();
            A.shrink_to_fit();
            B.clear();
            B.shrink_to_fit();
            pi.clear();
            pi.shrink_to_fit();
        }
};

class HMM_Island : public HMM{
    public:
        string method;
        int c;

        HMM_Island()
        {

        }

        HMM_Island(matrix A, matrix B, vector<float> pi, string method) : HMM(A, B, pi){
            this->method = method;
            if (method != "sqrt") 
                this->c = stoi(method);
            // cout << this->c;
        }

        int get_no_checkpoints(int l) {
            if (method == "sqrt") {
                return int(ceil(sqrt(l)));
            }
            return c;
        }

        void forward_backward(vector<int>& obs_seq) {
            vector<float> F0(nS);
            for (int i = 0; i < F0.size(); ++i) 
                F0[i] = pi[i] * B[obs_seq[0]][i];
            vector<float> Bn(nS);
            for (int i = 0; i < Bn.size(); i++)
                Bn[i] = 1.0f;
            _forward_backward(F0, Bn, obs_seq, 0, obs_seq.size()-1);
        }

        vector<float>* normalize(vector<float>& F, vector<float>& B) {
            float sum = 0.0f;
            for (int i=0; i < nS; ++i) 
                sum += F[i] * B[i];
            vector<float> res(nS);
            for (int i=0; i < nS; ++i)
                res[i] = F[i] * B[i] / sum;
            return new vector<float>(res);
        }
        
        void _forward_backward(vector<float>& F0, vector<float>& Bn, vector<int>& obs_seq, int low, int high) {
            int T = high - low + 1;
            // cout << T << endl;
            if (T == 1) {
                //compute and output P
                vector<float>* p = normalize(F0, Bn);
                // for (int i = 0; i < nS; ++i)
                    // cout << "P(x" << i << "|t=" << high << ")=" << p->at(i) << setprecision(5) << endl;
                delete p;
                max_mem_used = max(max_mem_used, getMem());
                return;
            }
            vector<vector<float> > Forw;
            vector<vector<float> > Back;
            
            Forw.push_back(F0);
            vector<float> prev_prob(F0);

            vector<int> idx;
            idx.push_back(low);
            vector<int> idx_2;

            int k = get_no_checkpoints(T);
            int chunk_size = int(ceil(1.0f * T / k));
            // cout << chunk_size <<endl;

            vector<float> alpha(this->nS);
            for (int t = low+1; t <= high; ++t) {
                for (int i = 0; i < nS; ++i) {
                    float sum = 0.0;
                    for (int j = 0; j < nS; ++j) 
                        sum += prev_prob[j] * A[j][i];
                    alpha[i] = sum * B[obs_seq[t]][i];
                }
                prev_prob = alpha;
                if (chunk_size == 0 || T==2) {
                    Forw.push_back(alpha);
                    idx.push_back(t);
                }
                else
                if ((t-low+1) % chunk_size == 0 && t != high) {
                    Forw.push_back(alpha);
                    idx.push_back(t);
                }
            }
            vector<float> beta(this->nS);
            prev_prob = Bn;
            for (int t = high-1; t >= low; --t) {
                for (int i = 0; i < nS; ++i) {
                    float sum = 0.0;
                    for (int j = 0; j < nS; ++j) 
                        sum += prev_prob[j] * A[i][j] * B[obs_seq[t+1]][j];
                    beta[i] = sum;
                }
                prev_prob = beta;
                if (chunk_size == 0 || T==2) {
                    Back.push_back(beta);
                    idx_2.push_back(t);
                }
                else
                if ((t-low+1) % chunk_size == 0 && t != low) {
                    Back.push_back(beta);
                    idx_2.push_back(t);
                }
            }
            reverse(Back.begin(), Back.end());
            reverse(idx_2.begin(), idx_2.end());
            Back.push_back(Bn);
            idx_2.push_back(high);
            //compute and output a_k * b_k
            for (int i = 0; i < Forw.size(); ++i)
            {
                // cout << idx[i] << " " << idx_2[i] << endl;
                _forward_backward(Forw[i], Back[i], obs_seq, idx[i], idx_2[i]);
            }
        }
};

// float island()

matrix generate_transition(int nS) {
    matrix A(nS);
    for (int i = 0; i < nS; ++i) {
        float norm = 0.0f;
        for (int j = 0; j < nS; ++j) {
            float v = 1.0f * (rand() % 2000 + 1);
            A[i].push_back(v);
            norm += v;
        }
        for (int j = 0; j < nS; ++j) {
            A[i][j] = A[i][j] / norm;
        }
    }
    return A;
}
matrix generate_emission(int nS, int nO) {
    matrix B(nO);
    for (int o = 0; o < nO; ++o) {
        for (int j = 0; j < nS; ++j) {
            B[o].push_back(1.0f * (rand() % 2000 + 1));
        }
    }
    for (int j = 0; j < nS; ++j) {
        float norm = 0.0f;
        for (int o = 0; o < nO; ++o) {
            norm += B[o][j];
        }
        for (int o = 0; o < nO; ++o) {
            B[o][j] /= norm;
        }
    }
    return B;
}

vector<float> generate_pi(int nS) {
    vector<float> pi(nS);
    float norm = 0.0f;
    for (int i = 0; i < nS; ++i) {
        float v = 1.0f * (rand() % 2000 + 1);
        pi[i] = v;
        norm += v;
    }
    for (int i = 0; i < nS; ++i) {
        pi[i] /= norm;
    }
    return pi;
}


int main(int argc, char* argv[]) {
    int T = stoi(argv[1]);
    int nS = stoi(argv[4]);
    int nO = stoi(argv[5]);
    matrix A = generate_transition(nS);
    matrix B = generate_emission(nS, nO);
    vector<float> pi = generate_pi(nS);

    
    vector<float> g;
    vector<int> y;
    for (int i = 0; i < T; ++i) 
    {
        y.push_back(rand() % nO);
    }

    int max_mem_used = 0;
    int opt = stoi(argv[2]);
    int trial = stoi(argv[3]);
    string method = "";
    HMM_Island normal_model = HMM_Island();
    auto start = std::chrono::high_resolution_clock::now();
    switch (opt)
    {
    case 0:
        method = "linear";
        normal_model = HMM_Island(A, B, pi, to_string(T));
        normal_model.forward_backward(y);
        max_mem_used = normal_model.max_mem_used;
        break;
    case 1:
        method = "island-log2";
        normal_model = HMM_Island(A, B, pi, "2");
        normal_model.forward_backward(y);
        max_mem_used = normal_model.max_mem_used;
        break;
    case 2:
        method = "island-sqrtN";
        normal_model = HMM_Island(A, B, pi, to_string(int(sqrt(T))));
        normal_model.forward_backward(y);
        max_mem_used = normal_model.max_mem_used;
        break;
    case 3:
        method = "sqrt";
        normal_model = HMM_Island(A, B, pi, "sqrt");
        normal_model.forward_backward(y);
        max_mem_used = normal_model.max_mem_used;
        break;

    default:
        break;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    ofstream fout;
    fout.open(argv[6], ios_base::app);
    fout << method << "," << to_string(T) << "," << to_string(max_mem_used) << "," << to_string(trial) << "," << to_string(duration.count()) <<  endl;
    fout.close();

    return 0;
}