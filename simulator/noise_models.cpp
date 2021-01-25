# include <iostream>
# include <random>
# include <vector>
# include <cstdlib>
# include <algorithm>
# include <assert.h>
# include <functional>

# include "functions.hpp"



void flip_if_exists(std::vector<int>& Z_vec, int i) {
    if (i >= 0 && i < Z_vec.size()) {
        Z_vec[i] *= -1;
    } else {
        // std::cout << "got a bad i="<<i<<std::endl;
    }
}

std::pair<int,int> unif_random_pauli_len2(std::function<double()> rand) {
    // returns a random vector of 2 paulis
    // other than II
    // there are 15 possibilities
    // Samples a random int in range [1, 2, ..., 14, 15]
    // then returns NUM % 4, int(NUM / 4)
    int num = int(15*rand());
    int p1 = num % 4;
    int p2 = int(num/4);
    assert(p1+p2*4 == num);
    return std::make_pair(p1, p2);
}

void SingleSite_err(const double p, const int u, std::vector<int>& Z_vec, std::vector<int>& X_vec, std::function<double()> rand) {
    if (u>=0 && u<Z_vec.size()) {
        double rp = rand();
        // if (rp <= 2*p/3) {
        //     // X error
        //     X_vec[u] *= -1;
        // }
        // if (rp > p/3 && rp <= p) {
        //     // Z error
        //     Z_vec[u] *= -1;
        // }

        // for measurement and initialization, either Z or I
        if (rp < p) {
            Z_vec[u] *= -1;
        }
    }
}

void CZ_err(const double p, const int u, const int v, std::vector<int>& Z_vec, std::vector<int>& X_vec, std::function<double()> rand) {
    if (u >= 0 && u < Z_vec.size() && v >= 0 && v < Z_vec.size()) {
        // propagate X errors past this CZ gate
        Z_vec[u] *= X_vec[v];
        Z_vec[v] *= X_vec[u];

        // new errors introduced due to the present CZ gate
        // u

        double rp = rand();
        if (rp <= p) {
            // some pauli other than II
            auto err =  unif_random_pauli_len2(rand);
            if (err.first % 2) {
                // X_1
                flip_if_exists(X_vec, u);
            }
            if (int(err.first/2)) {
                // Z_1
                flip_if_exists(Z_vec, u);
            }
            if (err.second % 2) {
                // X_1
                flip_if_exists(X_vec, v);
            }
            if (int(err.second/2)) {
                // Z_1
                flip_if_exists(Z_vec, v);
            }
        }
        // double rp11 = rand();
        // if (rp11 <= 2*p/3) {
        //     // X error on u qubit
        //     flip_if_exists(X_vec, u);
        // }
        // if (p/3 < rp11 && rp11 <= p) {
        //     // Z error on u
        //     flip_if_exists(Z_vec, u);
        // }
        //
        // // v
        // double rp12 = rand();
        // if (rp12 <= 2*p/3) {
        //     // X error on neighbor edge qubit
        //     flip_if_exists(X_vec, v);
        // }
        // if (p/3 < rp12 && rp12 <= p) {
        //     // Z error on neighbor edge qubit
        //     flip_if_exists(Z_vec, v);
        // }

        // note that all Z's commute straight to the end
    } else {
        // std::cerr << "invalid edge: (u,v)=("<<u<<","<<v<<")"<<std::endl;
    }
}

std::vector<std::vector<int>> get_sequences(int L, int M) {
    // directions
    int top = -L;
    int bot = +L;
    int right = +1;
    int left = -1;
    int back = +L*M;
    int front = -L*M;


    #if 0
    // Noah's arbitrary sequence
    std::vector<int> top_seq = {back, left, right, front};
    std::vector<int> bot_seq = top_seq;
    std::vector<int> front_seq = {right, top, bot, left};
    std::vector<int> back_seq = front_seq;
    std::vector<int> right_seq = {top, front, back, bot};
    std::vector<int> left_seq = right_seq;
    #endif

    // #if 0
    // RHG sequence
    std::vector<int> top_seq = {front, back, right, left};
    std::vector<int> bot_seq = top_seq;
    std::vector<int> front_seq = {right, left, top, bot};
    std::vector<int> back_seq = front_seq;
    std::vector<int> right_seq = {bot, top, back, front};
    std::vector<int> left_seq = right_seq;
    // #endif


    return {
        front_seq,
        back_seq,
        left_seq,
        right_seq,
        top_seq,
        bot_seq
    };
}


std::vector<int> noise_stephens_standard_gateerror(int L, int M, int N, double p, std::function<double()> rand) {
    std::vector<int> Z_vec(L*M*N, 1);
    std::vector<int> X_vec(L*M*N, 1);

    std::vector<int> Vxy = get_Vxy(L, M, N); // +/- y (+-1) +/- x
    std::vector<int> Vyz = get_Vyz(L, M, N);
    std::vector<int> Vzx = get_Vzx(L, M, N);


    int l = (L + 1) / 2;
    int m = (M + 1) / 2;
    int n = (N - 1) / 2;


    std::vector<std::vector<int>> front_back_left_right_top_bot_seq = get_sequences(L, M);
    for (int step=0; step<6; step++) {
        // keep track of when we visited a face, i.e. applied the noise to all gates
        std::vector<bool> visited(L*M*N, false);
        for (int c = 0; c < l*m*n; c++) { // iterate over all cubes; the depth-4 CZ gate procedure is defined at cube-level
            int x = (c % (l*m)) % l;
            int y = (c % (l*m)) / l;
            int z = c / (l*m);


            int front_face = 2*x + 2*y*L + 2*z*L*M; // qubit number of qubit on the front face (smallest z value) of cube c

            if (!visited[front_face] ) {
                visited[front_face] = true;

                if (step == 0 || step == 5) {
                    SingleSite_err(p, front_face, Z_vec, X_vec, rand);
                } else {
                    int dir = front_back_left_right_top_bot_seq.at(0).at(step-1);
                    CZ_err(p, front_face, front_face + dir, Z_vec, X_vec, rand);
                }
            }


            // back face
            int back_face = front_face + 2*L*M;
            if (!visited[back_face]) {
                visited[back_face] = true;
                if (step == 0 || step == 5) {
                    SingleSite_err(p, back_face, Z_vec, X_vec, rand);
                } else {
                    int dir = front_back_left_right_top_bot_seq.at(1).at(step-1);
                    CZ_err(p, back_face, back_face + dir, Z_vec, X_vec, rand);
                }
            }

            // left face
            if (x != 0) { // cubes at the left (x = 0) boundary are partial and don't have left faces
                int left_face = front_face - 1 + L*M;
                if(!visited[left_face]) {
                    visited[left_face] = true;
                    if (step == 0 || step == 5) {
                        SingleSite_err(p, left_face, Z_vec, X_vec, rand);
                    } else {
                        int dir = front_back_left_right_top_bot_seq.at(2).at(step-1);
                        CZ_err(p, left_face, left_face + dir, Z_vec, X_vec, rand);
                    }
                }
            }
            // right face
            if (x != l - 1) { // cubes at the right boundary don't have right faces
                int right_face = front_face + 1 + L*M;
                // TODO should compute X, Y, Z
                // if X in {0, L-1} you are at an X edge i.e. a YZ plane
                // if Y == {0, LM-1} you are at a Y edge i.e. a XZ plane
                // if Z == {0, LMN-1} you are at a Z edge i.e. a XY plane
                if(!visited[right_face]) {
                    visited[right_face] = true;
                    if (step == 0 || step == 5) {
                        SingleSite_err(p, right_face, Z_vec, X_vec, rand);
                    } else {
                        int dir = front_back_left_right_top_bot_seq.at(3).at(step-1);
                        CZ_err(p, right_face, right_face + dir, Z_vec, X_vec, rand);
                    }
                }
            }
            // top face
            if (y != 0) {
                int top_face = front_face - L + L*M;
                if(!visited[top_face]) {
                    visited[top_face] = true;
                    if (step == 0 || step == 5) {
                        SingleSite_err(p, top_face, Z_vec, X_vec, rand);
                    } else {
                        int dir = front_back_left_right_top_bot_seq.at(4).at(step-1);
                        CZ_err(p, top_face, top_face + dir, Z_vec, X_vec, rand);
                    }
                }
            }
            // bottom face
            if (y != m - 1) {
                int bottom_face = front_face + L + L*M;
                if(!visited[bottom_face]) {
                    visited[bottom_face] = true;
                    if (step == 0 || step == 5) {
                        SingleSite_err(p, bottom_face, Z_vec, X_vec, rand);
                    } else {
                        int dir = front_back_left_right_top_bot_seq.at(5).at(step-1);
                        CZ_err(p, bottom_face, bottom_face + dir, Z_vec, X_vec, rand);
                    }
                }
            }
        }
    }

    return Z_vec;
}


std::vector<int> noise_RHG_error_model1_local_depolarizing(int L, int M, int N, double p, std::function<double()> rand) {
    std::vector<int> Z_vec(L*M*N, 1);

    std::vector<int> Vxy = get_Vxy(L, M, N); // +/- y (+-1) +/- x
    std::vector<int> Vyz = get_Vyz(L, M, N);
    std::vector<int> Vzx = get_Vzx(L, M, N);

    // maybe off by a factor of 2/3 or 3/2
    // for (int v = 0; v < Vxy.size(); v++) {
    for (int i=0; i<L*M*N; i++) {
        // int i = Vxy[v]; // qubit number
        int X = (i % (L*M)) % L;
        int Y = (i % (L*M)) / L;
        int Z = i / (L*M);

        // see https://iopscience.iop.org/article/10.1088/1367-2630/9/6/199/meta
        // section 7.1
        double rp1 = rand();
        if (rp1 <= 2*p/3) {
            // Z error
            // flip this site
            Z_vec[i] *= -1;
        }

    }
    return Z_vec;
}

std::vector < int > noise_single_emitter(int L, int M, int N, double p, std::function<double()> rand) {
    std::vector < int > Z_vec(L * M * N, 1);

    std::vector < int > Vxy = get_Vxy(L, M, N);
    std::vector < int > Vyz = get_Vyz(L, M, N);
    std::vector < int > Vzx = get_Vzx(L, M, N);

    for (int v = 0; v < Vxy.size(); v++) {
        int i = Vxy[v]; // qubit number
        int X = (i % (L * M)) % L;
        int Y = (i % (L * M)) / L;
        int Z = i / (L * M);

        // after initialisation
        double rp0 = rand();
        if (rp0 < 2 * p / 3) {
            if ((i + L < L * M * N) && (Y != M - 1)) { // no error if i + L isn't actually a neighbour of i
                Z_vec[i + L] *= -1;
            }
        }

        // after Z_{Q,i-L} (which is applied only for i-L >= 0 AND (i - L, i) \in E, i.e., Y > 0)
        if ((i >= L) && (Y != 0)) {
            double rp1 = rand();

            if (1 - p + 3 * p / 15 < rp1 && rp1 < 1 - p + 11 * p / 15) { // X_Q
                Z_vec[i - L] *= -1;
                if (X != 0) {
                    Z_vec[i - 1] *= -1;
                }
            }

            if (1 - p + 7 * p / 15 < rp1 && rp1 < 1) { // Z_Q
                Z_vec[i] *= -1;
            }
            // *note: Y_Q if 1-p + 7*p/15 < rp1 && rp1 < 1-p + 11*p/15

            if ((1 - p < rp1 && rp1 < 1 - p + 2 * p / 15) ||
                (1 - p + 4 * p / 15 < rp1 && rp1 < 1 - p + 6 * p / 15) ||
                (1 - p + 8 * p / 15 < rp1 && rp1 < 1 - p + 10 * p / 15) ||
                (1 - p + 12 * p / 15 < rp1 && rp1 < 1 - p + 14 * p / 15)) { // X_{i-L}
                if (i - L + L * M < L * M * N) {
                    Z_vec[i - L + L * M] *= -1;
                }
            }

            if ((1 - p + p / 15 < rp1 && rp1 < 1 - p + 3 * p / 15) ||
                (1 - p + 5 * p / 15 < rp1 && rp1 < 1 - p + 7 * p / 15) ||
                (1 - p + 9 * p / 15 < rp1 && rp1 < 1 - p + 11 * p / 15) ||
                (1 - p + 13 * p / 15 < rp1 && rp1 < 1)) { // Z_{i-L}
                Z_vec[i - L] *= -1;
            }
            // *again, this includes Y_{i-L} errors
        }

        // after X_{Q,i}
        double rp2 = rand();

        if (1 - p + 3 * p / 15 < rp2 && rp2 < 11 * p / 15) { // X_Q
            if ((i + 1 < L * M * N) && (X != L - 1)) {
                Z_vec[i + 1] *= -1;
            }
        }

        if (1 - p + 7 * p / 15 < rp2 && rp2 < 1) { // Z_Q
            Z_vec[i] *= -1;
        }

        if ((1 - p < rp2 && rp2 < 1 - p + 2 * p / 15) ||
            (1 - p + 4 * p / 15 < rp2 && rp2 < 1 - p + 6 * p / 15) ||
            (1 - p + 8 * p / 15 < rp2 && rp2 < 1 - p + 10 * p / 15) ||
            (1 - p + 12 * p / 15 < rp2 && rp2 < 1 - p + 14 * p / 15)) { // X_i
            if ((i + L < L * M * N) && (Y != M - 1)) {
                Z_vec[i + L] *= -1;
            }
        }

        if ((1 - p + p / 15 < rp2 && rp2 < 1 - p + 3 * p / 15) ||
            (1 - p + 5 * p / 15 < rp2 && rp2 < 1 - p + 7 * p / 15) ||
            (1 - p + 9 * p / 15 < rp2 && rp2 < 1 - p + 11 * p / 15) ||
            (1 - p + 13 * p / 15 < rp2 && rp2 < 1)) { // Z_i
            Z_vec[i] *= -1;
        }

        // after H_Q
        double rp3 = rand();

        if (rp3 < 2 * p / 3) {
            Z_vec[i] *= -1;
        }

        if (p / 3 < rp3 && rp3 < p) {
            if ((i + 1 < L * M * N) && (X != L - 1)) {
                Z_vec[i + 1] *= -1;
            }
        }

    }

    for (int v = 0; v < Vyz.size(); v++) {
        int i = Vyz[v]; // qubit number
        int X = (i % (L * M)) % L;
        int Y = (i % (L * M)) / L;
        int Z = i / (L * M);

        // after initialisation
        double rp0 = rand();
        if (rp0 < 2 * p / 3) {
            if ((i + L < L * M * N) && (Y != M - 1)) {
                Z_vec[i + L] *= -1;
            }
            if (i + L * M < L * M * N) {
                Z_vec[i + L * M] *= -1;
            }
        }

        // after Z_{Q,i-LM} (which is only applied for i-LM >= 0 AND (i - LM, i) \in E)
        if (i >= L * M) { // if i >= L*M, then (i - LM, i) \in E (for i \in V^{yz})
            double rp1 = rand();

            if (1 - p + 3 * p / 15 < rp1 && rp1 < 1 - p + 11 * p / 15) { // X_Q
                Z_vec[i - L * M] *= -1;
            }

            if (1 - p + 7 * p / 15 < rp1 && rp1 < 1) { // Z_Q
                Z_vec[i] *= -1;
            }

            // if an X_{i-LM} error occurs after Z_{Q,i-LM}, the effective error is X_{i-LM} -- doesn't affect syndrome

            if ((1 - p + p / 15 < rp1 && rp1 < 1 - p + 3 * p / 15) ||
                (1 - p + 5 * p / 15 < rp1 && rp1 < 1 - p + 7 * p / 15) ||
                (1 - p + 9 * p / 15 < rp1 && rp1 < 1 - p + 11 * p / 15) ||
                (1 - p + 13 * p / 15 < rp1 && rp1 < 1)) { // Z_{i-LM}
                Z_vec[i - L * M] *= -1;
            }
        }

        // after Z_{Q,i-L} (which is only applied for i-L >= 0 AND (i - L, i) \in E)
        if ((i >= L) && (Y != 0)) {
            double rp2 = rand();

            if (1 - p + 3 * p / 15 < rp2 && rp2 < 1 - p + 11 * p / 15) { // X_Q
                if (i - L * M >= 0) {
                    Z_vec[i - L * M] *= -1;
                }
                Z_vec[i - L] *= -1;
            }

            if (1 - p + 7 * p / 15 < rp2 && rp2 < 1) { // Z_Q
                Z_vec[i] *= -1;
            }

            // if an X_{i-L} error occurs after Z_{Q,i-L}, the effective error is X_{i-L} -- doesn't affect syndrome

            if ((1 - p + p / 15 < rp2 && rp2 < 1 - p + 3 * p / 15) ||
                (1 - p + 5 * p / 15 < rp2 && rp2 < 1 - p + 7 * p / 15) ||
                (1 - p + 9 * p / 15 < rp2 && rp2 < 1 - p + 11 * p / 15) ||
                (1 - p + 13 * p / 15 < rp2 && rp2 < 1)) { // Z_{i-L}
                Z_vec[i - L] *= -1;
            }
        }

        // after X_{Q,i}
        double rp3 = rand();

        // if an X_Q error occurs after X_{Q,i} for i \in Vyz, no effective error *FOR THE MEASUREMENT PROTOCOL*

        if (1 - p + 7 * p / 15 < rp3 && rp3 < 1) { // Z_Q
            Z_vec[i] *= -1;
        }

        if ((1 - p < rp3 && rp3 < 1 - p + 2 * p / 15) ||
            (1 - p + 4 * p / 15 < rp3 && rp3 < 1 - p + 6 * p / 15) ||
            (1 - p + 8 * p / 15 < rp3 && rp3 < 1 - p + 10 * p / 15) ||
            (1 - p + 12 * p / 15 < rp3 && rp3 < 1 - p + 14 * p / 15)) { // X_i
            if ((i + L < L * M * N) && (Y != M - 1)) {
                Z_vec[i + L] *= -1;
            }
            if (i + L * M < L * M * N) {
                Z_vec[i + L * M] *= -1;
            }
        }

        if ((1 - p + p / 15 < rp3 && rp3 < 1 - p + 3 * p / 15) ||
            (1 - p + 5 * p / 15 < rp3 && rp3 < 1 - p + 7 * p / 15) ||
            (1 - p + 9 * p / 15 < rp3 && rp3 < 1 - p + 11 * p / 15) ||
            (1 - p + 13 * p / 15 < rp3 && rp3 < 1)) { // Z_i
            Z_vec[i] *= -1;
        }

        // after H_Q
        double rp4 = rand();
        if (rp4 < 2 * p / 3) {
            Z_vec[i] *= -1;
        }

        // before measurement of Q
        double rp5 = rand();
        if (rp5 < 2 * p / 3) {
            Z_vec[i] *= -1;
        }

        // after re-initialisation of Q
        double rp6 = rand();
        if (rp6 < 2 * p / 3) {
            if (i + 2 < L * M * N) {
                Z_vec[i + 2] *= -1;
            }
        }
    }

    for (int v = 0; v < Vzx.size(); v++) {
        int i = Vzx[v]; // qubit number
        int X = (i % (L * M)) % L;
        int Y = (i % (L * M)) / L;
        int Z = i / (L * M);

        // after initialisation
        double rp0 = rand();
        if (rp0 < 2 * p / 3) {
            if (i + L * M < L * M * N) {
                Z_vec[i + L * M] *= -1;
            }
        }

        // after Z_{Q,i-LM} (which is only applied for i-LM >= 0)
        if (i >= L * M) {
            double rp1 = rand();

            if (1 - p + 3 * p / 15 < rp1 && rp1 < 1 - p + 11 * p / 15) { // X_Q
                Z_vec[i - L * M] *= -1;
                if (X != 0) {
                    Z_vec[i - 1] *= -1;
                }
            }

            if (1 - p + 7 * p / 15 < rp1 && rp1 < 1) { // Z_Q
                Z_vec[i] *= -1;
            }

            // if an X_{i-L} error occurs after Z_{Q,i-L}, the effective error is X_{i-L} -- doesn't affect syndrome

            if ((1 - p + p / 15 < rp1 && rp1 < 1 - p + 3 * p / 15) ||
                (1 - p + 5 * p / 15 < rp1 && rp1 < 1 - p + 7 * p / 15) ||
                (1 - p + 9 * p / 15 < rp1 && rp1 < 1 - p + 11 * p / 15) ||
                (1 - p + 13 * p / 15 < rp1 && rp1 < 1)) { // Z_{i-LM}
                Z_vec[i - L * M] *= -1;
            }
        }

        // after X_{Q,i}
        double rp2 = rand();

        if (1 - p + 3 * p / 15 < rp2 && rp2 < 11 * p / 15) { // X_Q
            if ((i + 1 < L * M * N) && (X != L - 1)) {
                Z_vec[i + 1] *= -1;
            }
        }

        if (1 - p + 7 * p / 15 < rp2 && rp2 < 1) { // Z_Q
            Z_vec[i] *= -1;
        }

        if ((1 - p < rp2 && rp2 < 1 - p + 2 * p / 15) ||
            (1 - p + 4 * p / 15 < rp2 && rp2 < 1 - p + 6 * p / 15) ||
            (1 - p + 8 * p / 15 < rp2 && rp2 < 1 - p + 10 * p / 15) ||
            (1 - p + 12 * p / 15 < rp2 && rp2 < 1 - p + 14 * p / 15)) { // X_i
            if (i + L * M < L * M * N) {
                Z_vec[i + L * M] *= -1;
            }
        }

        if ((1 - p + p / 15 < rp2 && rp2 < 1 - p + 3 * p / 15) ||
            (1 - p + 5 * p / 15 < rp2 && rp2 < 1 - p + 7 * p / 15) ||
            (1 - p + 9 * p / 15 < rp2 && rp2 < 1 - p + 11 * p / 15) ||
            (1 - p + 13 * p / 15 < rp2 && rp2 < 1)) { // Z_i
            Z_vec[i] *= -1;
        }

        // after H_Q
        double rp3 = rand();

        if (rp3 < 2 * p / 3) {
            Z_vec[i] *= -1;
        }

        if (p / 3 < rp3 && rp3 < p) {
            if ((i + 1 < L * M * N) && (X != L - 1)) {
                Z_vec[i + 1] *= -1;
            }
        }
    }

    // syndrome measurement errors!
    for (int i = 0; i < Z_vec.size(); i++) {
        double rp = rand();
        if (rp < 2 * p / 3) {
            Z_vec[i] *= -1;
        }
    }

    // initialisation of Q at the very beginning; note that the first qubit (qubit 0) is in Vxy
    double rpi = rand();
    if (rpi < 2 * p / 3) {
        Z_vec[0] *= -1;
    }

    // measurement of Q at the very end (to remove it from the cluster)
    double rpf = rand();
    if (rpf < 2 * p / 3) {
        Z_vec[L * M * N - 1] *= -1;
    }

    return Z_vec;
}




std::vector<int> get_Vbcc(int L, int M, int N) { // returns indices in LARGER lattice of bcc lattice qubits
  std::vector<int> Vbcc;
  int LL = L + 2;
  int MM = M + 2;
  int NN = N + 2;

  for (int k = 0; k < NN/2; k++) {
    for (int j = 0; j < MM/2; j++) {
      for (int i = 0; i < LL/2; i++) {
        Vbcc.push_back((2*k + 1)*LL*MM + (2*j + 1)*LL + (2*i + 1));
      }
    }
  }
  for (int k = 0; k < NN/2 - 1; k++) {
    for (int j = 0; j < MM/2 - 1; j++) {
      for (int i = 0; i < LL/2 - 1; i++) {
        Vbcc.push_back((2*k + 2)*LL*MM + (2*j + 2)*LL + (2*i + 2));
      }
    }
  }
  for (int k = 0; k < NN/2; k++) {
    for (int j = 0; j < MM/2 - 1; j++) {
      for (int i = 0; i < LL/2; i++) {
        Vbcc.push_back((2*k + 1)*LL*MM + (2*j + 2)*LL + (2*i + 1));
      }
    }
  }
  for (int k = 0; k < NN/2 - 1; k++) {
    for (int j = 0; j < MM/2; j++) {
      for (int i = 0; i < LL/2 - 1; i++) {
        Vbcc.push_back((2*k + 2)*LL*MM + (2*j + 1)*LL + (2*i + 2));
      }
    }
  }
  for (int k = 0; k < NN/2; k++) {
    for (int j = 0; j < MM/2; j++) {
      for (int i = 0; i < LL/2 - 1; i++) {
        Vbcc.push_back((2*k + 1)*LL*MM + (2*j + 1)*LL + (2*i + 2));
      }
    }
  }
  for (int k = 0; k < NN/2 - 1; k++) {
    for (int j = 0; j < MM/2 - 1; j++) {
      for (int i = 0; i < LL/2; i++) {
        Vbcc.push_back((2*k + 2)*LL*MM + (2*j + 2)*LL + (2*i + 1));
      }
    }
  }

  std::sort(Vbcc.begin(), Vbcc.end());
  return Vbcc;
}


std::vector<int> get_Vextra(int L, int M, int N) {
  std::vector<int> Vextra;
  int LL = L + 2;
  int MM = M + 2;
  int NN = N + 2;

  for (int l = 0; l < LL*MM; l++) { // first layer (Z = 0)
    Vextra.push_back(l);
  }
  for (int l = LL*MM*NN - LL*MM; l < LL*MM*NN; l++) { // last layer (Z = NN - 1)
    Vextra.push_back(l);
  }

  for (int k = 1; k < NN - 1; k++) {  // rest of boundary
    for (int i = 0; i < LL; i++) {
      Vextra.push_back(k*LL*MM + i);
      Vextra.push_back(k*LL*MM + (MM - 1)*LL + i);
    }
    for (int j = 1; j < MM - 1; j++) {
      Vextra.push_back(k*LL*MM + j*LL);
      Vextra.push_back(k*LL*MM + j*LL + (LL - 1));
    }
  }

  for (int k = 0; k < NN/2; k++) {
    for (int j = 0; j < MM/2 - 1; j++) {
      for (int i = 0; i < LL/2 - 1; i++) {
        Vextra.push_back((2*k + 1)*LL*MM + (2*j + 2)*LL + (2*i + 2));
      }
    }
  }
  for (int k = 0; k < NN/2 - 1; k++) {
    for (int j = 0; j < MM/2; j++) {
      for (int i = 0; i < LL/2; i++) {
        Vextra.push_back((2*k + 2)*LL*MM + (2*j + 1)*LL + (2*i + 1));
      }
    }
  }

  std::sort(Vextra.begin(), Vextra.end());
  return Vextra;
}


std::vector<int> big_cubic_noise(int L, int M, int N, double p, std::function<double()> rand) { // here, L, M, N are the dimensions of the BIGGER lattice; this function is exactly the same as the one in cubic.cpp

  std::vector<int> Z_vec(L*M*N, 1);

  std::vector<int> Vbcc = get_Vbcc(L - 2, M - 2, N - 2);
  std::vector<int> Vextra = get_Vextra(L - 2, M - 2, N - 2);

  // initialisation of Q at the very beginning
  double rpi = rand();
  if (rpi < 2*p/3) {
    Z_vec[0] *= -1;
  }

  // measurement of Q at the very end
  double rpf = rand();
  if (rpf < 2*p/3) {
    Z_vec[L*M*N - 1] *= -1;
  }

  for (int i = 0; i < L*M*N; i++) {
    // after initialisation of i
    double rp0 = rand();
    if (rp0 < 2*p/3) {
      if (i + L < L*M*N) {
        Z_vec[i + L] *= -1;
      }
      if (i + L*M < L*M*N) {
        Z_vec[i + L*M] *= -1;
      }
    }


    // after Z_{Q,i-LM} (which is only applied for i-LM >= 0)
    if (i >= L*M) {
      double rp1 = rand();

      if (1-p + 3*p/15 < rp1 && rp1 < 1-p + 11*p/15) {  // X_Q
        Z_vec[i - L*M] *= -1;
        Z_vec[i - 1] *= -1;
      }

      if (1-p + 7*p/15 < rp1 && rp1 < 1) {  // Z_Q
        Z_vec[i] *= -1;
      }

      // if an X_{i-LM} error occurs after Z_{Q,i-LM}, the effective error is X_{i-LM} -- doesn't affect syndrome

      if ((1-p +    p/15 < rp1 && rp1 < 1-p +  3*p/15) ||
          (1-p +  5*p/15 < rp1 && rp1 < 1-p +  7*p/15) ||
          (1-p +  9*p/15 < rp1 && rp1 < 1-p + 11*p/15) ||
          (1-p + 13*p/15 < rp1 && rp1 < 1            )) { // Z_{i-LM}
            Z_vec[i - L*M] *= -1;
      }
    }

    // after Z_{Q,i-L} (which is only applied for i-L >= 0)
    if (i >= L) {
      double rp2 = rand();

      if (1-p + 3*p/15 < rp2 && rp2 < 1-p + 11*p/15) {  // X_Q
        if (i - L*M >= 0) {
          Z_vec[i - L*M] *= -1;
        }
        Z_vec[i - L] *= -1;
        Z_vec[i - 1] *= -1;
      }

      if (1-p + 7*p/15 < rp2 && rp2 < 1) {  // Z_Q
        Z_vec[i] *= -1;
      }

      if ((1-p           < rp2 && rp2 < 1-p +  2*p/15) ||
          (1-p +  4*p/15 < rp2 && rp2 < 1-p +  6*p/15) ||
          (1-p +  8*p/15 < rp2 && rp2 < 1-p + 10*p/15) ||
          (1-p + 12*p/15 < rp2 && rp2 < 1-p + 14*p/15)) { // X_{i-L}
            if (i - L + L*M < L*M*N) {
              Z_vec[i - L + L*M] *= -1;
            }
      }

      if ((1-p +    p/15 < rp2 && rp2 < 1-p +  3*p/15) ||
          (1-p +  5*p/15 < rp2 && rp2 < 1-p +  7*p/15) ||
          (1-p +  9*p/15 < rp2 && rp2 < 1-p + 11*p/15) ||
          (1-p + 13*p/15 < rp2 && rp2 < 1            )) { // Z_{i-L}
            Z_vec[i - L] *= -1;
      }
    }

    // after X_{Q,i}
    double rp3 = rand();

    if (1-p + 3*p/15 < rp3 && rp3 < 11*p/15) {  // X_Q
      if (i + 1 < L*M*N) {
        Z_vec[i + 1] *= -1;
      }
    }

    if (1-p + 7*p/15 < rp3 && rp3 < 1) {  // Z_Q
      Z_vec[i] *= -1;
    }

    if ((1-p           < rp3 && rp3 < 1-p +  2*p/15) ||
        (1-p +  4*p/15 < rp3 && rp3 < 1-p +  6*p/15) ||
        (1-p +  8*p/15 < rp3 && rp3 < 1-p + 10*p/15) ||
        (1-p + 12*p/15 < rp3 && rp3 < 1-p + 14*p/15)) { // X_i
          if (i + L < L*M*N) {
            Z_vec[i + L] *= -1;
          }
          if (i + L*M < L*M*N) {
            Z_vec[i + L*M] *= -1;
          }
    }

    if ((1-p +    p/15 < rp3 && rp3 < 1-p +  3*p/15) ||
        (1-p +  5*p/15 < rp3 && rp3 < 1-p +  7*p/15) ||
        (1-p +  9*p/15 < rp3 && rp3 < 1-p + 11*p/15) ||
        (1-p + 13*p/15 < rp3 && rp3 < 1            )) { // Z_i
          Z_vec[i] *= -1;
    }

    // after H_Q
    double rp4 = rand();
    if (rp4 < 2*p/3) {  // X_Q
      Z_vec[i] *= -1;
    }

    if (p/3 < rp4 && rp4 < p) { // Z_Q
      if (i + 1 < L*M*N) {
        Z_vec[i + 1] *= -1;
      }
    }
  }

  // Z measurements to remove qubits
  for (int v = 0; v < Vextra.size(); v++) {
    int i = Vextra[v];  // i = qubit number
    double rpZ = rand();
    if (rpZ < 2*p/3) {
      if (i - L*M >= 0) {
        Z_vec[i - L*M] *= -1;
      }
      if (i - L >= 0) {
        Z_vec[i - L] *= -1;
      }
      if (i - 1 >= 0) {
        Z_vec[i - 1] *= -1;
      }
      if (i + 1 < L*M*N) {
        Z_vec[i + 1] *= -1;
      }
      if (i + L < L*M*N) {
        Z_vec[i + L] *= -1;
      }
      if (i + L*M < L*M*N) {
        Z_vec[i + L*M] *= -1;
      }
    }
  }

  // X measurements to obtain syndrome
  for (int v = 0; v < Vbcc.size(); v++) {
    int i = Vbcc[v];  // i = qubit number
    double rpX = rand();
    if (rpX < 2*p/3) {
      Z_vec[i] *= -1;
    }
  }

  // // TODO put into noise_model so it's compatible with CLI
  // // first, cut out boundaries
  // std::vector<int> Z_vec_smaller;
  // L-=2;M-=2;N-=2;
  // int LL = L + 2;
  // int MM = M + 2;
  // int NN = N + 2;
  //
  // for (int k = 1; k < NN - 1; k++) {
  //   for (int j = 1; j < MM - 1; j++) {
  //     for (int i = 1; i < LL - 1; i++) {
  //       Z_vec_smaller.push_back(Z_vec[k*LL*MM + j*LL + i]);
  //     }
  //   }
  // }
  //
  // return Z_vec_smaller;

  return Z_vec;
}

vector<int> smallZ_vec(int L, int M, int N, double p, std::function<double()> rand) {  // here, L, M, N are dimensions of the SMALLER lattice
  // note: bigZ_vec has size (L + 2) * (M + 2) * (N + 2)
  int LL = L + 2;
  int MM = M + 2;
  int NN = N + 2;
  vector<int> bigZ_vec = big_cubic_noise(LL, MM, NN, p, rand);

  // cut out boundaries
  vector<int> smallZ_vec;
  for (int k = 1; k < NN - 1; k++) {
    for (int j = 1; j < MM - 1; j++) {
      for (int i = 1; i < LL - 1; i++) {
        smallZ_vec.push_back(bigZ_vec[k*LL*MM + j*LL + i]);
      }
    }
  }
  return smallZ_vec;
}


// returns a binary vector indicating lost qubits
vector<int> big_loss(int L, int M, int N, double p_loss, std::function<double()> rand) { // loss_vec[i] = 1 if qubit i (on the BIGGER lattice) is lost
  vector<int> loss_vec(L*M*N, 0);

  // qubits in the final BCC lattice
  vector<int> Vbcc = get_Vbcc(L - 2, M - 2, N - 2);
  // extra qubits in the bigger lattice which must be removed
  vector<int> Vextra = get_Vextra(L - 2, M - 2, N - 2);

  for (int v = 0; v < Vbcc.size(); v++) {
    int i = Vbcc[v];  // i = qubit number
    double rp0 = rand();
    if (rp0 < p_loss) {
      loss_vec[i] = 1;
    }
  }

  // if a qubit in Vextra is lost, all its neighbors must be lost too
  for (int v = 0; v < Vextra.size(); v++) {
    int i = Vextra[v];  // i = qubit number
    double rp1 = rand();
    if (rp1 < p_loss) {
      loss_vec[i] = 1;
      // also need to treat all the neighbours of i as being lost
      if (i - L*M >= 0) {
        loss_vec[i - L*M] = 1;
      }
      if (i - L >= 0) {
        loss_vec[i - L] = 1;
      }
      if (i - 1 >= 0) {
        loss_vec[i - 1] = 1;
      }
      if (i + 1 < L*M*N) {
        loss_vec[i + 1] = 1;
      }
      if (i + L < L*M*N) {
        loss_vec[i + L] = 1;
      }
      if (i + L*M < L*M*N) {
        loss_vec[i + L*M] = 1;
      }
    }
  }

  return loss_vec;
}

// trims out the fat -- removes the boundary qubits to just give a standard sized lost qubits indicator vector
vector<int> small_loss_vec(int L, int M, int N, double p_loss, std::function<double()> rand) {
  // note: big_loss_vec has size (L + 2) * (M + 2) * (N + 2)
  int LL = L + 2;
  int MM = M + 2;
  int NN = N + 2;
  vector<int> big_loss_vec = big_loss(LL, MM, NN, p_loss, rand);

  // cut out boundaries
  vector<int> small_loss_vec;
  for (int k = 1; k < NN - 1; k++) {
    for (int j = 1; j < MM - 1; j++) {
      for (int i = 1; i < LL - 1; i++) {
        small_loss_vec.push_back(big_loss_vec[k*LL*MM + j*LL + i]);
      }
    }
  }
  return small_loss_vec;
}


void standardNoise_bcc(int L, int M, int N, double p, vector<int> &Z_vec, std::function<double()> rand) {
  vector<int> Vxy = get_Vxy(L, M, N);
  vector<int> Vyz = get_Vyz(L, M, N);
  vector<int> Vzx = get_Vzx(L, M, N);

  for (int v = 0; v < Vxy.size(); v++) {
    int i = Vxy[v]; // qubit number
    int X = (i % (L*M)) % L;
    int Y = (i % (L*M)) / L;
    int Z = i / (L*M);

    // after initialisation
    double rp0 = rand();
    if (rp0 < 2*p/3) {
      if ((i + L < L*M*N) && (Y != M - 1)) {  // no error if i + L isn't actually a neighbour of i
        Z_vec[i + L] *= -1;
      }
    }


    // after Z_{Q,i-L} (which is applied only for i-L >= 0 AND (i - L, i) \in E, i.e., Y > 0)
    if ((i >= L) && (Y != 0)) {
      double rp1 = rand();

      if (1-p + 3*p/15 < rp1 && rp1 < 1-p + 11*p/15) {  // X_Q
        Z_vec[i - L] *= -1;
        if (X != 0) {
          Z_vec[i - 1] *= -1;
        }
      }

      if (1-p + 7*p/15 < rp1 && rp1 < 1) {  // Z_Q
        Z_vec[i] *= -1;
      }
      // *note: Y_Q if 1-p + 7*p/15 < rp1 && rp1 < 1-p + 11*p/15

      if ((1-p           < rp1 && rp1 < 1-p +  2*p/15) ||
          (1-p +  4*p/15 < rp1 && rp1 < 1-p +  6*p/15) ||
          (1-p +  8*p/15 < rp1 && rp1 < 1-p + 10*p/15) ||
          (1-p + 12*p/15 < rp1 && rp1 < 1-p + 14*p/15)) { // X_{i-L}
            if (i - L + L*M < L*M*N) {
              Z_vec[i - L + L*M] *= -1;
            }
      }

      if ((1-p +    p/15 < rp1 && rp1 < 1-p +  3*p/15) ||
          (1-p +  5*p/15 < rp1 && rp1 < 1-p +  7*p/15) ||
          (1-p +  9*p/15 < rp1 && rp1 < 1-p + 11*p/15) ||
          (1-p + 13*p/15 < rp1 && rp1 < 1            )) { // Z_{i-L}
            Z_vec[i - L] *= -1;
      }
      // *again, this includes Y_{i-L} errors
    }

    // after X_{Q,i}
    double rp2 = rand();

    if (1-p + 3*p/15 < rp2 && rp2 < 11*p/15) {  // X_Q
      if ((i + 1 < L*M*N) && (X != L - 1)) {
        Z_vec[i + 1] *= -1;
      }
    }

    if (1-p + 7*p/15 < rp2 && rp2 < 1) {  // Z_Q
      Z_vec[i] *= -1;
    }

    if ((1-p           < rp2 && rp2 < 1-p +  2*p/15) ||
        (1-p +  4*p/15 < rp2 && rp2 < 1-p +  6*p/15) ||
        (1-p +  8*p/15 < rp2 && rp2 < 1-p + 10*p/15) ||
        (1-p + 12*p/15 < rp2 && rp2 < 1-p + 14*p/15)) { // X_i
          if ((i + L < L*M*N) && (Y != M - 1)) {
            Z_vec[i + L] *= -1;
          }
    }

    if ((1-p +    p/15 < rp2 && rp2 < 1-p +  3*p/15) ||
        (1-p +  5*p/15 < rp2 && rp2 < 1-p +  7*p/15) ||
        (1-p +  9*p/15 < rp2 && rp2 < 1-p + 11*p/15) ||
        (1-p + 13*p/15 < rp2 && rp2 < 1            )) { // Z_i
          Z_vec[i] *= -1;
    }

    // after H_Q
    double rp3 = rand();

    if (rp3 < 2*p/3) {
      Z_vec[i] *= -1;
    }

    if (p/3 < rp3 && rp3 < p) {
      if ((i + 1 < L*M*N) && (X != L - 1)) {
        Z_vec[i + 1] *= -1;
      }
    }

  }

  for (int v = 0; v < Vyz.size(); v++) {
    int i = Vyz[v]; // qubit number
    int X = (i % (L*M)) % L;
    int Y = (i % (L*M)) / L;
    int Z = i / (L*M);

    // after initialisation
    double rp0 = rand();
    if (rp0 < 2*p/3) {
      if ((i + L < L*M*N) && (Y != M - 1)) {
        Z_vec[i + L] *= -1;
      }
      if (i + L*M < L*M*N) {
        Z_vec[i + L*M] *= -1;
      }
    }

    // after Z_{Q,i-LM} (which is only applied for i-LM >= 0 AND (i - LM, i) \in E)
    if (i >= L*M) { // if i >= L*M, then (i - LM, i) \in E (for i \in V^{yz})
      double rp1 = rand();

      if (1-p + 3*p/15 < rp1 && rp1 < 1-p + 11*p/15) {  // X_Q
        Z_vec[i - L*M] *= -1;
      }

      if (1-p + 7*p/15 < rp1 && rp1 < 1) {  // Z_Q
        Z_vec[i] *= -1;
      }

      // if an X_{i-LM} error occurs after Z_{Q,i-LM}, the effective error is X_{i-LM} -- doesn't affect syndrome

      if ((1-p +    p/15 < rp1 && rp1 < 1-p +  3*p/15) ||
          (1-p +  5*p/15 < rp1 && rp1 < 1-p +  7*p/15) ||
          (1-p +  9*p/15 < rp1 && rp1 < 1-p + 11*p/15) ||
          (1-p + 13*p/15 < rp1 && rp1 < 1            )) { // Z_{i-LM}
            Z_vec[i - L*M] *= -1;
      }
    }

    // after Z_{Q,i-L} (which is only applied for i-L >= 0 AND (i - L, i) \in E)
    if ((i >= L) && (Y != 0)) {
      double rp2 = rand();

      if (1-p + 3*p/15 < rp2 && rp2 < 1-p + 11*p/15) {  // X_Q
        if (i - L*M >= 0) {
          Z_vec[i - L*M] *= -1;
        }
        Z_vec[i - L] *= -1;
      }

      if (1-p + 7*p/15 < rp2 && rp2 < 1) {  // Z_Q
        Z_vec[i] *= -1;
      }

      // if an X_{i-L} error occurs after Z_{Q,i-L}, the effective error is X_{i-L} -- doesn't affect syndrome

      if ((1-p +    p/15 < rp2 && rp2 < 1-p +  3*p/15) ||
          (1-p +  5*p/15 < rp2 && rp2 < 1-p +  7*p/15) ||
          (1-p +  9*p/15 < rp2 && rp2 < 1-p + 11*p/15) ||
          (1-p + 13*p/15 < rp2 && rp2 < 1            )) { // Z_{i-L}
            Z_vec[i - L] *= -1;
      }
    }

    // after X_{Q,i}
    double rp3 = rand();

    // if an X_Q error occurs after X_{Q,i} for i \in Vyz, no effective error *FOR THE MEASUREMENT PROTOCOL*

    if (1-p + 7*p/15 < rp3 && rp3 < 1) {  // Z_Q
      Z_vec[i] *= -1;
    }

    if ((1-p           < rp3 && rp3 < 1-p +  2*p/15) ||
        (1-p +  4*p/15 < rp3 && rp3 < 1-p +  6*p/15) ||
        (1-p +  8*p/15 < rp3 && rp3 < 1-p + 10*p/15) ||
        (1-p + 12*p/15 < rp3 && rp3 < 1-p + 14*p/15)) { // X_i
          if ((i + L < L*M*N) && (Y != M - 1)) {
            Z_vec[i + L] *= -1;
          }
          if (i + L*M < L*M*N) {
            Z_vec[i + L*M] *= -1;
          }
    }

    if ((1-p +    p/15 < rp3 && rp3 < 1-p +  3*p/15) ||
        (1-p +  5*p/15 < rp3 && rp3 < 1-p +  7*p/15) ||
        (1-p +  9*p/15 < rp3 && rp3 < 1-p + 11*p/15) ||
        (1-p + 13*p/15 < rp3 && rp3 < 1            )) { // Z_i
          Z_vec[i] *= -1;
    }

    // after H_Q
    double rp4 = rand();
    if (rp4 < 2*p/3) {
      Z_vec[i] *= -1;
    }

    // before measurement of Q
    double rp5 = rand();
    if (rp5 < 2*p/3) {
      Z_vec[i] *= -1;
    }

    // after re-initialisation of Q
    double rp6 = rand();
    if (rp6 < 2*p/3) {
      if (i + 2 < L*M*N) {
        Z_vec[i + 2] *= -1;
      }
    }
  }

  for (int v = 0; v < Vzx.size(); v++) {
    int i = Vzx[v]; // qubit number
    int X = (i % (L*M)) % L;
    int Y = (i % (L*M)) / L;
    int Z = i / (L*M);

    // after initialisation
    double rp0 = rand();
    if (rp0 < 2*p/3) {
      if (i + L*M < L*M*N) {
        Z_vec[i + L*M] *= -1;
      }
    }

    // after Z_{Q,i-LM} (which is only applied for i-LM >= 0)
    if (i >= L*M) {
      double rp1 = rand();

      if (1-p + 3*p/15 < rp1 && rp1 < 1-p + 11*p/15) {  // X_Q
        Z_vec[i - L*M] *= -1;
        if (X != 0) {
          Z_vec[i - 1] *= -1;
        }
      }

      if (1-p + 7*p/15 < rp1 && rp1 < 1) {  // Z_Q
        Z_vec[i] *= -1;
      }

      // if an X_{i-L} error occurs after Z_{Q,i-L}, the effective error is X_{i-L} -- doesn't affect syndrome

      if ((1-p +    p/15 < rp1 && rp1 < 1-p +  3*p/15) ||
          (1-p +  5*p/15 < rp1 && rp1 < 1-p +  7*p/15) ||
          (1-p +  9*p/15 < rp1 && rp1 < 1-p + 11*p/15) ||
          (1-p + 13*p/15 < rp1 && rp1 < 1            )) { // Z_{i-LM}
            Z_vec[i - L*M] *= -1;
      }
    }

    // after X_{Q,i}
    double rp2 = rand();

    if (1-p + 3*p/15 < rp2 && rp2 < 11*p/15) {  // X_Q
      if ((i + 1 < L*M*N) && (X != L - 1)) {
        Z_vec[i + 1] *= -1;
      }
    }

    if (1-p + 7*p/15 < rp2 && rp2 < 1) {  // Z_Q
      Z_vec[i] *= -1;
    }

    if ((1-p           < rp2 && rp2 < 1-p +  2*p/15) ||
        (1-p +  4*p/15 < rp2 && rp2 < 1-p +  6*p/15) ||
        (1-p +  8*p/15 < rp2 && rp2 < 1-p + 10*p/15) ||
        (1-p + 12*p/15 < rp2 && rp2 < 1-p + 14*p/15)) { // X_i
          if (i + L*M < L*M*N) {
            Z_vec[i + L*M] *= -1;
          }
    }

    if ((1-p +    p/15 < rp2 && rp2 < 1-p +  3*p/15) ||
        (1-p +  5*p/15 < rp2 && rp2 < 1-p +  7*p/15) ||
        (1-p +  9*p/15 < rp2 && rp2 < 1-p + 11*p/15) ||
        (1-p + 13*p/15 < rp2 && rp2 < 1            )) { // Z_i
          Z_vec[i] *= -1;
    }

    // after H_Q
    double rp3 = rand();

    if (rp3 < 2*p/3) {
      Z_vec[i] *= -1;
    }

    if (p/3 < rp3 && rp3 < p) {
      if ((i + 1 < L*M*N) && (X != L - 1)) {
        Z_vec[i + 1] *= -1;
      }
    }
  }

  // syndrome measurement errors!
  for (int i = 0; i < Z_vec.size(); i++) {
    double rp = rand();
    if (rp < 2*p/3) {
      Z_vec[i] *= -1;
    }
  }

  // initialisation of Q at the very beginning; note that the first qubit (qubit 0) is in Vxy
  double rpi = rand();
  if (rpi < 2*p/3) {
    Z_vec[0] *= -1;
  }

  // measurement of Q at the very end (to remove it from the cluster)
  double rpf = rand();
  if (rpf < 2*p/3) {
    Z_vec[L*M*N - 1] *= -1;
  }

}

void quasistandardNoise_bcc(int L, int M, int N, double p, vector<int> &Z_vec, std::function<double()> rand) {
    // no initialisation errors (subsumed into CNOT errors), and no error for reinitialisation of Q
  vector<int> Vxy = get_Vxy(L, M, N);
  vector<int> Vyz = get_Vyz(L, M, N);
  vector<int> Vzx = get_Vzx(L, M, N);

  for (int v = 0; v < Vxy.size(); v++) {
    int i = Vxy[v]; // qubit number
    int X = (i % (L*M)) % L;
    int Y = (i % (L*M)) / L;
    int Z = i / (L*M);

    // // after initialisation
    // double rp0 = rand();
    // if (rp0 < 2*p/3) {
    //   if ((i + L < L*M*N) && (Y != M - 1)) {  // no error if i + L isn't actually a neighbour of i
    //     Z_vec[i + L] *= -1;
    //   }
    // }


    // after Z_{Q,i-L} (which is applied only for i-L >= 0 AND (i - L, i) \in E, i.e., Y > 0)
    if ((i >= L) && (Y != 0)) {
      double rp1 = rand();

      if (1-p + 3*p/15 < rp1 && rp1 < 1-p + 11*p/15) {  // X_Q
        Z_vec[i - L] *= -1;
        if (X != 0) {
          Z_vec[i - 1] *= -1;
        }
      }

      if (1-p + 7*p/15 < rp1 && rp1 < 1) {  // Z_Q
        Z_vec[i] *= -1;
      }
      // *note: Y_Q if 1-p + 7*p/15 < rp1 && rp1 < 1-p + 11*p/15

      if ((1-p           < rp1 && rp1 < 1-p +  2*p/15) ||
          (1-p +  4*p/15 < rp1 && rp1 < 1-p +  6*p/15) ||
          (1-p +  8*p/15 < rp1 && rp1 < 1-p + 10*p/15) ||
          (1-p + 12*p/15 < rp1 && rp1 < 1-p + 14*p/15)) { // X_{i-L}
            if (i - L + L*M < L*M*N) {
              Z_vec[i - L + L*M] *= -1;
            }
      }

      if ((1-p +    p/15 < rp1 && rp1 < 1-p +  3*p/15) ||
          (1-p +  5*p/15 < rp1 && rp1 < 1-p +  7*p/15) ||
          (1-p +  9*p/15 < rp1 && rp1 < 1-p + 11*p/15) ||
          (1-p + 13*p/15 < rp1 && rp1 < 1            )) { // Z_{i-L}
            Z_vec[i - L] *= -1;
      }
      // *again, this includes Y_{i-L} errors
    }

    // after X_{Q,i}
    double rp2 = rand();

    if (1-p + 3*p/15 < rp2 && rp2 < 11*p/15) {  // X_Q
      if ((i + 1 < L*M*N) && (X != L - 1)) {
        Z_vec[i + 1] *= -1;
      }
    }

    if (1-p + 7*p/15 < rp2 && rp2 < 1) {  // Z_Q
      Z_vec[i] *= -1;
    }

    if ((1-p           < rp2 && rp2 < 1-p +  2*p/15) ||
        (1-p +  4*p/15 < rp2 && rp2 < 1-p +  6*p/15) ||
        (1-p +  8*p/15 < rp2 && rp2 < 1-p + 10*p/15) ||
        (1-p + 12*p/15 < rp2 && rp2 < 1-p + 14*p/15)) { // X_i
          if ((i + L < L*M*N) && (Y != M - 1)) {
            Z_vec[i + L] *= -1;
          }
    }

    if ((1-p +    p/15 < rp2 && rp2 < 1-p +  3*p/15) ||
        (1-p +  5*p/15 < rp2 && rp2 < 1-p +  7*p/15) ||
        (1-p +  9*p/15 < rp2 && rp2 < 1-p + 11*p/15) ||
        (1-p + 13*p/15 < rp2 && rp2 < 1            )) { // Z_i
          Z_vec[i] *= -1;
    }

    // after H_Q
    double rp3 = rand();

    if (rp3 < 2*p/3) {
      Z_vec[i] *= -1;
    }

    if (p/3 < rp3 && rp3 < p) {
      if ((i + 1 < L*M*N) && (X != L - 1)) {
        Z_vec[i + 1] *= -1;
      }
    }

  }

  for (int v = 0; v < Vyz.size(); v++) {
    int i = Vyz[v]; // qubit number
    int X = (i % (L*M)) % L;
    int Y = (i % (L*M)) / L;
    int Z = i / (L*M);

    // // after initialisation
    // double rp0 = rand();
    // if (rp0 < 2*p/3) {
    //   if ((i + L < L*M*N) && (Y != M - 1)) {
    //     Z_vec[i + L] *= -1;
    //   }
    //   if (i + L*M < L*M*N) {
    //     Z_vec[i + L*M] *= -1;
    //   }
    // }

    // after Z_{Q,i-LM} (which is only applied for i-LM >= 0 AND (i - LM, i) \in E)
    if (i >= L*M) { // if i >= L*M, then (i - LM, i) \in E (for i \in V^{yz})
      double rp1 = rand();

      if (1-p + 3*p/15 < rp1 && rp1 < 1-p + 11*p/15) {  // X_Q
        Z_vec[i - L*M] *= -1;
      }

      if (1-p + 7*p/15 < rp1 && rp1 < 1) {  // Z_Q
        Z_vec[i] *= -1;
      }

      // if an X_{i-LM} error occurs after Z_{Q,i-LM}, the effective error is X_{i-LM} -- doesn't affect syndrome

      if ((1-p +    p/15 < rp1 && rp1 < 1-p +  3*p/15) ||
          (1-p +  5*p/15 < rp1 && rp1 < 1-p +  7*p/15) ||
          (1-p +  9*p/15 < rp1 && rp1 < 1-p + 11*p/15) ||
          (1-p + 13*p/15 < rp1 && rp1 < 1            )) { // Z_{i-LM}
            Z_vec[i - L*M] *= -1;
      }
    }

    // after Z_{Q,i-L} (which is only applied for i-L >= 0 AND (i - L, i) \in E)
    if ((i >= L) && (Y != 0)) {
      double rp2 = rand();

      if (1-p + 3*p/15 < rp2 && rp2 < 1-p + 11*p/15) {  // X_Q
        if (i - L*M >= 0) {
          Z_vec[i - L*M] *= -1;
        }
        Z_vec[i - L] *= -1;
      }

      if (1-p + 7*p/15 < rp2 && rp2 < 1) {  // Z_Q
        Z_vec[i] *= -1;
      }

      // if an X_{i-L} error occurs after Z_{Q,i-L}, the effective error is X_{i-L} -- doesn't affect syndrome

      if ((1-p +    p/15 < rp2 && rp2 < 1-p +  3*p/15) ||
          (1-p +  5*p/15 < rp2 && rp2 < 1-p +  7*p/15) ||
          (1-p +  9*p/15 < rp2 && rp2 < 1-p + 11*p/15) ||
          (1-p + 13*p/15 < rp2 && rp2 < 1            )) { // Z_{i-L}
            Z_vec[i - L] *= -1;
      }
    }

    // after X_{Q,i}
    double rp3 = rand();

    // if an X_Q error occurs after X_{Q,i} for i \in Vyz, no effective error *FOR THE MEASUREMENT PROTOCOL*

    if (1-p + 7*p/15 < rp3 && rp3 < 1) {  // Z_Q
      Z_vec[i] *= -1;
    }

    if ((1-p           < rp3 && rp3 < 1-p +  2*p/15) ||
        (1-p +  4*p/15 < rp3 && rp3 < 1-p +  6*p/15) ||
        (1-p +  8*p/15 < rp3 && rp3 < 1-p + 10*p/15) ||
        (1-p + 12*p/15 < rp3 && rp3 < 1-p + 14*p/15)) { // X_i
          if ((i + L < L*M*N) && (Y != M - 1)) {
            Z_vec[i + L] *= -1;
          }
          if (i + L*M < L*M*N) {
            Z_vec[i + L*M] *= -1;
          }
    }

    if ((1-p +    p/15 < rp3 && rp3 < 1-p +  3*p/15) ||
        (1-p +  5*p/15 < rp3 && rp3 < 1-p +  7*p/15) ||
        (1-p +  9*p/15 < rp3 && rp3 < 1-p + 11*p/15) ||
        (1-p + 13*p/15 < rp3 && rp3 < 1            )) { // Z_i
          Z_vec[i] *= -1;
    }

    // after H_Q
    double rp4 = rand();
    if (rp4 < 2*p/3) {
      Z_vec[i] *= -1;
    }

    // before measurement of Q
    double rp5 = rand();
    if (rp5 < 2*p/3) {
      Z_vec[i] *= -1;
    }

    // // after re-initialisation of Q
    // double rp6 = rand();
    // if (rp6 < 2*p/3) {
    //   if (i + 2 < L*M*N) {
    //     Z_vec[i + 2] *= -1;
    //   }
    // }
  }

  for (int v = 0; v < Vzx.size(); v++) {
    int i = Vzx[v]; // qubit number
    int X = (i % (L*M)) % L;
    int Y = (i % (L*M)) / L;
    int Z = i / (L*M);

    // // after initialisation
    // double rp0 = rand();
    // if (rp0 < 2*p/3) {
    //   if (i + L*M < L*M*N) {
    //     Z_vec[i + L*M] *= -1;
    //   }
    // }

    // after Z_{Q,i-LM} (which is only applied for i-LM >= 0)
    if (i >= L*M) {
      double rp1 = rand();

      if (1-p + 3*p/15 < rp1 && rp1 < 1-p + 11*p/15) {  // X_Q
        Z_vec[i - L*M] *= -1;
        if (X != 0) {
          Z_vec[i - 1] *= -1;
        }
      }

      if (1-p + 7*p/15 < rp1 && rp1 < 1) {  // Z_Q
        Z_vec[i] *= -1;
      }

      // if an X_{i-L} error occurs after Z_{Q,i-L}, the effective error is X_{i-L} -- doesn't affect syndrome

      if ((1-p +    p/15 < rp1 && rp1 < 1-p +  3*p/15) ||
          (1-p +  5*p/15 < rp1 && rp1 < 1-p +  7*p/15) ||
          (1-p +  9*p/15 < rp1 && rp1 < 1-p + 11*p/15) ||
          (1-p + 13*p/15 < rp1 && rp1 < 1            )) { // Z_{i-LM}
            Z_vec[i - L*M] *= -1;
      }
    }

    // after X_{Q,i}
    double rp2 = rand();

    if (1-p + 3*p/15 < rp2 && rp2 < 11*p/15) {  // X_Q
      if ((i + 1 < L*M*N) && (X != L - 1)) {
        Z_vec[i + 1] *= -1;
      }
    }

    if (1-p + 7*p/15 < rp2 && rp2 < 1) {  // Z_Q
      Z_vec[i] *= -1;
    }

    if ((1-p           < rp2 && rp2 < 1-p +  2*p/15) ||
        (1-p +  4*p/15 < rp2 && rp2 < 1-p +  6*p/15) ||
        (1-p +  8*p/15 < rp2 && rp2 < 1-p + 10*p/15) ||
        (1-p + 12*p/15 < rp2 && rp2 < 1-p + 14*p/15)) { // X_i
          if (i + L*M < L*M*N) {
            Z_vec[i + L*M] *= -1;
          }
    }

    if ((1-p +    p/15 < rp2 && rp2 < 1-p +  3*p/15) ||
        (1-p +  5*p/15 < rp2 && rp2 < 1-p +  7*p/15) ||
        (1-p +  9*p/15 < rp2 && rp2 < 1-p + 11*p/15) ||
        (1-p + 13*p/15 < rp2 && rp2 < 1            )) { // Z_i
          Z_vec[i] *= -1;
    }

    // after H_Q
    double rp3 = rand();

    if (rp3 < 2*p/3) {
      Z_vec[i] *= -1;
    }

    if (p/3 < rp3 && rp3 < p) {
      if ((i + 1 < L*M*N) && (X != L - 1)) {
        Z_vec[i + 1] *= -1;
      }
    }
  }

  // syndrome measurement errors!
  for (int i = 0; i < Z_vec.size(); i++) {
    double rp = rand();
    if (rp < 2*p/3) {
      Z_vec[i] *= -1;
    }
  }

  // // initialisation of Q at the very beginning; note that the first qubit (qubit 0) is in Vxy
  // double rpi = rand();
  // if (rpi < 2*p/3) {
  //   Z_vec[0] *= -1;
  // }

  // // measurement of Q at the very end (to remove it from the cluster)
  // double rpf = rand();
  // if (rpf < 2*p/3) {
  //   Z_vec[L*M*N - 1] *= -1;
  // }

}

void dephasingPerBin_bcc(int L, int M, int N, double p, vector<int> &Z_vec, std::function<double()> rand) {
  vector<int> Vxy = get_Vxy(L, M, N);
  vector<int> Vyz = get_Vyz(L, M, N);
  vector<int> Vzx = get_Vzx(L, M, N);

  for (int v = 0; v < Vxy.size(); v++) {  // iterate over qubits (instead of over gate blocks like in the noise.cpp code)
    int i = Vxy[v]; // qubit number
    int X = (i % (L*M)) % L;
    int Y = (i % (L*M)) / L;
    int Z = i / (L*M);

    // after X_{Q,i} -- L time bins
    for (int j = 0; j < L; j++) {
      double rp0 = rand();
      if (rp0 < p) { // Z_i
        Z_vec[i] *= -1;
      }
    }
  }

  for (int v = 0; v < Vzx.size(); v++) {
    int i = Vzx[v]; // qubit number
    int X = (i % (L*M)) % L;
    int Y = (i % (L*M)) / L;
    int Z = i / (L*M);

    // after X_{Q,i} -- LM time bins
    for (int j = 0; j < L*M; j++) {
      double rp0 = rand();
      if (rp0 < p) { // Z_i
        Z_vec[i] *= -1;
      }
    }
  }

  for (int v = 0; v < Vyz.size(); v++) {
    int i = Vyz[v]; // qubit number
    int X = (i % (L*M)) % L;
    int Y = (i % (L*M)) / L;
    int Z = i / (L*M);

    // after X_{Q,i} -- L time bins
    for (int j = 0; j < L; j++) {
      double rp0 = rand();
      if (rp0 < p) { // Z_i
        Z_vec[i] *= -1;
      }
    }

    // after Z_{Q,i} -- L time bins
    for (int j = 0; j < L*(M-1); j++) {
      double rp1 = rand();
      if (rp1 < p) { // Z_i
        Z_vec[i] *= -1;
      }
    }
  }
}


vector<vector<int> > supercheck_operators_loss_per_time(int L, int M, int N, double rate, std::function<double()> rand) {  // loss rate per time bin
  vector<int> Vxy_primal = get_Vxy_primal(L, M, N);
  vector<int> Vyz_primal = get_Vyz_primal(L, M, N);
  vector<int> Vzx_primal = get_Vzx_primal(L, M, N);

  int l = (L + 1) / 2;
  int m = (M + 1) / 2;
  int n = (N - 1) /2;

  // create adjacency list for connected components problem
  vector<vector<int> > adj(l*m*n + 2);  // front boundary has vertex number l*m*n, back boundary has vertex number l*m*n + 1
  int front_boundary_vertex = l*m*n;
  int back_boundary_vertex = l*m*n + 1;

  for (int i = 0; i < Vxy_primal.size(); i++) { // note: only need to consider PRIMAL qubits that are lost
    double rp = rand();
    if (rp < 1 - exp(-rate * L)) {  // lose qubit Vxy_primal[i] with probability 1 - exp(-rate * L), since Vxy only travels through the first delay line!!
      int q = Vxy_primal[i];
      // cout << q << ", ";
      int X = (q % (L*M)) % L;
      int Y = (q % (L*M)) / L;
      int Z = q / (L*M);
      if (Z == 0) { // if lost face is on the front boundary, add edge between back cube and front boundary vertex
        int back_cube = X/2 + (Y/2)*l + (Z/2)*l*m;
        adj[back_cube].push_back(front_boundary_vertex);
        adj[front_boundary_vertex].push_back(back_cube);
      }
      else if (Z == N - 1) {  // if lost face is on the back boundary, add edge between front cube and back boundary vertex
        int front_cube = X/2 + (Y/2)*l + (Z/2 - 1)*l*m;
        adj[front_cube].push_back(back_boundary_vertex);
        adj[back_boundary_vertex].push_back(front_cube);
      }
      else {  // if lost face is on neither boundary, add edge between front cube and back cube
        int front_cube = X/2 + (Y/2)*l + (Z/2 - 1)*l*m;
        int back_cube = X/2 + (Y/2)*l + (Z/2)*l*m;
        adj[front_cube].push_back(back_cube);
        adj[back_cube].push_back(front_cube);
      }
    }
  }

  for (int i = 0; i < Vyz_primal.size(); i++) {
    double rp = rand();
    if (rp < 1 - exp(-rate * L * M)) {  // lose qubit Vyz_primal[i] with probability 1 - exp(-rate * L * M), since Vyz goes through both delay lines
      int q = Vyz_primal[i];
      // cout << q << ", ";
      int X = (q % (L*M)) % L;
      int Y = (q % (L*M)) / L;
      int Z = q / (L*M);
      int left_cube = (X - 1)/2 + (Y/2)*l + ((Z-1)/2)*l*m;
      int right_cube = (X - 1)/2 + 1 + (Y/2)*l + ((Z-1)/2)*l*m;
      adj[left_cube].push_back(right_cube);
      adj[right_cube].push_back(left_cube);
    }
  }

  for (int i = 0; i < Vzx_primal.size(); i++) {
    double rp = rand();
    if (rp < 1 - exp(-rate * L * M)) {  // lose qubit Vzx_primal[i] with probability 1 - exp(-rate * L * M), since Vzx goes through both delay lines
      int q = Vzx_primal[i];
      // cout << q << ", ";
      int X = (q % (L*M)) % L;
      int Y = (q % (L*M)) / L;
      int Z = q / (L*M);
      int top_cube = X/2 + ((Y-1)/2)*l + ((Z-1)/2)*l*m;
      int bottom_cube = X/2 + ((Y-1)/2 + 1)*l + ((Z-1)/2)*l*m;
      adj[top_cube].push_back(bottom_cube);
      adj[bottom_cube].push_back(top_cube);
    }
  }

  // cout << "\n";
  vector<vector<int> > superchunks = connected_components(adj);

  sort(superchunks.begin(), superchunks.end()); // last vector in superchunks is the connected component containing the back boundary vertex
  // for (int i = 0; i < superchunks.size(); i++) {
    // print_vec(superchunks[i]);
  // }

  return superchunks;

}
