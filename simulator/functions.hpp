using namespace std;

void print_vec(const vector<int> &vec);

void print_vec_double(const vector<double> &vec);

vector<int> get_Vxy(int L, int M, int N);

vector<int> get_Vyz(int L, int M, int N);

vector<int> get_Vzx(int L, int M, int N);

vector<int> get_Vxy_primal(int L, int M, int N);

vector<int> get_Vxy_dual(int L, int M, int N);

vector<int> get_Vyz_primal(int L, int M, int N);

vector<int> get_Vyz_dual(int L, int M, int N);

vector<int> get_Vzx_primal(int L, int M, int N);

vector<int> get_Vzx_dual(int L, int M, int N);

vector<int> get_Vprimal(int L, int M, int N);

vector<int> get_faces(int c, int L, int M, int N);

void DFS(int v, const vector<vector<int> > &adj, bool visited[], vector<int> &components);

vector<vector<int> > connected_components(const vector<vector<int> > &adj);

int taxicab_distance(int c1, int c2, int l, int m, int n);

int shortest_distance(const vector<int> &chunk1, const vector<int> &chunk2, int l, int m, int n);

int decode(int L, int M, int N, const std::vector<int>& Z_vec, double p, std::function<double()> rand);

int loss_decode(int L, int M, int N, const vector<int> &Z_vec, const vector<vector<int> > &superchunks, std::function<double()> rand);

int weighted_decode(int L, int M, int N, const vector<int> &Z_vec, double p, std::function<double()> rand);

void phenom_noise(double p, vector<int> &Z_vec, vector<int> &loss_vec);
