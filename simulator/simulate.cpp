# include <iostream>
# include <random>
# include <vector>
# include <cstdlib>
# include <functional>
# include <algorithm>
# include <fstream>
# include "cxxopts.hpp"

# include <assert.h>

using namespace std;
// todo: change to for(auto: ...) iteration where possible

# include "functions.hpp"
# include "noise_models.hpp"

// loss errors + decoding
# define demand(cond, str) if(!(cond)){std::cout << str << std::endl;exit(1);}


vector<vector<int> > supercheck_operators(int L, int M, int N, std::vector<int> loss_vec) {
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

  // TODO: refactor -- take in a loss_vec and for each of the following q's, don't flip a coin but instead use the loss_vec[q] to tell us whether there was an error
  for (int i = 0; i < Vxy_primal.size(); i++) { // note: only need to consider PRIMAL qubits that are lost
    // double rp = rand();
    // if (rp < p_loss) {  // lose qubit Vxy_primal[i] with probability p_loss
    int q = Vxy_primal[i];
    if (loss_vec[q]) {
      //      cout << q << ", ";
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
    int q = Vyz_primal[i];
    if (loss_vec[q]) {
    // double rp = rand();
    // if (rp < p_loss) {  // lose qubit Vyz_primal[i] with probability p_loss
      // int q = Vyz_primal[i];
      //      cout << q << ", ";
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
    // double rp = rand();
    int q = Vzx_primal[i];
    if (loss_vec[q]) {
    // if (rp < p_loss) {  // lose qubit Vzx_primal[i] with probability p_loss
      // int q = Vzx_primal[i];
      //      cout << q << ", ";
      int X = (q % (L*M)) % L;
      int Y = (q % (L*M)) / L;
      int Z = q / (L*M);
      int top_cube = X/2 + ((Y-1)/2)*l + ((Z-1)/2)*l*m;
      int bottom_cube = X/2 + ((Y-1)/2 + 1)*l + ((Z-1)/2)*l*m;
      adj[top_cube].push_back(bottom_cube);
      adj[bottom_cube].push_back(top_cube);
    }
  }

  //  cout << "\n";
  vector<vector<int> > superchunks = connected_components(adj);

  sort(superchunks.begin(), superchunks.end()); // last vector in superchunks is the connected component containing the back boundary vertex
  /*  for (int i = 0; i < superchunks.size(); i++) {
    print_vec(superchunks[i]);
    }*/

  return superchunks;

}

vector<vector<int> > supercheck_operators(int L, int M, int N, double p_loss, std::function<double()> rand) {
    std::vector<int> loss_vec(L*M*N, 0);
    for (int q=0; q<L*M*N; q++) {
        double rp = rand();
        if (rp<p_loss)
            loss_vec[q] = 1;
    }
    return supercheck_operators(L, M, N, loss_vec);
}

// typedef std::vector <int> (*NoiseModel)(int, int, int, double, std::function <double()>);
// NoiseModel noise_model;
typedef std::function<std::vector<int>(int, int, int, double, std::function <double()>)> NoiseModel;

void run_sims(
    NoiseModel noise,
    int decoder (int, int, int, const std::vector<int>&, double, std::function<double()>),
    std::vector<int> Ls,
    double p1, double p2,
    int points, int trials,
    const char* fname,
    std::function<bool(int,int)> stopping_condition,
    std::function<double()> rand
) {
    std::ofstream outfile;
    outfile.open(fname);
    outfile << "L,p_error,num_success,num_fail\n";

    std::vector<double> ps;
    for (int j = 0; j <= points; j++) {
      double p = p1 + j * (p2 - p1) / points;
      ps.push_back(p);
    }

    std::cout << "~ps~ \n";
    print_vec_double(ps);
    std::cout << "\n";

    for (auto L: Ls) {
        for (auto p: ps) {
            std::cout << "L=" << L << " p=" << p << "\n";
            int num_success = 0;
            int k=0;
            while (1) {
                std::vector<int> Z_vec = noise(L, L, L, p, rand);
                num_success += decoder(L, L, L, Z_vec, p, rand);
                k++;
                if (stopping_condition(k, k-num_success))
                    break;
            }
            outfile << L << "," << p << "," << num_success << "," <<
                k-num_success << "\n";
        }
    }
    outfile.close();
}



void run_sims(
    // std::vector<int> noise (int, int, int, double, std::function<double()>),
    NoiseModel noise,
    std::vector<std::vector<int> > loss_superchecks(int, int, int, double, std::function<double()>),
    int decoder (int, int, int, const std::vector<int>&, const vector<vector<int> > &, std::function<double()>),
    std::vector<int> Ls,
    double p1, double p2,
    double pl1, double pl2,
    int points, int points_loss, int trials,
    const char* fname,
    std::function<bool(int,int)> stopping_condition,
    std::function<double()> rand
) {
    std::ofstream outfile;
    outfile.open(fname);
    outfile << "L,p_error,p_loss,num_success,num_fail\n";

    std::vector<double> ps;
    for (int j = 0; j <= points; j++) {
      double p = p1 + j * (p2 - p1) / points;
      ps.push_back(p);
    }

    std::vector<double> pls;
    for (int j = 0; j <= points_loss; j++) {
      double pl = pl1 + j * (pl2 - pl1) / points_loss;
      pls.push_back(pl);
    }


    std::cout << "~ps~ \n";
    print_vec_double(ps);
    std::cout << "\n";

    std::cout << "~pls~ \n";
    print_vec_double(pls);
    std::cout << "\n";

    for (auto L: Ls) {
        for (auto p: ps) {
            for (auto pl: pls) {
                std::cout << "L=" << L << " p=" << p << " pl=" << pl << "\n";
                int num_success = 0;
                int k=0;
                while (1) {
                    std::vector<int> Z_vec = noise(L, L, L, p, rand);
                    vector<vector<int> > superchunks = loss_superchecks(L, L, L, pl, rand);
                    num_success += decoder(L, L, L, Z_vec, superchunks, rand);
                    k++;
                    if (stopping_condition(k, k-num_success))
                        break;
                }
                outfile << L << "," << p << "," << pl << "," << num_success << "," <<
                    k-num_success << "\n";
            }
        }
    }
    outfile.close();
}


int main(int argc, const char * argv[]) {
    // REMOVE RANDOM SEED AFTER DEBUGGING
    random_device rd;
    auto randseed = rd();
    mt19937 gen(randseed); // optimise: maybe seed this with a number, outside the function
    uniform_real_distribution<> dis(0.0, 1.0);

    // if (argc != 11) {
    //     std::cout << "ERROR: invalid arguments -- usage:\n./_ <noise_model: stephens_standard / single_emitter / single_emitter_cubic / local_depolarizing> <p_min> <p_max> <p_loss_min> <p_loss_max> <points> <points_loss> <trials> <randseed> <fname>" << std::endl;
    //     exit(1);
    // }

    std::vector <int> Ls;

    for(auto L: Ls) {
        if (L%2==0) {
            std::cerr << "error: boundary conditions require L to be even" << std::endl;
            exit(1);
        }
    }

    // int argi = 1;
    std::string noise_model_name, decoder_name, fname, Ls_string;// = std::string(argv[argi]);
    double p1, p2, pl1, pl2;
    int points, points_loss, trials, min_errors, seed, max_trials;
    double noise_per_unit_time;
    bool sim_loss = false;
    cxxopts::Options options(
        argv[0],
        "Simulator for fault-tolerant measurement-based quantum "
        "computation on foliated surface code cluster states"
    );
    options.add_options()
    (
        "m,noise-model", "Noise Model Name",
        cxxopts::value(noise_model_name)->default_value("single_emitter")
    )
    (
        "decoder", "Decoder (leave blank to use "
        "standard decoder. Options: (standard/weighted)",
        cxxopts::value(decoder_name)->default_value("standard")
    )
    (
        "p1", "Physical Error Probability Lower Bound",
        cxxopts::value(p1)->default_value("0")
    )
    (
        "p2", "Physical Error Probability Upper Bound",
        cxxopts::value(p2)->default_value("1")
    )
    (
        "loss", "Simulate loss",
        cxxopts::value(sim_loss)->default_value("false")
    )
    (
        "pl1", "Physical Loss Probability Lower Bound",
        cxxopts::value(pl1)->default_value("0")
    )
    (
        "pl2", "Physical Loss Probability Upper Bound",
        cxxopts::value(pl2)->default_value("1")
    )
    (
        "noise-rate", "Noise Rate per Unit Time (only applies for quasistandard error models)",
        cxxopts::value(noise_per_unit_time)->default_value("0.001")
    )
    (
        "points", "Number of grid point values along the "
        "physical error axis",
        cxxopts::value(points)->default_value("3")
    )
    (
        "points-loss", "Number of grid point values along the "
        "physical loss axis",
        cxxopts::value(points_loss)->default_value("3")
    )
    (
        "trials", "Number of trials per grid point",
        cxxopts::value(trials)->default_value("-1")
    )
    (
        "max-trials", "Maximum number of trials per grid point",
        cxxopts::value(max_trials)->default_value("-1")
    )
    (
        "min-errors", "Minimum number of errors per grid point",
        cxxopts::value(min_errors)->default_value("-1")
    )
    (
        "Ls", "comma-delimited list of integer values for L",
        cxxopts::value(Ls_string)->default_value("11,15,21,25,29,31")
    )
    (
        "seed", "Random seed (leave blank to use "
        "hardware-generated random seed)",
        cxxopts::value(seed)->default_value("-1")
    )
    (
        "fname", "Output file name to write simulation results to",
        cxxopts::value(fname)
    );

    auto result = options.parse(argc, argv);
    demand(result.count("fname")==1, options.help());

    demand((min_errors!=-1)^(trials!=-1), "must supply exactly one of --trials --min-errors");
    // typedef bool (*StoppingCondition)(int, int);
    // StoppingCondition stopping_condition;
    std::function<bool(int,int)> stopping_condition;
    // auto stopping_condition = [](){return true;};
    if (min_errors==-1) {
        stopping_condition = [&](int num_trials_run, int num_errors_obs) {
            // std::cout << "num_trials_run = " << num_trials_run << " num_errors_obs = " << num_errors_obs << std::endl;
            if (num_trials_run >= trials)
                return true;

            return false;
        };
    } else {
        assert(trials==-1);
        stopping_condition = [&](int num_trials_run, int num_errors_obs) {
            // std::cout << "num_trials_run = " << num_trials_run << " num_errors_obs = " << num_errors_obs << std::endl;
            if (num_errors_obs >= min_errors)
                return true;

            if (max_trials != -1 && num_trials_run >= max_trials)
                return true;

            return false;
        };
    }

    if (seed == -1) {
        std::cout << "using HARDWARE-GENERATED random seed " << randseed << std::endl;
    } else {
        std::cout << "using USER-SUPPLIED random seed " << seed << std::endl;
        gen = std::mt19937(seed);
    }

    std::stringstream Ls_stringstream = std::stringstream(Ls_string);
    std::string placeholder;
    while(std::getline(Ls_stringstream, placeholder, ','))
    {
        int L = std::stoi(placeholder);
        Ls.push_back(L);
        std::cout << "simulating L=" << L << "\n";
    }

    NoiseModel noise_model;
    typedef std::vector<std::vector<int> > (*LossModel)(int, int, int, double, std::function<double()>);
    LossModel loss_model;
    typedef int (*Decoder)(int, int, int, const std::vector<int>&, double, std::function<double()>);
    Decoder decoder;
    typedef int (*LossDecoder)(int, int, int, const std::vector<int>&, const vector<vector<int> > &, std::function<double()>);
    LossDecoder loss_decoder;

    if (noise_model_name == "stephens_standard") {
        noise_model = noise_stephens_standard_gateerror;
    } else if (noise_model_name == "local_depolarizing") {
        noise_model = noise_RHG_error_model1_local_depolarizing;
    } else if (noise_model_name == "single_emitter") {
        noise_model = noise_single_emitter;
    } else if (noise_model_name == "single_emitter_cubic") {
        noise_model = smallZ_vec;
    } else if (
        noise_model_name == "bcc_quasistandard"
    ) {
        std::cout << "noise_per_unit_time = " << noise_per_unit_time << std::endl;
        noise_model = [&](int L, int M, int N, double p, std::function <double()> rand) {
            std::vector<int> Z_vec(L*M*N, 1);
                                         // p=1e-3
            quasistandardNoise_bcc(L, M, N, noise_per_unit_time, Z_vec, rand);
            return Z_vec;
        };
    } else if (noise_model_name == "bcc_quasistandard_dephasing") {
        std::cout << "noise_per_unit_time = " << noise_per_unit_time << std::endl;
        noise_model = [&](int L, int M, int N, double p, std::function <double()> rand) {
            std::vector<int> Z_vec(L*M*N, 1);
                                         // p=1e-3
            quasistandardNoise_bcc(L, M, N, noise_per_unit_time, Z_vec, rand);
            dephasingPerBin_bcc(L, M, N, p, Z_vec, rand);
            return Z_vec;
        };
    } else {
        std::cerr << "unsupported noise_model: " << noise_model_name << std::endl;
        exit(1);
    }
    if (decoder_name == "standard") {
        std::cout << "using standard decoder\n";
        loss_decoder = loss_decode;
        decoder = decode;
    } else if (decoder_name == "weighted") {
        demand(!sim_loss, "no weighted loss decoder implementation exists");
        std::cout << "using weighted decoder\n";
        decoder = weighted_decode;
    } else {
        std::cerr<<"invalid decoder: "<<decoder_name<<"\n";
        exit(1);
    }

    if (noise_model_name == "single_emitter_cubic") {
        std::cout << "using cubic loss model (superchecks)"<<std::endl;
        loss_model = [](int L, int M, int N, double p_loss, std::function <double()> rand){
            auto loss_vec = small_loss_vec(L, M, N, p_loss, rand);
            return supercheck_operators(L, M, N, loss_vec);
        };
    } else if (
        (noise_model_name == "bcc_quasistandard") ||
        (noise_model_name == "bcc_quasistandard_dephasing")) {
        std::cout << "using BCC gate and delay-line loss model (superchecks)"<<std::endl;
        loss_model = supercheck_operators_loss_per_time;
    } else {
        std::cout << "using non-cubic loss model (superchecks)"<<std::endl;
        loss_model = supercheck_operators;
        // [](int L, int M, int N, std::function <double()> rand){
        //     supercheck_operators(L, M, N, p_loss, rand)
        // };
    }

    std::cout << "p1 = " << p1 <<
        "\np2 = " << p2 <<
        "\nsim_loss = " << sim_loss <<
        "\npl1 = " << pl1 <<
        "\npl2 = " << pl2 <<
        "\npoints = " << points <<
        "\ntrials = " << trials <<
        "\nmin_errors = " << min_errors <<
        "\nout fname = " << fname << std::endl;

    if (sim_loss) {
        std::cout << "simulating loss" << std::endl;
        run_sims(
            noise_model, loss_model, loss_decoder,
            Ls, p1, p2, pl1, pl2,
            points, points_loss, trials, fname.c_str(),
            stopping_condition,
            [ & gen, & dis]() {
                return dis(gen);
            }
        );
    } else {
        std::cout << "not simulating loss" << std::endl;
        run_sims(
            noise_model, decoder, Ls, p1, p2,
            points, trials, fname.c_str(),
            stopping_condition,
            [ & gen, & dis]() {
                return dis(gen);
            }
        );
    }

}
