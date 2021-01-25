# include <iostream>
# include <random>
# include <vector>
# include <cstdlib>
# include <algorithm>
# include <functional>

using namespace std;

# include "functions.hpp"

#include "PerfectMatching.h"
#include "GEOM/GeomPerfectMatching.h"

// REMOVE RANDOM SEED AFTER DEBUGGING?
// random_device rd;
// mt19937 gen(rd()); // optimise: maybe seed this with a number, outside the function
// uniform_real_distribution<> dis(0.0, 1.0);

void print_vec(const vector<int> &vec) {
  for (int i = 0; i < vec.size(); i++) {
    cout << vec[i] << ", ";
  }
  cout << "\n\n";
}

void print_vec_double(const vector<double> &vec) {
  for (int i = 0; i < vec.size(); i++) {
    cout << vec[i] << ", ";
  }
  cout << "\n\n";
}

vector<int> get_Vxy(int L, int M, int N) {
  vector<int> Vxy;
  for (int k = 0; k < N/2 + 1; k++) {
    for (int j = 0; j < M/2 + 1; j++) {
      for (int i = 0; i < L/2 + 1; i++) {
        Vxy.push_back(2*k*L*M + 2*j*L + 2*i);
      }
    }
  }
  for (int k = 0; k < N/2; k++) {
    for (int j = 0; j < M/2; j++) {
      for (int i = 0; i < L/2; i++) {
        Vxy.push_back((2*k + 1)*L*M + (2*j + 1)*L + (2*i + 1));
      }
    }
  }
  sort(Vxy.begin(), Vxy.end());
  return Vxy;
}

vector<int> get_Vyz(int L, int M, int N) {
  vector<int> Vyz;
  for (int k = 0; k < N/2 + 1; k++) {
    for (int j = 0; j < M/2; j++) {
      for (int i = 0; i < L/2 + 1; i++) {
        Vyz.push_back(2*k*L*M + (2*j + 1)*L + 2*i);
      }
    }
  }
  for (int k = 0; k < N/2; k++) {
    for (int j = 0; j < M/2 + 1; j++) {
      for (int i = 0; i < L/2; i++) {
        Vyz.push_back((2*k + 1)*L*M + 2*j*L + (2*i + 1));
      }
    }
  }
  sort(Vyz.begin(), Vyz.end());
  return Vyz;
}

vector<int> get_Vzx(int L, int M, int N) {
  vector<int> Vzx;
  for (int k = 0; k < N/2 + 1; k++) {
    for (int j = 0; j < M/2 + 1; j++) {
      for (int i = 0; i < L/2; i++) {
        Vzx.push_back(2*k*L*M + 2*j*L + (2*i + 1));
      }
    }
  }
  for (int k = 0; k < N/2; k++) {
    for (int j = 0; j < M/2; j++) {
      for (int i = 0; i < L/2 + 1; i++) {
        Vzx.push_back((2*k + 1)*L*M + (2*j + 1)*L + 2*i);
      }
    }
  }
  sort(Vzx.begin(), Vzx.end());
  return Vzx;
}

vector<int> get_Vxy_primal(int L, int M, int N) {
  vector<int> Vxy_primal;
  for (int k = 0; k < N/2 + 1; k++) {
    for (int j = 0; j < M/2 + 1; j++) {
      for (int i = 0; i < L/2 + 1; i++) {
        Vxy_primal.push_back(2*k*L*M + 2*j*L + 2*i);
      }
    }
  }
  sort(Vxy_primal.begin(), Vxy_primal.end());
  return Vxy_primal;
}

vector<int> get_Vxy_dual(int L, int M, int N) {
  vector<int> Vxy_dual;
  for (int k = 0; k < N/2; k++) {
    for (int j = 0; j < M/2; j++) {
      for (int i = 0; i < L/2; i++) {
        Vxy_dual.push_back((2*k + 1)*L*M + (2*j + 1)*L + (2*i + 1));
      }
    }
  }
  sort(Vxy_dual.begin(), Vxy_dual.end());
  return Vxy_dual;
}

vector<int> get_Vyz_primal(int L, int M, int N) {
  vector<int> Vyz_primal;
  for (int k = 0; k < N/2; k++) {
    for (int j = 0; j < M/2 + 1; j++) {
      for (int i = 0; i < L/2; i++) {
        Vyz_primal.push_back((2*k + 1)*L*M + 2*j*L + (2*i + 1));
      }
    }
  }
  sort(Vyz_primal.begin(), Vyz_primal.end());
  return Vyz_primal;
}

vector<int> get_Vyz_dual(int L, int M, int N) {
  vector<int> Vyz_dual;
  for (int k = 0; k < N/2 + 1; k++) {
    for (int j = 0; j < M/2; j++) {
      for (int i = 0; i < L/2 + 1; i++) {
        Vyz_dual.push_back(2*k*L*M + (2*j + 1)*L + 2*i);
      }
    }
  }
  sort(Vyz_dual.begin(), Vyz_dual.end());
  return Vyz_dual;
}

vector<int> get_Vzx_primal(int L, int M, int N) {
  vector<int> Vzx_primal;
  for (int k = 0; k < N/2; k++) {
    for (int j = 0; j < M/2; j++) {
      for (int i = 0; i < L/2 + 1; i++) {
        Vzx_primal.push_back((2*k + 1)*L*M + (2*j + 1)*L + 2*i);
      }
    }
  }
  sort(Vzx_primal.begin(), Vzx_primal.end());
  return Vzx_primal;
}

vector<int> get_Vzx_dual(int L, int M, int N) {
  vector<int> Vzx_dual;
  for (int k = 0; k < N/2 + 1; k++) {
    for (int j = 0; j < M/2 + 1; j++) {
      for (int i = 0; i < L/2; i++) {
        Vzx_dual.push_back(2*k*L*M + 2*j*L + (2*i + 1));
      }
    }
  }
  sort(Vzx_dual.begin(), Vzx_dual.end());
  return Vzx_dual;
}

vector<int> get_Vprimal(int L, int M, int N) {
  vector<int> Vprimal;
  for (int k = 0; k < N/2 + 1; k++) {
    for (int j = 0; j < M/2 + 1; j++) {
      for (int i = 0; i < L/2 + 1; i++) {
        Vprimal.push_back(2*k*L*M + 2*j*L + 2*i);
      }
    }
  }

  for (int k = 0; k < N/2; k++) {
    for (int j = 0; j < M/2 + 1; j++) {
      for (int i = 0; i < L/2; i++) {
        Vprimal.push_back((2*k + 1)*L*M + 2*j*L + (2*i + 1));
      }
    }
  }
  return Vprimal;
}

vector<int> get_faces(int c, int L, int M, int N) {  // returns qubit numbers of qubits on the faces of cube with cube number c;
  vector<int> faces;

  int l = (L + 1) / 2;
  int m = (M + 1) / 2;
  int n = (N - 1) /2;

  int x = (c % (l*m)) % l;
  int y = (c % (l*m)) / l;
  int z = c / (l*m);

  int f = 2*x + 2*y*L + 2*z*L*M;  // qubit number of qubit on the front face (smallest z value) of cube c
  faces.push_back(f);
  faces.push_back(f + 2*L*M);
  if (x != 0) {
    faces.push_back(f - 1 + L*M);
  }
  if (x != l - 1) {
    faces.push_back(f + 1 + L*M);
  }
  if (y != 0) {
    faces.push_back(f - L + L*M);
  }
  if (y != m - 1) {
    faces.push_back(f + L + L*M);
  }
  return faces;
}

void DFS(int v, const vector<vector<int> > &adj, bool visited[], vector<int> &components) {
  visited[v] = true;
  components.push_back(v);

  for (int i = 0; i < adj[v].size(); i++) {
    int w = adj[v][i];
    if (!visited[w]) {
      DFS(w, adj, visited, components);
    }
  }
}



vector<vector<int> > connected_components(const vector<vector<int> > &adj) {
  vector<vector<int> > components_vec;
  int num_vertices = adj.size();
  bool *visited = new bool[adj.size()];

  for (int v = 0; v < num_vertices; v++) {
    visited[v] = false;
  }

  for (int v = 0; v < num_vertices; v++) {
    if (visited[v] == false) {
      vector<int> components;
      DFS(v, adj, visited, components);
      sort(components.begin(), components.end(), greater<int>());
      components_vec.push_back(components);
    }
  }

  return components_vec;
}

int taxicab_distance(int c1, int c2, int l, int m, int n) {  // compute taxicab distance between two cubes, given their cube numbers
  int x1 = (c1 % (l*m)) % l;
  int y1 = (c1 % (l*m)) / l;
  int z1 = c1 / (l*m);
  int x2 = (c2 % (l*m)) % l;
  int y2 = (c2 % (l*m)) / l;
  int z2 = c2 / (l*m);
  int distance = abs(x1 - x2) + abs(y1 - y2) + abs(z1 - z2);
  return distance;
}

int shortest_distance(const vector<int> &chunk1, const vector<int> &chunk2, int l, int m, int n) {
  int d_min = l + m + n;  // distance can't be any larger than this
  for (int i = 0; i < chunk1.size(); i++) {
    int c1 = chunk1[i];
    for (int j = 0; j < chunk2.size(); j++) {
      int c2 = chunk2[j];

      int d = taxicab_distance(c1, c2, l, m, n);
      if (d < d_min) {
        d_min = d;
      }
    }
  }
  return d_min;
}
//
// int main() {
//   vector<int> chunk1;
//   vector<int> chunk2;
//
//   chunk1.push_back(11);
//   chunk1.push_back(15);
//   chunk1.push_back(14);
//   chunk1.push_back(31);
//
//   chunk2.push_back(20);
//
//   cout << shortest_distance(chunk1, chunk2, 4,5,1);
//   cout << min(12, 13);
// }




int decode(int L, int M, int N, const std::vector<int>& Z_vec, double p, std::function<double()> rand) {
  int l = (L + 1) / 2;
  int m = (M + 1) / 2;
  int n = (N - 1) / 2;

  std::vector<int> red_cubes;

  for (int c = 0; c < l*m*n; c++) { // iterate over all cubes, figure out which ones are "red", i.e., X on faces multiply to -1, add these to red_cubes
    int x = (c % (l*m)) % l;
    int y = (c % (l*m)) / l;
    int z = c / (l*m);

    int syndrome = 1;

    int front_face = 2*x + 2*y*L + 2*z*L*M; // qubit number of qubit on the front face (smallest z value) of cube c

    // front_face
    syndrome *= Z_vec[front_face];
    // back face
    syndrome *= Z_vec[front_face + 2*L*M];
    // left face
    if (x != 0) { // cubes at the left (x = 0) boundary are partial and don't have left faces
      syndrome *= Z_vec[front_face - 1 + L*M];
    }
    // right face
    if (x != l - 1) { // cubes at the right boundary don't have right faces
      syndrome *= Z_vec[front_face + 1 + L*M];
    }
    // top face
    if (y != 0) {
      syndrome *= Z_vec[front_face - L + L*M];
    }
    // bottom face
    if (y != m - 1) {
      syndrome *= Z_vec[front_face + L + L*M];
    }

    if (syndrome == -1) {
      red_cubes.push_back(c);
    }
  }

  // print_vec(red_cubes);
  int num_red = red_cubes.size();

  /*
  for each node (i.e., red cube), need to add closest "node" on that's on one of the two boundaries that are parallel to the xy plane
  these should be thought of as extra cubes, with z = -1 or z = n
  */

  /* NOTE: each red cube has a cube number that we'll denote by c (c \in [0, l*m*n - 1]);
  it will also have an index i \in [0, num_red - 1] -- these indices label the nodes in the perfect matching;
  c and i are related by c = red_cubes[i] */



  // for each red cube, compute distance to boundary and determine which plane the closest node is on
  std::vector<int> bound_dist(num_red, 0); // distance of each interior node (red cube) to its closest boundary node
  std::vector<int> which_plane(num_red, 0); // will have which_plane[i] = 0 if the i-th red cube is closer to the Z = 0 plane, and which_plane[i] = 1 if it's closer to the Z = n plane

  for (int i = 0; i < num_red; i++) {
    int c = red_cubes[i];
    int z = c / (l*m);

    int d0 = z + 1;
    int d1 = n - z;

    if (d0 < d1) {
      bound_dist[i] = d0;
    }
    else if (d1 < d0) {
      bound_dist[i] = d1;
      which_plane[i] = 1;
    }
    else {  // if the cube is the same distance from both planes, randomly choose a plane [check that this is okay]
      bound_dist[i] = d0;
      double random_p = rand();
      if (random_p < 0.5) {
        which_plane[i] = 1;
      }
    }
  }

  // for each pair of red cubes, add an edge between them only if the distance between them is no larger than the sum of the distances from each cube to the boundary
  std::vector<int> nodes1, nodes2;
  std::vector<int> int_weights; // weights of edges between red cubes ("interior" nodes) -- note: not all-to-all (see arXiv:0803.0272)

  for (int i = 0; i < num_red; i++) {
    int ci = red_cubes[i];
    int xi = (ci % (l*m)) % l; // I guess storing these in the syndrome for loop might make things faster
    int yi = (ci % (l*m)) / l;
    int zi = ci / (l*m);

    for (int j = 0; j < i; j++) {
      int cj = red_cubes[j];
      int xj = (cj % (l*m)) % l;
      int yj = (cj % (l*m)) / l;
      int zj = cj / (l*m);

      int distance = abs(xi - xj) + abs(yi - yj) + abs(zi - zj);  // taxicab distance

      if (distance <= bound_dist[i] + bound_dist[j]) {
        nodes1.push_back(i);
        nodes2.push_back(j);
        int_weights.push_back(distance);
      }
    }
  }

  int node_num = 2 * num_red;  // # of nodes in perfect matching problem -- factor of 2 since we add in a boundary node for each red cube
  int edge_num = num_red + int_weights.size() + num_red * (num_red - 1) / 2;
  /*
    # of edges between red cubes and boundary nodes = n;
    # of edges between red cubes = int_weights.size();
    # of edges between boundary nodes = (num_red)Choose2
  */

  PerfectMatching *pm = new PerfectMatching(node_num, edge_num);
  /*
    nodes i = 0 to num_red - 1 in PerfectMatching correspond to red cubes; red_cube[i] = cube # in my numbering system
    nodes j = num_red to 2*num_red - 1 correponds to boundary nodes; j corresponds to the closest boundary node to red_cube[j - red_num]
  */

  // add edges between red cubes
  for (int k = 0; k < int_weights.size(); k++) {
    pm->AddEdge(nodes1[k], nodes2[k], int_weights[k]);
  }

  // add edges between each red cube and its closest boundary node
  for (int i = 0; i < num_red; i++) {
    pm->AddEdge(i, i + num_red, bound_dist[i]);
  }

  // add an edge of weight 0 between each pair of boundary nodes
  for (int i = num_red; i < 2*num_red; i++) {
    for (int j = num_red; j < i; j++) {
      pm->AddEdge(i, j, 0);
    }
  }

  struct PerfectMatching::Options options;
  options.verbose = false;
  pm->options = options;
  pm->Solve();

  // determine whether decoding introduced a logical error
  int parity = 1;

  // initial Z errors on front boundary (faces in z = 0 plane)
  for (int c = 0; c < l*m; c++) {
    int x = (c % (l*m)) % l;
    int y = (c % (l*m)) / l;
    int z = c / (l*m);
    parity *= Z_vec[2*x + 2*y*L + 2*z*L*M];
  }

  /* Z operators applied by decoding to the front boundary -- count the number of red cubes that were matched to boundary nodes on one of the two boundaries
  (by checking whether which_plane = 0, we're choosing to check the  Z = 0 boundary) */
  for (int i = 0; i < num_red; i++) {
    int j = pm->GetMatch(i);
    if (j == i + num_red && which_plane[i] == 0) {
      parity *= -1;
    }
  }


  /*
  // the code below is to check both front and back planes -- these should have the same parity!
  int parity2 = 1;

  for (int c = l*m*n - l*m; c < l*m*n; c++) {
    int x = (c % (l*m)) % l;
    int y = (c % (l*m)) / l;
    int z = c / (l*m);
    parity2 *= Z_vec[2*x + 2*y*L + (2*z + 2)*L*M];
  }

  for (int i = 0; i < num_red; i++) {
    int j = pm->GetMatch(i);
    if (j == i + num_red && which_plane[i] == 1) {
      parity2 *= -1;
    }
  }

  std::cout << parity << "  " << parity2 << "\n";
  */
  // std::cout << parity << "   " << (parity + 1) / 2 << "\n";

  delete pm;

  return (parity + 1) / 2;  // returns 1 if no logical error, and 0 if there is a logical error

}

int loss_decode(
    int L, int M, int N, const vector<int> &Z_vec,
    const vector<vector<int> > &superchunks,
    std::function<double()> rand
) {
  int l = (L + 1) / 2;
  int m = (M + 1) / 2;
  int n = (N - 1) / 2;

  vector<int> back_boundary = superchunks[superchunks.size() - 1];

  // first, check whether front boundary vertex and back boundary vertex are in the same connected component -- if so, failure
  if (find(back_boundary.begin(), back_boundary.end(), l*m*n) != back_boundary.end()) {
    // cout << "lost too many qubits \n";
    return 0; // failure
  }
  else {  // otherwise, haven't percolated
    vector<int> front_boundary = superchunks[superchunks.size() - 2];
    vector<int> fbc = vector<int>(front_boundary.begin() + 1, front_boundary.end());  // same vector as front_boundary, except without the element l*m*n
    vector<int> bbc = vector<int>(back_boundary.begin() + 1, back_boundary.end());  // same vector as back_boundary, except without the element l*m*n + 1

    // find syndromes
    vector<vector<int> > red_chunks; // check and supercheck operators with parity -1
    for (int i = 0; i < superchunks.size() - 2; i++) {  // iterate over all check and supercheck operators -- all vectors in superchunks except the last two, which correspond to the boundaries
      vector<int> chunk = superchunks[i];
      int syndrome = 1;
      for (int j = 0; j < chunk.size(); j++) {  // iterate over all cubes in chunk, multiply together measurement results on all faces
        int c = chunk[j];
        vector<int> faces = get_faces(c, L, M, N);
        for (int k = 0; k < faces.size(); k++) {
          syndrome *= Z_vec[faces[k]];
        }
      }

      if (syndrome == -1) {
        red_chunks.push_back(chunk);
      }
    }
    // for (int i = 0; i < red_chunks.size(); i++) {
    //   print_vec(red_chunks[i]);
    // }

    int num_red = red_chunks.size();

    // for each red chunk, compute the shortest distance to each (possibly deformed) boundary, and determine which boundary it's closest to
    vector<int> bound_dist(num_red, 0); // distance of each interior node (red cube) to its closest boundary node (in the MWPM problem)
    vector<int> which_plane(num_red, 0); // will have which_plane[i] = 0 if the i-th red cube is closer to the front boundary, and which_plane[i] = 1 if it's closer to the back boundary
    for (int i = 0; i < num_red; i++) {
      vector<int> chunk = red_chunks[i];

      // find shortest distance between chunk and Z = 0 plane, and between chunk and Z = N - 1 plane
      int z_min0 = n;
      int z_min1 = n;
      for (int j = 0; j < chunk.size(); j++) {
        int c = chunk[j];
        int z = c / (l*m);
        if (z + 1 < z_min0) {
          z_min0 = z + 1;
        }
        if (n - z < z_min1) {
          z_min1 = n - z;
        }
      }

      // find shortest distance between chunk and any cubes deforming the Z = 0 plane (if any -- if none, this just gives d_min0 = l*m*n)
      int d_min0 = shortest_distance(chunk, fbc, l, m, n);
      // find shortest distance between chunk and any cubes deforming the Z = N - 1 plane
      int d_min1 = shortest_distance(chunk, bbc, l, m, n);

      // find shortest distance between chunk and front boundary (front boundary = (Z = 0 plane) + any cubes in fbc)
      int d0 = min(z_min0, d_min0);
      // find shortest distance between chunk and back boundary
      int d1 = min(z_min1, d_min1);

      if (d0 < d1) {
        bound_dist[i] = d0;
      }
      else if (d1 < d0) {
        bound_dist[i] = d1;
        which_plane[i] = 1;
      }
      else {  // if the cube is the same distance from both boundaries, randomly choose one of the boundaries
        bound_dist[i] = d0;
        double rp = rand();
        if (rp < 0.5) {
          which_plane[i] = 1;
        }
      }
    }

    // for each pair of red chunks, add an edge between them only if the distance between them is no larger than the sum of the distances from each chunk to the boundary
    vector<int> nodes1, nodes2;
    vector<int> int_weights;  // weights of edges between red chunks ("interior" nodes) -- note: not all-to-all (see arXiv:0803.0272)

    for (int i = 0; i < num_red; i++) {
      vector<int> chunk1 = red_chunks[i];
      for (int j = 0; j < i; j++) {
        vector<int> chunk2 = red_chunks[j];
        int distance = shortest_distance(chunk1, chunk2, l, m, n);

        if (distance <= bound_dist[i] + bound_dist[j]) {
          nodes1.push_back(i);
          nodes2.push_back(j);
          int_weights.push_back(distance);
        }
      }
    }

    int node_num = 2 * num_red;  // # of nodes in perfect matching problem -- factor of 2 since we add in a boundary node for each red chunk
    int edge_num = num_red + int_weights.size() + num_red * (num_red - 1) / 2;
    /*
      # of edges between red chunks and boundary nodes = n;
      # of edges between red chunks = int_weights.size();
      # of edges between boundary nodes = (num_red)Choose2
    */

    PerfectMatching *pm = new PerfectMatching(node_num, edge_num);
    /*
    nodes i = 0 to num_red - 1 in PerfectMatching correspond to red chunks
    nodes j = num_red to 2*num_red - 1 correponds to boundary nodes; j corresponds to the closest boundary node to red_chunks[j - num_red]
    */

    // add edges between red chunks
    for (int k = 0; k < int_weights.size(); k++) {
      pm->AddEdge(nodes1[k], nodes2[k], int_weights[k]);
    }

    // add edges between each red cube and its closest boundary node
    for (int i = 0; i < num_red; i++) {
      pm->AddEdge(i, i + num_red, bound_dist[i]);
    }

    // add an edge of weight 0 between each pair of boundary nodes
    for (int i = num_red; i < 2*num_red; i++) {
      for (int j = num_red; j < i; j++) {
        pm->AddEdge(i, j, 0);
      }
    }

    struct PerfectMatching::Options options;
    options.verbose = false;
    pm->options = options;
    pm->Solve();

    // determine whether decoding introduced a logical error
    int parity = 1;

    // initial Z errors on front boundary (multiply syndrome on faces in Z = 0 plane and on deforming cubes (in fbc))
    // faces in Z = 0 plane
    for (int c = 0; c < l*m; c++) {
        int x = (c % (l*m)) % l;
        int y = (c % (l*m)) / l;
        int z = c / (l*m);
        parity *= Z_vec[2*x + 2*y*L + 2*z*L*M];
    }
    // faces of deforming cubes
    for (int j = 0; j < fbc.size(); j++) {
      int c = fbc[j];
      vector<int> faces = get_faces(c, L, M, N);
      for (int k = 0; k < faces.size(); k++) {
        parity *= Z_vec[faces[k]];
      }
    }

    /* Z operators applied by decoding to the front boundary -- count the number of red chunks that were matched to boundary nodes on one of the two boundaries
    (by checking whether which_plane = 0, we're choosing to check the  Z = 0 boundary) */
    for (int i = 0; i < num_red; i++) {
      int j = pm->GetMatch(i);
      if (j == i + num_red && which_plane[i] == 0) {
        parity *= -1;
      }
    }

    delete pm;

    return (parity + 1) / 2;//////

  }
}


int weighted_decode(int L, int M, int N, const vector<int> &Z_vec, double p, std::function<double()> rand) {
    // for dephasing per bin, where Vxy qubits have lower total dephasing probability
  int l = (L + 1) / 2;
  int m = (M + 1) / 2;
  int n = (N - 1) / 2;

  vector<int> red_cubes;

  double pL = (1-pow((1-2*p),L))/2;
  double pLM = (1-pow((1-2*p),L*M))/2;
  double w = log(pLM) / log(pL); // weight corresponding to the qubits with higher error probability

  for (int c = 0; c < l*m*n; c++) { // iterate over all cubes, figure out which ones are "red", i.e., X on faces multiply to -1, add these to red_cubes
    int x = (c % (l*m)) % l;
    int y = (c % (l*m)) / l;
    int z = c / (l*m);

    int syndrome = 1;

    int front_face = 2*x + 2*y*L + 2*z*L*M; // qubit number of qubit on the front face (smallest z value) of cube c

    // front_face
    syndrome *= Z_vec[front_face];
    // back face
    syndrome *= Z_vec[front_face + 2*L*M];
    // left face
    if (x != 0) { // cubes at the left (x = 0) boundary are partial and don't have left faces
      syndrome *= Z_vec[front_face - 1 + L*M];
    }
    // right face
    if (x != l - 1) { // cubes at the right boundary don't have right faces
      syndrome *= Z_vec[front_face + 1 + L*M];
    }
    // top face
    if (y != 0) {
      syndrome *= Z_vec[front_face - L + L*M];
    }
    // bottom face
    if (y != m - 1) {
      syndrome *= Z_vec[front_face + L + L*M];
    }

    if (syndrome == -1) {
      red_cubes.push_back(c);
    }
  }

  // print_vec(red_cubes);
  int num_red = red_cubes.size();

  /*
  for each node (i.e., red cube), need to add closest "node" on that's on one of the two boundaries that are parallel to the xy plane
  these should be thought of as extra cubes, with z = -1 or z = n
  */

  /* NOTE: each red cube has a cube number that we'll denote by c (c \in [0, l*m*n - 1]);
  it will also have an index i \in [0, num_red - 1] -- these indices label the nodes in the perfect matching;
  c and i are related by c = red_cubes[i] */



  // for each red cube, compute distance to boundary and determine which boundary (Z = 0 plane or Z = N - 1 plane) the closest node is on
  vector<double> bound_dist(num_red, 0); // distance of each interior node (red cube) to its closest boundary node
  vector<int> which_plane(num_red, 0); // will have which_plane[i] = 0 if the i-th red cube is closer to the Z = 0 plane, and which_plane[i] = 1 if it's closer to the Z = N - 1 plane

  for (int i = 0; i < num_red; i++) {
    int c = red_cubes[i];
    int z = c / (l*m);

    int d0 = z + 1;
    int d1 = n - z;

    if (d0 < d1) {
      bound_dist[i] = double(d0);
    }
    else if (d1 < d0) {
      bound_dist[i] = double(d1);
      which_plane[i] = 1;
    }
    else {  // if the cube is the same distance from both planes, randomly choose a plane [check that this is okay]
      bound_dist[i] = double(d0);
      double random_p = rand();
      if (random_p < 0.5) {
        which_plane[i] = 1;
      }
    }
  }

  // for each pair of red cubes, add an edge between them only if the distance between them is no larger than the sum of the distances from each cube to the boundary
  vector<int> nodes1, nodes2;
  vector<double> int_weights; // weights of edges between red cubes ("interior" nodes) -- note: not all-to-all (see arXiv:0803.0272)

  for (int i = 0; i < num_red; i++) {
    int ci = red_cubes[i];
    int xi = (ci % (l*m)) % l; // I guess storing these in the syndrome for loop might make things faster
    int yi = (ci % (l*m)) / l;
    int zi = ci / (l*m);

    for (int j = 0; j < i; j++) {
      int cj = red_cubes[j];
      int xj = (cj % (l*m)) % l;
      int yj = (cj % (l*m)) / l;
      int zj = cj / (l*m);

      double distance = w * (abs(xi - xj) + abs(yi - yj)) + abs(zi - zj);  // edges along z direction have larger weight since Vxy qubits have lower total error probability

      if (distance <= bound_dist[i] + bound_dist[j]) {
        nodes1.push_back(i);
        nodes2.push_back(j);
        int_weights.push_back(distance);
      }
    }
  }

  int node_num = 2 * num_red;  // # of nodes in perfect matching problem -- factor of 2 since we add in a boundary node for each red cube
  int edge_num = num_red + int_weights.size() + num_red * (num_red - 1) / 2;
  /*
    # of edges between red cubes and boundary nodes = n;
    # of edges between red cubes = int_weights.size();
    # of edges between boundary nodes = (num_red)Choose2
  */

  PerfectMatching *pm = new PerfectMatching(node_num, edge_num);
  /*
    nodes i = 0 to num_red - 1 in PerfectMatching correspond to red cubes; red_cubes[i] = cube # in my numbering system
    nodes j = num_red to 2*num_red - 1 correponds to boundary nodes; j corresponds to the closest boundary node to red_cubes[j - red_num]
  */

  // add edges between red cubes
  for (int k = 0; k < int_weights.size(); k++) {
    pm->AddEdge(nodes1[k], nodes2[k], int_weights[k]);
  }

  // add edges between each red cube and its closest boundary node
  for (int i = 0; i < num_red; i++) {
    pm->AddEdge(i, i + num_red, bound_dist[i]);
  }

  // add an edge of weight 0 between each pair of boundary nodes
  for (int i = num_red; i < 2*num_red; i++) {
    for (int j = num_red; j < i; j++) {
      pm->AddEdge(i, j, 0);
    }
  }

  struct PerfectMatching::Options options;
  options.verbose = false;
  pm->options = options;
  pm->Solve();

  // determine whether decoding introduced a logical error
  int parity = 1;

  // initial Z errors on front boundary (faces in z = 0 plane)
  for (int c = 0; c < l*m; c++) {
    int x = (c % (l*m)) % l;
    int y = (c % (l*m)) / l;
    int z = c / (l*m);
    parity *= Z_vec[2*x + 2*y*L + 2*z*L*M];
  }

  /* Z operators applied by decoding to the front boundary -- count the number of red cubes that were matched to boundary nodes on one of the two boundaries
  (by checking whether which_plane = 0, we're choosing to check the  Z = 0 boundary) */
  for (int i = 0; i < num_red; i++) {
    int j = pm->GetMatch(i);
    if (j == i + num_red && which_plane[i] == 0) {
      parity *= -1;
    }
  }


  /*
  // the code below is to check both front and back planes -- these should have the same parity!
  int parity2 = 1;

  for (int c = l*m*n - l*m; c < l*m*n; c++) {
    int x = (c % (l*m)) % l;
    int y = (c % (l*m)) / l;
    int z = c / (l*m);
    parity2 *= Z_vec[2*x + 2*y*L + (2*z + 2)*L*M];
  }

  for (int i = 0; i < num_red; i++) {
    int j = pm->GetMatch(i);
    if (j == i + num_red && which_plane[i] == 1) {
      parity2 *= -1;
    }
  }

  cout << parity << "  " << parity2 << "\n";
  */
  // cout << parity << "   " << (parity + 1) / 2 << "\n";

  delete pm;

  return (parity + 1) / 2;  // returns 1 if no logical error, and 0 if there is a logical error

}
