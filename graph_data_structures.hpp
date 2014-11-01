#ifndef GRAPH_DATA_STRUCTURES_HPP
#define GRAPH_DATA_STRUCTURES_HPP
#include <vector>
#include <cstdint>
#include <utility>
#include <Eigen/SparseCore>

namespace csc715
{
  #define HASH_KEY_BYTES 16
  #define HASH_KEY_BITS (HASH_KEY_BYTES*8)
  struct HashKey
  {
    uint8_t val[HASH_KEY_BYTES];
  };  
  
  typedef uint32_t node_tp;
  typedef Eigen::SparseMatrix<double> AdjacencyMatrix;
  typedef std::vector<std::tuple<node_tp,node_tp>> EdgeList;
  typedef std::vector<std::vector<node_tp>> AdjacencyList;

  EdgeList read_edge_list(const char* filename);
  
  AdjacencyList edge_list2adjacency_list(const EdgeList& e);
  
  std::vector<double> page_rank(const AdjacencyList& al, const double d = 1.0,
  const int max_iter = -1);
  
  AdjacencyList reverse_adjacency_list(const AdjacencyList& al);
  
  AdjacencyMatrix adjacency_list_to_adjacency_matrix(const AdjacencyList& E);
  
  typedef std::vector<std::tuple<HashKey,double>> PageFeatures;
  
  PageFeatures generateFeatures(const AdjacencyList& al);


  HashKey simhash(const PageFeatures& pf);
  size_t hamming_distance(const HashKey& x, const HashKey& y);

}

#endif
