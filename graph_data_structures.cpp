#include "graph_data_structures.hpp"
#include <fstream>
#include <stdlib.h>
#include <openssl/sha.h>
#include <openssl/md5.h>
#include <iostream>

namespace csc715
{
  std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) 
  {
      std::stringstream ss(s);
      std::string item;
      while (std::getline(ss, item, delim)) {
          elems.push_back(item);
      }
      return elems;
  }

  const double* get_data(const Eigen::MatrixXd& xx) {return xx.data();}
  std::vector<std::string> split(const std::string &s, char delim)
  {
      std::vector<std::string> elems;
      split(s, delim, elems);
      return elems;
  }

  EdgeList read_edge_list(const char* filename)
  {
    std::ifstream file(filename);
    
    std::string line;
    EdgeList ret;
    
    while (getline(file,line))
    {
      auto splitted = split(line,' ');
      char* end;
      ret.push_back(std::make_tuple(std::strtoull(splitted[0].c_str(),&end,0),
        std::strtoull(splitted[1].c_str(),&end,0)));
    }
    return ret;
  }
  AdjacencyList edge_list2adjacency_list(const EdgeList& e)
  {
    //find max
    node_tp max_node = 0;
    for (auto i = e.cbegin(); i!=e.cend(); ++i)
    {
      max_node = std::max(max_node,std::max(std::get<0>(*i),std::get<1>(*i)));
    }
    
    AdjacencyList ret;
    ret.resize(max_node+1);
    
    for (auto i = e.cbegin(); i!=e.cend(); ++i)
    {
      ret[std::get<0>(*i)].push_back(std::get<1>(*i));
    }
    return ret;
  }
  
  std::vector<double> page_rank(const AdjacencyList& al, double d, int max_iter)
  {
    const size_t n = al.size();
    std::vector<double> out_degree(n);
    
    for (int i = 0; i < n; i++)
    {
      out_degree[i] = al[i].size();
    }
    
    //remove dangling nodes, we'll compute them later
    size_t gn = 0;
    std::vector<size_t> oindex(n);
    size_t ii = 0;
    std::vector<double> gout_degree; gout_degree.reserve(n);
    
    AdjacencyList gE;

    for (size_t i = 0; i < n; i++)
    {
       const size_t outdeg = out_degree[i];
       if (outdeg)
       {
         gout_degree.push_back(outdeg);
         
         std::vector<node_tp> ls;
         for (auto j = al[i].begin(); j != al[i].end(); ++j)
         {
           if (out_degree[*j])
           {
             ls.push_back(*j);
           }
         }
         
         gE.push_back(std::move(ls));
         gn++;
       }
       oindex[i] = gn-1;
    };
    
    for (size_t i = 0; i < gn; i++)
    {
      for (size_t j = 0; j < gE[i].size(); j++)
      {
        gE[i][j] = oindex[gE[i][j]];
      }
    }
    
    auto gW = Eigen::SparseMatrix<double>(adjacency_list_to_adjacency_matrix(gE).transpose());
    //scale all the columns
    for (size_t i = 0; i < gn; i++)
    {
      gW.col(i) /= gout_degree[i];
    }
    
    Eigen::MatrixXd xkv_1 = Eigen::MatrixXd::Constant(gn,1,1.0/(double)n);
    Eigen::MatrixXd xkv_2 = Eigen::MatrixXd(gn,1);
    
    Eigen::MatrixXd& xk = xkv_1;
    Eigen::MatrixXd& xk_1 = xkv_2;
    
    double diff;
    size_t k = 0;
    do
    {
      k++;
      xk_1 = (1-d)+Eigen::ArrayXXd(d*gW*xk);
      diff = (xk-xk_1).norm();
      std::swap(xk,xk_1);
    }
    while(diff > 1e-6 && (max_iter < 0 || k < max_iter));
    
    std::vector<double> ret(n);
    auto E_reverse = reverse_adjacency_list(al);
    
    for (int i = 0; i < n; i++)
    {
      if (out_degree[i])
      {
        ret[i] = xk(oindex[i]);
      }
      else
      {
        double sum = 0;
        
        for (auto j = E_reverse[i].cbegin();
          j!=E_reverse[i].cend();
          ++j)
        {
          sum += xk(oindex[*j])/static_cast<double>(out_degree[*j]);
        }
        ret[i] = (1-d)+d*sum;
      }
    }
    
    return ret;
  }
  
    AdjacencyMatrix adjacency_list_to_adjacency_matrix(const AdjacencyList& E)
    {
       const size_t n = E.size();
       AdjacencyMatrix ret = AdjacencyMatrix(n,n);
       
       for (int i = 0; i < n; i++)
       {
          for (auto j = E[i].cbegin(); j != E[i].cend(); j++)
          {
            ret.insert(i,*j)=1;
          }
       }
       
       return ret;
    }
    
    AdjacencyList reverse_adjacency_list(const AdjacencyList& al)
    {
       AdjacencyList ret(al.size());
       
       for (node_tp i = 0; i<al.size(); ++i)
       {
         for (auto j = al[i].cbegin(); j != al[i].end(); ++j)
         {
           ret[*j].push_back(i);
         }
       }
       
       return ret;
    }
    

    #define HASH_SIZE 16
   template <typename T>
    HashKey get_hash(const T& val)
    {
      uint8_t full_hash[HASH_SIZE];
      HashKey hash_key;
      MD5(reinterpret_cast<const uint8_t*>(&val), sizeof(val), full_hash);

//      SHA256(reinterpret_cast<const uint8_t*>(&val), sizeof(val), full_hash);
      for (int i = 0; i < HASH_KEY_BYTES; i++)
      {
         hash_key.val[i] = full_hash[i];
      }
      //SHA1(reinterpret_cast<const uint8_t*>(&val), sizeof(val), hash_key.val);
      //*reinterpret_cast<uint64_t*>(HashKey.val) = vakl
      return hash_key;
    }
    
    PageFeatures generateFeatures(const AdjacencyList& al)
    {
      auto pr = page_rank(al,0.85);
      //count the number of edges in list
      size_t num_keys = 0;
      for (int i = 0; i < al.size(); i++)
      {
        num_keys+=al[i].size();
        if (al[i].size()) num_keys++;
      }
      PageFeatures ret;
      ret.reserve(num_keys);
      for (uint32_t i = 0; i < al.size(); i++)
      {
        size_t n = al[i].size();
        if (n)
        {
          uint64_t key_start = static_cast<uint64_t>(i)<<static_cast<uint64_t>(32);
          ret.emplace_back(get_hash(key_start|static_cast<uint64_t>(i)),pr[i]);
          
          for (auto j = al[i].cbegin(); j!=al[i].cend();++j)
          {
            ret.emplace_back(get_hash(key_start|static_cast<uint64_t>(*j)),pr[i]/static_cast<double>(n));
          }
        }
      }
      return ret;
    }
    

    HashKey simhash(const PageFeatures& pf)
    {
      const uint64_t num_bits = 8*HASH_KEY_BYTES;
      std::vector<double> sums(num_bits);
      
      for (auto iter = pf.cbegin(); iter != pf.cend(); ++iter)
      {
        HashKey key = std::get<0>(*iter);
        double value = std::get<1>(*iter);
        
        for (uint64_t j = 0; j < num_bits; j++)
        {
          uint64_t byte_num = (j/8);
          uint64_t bit_num = (j%8);
          int passed = (key.val[byte_num] & (1<<bit_num))>0;
          int multiplyer = 2*passed-1;
          sums[j]+=(multiplyer*value);
        }
      }
      
      HashKey ret;
      memset(ret.val,0,HASH_KEY_BYTES);
      
      for (int i = 0; i < num_bits; i++)
      {
        uint64_t byte_num = (i/8);
        uint64_t bit_num = (i%8);
        ret.val[byte_num] |= (sums[i]>0)<<bit_num;
      }
      
      return ret;
    }
    size_t hamming_distance(const HashKey& x, const HashKey& y)
    {
      int sum = 0;
      int lim = HASH_KEY_BYTES;
      for (int i = 0; i < lim; i++)
      {
         sum+=__builtin_popcountll(x.val[i]^y.val[i]);
      }
      return sum;
    }
}
