#include <iostream>
#include "graph_data_structures.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>
#include <algorithm>
#include <fstream>
#include <cstdio>

int getdir (const std::string& dir, std::vector<std::string> &files)
{
    using namespace std;
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
      struct stat st;
      std::string fullname = dir+string(dirp->d_name);
      stat(fullname.c_str(), &st);
      if (!S_ISDIR(st.st_mode))
      {
          files.push_back(std::move(fullname));
      }
    }
    closedir(dp);
    return 0;
}

timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

void print_help()
{
  std::cout << 
  "simsimhash: Utility for running the algortihm described in Papadimitiou et. al" << std::endl <<
  "Usage" << std::endl << std::endl <<
  "  simsimhash <folder> <hash-difference-file> <anomaly_file>" << std::endl;
};

struct p4_dir_comp
{
  bool operator() (const std::string& lhs, const std::string& rhs)
  {
    int i=0;
    int j=0;
    sscanf(lhs.c_str(),"%i_",&i);
    sscanf(rhs.c_str(),"%i_",&j);
    return i < j;
  }
};
int main(int argc,const char** argv)
{
  using namespace csc715;
  
  if (argc < 4)
  {
    std::cerr << "Error, you must have 4 arguments. See the following usage" << std::endl;
    print_help();
    return -1;
  }
  std::string infolder = argv[1];
  std::string hash_values = argv[2];
  std::string anomalys = argv[3];
  
  if (infolder.back() != '/')
  {
    infolder+='/';
  }
  std::vector<std::string> files;
  getdir(infolder,files);
  p4_dir_comp comp;
  std::sort(files.begin(), files.end(), comp);
  
  uint64_t last_data;
  struct timespec start,stop;
  std::vector<double> diffs;

  for (int i = 0; i < files.size(); i++)
  {

    const std::string& infile = files[i];  
    
    
    //std::cout << "Reading File " << infile << "..." << std::endl; 
    EdgeList el = read_edge_list(infile.c_str());

    clock_gettime(CLOCK_REALTIME,&start);
    //std::cout << "Converting to Adjacency List..." << std::endl;
    AdjacencyList al = edge_list2adjacency_list(el);

    //std::cout << "Calculating page feautres..." << std::endl;
    auto features = generateFeatures(al);
    uint64_t curr_hash = simhash(features);
    

    double d = 0;
    if (i > 0)
    {
      d=1.0-(hamming_distance(last_data,curr_hash)/64.0);
      diffs.push_back(d);
    }
    clock_gettime(CLOCK_REALTIME,&stop);
    auto td = diff(start,stop);
    double ms = 1000.0*td.tv_sec+1.0e-6*td.tv_nsec;
    std::cout << i << ':' << std::hex << curr_hash << std::dec;
    std::cout << '=' << d << '@' << ms;
    std::cout << std::endl;
    last_data=curr_hash;
  }
  
  //compute median
  double median;
  
  {
    std::vector<double> diffs2(diffs);
    size_t medposn = diffs2.size()/2;
    std::nth_element(diffs2.begin(),diffs2.begin()+medposn,diffs2.end());
    median = diffs2[medposn];
  }
  
  //calculate   MR
  double MR = 0;
  for (int i = 1; i < diffs.size(); i++)
  {
    MR+=std::abs(diffs[i]-diffs[i-1]);
  }
  MR /= (diffs.size()-1);
  
  double lower_bound = median-3*MR;
  std::cout << "Lower Bound: " << lower_bound;
  
  //find anomalous points
  std::vector<size_t> anomolous_points;
  anomolous_points.reserve(diffs.size());
  bool prev_anomolouse = false;
  for (int i = 0; i < diffs.size(); i++)
  {
    if (diffs[i]>lower_bound)
    {
       prev_anomolouse = false;
    }
    else if (!prev_anomolouse)
    {
      prev_anomolouse = true;
      anomolous_points.push_back(i+1);
    }
  }

  std::cout << "," << anomolous_points.size() << " found" << std::endl;
  for (int i = 0; i < anomolous_points.size(); i++)
  {
    std::cout << anomolous_points[i] << std::endl;
  }
  
  std::cout << "Writing Files..." << std::endl;
  
  std::ofstream ofs(hash_values,std::ios_base::out | std::ios_base::trunc);
  ofs<<lower_bound << std::endl;
  
  for (int i = 0; i < diffs.size(); i++)
  {
    ofs << diffs[i] << std::endl;
  }
  ofs.close();
  ofs.open(anomalys,std::ios_base::out | std::ios_base::trunc);
  ofs << anomolous_points.size() << std::endl;
  
  for (int i = 0; i < anomolous_points.size(); i++)
  {
    ofs << anomolous_points[i] << std::endl;
  }
  ofs.close();
  return 0;
}
