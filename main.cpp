#include <iostream>
#include "graph_data_structures.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>

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

int main(int argc,const char** argv)
{
  using namespace csc715;
  
  std::string infolder;
  if (argc < 2)
  {
    infolder = "/home/njclimer/data/p4/as-733";
  }
  else
  {
    infolder = argv[1];
  }
  
  if (infolder.back() != '/')
  {
    infolder+='/';
  }
  std::vector<std::string> files;
  getdir(infolder,files);
  
  uint64_t last_data;
  struct timespec start,stop;
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
    }
    clock_gettime(CLOCK_REALTIME,&stop);
    auto td = diff(start,stop);
    double ms = 1000.0*td.tv_sec+1.0e-6*td.tv_nsec;
    std::cout << i << ':' << std::hex << curr_hash << std::dec;
    std::cout << '=' << d << '@' << ms;
    std::cout << std::endl;
    last_data=curr_hash;
  }
}
