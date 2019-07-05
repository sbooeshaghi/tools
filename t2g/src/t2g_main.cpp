#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>
#include <unistd.h> 

#include <iostream>
#include <fstream>
#include <thread>
#include <algorithm>
// #include <atomic>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <sstream>

#include "Common.hpp"


std::vector<std::string> split(const std::string& str, const std::string& delim)
{
    std::vector<std::string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == std::string::npos) pos = str.length();
        std::string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

bool checkFileExists(const std::string &fn) {
    struct stat stFileInfo;
    auto intStat = stat(fn.c_str(), &stFileInfo);
    return intStat == 0;
}

void parse_ProgramOptions(int argc, char **argv, t2g_opt& opt) {

    const char* opt_string = "o:p";
    int version_flag = 0;
    static struct option long_options[] = {
        {"output",          required_argument,  0, 'o'},
        {"pipe",            no_argument, 0, 'p'},
        {"version",         no_argument, &version_flag, 1},
        {0,                 0,                  0,  0 }
    };

    int option_index = 0, c;

    while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

        switch (c) {
            case 'o':
                opt.output = optarg;
                break;
            case 'p':
                opt.stream_out = true;
                break;
            default:
                break;
        }
    }
    if (version_flag) {
        opt.version = true;
    }

    // all other arguments are fast[a/q] files to be read
    while (optind < argc) opt.files.push_back(argv[optind++]);

    if (opt.files.size() == 1 && opt.files[0] == "-") {
        opt.stream_in = true;
    }
}

bool check_ProgramOptions(t2g_opt& opt) {
    bool ret = true;

  if (opt.output.empty() && !opt.stream_in) {
    std::cerr << "Error: Missing output file" << std::endl;
    ret = false;
  } 


  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing GTF input files" << std::endl;
    ret = false;
  } else {
    if (!opt.stream_in) {
      for (const auto& it : opt.files) {  
        if (!checkFileExists(it)) {
          std::cerr << "Error: File not found, " << it << std::endl;
          ret = false;
        }
      }
    }
  }
    return ret;
}


void t2g_Usage() {
    std::cout << "Usage: t2g [options] GTF-file" << std::endl << std::endl
        << "Options: " << std::endl
        << "-o, --output        File for t2g output" << std::endl
        << "-p, --pipe          Write to standard output" << std::endl
        << "--version           Use version numbers" << std::endl
        << std::endl;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        t2g_Usage();
        exit(1);
    } else {
        bool disp_help = argc == 2;
        t2g_opt opt;
        std::string cmd(argv[1]);
        if (cmd == "make") {
            if (disp_help) {
                t2g_Usage();
                exit(0);
            }
            parse_ProgramOptions(argc-1, argv+1, opt);
            if (check_ProgramOptions(opt)) { // program options valid 
                // t2g(opt);
                // std::cout  << "Everything works." << std::endl;
                std::string line;
                std::streambuf *buf = nullptr;
                std::ofstream of;

                if (!opt.stream_out) {
                    of.open(opt.output);
                    buf = of.rdbuf();
                } else {
                    buf = std::cout.rdbuf();
                }
                std::ostream o(buf);

                std::string tmp;
                std::vector<std::string> token;
                char delim_line = '\t';
                std::string delim_attr = " ";
                std::string gid, gid_version, tid, tid_version, gene;
                while(std::getline(std::cin, line)) {
                    if (line.length() == 0 || line[0] == '#') {
                        continue;
                    }
                    std::stringstream ss(line);
                    while (std::getline(ss, tmp, delim_line)){
                        token.push_back(tmp);
                    }
                    if (token[2] == "transcript") {
                        token[8].erase(std::remove(token[8].begin(), token[8].end(), ';'), token[8].end());
                        std::vector<std::string> v = split(token[8], delim_attr);
                        for (int i=0; i<v.size();i++) {
                            if (v[i]=="gene_id"){
                                gid = v[i+1];
                                gid.erase(std::remove(gid.begin(), gid.end(), '"'), gid.end());
                            } else if (v[i] == "gene_version"){
                                gid_version = v[i+1];
                                gid_version.erase(std::remove(gid_version.begin(), gid_version.end(), '"'), gid_version.end());

                            } else if (v[i] == "transcript_id"){
                                tid = v[i+1];
                                tid.erase(std::remove(tid.begin(), tid.end(), '"'), tid.end());
                            } else if (v[i] == "transcript_version"){
                                tid_version = v[i+1];
                                tid_version.erase(std::remove(tid_version.begin(), tid_version.end(), '"'), tid_version.end());

                            } else if (v[i] == "gene_name"){
                                gene = v[i+1];
                                gene.erase(std::remove(gene.begin(), gene.end(), '"'), gene.end());
                            }
                        }
                        if (!opt.version) {
                            o << tid << "\t" << gid << "\t" << gene << "\n";
                        } else {
                            o << tid << "." << tid_version << "\t" << gid << "." << gid_version << "\t" << gene << "\n";
                        }
                        token.clear();
                    } else {
                        token.clear();
                    }
                }

            } else {
                t2g_Usage();
                exit(1);
            }
        } else {
            t2g_Usage();
            exit(0);
        }
    }
}
