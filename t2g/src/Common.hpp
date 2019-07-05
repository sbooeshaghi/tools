#ifndef T2G_COMMON_HPP
#define T2G_COMMON_HPP

#include <cassert>
#include <algorithm>
#include <stdint.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>


#define T2G_VERSION "0.24.0"


struct t2g_opt {
    std::string output;
    std::vector<std::string> files;

    bool stream_in = false;
    bool stream_out = false;
    bool version = false;
};


#endif
