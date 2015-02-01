#ifndef CHEMREAC_UTIL_H__AZQIEGXT55HCHGTXM3T47G6BRU
#define CHEMREAC_UTIL_H__AZQIEGXT55HCHGTXM3T47G6BRU

#include <cstddef>
#include <string>
#include <iostream>
#include <fstream>

namespace chemreac {

    template <typename T>
    inline bool save_array(const T * const data, std::size_t length, const std::string& path)
    {
        std::ofstream ofs(path.c_str(), std::ios::out | std::ios::binary);
        if (!ofs.is_open())
            return false;
        ofs.write(reinterpret_cast<const char * const>(data),
                  std::streamsize(length*sizeof(T)));
        return true;
    }

    template <typename T>
    inline bool load_array(T * const data, std::size_t length, const std::string& path)
    {
        std::ifstream ifs(path.c_str(), std::ios::in | std::ios::binary);
        if (!ifs.is_open())
            return false;
        ifs.read(reinterpret_cast<char * const>(data),
                 std::streamsize(length*sizeof(T)));
        return true;
    }

}
#endif /* CHEMREAC_UTIL_H__AZQIEGXT55HCHGTXM3T47G6BRU */
