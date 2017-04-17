#ifndef PathHelpers_hpp
#define PathHelpers_hpp

#include <string>

class PathHelpers
{
public:
    static std::string getParentDirectory (std::string pathName);
    static void ensureDirectoryExists (std::string pathName);
    static void ensureParentDirectoryExists (std::string dirName);
};


#endif
