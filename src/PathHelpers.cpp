
#include <sys/stat.h>
#include "PathHelpers.hpp"


std::string PathHelpers::getParentDirectory (std::string pathName)
{
	std::string::size_type lastSlash = pathName.find_last_of ("/");
    return pathName.substr (0, lastSlash);
}

void PathHelpers::ensureDirectoryExists (std::string dirName)
{
    mkdir (dirName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

void PathHelpers::ensureParentDirectoryExists (std::string pathName)
{
	std::string parentDir = getParentDirectory (pathName);
    mkdir (parentDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}
