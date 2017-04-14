
#include <sys/stat.h>
#include "PathHelpers.hpp"


std::string PathHelpers::getParentDirectory (std::string filePath)
{
	std::size_t lastSlash = filePath.find_last_of ("/");
    return filePath.substr (0, lastSlash);
}

void PathHelpers::ensureParentDirectoryExists (std::string pathName)
{
	std::string parentDir = getParentDirectory (pathName);
    mkdir (parentDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}
