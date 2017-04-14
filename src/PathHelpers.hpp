#ifndef PathHelpers_hpp
#define PathHelpers_hpp

#include <string>

class PathHelpers
{
public:
	static std::string getParentDirectory (std::string filePath);
	static void ensureParentDirectoryExists (std::string filePath);
};


#endif
