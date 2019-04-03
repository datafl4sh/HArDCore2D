#ifndef UTILS_HPP
#define UTILS_HPP
#include <iostream>
using namespace HArDCore2D;
// macro for debugging, makes python like print function
#ifndef DEBUG
#ifndef print
#define print(message) std::cout << message << std::endl;
#endif
#else
#ifndef print
#define print(message)
#endif
#endif
#endif /* UTILS_HPP */
