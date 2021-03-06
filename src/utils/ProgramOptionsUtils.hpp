/*
 * ProgramOptionsUtils.hpp
 *
 *  Created on: Jan 17, 2013
 *      Author: ethan
 */

#ifndef PROGRAMOPTIONSUTILS_HPP_
#define PROGRAMOPTIONSUTILS_HPP_

#include "../globals.hpp"

// =================================================
// For Boost Options since their required function
// is not supported by older versions on brimstone
// =================================================
template <class T>
void optionRequired (string paramName, T value, T unset) {
 if (value == unset) {
   std::ostringstream out;
   out << paramName << " is required.";
   throw std::runtime_error (out.str() );
 }
}
#endif /* PROGRAMOPTIONSUTILS_HPP_ */
