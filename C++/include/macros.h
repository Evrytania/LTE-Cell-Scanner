// Copyright 2012 Evrytania LLC (http://www.evrytania.com)
//
// Written by James Peroulas <james@evrytania.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef HAVE_MACROS_H
#define HAVE_MACROS_H

// Used for debugging. Simply prints an "I am here!" statement to cout.
#define MARK std::cout << "Program execution has reached " << __FILE__ << " line: " << __LINE__ << std::endl

// Improved assert() macro. Takes one or two arguments.
// Ex: ASSERT(n_r==3);
// Ex: ASSERT(n_r==3,"number of rows must be equal to 3!");
#ifndef NDEBUG
#   define ASSERT_CORE(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion (" #condition ") failed in '" << __FILE__ \
                      << "' line " << __LINE__ << ": " << message << std::endl; \
            std::exit(EXIT_FAILURE); \
        } \
    } while (false)
#else
#   define ASSERT_CORE(condition, message) do { } while (false)
#endif
#define ASSERT_1(condition) ASSERT_CORE(condition,"no further information available")
#define ASSERT_2(condition,message) ASSERT_CORE(condition,message)
#define ASSERT_X(x,A,B,FUNC,...) FUNC
#define ASSERT(...) ASSERT_X(,##__VA_ARGS__,ASSERT_2(__VA_ARGS__),ASSERT_1(__VA_ARGS__))

// max/min
#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)<(B)?(A):(B))

// wrap/rail
#define WRAP(x,sm,lg) (itpp_ext::matlab_mod((x)-(sm),(lg)-(sm))+(sm))
#define RAIL(x,min,max) (((x)<(min))?(min):(((x)>(max))?(max):(x)))

// Is this defined in ITPP somewhere?
#define J (complex<double>(0,1))

// Used to easily save variables in the itpp debug file.
// Ex: ITPP_DEBUG_EXPORT(myvar);
//   will save the variable myvar in the ITPP_DEBUG file as ITPP_myvar.
// Ex: ITPP_DEBUG_EXPORT(myvar,newname);
//   will save the variable myvar in the ITPP_DEBUG file as ITPP_newname.
#define ITPP_DEBUG_EXPORT_CORE(CPP_NAME,ITPP_NAME) \
do { \
  ITPP_DEBUG << Name("ITPP_" #ITPP_NAME) << CPP_NAME; \
  ITPP_DEBUG.flush(); \
} while (false)
#define ITPP_DEBUG_EXPORT_1(CPP_NAME) ITPP_DEBUG_EXPORT_CORE(CPP_NAME,CPP_NAME)
#define ITPP_DEBUG_EXPORT_2(CPP_NAME,ITPP_NAME) ITPP_DEBUG_EXPORT_CORE(CPP_NAME,ITPP_NAME)
#define ITPP_DEBUG_EXPORT_X(x,A,B,FUNC,...) FUNC
#define ITPP_DEBUG_EXPORT(...) ITPP_DEBUG_EXPORT_X(,##__VA_ARGS__,ITPP_DEBUG_EXPORT_2(__VA_ARGS__),ITPP_DEBUG_EXPORT_1(__VA_ARGS__))

#ifndef NDEBUG
extern itpp::it_file ITPP_DEBUG;
#endif

#endif

