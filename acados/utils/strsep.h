/* strsep.h

  Provides the 4.4BSD strsep(3) function for those that don't have it.

  Copyright 2011 Michael Thomas Greer
  Distributed under the Boost Software License, Version 1.0.
  ( See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt )

  Including this file modifies the std namespace in C++.

  Don't include this file if your compiler provides the strsep function in <string.h>.
  Make sure your build process tests for this and behaves accordingly!

*/
#ifdef __cplusplus
  #include <cstring>
  #define STD(x) std::x
  namespace std {
#else
  #include <string.h>
  #define STD(x) x
#endif

char* strsep( char** stringp, const char* delim )
  {
  char* result;

  if ((stringp == NULL) || (*stringp == NULL)) return NULL;

  result = *stringp;

  while (**stringp && !STD(strchr)( delim, **stringp )) ++*stringp;

  if (**stringp) *(*stringp)++ = '\0';
  else             *stringp    = NULL;

  return result;
  }

#ifdef __cplusplus
  } // namespace std;
#endif

#undef STD 