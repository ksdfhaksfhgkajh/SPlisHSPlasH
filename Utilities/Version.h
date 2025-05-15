#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "36648620eb68e2c529ff52557e369948ebf5a9b2"
#define GIT_REFSPEC "refs/heads/viscosity_predict"
#define GIT_LOCAL_STATUS "CLEAN"

#define SPLISHSPLASH_VERSION "2.13.1"

#ifdef DL_OUTPUT

#endif

#endif
