#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "84bced16de1faa31579125b0e736b138a7e59c9b"
#define GIT_REFSPEC "refs/heads/particle_project"
#define GIT_LOCAL_STATUS "CLEAN"

#define SPLISHSPLASH_VERSION "2.13.1"

#ifdef DL_OUTPUT

#endif

#endif
