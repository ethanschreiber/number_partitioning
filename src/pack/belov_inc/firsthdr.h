#ifdef  __FIRSTHDR_H
#error should be included only once
#endif
#define __FIRSTHDR_H


/* Author: Gleb Belov <belov@web.de> */


#define COMMON_NAMESPACE__ bgn
#define BEGIN_COMMON_NAMESPACE__ namespace bgn {
#define END_COMMON_NAMESPACE__   }


BEGIN_COMMON_NAMESPACE__
//using namespace std;
END_COMMON_NAMESPACE__


//#define DBG_ON


#include "mydebug.h"
#include "mytools.h"
#include "myopt.h"
#include "mystat.h"


//#include "spxdefines.h"
#include "timer.h"
#include "random.h"
BEGIN_COMMON_NAMESPACE__
typedef soplex::Timer Timer;
typedef soplex::Random Random;
END_COMMON_NAMESPACE__


#include "sshdrs.h"
#include "mymath.h"


