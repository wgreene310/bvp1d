#ifndef _mexDebug_h_
#define _mexDebug_h_

#ifdef WIN32

#include <windows.h>

#include <mex.h>

namespace {
  void checkHeap(const char *msg) {
    HANDLE heapHndl = GetProcessHeap();
    BOOL heapOK = HeapValidate(heapHndl, 0, 0);
    if (heapOK)
      mexPrintf("### %s: Heap OK.\n", msg);
    else
      mexPrintf("### %s: Heap Corrupted.\n", msg);

  }
}
#endif

#endif