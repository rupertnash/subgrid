#ifndef XArray_H
#define XArray_H

#define XArray_ERROR 255

#include <Python.h>

typedef struct XArrayCtx {
  PyObject* self;
  size_t nKeys, nCols, nRows;
  char* table;
  double* data;
} XArrayCtx;

typedef struct XArrayMember {
  XArrayCtx* ctx;
  size_t start;
} XArrayMember;

XArrayCtx* XArray_initCtx(PyObject* self);
void XArray_delCtx(XArrayCtx* ctx);

int XArray_hasKey(XArrayCtx* ctx, char k);
size_t XArray_getStart(XArrayCtx* ctx, char k);
size_t XArray_getLen(XArrayCtx* ctx, char k);

void XArray_getMember(XArrayCtx* ctx, char k, XArrayMember* mem);

double* XArray_getItem(XArrayMember* mem, size_t i);

#endif
