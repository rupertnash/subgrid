#include "XArray.h"
#define NO_IMPORT_ARRAY
#include "DQarrayobject.h"

/* cTable has 4 columns: key, startPos, length, reserved */
#define XARRAY_TABLE_COLS 4

typedef struct XArrayCtxFull {
  PyObject* self;
  size_t nKeys, nCols, nRows;
  char* table;
  double* data;
  PyArrayObject* tableObj;
  PyArrayObject* dataObj;
} XArrayCtxFull;



XArrayCtx* XArray_initCtx(PyObject* self) {
  XArrayCtxFull* ctx = (XArrayCtxFull*)malloc(sizeof(XArrayCtxFull));
  /* we haven't so extract some C values from the Python structures -
     note these are NEW REFERENCES - hang onto them for the lifetime
     of the context to ensure the PyObject isn't killed.
  */
  ctx->self = self;
  /* get the lookup table for columns */
  ctx->tableObj = (PyArrayObject *)PyObject_GetAttrString(self, "cTable");
  /* Number of keys */
  ctx->nKeys = PyArray_DIM(ctx->tableObj, 0);
  /* Raw pointer to table */
  ctx->table = PyArray_DATA(ctx->tableObj);
  
  /* Now get the actual data numpy array */
  ctx->dataObj = (PyArrayObject *)PyObject_GetAttrString(self, "data");
  /* Store its shape */
  ctx->nRows = PyArray_DIM(ctx->dataObj, 0);
  ctx->nCols = PyArray_DIM(ctx->dataObj, 1);
  /* Raw pointer to data */
  ctx->data = PyArray_DATA(ctx->dataObj);
  
  return (XArrayCtx*)ctx;
}
void XArray_delCtx(XArrayCtx* ctx) {
  XArrayCtxFull* fullCtx = (XArrayCtxFull*)ctx;
  Py_DECREF(fullCtx->dataObj);
  Py_DECREF(fullCtx->tableObj);
  free(fullCtx);
}


int XArray_hasKey(XArrayCtx* ctx, char k) {
  int i;
  for (i=0; i<ctx->nKeys; i++) {
    if (ctx->table[i*XARRAY_TABLE_COLS + 0] == k) {
      /* found the right row */
      return 1;
    }
  }
  return 0;
}

size_t XArray_getStart(XArrayCtx* ctx, char k) {
  int i;
  
  for (i=0; i<ctx->nKeys; i++) {
    if (ctx->table[i*XARRAY_TABLE_COLS + 0] == k) {
      /* found the right row */
      return ctx->table[i*XARRAY_TABLE_COLS + 1];
    }
  }
  PyErr_Format(PyExc_AttributeError, "XArray has no attribute '%c'", k);
  return XArray_ERROR;
}

size_t XArray_getLen(XArrayCtx* ctx, char k) {
  int i;
  
  for (i=0; i<ctx->nKeys; i++) {
    if (ctx->table[i*XARRAY_TABLE_COLS + 0] == k) {
      /* found the right row */
      return ctx->table[i*XARRAY_TABLE_COLS + 1];
    }
  }
  PyErr_Format(PyExc_AttributeError, "XArray has no attribute '%c'", k);
  return XArray_ERROR;
}

void XArray_getMember(XArrayCtx* ctx, char k, XArrayMember* mem) {
  mem->ctx = ctx;
  mem->start = XArray_getStart(ctx, k);
}

double* XArray_getItem(XArrayMember* mem, size_t i) {
  return &(mem->ctx->data[mem->ctx->nCols*i + mem->start]);
}
