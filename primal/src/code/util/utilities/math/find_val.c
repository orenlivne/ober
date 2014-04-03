#include "Python.h"
#include "Numeric/arrayobject.h"

#define SEARCH_LOOP(OP, INC) { \
for (i = start_addr; i != end_addr; INC) { \
    if (*i OP value) { \
        return (void*)i; \
    } \
} \
return (void*)0; }

#define SWITCH_ON_OP(TYPE, INC) { \
TYPE* i; \
switch (op) { \
    case EQ: \
        SEARCH_LOOP(==, INC) \
    case GT: \
        SEARCH_LOOP(>, INC) \
    case GTE: \
        SEARCH_LOOP(>=, INC) \
    case LT: \
        SEARCH_LOOP(<, INC) \
    case LTE: \
        SEARCH_LOOP(<=, INC) \
    case NEQ: \
        SEARCH_LOOP(!=, INC) \
    default: \
        return (void*)0; \
    } \
}

#define SWITCH_ON_STRIDE(TYPE) { \
    if (stride > 0) \
    { \
        SWITCH_ON_OP(TYPE, i++) \
    } \
    else \
    { \
        SWITCH_ON_OP(TYPE, i--) \
    } \
}

enum comp_operators {
    EQ,
    GT,
    GTE,
    LT,
    LTE,
    NEQ
};

static void* array_find_byte(PyObject* searchval, void* start_addr, void* end_addr, int stride, int op) {
    char value = (char)(PyInt_AsLong(searchval));
    SWITCH_ON_STRIDE(char)
}

static void* array_find_short(PyObject* searchval, void* start_addr, void* end_addr, int stride, int op) {
    short value = (short)(PyInt_AsLong(searchval));
    SWITCH_ON_STRIDE(short)
}

static void* array_find_long(PyObject* searchval, void* start_addr, void* end_addr, int stride, int op) {
    long value = PyInt_AsLong(searchval);
    SWITCH_ON_STRIDE(long)
}

static void* array_find_longlong(PyObject* searchval, void* start_addr, void* end_addr, int stride, int op) {
    long long value = PyLong_AsLongLong(searchval);
    SWITCH_ON_STRIDE(long long)
}

static void* array_find_ubyte(PyObject* searchval, void* start_addr, void* end_addr, int stride, int op) {
    unsigned char value = (unsigned char)(PyInt_AsUnsignedLongMask(searchval));
    SWITCH_ON_STRIDE(unsigned char)
}

static void* array_find_ushort(PyObject* searchval, void* start_addr, void* end_addr, int stride, int op) {
    unsigned short value = (unsigned short)(PyInt_AsUnsignedLongMask(searchval));
    SWITCH_ON_STRIDE(unsigned short)
}

static void* array_find_ulong(PyObject* searchval, void* start_addr, void* end_addr, int stride, int op) {
    unsigned long value = PyInt_AsUnsignedLongMask(searchval);
    SWITCH_ON_STRIDE(unsigned long)
}

static void* array_find_ulonglong(PyObject* searchval, void* start_addr, void* end_addr, int stride, int op) {
    unsigned long long value = PyLong_AsUnsignedLongLong(searchval);
    SWITCH_ON_STRIDE(unsigned long long)
}

static void* array_find_float(PyObject* searchval, void* start_addr, void* end_addr, int stride, int op) {
    float value = (float)(PyFloat_AsDouble(searchval));
    SWITCH_ON_STRIDE(float)
}

static void* array_find_double(PyObject* searchval, void* start_addr, void* end_addr, int stride, int op) {
    double value = PyFloat_AsDouble(searchval);
    SWITCH_ON_STRIDE(double)
}

static void* array_find_eq_string(char* value, void* start_addr, void* end_addr, int stride, int elsize)
{
    char dest[elsize+1];
    dest[elsize] = '\0';
    void* i;
    for (i = start_addr; i != end_addr; i += stride)
    {
        memcpy((void*)dest, i, elsize);
        if (strcmp(dest, value) == 0)
        {
            return i;
        }
    }
    return 0;
}

static void* array_find_neq_string(char* value, void* start_addr, void* end_addr, int stride, int elsize)
{
    char dest[elsize+1];
    dest[elsize] = '\0';
    void* i;
    for (i = start_addr; i != end_addr; i += stride)
    {
        memcpy((void*)dest, i, elsize);
        if (strcmp(dest, value) != 0)
        {
            return i;
        }
    }
    return 0;
}

static void* array_find_string(PyObject* searchval, void* start_addr, void* end_addr, int stride, int op, int elsize) {
    char* value = PyString_AsString(searchval);
    if (op == NEQ)
    {
        return array_find_neq_string(value, start_addr, end_addr, stride, elsize);
    }
    else
    {
        return array_find_eq_string(value, start_addr, end_addr, stride, elsize);
    }
}

static PyObject* multi_find(PyObject* self, PyObject* args, short direction)
{
    PyArrayObject* arr;
    PyObject* value;
    short eqtest;
    int start, end;
    if (!PyArg_ParseTuple(args, "O!Ohii", &PyArray_Type, &arr, &value, &eqtest, &start, &end))
    {
        return NULL;
    }
    char type = arr->descr->type;
    int elsize = arr->descr->elsize;
    int stride = arr->strides[0];

    void* base_addr = PyArray_DATA(arr);
    void* start_addr;
    void* end_addr;
    int dirstride = stride;
    void* temp;

    if (stride > 0)
    {
        start_addr = base_addr + stride * start;
        end_addr = base_addr + stride * end;
    }
    else
    {
        start_addr = base_addr + stride * start;
        end_addr = base_addr + stride * end;
    }

    if (!direction)
    {
        temp = start_addr;
        //Shift by one place so we read the same set of data backwards
        start_addr = end_addr - stride;
        end_addr = temp - stride;
        dirstride = -1 * stride;
    }

    void* dest_addr;
    switch (type)
    {
        case 'b':
            dest_addr = array_find_byte(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'h':
            dest_addr = array_find_short(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'i':
            dest_addr = array_find_long(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'l':
            dest_addr = array_find_long(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'q':
            dest_addr = array_find_longlong(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'B':
            dest_addr = array_find_ubyte(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'H':
            dest_addr = array_find_ushort(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'I':
            dest_addr = array_find_ulong(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'L':
            dest_addr = array_find_ulong(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'Q':
            dest_addr = array_find_ulonglong(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'f':
            dest_addr = array_find_float(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'd':
            dest_addr = array_find_double(value, start_addr, end_addr, dirstride, eqtest);
            break;
        case 'S':
            dest_addr = array_find_string(value, start_addr, end_addr, dirstride, eqtest, elsize);
            break;
        default:
            dest_addr = 0;
    }

    // printf("%s\t%s\t%d\t%d\t%p\t%p\t%p\t%p\n", &(arr->descr->kind), &type, elsize, dirstride, 
    //                                            base_addr, start_addr, end_addr, dest_addr);

    int dest_index;
    if (dest_addr > 0)
    {
        dest_index = (dest_addr - base_addr) / stride;
    }
    else
    {
        dest_index = -1;
    }
    return Py_BuildValue("i", dest_index);
}

static char py_find_doc[] = "Finds the first index of the occurrence of a value in an array.";
static PyObject* py_find(PyObject* self, PyObject* args)
{
    return multi_find(self, args, 1);
}

static char py_rfind_doc[] = "Finds the last index of the occurrence of a value in an array.";
static PyObject* py_rfind(PyObject* self, PyObject* args)
{
    return multi_find(self, args, 0);
}

static PyMethodDef _npfindmethods[] = {
    {"np_find", py_find, METH_VARARGS, py_find_doc},
    {"np_rfind", py_rfind, METH_VARARGS, py_rfind_doc},
    {NULL, NULL, 0, NULL}
};

void init_npfind(void)
{
    Py_InitModule("_npfind", _npfindmethods);
    import_array();
}
