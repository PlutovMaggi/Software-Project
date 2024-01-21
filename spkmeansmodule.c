#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"


/*fit function - kmeans algorithm from hw2*/
/*get arguments, apply kmeans algorithm and return final results*/
static PyObject* kmeans(PyObject *self, PyObject *args){
    /*variables declared*/
    int n,d,k,max_iter,i,j;
    double epsilon;
    PyObject *Pcentroids;
    PyObject *Pdatapoints; /*vectors*/
    PyObject *Pcentriod;
    PyObject *Pinner_centriod;
    PyObject *Pdatapoint;
    PyObject *Pinner_datapoint;
    PyObject *Presults;
    PyObject *Presult;
    PyObject *Pinner_result;
    double **centroids;
    double **datapoints;

    /*get arguments from python*/
    if(!PyArg_ParseTuple(args,"iiiidOO",&n,&d,&k,&max_iter,&epsilon,&Pcentroids,&Pdatapoints)) /*in case of an error*/
    {
        return NULL;
    }

    /*initialize datapoints(nXd) and centroids(kXd)*/
    datapoints = allocate_mat(n,d); /*allocating n vectors,each with double array in it*/

    if(datapoints==NULL) /*in case of an error*/
    {
        return NULL;
    }

    centroids = allocate_mat(k,d); /*allocating k centroids,each with double array in it*/

    if(centroids==NULL) /*in case of an error*/
    {
        free_mat(datapoints,n);
        return NULL;
    }

    /*get data from Pdatapoints*/
    for (i = 0;i<n;i++)
    {
        Pdatapoint = PyList_GetItem(Pdatapoints,i);

        if(!PyList_Check(Pdatapoint)) /*in case of an error*/
        {   
            free_mat(datapoints,n);
            return NULL;
        }

        for (j = 0;j<d;j++)
        {
            Pinner_datapoint = PyList_GetItem(Pdatapoint ,j);
            if (!PyFloat_Check(Pinner_datapoint)) /*in case of an error*/
            {
                free_mat(datapoints,n);
                return NULL;
            }
            datapoints[i][j] = PyFloat_AsDouble(Pinner_datapoint);
        }
    }    

    /*get data from Pcentroids*/
    for (i = 0;i<k;i++)
    {
        Pcentriod = PyList_GetItem(Pcentroids,i);

        if(!PyList_Check(Pcentriod)) /*in case of an error*/
        {   
            free_mat(centroids,k);
            return NULL;
        }

        for (j = 0;j<d;j++)
        {
            Pinner_centriod = PyList_GetItem(Pcentriod ,j);
            if (!PyFloat_Check(Pinner_centriod)) /*in case of an error*/
            {
                free_mat(centroids,k);
                return NULL;
            }
            centroids[i][j] = PyFloat_AsDouble(Pinner_centriod);
        }
    }

    /*running kmeans algorithm*/
    if (kmeans_calc(datapoints,centroids,n,d,k,max_iter,epsilon) != SUCCESSFULL_RUN)
    {
        free_mat(datapoints,n);
        free_mat(centroids,k);
        return NULL;
    }

    free_mat(datapoints,n); /*deallocate datapoints*/

    /*return results*/
    Presults = PyList_New(k);
    if (Presults == NULL) /*in case of an error*/
    {   
        free_mat(centroids,k);
        return NULL;
    }

    for (i = 0;i<k;i++)
    {
        Presult = PyList_New(d);
        if (Presult == NULL) /*in case of an error*/
        {
            free_mat(centroids,k);
            return NULL;
        }

        for( j = 0; j<d;j++)
        {
            Pinner_result = PyFloat_FromDouble(centroids[i][j]);
            if (Pinner_result == NULL) /*in case of an error*/
            {
                free_mat(centroids,k);
                return NULL;
            }
            PyList_SET_ITEM(Presult,j,Pinner_result);   
        }
        PyList_SET_ITEM(Presults,i,Presult);  
    }

    free_mat(centroids,k); /*deallocate centeroids*/

    /*return final centroids*/
    return Py_BuildValue("O", Presults);
}

/*python cases reader*/
static PyObject* pythonReader(PyObject *self, PyObject *args)
{
    char *Pygoal;
    double **datapoints;
    double **finalresult;
    int n,k,d,i,j;    
    PyObject *Pdatapoints;
    PyObject *Pdatapoint;
    PyObject *Pinner_datapoint;
    PyObject *Presults;
    PyObject *Presult;
    PyObject *Pinner_result;
    
    if(!PyArg_ParseTuple(args,"iiisO",&n,&k,&d,&Pygoal,&Pdatapoints)) /*in case of an error*/
    {
        return NULL;
    }

    /*initialize datapoints(nXd) and centroids(kXd)*/
    datapoints = allocate_mat(n,d); /*allocating n vectors,each with double array in it*/
    
    if(datapoints==NULL) /*in case of an error*/
    {
        return NULL;
    }

    for (i = 0;i<n;i++)
    {
        Pdatapoint = PyList_GetItem(Pdatapoints,i);

        if(!PyList_Check(Pdatapoint)) /*in case of an error*/
        {   
            free_mat(datapoints,n);
            return NULL;
        }

        for (j = 0;j<d;j++)
        {
            Pinner_datapoint = PyList_GetItem(Pdatapoint ,j);
            if (!PyFloat_Check(Pinner_datapoint)) /*in case of an error*/
            {
                free_mat(datapoints,n);
                return NULL;
            }
            datapoints[i][j] = PyFloat_AsDouble(Pinner_datapoint);
        }
    }
    
    /*c help function*/
    finalresult = pythonCases(datapoints,n,k,d,Pygoal);

    /*return matrix - according to each case*/
    if(strcmp(Pygoal,"spk") == 0)
    {
        /*update k to the calculated value if needed*/
        if (k == 0) 
        {
            k = update_k();
        }
        Presults = PyList_New(n);
        if (Presults == NULL) /*in case of an error*/
        {   
            free_mat(finalresult,n);
            return NULL;
        }

        for (i = 0;i<n;i++)
        {
            Presult = PyList_New(k);
            if (Presult == NULL) /*in case of an error*/
            {
                free_mat(finalresult,n);
                return NULL;
            }

            for(j = 0; j<k;j++)
            {
                Pinner_result = PyFloat_FromDouble(finalresult[i][j]);
                if (Pinner_result == NULL) /*in case of an error*/
                {
                    free_mat(finalresult,n);
                    return NULL;
                }
                PyList_SET_ITEM(Presult,j,Pinner_result);   
            }
            PyList_SET_ITEM(Presults,i,Presult);
        }
        free_mat(finalresult,n); 
    }
    
    else
    {
        if ((strcmp(Pygoal,"wam") == 0) || (strcmp(Pygoal,"ddg") == 0) || (strcmp(Pygoal,"lnorm") == 0))
        {
            Presults = PyList_New(n);
            if (Presults == NULL) /*in case of an error*/
            {   
                free_mat(finalresult,n);
                return NULL;
            }
        
            for (i=0 ;i<n;i++)
            {
                Presult = PyList_New(n);
                if (Presult == NULL) /*in case of an error*/
                {
                    free_mat(finalresult,n);
                    return NULL;
                }

                for(j = 0; j<n;j++)
                {
                    Pinner_result = PyFloat_FromDouble(finalresult[i][j]);
                    if (Pinner_result == NULL) /*in case of an error*/
                    {
                        free_mat(finalresult,n);
                        return NULL;
                    }
                    PyList_SET_ITEM(Presult,j,Pinner_result);   
                }
                PyList_SET_ITEM(Presults,i,Presult);
            }
            free_mat(finalresult,n);  
        }

        else /*goal = jacobi */
        {
            Presults = PyList_New(n+1);
            if (Presults == NULL) /*in case of an error*/
            {   
                free_mat(finalresult,n);
                return NULL;
            }

            for (i=0 ;i<n+1;i++)
            {
                Presult = PyList_New(n);
                if (Presult == NULL) /*in case of an error*/
                {
                    free_mat(finalresult,n);
                    return NULL;
                }

                for(j = 0; j<n;j++)
                {
                    Pinner_result = PyFloat_FromDouble(finalresult[i][j]);
                    if (Pinner_result == NULL) /*in case of an error*/
                    {
                        free_mat(finalresult,n);
                        return NULL;
                    }
                    PyList_SET_ITEM(Presult,j,Pinner_result);   
                }
                PyList_SET_ITEM(Presults,i,Presult);
            }
            /*free final result*/
            free_mat(finalresult,n+1);
        }
    } 

    /*return final result*/
    return Py_BuildValue("O", Presults);
}

/*for python*/
static PyMethodDef spkmeansMethods[] = {
    {"cases", (PyCFunction) pythonReader, METH_VARARGS, PyDoc_STR("all cases in C")},
    {"fit", (PyCFunction) kmeans, METH_VARARGS, PyDoc_STR("return centroids")},
    {NULL,NULL,0,NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, 
    "spkmeansmodule",
    NULL,
    -1,
    spkmeansMethods
};

PyMODINIT_FUNC PyInit_spkmeansmodule (void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if(!m) {
        return NULL;
    }
    return m;
}