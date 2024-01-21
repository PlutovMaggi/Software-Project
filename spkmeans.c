#include "spkmeans.h"

/*help functions declared in header file*/

/*global variables declered*/
int k,n,d = 1; /* k - number of clusters, n - number of data points, d - dimension of data point*/
double **datapoints; /*vectors*/
double **weights; /*WAM*/
double **diagonal; /*DDG*/
double **diagonalMinusHalf;/*helper for lnorm*/
double **laplacian;/*LNORM*/
double *eigenvalues;/*jacobi eigenvalues*/
double **eigenvectors;/*jacobi eigenvectors*/
double **finalresult; /*pointer to return value*/

/*main (C) function:*/
/*recieving arguments from cmd - intended for C program direct use*/
int main(int argc, char *argv[])
{   
    /*declerations*/
    char *input_file_name;
    char *cur_goal;

    /*get 2 argument*/
    if (argc != 3)
    {		
		return exitCodeInvalidInput();
    }
    cur_goal = argv[1];
    input_file_name = argv[2];

    /*asssert input argument*/
    if (assert_input(cur_goal,input_file_name)!=SUCCESSFULL_RUN)
    {
        return exitCodeInvalidInput();
    }

    /*initialization*/
    if (initialize_data(input_file_name)!=SUCCESSFULL_RUN) /*in case of an error*/
    {
        return EXIT_CODE_ERROR; 
    }

    /*check goal cases - will treat each goal differently according to its requirements*/
    if(kmeans_switcher(cur_goal)!=SUCCESSFULL_RUN)/*in case of an error*/
    {
        return exitCodeError();
    }

    /*print and free allocations*/
    print_and_free(cur_goal);

    return SUCCESSFULL_RUN;  
} 


/*help-functions:*/

/*python functions*/
/*recieving arguments from python - intended for python program call*/
double** pythonCases(double **Pdatapoints,int py_n,int py_k, int py_d,char *Pygoal)
{
    double **Ueigen;
    double **Teigen;
    int* index_array;
    
    datapoints = Pdatapoints;
    n = py_n;
    k = py_k;
    d = py_d;

    if(strcmp(Pygoal,"spk") == 0)
    {
        /*steps 1+2 - form W & compute Lnorm*/
        if(compute_all() !=SUCCESSFULL_RUN) /*in case of an error*/
        {
            return NULL;
        } 

        /*create indices array sorted by decreasing eigenvalues order*/
        index_array = index_sort();

        if (index_array == NULL) /*in case of an error*/
        {
            free_mat(eigenvectors,n);
            free(eigenvalues);
            return NULL;
        }

        /*step 3 - determine global variable k*/
        if(k==0)
        {
            /*update global variable k*/
            determine_k(index_array); 
        }
        
        /*step 4 - form U matrix*/
        Ueigen = form_U(index_array);
        
        if(Ueigen==NULL) /*in case of an error*/
        {
            free_mat(eigenvectors,n);
            free(eigenvalues);
            return NULL;
        }

        free_mat(eigenvectors,n);
        free(eigenvalues);
        free(index_array);

        /*step 5 - form T matrix*/
        Teigen = form_T(Ueigen);
        
        if(Teigen==NULL) /*in case of an error*/
        {
            free_mat(Ueigen,n);
            return NULL;
        }

        free_mat(Ueigen,n);
        finalresult = Teigen;
    }

    /*if goal != spk*/
    else
    {
        kmeans_switcher(Pygoal);
        
        /*jacobi - special case*/
        if(strcmp(Pygoal,"jacobi") == 0) 
        {
            /*copy eigenvlues & eigenvectors to final result matrix*/
            if (jacobi_finalresult() !=SUCCESSFULL_RUN) /*in case of an error*/
            {
                free_mat(eigenvectors,n);
                free(eigenvalues);
                return NULL; 
            }
            free(eigenvalues);
            free_mat(eigenvectors,n);
        }    
    }

    /*return final result*/
    return finalresult;
}

/*return updated value of k*/
int update_k(void)
{
    return k;
}


/*cases function*/
/*switch between goal cases*/
int kmeans_switcher(char *cur_goal)
{ 
    /*calculate result for each case -for some has to do previous cases*/
    /*wam case*/
    if (strcmp(cur_goal,"wam") == 0)
    {
        if (weightMat_calc()!=SUCCESSFULL_RUN) /*in case of an error*/
        {
            free_mat(datapoints,n);
            return EXIT_CODE_ERROR;
        }
        free_mat(datapoints,n);
        finalresult = weights;
    }
    /*ddg case*/
    if (strcmp(cur_goal,"ddg") == 0)
    {
        if (weightMat_calc()!=SUCCESSFULL_RUN) /*in case of an error*/
        {
            free_mat(datapoints,n);
            return EXIT_CODE_ERROR;
        }
        if (diagonalMat_calc()!=SUCCESSFULL_RUN) /*in case of an error*/
        {
            free_mat(datapoints,n);
            free_mat(weights,n);
            return EXIT_CODE_ERROR;
        }
        free_mat(datapoints,n);
        free_mat(weights,n);
        free_mat(diagonalMinusHalf,n);
        finalresult = diagonal;
    }
    /*lnorm case*/
    if (strcmp(cur_goal,"lnorm") == 0)
    {
        if (weightMat_calc()!=SUCCESSFULL_RUN) /*in case of an error*/
        {
            free_mat(datapoints,n);
            return EXIT_CODE_ERROR;
        }
        if (diagonalMat_calc()!=SUCCESSFULL_RUN) /*in case of an error*/
        {
            free_mat(datapoints,n);
            free_mat(weights,n);
            return EXIT_CODE_ERROR;
        }
        if (lnormMat_calc()!=SUCCESSFULL_RUN) /*in case of an error*/
        {
            free_mat(datapoints,n);
            free_mat(weights,n);
            free_mat(diagonal,n);
            free_mat(diagonalMinusHalf,n);
            return EXIT_CODE_ERROR;
        }
        free_mat(datapoints,n);
        free_mat(diagonal,n);
        free_mat(diagonalMinusHalf,n);
        free_mat(weights,n);
        finalresult = laplacian;
    }
    /*jacobi case*/
    if (strcmp(cur_goal, "jacobi") == 0)
    {
        /*calculate jacobi*/
        if (jacobi_calc(datapoints)!=SUCCESSFULL_RUN) /*in case of an error*/
        {
            free_mat(datapoints,n);
            return EXIT_CODE_ERROR;
        }
        free_mat(datapoints,n);
    }
        
    return SUCCESSFULL_RUN;  
}


/*errors functions*/
int exitCodeError(void)
{
    printf("An Error Has Occured\n");
    return EXIT_CODE_ERROR;
}

int exitCodeInvalidInput(void)
{
    printf("Invalid Input!\n");
    return EXIT_CODE_INVALID_INPUT;
}


/*assert & initialize data functions (for reciving arguments from user)*/
/*reads N data points - final data point will be placed in the global variable datapoints*/
int initialize_data(char *input_file_name)
{
    /*determine datapoints dimension - nXd*/
    if(determine_dim(input_file_name)!=SUCCESSFULL_RUN) /*in case of an error*/
    {
        return EXIT_CODE_ERROR;
    }

    /*initialize datapoints to global variable*/
    if(read_data(input_file_name)!=SUCCESSFULL_RUN) /*in case of an error*/
    {   
        return EXIT_CODE_ERROR;
    }

    return SUCCESSFULL_RUN;
}

/*assert input argiments - goal + input file name*/
int assert_input(char *cur_goal,char *input_file_name)
{
    /*asssert goal argument*/
    if ((strcmp(cur_goal,"wam") != 0) && (strcmp(cur_goal,"ddg") != 0) && (strcmp(cur_goal,"lnorm") != 0 ) && (strcmp(cur_goal,"jacobi") != 0 ))
    {
        return EXIT_CODE_INVALID_INPUT;
    }

    /*assert input file argument*/
    if(strlen(input_file_name) < 5 || (strcmp(input_file_name + strlen(input_file_name) - 4, ".txt") !=0 && strcmp(input_file_name + strlen(input_file_name) - 4, ".csv") !=0))
    {    
        return EXIT_CODE_INVALID_INPUT;
    }

    return SUCCESSFULL_RUN;
}

/*determine global variables n and d - datapoint dimensions. internal prints*/
int determine_dim(char *input_file_name)
{
    int curr_char,prev_char,flag;
    FILE *input;

    /*read file*/
    input = fopen(input_file_name, "r");

    if(input==NULL) /*in case of an error*/
    {
        return exitCodeInvalidInput();
    }

    curr_char = '\0';
    prev_char = '\n';
    flag = 1; /*indicates the first row - to count d*/

    /*find size of vectors (d) and number of vectors (n)*/
    while ((curr_char = fgetc(input)) && (curr_char != EOF))
    {
        if (curr_char == ',' && flag)
        {
            d++;
        }
        if (curr_char == '\n'  &&  prev_char != '\n')
        {
            n++;
            flag = 0;
        }
        prev_char = curr_char;
    }
    /*find n*/
    if (prev_char != '\n')
    {
        n++;
    }
    /*close file*/
    if(fclose(input)!=0) /*in case of an error*/
    {
        return exitCodeError();
    }

    return SUCCESSFULL_RUN;
}

/*read data from file & initialize global variable datapoints. internal prints*/
int read_data(char *input_file_name)
{
    int i,j;
    FILE *input;

    input = fopen(input_file_name, "r");

    if(input==NULL) /*in case of an error*/
    {   
        return exitCodeInvalidInput();
    }

    /*allocating n vectors,each with double array in it*/
    datapoints = allocate_mat(n,d);
    /*assert allocation*/
    if(datapoints==NULL) 
    {
        return exitCodeError();
    }

    for (i=0;i<n;i++)
    { 
        for (j=0; j<d; j++)
        {
            if (fscanf(input, "%lf%*c", &datapoints[i][j]) <= 0)
            {
                free_mat(datapoints,n);
                return exitCodeError();
            }
        }
    }

    if(fclose(input)!=0) /*in case of an error*/
    {
        free_mat(datapoints,n);
        return exitCodeError();
    }

    return SUCCESSFULL_RUN;
}


/*SPK case functions*/
/*compare function - comperator for decreasing order of eigenvalues*/
/*assumes successful initialization of necessary global variables(eigenvalues)*/
int compare_val(const void *p1,const void *p2)
{
    double result;

    result = eigenvalues[*(int *)p1] - eigenvalues[*(int *)p2];
    if (result > 0)
    {
        return -1;
    }
    else
    {
        if (result < 0)
        {
            return 1;
        }  
        else /* result == 0 */
        {
            return 0;
        }
    } 
}

/*calculate all needed matrices (steps 1+2 - from W to Laplacian (+JACOBI))*/
/*assumes successful initialization of necessary global variables(datapoints)*/
int compute_all(void)
{
    /*1 - calculate weighted adjacency matrix*/
    if (weightMat_calc()!=SUCCESSFULL_RUN) /*in case of an error*/
    {
        free_mat(datapoints,n);
        return EXIT_CODE_ERROR;
    }
    free_mat(datapoints,n);
    /*calculate diaginal degree matrix*/
    if (diagonalMat_calc()!=SUCCESSFULL_RUN) /*in case of an error*/
    {
        free_mat(weights,n);
        return EXIT_CODE_ERROR;
    }
    /*2 - calculate normalized graph laplacian*/
    if (lnormMat_calc()!=SUCCESSFULL_RUN) /*in case of an error*/
    {
        free_mat(diagonal,n);
        free_mat(diagonalMinusHalf,n);
        return EXIT_CODE_ERROR;
    }
    free_mat(weights,n);
    free_mat(diagonal,n);
    free_mat(diagonalMinusHalf,n);
    /*3 - calculate eigenvalues & eigenvectors*/
    if (jacobi_calc(laplacian)!=SUCCESSFULL_RUN) /*in case of an error*/
    {
        free_mat(laplacian,n);
        return EXIT_CODE_ERROR;
    }
    free_mat(laplacian,n);

    return SUCCESSFULL_RUN;
}

/*if k==0 - determine global variable k*/
/*assumes successful initialization of necessary global variables(eigenvalues&eigenvectors)*/
void determine_k(int *index_array)
{
    int max_index = 0,i;
    double max_gap=-1;

    /*find k*/
    for (i=0;i<(n/2);i++)
    {
        if(fabs(eigenvalues[index_array[i]]-eigenvalues[index_array[i+1]])>max_gap)
        {
            max_gap = fabs(eigenvalues[index_array[i]]-eigenvalues[index_array[i+1]]);
            max_index = i;
        }
    }
    k = max_index+1;
}

/*form & return calculated nXk U matrix*/
/*assumes successful initialization of necessary global variables(eigenvalues&eigenvectors)*/
double** form_U(int *index_array)
{
    int i,j;
    double **Ueigen;

    /*alocate nXk U matrix*/
    Ueigen = allocate_mat(n,k);
    
    if(Ueigen==NULL) /*in case of an error*/
    {
        return NULL;
    }

    /*calculate U*/
    for (i=0;i<k;i++)
    {
        for(j=0;j<n;j++)
        {
            Ueigen[j][i]=eigenvectors[j][index_array[i]];
        }
    }

    return Ueigen;
}

/*form & return calculated nXk T matrix*/
double** form_T(double **Ueigen)
{
    int i,j;
    double row_sum;
    double **Teigen;

    /*allocate T*/
    Teigen = allocate_mat(n,k);
    
    if(Teigen==NULL) /*in case of an error*/
    {
        return NULL;
    }

    /*calculate T*/
    for (i=0;i<n;i++)
    {
        row_sum = 0;
        for (j=0;j<k;j++)
        {
            row_sum += pow(Ueigen[i][j],2);
        }
        row_sum= pow(row_sum,0.5);
        for (j=0;j<k;j++)
        {
            if(row_sum == 0)
            {
                Teigen[i][j] = 0;
            }
            else
            {   
                Teigen[i][j] = Ueigen[i][j]/row_sum;
            }
        }    
    }

    return Teigen;
}

/*returns the order of indices sorted by decreasing eigenvalues*/
int* index_sort(void)
{
    int i;
    int *index_array;

    index_array = (int *) malloc((int)(sizeof(int))*n); /*allocating nX1 vector*/

    if (index_array == NULL)
    {
        return NULL;
    }
    /*index_array = [0,1,...,n-1]*/
    for (i = 0;i<n;i++)
    {
        index_array[i] = i;
    }

    qsort(index_array,n,sizeof(int),compare_val);

    return index_array;
}


/*WAM case functions*/
/*calculate W=weighted adjacency matrix*/
/*assumes successful initialization of necessary global variables(datapoints)*/
int weightMat_calc(void)
{
    int i,j;
    double result;

    /*allocating W*/
    weights = allocate_mat(n,n); 

    /*assert allocation*/
    if(weights==NULL)
    {
        return EXIT_CODE_ERROR;
    }

    /*calculate W*/
    for (i=0;i<n;i++)
    {
        for (j=i+1;j<n;j++)
        {
            result = exp(euclidean_norm(datapoints[i],datapoints[j])*(-0.5));
            if (result > 0)
            {
                weights[i][j] = result;
                weights[j][i] = weights[i][j] ; /*symmetry*/
            }
        }
    }
    return SUCCESSFULL_RUN;
}

/*euclidian norm calculation between two 1Xd vectors*/
double euclidean_norm(double *vec1,double *vec2 )
{
    int i;
    double delta;

    delta = 0;
    for(i=0;i<d;i++)
    {
        delta+=pow((vec1[i]-vec2[i]),2);
    }
    delta = pow(delta,0.5);

    return delta;
}


/*DDG case functions*/
/*calculate D=diagonal degree matrix (& D^(-0.5)) */
/*assumes successful initialization of necessary global variables(weights)*/
int diagonalMat_calc(void)
{
    int i,j;
    double result;

    /*allocating D*/
    diagonal = allocate_mat(n,n); 
    
    /*assert allocation*/
    if(diagonal==NULL) 
    {
        return EXIT_CODE_ERROR;
    }

    /*allocating D ^(-0.5) */
    diagonalMinusHalf = allocate_mat(n,n); 

    /*assert allocation*/
    if(diagonalMinusHalf==NULL) 
    {
        free_mat(diagonal,n);
        return EXIT_CODE_ERROR;
    }

    /*calculate D*/
    for (i=0;i<n;i++)
    {
        result = 0;
        for (j=0;j<n;j++)
        {
            result += weights[i][j];
        }
        diagonal[i][i] = result;
        diagonalMinusHalf[i][i] = 1/(pow(result,0.5)); 
    }

    return SUCCESSFULL_RUN;
}


/*LNORM case functions*/
/*calculate normalized graph laplacian = I-D(^-0.5)*W*D(^-0.5)*/
/*based on D(^-0.5) being diagonal matrix*/
/*assumes successful initialization of necessary global variables(weights & diagonalMinusHalf)*/
int lnormMat_calc(void) 
{
    int i,j;

    laplacian = allocate_mat(n,n); /*allocating matrix*/

    /*assert allocation*/
    if(laplacian==NULL) 
    {
        return EXIT_CODE_ERROR;
    }

    for(i=0;i<n;i++)    
    {
        for(j=0;j<n;j++)    
        {       
            laplacian[i][j] = diagonalMinusHalf[i][i] * weights[i][j] * diagonalMinusHalf[j][j];

            if(i==j)
            {
                laplacian[i][j] = 1-laplacian[i][j];
            }
            else
            {
                laplacian[i][j]= 0-laplacian[i][j];
            }   
        }    
    }

    return SUCCESSFULL_RUN;
}


/*JACOBI case functions*/
/*calculate eigenvalues and eigenvectors by jacobi algorithm*/
/*assumes given matrix is not diagonal(in particular dimensions != 1X1)*/
int jacobi_calc(double **A)
{
    int iter = 1, max_i = -1,max_j = -1,i,j;
    double epsilon = 0.00001,teta,s,t,c,result;
    double max_val = 0; /*max_val = max A[i][j] in pivot step*/
    double **A_tag; /* A' */
    double **mult_P; /*final eigenvectors - V=P1*P2.., P - rotation matrix  */
    
    /*allocate initialize mult_P & A_tag*/
    mult_P = allocate_mat(n,n); 
    if(mult_P==NULL) /*in case of an error*/
    {
        return EXIT_CODE_ERROR;
    }

    A_tag = allocate_mat(n,n);
    if(A_tag==NULL) /*in case of an error*/
    {
        free_mat(mult_P,n);
        return EXIT_CODE_ERROR;
    }

    /*initialize mult_P = I(identity matrix)*/
    init_I(mult_P);

    do
    {
        max_val = 0;
        /* pivot step - find Aij (off diagonal element with largest absolute value*/
        for(i=0;i<n;i++)
        {
            for (j=i+1;j<n;j++)
            {
                if (fabs(A[i][j])>max_val)
                {
                    max_i = i;
                    max_j = j;
                    max_val = fabs(A[i][j]);
                }
            }
        }

        /*in case of diagonal matrix - special case(discussed in the forum)*/
        if (max_val == 0)
        {
            break;
        }

        /*obtain c,t*/
        teta = (A[max_j][max_j]-A[max_i][max_i])/(2*A[max_i][max_j]);
        if (teta>=0)
        {
            t = 1/(teta+pow((pow(teta,2)+1),0.5));
        }
        else
        {
            t = -1/((-1*teta)+pow((pow(teta,2)+1),0.5));
        }  
        c = 1/pow(pow(t,2)+1,0.5);
        s = t*c;

        /*build A_tag*/
        update_A_tag(A_tag,A,max_i,max_j,c,s);

        /*update mult_P=V=P1*P2.. (multiply from the right with new rotation matrix P*/
        update_mult_P(mult_P,max_i,max_j,s,c);

        iter +=1;
        result = off(A)-off(A_tag);

        /*copying values*/
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                A[i][j] = A_tag[i][j];
            }
        }
    } while ((result>epsilon) && (iter<=100));

    /*free A_tag (keep in global variable eigenvectors = mult_p)*/
    /*A = data points - memory released in outer function*/
    free_mat(A_tag,n);

    /*update global eigenvectors & eigenvalues */
    eigenvectors = mult_P;
    eigenvalues = (double *) malloc((int)(sizeof(double))*n); /*allocating nX1 vector*/

    if(eigenvalues==NULL) /*in case of an error*/
    {
        free_mat(mult_P,n);
        return EXIT_CODE_ERROR;
    }

    for(i=0; i<n; i++)
    {
        eigenvalues[i] = A[i][i];
    }
    
    return SUCCESSFULL_RUN;
}

/*calculate off value of matrix*/
double off(double **mat)
{
    int i,j;
    double result = 0;
    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {
            if (j!=i)
            {
                result += pow(mat[i][j],2);
            }
        }
    }
    return result;
}

/*initialize nXn matrix to I-identity matrix*/
void init_I(double **mat)
{
    int i,j;

    for(i=0; i<n; i++)
    {
        for(j=0;j<n;j++)
        {
            if(i==j)
            {
                mat[i][j] = 1;
            }
            else
            {
                mat[i][j] = 0;
            }
        }
    }
}

/*build nXn A' matrix - transformation of matrix A*/
void update_A_tag(double **A_tag,double **A,int max_i,int max_j,double c,double s)
{
    int i,j;

    /*build A_tag*/
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if (i==max_i && j==max_i)/*A'[i][i]*/
            {
                A_tag[i][j]= pow(c,2)*A[max_i][max_i] + pow(s,2)*A[max_j][max_j] - 2*s*c*A[max_i][max_j];
            }
            if (i==max_j && j==max_j)/*A'[j][j]*/
            {
                A_tag[i][j]= pow(s,2)*A[max_i][max_i] + pow(c,2)*A[max_j][max_j] + 2*s*c*A[max_i][max_j];
            }
            if (i==max_i && j==max_j)/*A'[i][j] = A'[j][i] = 0*/
            {
                A_tag[i][j] = 0;
                A_tag[j][i] = 0;
            }
            if (j==max_i && i!=max_i && i!=max_j)/*A'[r][i] = A'[i][r] */
            {
                A_tag[i][j] = c*A[i][max_i] - s*A[i][max_j];
                A_tag[j][i] = A_tag[i][j];
            }
            if (j==max_j && i!=max_i && i!=max_j) /*A'[r][j] = A'[j][r]*/
            {
                A_tag[i][j] = c*A[i][max_j] + s*A[i][max_i];
                A_tag[j][i] = A_tag[i][j];
            }
            /*other rows & columns*/
            if (i!=max_i && j!=max_j && i!= max_j && j!=max_i)
            {
                A_tag[i][j] = A[i][j];
            }
        }
    }
}

/*allocate global variable finalresult for jacobi case (special case)*/
/*assumes successful initialization of necessary global variables(eigenvalues&eigenvectors)*/
int jacobi_finalresult(void)
{
    int i,j;

    /*allocate finalresult*/
    finalresult = allocate_mat(n+1,n);

    if (finalresult == NULL) /*in case of an error*/
    {
        return EXIT_CODE_ERROR; 
    }

    /*copy eigenvlues & eigenvectors to final result matrix*/
    for (i=0;i<n+1;i++)
    {
        for (j=0;j<n;j++)
        {
            if(i==0) /*first row - eigenvalues*/
            {
                finalresult[i][j] = eigenvalues[j];
            }
            else /*other - eigenvectors*/
            {
                finalresult[i][j] = eigenvectors[i-1][j];
            }   
        }
    }

    return SUCCESSFULL_RUN;
}

/*update V - multiply rotaion matrices to calculate V = P1*P2...*/
/*each time multiply current V with Pi from the right, max_i < max_j*/
void update_mult_P(double **mult_P,int max_i,int max_j,double s, double c)
{
    int i;
    double a,b; /*keeps values in max_i & max_j columns respectively*/

    /*update based on P definition - update only max_i & max_j columns*/
    for(i=0;i<n;i++)    
    {
        a = mult_P[i][max_i];
        b = mult_P[i][max_j];
        mult_P[i][max_i] = c*a - s*b;
        mult_P[i][max_j] = s*a + c*b; 
    } 
} 


/*printing & (de)allocating memory functions*/
/*print matrix nXn*/
void print_mat(double **mat)
{
    int i,j;
    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {
            printf("%.4f",mat[i][j]);
            if (j != (n-1))
            {
                printf(",");
            }
            else
            {
                printf("\n");
            }
        }
    }  
}

/*deallocate memory for matrix with m rows*/
void free_mat(double **mat,int m)
{
    int i;
    for (i=0;i<m;i++)
    {
        free(mat[i]);
    }
    free(mat);
}

/*print eigenvalues and eigenvectors*/
void print_eigen(void)
{
    int i;
    /*print eigenvalues*/
    for (i=0;i<n;i++)
    {
        printf("%.4f",eigenvalues[i]);
        if (i != (n-1))
        {
            printf(",");
        }
        else
        {
            printf("\n");
        }
    }
    
    /*print eigenvectors*/
    print_mat(eigenvectors);
}

/*print final result and free allocated mamory*/
void print_and_free(char *cur_goal)
{
    /*for jacobi case - print eigenvalues & eigenvectors separately*/ 
    if (strcmp(cur_goal,"jacobi") == 0) 
    {
        print_eigen();
        free(eigenvalues);
        free_mat(eigenvectors,n);
    }
    /*for other cases*/
    else
    {
        print_mat(finalresult);
        free_mat(finalresult,n);
    }
}

/*allocate nXm double matrix*/
double** allocate_mat(int n, int m)
{
    int i;
    double** mat;

    mat = (double **) malloc((int)(sizeof(double *))*n); 
    /*assert allocation*/
    if(mat==NULL) /*in case of an error*/
    {
        return NULL;
    }

    for (i=0;i<n;i++)
    { 
        mat[i] = (double*)calloc(m,sizeof(double)); /*initialized to 0's*/
        /*assert allocation*/
        if(mat[i]==NULL) /*in case of an error*/
        {   
            free_mat(mat,i);
            return NULL;
        }
    }
    return mat;
}


/*kmeans algorithm functions*/
/*all functions used in spkmeansmodule*/
/*calculates centroids - code from HW2*/
int kmeans_calc(double **datapoints,double **centroids,int n,int d,int k,int max_iter,double epsilon)
{
    /*variables declared*/
    int numofiter,flag;
    double maxeuclid;
    double **nowcentroids; /*remains current sum*/
    int *clusters;/*counts how many in each bucket*/

    /*running kmeans algorithm*/
    flag = 1;
    numofiter = 0;
    while (flag && numofiter<max_iter)
    {
        numofiter +=1;

        /*allocate kXd nowcentroids help matrix*/
        nowcentroids = allocate_mat(k,d); /*remains the sum till now, we need k of this*/
        if(nowcentroids==NULL) /*in case of an error*/
        {   
            return EXIT_CODE_ERROR;
        }

        /*allocate 1Xk clusters help vector*/
        clusters = calloc(k,sizeof(int)); 
        if(clusters ==NULL) /*in case of an error*/
        {
            free_mat(nowcentroids,k);
            return EXIT_CODE_ERROR;
        }

        /*assign vectors*/
        assign_vec(datapoints,centroids,nowcentroids,clusters,n,k,d);

        /*update centroids & check convergence*/
        maxeuclid = update_centroids(centroids,nowcentroids,clusters,k,d);
        if (maxeuclid < epsilon)
        {
            flag = 0;
        }

        /*deallocate clusters & nowcentroids*/
        free(clusters);
        free_mat(nowcentroids,k);
    }

    return SUCCESSFULL_RUN;
}

/*assign vectors to the closest cluster*/
void assign_vec(double **datapoints,double **centroids,double **nowcentroids,int *clusters,int n,int k,int d)
{
    int i,j,m,minCentroid;
    double mindelta,delta;
    /*going through all datapoints*/
    for(i=0;i<n;i++)
    {
        minCentroid = -1;
        mindelta =  __DBL_MAX__;
        for(j=0;j<k;j++)
        {
            delta = 0;
            for(m=0;m<d;m++)
            {
                delta+=pow((datapoints[i][m]-centroids[j][m]),2);
            }
            delta = pow(delta,0.5);
            if (delta < mindelta) /*new min*/
            {
                mindelta = delta;
                minCentroid = j;
            }
        }
        
        clusters[minCentroid] += 1;

        for(j=0;j<d;j++)
        {
            nowcentroids[minCentroid][j] += datapoints[i][j];
        }
    }
}

/*update the centroids and return the maximum euclidean norm*/
int update_centroids(double **centroids,double **nowcentroids,int *clusters,int k,int d)
{
    int i,j;
    double euclideanNorm,maxeuclid;

    maxeuclid = -1;
    for(i=0;i<k;i++)
    {
        euclideanNorm=0;
        if (clusters[i] != 0)
        {
            for(j=0;j<d;j++)
            {
                nowcentroids[i][j] = nowcentroids[i][j]/(clusters[i]);
                euclideanNorm += pow((nowcentroids[i][j]-centroids[i][j]),2);
                centroids[i][j] = nowcentroids[i][j];
            }
            if (euclideanNorm>maxeuclid)
            {
                maxeuclid = euclideanNorm;
            }
        }
    }
    maxeuclid = pow(maxeuclid,0.5);
    
    return maxeuclid;
}





