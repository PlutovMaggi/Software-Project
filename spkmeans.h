#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define SUCCESSFULL_RUN         (0)
#define EXIT_CODE_INVALID_INPUT (1)
#define EXIT_CODE_ERROR         (1)


#ifndef HEADER_FILE
#define HEADER_FILE

/* help functions declared */

/*python functions*/
double** pythonCases(double **Pdatapoints,int py_n,int py_k, int py_d,char *Pygoal); /*switch between goal cases recieved from python program*/
int update_k(void); /*return updated value of k*/
/*cases function*/
int kmeans_switcher(char *cur_goal); /*switch between goal cases to creat final results*/
/*errors functions*/
int exitCodeError(void); /*prints "An Error Has Accured" in case of an unexpected error*/
int exitCodeInvalidInput(void); /*prints "Invalid Input!" in case of an input error*/
/*assert & initialize data functions (for reciving arguments from user)*/
int initialize_data(char *input_file_name); /*read N data points*/
int assert_input(char *cur_goal,char *input_file_name); /*assert input argiments*/
int determine_dim(char *input_file_name); /*determine datapoints dimension - nXd*/
int read_data(char *input_file_name); /*read data & initialize datapoints*/
/*SPK case functions*/
int compare_val(const void *p1,const void *p2); /*comperator for decreasing order of eigenvalues*/
int compute_all(void); /*calculate all needed matrices (steps 1+2 - from W to Laplacian (+JACOBI))*/
void determine_k(int *index_array); /*if k==0 - determine global variable k*/
double** form_U(int *index_array); /*form & return calculated nXk U matrix*/
double** form_T(double **Ueigen); /*form & return calculated nXk T matrix*/
int* index_sort(void); /*returns the order of indices sorted by decreasing eigenvalues*/
/*WAM case functions*/
int weightMat_calc(void); /*calculate W - weighted adjacency matrix*/
double euclidean_norm(double *vec1,double *vec2 ); /*euclidian norm calculation*/
/*DDG case functions*/
int diagonalMat_calc(void); /*calculate D - diagonal degree matrix*/
/*LNORM case functions*/
int lnormMat_calc(void); /*calculate lnorm - normalized graph laplacian*/
/*JACOBI case functions*/
int jacobi_calc(double **A); /*calculate eigenvalues and eigenvectors*/
double off(double **mat); /*calculate off value of matrix*/
void init_I(double **mat);/*initialize nXn matrix to I-identity matrix*/
void update_A_tag(double **A_tag,double **A,int max_i,int max_j,double c,double s); /*build nXn A' matrix - transformation of matrix A for jacobi algorithm*/
int jacobi_finalresult(void);/*allocate global variable finalresult for jacobi case (special case)*/
void update_mult_P(double **mult_P,int max_i,int max_j,double s, double c); /*update V - multiply rotaion matrices to calculate V = P1*P2...*/
/*printing & (de)allocating memory functions*/
void print_mat(double **mat); /*print matrix nXn*/
void free_mat(double **mat,int m); /*deallocate memory for matrix with m rows*/
void print_eigen(void); /*print eigenvalues and eigenvectors*/
void print_and_free(char *cur_goal); /*print final result and free allocated mamory*/
double** allocate_mat(int n, int m); /*allocate double matrix with n rows & m collumns*/
/*kmeans algorithm functions*/
int kmeans_calc(double **datapoints,double **centroids,int n,int d,int k,int max_iter,double epsilon);/*calculates centroids - code from HW2*/
void assign_vec(double **datapoints,double **centroids,double **nowcentroids,int *clusters,int n,int k,int d);/*assign vectors to the closest cluster*/
int update_centroids(double **centroids,double **nowcentroids,int *clusters,int k,int d);/*update the centroids and return the maximum euclidean norm*/

#endif