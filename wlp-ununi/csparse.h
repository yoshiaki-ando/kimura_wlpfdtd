#include <stdio.h>
#include <stdlib.h>

#ifndef _CS_H
#define _CS_H

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#define CS_VER 1		    /* CSparse Version 1.2.0 */
#define CS_SUBVER 2
#define CS_SUBSUB 0
#define CS_DATE "Mar 6, 2006"	    /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006"

/* --- primary CSparse routines and data structures ------------------------- */
typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    long long nzmax ;	    /* maximum number of entries */
    long long m ;	    /* number of rows */
    long long n ;	    /* number of columns */
    long long *p ;	    /* column pointers (size n+1) or col indices (size nzmax) */
    long long *i ;	    /* row indices, size nzmax */
    double *x ;	    /* numerical values, size nzmax */
    long long nz ;	    /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

cs *cs_add (const cs *A, const cs *B, double alpha, double beta) ;
long long cs_cholsol (const cs *A, double *b, long long order) ;
long long cs_dupl (cs *A) ;
long long cs_entry (cs *T, long long i, long long j, double x) ;
long long cs_lusol (const cs *A, double *b, long long order, double tol) ;
long long cs_gaxpy (const cs *A, const double *x, double *y) ;
cs *cs_multiply (const cs *A, const cs *B) ;
long long cs_qrsol (const cs *A, double *b, long long order) ;
cs *cs_transpose (const cs *A, long long values) ;
cs *cs_triplet (const cs *T) ;
double cs_norm (const cs *A) ;
long long cs_print (const cs *A, long long brief) ;
cs *cs_load (FILE *f) ;
/* utilities */
void *cs_calloc (long long n, size_t size) ;
void *cs_free (void *p) ;
void *cs_realloc (void *p, long long n, size_t size, long long *ok) ;
cs *cs_spalloc (long long m, long long n, long long nzmax, long long values, long long triplet) ;
cs *cs_spfree (cs *A) ;
long long cs_sprealloc (cs *A, long long nzmax) ;
void *cs_malloc (long long n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */
typedef struct cs_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    long long *Pinv ;	    /* inverse row perm. for QR, fill red. perm for Chol */
    long long *Q ;	    /* fill-reducing column permutation for LU and QR */
    long long *parent ;   /* elimination tree for Cholesky and QR */
    long long *cp ;	    /* column pointers for Cholesky, row counts for QR */
    long long m2 ;	    /* # of rows for QR, after adding fictitious rows */
    long long lnz ;	    /* # entries in L for LU or Cholesky; in V for QR */
    long long unz ;	    /* # entries in U for LU; in R for QR */
} css ;

typedef struct cs_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    cs *L ;	    /* L for LU and Cholesky, V for QR */
    cs *U ;	    /* U for LU, R for QR, not used for Cholesky */
    long long *Pinv ;	    /* partial pivoting for LU */
    double *B ;	    /* beta [0..n-1] for QR */
} csn ;

typedef struct cs_dmperm_results    /* cs_dmperm or cs_scc output */
{
    long long *P ;	    /* size m, row permutation */
    long long *Q ;	    /* size n, column permutation */
    long long *R ;	    /* size nb+1, block k is rows R[k] to R[k+1]-1 in A(P,Q) */
    long long *S ;	    /* size nb+1, block k is cols S[k] to S[k+1]-1 in A(P,Q) */
    long long nb ;	    /* # of blocks in fine dmperm decomposition */
    long long rr [5] ;    /* coarse row decomposition */
    long long cc [5] ;    /* coarse column decomposition */
} csd ;

long long *cs_amd (const cs *A, long long order) ;
csn *cs_chol (const cs *A, const css *S) ;
csd *cs_dmperm (const cs *A) ;
long long cs_droptol (cs *A, double tol) ;
long long cs_dropzeros (cs *A) ;
long long cs_happly (const cs *V, long long i, double beta, double *x) ;
long long cs_ipvec (long long n, const long long *P, const double *b, double *x) ;
long long cs_lsolve (const cs *L, double *x) ;
long long cs_ltsolve (const cs *L, double *x) ;
csn *cs_lu (const cs *A, const css *S, double tol) ;
cs *cs_permute (const cs *A, const long long *P, const long long *Q, long long values) ;
long long *cs_pinv (const long long *P, long long n) ;
long long cs_pvec (long long n, const long long *P, const double *b, double *x) ;
csn *cs_qr (const cs *A, const css *S) ;
css *cs_schol (const cs *A, long long order) ;
css *cs_sqr (const cs *A, long long order, long long qr) ;
cs *cs_symperm (const cs *A, const long long *Pinv, long long values) ;
long long cs_usolve (const cs *U, double *x) ;
long long cs_utsolve (const cs *U, double *x) ;
long long cs_updown (cs *L, long long sigma, const cs *C, const long long *parent) ;
/* utilities */
css *cs_sfree (css *S) ;
csn *cs_nfree (csn *N) ;
csd *cs_dfree (csd *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */
long long *cs_counts (const cs *A, const long long *parent, const long long *post, long long ata) ;
long long cs_cumsum (long long *p, long long *c, long long n) ;
long long cs_dfs (long long j, cs *L, long long top, long long *xi, long long *pstack, const long long *Pinv) ;
long long *cs_etree (const cs *A, long long ata) ;
long long cs_fkeep (cs *A, long long (*fkeep) (long long, long long, double, void *), void *other) ;
double cs_house (double *x, double *beta, long long n) ;
long long *cs_maxtrans (const cs *A) ;
long long *cs_post (long long n, const long long *parent) ;
long long cs_reach (cs *L, const cs *B, long long k, long long *xi, const long long *Pinv) ;
csd *cs_scc (cs *A) ;
long long cs_scatter (const cs *A, long long j, double beta, long long *w, double *x, long long mark,
    cs *C, long long nz) ;
long long cs_splsolve (cs *L, const cs *B, long long k, long long *xi, double *x,
    const long long *Pinv) ;
long long cs_tdfs (long long j, long long k, long long *head, const long long *next, long long *post,
    long long *stack) ;
/* utilities */
csd *cs_dalloc (long long m, long long n) ;
cs *cs_done (cs *C, void *w, void *x, long long ok) ;
long long *cs_idone (long long *p, cs *C, void *w, long long ok) ;
csn *cs_ndone (csn *N, cs *C, void *w, void *x, long long ok) ;
csd *cs_ddone (csd *D, cs *C, void *w, long long ok) ;

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > LONG_MAX / (long long) size)
#endif
