#include <assert.h>
#include <string.h>

#include "cs.h"
#include "mmio.h"
#include "cs_misc.h"

#define CS_MISC_DEBUG  1

// cf: https://math.nist.gov/MatrixMarket/mmio/c/example_read.c
cs* CS_load_mm(const char* path)
{
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    csi i, *I, *J;
    double *val;
    int tmp;

    if ((f = fopen(path, "r")) == NULL) {
        printf("file %s not founc!", path);
        return NULL;
    }

    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        return NULL;
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        return NULL;
    }

    /* find out size of sparse matrix .... */

    if (0 != mm_read_mtx_crd_size(f, &M, &N, &nz)) {
        printf("Sorry, mm_read_mtx_crd_size() failed!");
        return NULL;
    }

    printf("(M, N, nz)=(%d, %d, %d)\n", M, N, nz);

    /* reseve memory for matrices */
    I = (csi*) malloc(nz * sizeof(csi));
    J = (csi*) malloc(nz * sizeof(csi));
    val = (double *) malloc(nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i=0; i<nz; i++)
    {
        tmp = fscanf(f, "%ld %ld %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    fclose(f);

    cs* A = (cs*)malloc(sizeof(cs));
    A->m = M;
    A->n = N;
    A->nzmax = nz;
    A->nz = nz;  // not compressed
    A->p = J;
    A->i = I;
    A->x = val;


    return A;
}

cs* CS_unzip(const cs* A)
{
    if (!A || (A->nz != -1)) {
        assert( 0 && "input is either NULL, or it's already triplet-form");
        return NULL;
    }

    // allocate a triplet-form matrix, initially with 1 value item only
    cs* t = cs_spalloc(A->m, A->n, A->nzmax, 1, 1);
    assert(t && "cs_spalloc() failure");

    // travese the elements
    int ret;
    for (int j = 0 ; j < A->n ; ++j) {
        for (int p = A->p[j] ; p < A->p[j+1] ; p++) {
            int row = A->i[p]; // col = j
            ret = cs_entry(t, row, j, A->x[p]);
            assert((ret == 1) && "cs_entry() failure!\n");
        }
    }

    return t;
}

void CS_dump(const cs* M, const char* text)
{
    printf("\n%s (%s): (m, n, nzmax, nz) = (%ld, %ld, %ld, %ld)\n", text, (M->nz==-1)? "compressed": "triplet", M->m, M->n, M->nzmax, M->nz);
    if (M->nz != - 1) {  // triplet
        printf("p[]=[");
        for (int i = 0; i < M->nz; ++ i) printf("%ld, ", M->p[i]);
        printf("]  <- col start indices\n");

        printf("i[]=[");
        for (int i = 0; i < M->nz; ++ i) printf("%ld, ", M->i[i]);
        printf("] <- row indices\n");

        printf("x[]=[");
        for (int i = 0; i < M->nz; ++ i) printf("%lf, ", M->x[i]);
        printf("]\n");
    } else {  // compressed
        printf("p[]=[");
        for (int i = 0; i < M->n; ++ i) printf("%ld, ", M->p[i]);
        printf("] <- col start indices\n");

        printf("i[]=[");
        for (int i = 0; i < M->nzmax; ++ i) printf("%ld, ", M->i[i]);
        printf("] <- row indices\n");

        printf("x[]=[");
        for (int i = 0; i < M->nzmax; ++ i) printf("%lf, ", M->x[i]);
        printf("]\n");
    }
}

cs* CS_diff(const cs* A, const cs* B)
{
    return cs_add(A, B, 1, -1);
}

cs* CS_copy(const cs* src)
{
    assert(src);

    cs* c = (cs*)malloc(sizeof(cs));
    assert(c && "malloc() failure");

    c->nzmax = src->nzmax;
    c->nz = src->nz;
    c->m = src->m;
    c->n = src->n;

    if (src->i) {
        c->i = (csi*)calloc(src->nzmax, sizeof(csi));
        assert(c->i);
        memcpy(c->i, src->i, sizeof(csi) * src->nzmax);
    } else {
        c->i = NULL;
    }

    if (src->x) {
        c->x = (double*)calloc(src->nzmax, sizeof(double));
        assert(c->x);
        memcpy(c->x, src->x, sizeof(double) * src->nzmax);
    } else {
        c->x = NULL;
    }

    if (src->p) {
        int len = 0;  // #elements of p[]
        if (src->nz == -1) { // compressed
            len = src->n + 1;
        } else { // triplet
            len = src->nzmax;
        }

        c->p = (csi*)malloc(sizeof(csi)*len);
        assert(c->p);
        memcpy(c->p, src->p, sizeof(csi)*len);
    } else {
        c->p = NULL;
    }

    return c;
}



cs* CS_clip(cs* A, csi ri, csi ci, csi m, csi n)
{
    int ret;

    // sanity checks
    if (!A || m < 1 || n < 1 || ri < 0 || ci < 0 ||
        (ri >= A->m) || (ci >= A->n) ||
        ((ri + m) > A->m) || ((ci + n) > A->n)) {
        assert(0 && "CS_submatrix() error: wrong parameters!");
        return NULL;
    }

    cs* t = CS_skeleton(m, n, 1);
    //CS_dump(t, "clip");

    if (A->nz != -1) { // source is triplet-form
        for(int i = 0; i < A->nzmax; ++ i) {
            if(A->p[i] >= ci && A->p[i] < (ci + n) &&
               A->i[i] >= ri && A->i[i] < (ri + m)) {
                ret = cs_entry(t, A->i[i] - ri, A->p[i] - ci, A->x[i]);
                assert((ret == 1) && "cs_entry() failure!\n");
            }
        }
    } else { // source is compressed
        for (int j = 0 ; j < A->n ; ++j) {
            for (int p = A->p[j] ; p < A->p[j+1] ; p++) {
                int row = A->i[p]; // col = j
                if(row >= ri &&
                   row < (ri + m) &&
                   j >= ci &&
                   j < (ci + n)) {
                    //printf("(row, col, value)=%d, %d, %g\n", row, j, A->x[p]);
                    //printf("new (row, col)=%d, %d\n", row-ri, j - ci);
                    ret = cs_entry(t, row - ri, j - ci, A->x[p]);
                    assert((ret == 1) && "cs_entry() failure!\n");
                }
            }
        }
    }

    cs* c = cs_compress(t);
    cs_spfree(t);

    return c;
}


cs* CS_skeleton(csi m, csi n, int triplet)
{
    cs* A = cs_spalloc(m, n,
                       1,  // nzmax
                       1,  // values
                       triplet);

    assert(A && "cs_spalloc() fail");

    if (triplet) {
        A->nz = 0;
    } else {
        // set all Ap[] to zero
        for (int i = 0; i <= n; i ++)
            A->p[i] = 0;
    }



    return A;
}

cs* CS_patch(const cs* A, const cs* P, csi ri, csi ci, double k)
{
    int ret;

    if (!A || !P ||
        !CS_CSC (A) ||!CS_CSC (P) ||
        ri < 0 || ci < 0 ||
        ri >= A->m || ci >= A->n ||
        ((ri + P->m) > A->m) ||
        ((ci + P->n) > A->n)) {
        assert (0 && "input parameter wrong");
        return 0;
    }

    // in case A & P are of the same size and completely overlap
    if (ri == 0 && ci == 0 && A->m == P->m && A->n == P->n) {
        return cs_add(A, P, 1, k);
    }

    // extend P to the same size of A, in triplet-form
    cs* PP = cs_spalloc(A->m, A->n, P->nzmax,
                             1, // allocate x[]: this necessary!
                             1); // triplet-form

    // populate PP
    for (int j = 0 ; j < P->n ; ++j) {
        for (int p = P->p[j] ; p < P->p[j+1] ; p++) {
            int row = P->i[p]; // col = j
            ret = cs_entry(PP, row + ri, j + ci, P->x[p]);
            assert((ret == 1) && "cs_entry() failure!\n");
        }
    }

    // compress PP
    cs* C = cs_compress(PP);

    // patch by cs_add()
    cs* AP = cs_add(A, C, 1, k);

    cs_spfree(PP);
    cs_spfree(C);

    return AP;
}
