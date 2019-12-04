#include <assert.h>
#include <string.h>

#include "cs.h"
#include "cs_misc.h"

///////////////////////////////////////////////////////////////////////////////
// a global test 3x3 matrix: A = LU
// | 2  0  3 |   | 1               || 2  0    3   |
// | 3  6  0 | = | 1.5   1         ||    6   -4.5 |
// | 0  1  4 |   | 0     1.667   1 ||         4.75|

// triplet form
csi    g_p[] = {0,  0,  0,  1,  1,  1,  2,  2,  2};   // col idx
csi    g_i[] = {0,  1,  2,  0,  1,  2,  0,  1,  2};   // row idx
double g_x[] = {2,  3,  0,  0,  6,  1,  3,  0,  4};
cs     g_m = {.nzmax = 9, .nz = 9, .m=3, .n=3, .p = g_p, .i = g_i, .x = g_x};
///////////////////////////////////////////////////////////////////////////////


//const char g_path[] = "../../KLU/Matrix/impcol_a.mtx"; // pivoting failed
//const char g_path[] = "./fpga_dcop_01.mtx";
//const char g_path[] = "./rajat19.mtx"; // https://www.cise.ufl.edu/research/sparse/matrices/Rajat/rajat19.html
const char g_path[] = "../../SPQR/Matrix/lfat5b.mtx";


int test(void);

int main(void)
{
    return test();
}

int test(void)
{

    cs* A = NULL;

#if (1)
    // triplet-form
    cs* T = CS_load_mm(g_path);
    if (!T) {
        printf("CS_load_mm() failed\n");
        exit(1);
    }
#else
    cs* T = CS_copy(&g_m);
    //CS_dump(T, "copy");
#endif

    // test CS_id()
    cs* id = CS_id(5);
    CS_dump(id, "id");

    // test CS_unzip()
    cs* A1 = cs_compress(T) ;
    cs* A2 = CS_unzip(A1);
    cs* A3 = cs_compress(A2) ;


    // test CS_clip() and CS_patch()
    cs* A4 = CS_skeleton(A3->m, A3->n, 0);
    cs* A5 = CS_patch(A4, A3, 0, 0, 1.);


#if (1)
    cs* D = CS_clip(A5, 0, 0, 2, 2);
    CS_dump(D, "D");
    printf("-----------\n");
    cs* E = CS_clip(A5, 2, 2, A5->m-2, A5->n-2);
    cs* R = CS_clip(A5, 2, 0, A5->m-2, 2);
    cs* C = CS_clip(A5, 0, 2, 2, A5->n-2);
    cs* A6 = CS_skeleton(A5->m, A5->n, 0);
    cs* A7 = CS_patch(A6, D, 0, 0, 1.);
    cs* A8 = CS_patch(A7, E, 2, 2, 1.);
    cs* A9 = CS_patch(A8, R, 2, 0, 1.);
    cs* A10 = CS_patch(A9, C, 0, 2, 1.);
#endif

    A = A10;
    // symbolic LU ordering and analysis
    int order = 0; // 0: natural ordering, 1: A+A'
    int is_qr = 0; // lu decomp
    css* sym = cs_sqr(order, A, is_qr);

    if (!sym) {
        printf("ERROR: cs_sqr() failed\n");
        exit(1);
    }

    if (sym->pinv) {
        printf("sym->pinv != NULL\n");
    } else {
        printf("sym->pinv = NULL\n");
    }

    if (sym->q) {
        printf("sym->q[] = [");
        for (int i = 0; i < A->n; ++i) {
            printf("%ld ", sym->q[i]);
        }
        printf("]\n");
    } else {
        printf("sym->q = NULL\n");
    }

    if (sym->parent) {
        printf("sym->parent != NULL\n");
    } else {
        printf("sym->parent = NULL\n");
    }

    if (sym->cp) {
        printf("sym->cp != NULL\n");
    } else {
        printf("sym->cp = NULL\n");
    }

    if (sym->leftmost) {
        printf("sym->leftmost != NULL\n");
    } else {
        printf("sym->leftmost = NULL\n");
    }

    printf("sym->m2=%ld\n", sym->m2);
    printf("sym->lnz=%lf\nsym->unz=%lf\n", sym->lnz, sym->unz);

    // LU factor
    double tol = 0.;
    csn* num = cs_lu(A, sym, tol); // sym can't be NULL

    if (!num) {
        printf("ERROR: cs_lu() failed\n");
        exit(1);
    }

    // verify that L*U=A
    cs* LU = cs_multiply(num->L, num->U);
    assert(LU && "cs_multiply() failure");
    cs* diff = CS_diff(A, LU);
    assert(diff && "cs_multiply() failure");
    printf("|A-LU|=%g\n", cs_norm(diff));
    cs_spfree(LU);
    cs_spfree(diff);

    // inverse of L
    cs* id2 = CS_id(num->L->n);
    cs* Linv = CS_tri_inv(num->L, 1);
    //CS_dump(Linv, "Linv");
    //cs_print(Linv, 1);
    cs* mul = cs_multiply(num->L, Linv);
    //cs_print(mul, 1);
    cs* Ldiff= cs_add(id2, mul, 1, -1);
    printf("|I-L*Linv|=%g\n", cs_norm(Ldiff));

    // inverse of U
    cs* Uinv = CS_tri_inv(num->U, 0);
    cs* Umul = cs_multiply(num->U, Uinv);
    cs* Udiff= cs_add(id2, Umul, 1, -1);
    printf("|I-U*Uinv|=%g\n", cs_norm(Udiff));

#if (0)
    // L/U
    int brief = 0;
    CS_dump(num->L, "L");
    cs_print(num->L, brief);
    CS_dump(num->U, "U");
    cs_print(num->U, brief);

    if (num->pinv) {
        printf("num->pinv[]=[");
        for (int i = 0; i < A->m; ++ i) printf("%ld, ", num->pinv[i]);
        printf("]\n");
    } else {
        printf("num->pinv = NULL\n");
    }
#endif

    // tear down
    cs_nfree(num);
    cs_sfree(sym);
    cs_spfree (A) ;
    cs_spfree (T) ;

    return 0;
}

