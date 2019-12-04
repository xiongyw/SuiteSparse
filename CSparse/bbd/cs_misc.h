#ifndef __CS_MISC__
#define __CS_MISC__


/* read mtx format file into triplet-form
 * return NULL for failure
 */
cs* CS_load_mm(const char* path);

void CS_dump(const cs* M, const char* text);

/*
 * convert from CSC to triplet-form
 * return NULL for failure (or A is already triplet)
 */
cs* CS_unzip(const cs* A);

/*
 * make an indentical copy of a matrix
 */
cs* CS_copy(const cs* A);

/*
 * return A-B; a wrapper of cs_add()
 */
cs* CS_diff(const cs* A, const cs* B);

/*
 * allocate an empty matrix
 * triplet: 0 for CSC, otherwise triplet-form
 */
cs* CS_skeleton(csi m, csi n, int triplet);

/*
 * get an identity matrix, compressed-form
 * n: dimension
 */
cs* CS_id(csi n);

/*
 * return a submatrix (compressed) from matrix A
 * ri: the row index (in A) for the 1st row
 * ci: the col index (in A) for the 1st col
 * m: # rows
 * n: # cols
 *
 * return NULL for failures
 */
cs* CS_clip(cs* A, csi ri, csi ci, csi m, csi n);


/*
 * patch a submatrix of A: A += k * P. both A & P should be compressed-form.
 * it's required that dim of A is bigger enough to accommendate P.
 *
 * A: input matrix
 * P: input submatrix
 * ri, ci: row and col index (in A) for P[0][0]
 * k: a scalar
 *
 * return NULL for failures
 */
cs* CS_patch(const cs* A, const cs* P, csi ri, csi ci, double k);

/*
 * invert a triangular matrix
 *
 * G: triangular matrix
 * lo: L or U
 *
 * return NULL for failures
 */
cs* CS_tri_inv(cs* t, int low);

#endif // __CS_MISC__
