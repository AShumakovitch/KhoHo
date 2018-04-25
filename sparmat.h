/*
 *    sparmat.h --- computation library for working with sparse matrices.
 *
 * Copyright (C) 2002-2007 Alexander Shumakovitch <Shurik@gwu.edu>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program  is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; see COPYING.gz. If not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/*
 * Types for indices and values of matrix entries.
 * 4 bytes (i.e. up to 2^32) is enough for the moment.
 */
typedef int SM_index_t;
typedef int SM_value_t;

/*
 * Function to compute the absolute value (depends on the type of entries)
 */
#define ABSFUNC abs

/*
 * Maximal entry allowed (in absolute value).
 */
#define ENTRY_MAX (INT_MAX / 2)

/*
 * Maximal possible entry value.
 */
#define ERR_MVAL INT_MAX

extern char *ERR_MESSAGE;

/*
 * Entry of a sparse vector. Members are:
 *   index of the entry in the vector (starting with 1),
 *   pointer to the next entry (NULL if none).
 */
typedef struct sparse_entry {
	SM_index_t index;
	SM_value_t value;
	struct sparse_entry *next;
} SparseEntry;

/*
 * Root of a sparse vector (either a row or a column). Members are:
 *   number of (non-zero) entries, 0 ==> no entries, -1 ==> vector is deleted,
 *   pointer to the first entry (NULL if none).
 *
 * Deleted vectors are assumed to be not modifiable.
 */
typedef struct sparse_vector {
	SM_index_t num_entries;
	SparseEntry *entries;
} SparseVector;

/*
 * Root of a sparse matrix. Members are:
 *   number of rows and columns in the matrix,
 *   pointers to arrays of sparse vector roots representing rows and columns.
 *
 * Arrays must be of length num_rows and num_cols, respectively.
 */
typedef struct sparse_matrix {
	SM_index_t num_rows, num_cols;
	SparseVector *rows, *columns;
} SparseMatrix;

/*
 * Return index of the first invertible entry in a sparse vector (that is,
 * the one whose absolute value equals 1) or 0 if no such entry exists.
 * If val is not NULL, use it to store the value of the entry deleted.
 * Return -1 and set ERR_MESSAGE if the vector is already deleted.
 */
SM_index_t find_v_unit(SparseVector *vec, SM_value_t *val);

/*
 * Return the entry's value from a sparse matrix (or 0 if there is none).
 * If row and column entries differ, set an error message and return ERR_MVAL.
 */
SM_value_t get_m_entry(SparseMatrix *matr, SM_index_t row, SM_index_t col);

/*
 * Remove entry given by its indices from a sparse matrix.
 * Return the value of the entry deleted (or 0 if there is none).
 * If row and column entries differ, set an error message and return ERR_MVAL.
 */
SM_value_t remove_m_entry(SparseMatrix *matr, SM_index_t row, SM_index_t col);

/*
 * Add entry given by its indices and value to a sparse matrix.
 * Return 0 on success and -1 otherwise.
 */
int add_m_entry(SparseMatrix *matr,
		SM_index_t row, SM_index_t col, SM_value_t val);

/*
 * Erase all entries in a given row of a sparse matrix.
 * If do_del is set, mark the corresponding sparse vector as being deleted.
 * Return 0 on success and -1 otherwise.
 */
int erase_m_row(SparseMatrix *matr, SM_index_t row, int do_del);

/*
 * Erase all entries in a given column of a sparse matrix.
 * If do_del is set, mark the corresponding sparse vector as being deleted.
 * Return 0 on success and -1 otherwise.
 */
int erase_m_column(SparseMatrix *matr, SM_index_t col, int do_del);

/*
 * Add (with a scalar) two rows in a given sparse matrix.
 * Return the maximal absolute value of the new entries and -1 on failure.
 */
SM_value_t add_m_rows(SparseMatrix *matr,
		SM_index_t row1, SM_index_t row2, SM_value_t scalar);

/*
 * Add (with a scalar) two columns in a given sparse matrix.
 * Return the maximal absolute value of the new entries and -1 on failure.
 */
SM_value_t add_m_cols(SparseMatrix *matr,
		SM_index_t col1, SM_index_t col2, SM_value_t scalar);

/*
 * Allocate memory for rows and columns of a new sparse matrix.
 * Return 0 on success and -1 otherwise.
 */
int init_s_matrix(SparseMatrix *matr,
		SM_index_t n_rows, SM_index_t n_cols);

/*
 * Free the memory allocated for a sparse matrix.
 */
void kill_s_matrix(SparseMatrix *matr);

/*
 * For testing only: print the content of a sparse vector to stdout.
 * Deleted vectors are ignored.
 */
int print_s_vector(SparseVector *vec);

/*
 * For testing only: print the content of a sparse matrix to stdout.
 */
int print_s_matrix(SparseMatrix *matr);

/*
 * For testing only: check that the vector data are consistent.
 * If the third argument is not NULL, compare the content with
 * the ``orthogonal'' family of vectors.
 */
int check_v_data(SparseVector *vec, SM_index_t max_index, SM_index_t v_ind,
							SparseVector *others);

/*
 * For testing only: check that the matrix data are consistent.
 */
int check_m_data(SparseMatrix *matr);
