/*
 *    sparmat.c --- computation library for working with sparse matrices.
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

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "sparmat.h"

/* perform some consistency tests */
#define SPARMAT_DEBUG

char *ERR_MESSAGE;
#define ERR_RET(msg, val) { ERR_MESSAGE = (msg); return (val); }
#define ERRET_1(msg) ERR_RET((msg), -1)
#define ERRET_M(msg) ERR_RET((msg), ERR_MVAL)

/*
 * Return the entry's value from a sparse vector (or 0 if there is none).
 */
static SM_value_t get_v_entry(SparseVector *vec, SM_index_t ind)
{
	SparseEntry *eptr;

	/* find the desired entry (if it exists) */
	for (eptr = vec->entries; eptr != NULL; eptr = eptr->next)
		if (eptr->index >= ind) break;

	if (eptr != NULL && eptr->index == ind) return eptr->value;

	return 0;
}

/*
 * Return index of the first invertible entry in a sparse vector (that is,
 * the one whose absolute value equals 1) or 0 if no such entry exists.
 * If val is not NULL, use it to store the value of the entry deleted.
 * Return -1 and set ERR_MESSAGE if the vector is already deleted.
 */
SM_index_t find_v_unit(SparseVector *vec, SM_value_t *val)
{
	SparseEntry *eptr;

	if (vec->num_entries == -1)
		ERRET_1("find_v_unit: vector is already deleted");

	/* find the desired entry (if it exists) */
	for (eptr = vec->entries; eptr != NULL; eptr = eptr->next) {
		if (ABSFUNC(eptr->value) == 1) {
			if (val != NULL) *val = eptr->value;
			return eptr->index;
		}
	}

	return 0;
}

/*
 * Remove an entry given by its index from a sparse vector.
 * Return the value of the entry deleted (or 0 if there is none).
 * Return ERR_MVAL and set ERR_MESSAGE if the vector is already deleted.
 */
static SM_value_t remove_v_entry(SparseVector *vec, SM_index_t ind)
{
	SM_value_t val;
	SparseEntry *eptr, *prev = NULL;

	if (vec->num_entries == -1)
		ERRET_M("remove_v_entry: vector is already deleted");

	/* find the entry to remove (if it exists) */
	for (eptr = vec->entries; eptr != NULL; prev = eptr, eptr = eptr->next)
		if (eptr->index >= ind) break;

	if (eptr != NULL && eptr->index == ind) {
		val = eptr->value;
		if (prev != NULL)
			prev->next = eptr->next;
		else
			vec->entries = eptr->next;
		free(eptr);

		/* do we want to check that the number of entries >= 0 ?? */
		vec->num_entries--;
	} else return 0;

	return val;
}

/*
 * Add an entry given by its index and value to a sparse vector.
 * Return 0 on success and -1 otherwise (also if the vector is deleted).
 */
static int add_v_entry(SparseVector *vec, SM_index_t ind, SM_value_t val)
{
	SparseEntry *eptr, *prev = NULL;

	if (vec->num_entries == -1)
		ERRET_1("add_v_entry: vector is already deleted");

	/* zero entries don't exist. Not having an entry to remove is OK */
	if (val == 0) { remove_v_entry(vec, ind); return 0; }

	/* find an entry in the vector the new one should be added after */
	for (eptr = vec->entries; eptr != NULL; prev = eptr, eptr = eptr->next)
		if (eptr->index >= ind) break;

	if (eptr != NULL && eptr->index == ind)
		/* the entry already exists, only change its value than */
		eptr->value = val;
	else {
		if ((eptr = malloc(sizeof(SparseEntry))) == NULL)
			ERRET_1("add_v_entry: not enough memory");

		eptr->index = ind;
		eptr->value = val;

		if (prev != NULL) {
			eptr->next = prev->next;
			prev->next = eptr;
		} else {
			eptr->next = vec->entries;
			vec->entries = eptr;
		}

		/* do we want to check the number of entries is not too big? */
		vec->num_entries++;
	}

	return 0;
}

/*
 * Check that row and column indices are within defined boundaries.
 */
static int check_m_indices(SparseMatrix *matr, SM_index_t row, SM_index_t col)
{
	if (row < 1 || row > matr->num_rows || col < 1 || col > matr->num_cols)
		ERRET_1("check_m_indices: wrong matrix indices");

	return 0;
}

/*
 * Return the entry's value from a sparse matrix (or 0 if there is none).
 * If row and column entries differ, set an error message and return ERR_MVAL.
 */
SM_value_t get_m_entry(SparseMatrix *matr, SM_index_t row, SM_index_t col)
{
	SM_value_t val;

	if (check_m_indices(matr, row, col) == -1) return ERR_MVAL;

	val = get_v_entry(matr->rows + row - 1, col);

#ifdef SPARMAT_DEBUG
	if (val != get_v_entry(matr->columns + col - 1, row))
		ERRET_M("get_m_entry: row and column entries don't match");
#endif

	return val;
}

/*
 * Remove an entry given by its indices from a sparse matrix.
 * Return the value of the entry deleted (or 0 if there is none).
 * If row and column entries differ, set an error message and return ERR_MVAL.
 */
SM_value_t remove_m_entry(SparseMatrix *matr, SM_index_t row, SM_index_t col)
{
	SM_value_t valr, valc;

	if (check_m_indices(matr, row, col) == -1) return ERR_MVAL;

	/* the row could be already deleted */
	if ((valr = remove_v_entry(matr->rows + row - 1, col)) == ERR_MVAL)
		return ERR_MVAL;

	/* the column could be already deleted */
	if ((valc = remove_v_entry(matr->columns + col - 1, row)) == ERR_MVAL)
		return ERR_MVAL;

	if (valr != valc)
		ERRET_M("remove_m_entry: row and column entries don't match");

	return valr;
}

/*
 * Add entry given by its indices and value to a sparse matrix.
 * Return 0 on success and -1 otherwise.
 */
int add_m_entry(SparseMatrix *matr,
		SM_index_t row, SM_index_t col, SM_value_t val)
{
	if (check_m_indices(matr, row, col) == -1) return -1;

	if (val > ENTRY_MAX || val < -ENTRY_MAX)
		ERRET_1("add_m_entry: entry's value is too big");

	/* zero entries don't exist. Not having an entry to remove is OK */
	if (val == 0) {
		if (remove_m_entry(matr, row, col) == ERR_MVAL) return -1;
		else return 0;
	}

	if (add_v_entry(matr->rows + row - 1, col, val) == -1) return -1;

	return add_v_entry(matr->columns + col - 1, row, val);
}

/*
 * Erase all entries in a given row or column of a sparse matrix and clean up
 * the corresponding entries in the ``orthogonal'' family of vectors.
 * If do_del is set, mark the sparse vector as being deleted.
 * Return 0 on success and -1 otherwise.
 */
static int erase_m_colrow(SparseVector *cr_vec, SM_index_t cr_ind,
					SparseVector *others, int do_del)
{
	SM_value_t val;
	SparseEntry *eptr = cr_vec->entries;

	if (cr_vec->num_entries == -1)
		ERRET_1("erase_m_colrow: vector is already deleted");

	while (eptr != NULL) {
		val = remove_v_entry(others + eptr->index - 1, cr_ind);
		if (val == ERR_MVAL) return -1;

#ifdef SPARMAT_DEBUG
		if (val != eptr->value)
			ERRET_1 \
			("erase_m_colrow: row and column entries don't match");
#endif

		/* this way the vector stays always sane */
		cr_vec->entries = eptr->next;
		free(eptr);
		cr_vec->num_entries--;
		eptr = cr_vec->entries;
	}

#ifdef SPARMAT_DEBUG
	if (cr_vec->num_entries != 0)
		ERRET_1("erase_m_colrow: num_entries is not 0 after erasing");
#endif

	if (do_del) cr_vec->num_entries = -1;
	return 0;
}

/*
 * Erase all entries in a given row of a sparse matrix.
 * If do_del is set, mark the corresponding sparse vector as being deleted.
 * Return 0 on success and -1 otherwise.
 */
int erase_m_row(SparseMatrix *matr, SM_index_t row, int do_del)
{
	if (check_m_indices(matr, row, 1) == -1) return -1;

	return erase_m_colrow(matr->rows + row - 1, row, matr->columns, do_del);
}

/*
 * Erase all entries in a given column of a sparse matrix.
 * If do_del is set, mark the corresponding sparse vector as being deleted.
 * Return 0 on success and -1 otherwise.
 */
int erase_m_column(SparseMatrix *matr, SM_index_t col, int do_del)
{
	if (check_m_indices(matr, 1, col) == -1) return -1;

	return erase_m_colrow(matr->columns + col - 1, col, matr->rows, do_del);
}

/*
 * Add (with a scalar) two rows or columns of a sparse matrix, and assign the
 * corresponding entries in the ``orthogonal'' family of vectors appropriately.
 * Return the maximal absolute value of the new entries and -1 on failure.
 */
static SM_value_t add_m_colrows(SparseVector *cr_vec1, SM_index_t cr_ind1,
	SparseVector *cr_vec2, SparseVector *others, SM_index_t scalar)
{
	SM_value_t maxval = 0;
	SparseEntry *prev = NULL, *new, *old;
	SparseEntry *eptr1 = cr_vec1->entries, *eptr2 = cr_vec2->entries;

	if (cr_vec1->num_entries == -1 || cr_vec2->num_entries == -1)
		ERRET_1("erase_m_colrow: vector is already deleted");

	while (eptr2 != NULL) {
		/* there is an unmatched entries in the first vector */
		if (eptr1 != NULL && eptr1->index < eptr2->index) {
			prev = eptr1;
			eptr1 = eptr1->next;
			continue;
		}

		/* a trick to deal with newly created zero entries */
		old = prev;

		/* there is an unmatched entry in the second vector */
		if (eptr1 == NULL || eptr1->index > eptr2->index) {
			if ((new = malloc(sizeof(SparseEntry))) == NULL)
				ERRET_1("add_m_colrows: not enough memory");

			new->index = eptr2->index;
			/* need to check for admissible values here !!! */
			new->value = scalar * eptr2->value;
			new->next = eptr1;
			if (prev != NULL) {
				prev->next = new;
			} else {
				cr_vec1->entries = new;
			}
			cr_vec1->num_entries++;

			if (ABSFUNC(new->value) > maxval)
				maxval = ABSFUNC(new->value);

			prev = new;
			eptr2 = eptr2->next;
		} else {
			/* entries are matched: both indices are the same */
			eptr1->value += scalar * eptr2->value;

			if (ABSFUNC(eptr1->value) > maxval)
				maxval = ABSFUNC(eptr1->value);

			prev = eptr1;
			eptr1 = eptr1->next;
			eptr2 = eptr2->next;
		}
		if (ABSFUNC(prev->value) > ENTRY_MAX)
			ERRET_1("add_m_colrows: entry's value is too big");

		if (add_v_entry(others + prev->index - 1,
					cr_ind1, prev->value) == -1)
			return -1;

		/* if the new entry's value is 0, remove the entry */
		if (prev->value == 0) {
			free(prev);
			prev = old;

			if (prev != NULL)
				prev->next = eptr1;
			else
				cr_vec1->entries = eptr1;
			cr_vec1->num_entries--;
		}
	}

	return maxval;
}

/*
 * Add (with a scalar) two rows in a given sparse matrix.
 * Return the maximal absolute value of the new entries and -1 on failure.
 */
SM_value_t add_m_rows(SparseMatrix *matr,
		SM_index_t row1, SM_index_t row2, SM_value_t scalar)
{
	if (check_m_indices(matr, row1, 1) == -1) return -1;
	if (check_m_indices(matr, row2, 1) == -1) return -1;

	return add_m_colrows(matr->rows + row1 - 1, row1,
				matr->rows + row2 - 1, matr->columns, scalar);
}

/*
 * Add (with a scalar) two columns in a given sparse matrix.
 * Return the maximal absolute value of the new entries and -1 on failure.
 */
SM_value_t add_m_cols(SparseMatrix *matr,
		SM_index_t col1, SM_index_t col2, SM_value_t scalar)
{
	if (check_m_indices(matr, 1, col1) == -1) return -1;
	if (check_m_indices(matr, 1, col2) == -1) return -1;

	return add_m_colrows(matr->columns + col1 - 1, col1,
				matr->columns + col2 - 1, matr->rows, scalar);
}

/*
 * Allocate memory for rows and columns of a new sparse matrix.
 * Return 0 on success and -1 otherwise.
 */
int init_s_matrix(SparseMatrix *matr, SM_index_t n_rows, SM_index_t n_cols)
{
	SM_index_t i;
	SparseVector *rvec, *cvec, *vptr;

	if (n_rows < 1 || n_cols < 1)
		ERRET_1("init_s_matrix: number of rows or columns is too small");

	if ((rvec = (SparseVector *) malloc(n_rows * sizeof(SparseVector)))
			== NULL)
		ERRET_1("init_s_matrix: not enough memory");

	if ((cvec = (SparseVector *) malloc(n_cols * sizeof(SparseVector)))
			== NULL) {
		free(rvec);
		ERRET_1("init_s_matrix: not enough memory");
	}

	matr->num_rows = n_rows;
	matr->num_cols = n_cols;
	matr->rows = rvec;
	matr->columns = cvec;

	for (i = 0, vptr = rvec; i < n_rows; i++, vptr++) {
		vptr->num_entries = 0;
		vptr->entries = NULL;
	}
	for (i = 0, vptr = cvec; i < n_cols; i++, vptr++) {
		vptr->num_entries = 0;
		vptr->entries = NULL;
	}

	return 0;
}

/*
 * Free the memory allocated for a sparse vector (but don't touch the root).
 */
static inline void kill_s_vector(SparseVector *vec)
{
	SparseEntry *eptr = vec->entries;

	while (eptr != NULL) {
		/* this way the vector stays always sane */
		vec->entries = eptr->next;
		free(eptr);
		vec->num_entries--;
		eptr = vec->entries;
	}
}

/*
 * Free the memory allocated for a sparse matrix.
 */
void kill_s_matrix(SparseMatrix *matr)
{
	SM_index_t i;
	SparseVector *vptr;

	vptr = matr->rows;
	/* check for NULL in case we kill the matrix before initialization */
	if (vptr != NULL) {
		for (i = 0; i < matr->num_rows; i++) kill_s_vector(vptr++);
			free(matr->rows);
	}

	vptr = matr->columns;
	if (vptr != NULL) {
		for (i = 0; i < matr->num_cols; i++) kill_s_vector(vptr++);
			free(matr->columns);
	}
}

/*
 * For testing only: print the content of a sparse vector to stdout.
 * Deleted vectors are ignored.
 */
int print_s_vector(SparseVector *vec)
{
	SM_index_t i, num_e = vec->num_entries;
	SparseEntry *eptr = vec->entries;

	if (num_e == -1) {
		printf("vector is deleted\n");
		return 0;
	}

	printf("%d entries: ", num_e);
	for (i = 0; i < num_e; i++, eptr = eptr->next) {
		if (eptr == NULL) ERRET_1("print_vector: vector is corrupt");

		if (i) printf("; ");
		printf("%d, %d", eptr->index, eptr->value);
	}

	printf(".\n");
	return 0;
}

/*
 * For testing only: print the content of a sparse matrix to stdout.
 */
int print_s_matrix(SparseMatrix *matr)
{
	SM_index_t i;
	SparseVector *vec;

	printf("%d rows and %d columns:\n", matr->num_rows, matr->num_cols);
	for (i = 0, vec = matr->rows; i < matr->num_rows; i++, vec++) {
		printf("  The row number %d, ", i + 1);
		if (print_s_vector(vec) == -1) return -1;
	}
	printf("\n");
	for (i = 0, vec = matr->columns; i < matr->num_cols; i++, vec++) {
		printf("  The column number %d, ", i + 1);
		if (print_s_vector(vec) == -1) return -1;
	}

	return 0;
}

/*
 * For testing only: check that the vector data are consistent.
 * If the third argument is not NULL, compare the content with
 * the ``orthogonal'' family of vectors.
 */
int check_v_data(SparseVector *vec, SM_index_t max_index, SM_index_t v_ind,
							SparseVector *others)
{
	SM_index_t i, oldind = -1, n_entries;
	SM_value_t val;
	SparseEntry *eptr;

	n_entries = vec->num_entries;
	if (n_entries > max_index)
		ERRET_1("check_v_data: number of entries is too big");
	if (n_entries < -1 )
		ERRET_1("check_v_data: number of entries is negative");
	if (n_entries == -1 && vec->entries != NULL)
		ERRET_1("check_v_data: deleted vector is not empty");

	if (n_entries == -1) n_entries = 0;

	for (eptr = vec->entries, i = 1; eptr != NULL; eptr = eptr->next, i++) {
		if (eptr->index < 1)
			ERRET_1("check_v_data: index is not positive");
		if (eptr->index > max_index)
			ERRET_1("check_v_data: index is too big");
		if (eptr->index <= oldind)
			ERRET_1("check_v_data: index is not increasing");

		if (eptr->value == 0)
			ERRET_1("check_v_data: value is 0");

		oldind = eptr->index;
		if (others == NULL) continue;

		val = get_v_entry(others + eptr->index - 1, v_ind);
		if (val != eptr->value)
			ERRET_1("check_v_data: rows and columns don't match");

	}
	if (i != n_entries + 1)
		ERRET_1("check_v_data: wrong number of entries");

	return 0;
}

/*
 * For testing only: check that the matrix data are consistent.
 */
int check_m_data(SparseMatrix *matr)
{
	SM_index_t i;
	SM_index_t n_rows = matr->num_rows, n_cols = matr->num_cols;
	SparseVector *vec;

	for (i = 0, vec = matr->rows; i < n_rows; i++, vec++)
		if (check_v_data(vec, n_rows, i + 1, matr->columns) == -1)
			return -1;

	for (i = 0, vec = matr->columns; i < n_cols; i++, vec++)
		if (check_v_data(vec, n_rows, i + 1, matr->rows) == -1)
			return -1;

	return 0;
}
