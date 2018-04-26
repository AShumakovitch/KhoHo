/*
 *    sparreduce-U.c --- reduce a chain complex with only free Abelian chain
 *                       groups as far as possible using a sequence of
 *                       elementary collapses and merging of cells.
 *                       Matrices of differentials are presented in the
 *                       PARI/GP implementation of a sparse format, and all
 *                       computations are done using sparmat library.
 *
 * Copyright (C) 2002--2018 Alexander Shumakovitch <Shurik@gwu.edu>
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
 *
 * To load from PARI/GP:
 *    install(reduce_s_complex_U, "LGGG", reduce_s_complex, "./sparreduce-U.so")
 */

#include <stdlib.h>
#include <limits.h>
#include <pari/pari.h>

#if PARI_VERSION_CODE > PARI_VERSION(2,7,0)
#  define talker e_MISC
#endif

#include "sparmat-U.h"

/* print some reduction statistic */
// #define PRINT_REDSTAT

/* print some debugging messages */
// #define PRINT_DEBUG

#define ERR_BAIL(msg) { ERR_MESSAGE = (msg); bailout(); }

/*
 * Type for complex sizes. 4 bytes (i.e. up to 2^32) is enough.
 */
typedef int SM_complex_t;

static SM_complex_t cplx_size = 0;		// size of the chain complex
static SM_complex_t first_group = 0;		// first non-empty chain group
static SM_complex_t last_group = 0;		// last non-empty chain group
static SparseMatrix *cplx_matrices = NULL;	// matrices of differentials
static SM_index_t *cplx_group_ranks = NULL;	// ranks of the chain groups
static SM_index_t *num_generators = NULL;	// current number of generators
static GEN pari_matrices = NULL;		// matrices in Pari format
static GEN num_entries = NULL;  		// number of matrix entries

/*
 * Free all the memory allocated in the process.
 */
static void cleanup(void)
{
	SM_complex_t i;

	if (cplx_matrices != NULL) {
		for (i = 0; i < cplx_size - 1; i++)
			kill_s_matrix(cplx_matrices + i);
		free(cplx_matrices);
	}

	if (cplx_group_ranks != NULL) free(cplx_group_ranks);
	if (num_generators != NULL) free(num_generators);

	cplx_size = 0;
	first_group = 0;
	last_group = 0;
	cplx_matrices = NULL;
	cplx_group_ranks = NULL;
	num_generators = NULL;
	pari_matrices = NULL;
	num_entries = NULL;
}

/*
 * Bailout on error.
 */
static void bailout(void)
{
	cleanup();
	pari_err(talker, ERR_MESSAGE);
}

/*
 * Allocate memory for main arrays.
 */
static void malloc_arrays(SM_complex_t c_size)
{
	SM_complex_t i;
	char *mem_error = "malloc_arrays: not enough memory";

	cplx_size = c_size;

	cplx_matrices = (SparseMatrix *)
				malloc((cplx_size - 1) * sizeof(SparseMatrix));
	/* check for problems early, to be able to kill cplx_matrices later */
	if (cplx_matrices == NULL) ERR_BAIL(mem_error)
	for (i = 0; i < cplx_size - 1; i++) {
		cplx_matrices[i].num_rows = 0;
		cplx_matrices[i].num_cols = 0;
		cplx_matrices[i].rows = NULL;
		cplx_matrices[i].columns = NULL;
	}

	cplx_group_ranks = (SM_index_t *) malloc(cplx_size * sizeof(SM_index_t));
	num_generators = (SM_index_t *) malloc(cplx_size * sizeof(SM_index_t));
	if (cplx_group_ranks == NULL || num_generators == NULL)
		ERR_BAIL(mem_error)
}

#ifdef PRINT_DEBUG

/* this is BROKEN for the unified homology, really BROKEN */

/*
 * For testing only: print sums of squared incidence numbers for every
 * generator in the group.
 */
static void print_inums(SM_complex_t group)
{
	SM_index_t i;
	SM_value_t tmp;
	SparseVector *vec;
	SparseEntry *eptr;

	if (group < first_group || group > last_group) return;

	printf("  Group %d has %d generator(s)", group,
						num_generators[group]);
	if (group == first_group) { printf(".\n"); return; }
	else printf(": [");

 	vec = cplx_matrices[group - 1].rows;

	for (i = 0; i < cplx_group_ranks[group]; i++, vec++) {
		tmp = 0;
		for (eptr = vec -> entries; eptr != NULL; eptr = eptr->next)
			tmp += eptr->value * eptr->value;
		if (vec -> num_entries == -1) tmp = -1;
		printf("%d, ", tmp);
	}
	printf("]\n");
}

/*
 * For testing only: do print_inums for every non-empty group.
 */
static void print_all_inums(void)
{
	SM_complex_t group;

	printf("\n");
	for (group = first_group; group <= last_group; group++)
		print_inums(group);
}

#endif // #ifdef PRINT_DEBUG

/*
 * Given a matrix in PARI's sparse format, translate it into the internal one.
 * Four formats are supported:
 *   shrunk:   the entry is value * (row * 2^32 + column) with value = \pm1
 *             row and column are assumed to be not bigger than 2^31
 *   packed:   on a 64-bit architecture the same as shrunk
 *             on a 32-bit architecture every matrix entry occupies 2 places
 *             in the VECSMALL vector: ..., row, value * column, ...
 *             row and column are assumed to be not bigger than 2^31
 */
static void assign_matrix(SparseMatrix *matr, GEN entries_list, long list_len)
{
	GEN GEN_ptr = entries_list + 1;
	long i;
#ifdef LONG_IS_64BIT
	long tmp_val;
#endif
	SM_index_t row, column;
	SM_value_t value[2];
	int is_val_neg, is_odd_var;

	if (typ(entries_list) != t_VECSMALL) bailout();

	for (i = 0; i < list_len; i++, GEN_ptr++) {
#ifdef LONG_IS_64BIT
		tmp_val = (long) *GEN_ptr;
		is_val_neg = tmp_val < 0;
		tmp_val = labs(tmp_val);
		is_odd_var = (tmp_val & (1L << 31)) != 0;

		row = tmp_val >> 32;
		column = tmp_val & ((1L << 31) - 1);
#else
		row = (SM_index_t) *(GEN_ptr++);
		column = (SM_index_t) *GEN_ptr;
		is_val_neg = column < 0;
		is_odd_var = row < 0;

		column = abs(column);
		row = abs(row);
#endif
		value[0] = value[1] = 0;
		value[is_odd_var] = is_val_neg ? -1 : 1;

		if (add_m_entry(matr, row, column, value) == -1) bailout();
	}

	return;
}

static void init_diff_matrix(SM_complex_t matrix)
{
	/* only differentials between non-empty groups are interesting */
	if (matrix < first_group || matrix >= last_group) return;

#ifdef PRINT_DEBUG
	printf("Initializing matrix number %d.", matrix);
#endif

	if (init_s_matrix(cplx_matrices + matrix, cplx_group_ranks[matrix + 1],
				cplx_group_ranks[matrix]) == -1) bailout();
	assign_matrix(cplx_matrices + matrix, (GEN) pari_matrices[matrix + 1],
					itos((GEN) num_entries[matrix + 1]));

#ifdef PRINT_DEBUG
	printf("Matrix number %d has ", matrix);
	print_s_matrix(cplx_matrices + matrix);
#endif
}

/*
 * Assign values to the main internal variables.
 */
static int init_ranks(GEN c_ranks)
{
	SM_complex_t i;

	first_group = last_group = -1;

	for (i = 0; i < cplx_size; i++) {
		num_generators[i] = cplx_group_ranks[i] =
					(SM_index_t) itos((GEN) c_ranks[i + 1]);

		if (cplx_group_ranks[i] > 0) {
			if (first_group < 0) first_group = i;
			last_group = i;
		}
	}

#ifdef PRINT_DEBUG
	printf("\n first: %d  last: %d \n", first_group, last_group);
#endif

	/* the complex is empty */
	if (first_group < 0) return -1;

	return 0;
}

/*
 * Translate a matrix in the internal format into a PARI's matrix.
 */
static void matr2pari(SM_complex_t group, GEN pari_matr[2])
{
	SM_index_t i, j, row, col;
	SM_value_t val[2];
	SparseMatrix *matr = cplx_matrices + group;
	SparseVector *matr_cols = matr->columns, *matr_rows = matr->rows;
	SM_index_t n_rows = num_generators[group + 1];
	SM_index_t n_cols = num_generators[group];
	SM_index_t n_m_rows = matr->num_rows, n_m_cols = matr->num_cols;
	GEN pari_vec[2];
	char *matr_error = "matr2pari: matrix is corrupt";

	pari_matr[0] = cgetg(n_cols + 1, t_MAT);
	pari_matr[1] = cgetg(n_cols + 1, t_MAT);

	for (i = 1, col = 1; i <= n_cols; i++, col++) {
		while(matr_cols[col - 1].num_entries == -1 && col <= n_m_cols)
			col++;
		if (col > n_m_cols) ERR_BAIL(matr_error)

		pari_vec[0] = cgetg(n_rows + 1, t_COL);
		pari_vec[1] = cgetg(n_rows + 1, t_COL);
		for (j = 1, row = 1; j <= n_rows; j++, row++) {
			while(matr_rows[row - 1].num_entries == -1 &&
					row <= n_m_rows) row++;
			if (row > n_m_rows) ERR_BAIL(matr_error)

			/* remove_m_entry checks much more than get_m_entry */
			if (remove_m_entry(matr, row, col, val) == -1)
				bailout();

			pari_vec[0][j] = (long)stoi(val[0]);
			pari_vec[1][j] = (long)stoi(val[1]);
		}
		pari_matr[0][i] = (long) pari_vec[0];
		pari_matr[1][i] = (long) pari_vec[1];
	}
}

/*
 * Prepare the result to be sent back to PARI.
 */
static GEN feed2pari()
{
	SM_complex_t i, group;
	GEN main_vec = cgetg(4, t_VEC);
	GEN pari_matr[2];
	GEN numgen_vec = cgetg(cplx_size + 1, t_VEC);
	GEN matrices_vec[] = {cgetg(cplx_size, t_VEC), cgetg(cplx_size, t_VEC)};

	for (i = 1; i <= cplx_size; i++) numgen_vec[i] = (long)gen_0;
	main_vec[1] = (long) numgen_vec;

	for (i = 1; i < cplx_size; i++) matrices_vec[0][i] = (long)gen_0;
	main_vec[2] = (long) matrices_vec[0];

	for (i = 1; i < cplx_size; i++) matrices_vec[1][i] = (long)gen_0;
	main_vec[3] = (long) matrices_vec[1];

	if (first_group < 0) return main_vec;

	for (group = first_group; group <= last_group; group++) {
		/* nothing to do with an empty group */
		if (num_generators[group] == 0) continue;
		numgen_vec[group + 1] = (long)stoi(num_generators[group]);

		/* no matrices after the last group */
		if (group == last_group) continue;

		/* no matrices with zero size */
		if (num_generators[group + 1] == 0) continue;

		matr2pari(group, pari_matr);
		matrices_vec[0][group + 1] = (long)(pari_matr[0]);
		matrices_vec[1][group + 1] = (long)(pari_matr[1]);
	}

	/*
	for (group = first_group; group < last_group; group++)
		if (!is_matr_empty(cplx_matrices + group))
			ERR_BAIL("feed2pari: not all matrices are returned")
	*/

	return main_vec;
}

/*
 * Kill a generator gen_num in the group.
 */
static void kill_gen(SM_complex_t group, SM_index_t gen_num)
{
	if (group > first_group)
		if (erase_m_row(cplx_matrices + group - 1, gen_num, 1) == -1)
			bailout();

	if (group < last_group)
		if (erase_m_column(cplx_matrices + group, gen_num, 1) == -1)
			bailout();

	num_generators[group]--;
}

/*
 * Eliminate as many generators as possible in a given group.
 * Return 1 if some elimination was done and 0 otherwise.
 * If do_short is set, only consider generators with at most 2 incident ones.
 */
static int eliminate_gens(SM_complex_t group, int do_short)
{
	SM_index_t elim_cnt = 0, gen, inc_gen;
	SM_value_t gen_coeff[2], add_coeff[2];
	SparseEntry *cur_entry, *next_entry;
	int isfound = 0;

	SparseVector *inum_vectors;
	char *gen_error = "eliminate_gens: generator is not killed cleanly";

	/* check that the differential matrices are already initialized */
	if (group > first_group + 1 && cplx_matrices[group - 2].rows == 0)
		init_diff_matrix(group - 2);
	if (group > first_group && cplx_matrices[group - 1].rows == 0)
		init_diff_matrix(group - 1);
	if (group < last_group && cplx_matrices[group].rows == 0)
		init_diff_matrix(group);

	/* searching for invertible entries across rows is for some reason
	 * _much_ faster than down columns, especially when do_short is set */
	inum_vectors = cplx_matrices[group - 1].rows;

	for (gen = 1; gen <= cplx_group_ranks[group]; gen++, inum_vectors++) {
		if (inum_vectors->num_entries == -1) continue; // already gone
		if (do_short && (inum_vectors->num_entries > 2)) continue;

		inc_gen = find_v_unit(inum_vectors, gen_coeff);
		/* should be impossible after the previous check */
		// if (inc_gen == ERR_MVAL) bailout ();

		if (inc_gen == 0) continue; // no invertible incidence numbers

		elim_cnt++;
		isfound = 1;

		/* gen_coeff^2 == 1 and we need to use it for subtraction */
		gen_coeff[0] = -gen_coeff[0];
		gen_coeff[1] = -gen_coeff[1];

		/* entries in this column are being erased as the elimination
		 * is taking place, so a simple 'for' loop is not enough */
		cur_entry = inum_vectors->entries;
		while (cur_entry != NULL) {
			/* the current entry is going to be erased,
			 * so we need to remember where it is pointing to */
			next_entry = cur_entry->next;
			mult_Uvals(cur_entry->value, gen_coeff, add_coeff);
			if (cur_entry->index != inc_gen) {
				if (add_m_cols(cplx_matrices + group - 1,
					cur_entry->index, inc_gen, add_coeff) == -1)
						bailout();
			}
			cur_entry = next_entry;
		}

		/* a single entry should remain in this column by now ... */
		if (inum_vectors->num_entries != 1) ERR_BAIL(gen_error);
		kill_gen(group - 1, inc_gen);

		/* ... and now it has to be gone too */
		if (inum_vectors->num_entries != 0) ERR_BAIL(gen_error);
		kill_gen(group, gen);
	}

#ifdef PRINT_REDSTAT
	if (isfound) printf("%d | ", elim_cnt);
#endif

	return isfound;
}

/*
 * Reduce a chain complex with only free Abelian chain groups as far as
 * possible using a sequence of elementary collapses and merging of cell.
 *
 * Arguments are:
 *   size (length) of the complex probably without zero groups at the ends
 *   ranks of the chain groups
 *   matrices of chain differentials in the sparse format
 *   lengths of arrays representing the matrices
 *
 * The return value is a 2-component vector which contain
 *   ranks of the chain groups after the reduction
 *   matrices of chain differentials after the reduction
 *     (matrices of size 0 are substituted with 0 for better visualization)
 */
GEN reduce_s_complex_U(long c_size, GEN c_ranks, GEN d_matrices, GEN matr_lengths)
{
	GEN answer;
	SM_complex_t group;
	int cnt_short, cnt_full;

	malloc_arrays((SM_complex_t) c_size);
	pari_matrices = d_matrices;
	num_entries = matr_lengths;

	if (init_ranks(c_ranks) == -1) {
		/* complex is empty */
		answer = feed2pari();
		cleanup();
		return answer;
	}

#ifdef PRINT_REDSTAT
	printf("\n   ");
#endif
	for (group = first_group + 1; group <= last_group; group++) {
		/* the number of successful iterations for each group */
		cnt_short = cnt_full = 0;
#ifdef PRINT_REDSTAT
		printf("%d: ", group);
#endif
		/* eliminate as many generators in this group as possible,
		 * repeat the procedure if something was eliminated */
		while (eliminate_gens(group, 1)) cnt_short++;
		while (eliminate_gens(group, 0)) cnt_full++;
#ifdef PRINT_REDSTAT
		printf("%d+%d;  ", cnt_short, cnt_full);
#endif
	}

	answer = feed2pari();
	cleanup();
	return answer;
}
