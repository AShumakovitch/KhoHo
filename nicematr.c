/*
 *    nicematr.c --- simple program to print a PARI matrix in a nice way
 *                   (with aligned columns).
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
 * 	install(nicematr, vGL, nicematr, "./nicematr.so")
 */

#include <string.h>
#include <stdio.h>
#include <pari/pari.h>

#if PARI_VERSION_CODE > PARI_VERSION(2,7,0)
#  define talker e_MISC
#endif

/*
 * Compute maximal length of entries in every column of matr. The result
 * is put into an array of the appropriate length pointed to by maxlen.
 */
static void maxlength(GEN matr, int *maxlen)
{
	long i, j, max_i, max_j;
	int curlength;
	char *buf;
	GEN entry;

	max_j = lg(matr);
	max_i = lg((GEN) matr[1]);

	for (j = 1; j < max_j; j++) maxlen[j] = -1;

	for (i = 1; i < max_i; i++) {
		for (j = 1; j < max_j; j++) {
			entry = gcoeff(matr, i, j);

			buf = GENtostr(entry);
			curlength = strlen(buf);
			if (typ(entry) == t_STR) curlength += 2;
			if (curlength > maxlen[j]) maxlen[j] = curlength;

			free(buf);
		}
	}
}

/*
 * Print a PARI matrix in a nice way (with aligned columns).
 * It's done by computing the maximal length of entries in every column first.
 * If replace_empty is true, replace all empty strings with '.'
 */
void nicematr(GEN matr, long replace_empty)
{
	long i, j, max_i, max_j;
	int maxlen[lg(matr)], sp, matsize;
	char *buf;
	GEN entry;

	if (typ(matr) != t_MAT) pari_err(talker, "the argument is not a matrix");

	max_j = lg(matr);
	if (max_j == 1 || lg((GEN) matr[1]) == 1) {
		pariputs("[;]\n");
		return;
	}
	max_i = lg((GEN) matr[1]);

	maxlength(matr, maxlen);
	for (j = 1, matsize = 0; j < max_j; j++) matsize += maxlen[j];
	matsize += 2 * (max_j - 2);

	pariputc('\n');
	for (i = 1; i < max_i; i++) {
		pariputc('[');
		for (j = 1; j < max_j; j++) {
			entry = gcoeff(matr, i, j);
			buf = GENtostr(entry);

			sp = maxlen[j] - strlen(buf);
			if (typ(entry) == t_STR) sp -= 2;
			for (; sp > 0; sp--) pariputc(' ');

			if (replace_empty && !strlen(buf))
				pariputs(" .");
			else {
				if (typ(entry) == t_STR) pariputc('"');
				pariputs(buf);
				if (typ(entry) == t_STR) pariputc('"');
			}

			free(buf);
			if (j < max_j - 1) pariputs("  ");
		}
		pariputs("]\n");
		if (i < max_i - 1) {
			pariputc('[');
			for (j = 0; j < matsize; j++) pariputs(" ");
			pariputs("]\n");
		}
	}
}
