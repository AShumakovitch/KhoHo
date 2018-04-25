/*
 *    print_ranks.c --- 1. combine the ranks of chain groups and homology
 *                         together with ranks of chain differentials into
 *                         a single matrix in the TeX format.
 *                      2. combine the ranks of standard and reduced homology
 *                         (including optional torsion) into a single table
 *                         in the TeX format.
 *
 * Copyright (C) 2002, 2003, 2004 Alexander Shumakovitch <Shurik@Dartmouth.edu>
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
 *	install(print_ranks, vGGGGLLL, print_ranks, "./print_ranks.so")
 *	install(print_homology, vGGLGGLLL, print_homology, "./print_ranks.so")
 *	install(print_both_homology, vGGGGLL, print_both_homology,
 *							"./print_ranks.so")
 */

#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <pari/pari.h>

#ifdef min
#  undef min
#endif
#ifdef max
#  undef max
#endif
#define min(a,b)  ((a)<(b)?(a):(b))
#define max(a,b)  ((a)>(b)?(a):(b))
#define M_ENTR(mat, j, i)  ((int) (itos(gcoeff((mat), (j) + 1, (i) + 1))))

#if PARI_VERSION_CODE > PARI_VERSION(2,7,0)
#  define talker e_MISC
#endif

/*
 * Figure out whether there are non-zero torsion ranks in a torsion matrix
 * entry. Both old (the entry is an integer) and new (the entry is a vector)
 * formats are accepted.
 */
static int got_torsion(GEN tors_entry)
{
	int len, i;

	if (tors_entry == NULL) return 0;

	if (typ(tors_entry) == t_VEC || typ(tors_entry) == t_COL) {
		/* entry is a vector --- the new multiple torsion format;
		 * two formats are now supported:
		 *   list of ranks and 
		 *   list of pairs: order of the factors and their ranks */
		len = glength(tors_entry);
		if (len == 0) return 0;
		if (typ((GEN) tors_entry[1]) == t_VEC) return 1;
		for (i = 1; i <= len; i ++) {
			if (signe((GEN) tors_entry[i])) return 1;
		}
		return 0;
	} else {
		/* entry is an integer --- the old single torsion format */
		return (int) signe(tors_entry);
	}
}

/*
 * Print all torsion ranks separated by commas.  Both old (the entry is an
 * integer) and new (the entry is a vector) formats are accepted.
 */
static void print_torsion(FILE *fd, GEN tors_entry, char *preamble)
{
	int len, i, t_order;
	GEN tor1;

	if (tors_entry == NULL) return;

	if (typ(tors_entry) == t_VEC || typ(tors_entry) == t_COL) {
		/* entry is a vector --- the new multiple torsion format;
		 * two formats are now supported:
		 *   list of ranks (print in square brackets) and 
		 *   list of pairs: order of the factors and their ranks */
		len = glength(tors_entry);
		if (len == 0) return;
		if (typ((GEN) tors_entry[1]) == t_VEC) {
			fputs(preamble, fd);
			for (i = 1; i <= len; i ++) {
				tor1 = (GEN) tors_entry[i];
				t_order = (int) itos((GEN) tor1[1]);
				
				fprintf(fd, t_order > 2 ? 
					"\\torframe{$\\fam6{%d}_{%d}$}" :
								"{%d}_{%d}",
					(int) itos((GEN) tor1[2]), t_order);
			}
			fputs("$", fd);
		} else {
			fputs("\\lbr", fd);
			for (i = 1; i <= len; i ++) {
				fprintf(fd, "%d",
					(int) itos((GEN) tors_entry[i]));
				if (i < len) fprintf(fd, ",");
			}
			fputs("\\rbr", fd);
		}
	} else {
		/* entry is an integer --- the old single torsion format */
		fprintf(fd, "\\lbr%d\\rbr", (int) itos(tors_entry));
	}
}

/*
 * Print a repeated pattern of strings <l_half, r_half> (n_str times), but
 * finish with final instead of r_half.
 */
static
void print_h_line(FILE *fd, int n_str, char *l_half, char *r_half, char *final)
{
	int i;

	for (i = 0; i < n_str - 1; i++) {
		fputs(l_half, fd); fputs(r_half, fd); fputc('\n', fd);
	}
	fputs(l_half, fd); fputs(final, fd);
}

/*
 * Print the TeX preamble.
 */
static void print_TeX_common_preamble(FILE *fd)
{
	fputs("\
\\advance\\hsize 1.4truein \\advance\\hoffset -0.7truein\n\
\\advance\\vsize 1.0truein \\advance\\voffset -0.5truein\n\
\n\
\\def\\Cal#1{{\\fam2#1}}\n\
\n\
\\let\\TSp\\thinspace\n\
\\let\\NSp\\negthinspace\n\
\\def\\DSp{\\thinspace\\thinspace}\n\
\\def\\QSp{\\thinspace\\thinspace\\thinspace\\thinspace}\n\
\\def\\lbr{\\raise 0.5pt\\hbox{[}}\n\
\\def\\rbr{\\raise 0.5pt\\hbox{]}}\n\
\n\
\\nopagenumbers\n\
\\offinterlineskip\n\
\n\
\\def\\gobble#1{}\n\
\\def\\hline{\\noalign{\\hrule}}\n\
\\def\\torframe#1{\\vtop{\\vbox{\\hrule\\hbox{\\vrule\\strut #1\\vrule}}\\hrule}}\n\
", fd);
}

static void print_TeX_ranks_preamble(FILE *fd, int n_str)
{
	fputs("\\def\\dblhline{height 0.16667em\\gobble&&\n", fd);

	print_h_line(fd, n_str, "", "&\\vrule&", "\\cr\n");

	fputs("\
\\noalign{\\hrule}}\n\
\n\
\\def\\group#1#2#3#4{${\\raise 1pt\\vbox to 6.5pt{\\vss%\n\
% \\empty is repeated twice to account for the case of #2 being, well, empty.\n\
\\hbox{#4#1\\ifx\\empty#2\\empty\\else #2\\fi}}\\over\\lower 3pt\\hbox{#4#3}}$}\n\
\\def\\putdown#1{\\smash{\\vtop{\\null\\hbox{#1}}}}\n\
\n\
", fd);
}

static void print_TeX_homol_preamble(FILE *fd, int n_str)
{
	int i;

	fputs("\\def\\dblhline{\\hline height 0.16667em\\gobble&&\n", fd);

	for (i = 0; i < n_str; i++) fputc('&', fd);

	fputs("\\cr\n\
\\noalign{\\hrule}}\n\
\n\
\\def\\inline#1{\\leaders\\hrule\\hskip 0.33em plus 1fill%\n\
\\hbox{\\vbox to 0pt{\\vss\\hbox{\\TSp #1\\TSp}\\vss}}%\n\
\\leaders\\hrule\\hskip 0.33em plus 1fill\\vrule}\n\
\n\
", fd);
}

/*
 * Print a TeX macro for scaling a box to fit into the page.
 */
static void print_fitbox(FILE *fd)
{
	fputs("\
\\newdimen\\tabheight\n\
\\newcount\\sfactor\n\
\\newcount\\fittmp\n\
% scale the box #1 s.t. its width fits into #2; inspired by rotate.tex\n\
\\def\\fitbox#1#2{%\n\
% box's height shouldn't be too precise, 0.1pt error is fine\n\
% the sizes lie in the range of 2^26sp=2^10pt\\approx 36cm, and half\n\
% of this precision should suffice; \\sfactor is further multiplied by 16\n\
\\fittmp\\wd#1\\divide\\fittmp 8192%\n\
\\sfactor #2\\multiply\\sfactor 16\\divide\\sfactor\\fittmp%\n\
\\tabheight=\\ht#1\\advance\\tabheight\\dp#1%\n\
\\divide\\tabheight 8192\\multiply\\tabheight\\sfactor\\divide\\tabheight 16%\n\
\\hbox to#2{\\vbox to\\tabheight{%\n\
% no problems with precision of PostScript though\n\
\\rotstart{\\number #2\\space\\number\\wd #1\\space div dup scale}\\box#1\\rotfinish\n\
\\vss}\\hss}}\n\
", fd);
}

/*
 * Combine the ranks of chain groups and homology together with ranks
 * of chain differentials into a single table in the TeX format.
 *
 * Arguments are:
 *   an array with 4 elements: matrix size and gradings in the top left corner
 *   an array with 4 elements:
 *     1) matrix with chain group ranks
 *     2) matrix with chain differential ranks
 *     3) matrix with homology ranks (Betti numbers)
 *     4) matrix with homology torsion ranks (optional)
 *   name of the knot to be printed below the table
 *   name of the output file
 *   reduced flag: 1 ==> homology are assumed to be computed for the reduced
 *                       Khovanov chain complex (only the caption is changed).
 *   landscape flag: 1 ==> print the table rotated by 90 degrees
 *   fit-to-page flag: 1 ==> make the table as wide as the page;
 *                     2 ==> shrink the table to fit the page if it's too wide
 */
void print_ranks(GEN matr_sizes, GEN all_ranks, GEN knot_name, GEN outfile,
		long is_reduced, long is_landscape, long fit_table)
{
	int do_tors, tors2do, i_size, j_size, first_i, first_j;
	int left, right, top, bottom;
	int i, j, tmpi;
	GEN chain_ranks, D_ranks, Betti, tors_ranks, tors_entry;
	char *tmps;
	FILE *fd;

	if (typ(matr_sizes) != t_VEC || lg(matr_sizes) != 5)
		pari_err(talker, "the first argument is not a vector of length 4");
	if (typ(all_ranks) != t_VEC || lg(all_ranks) != 5)
		pari_err(talker, "the second argument is not a vector of length 4");

	chain_ranks = (GEN) all_ranks[1];
	D_ranks = (GEN) all_ranks[2];
	Betti = (GEN) all_ranks[3];
	tors_ranks = (GEN) all_ranks[4];

	if (typ(chain_ranks) != t_MAT)
		pari_err(talker, "the second argument's 1st entry is not a matrix");
	if (typ(D_ranks) != t_MAT)
		pari_err(talker, "the second argument's 2nd entry is not a matrix");
	if (typ(Betti) != t_MAT)
		pari_err(talker, "the second argument's 3rd entry is not a matrix");

	if (typ(knot_name) != t_STR)
		pari_err(talker, "the knot name is not a string");
	if (typ(outfile) != t_STR)
		pari_err(talker, "the output file name is not a string");

	if ((fd = fopen(GSTR(outfile), "w")) == NULL)
		pari_err(talker, "file cannot be opened");

	do_tors = typ(tors_ranks) == t_MAT;

	i_size = itos((GEN) matr_sizes[1]);
	j_size = itos((GEN) matr_sizes[2]);
	first_i = itos((GEN) matr_sizes[3]);
	first_j = itos((GEN) matr_sizes[4]);

	// find boundaries of a submatrix with non-zero chain group ranks
	top = left = INT_MAX;
	bottom = right = INT_MIN;
	for (i = 0; i < i_size; i++)
		for (j = 0; j < j_size; j++)
			if (M_ENTR(chain_ranks, j, i)) {
				top = min(top, j);
				bottom = max(bottom, j);
				left = min(left, i);
				right = max(right, i);
			}

	// matrix is apparently empty
	if (top == INT_MAX) {
		fputs("\\bye\n", fd);
		fclose(fd);
		return;
	}

	// some debugging
	// printf("l: %d, r: %d, t: %d, b: %d", left, right, top, bottom);

	print_TeX_common_preamble(fd);
	print_TeX_ranks_preamble(fd, right - left + 1);

	// load rotate.tex if rotating or scaling of the table is required
	if (is_landscape || fit_table)
		fputs("\\input rotate.tex\n\n", fd);

	fputs("\\newbox\\tablebox\n", fd);
	if (fit_table) print_fitbox(fd);

	// start of ialign and preamble
	fputs("\n% main table is put into \\tablebox\n", fd);
	fputs("\\setbox\\tablebox\\vbox{\\ialign{\%\n", fd);
	fputs("\\vrule\\TSp\\vrule #\\strut&", fd);
	fputs("\\DSp\\hfil #\\DSp\\vrule\\TSp\\vrule&\\kern.75em\n", fd);
	print_h_line(fd, right - left + 1,
			"\\TSp\\hfil #\\hfil\\TSp",
			"&\\hfil #\\hfil&",
			"\\kern.75em\\vrule\\TSp\\vrule\\cr\n\%\n");

	// top of the table
	fputs("\\hline\\dblhline\n", fd);

	// first line: primary grading
	fputs("height 11pt depth 4pt&&", fd);
	for (i = left; i <= right; i++) {
		fprintf(fd, "\n%d", first_i + i);
		if (i < right)
			fputs("&\\vrule&", fd);
		else
			fputs("\\cr\\hline\\dblhline\n\%\n", fd);
	}

	// main loop: ranks and Betti numbers
	for (j = top; j <= bottom; j++) {
		// vertical lines above the ranks of matrices
		fputs("height 0.2em depth 0.2em\\gobble&&\n", fd);
		print_h_line(fd, right - left + 1, "",
				"&\\vrule depth0pt&",
				"\\cr\n");

		// ranks of differential matrices.
		fputs("height 0pt\\gobble&&", fd);
		for (i = left; i < right; i++) {
			if (M_ENTR(chain_ranks, j, i) &&
					M_ENTR(chain_ranks, j, i + 1)) {
				tmpi = M_ENTR(D_ranks, j, i);
				if (tmpi < 10) tmps="\\QSp";
				else {
					if (tmpi < 100) tmps="\\DSp";
					else tmps="\\TSp";
				}
				fprintf(fd, "\n&%s\\putdown{%d}%s&",
						tmps, tmpi, tmps);
			} else
				fputs("\n&&", fd);
		}
		fputs("\\cr\n", fd);

		// secondary grading
		fprintf(fd, "&%d&", first_j - 2 * j);
		// Betti numbers [ torsion ranks ] / ranks of the chain groups
		for (i = left; i <= right; i++) {
			if (M_ENTR(chain_ranks, j, i)) {
				tors_entry = do_tors ? gcoeff(tors_ranks,
					       	j + 1, i + 1) : NULL;
				tors2do = got_torsion(tors_entry);
				tmps = (M_ENTR(Betti, j, i) || tors2do) ?
						"\\bf" : "";
				fprintf(fd, "\n\\group{%d}{",
					M_ENTR(Betti, j, i));
				if (tors2do)
					print_torsion(fd, tors_entry,
								", $\\fam6");
				fprintf(fd, "}{%d}{%s}",
					M_ENTR(chain_ranks, j, i), tmps);
			} else
				fputc('\n', fd);
			if (i < right)
				if (M_ENTR(chain_ranks, j, i) &&
						M_ENTR(chain_ranks, j, i + 1))
					fputs("&\\rightarrowfill&", fd);
				else
					fputs("&\\smash{\\vrule height 15pt}&",
							fd);
			else
				fputs("\\cr\n", fd);
		}

		// vertical lines below the ranks of matrices
		fputs("height 0.3em\\gobble&&\n", fd);
		print_h_line(fd, right - left + 1, "",
				"&\\smash{\\vrule height 10pt}&",
				"\\cr\n");

		fputs("\\hline\n\%\n", fd);
	}

	fputs("\\dblhline\n}}\n\n", fd);

	// print in the landscape mode, if necessary
	if (is_landscape) {
		fputs("\\setbox\\tablebox\\hbox{\\vbox to\\hsize{", fd);
		fputs("\\hsize=\\vsize{\%\n", fd);
	}
	fputs("\\null\\vfill\n$$", fd);

	switch (fit_table) {
		case 0:
			fputs("\\box\\tablebox", fd);
			break;
		case 1:
			fputs("\\fitbox\\tablebox\\hsize", fd);
			break;
		default:
			fputs("\\ifdim\\wd\\tablebox>\\hsize", fd);
			fputs("\\fitbox\\tablebox\\hsize", fd);
			fputs("\\else\\box\\tablebox\\fi", fd);
	}

	fputs("$$\n\n\\medskip\n\\centerline{", fd);
	fprintf(fd, "Ranks of the chain groups $%s\\Cal{C}^{i,j}$,\n",
			is_reduced ? "\\widetilde" : "");
	fprintf(fd, "differentials, and homology $%s\\Cal{H}^{i,j}$}\n",
			is_reduced ? "\\widetilde" : "");
	fprintf(fd, "\\medskip\n\\centerline{of the %sKhovanov chain complex\n",
			is_reduced ? "reduced " : "");
	fprintf(fd, "for %s}\n\\vfill\n\n", GSTR(knot_name));

	if (is_landscape)
		fputs("}}}\\rotl\\tablebox\n\n", fd);

	fputs("\\bye\n", fd);
	fclose(fd);
	return;
}

/*
 * Print the ranks of Khovanov homology (including optional torsion) as a
 * table in the TeX format.
 *
 * Arguments are:
 *   an array with 4 elements: matrix size and gradings in the top left corner
 *   an array with 2 elements:
 *     1) matrix with homology ranks (Betti numbers)
 *     2) matrix with homology torsion ranks (optional)
 *   Z_2 flag: 1 ==> print Z_2-homology ranks (the torsion must be present)
 *   name of the knot to be printed below the table
 *   name of the output file
 *   reduced flag: 1 ==> homology are assumed to be computed for the reduced
 *                       Khovanov chain complex (only the caption is changed).
 *   landscape flag: 1 ==> print the table rotated by 90 degrees
 */
void print_homology(GEN matr_sizes, GEN all_homology, long is_Z2_coeff,
		GEN knot_name, GEN outfile, long is_reduced, long is_landscape,
		long fit_table)
{
	int do_tors, tors2do, i_size, j_size, first_i, first_j;
	int left, right, top, bottom, i, j;
	GEN Betti, tors_ranks, tors_entry;
	char *tmps;
	FILE *fd;

	if (typ(matr_sizes) != t_VEC || lg(matr_sizes) != 5)
		pari_err(talker, "the first argument is not a vector of length 4");
	if (typ(all_homology) != t_VEC || lg(all_homology) != 3)
		pari_err(talker, "the second argument is not a vector of length 2");

	Betti = (GEN) all_homology[1];
	tors_ranks = (GEN) all_homology[2];

	if (typ(Betti) != t_MAT)
		pari_err(talker, "the second argument's 1st entry is not a matrix");

	if (typ(knot_name) != t_STR)
		pari_err(talker, "the knot name is not a string");
	if (typ(outfile) != t_STR)
		pari_err(talker, "the output file name is not a string");

	if ((fd = fopen(GSTR(outfile), "w")) == NULL)
		pari_err(talker, "file cannot be opened");

	do_tors = typ(tors_ranks) == t_MAT;
	is_Z2_coeff = is_Z2_coeff && do_tors;

	i_size = itos((GEN) matr_sizes[1]);
	j_size = itos((GEN) matr_sizes[2]);
	first_i = itos((GEN) matr_sizes[3]);
	first_j = itos((GEN) matr_sizes[4]);

	// find boundaries of a submatrix with non-zero homology
	top = left = INT_MAX;
	bottom = right = INT_MIN;
	for (i = 0; i < i_size; i++)
		for (j = 0; j < j_size; j++) {
			tors_entry = do_tors ?
				gcoeff(tors_ranks, j + 1, i + 1) : NULL;
			tors2do = got_torsion(tors_entry);

			if (M_ENTR(Betti, j, i) || tors2do)
			{
				top = min(top, j);
				bottom = max(bottom, j);
				left = min(left, i);
				right = max(right, i);
			}
		}

	// matrices are apparently empty
	if (top == INT_MAX) {
		fputs("\\bye\n", fd);
		fclose(fd);
		return;
	}

	// some debugging
	// printf("l: %d, r: %d, t: %d, b: %d", left, right, top, bottom);

	print_TeX_common_preamble(fd);
	print_TeX_homol_preamble(fd, right - left + 1);

	// load rotate.tex if rotating or scaling of the table is required
	if (is_landscape || fit_table)
		fputs("\\input rotate.tex\n\n", fd);

	fputs("\\newbox\\tablebox\n", fd);
	if (fit_table) print_fitbox(fd);

	// start of ialign and preamble
	fputs("\n% main table is put into \\tablebox\n", fd);
	fputs("\\setbox\\tablebox\\vbox{\\ialign{\%\n", fd);
	fputs("\\vrule\\TSp\\vrule #\\strut&", fd);
	fputs("\\DSp\\hfil #\\DSp\\vrule\\TSp\\vrule&\n", fd);
	print_h_line(fd, right - left + 1,
			"\\DSp\\hfil\\bf #\\hfil\\DSp\\vrule&",
			"",
			"#\\TSp\\vrule\\cr\n\%\n");

	// top of the table
	fputs("\\dblhline\n", fd);

	// first line: primary grading
	fputs("height 11pt depth 4pt&&", fd);
	for (i = left; i <= right; i++) {
		fprintf(fd, "\n\\rm\\DSp%d\\DSp&", first_i + i);
	}
	fputs("\\cr\\dblhline\n\%\n", fd);

	// main loop: Betti numbers and torsion
	for (j = top; j <= bottom; j++) {
		// secondary grading
		fprintf(fd, "height 13pt depth 5pt&%d&\n", first_j - 2 * j);
		// Betti numbers [ torsion ranks ] of the homology
		for (i = left; i <= right; i++) {
			tors_entry = do_tors ?
				gcoeff(tors_ranks, j + 1, i + 1) : NULL;
			tors2do = got_torsion(tors_entry);
			tmps = tors2do ? "" : "\\DSp";
			if (M_ENTR(Betti, j, i)) {
				fprintf(fd, "%s%d%s", tmps,
					M_ENTR(Betti, j, i), tmps);
				tmps = ", $\\fam6";
			} else tmps = "$\\fam6";
			if (tors2do) print_torsion(fd, tors_entry, tmps);
			fputs("&\n", fd);
		}
		fputs("\\cr", fd);
		if (j < bottom) fputs("\\hline", fd);
		fputs("\n%\n", fd);
	}

	fputs("\\dblhline\n}}\n\n", fd);

	// print in the landscape mode, if necessary
	if (is_landscape) {
		fputs("\\setbox\\tablebox\\hbox{\\vbox to\\hsize{", fd);
		fputs("\\hsize=\\vsize{\%\n", fd);
	}
	fputs("\\null\\vfill\n$$", fd);

	switch (fit_table) {
		case 0:
			fputs("\\box\\tablebox", fd);
			break;
		case 1:
			fputs("\\fitbox\\tablebox\\hsize", fd);
			break;
		default:
			fputs("\\ifdim\\wd\\tablebox>\\hsize", fd);
			fputs("\\fitbox\\tablebox\\hsize", fd);
			fputs("\\else\\box\\tablebox\\fi", fd);
	}

	fputs("$$\n\n\\medskip\n\\centerline{", fd);
	fprintf(fd, "Ranks and torsions of the homology $%s\\Cal{H}^{i,j}$}\n",
			is_reduced ? "\\widetilde" : "");
	fprintf(fd, "\\medskip\n\\centerline{of the %sKhovanov chain complex\n",
			is_reduced ? "reduced " : "");
	fprintf(fd, "for %s}\n\\vfill\n\n", GSTR(knot_name));

	if (is_landscape)
		fputs("}}}\\rotl\\tablebox\n\n", fd);

	fputs("\\bye\n", fd);
	fclose(fd);
	return;
}

/*
 * Combine the ranks of standard and reduced homology (including
 * optional torsion) into a single table in the TeX format.
 *
 * Arguments are:
 *   an array with 4 elements: matrix size and gradings in the top left corner
 *   an array with 4 elements:
 *     1) matrix with homology ranks (Betti numbers)
 *     2) matrix with homology torsion ranks (optional)
 *     3) matrix with reduced homology ranks
 *     4) matrix with reduced homology torsion ranks (optional)
 *   name of the knot to be printed below the table
 *   name of the output file
 *   landscape flag: 1 ==> print the table rotated by 90 degrees
 */
void print_both_homology(GEN matr_sizes, GEN all_homology,
		GEN knot_name, GEN outfile, long is_landscape, long fit_table)
{
	int do_tors, do_red_tors, tors2do, i_size, j_size, first_i, first_j;
	int left, right, top, bottom, i, j;
	GEN Betti, red_Betti, tors_ranks, red_tors_ranks, tors_entry;
	char *tmps;
	FILE *fd;

	if (typ(matr_sizes) != t_VEC || lg(matr_sizes) != 5)
		pari_err(talker, "the first argument is not a vector of length 4");
	if (typ(all_homology) != t_VEC || lg(all_homology) != 5)
		pari_err(talker, "the second argument is not a vector of length 4");

	Betti = (GEN) all_homology[1];
	tors_ranks = (GEN) all_homology[2];
	red_Betti = (GEN) all_homology[3];
	red_tors_ranks = (GEN) all_homology[4];

	if (typ(Betti) != t_MAT)
		pari_err(talker, "the second argument's 1st entry is not a matrix");
	if (typ(red_Betti) != t_MAT)
		pari_err(talker, "the second argument's 3rd entry is not a matrix");

	if (typ(knot_name) != t_STR)
		pari_err(talker, "the knot name is not a string");
	if (typ(outfile) != t_STR)
		pari_err(talker, "the output file name is not a string");

	if ((fd = fopen(GSTR(outfile), "w")) == NULL)
		pari_err(talker, "file cannot be opened");

	do_tors = typ(tors_ranks) == t_MAT;
	do_red_tors = typ(red_tors_ranks) == t_MAT;

	i_size = itos((GEN) matr_sizes[1]);
	j_size = itos((GEN) matr_sizes[2]);
	first_i = itos((GEN) matr_sizes[3]);
	first_j = itos((GEN) matr_sizes[4]);

	// find boundaries of a submatrix with non-zero homology
	top = left = INT_MAX;
	bottom = right = INT_MIN;
	for (i = 0; i < i_size; i++)
		for (j = 0; j < j_size; j++) {
			tors_entry = do_tors ?
				gcoeff(tors_ranks, j + 1, i + 1) : NULL;
			tors2do = got_torsion(tors_entry);
			tors_entry = do_red_tors ?
				gcoeff(red_tors_ranks, j + 1, i + 1) : NULL;
			tors2do = tors2do || got_torsion(tors_entry);

			if (M_ENTR(Betti, j, i) || M_ENTR(red_Betti, j, i) ||
					tors2do)
			{
				top = min(top, j);
				bottom = max(bottom, j);
				left = min(left, i);
				right = max(right, i);
			}
		}

	// matrices are apparently empty
	if (top == INT_MAX) {
		fputs("\\bye\n", fd);
		fclose(fd);
		return;
	}

	// some debugging
	// printf("l: %d, r: %d, t: %d, b: %d", left, right, top, bottom);

	print_TeX_common_preamble(fd);
	print_TeX_homol_preamble(fd, right - left + 1);

	// load rotate.tex if rotating or scaling of the table is required
	if (is_landscape || fit_table)
		fputs("\\input rotate.tex\n\n", fd);

	fputs("\\newbox\\tablebox\n", fd);
	if (fit_table) print_fitbox(fd);

	// start of ialign and preamble
	fputs("\n% main table is put into \\tablebox\n", fd);
	fputs("\\setbox\\tablebox\\vbox{\\ialign{\%\n", fd);
	fputs("\\vrule\\TSp\\vrule #\\strut&", fd);
	fputs("\\DSp\\hfil #\\DSp\\vrule\\TSp\\vrule&\n", fd);
	print_h_line(fd, right - left + 1,
			"\\DSp\\hfil\\bf #\\hfil\\DSp\\vrule&",
			"",
			"#\\TSp\\vrule\\cr\n\%\n");

	// top of the table
	fputs("\\dblhline\n", fd);

	// first line: primary grading
	fputs("height 11pt depth 4pt&&", fd);
	for (i = left; i <= right; i++) {
		fprintf(fd, "\n\\rm\\DSp%d\\DSp&", first_i + i);
	}
	fputs("\\cr\\dblhline\n\%\n", fd);

	// main loop: Betti numbers and torsion
	for (j = top; j <= bottom; j++) {
		// secondary grading
		fprintf(fd, "height 15pt depth 9pt&%d&\n", first_j - 2 * j);
		// Betti numbers [ torsion ranks ] of the standard homology
		for (i = left; i <= right; i++) {
			tors_entry = do_tors ?
				gcoeff(tors_ranks, j + 1, i + 1) : NULL;
			tors2do = got_torsion(tors_entry);
			tmps = tors2do ? "" : "\\DSp";
			if (M_ENTR(Betti, j, i)) {
				fprintf(fd, "%s%d%s", tmps,
					M_ENTR(Betti, j, i), tmps);
				tmps = ", $\\fam6";
			} else tmps = "$\\fam6";
			if (tors2do) print_torsion(fd, tors_entry, tmps);
			fputs("&\n", fd);
		}
		fputs("\\cr\n%\n", fd);

		// Betti numbers [ torsion ranks ] of the reduced homology
		// They must occupy one row less than the standard ones
		if (j == bottom) continue;

		fputs("\\omit\\span\\omit\\leaders\\hrule\\hfill&", fd);
		for (i = left; i <= right; i++) {
			tors_entry = do_red_tors ?
				gcoeff(red_tors_ranks, j + 2, i + 1) : NULL;
			tors2do = got_torsion(tors_entry);
			if (M_ENTR(red_Betti, j + 1, i) || tors2do) {
				fputs("\\omit\\inline{", fd);
				if (M_ENTR(red_Betti, j + 1, i)) {
					fprintf(fd, "%d",
						M_ENTR(red_Betti, j + 1, i));
					tmps = ", $";
				} else tmps = "$";
				if (tors2do)
					print_torsion(fd, tors_entry, tmps);
				fputs("}&\n", fd);
			} else {
				fputs("\\omit\\leaders\\hrule\\hfill&\n", fd);
			}
		}
		fputs("\\omit\\leaders\\hrule\\hfill\\vrule\\cr\n%\n", fd);
	}

	fputs("\\dblhline\n}}\n\n", fd);

	// print in the landscape mode, if necessary
	if (is_landscape) {
		fputs("\\setbox\\tablebox\\hbox{\\vbox to\\hsize{", fd);
		fputs("\\hsize=\\vsize{\%\n", fd);
	}
	fputs("\\null\\vfill\n$$", fd);

	switch (fit_table) {
		case 0:
			fputs("\\box\\tablebox", fd);
			break;
		case 1:
			fputs("\\fitbox\\tablebox\\hsize", fd);
			break;
		default:
			fputs("\\ifdim\\wd\\tablebox>\\hsize", fd);
			fputs("\\fitbox\\tablebox\\hsize", fd);
			fputs("\\else\\box\\tablebox\\fi", fd);
	}

	fputs("$$\n\n\\medskip\n\\centerline{", fd);
	fputs("Ranks of the standard and reduced homology\n", fd);
	fputs("$\\Cal{H}^{i,j}$ and $\\widetilde\\Cal{H}^{i,j}$", fd);
	fputs(" as well as their torsion}\n", fd);
	fputs("\\medskip\n\\centerline{of the Khovanov chain complex\n", fd);
	fprintf(fd, "for %s}\n\\vfill\n\n", GSTR(knot_name));

	if (is_landscape)
		fputs("}}}\\rotl\\tablebox\n\n", fd);

	fputs("\\bye\n", fd);
	fclose(fd);
	return;
}
