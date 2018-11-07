/*
##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################
 */

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#ifndef _GUSbase
#define _GUSbase


SEXP pest_c(SEXP p, SEXP ep, SEXP v, SEXP ref, SEXP alt, SEXP nInd, SEXP nSnps);
SEXP pest_em_c(SEXP pinit, SEXP epinit, SEXP ref, SEXP alt, SEXP nInd, SEXP nSnps, SEXP nThreads, SEXP EMpara);
// SEXP gest_c(SEXP p, SEXP ep, SEXP v, SEXP ref, SEXP alt, SEXP nInd, SEXP nSnps);

#endif
