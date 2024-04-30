/* Copyright (C) 2000 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef SPARSEFFT3_H
#define SPARSEFFT3_H

#ifdef __cplusplus
extern "C" {
#endif                          /* __cplusplus */

#include <fftw.h>

struct sparsefft3_range_struct_ {
     int min, max;
     struct sparsefft3_range_struct_ *next;
};
typedef struct sparsefft3_range_struct_ sparsefft3_range_struct;
typedef sparsefft3_range_struct *sparsefft3_range;

typedef enum {
     SPARSEFFT3_SPARSEINPUT,
     SPARSEFFT3_SPARSEOUTPUT
} sparsefft3_sparsedir;

typedef struct {
     sparsefft3_sparsedir sparsedir;
     int n[3], nafter[3];
     fftw_plan p[3];
     sparsefft3_range *range1[3][3], range2[3];
     fftw_complex *scratch;
     int fft2_first, first_dim, second_dim, third_dim;
} sparsefft3_plan_struct;
typedef sparsefft3_plan_struct *sparsefft3_plan;

typedef int (*sparsefft3_nonzero_func) (int x[3], void *data);

extern sparsefft3_plan sparsefft3_create_plan(int nx, int ny, int nz,
					      int dir, int flags,
					      sparsefft3_sparsedir sparsedir,
					      sparsefft3_nonzero_func nonzero,
					      void *nonzero_data);

extern sparsefft3_plan sparsefft3_create_plan_specific(int nx, int ny, int nz,
					      int dir, 
                               int flags,
                               sparsefft3_sparsedir sparsedir,
					      sparsefft3_nonzero_func nonzero,
					      void *nonzero_data,
					      fftw_complex *data);

extern void sparsefft3_destroy_plan(sparsefft3_plan p);

extern void sparsefft3(sparsefft3_plan p, fftw_complex *data);

#ifdef __cplusplus
}                               /* extern "C" */
#endif                          /* __cplusplus */

#endif /* SPARSEFFT3_H */
