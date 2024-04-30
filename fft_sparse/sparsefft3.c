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

#include "sparsefft3.h"

static void do_fft(fftw_plan p,
		   int howmany, fftw_complex *data, int stride, int dist,
		   fftw_complex *scratch)
{
	fftw(p, howmany, data, stride, dist, scratch, 1, 0);
}

/***************************************************************************/
/* The following are various routines to transform one or two dimensions
   of a 3d array at a time, taking advantage of the sparsity pattern that
   was passed to sparsefft3_create_plan.   For some of them (the ones
   that do two FFTs at a time) we also supply an inverse routine. */

/* Given the 3d array data conforming to the sparsity pattern,
   transform first_dim and then second_dim, arranged in planes
   along plane_dim. */
static void fft2_planes(sparsefft3_plan p,
			int plane_dim, int first_dim, int second_dim,
			fftw_complex *data)
{
     int iplane;
     for (iplane = 0; iplane < p->n[plane_dim]; ++iplane) {
	  fftw_complex *dataplane = data + iplane * p->nafter[plane_dim];
	  sparsefft3_range r;

	  r = p->range1[plane_dim][second_dim][iplane];
	  while (r) {
	       do_fft(p->p[first_dim], r->max - r->min + 1,
		      dataplane + r->min * p->nafter[second_dim],
		      p->nafter[first_dim], p->nafter[second_dim],
		      p->scratch);
	       r = r->next;
	  }
	  
	  if (p->range1[plane_dim][second_dim][iplane])
	       do_fft(p->p[second_dim], p->n[first_dim],
		      dataplane,
		      p->nafter[second_dim], p->nafter[first_dim],
		      p->scratch);
     }
}

/* inverse of fft2_planes. */
static void ifft2_planes(sparsefft3_plan p,
			int plane_dim, int first_dim, int second_dim,
			fftw_complex *data)
{
     int iplane;
     for (iplane = 0; iplane < p->n[plane_dim]; ++iplane) {
	  fftw_complex *dataplane = data + iplane * p->nafter[plane_dim];
	  sparsefft3_range r;

	  r = p->range1[plane_dim][second_dim][iplane];

	  if (r)
	       do_fft(p->p[second_dim], p->n[first_dim],
		      dataplane,
		      p->nafter[second_dim], p->nafter[first_dim],
		      p->scratch);

	  r = p->range1[plane_dim][second_dim][iplane];
	  while (r) {
	       do_fft(p->p[first_dim], r->max - r->min + 1,
		      dataplane + r->min * p->nafter[second_dim],
		      p->nafter[first_dim], p->nafter[second_dim],
		      p->scratch);
	       r = r->next;
	  }
     }
}

/* As fft2_planes, except that we only transform first_dim and leave
   second_dim untransformed. */
static void fft1_planes(sparsefft3_plan p,
			int plane_dim, int first_dim, int second_dim,
			fftw_complex *data)
{
     int iplane;
     for (iplane = 0; iplane < p->n[plane_dim]; ++iplane) {
	  fftw_complex *dataplane = data + iplane * p->nafter[plane_dim];
	  sparsefft3_range r;

	  r = p->range1[plane_dim][second_dim][iplane];
	  while (r) {
	       do_fft(p->p[first_dim], r->max - r->min + 1,
		      dataplane + r->min * p->nafter[second_dim],
		      p->nafter[first_dim], p->nafter[second_dim],
		      p->scratch);
	       r = r->next;
	  }
     }
}

/* Given the 3d array data conforming to the sparsity pattern, but
   after the first dim (plane_dim) has been transformed, transform
   second_dim then third_dim, arranged in planes along plane_dim. */
static void fft2_planes2(sparsefft3_plan p,
			 int plane_dim, int second_dim, int third_dim,
			 fftw_complex *data)
{
     int iplane;
     for (iplane = 0; iplane < p->n[plane_dim]; ++iplane) {
          fftw_complex *dataplane = data + iplane * p->nafter[plane_dim];
          sparsefft3_range r;
	  
	  r = p->range2[third_dim];
	  while (r) {
               do_fft(p->p[second_dim], r->max - r->min + 1,
                      dataplane + r->min * p->nafter[third_dim],
                      p->nafter[second_dim], p->nafter[third_dim],
                      p->scratch);
               r = r->next;
          }

	  do_fft(p->p[third_dim], p->n[second_dim],
		 dataplane,
		 p->nafter[third_dim], p->nafter[second_dim],
		 p->scratch);
     }
}

/* inverse of fft2_planes2 */
static void ifft2_planes2(sparsefft3_plan p,
			  int plane_dim, int second_dim, int third_dim,
			  fftw_complex *data)
{
     int iplane;
     for (iplane = 0; iplane < p->n[plane_dim]; ++iplane) {
          fftw_complex *dataplane = data + iplane * p->nafter[plane_dim];
          sparsefft3_range r;
	  
	  do_fft(p->p[third_dim], p->n[second_dim],
		 dataplane,
		 p->nafter[third_dim], p->nafter[second_dim],
		 p->scratch);

	  r = p->range2[third_dim];
	  while (r) {
               do_fft(p->p[second_dim], r->max - r->min + 1,
                      dataplane + r->min * p->nafter[third_dim],
                      p->nafter[second_dim], p->nafter[third_dim],
                      p->scratch);
               r = r->next;
          }
     }
}

#define TEST_FFT1_LOOP 0

/* Given the 3d array data, completely transform the dimension which_dim,
   assuming no sparsity (i.e. all sparsity has been destroyed already). */
static void fft1_all(sparsefft3_plan p, int which_dim, fftw_complex *data)
{
#if ! (defined(USE_CRAYFFT) || TEST_FFT1_LOOP)
     /* If we are doing the first or the last dimension, we can
	process the whole thing in one do_fft call.  We can't do
	this on the Cray, though, since on that machine we have
	only allocated enough scratch space to do two dimensions at a time. */
     if (which_dim == 0)
	  do_fft(p->p[0], p->nafter[0], data, p->nafter[0], 1, p->scratch);
     else if (which_dim == 2)
	  do_fft(p->p[2], p->n[0] * p->n[1], data, 1, p->n[2], p->scratch);
     else
#endif
     {
	  int first_dim, second_dim, ifirst, isecond;
	  if ((which_dim + 1) % 3 < (which_dim + 2) % 3) {
	       first_dim = (which_dim + 1) % 3;
	       second_dim = (which_dim + 2) % 3;
	  }
	  else {
	       first_dim = (which_dim + 2) % 3;
	       second_dim = (which_dim + 1) % 3;
	  }
	  for (ifirst = 0; ifirst < p->n[first_dim]; ++ifirst)
	       do_fft(p->p[which_dim], p->n[second_dim],
		      data + ifirst * p->nafter[first_dim],
		      p->nafter[which_dim], p->nafter[second_dim], p->scratch);
     }
}

/***************************************************************************/

void sparsefft3(sparsefft3_plan p, fftw_complex *data)
{
     if (p->sparsedir == SPARSEFFT3_SPARSEINPUT) {
	  if (p->fft2_first) {
	       fft2_planes(p, p->third_dim,p->first_dim,p->second_dim, data);
	       fft1_all(p, p->third_dim, data);
	  }
	  else {
	       fft1_planes(p, p->third_dim,p->first_dim,p->second_dim, data);
	       fft2_planes2(p, p->first_dim,p->second_dim,p->third_dim, data);
	  }
     }
     else {
	  if (p->fft2_first) {
	       fft1_all(p, p->third_dim, data);
	       ifft2_planes(p, p->third_dim,p->first_dim,p->second_dim, data);
	  }
	  else {
	       ifft2_planes2(p, p->first_dim,p->second_dim,p->third_dim, data);
	       fft1_planes(p, p->third_dim,p->first_dim,p->second_dim, data);
	  }
     }
}
