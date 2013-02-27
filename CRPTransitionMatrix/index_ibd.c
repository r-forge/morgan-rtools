/* MORGAN: Copyright 1991-2005: Univ. of Washington. All rights reserved. */
/* You may copy, use, modify, but NOT redistribute, modified or otherwise.*/
/* See: http://www.stat.washington.edu/thompson/Genepi/license.shtml      */

/********************************************************************
 *   A set of routines for scoring ibd, and related stuff        *
 *       Updated  EAT, March 2000                                *
 *******************************************************************/

#include "markers.h"

static int *cst = NULL;

/*  routine to combine l-sampler and m-sampler scoring calls
   EAT  March 2000           */

/*  version for multiple gamete sets;            EAT  Mar 2001           */
/*  version for malloced arrays                  MAJ  Aug 2001           */

int score_genes (int jj, int np, int *gaga, int *inda, int *list,
                 int *proba, int *gama, double *nscore, int nps, int *list2)
{
   int   k, p, bb;

   bb = 0;
   get_genes (np, jj, gaga, proba, gama);
   p = score_ibd (np, gaga);

   for (k = 0; k < nps; ++k)
      if (p == list2[k])
      {
         bb = 1;
         break;
      }

   ++nscore[list[inda[p]]];             /* list now indexed by n gametes */
   return (bb);
}

/* version for scoring exact-probabilities  */
void score_xact (int np, int *gaga, int *inda, int *list,
                 int *proba, int *gama, double *xpscore, double xp)
{
   int   k, p;

   get_genes (np, 0, gaga, proba, gama);
   p = score_ibd (np, gaga);
   k = list[inda[p]];
   xpscore[k] += xp;                    /* list now indexed by n gametes */
   return;
}

void score_windows (double **wscore, int winsz, int *wloc, int nloc)
{
   int   i, pt, bfac;

   bfac = 1;
   for (i = 0; i < winsz - 1; ++i)
      bfac *= 2;
   pt = wloc[0];
   for (i = 1; i < winsz; ++i)
      pt = pt * 2 + wloc[i];

   for (i = 0; i <= (nloc - winsz + 1); ++i)
   {
      ++wscore[i][pt];
      pt = pt - (pt / bfac) * bfac;
      pt = 2 * pt + wloc[i + winsz];
   }
   return;
}

/******************************************************************

/* program to index the possible ibd states for nn genes */

void invert_label (int clab, int nn, int *genes)
{
   int   i, j, cc, gg, flag;

   flag = 0;
   i = nn - 1;
   cc = clab;
   while (i > 0)
   {
      if (cc == 0)
         break;
      gg = cc - ((cc / i) * i);
/*    if ((gg != 0) || (cc < max[i]))           */
      if ((gg != 0) || (cc < CDgg[i + 1]))
      {
         genes[i] = gg + 1;
         cc = cc / i;
      }
      else
      {
         flag = 1;
         break;
      }
      i--;
   }
   if (flag == 0)
      for (j = 0; j <= i; ++j)
         genes[j] = 1;
   if (flag == 1)
      for (j = 0; j <= i; ++j)
         genes[j] = j + 1;
   return;
}
/*********************************************************************


/* next state in cycle  */

int next_state (int nn, int *ng, int *state)
{
   int   i, j;

   i = nn - 1;
   while ((ng[i - 1] < ng[i]) && (i > 0))
      --i;
   if (i == 0)
      return (0);

   ++state[i];
   if (ng[i] < state[i])
      ng[i] = state[i];
   for (j = i + 1; j < nn; ++j)
      state[j] = 1;
   for (j = i + 1; j < nn; ++j)
      ng[j] = ng[i];

   return (1);
}

int score_either (int g1, int g2, int g3)
{
   if ((g2 == g1) | (g3 == g1))
      return (1);
   else
      return (0);
}

int score_whichever (int g1, int g2, int g3, int g4)
{
   if ((g1 == g3) || (g1 == g4) || (g2 == g3) || (g2 == g4))
      return (1);
   else
      return (0);
}

int list_states (int jj, int nn, int *state_list)
{
   int   i;
   for (i = 1; i <= nn; ++i)
      if ((jj + 1) == state_list[i])
         return (i);
   return (0);
}

int sibs_states (int nn, int *states)
{
   int   i, s1, s2;

   for (i = 0; i < nn; i = i + 2)
      if (states[i] == states[i + 1])
         return (0);

   s1 = score_whichever (states[0], states[1], states[2], states[3]);
   s2 = score_whichever (states[0], states[1], states[4], states[5]);
   if (s1 * s2 == 1)
      return (1);
   else if (s1 == 1)
      return (2);
   else if (s2 == 1)
      return (3);
   else
      return (4);
}

int reduce_states_jv (int *genes, int *index)
{
   int   nrd;
   int   gaga[4];
   int   ncd;
   nrd = 4;
   gaga[0] = genes[0];
   gaga[1] = genes[1];
   if (gaga[0] == gaga[1])
      ncd = 1;
   else
      ncd = 2;
   if (score_either (genes[0], genes[2], genes[3]))
      gaga[2] = gaga[0];
   else if (score_either (genes[1], genes[2], genes[3]))
      gaga[2] = gaga[1];
   else
   {
      ++ncd;
      gaga[2] = ncd;
   }

   if (score_either (genes[0], genes[4], genes[5]))
      gaga[3] = gaga[0];
   else if (score_either (genes[1], genes[4], genes[5]))
      gaga[3] = gaga[1];
   else if (score_whichever (genes[2], genes[3], genes[4], genes[5]))
      gaga[3] = gaga[2];
   else
   {
      ++ncd;
      gaga[3] = ncd;
   }

   return (index[score_ibd (nrd, gaga)]);
}

/* routine to read a fixed set of states to be scored */

int read_states (FILE * fp, int *state_list)
{
   int   i, nn;
   fscanf (fp, "%d", &nn);
   for (i = 0; i < nn; ++i)
      fscanf (fp, "%d", state_list + i);
   return (nn);
}

int index_ibd (int nn, int rd, int nrd, int *index, int *label, int *list,
               FILE * fp)
/* rd is flag, for reducing list to nrd genes */
{
   /* compute a label and set indices for each state */
   /* (for nn = 6, there are 203 states) */

   int   i, jj, clab;
/* int   states[8], ng[8];   */
   int  *states, *ng;
   int   rd_label[30], rd_index[40];    /* max nrd=4, currently */
   int   rdpt;

   states = malloc (nn * sizeof (int));
   alloc_chk (WHERE, states);
   ng = malloc (nn * sizeof (int));
   alloc_chk (WHERE, ng);

   rdpt = 0;
   if (rd == 1)
      rdpt = index_ibd (nrd, 0, 0, rd_index, rd_label, rd_index, fp);
   if (rd == -2)
      rdpt = read_states (fp, rd_index);
   if (rd == 2)
      rdpt = 4;

   jj = -1;
   for (i = 0; i < nn; ++i)
   {
      ng[i] = 1;
      states[i] = 1;
   }
   ng[nn - 1] = 0;
   states[nn - 1] = 0;

   while (next_state (nn, ng, states) > 0)
   {
      clab = 0;
      ++jj;
      for (i = 0; i < nn; ++i)          /* labels state number jj+1 */
         clab = clab * i + states[i] - 1;

      index[clab] = jj;
      label[jj] = clab;
      if (rd == 1)
         list[jj] = reduce_states_jv (states, rd_index);
      else if (rd == 0)
         list[jj] = jj + 1;
      else if (rd == -2)
         list[jj] = list_states (jj, rdpt, rd_index);
      else if (rd == 2)
         list[jj] = sibs_states (nn, states);
/*
   printf ("\n chk state %d %d %d   ", jj, clab, list[jj] );
   for (i=0; i<nn; ++i) printf (" %d", states[i]);
   invert_label(clab,nn,states);  printf(" chk invert ");
   for (i=0; i<nn; ++i) printf (" %d", states[i]);
 */
   }

   free (states);
   free (ng);

/*  printf ("\n number of patterns = %d", jj+1);                    */

   if (rdpt == 0)
      return (jj + 1);
   else
      return (rdpt);
}

/* program to score ibd states for nn genes */

int score_ibd (int nn, int *genes)
{
   int   i, j, ng, clab;                /* j = current gene(s) to consider
                                         * ng = distinct genes so far      */
   if (!cst)
   {
      cst = malloc (MxGam * sizeof (int));
      alloc_chk (WHERE, cst);
   }
   cst[0] = 1;
   clab = 0;
   ng = 1;
   for (j = 1; j < nn; ++j)
   {
      cst[j] = 0;
      for (i = 0; i < j; ++i)
      {
         if (cst[j] == 0)
            if (genes[j] == genes[i])
               cst[j] = cst[i];
      }
      if (cst[j] == 0)
      {
         ++ng;
         cst[j] = ng;
      }
      clab = j * clab + cst[j] - 1;
   }
   return (clab);
}

void get_genes (int nn, int j, int *genes, int *proband, int *gamete)
{
   int   i;
   for (i = 0; i < nn; ++i)
   {
      if (gamete[i] == 0)
         genes[i] = gg0[j][proband[i]];
      else
         genes[i] = gg1[j][proband[i]];
   }
   return;
}
void free_cst (void)
{
   if (cst)
   {
      free (cst);
      cst = NULL;
   }
   return;
}
