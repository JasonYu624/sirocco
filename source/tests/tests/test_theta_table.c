/** ********************************************************************************************************************
 *
 *  @file test_theta_table.c
 *  @date November 2026
 *
 *  @brief Unit tests for the theta_table central source (agn.c, Central_object.geometry_for_source = theta_table)
 *
 *  A three-bin, two-group table with P1 anisotropies +0.3 / 0 / -0.3 is written to a temporary file and read back.
 *  Checked: the total; band luminosities against hand values; and, over 2e5 sampled packets, the polar-bin shares,
 *  positions on the sphere inside the chosen bin, frequencies inside the bin's groups and the band, unit directions in
 *  the outward hemisphere, and the mean cosine to the normal against the analytic (1/3 + 3f/4) / (1/2 + f) of the
 *  emergent P1 law mu (1 + 3 f mu).  Statistical checks use 5 sigma with a fixed seed.
 *
 * ****************************************************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <CUnit/CUnit.h>

#include "../../atomic.h"
#include "../../sirocco.h"

#define NSAMPLE 200000
#define RSPHERE 1.0e13

static char table_path[] = "/tmp/sirocco_theta_table_testXXXXXX";
static const double FK[3] = { 0.3, 0.0, -0.3 };

/* L[k][g] in units of 1e40 erg/s: bin 0 both groups, bin 1 group 0 only, bin 2 group 1 only */
static const double LKG[3][2] = { {1.0, 1.0}, {2.0, 0.0}, {0.0, 6.0} };

static int
suite_init (void)
{
  int fd = mkstemp (table_path);
  FILE *fp;
  if (fd < 0)
    return -1;
  fp = fdopen (fd, "w");
  fprintf (fp, "# test table\n3 2\n0 30 60 90\n1e15 2e15 4e15\n%g %g %g\n", FK[0], FK[1], FK[2]);
  for (int k = 0; k < 3; k++)
    fprintf (fp, "%.1fe40 %.1fe40\n", LKG[k][0], LKG[k][1]);
  fclose (fp);
  init_rand (20261114);
  return 0;
}

static int
suite_teardown (void)
{
  unlink (table_path);
  return 0;
}

static void
test_read_total (void)
{
  double total = theta_table_read (table_path);
  CU_ASSERT_DOUBLE_EQUAL (total / 1.0e41, 1.0, 1e-12);
}

static void
test_band_luminosity (void)
{
  CU_ASSERT_DOUBLE_EQUAL (theta_table_band_lum (1e15, 4e15) / 1e41, 1.0, 1e-12);
  /* half of group 0: (1 + 2 + 0) / 2 */
  CU_ASSERT_DOUBLE_EQUAL (theta_table_band_lum (1e15, 1.5e15) / 1e40, 1.5, 1e-12);
  /* half of group 1: (1 + 0 + 6) / 2 */
  CU_ASSERT_DOUBLE_EQUAL (theta_table_band_lum (3e15, 4e15) / 1e40, 3.5, 1e-12);
  /* outside the table */
  CU_ASSERT_DOUBLE_EQUAL (theta_table_band_lum (5e15, 6e15), 0.0, 0.0);
}

/* sample n packets in [f1, f2]; check every packet; return per-bin counts and mean mu' */
static void
sample_and_check (double f1, double f2, int n, int count[3], double mu_sum[3], int *npos)
{
  struct photon p;
  double th_edge[4] = { 0.0, PI / 6.0, PI / 3.0, PI / 2.0 };
  double nu_edge[3] = { 1e15, 2e15, 4e15 };
  int bad_r = 0, bad_theta = 0, bad_freq = 0, bad_dir = 0;

  for (int k = 0; k < 3; k++)
  {
    count[k] = 0;
    mu_sum[k] = 0.0;
  }
  *npos = 0;
  theta_table_cdf (f1, f2);
  for (int i = 0; i < n; i++)
  {
    int k = theta_table_sample (&p, RSPHERE, f1, f2);
    double r = length (p.x);
    double ct = fabs (p.x[2]) / r;
    double mu = (p.lmn[0] * p.x[0] + p.lmn[1] * p.x[1] + p.lmn[2] * p.x[2]) / r;
    int g_ok = 0;

    if (k < 0 || k > 2)
    {
      CU_FAIL ("bin index out of range");
      return;
    }
    count[k]++;
    mu_sum[k] += mu;
    if (p.x[2] > 0)
      (*npos)++;
    if (fabs (r / RSPHERE - 1.0) > 1e-12)
      bad_r++;
    if (ct > cos (th_edge[k]) + 1e-12 || ct < cos (th_edge[k + 1]) - 1e-12)
      bad_theta++;
    for (int g = 0; g < 2; g++)
      if (LKG[k][g] > 0 && p.freq >= fmax (nu_edge[g], f1) && p.freq <= fmin (nu_edge[g + 1], f2))
        g_ok = 1;
    if (!g_ok)
      bad_freq++;
    if (fabs (length (p.lmn) - 1.0) > 1e-12 || mu < -1e-12 || mu > 1.0 + 1e-12)
      bad_dir++;
  }
  CU_ASSERT_EQUAL (bad_r, 0);
  CU_ASSERT_EQUAL (bad_theta, 0);
  CU_ASSERT_EQUAL (bad_freq, 0);
  CU_ASSERT_EQUAL (bad_dir, 0);
}

static void
test_sample_full_band (void)
{
  int count[3], npos;
  double mu_sum[3];
  double share[3] = { 0.2, 0.2, 0.6 };

  sample_and_check (1e15, 4e15, NSAMPLE, count, mu_sum, &npos);
  for (int k = 0; k < 3; k++)
  {
    double expect = NSAMPLE * share[k];
    double sigma = sqrt (NSAMPLE * share[k] * (1.0 - share[k]));
    double mu_mean = mu_sum[k] / count[k];
    double mu_expect = (1.0 / 3.0 + 0.75 * FK[k]) / (0.5 + FK[k]);
    /* the variance of mu' is below 1/12 + small; 0.3 / sqrt(n) is a safe sigma */
    double mu_sigma = 0.3 / sqrt ((double) count[k]);
    CU_ASSERT (fabs (count[k] - expect) < 5.0 * sigma);
    CU_ASSERT (fabs (mu_mean - mu_expect) < 5.0 * mu_sigma);
  }
  CU_ASSERT (fabs (npos - 0.5 * NSAMPLE) < 5.0 * sqrt (0.25 * NSAMPLE));
}

static void
test_sample_one_group_band (void)
{
  int count[3], npos;
  double mu_sum[3];

  /* group 0 only: bins 0 and 1 carry 1 and 2 (x 1e40), bin 2 none */
  sample_and_check (1e15, 2e15, NSAMPLE, count, mu_sum, &npos);
  CU_ASSERT_EQUAL (count[2], 0);
  CU_ASSERT (fabs (count[0] - NSAMPLE / 3.0) < 5.0 * sqrt (NSAMPLE * (1.0 / 3.0) * (2.0 / 3.0)));
}

void
create_theta_table_test_suite (void)
{
  CU_pSuite suite = CU_add_suite ("theta_table central source", suite_init, suite_teardown);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create theta_table test suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  if ((CU_add_test (suite, "read: total luminosity", test_read_total) == NULL) ||
      (CU_add_test (suite, "band luminosity", test_band_luminosity) == NULL) ||
      (CU_add_test (suite, "sampling, full band", test_sample_full_band) == NULL) ||
      (CU_add_test (suite, "sampling, one-group band", test_sample_one_group_band) == NULL))
  {
    fprintf (stderr, "Failed to add tests to theta_table suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
