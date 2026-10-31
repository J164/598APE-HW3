#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

static void die(const char *s) {
  fprintf(stderr, "%s\n", s);
  exit(1);
}

static size_t grid_index(size_t i, size_t j, size_t N) {
  return (i + 1) * (N + 2) + (j + 1);
}

// you can also use this function to debug print your sandpile
static void render(char *restrict filename, int8_t *restrict grid, size_t N) {
  size_t chan = 3;
  unsigned char *pixels = (unsigned char *)malloc(N * N * chan * sizeof(char));
  size_t stride = N * chan;
  if (!pixels)
    die("malloc failed for PNG pixel buffer");
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      unsigned char r = 0, g = 0, b = 0;
      switch (grid[grid_index(i, j, N)]) {
      case 0: // obsidian
        r = 15;
        g = 15;
        b = 25;
        break;
      case 1: // electric cyan
        r = 0;
        g = 210;
        b = 255;
        break;
      case 2: // vibrant magenta
        r = 230;
        g = 30;
        b = 130;
        break;
      case 3: // bright gold
        r = 255;
        g = 200;
        b = 0;
        break;
      default: // white
        r = 0;
        g = 0;
        b = 0;
      }

      size_t pid = (i * N + j) * chan;
      pixels[pid] = r;
      pixels[pid + 1] = g;
      pixels[pid + 2] = b;
    }
  }
  if (stbi_write_png(filename, (int)N, (int)N, (int)chan, pixels, (int)stride)) {
    printf("Successfully saved sandpile identity to: %s\n", filename);
  } else {
    die("Error: stb_image_write failed to save png file\n");
  }
  free(pixels);
}

static void cg_solve(size_t N, const size_t *restrict col_start, const size_t *restrict col_end, const double *restrict b, double *restrict x, double tol, double *restrict r, double *restrict p, double *restrict Ap) {
  double rr = 0.0;

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = col_start[i]; j <= col_end[i]; ++j) {
      size_t index = grid_index(i, j, N);

      r[index] = b[index];
      p[index] = b[index];
      x[index] = 0.0;
      rr += b[index] * b[index];
    }
  }

  double stop = tol * tol * rr;

  while (rr > stop) {
    double pAp = 0.0;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = col_start[i]; j <= col_end[i]; ++j) {
        size_t index = grid_index(i, j, N);

        Ap[index] = 4 * p[index] - p[grid_index(i - 1, j, N)] - p[grid_index(i, j - 1, N)] - p[grid_index(i + 1, j, N)] - p[grid_index(i, j + 1, N)];
        pAp += p[index] * Ap[index];
      }
    }

    double alpha = rr / pAp;
    double rr_new = 0.0;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = col_start[i]; j <= col_end[i]; ++j) {
        size_t index = grid_index(i, j, N);

        x[index] += alpha * p[index];
        r[index] -= alpha * Ap[index];
        rr_new += r[index] * r[index];
      }
    }

    double beta = rr_new / rr;
    rr = rr_new;

    for (size_t i = 0; i < N; ++i) {
      for (size_t j = col_start[i]; j <= col_end[i]; ++j) {
        size_t index = grid_index(i, j, N);

        p[index] = r[index] + beta * p[index];
      }
    }
  }
}

static void stabilize_by_toppling(int8_t *restrict grid, size_t N, const size_t *restrict col_start, const size_t *restrict col_end, int8_t *restrict q) {
  int changed = 0;
  do {
    changed = 0;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = col_start[i]; j <= col_end[i]; ++j) {
        q[grid_index(i, j, N)] = grid[grid_index(i, j, N)] >> 2;
      }
    }

    for (size_t i = 0; i < N; ++i) {
      for (size_t j = col_start[i]; j <= col_end[i]; ++j) {
        int new = (grid[grid_index(i, j, N)] & 3) + q[grid_index(i - 1, j, N)] + q[grid_index(i, j - 1, N)] + q[grid_index(i + 1, j, N)] + q[grid_index(i, j + 1, N)];
        grid[grid_index(i, j, N)] = new;
        changed |= (new > 3);
      }
    }
  } while (changed);
}

static void stabilize(int8_t *restrict grid, size_t N, const size_t *restrict col_start, const size_t *restrict col_end, int8_t *restrict q, int *restrict t, double *restrict b, double *restrict x, double *restrict r, double *restrict p, double *restrict Ap) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = col_start[i]; j <= col_end[i]; ++j) {
      size_t index = grid_index(i, j, N);

      b[index] = (double)grid[index] - 3.0;
    }
  }

  cg_solve(N, col_start, col_end, b, x, 1e-6, r, p, Ap);

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = col_start[i]; j <= col_end[i]; ++j) {
      size_t index = grid_index(i, j, N);

      t[index] = MAX(0, (int)floor(x[index]));
    }
  }

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = col_start[i]; j <= col_end[i]; ++j) {
      int new = grid[grid_index(i, j, N)] - 4 * t[grid_index(i, j, N)] + t[grid_index(i - 1, j, N)] + t[grid_index(i, j - 1, N)] + t[grid_index(i + 1, j, N)] + t[grid_index(i, j + 1, N)];
      grid[grid_index(i, j, N)] = (int8_t)new;
    }
  }

  stabilize_by_toppling(grid, N, col_start, col_end, q);
}

void find_identity(int8_t *restrict grid, size_t N, const size_t *restrict col_start, const size_t *restrict col_end, int8_t *restrict q, int *restrict t, double *restrict b, double *restrict x, double *restrict r, double *restrict p, double *restrict Ap) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      grid[grid_index(i, j, N)] = 6;
    }
  }
  stabilize(grid, N, col_start, col_end, q, t, b, x, r, p, Ap);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      grid[grid_index(i, j, N)] = 6 - grid[grid_index(i, j, N)];
    }
  }
  stabilize(grid, N, col_start, col_end, q, t, b, x, r, p, Ap);
}

static void cond_fix_sink(size_t i, size_t j, size_t N, const size_t *restrict col_start, const size_t *restrict col_end, int8_t *restrict grid) {
  size_t index = grid_index(i, j, N);

  if (i < N && j < N && (j < col_start[i] || j > col_end[i])) { 
    grid[index] = 4;
  }
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 2) {
    fprintf(stderr, "Usage: %s N\n", argv[0]);
    return EXIT_FAILURE;
  }
  char *end;
  uintmax_t parsed = strtoumax(argv[1], &end, 10);
  if (errno || argv[1][0] < '0' || argv[1][0] > '9' || *end || !parsed ||
      parsed > INT32_MAX)
    die("grid size must be a positive int");
  size_t N = (size_t)parsed;

  ptrdiff_t n = (ptrdiff_t)N;
  ptrdiff_t rad = n / 2;

  size_t *col_start = (size_t *)malloc(N * sizeof(size_t));
  size_t *col_end = (size_t *)malloc(N * sizeof(size_t));

  for (ptrdiff_t i = 0; i < n; ++i) {
    ptrdiff_t width_2 = rad * rad - (i - rad) * (i - rad);
    if (width_2 < 0) {
      col_start[i] = 1;
      col_end[i] = 0;
      continue;
    }

    ptrdiff_t width = (ptrdiff_t)sqrt((double)width_2);

    col_start[i] = (size_t)MAX(0, rad - width);
    col_end[i] = (size_t)MIN(n - 1, MIN(rad + width, MAX(i, n - 1 - i)));
  }

  int8_t *grid = (int8_t *)calloc((N + 2) * (N + 2), sizeof(int8_t));
  int8_t *q = (int8_t *)calloc((N + 2) * (N + 2), sizeof(int8_t));
  int *t = (int *)calloc((N + 2) * (N + 2), sizeof(int));
  double *b = (double *)calloc((N + 2) * (N + 2), sizeof(double));
  double *x = (double *)calloc((N + 2) * (N + 2), sizeof(double));
  double *r = (double *)calloc((N + 2) * (N + 2), sizeof(double));
  double *p = (double *)calloc((N + 2) * (N + 2), sizeof(double));
  double *Ap = (double *)calloc((N + 2) * (N + 2), sizeof(double));

  find_identity(grid, N, col_start, col_end, q, t, b, x, r, p, Ap);

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = col_start[i]; j <= col_end[i]; ++j) {
      cond_fix_sink(i - 1, j, N, col_start, col_end, grid);
      cond_fix_sink(i, j - 1, N, col_start, col_end, grid);
      cond_fix_sink(i + 1, j, N, col_start, col_end, grid);
      cond_fix_sink(i, j + 1, N, col_start, col_end, grid);
    }
  }

  free(col_start);
  free(col_end);

  free(q);
  free(t);
  free(b);
  free(x);
  free(r);
  free(p);
  free(Ap);

  render("identity.png", grid, N);

  free(grid);

  return EXIT_SUCCESS;
}
