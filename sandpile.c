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
static void render(char *restrict filename, int *restrict grid, size_t N) {
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

static void stabilize(int *restrict grid, size_t N, const size_t *restrict col_start, const size_t *restrict col_end, int *restrict q) {
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

void find_identity(int *restrict grid, size_t N, const size_t *restrict col_start, const size_t *restrict col_end, int *restrict q) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      grid[grid_index(i, j, N)] = 6;
    }
  }
  stabilize(grid, N, col_start, col_end, q);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      grid[grid_index(i, j, N)] = 6 - grid[grid_index(i, j, N)];
    }
  }
  stabilize(grid, N, col_start, col_end, q);
}

static void cond_fix_sink(size_t i, size_t j, size_t N, const size_t *restrict col_start, const size_t *restrict col_end, int *restrict grid) {
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

  int *grid = (int *)calloc((N + 2) * (N + 2), sizeof(int));
  int *q = (int *)calloc((N + 2) * (N + 2), sizeof(int));

  ptrdiff_t n = (ptrdiff_t)N;
  ptrdiff_t r = n / 2;

  size_t *col_start = (size_t *)malloc(N * sizeof(size_t));
  size_t *col_end = (size_t *)malloc(N * sizeof(size_t));

  for (ptrdiff_t i = 0; i < n; ++i) {
    ptrdiff_t width_2 = r * r - (i - r) * (i - r);
    if (width_2 < 0) {
      col_start[i] = 1;
      col_end[i] = 0;
      continue;
    }

    ptrdiff_t width = (ptrdiff_t)sqrt((double)width_2);

    col_start[i] = (size_t)MAX(0, r - width);
    col_end[i] = (size_t)MIN(n - 1, MIN(r + width, MAX(i, n - 1 - i)));
  }

  find_identity(grid, N, col_start, col_end, q);

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = col_start[i]; j <= col_end[i]; ++j) {
      cond_fix_sink(i - 1, j, N, col_start, col_end, grid);
      cond_fix_sink(i, j - 1, N, col_start, col_end, grid);
      cond_fix_sink(i + 1, j, N, col_start, col_end, grid);
      cond_fix_sink(i, j + 1, N, col_start, col_end, grid);
    }
  }

  render("identity.png", grid, N);

  free(grid);
  free(q);

  free(col_start);
  free(col_end);

  return EXIT_SUCCESS;
}
