/* Read the MatrixMarket count matrix of alevin-fry (rows = cells,
 * columns = features) into the slots of a features x cells dgCMatrix.
 *
 * The function reads the file two times. Pass 1 counts the entries of each
 * cell, which gives the slot p. Pass 2 puts each entry at its final position
 * in the slots i and x. Thus p, i and x are allocated one time, and the
 * memory is approximately 12 bytes for each entry.
 *
 * In alevin-fry files, the entries of a cell can be in more than one run,
 * and the runs are not in cell order. Within a run, the feature indices
 * increase. After pass 2, a cell with features that do not increase is
 * sorted. Two entries with the same cell and feature give an error.
 *
 * The C code does not open the file. It calls an R function that reads the
 * next block of bytes from an R connection (gzfile(), so that a compressed
 * file is also read), and it keeps the incomplete last line of a block for
 * the next block. */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define MTX_LINE 1024

typedef struct {
  const char *path;    /* for the messages only */
  SEXP call_open;      /* read(TRUE): open the file again, return the first block */
  SEXP call_next;      /* read(FALSE): return the next block, raw(0) at the end */
  PROTECT_INDEX ipx;   /* protects the current block */
  const char *data;    /* the bytes of the current block */
  R_xlen_t n, pos;     /* the length of the block, the position of the next byte */
  int eof;
  long long line;      /* the number of the line in buf */
  char buf[MTX_LINE];
} mtx_file;

typedef struct {
  int i;
  double x;
} mtx_pair;

/* Evaluate 'call', which returns the next block of the file. Return 0 at the
 * end of the file. */
static int get_block(mtx_file *m, SEXP call) {
  SEXP b = Rf_eval(call, R_GlobalEnv);
  REPROTECT(b, m->ipx);
  if (TYPEOF(b) != RAWSXP)
    Rf_error("the block reader of %s did not return a raw vector", m->path);
  m->data = (const char *) RAW(b);
  m->n = XLENGTH(b);
  m->pos = 0;
  m->eof = (m->n == 0);
  return !m->eof;
}

/* Open the file again and go to its start. */
static void start_pass(mtx_file *m) {
  get_block(m, m->call_open);
  m->line = 0;
}

/* Read one line into m->buf. Return 0 at the end of the file.
 * A data line must fit in MTX_LINE - 1 characters and must not have a NUL
 * character. For a comment line (it starts with '%'), the remainder of a
 * long line is ignored. */
static int read_line(mtx_file *m) {
  size_t len = 0;
  int any = 0, skip = 0;
  for (;;) {
    if (m->pos >= m->n && (m->eof || !get_block(m, m->call_next))) break;
    const char *s = m->data + m->pos;
    size_t avail = (size_t) (m->n - m->pos);
    const char *nl = (const char *) memchr(s, '\n', avail);
    size_t take = nl != NULL ? (size_t) (nl - s) + 1 : avail;
    any = 1;
    if (!skip) {
      if (len + take > MTX_LINE - 1) {
        if ((len > 0 ? m->buf[0] : s[0]) != '%')
          Rf_error("%s: line %lld is too long", m->path, m->line + 1);
        m->buf[0] = '%';
        len = 1;
        skip = 1;
      } else {
        memcpy(m->buf + len, s, take);
        len += take;
      }
    }
    m->pos += (R_xlen_t) take;
    if (nl != NULL) break;
  }
  if (!any) return 0;
  m->line++;
  m->buf[len] = '\0';
  if (m->buf[0] != '%' && strlen(m->buf) != len)
    Rf_error("%s: line %lld has a NUL character", m->path, m->line);
  return 1;
}

static const char *skip_space(const char *c) {
  while (*c == ' ' || *c == '\t' || *c == '\r' || *c == '\n') c++;
  return c;
}

static int is_end(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\0';
}

static int is_blank(const char *c) {
  return *skip_space(c) == '\0';
}

/* Remove the line end from m->buf, for an error message. */
static const char *line_text(mtx_file *m) {
  m->buf[strcspn(m->buf, "\r\n")] = '\0';
  return m->buf;
}

/* Parse a non-negative integer of the size line. Return 0 if the token is
 * not an integer. A value larger than LLONG_MAX / 10 stays larger than
 * INT_MAX, which the caller refuses. */
static int parse_count(const char **pc, long long *out) {
  const char *c = skip_space(*pc);
  if (*c < '0' || *c > '9') return 0;
  long long v = 0;
  while (*c >= '0' && *c <= '9') {
    if (v < LLONG_MAX / 10) v = v * 10 + (*c - '0');
    c++;
  }
  if (!is_end(*c)) return 0;
  *out = v;
  *pc = c;
  return 1;
}

/* Parse an index. Stop with an error if it is not an integer in 1..max. */
static int parse_index(mtx_file *m, const char **pc, int max, const char *what) {
  const char *c = skip_space(*pc);
  int neg = 0;
  if (*c == '-' || *c == '+') {
    neg = (*c == '-');
    c++;
  }
  if (*c < '0' || *c > '9')
    Rf_error("%s: line %lld does not have 3 numbers", m->path, m->line);
  long long v = 0;
  while (*c >= '0' && *c <= '9') {
    if (v <= INT_MAX) v = v * 10 + (*c - '0');
    c++;
  }
  if (!is_end(*c))
    Rf_error("%s: line %lld does not have 3 numbers", m->path, m->line);
  *pc = c;
  if (neg || v < 1 || v > max)
    Rf_error("%s has a %s index that is not in 1..%d (line %lld)",
             m->path, what, max, m->line);
  return (int) v;
}

/* Parse a value as scan() does. A token of at most 15 digits is converted
 * directly, because this conversion is exact. Other tokens go to R_strtod,
 * the parser of scan(). */
static double parse_value(mtx_file *m, const char **pc) {
  const char *c = skip_space(*pc);
  const char *s = c;
  long long v = 0;
  while (*c >= '0' && *c <= '9' && c - s < 15) {
    v = v * 10 + (*c - '0');
    c++;
  }
  if (c > s && is_end(*c)) {
    *pc = c;
    return (double) v;
  }
  char *end;
  double d = R_strtod(s, &end);
  if (end == s || !is_end(*end))
    Rf_error("%s: line %lld does not have 3 numbers", m->path, m->line);
  *pc = end;
  return d;
}

/* Read the banner, the comment lines and the size line. */
static void read_header(mtx_file *m, int *nr, int *nc, int *nnz) {
  if (!read_line(m))
    Rf_error("%s is empty. It is not a MatrixMarket file.", m->path);
  char banner[MTX_LINE];
  strcpy(banner, m->buf);
  char *t[6];
  int nt = 0;
  for (char *tok = strtok(banner, " \t\r\n"); tok != NULL && nt < 6;
       tok = strtok(NULL, " \t\r\n"))
    t[nt++] = tok;
  /* as in Matrix::readMM, the first token must be exact; the other tokens
   * can have upper case letters */
  for (int k = 1; k < nt; k++)
    for (char *c = t[k]; *c; c++)
      if (*c >= 'A' && *c <= 'Z') *c = (char) (*c - 'A' + 'a');
  if (nt != 5 || strcmp(t[0], "%%MatrixMarket") || strcmp(t[1], "matrix") ||
      strcmp(t[2], "coordinate") || (strcmp(t[3], "real") && strcmp(t[3], "integer")) ||
      strcmp(t[4], "general"))
    Rf_error("%s is not a MatrixMarket 'coordinate real general' file. The first line is: %s",
             m->path, line_text(m));
  for (;;) {
    if (!read_line(m))
      Rf_error("%s does not have a MatrixMarket size line", m->path);
    if (m->buf[0] != '%' && !is_blank(m->buf)) break;
  }
  const char *c = m->buf;
  long long d[3];
  for (int k = 0; k < 3; k++)
    if (!parse_count(&c, &d[k]))
      Rf_error("%s does not have a correct MatrixMarket size line: %s", m->path, line_text(m));
  if (!is_blank(c) || d[0] > INT_MAX || d[1] > INT_MAX)
    Rf_error("%s does not have a correct MatrixMarket size line: %s", m->path, line_text(m));
  if (d[2] > INT_MAX)
    Rf_error("%s has %lld entries. A dgCMatrix can have at most %d entries.",
             m->path, d[2], INT_MAX);
  *nr = (int) d[0];
  *nc = (int) d[1];
  *nnz = (int) d[2];
}

static int cmp_pair(const void *a, const void *b) {
  int ia = ((const mtx_pair *) a)->i, ib = ((const mtx_pair *) b)->i;
  return (ia > ib) - (ia < ib);
}

static SEXP read_mtx(mtx_file *m) {
  int nr, nc, nnz;
  start_pass(m);
  read_header(m, &nr, &nc, &nnz);

  /* pass 1: the number of entries of each cell, and the checks of each line */
  int *fill = (int *) R_alloc((size_t) nr + 1, sizeof(int));
  memset(fill, 0, ((size_t) nr + 1) * sizeof(int));
  long long n = 0;
  while (read_line(m)) {
    if (is_blank(m->buf)) continue;
    if ((m->line & 0xFFFFF) == 0) R_CheckUserInterrupt();
    const char *c = m->buf;
    int r = parse_index(m, &c, nr, "row");
    parse_index(m, &c, nc, "column");
    parse_value(m, &c);
    if (!is_blank(c))
      Rf_error("%s: line %lld does not have 3 numbers", m->path, m->line);
    if (++n > nnz)
      Rf_error("%s has more entries than its size line gives (%d)", m->path, nnz);
    fill[r - 1]++;
  }
  if (n != nnz)
    Rf_error("%s has %lld entries, but its size line gives %d entries", m->path, n, nnz);

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP sp = Rf_allocVector(INTSXP, (R_xlen_t) nr + 1);
  SET_VECTOR_ELT(res, 1, sp);
  SEXP si = Rf_allocVector(INTSXP, nnz);
  SET_VECTOR_ELT(res, 0, si);
  SEXP sx = Rf_allocVector(REALSXP, nnz);
  SET_VECTOR_ELT(res, 2, sx);
  SEXP sd = Rf_allocVector(INTSXP, 2);
  SET_VECTOR_ELT(res, 3, sd);
  int *p = INTEGER(sp), *ip = INTEGER(si);
  double *xp = REAL(sx);
  INTEGER(sd)[0] = nc;
  INTEGER(sd)[1] = nr;
  p[0] = 0;
  for (int r = 0; r < nr; r++) {
    p[r + 1] = p[r] + fill[r];
    fill[r] = p[r];
  }

  /* pass 2: put each entry at its final position */
  start_pass(m);
  int nr2, nc2, nnz2;
  read_header(m, &nr2, &nc2, &nnz2);
  if (nr2 != nr || nc2 != nc || nnz2 != nnz)
    Rf_error("%s changed while it was read", m->path);
  n = 0;
  while (read_line(m)) {
    if (is_blank(m->buf)) continue;
    if ((m->line & 0xFFFFF) == 0) R_CheckUserInterrupt();
    if (++n > nnz)
      Rf_error("%s changed while it was read", m->path);
    const char *c = m->buf;
    int r = parse_index(m, &c, nr, "row");
    int g = parse_index(m, &c, nc, "column");
    double v = parse_value(m, &c);
    int k = fill[r - 1]++;
    if (k >= p[r])
      Rf_error("%s changed while it was read", m->path);
    ip[k] = g - 1;
    xp[k] = v;
  }
  if (n != nnz)
    Rf_error("%s changed while it was read", m->path);

  /* in each cell, the feature indices must increase; sort a cell only if
   * necessary */
  mtx_pair *tmp = NULL;
  for (int r = 0; r < nr; r++) {
    int a = p[r], b = p[r + 1];
    int sorted = 1;
    for (int k = a + 1; k < b; k++) {
      if (ip[k] <= ip[k - 1]) {
        sorted = 0;
        break;
      }
    }
    if (sorted) continue;
    if (tmp == NULL) {
      int maxlen = 0;
      for (int q = 0; q < nr; q++)
        if (p[q + 1] - p[q] > maxlen) maxlen = p[q + 1] - p[q];
      tmp = (mtx_pair *) R_alloc((size_t) maxlen, sizeof(mtx_pair));
    }
    for (int k = a; k < b; k++) {
      tmp[k - a].i = ip[k];
      tmp[k - a].x = xp[k];
    }
    qsort(tmp, (size_t) (b - a), sizeof(mtx_pair), cmp_pair);
    for (int k = a; k < b; k++) {
      ip[k] = tmp[k - a].i;
      xp[k] = tmp[k - a].x;
      if (k > a && ip[k] == ip[k - 1])
        Rf_error("%s has more than one entry for row %d, column %d",
                 m->path, r + 1, ip[k] + 1);
    }
  }
  UNPROTECT(1);
  return res;
}

/* .Call entry point. 'read' is an R function: read(TRUE) opens the file
 * (again) and returns its first block of bytes, read(FALSE) returns the next
 * block, and raw(0) at the end of the file. 'path' is the file name for the
 * messages. Return list(i, p, x, Dim) of the features x cells dgCMatrix. */
SEXP readFryMtx(SEXP read, SEXP path) {
  if (!Rf_isFunction(read))
    Rf_error("'read' must be a function");
  if (!Rf_isString(path) || Rf_xlength(path) != 1 || STRING_ELT(path, 0) == NA_STRING)
    Rf_error("'path' must be one file name");
  mtx_file *m = (mtx_file *) R_alloc(1, sizeof(mtx_file));
  m->path = Rf_translateChar(STRING_ELT(path, 0));
  m->call_open = PROTECT(Rf_lang2(read, Rf_ScalarLogical(TRUE)));
  m->call_next = PROTECT(Rf_lang2(read, Rf_ScalarLogical(FALSE)));
  PROTECT_WITH_INDEX(R_NilValue, &m->ipx);
  m->data = NULL;
  m->n = m->pos = 0;
  m->eof = 0;
  m->line = 0;
  SEXP res = read_mtx(m);
  UNPROTECT(3);
  return res;
}

static const R_CallMethodDef callMethods[] = {
  {"readFryMtx", (DL_FUNC) &readFryMtx, 2},
  {NULL, NULL, 0}
};

void R_init_fishpond(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
