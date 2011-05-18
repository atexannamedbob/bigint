#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "bigint.h"

typedef bigint_t *(*bz_binop_t)(bigint_t *, bigint_t *);
typedef void (*bz_binop_inplace_t)(bigint_t *, bigint_t *, bigint_t *);
typedef bigint_t *(*bz_shift_t)(bigint_t *, uint32_t);
typedef int (*bz_compare_t)(bigint_t *, bigint_t *);

static const uint64_t upper_mask = ~0ull << RADIXBITS;
static const uint64_t lower_mask = ~0ull >> RADIXBITS;

/* Trim length to be actual length, not counting leading zeros. */
static bigint_t *trim(bigint_t *x) {
  uint32_t i = x->len - 1;

  do {
    if (x->a[i] == 0)
      x->len--;
    else
      break;
  } while (i-- != 0);

  return x;
}

bigint_t *bz_new(uint32_t cap) {
  bigint_t *x;
  uint32_t i;

  x = (bigint_t *)malloc(sizeof(bigint_t));
  assert(x != NULL);
  x->a = (digit_t *)malloc(cap * sizeof(uint32_t));
  assert(x->a != NULL);

  for (i = 0; i < cap; ++i)
    x->a[i] = 0;
  x->len = 0;
  x->cap = cap;
  x->neg = false;

  return x;
}

void bz_free(bigint_t *x) {
  if (x == NULL)
    return;

  free(x->a);
  free(x);
}

static void bz_copy(bigint_t *from, bigint_t *to) {
  uint32_t i;

  assert(to->cap >= from->cap);
  to->len = from->len;
  to->neg = from->neg;

  for (i = 0; i < from->len; ++i)
    to->a[i] = from->a[i];
}

bigint_t *bz_clone(bigint_t *x) {
  bigint_t *y;

  y = bz_new(x->cap);
  bz_copy(x, y);

  return y;
}

void bz_print(bigint_t *x) {
  uint32_t i = x->len - 1;

  if (x->len == 0)
    return;

  do {
    if (i != x->len - 1)
      printf(" ");
    printf("%08X", x->a[i]);
  } while (i-- != 0);
}

/* -1 = less than, 0 = equal, 1 = greater than. */
static int bz_scompare_abs(bigint_t *x, bigint_t *y) {
  uint32_t i, xi, yi;

  if (x->len < y->len)
    return -1;
  if (x->len > y->len)
    return 1;
  
  for (i = 0; i < x->len; ++i) {
    xi = x->a[i];
    yi = y->a[i];
    if (xi < yi)
      return 1;
    if (xi > yi)
      return -1;
  }

  return 0;
}

int bz_scompare(bigint_t *x, bigint_t *y) {
  if (x->neg && y->neg)
    return bz_scompare_abs(y, x);
  if (x->neg && !y->neg)
    return -1;
  if (!x->neg && y->neg)
    return 1;
  return bz_scompare_abs(x, y);
}

/* Returns the product of x * y. */
static void bz_smult_(bigint_t *x, bigint_t *y, bigint_t *p) {
  uint32_t i, j, n, t;
  digit_t c, u, v;
  ddigit_t inner;

  n = x->len;
  t = y->len;
  p->neg = x->neg != y->neg;
  p->len = n + t;
  assert(p->cap >= n+t);

  for (i = 0; i < t; ++i) {
    c = 0;
    for (j = 0; j < n; ++j) {
      /* Calculate the inner product. */
      inner = (ddigit_t)p->a[i+j] +
	((ddigit_t)x->a[j] * (ddigit_t)y->a[i]) +
	(ddigit_t)c;

      /* Get out the two 32-bit digits and update the product and carry. */
      u = (digit_t)((inner & upper_mask) >> RADIXBITS);
      v = (digit_t)(inner & lower_mask);

      p->a[i+j] = v;
      p->a[i+n] = u;
      c = u;
    }
  }

  trim(p);
}

bigint_t *bz_smult(bigint_t *x, bigint_t *y) {
  bigint_t *p;

  p = bz_new(x->len + y->len);
  bz_smult_(x, y, p);
  
  return p;
}

bigint_t *bz_pmult(bigint_t *x, bigint_t *y) {
  int i, nthreads;
  uint32_t j, n, t, d, l;
  bigint_t **partials, *p, *partial, *tmp;

  n = x->len;
  t = y->len;

  nthreads = omp_get_max_threads();
  omp_set_num_threads(nthreads);

  d = n / nthreads;
  partials = (bigint_t **)malloc(nthreads * sizeof(bigint_t));

#pragma omp parallel for private(tmp, partial, j) shared(partials)
  /* To forego some overhead we divide up the first multiplicand into
     as many chunks of digits as we have threads. */
  for (i = 0; i < nthreads; ++i) {
    /* The last chunk of digits gets the remainder. */
    l = i == nthreads-1 ? d + (n % nthreads) : d;

    tmp = bz_new(l);
    for (j = 0; j < l; ++j)
      tmp->a[j] = x->a[i*d+j];
    partial = bz_smult(tmp, y);
    partials[i] = bz_pshiftl(partial, i*l*RADIXBITS);
    partials[i]->neg = false;

    bz_free(partial);
    bz_free(tmp);
  }

  p = bz_new(n+t);
  p->neg = x->neg != y->neg;
  p->len = n + t;

  tmp = bz_new(n+t);
  for (i = 0; i < nthreads; ++i) {
    bz_padd_(p, partials[i], tmp);
    bz_copy(tmp, p);
    bz_free(partials[i]);
  }
  bz_free(tmp);
  free(partials);

  return trim(p);
}

/* Returns the quotient of x / y and stores the remainder into r. */
static void bz_div_(bigint_t *x, bigint_t *y, bigint_t *q, bigint_t *r,
		    bz_binop_t mult, bz_binop_inplace_t mult_,
		    bz_binop_t add,
		    bz_binop_t minus, bz_binop_inplace_t minus_,
		    bz_shift_t shiftl, bz_compare_t comp) {
  bigint_t *yshifted, *qit, *ytz, *x3, *yb, *qyb;
  bigint_t *tmp;
  uint32_t i, n, t;
  digit_t yt1, yt2;

  /* If x < y, then quotient is 0 with x as the remainder. */
  if ((*comp)(x, y) < 0) {
    if (r)
      bz_copy(x, r);

    assert(q->cap >= 1);
    q->len = 1;
    q->a[0] = 0;

    return;
  }

  /* Clone x because we're going to be updating it in-place. */
  x = bz_clone(x);
  n = x->len;
  t = y->len;
  q->neg = x->neg != y->neg;
  q->len = n-t-1;
  assert(q->cap >= n-t-1);
  if (r) {
    r->neg = false;
    r->len = t;
    assert(r->cap >= t);
  }

  /* 2 */
  yshifted = (*shiftl)(y, (n-t-1)*RADIXBITS);
  tmp = bz_new(x->len);
  while ((*comp)(x, yshifted) > 0) {
    q->a[n-t-1] += 1;
    (*minus_)(x, yshifted, tmp);
    bz_copy(tmp, x);
  }
  bz_free(tmp);
  bz_free(yshifted);

  /* 3 */
  yt1 = y->a[t-1]; yt2 = y->a[t-2];
  qit = bz_new(1); qit->len = 1;
  ytz = bz_new(2); ytz->len = 2;
  x3  = bz_new(3); x3->len = 3;

  for (i = n-1; i >= t; --i) {
    /* 3.1 */
    if (x->a[i] == yt1) {
      q->a[i-t-2] = ~0 - 1;
    } else {
      q->a[i-t-2] = ((ddigit_t)(x->a[i] << RADIXBITS) +
		     (ddigit_t)x->a[i-1]) / (ddigit_t)yt1;
    }

    /* 3.2 */
    qit->a[0] = q->a[i-t-2];
    ytz->a[1] = yt1;
    ytz->a[0] = yt2;
    x3->a[2] = x->a[i];
    x3->a[1] = x->a[i-1];
    x3->a[0] = x->a[i-2];

    tmp = (*mult)(qit, ytz);
    while ((*comp)(tmp, x3) > 0) {
      q->a[i-t-2] -= 1;
      qit->a[0] = q->a[i-t-2];
      (*mult_)(qit, ytz, tmp);
    }
    bz_free(tmp);

    /* 3.3 */
    qit->a[0] = q->a[i-t-2];
    yb = (*shiftl)(y, (i-t-2)*RADIXBITS);
    qyb = (*mult)(qit, tmp);
    tmp = bz_new(x->len+1);
    tmp = (*minus)(x, qyb);
    bz_copy(tmp, x);
    bz_free(qyb);
    bz_free(tmp);

    /* 3.4 */
    if (x->neg) {
      tmp = (*add)(x, yb);
      bz_copy(tmp, x);
      q->a[i-t-2] -= 1;
      bz_free(tmp);
    }

    bz_free(yb);
  }

  bz_free(qit);
  bz_free(ytz);
  bz_free(x3);

  if (r)
    bz_copy(x, r);

  trim(q);
  trim(r);
  bz_free(x);
}

#if 0
bigint_t *bz_sdiv(bigint_t *x, bigint_t *y, bigint_t **r) {
  bigint_t *q;

  q = bz_new(x->len - y->len - 1);
  *r = bz_new(y->len);
  bz_sdiv_(x, y, q, *r,
	   &bz_smult, &bz_smult_, &bz_sadd, &bz_sminus, &bz_sminus_,
	   &bz_sshiftl, &bz_scompare);

  return q;
}


bigint_t *bz_pdiv(bigint_t *x, bigint_t *y, bigint_t **r) {
  bigint_t *q;

  q = bz_new(x->len - y->len - 1);
  *r = bz_new(y->len);
  bz_sdiv_(x, y, q, *r,
	   &bz_pmult, &bz_pmult_, &bz_padd, &bz_pminus, &bz_pminus_,
	   &bz_pshiftl, &bz_pcompare);

  return q;
}

int main() {
  bigint_t *x, *y, *p;
  x = bz_new(1);
  y = bz_new(1);
  x->len = 1;
  y->len = 1;
  x->a[0] = UINT32_MAX;
  y->a[0] = UINT32_MAX;

  p = bz_smult(x, y);
  bz_print(p);
  printf("\n");
  printf("%016llX\n", (uint64_t)((uint64_t)UINT32_MAX * (uint64_t)UINT32_MAX));
  
  return 0;
}
#endif
