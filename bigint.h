#include <stdint.h>
#include <stdbool.h>

/* We use radix 2^32 digits. */
#define RADIXBITS 32
typedef uint32_t digit_t;
typedef uint64_t ddigit_t;

typedef struct {
  digit_t *a;
  uint32_t len;
  uint32_t cap;
  bool neg;
} bigint_t;

bigint_t *bz_new();
void bz_free(bigint_t *);
void bz_print(bigint_t *);

int bz_scompare(bigint_t *, bigint_t *);
int bz_pcompare(bigint_t *, bigint_t *);

bigint_t *bz_sadd(bigint_t *, bigint_t *);
bigint_t *bz_sminus(bigint_t *, bigint_t *);
bigint_t *bz_smult(bigint_t *, bigint_t *);
bigint_t *bz_sdiv(bigint_t *, bigint_t *, bigint_t **);
bigint_t *bz_sshiftl(bigint_t *, uint32_t);
bigint_t *bz_sshiftr(bigint_t *, uint32_t);
bigint_t *bz_sand(bigint_t *, bigint_t *);
bigint_t *bz_sor(bigint_t *, bigint_t *);
bigint_t *bz_sxor(bigint_t *, bigint_t *);
bigint_t *bz_snot(bigint_t *, bigint_t *);

bigint_t *bz_padd(bigint_t *, bigint_t *);
bigint_t *bz_pminus(bigint_t *, bigint_t *);
bigint_t *bz_pmult(bigint_t *, bigint_t *);
bigint_t *bz_pdiv(bigint_t *, bigint_t *, bigint_t **);
bigint_t *bz_pshiftl(bigint_t *, uint32_t);
bigint_t *bz_pshiftr(bigint_t *, uint32_t);
bigint_t *bz_pand(bigint_t *, bigint_t *);
bigint_t *bz_por(bigint_t *, bigint_t *);
bigint_t *bz_pxor(bigint_t *, bigint_t *);
bigint_t *bz_pnot(bigint_t *, bigint_t *);
