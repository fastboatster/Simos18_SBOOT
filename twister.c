#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

// This is the ``Mersenne Twister'' random number generator MT19937
// http://www.math.keio.ac.jp/~matumoto/ver980409.html
// with inspiration from code attributed to Cokus@math.washington.edu

#define N (624)                              // length of state vector
#define M (397)                              // a period parameter
#define K (0x9908B0DFU)                      // a magic constant
#define hiBit(u) ((u)&0x80000000U)           // mask all but highest   bit of u
#define loBit(u) ((u)&0x00000001U)           // mask all but lowest    bit of u
#define loBits(u) ((u)&0x7FFFFFFFU)          // mask     the highest   bit of u
#define mixBits(u, v) (hiBit(u) | loBits(v)) // move hi bit of u to hi bit of v

static inline uint32_t temperMT(uint32_t s1)
{
  s1 ^= (s1 >> 11);
  s1 ^= (s1 << 7) & 0x9D2C5680U;
  s1 ^= (s1 << 15) & 0xEFC60000U;
  return (s1 ^ (s1 >> 18));
}

void seedMT(uint32_t seed, uint32_t *output)
{
  uint32_t state[N + 1];
  uint32_t *next;

  uint32_t x = (seed | 1U) & 0xFFFFFFFFU, *s = state;
  int j;

  for (*s++ = x, j = N; --j;
       *s++ = (x *= 69069U) & 0xFFFFFFFFU)
    ;
  uint32_t *p0 = state, *p2 = state + 2, *pM = state + M, s0, s1;

  next = state + 1;

  for (s0 = state[0], s1 = state[1], j = N - M + 1; --j; s0 = s1, s1 = *p2++)
    *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

  for (pM = state, j = M; --j; s0 = s1, s1 = *p2++)
    *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

  s1 = state[0], *p0 = *pM ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
  output[0] = temperMT(s1);
  uint32_t y;
  // generate 32 uints:
  for (int i = 1; i < 32; i++)
  {
    y = *next++;
    output[i] = temperMT(y);
  }
}

static inline uint32_t bswap32(uint32_t x)
{
  return (((x & 0xff000000u) >> 24) |
          ((x & 0x00ff0000u) >> 8) |
          ((x & 0x0000ff00u) << 8) |
          ((x & 0x000000ffu) << 24));
}

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    return 1;
  }
  uint32_t seed = strtol(argv[1], NULL, 16);
  uint64_t match = strtoll(argv[2], NULL, 16);
  int suppress_output = 0;
  if (argc > 3)
  {
    suppress_output = 1;
  }

  mpz_t n;
  // init RSA 1024 public key with id 0x96
  // hopefully, the endianness is correct here
  // and it's not 099122370e5bb81c8ee2381ccda9b10f5c48d50d6c1e275620fc6b1ddb4cb45b08fa041cc365367c436732d26e8ce291de920c84bb2b32aca047f6bc8a56059614f91cf588167e5dcd4e31d0a7a2cda76845061df7c490363e2acb0503ab68bc74744dee34ab587829eac8f71dbd2cb9f4b0aa50ca6dd2edbc312c4f13fed8e1
  mpz_init_set_str(n, "e1d8fe134f2c31bcedd26dca50aab0f4b92cbd1df7c8ea297858ab34ee4d7474bc68ab0305cb2a3e3690c4f71d064568a7cda2a7d0314ecd5d7e1688f51cf9149605568abcf647a0ac322bbb840c92de91e28c6ed23267437c3665c31c04fa085bb44cdb1d6bfc2056271e6c0dd5485c0fb1a9cd1c38e28e1cb85b0e37229109", 16);
  mpz_t e;
  mpz_init_set_ui(e, 65537U);

  int done = 0;
#pragma omp parallel shared(seed, done, match, n, e, suppress_output)
  while (!done)
  {
    int j = 0;
    uint32_t current_seed = 0;
#pragma omp atomic capture
    {
      current_seed = seed;
      seed = seed + 2;
    }
    // an arr of 32 32-bit uints
    uint32_t rand_data[32];
    seedMT(current_seed, rand_data);
    unsigned char *rand_data_bytes = (unsigned char *)rand_data;

    // The last int has the high bits replaced with 0x200, presumably to make the resulting bigint valid within the RSA parameters.
    rand_data[31] = bswap32(bswap32(rand_data[31] & 0xFFFF) + 0x0200);

    // This byte is set to 0x00 to serve as a PKCS v1.5. separator between the message and pseudorandom padding
    // 0x75 = 117
    rand_data_bytes[117] = 0;
    // DONE: check that all the data bytes in the "padding" area are not 0, otherwise replace with 0x01.
    //  not sure how important this really is, perhaps this happens very infrequently
    int k;
    // go over the padding bytes and replace any 0 bytes with 0x01
    for (k=118; k < 126; k++ )
      {
        // check if any bytes in padding area are 0:
          if (rand_data_bytes[k] == 0) {
            //
              if (!suppress_output)
              {
                  printf("**** Detected 0x00 byte in the padding area, replacing with 0x01 ****\n");
              };
              rand_data_bytes[k] = 0x01
          };
      };

    mpz_t data_num;
    mpz_init(data_num);
    mpz_import(data_num, 32, -1, 4, -1, 0, rand_data_bytes);
    mpz_t output;
    mpz_init(output);
    mpz_powm(output, data_num, e, n);
    unsigned char rsa_output[128];
    mpz_export(rsa_output, NULL, -1, 4, -1, 0, output);
    mpz_clear(output);
    mpz_clear(data_num);

    uint32_t *rsa_output_ints = (uint32_t *)rsa_output;
    if (rsa_output_ints[0] == match)
    {
      done = 1;
      if (!suppress_output)
      {
        printf("**** FOUND ****\n");
        printf("Seed: %08X\n", current_seed);
        printf("\nKey Data: \n");
      }
      for (j = 0; j < 32; j++)
      {
        printf("%08X", rand_data[j]);
      }
      if (!suppress_output)
      {

        printf("\nSeed Data: \n");
        for (j = 0; j < 32; j++)
        {
          printf(" %08X%s", rsa_output_ints[j], (j % 4) == 4 ? " " : "");
        }
        printf("\n");
      }
    }
  }
  mpz_clear(e);
  mpz_clear(n);

  return (EXIT_SUCCESS);
}
