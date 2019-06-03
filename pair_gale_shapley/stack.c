#include <stdio.h>
int power(int base, int exp);
int main()
{
  int i;
#pragma omp parallel for
  for (i=0; i<10; ++i)
    printf("%d %d %d\n", i, power(2,i), power(-3,i));

  return 0;
}

int power(int base, int exp)
{
  int result = 1;
  while (exp)
    {
      if (exp & 1)
	result *= base;
      exp >>= 1;
      base *= base;
    }

  return result;
}
