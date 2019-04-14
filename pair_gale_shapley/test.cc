#include <stdio.h>
#include <vector>
#include <math.h>

int main(){
  std::vector<int> test_vec;

  for (int i = 0; i < 1000000000; i++){
    float j = log(log(i));
    if (i%10000000 ==0)
      fprintf(stderr,"%f\n",j);
  }
  return 0;
}
