#include <gtest/gtest.h>
#include <fstream>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <thrust/system/omp/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/for_each.h>

TEST(OMP,NumThreads){
  /* Fork a team of threads giving them their own copies of variables */
  #pragma omp parallel
  {
    /* Obtain thread number */
    int tid = omp_get_thread_num();
    printf("Hello World from thread = %d\n", tid);

    /* Only master thread does this */
    if (tid == 0)
      printf("Number of threads = %d\n", omp_get_num_threads());
      EXPECT_NE(omp_get_num_threads(), 1) << "Only one thread, not running in parallel environment";

  }  /* All threads join master thread and disband */
}


TEST(OMP, Thrust){
  auto printer = [](int i) {
    printf("Hello on element %2d from thread %2d\n", i, omp_get_thread_num());
    EXPECT_NE(omp_get_num_threads(), 1) << "Only one thread, not running in parallel environment";
  };
  using enumerator = thrust::counting_iterator<int>;
  thrust::for_each(thrust::omp::par, enumerator(0), enumerator(100), printer);
  EXPECT_EQ(omp_get_num_threads(), 1) << "More theads than expected"; //number of threads is one
}