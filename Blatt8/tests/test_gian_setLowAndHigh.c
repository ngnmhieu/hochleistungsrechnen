#include <stdio.h>
void
setLowAndHigh(int num_rows, int g_num_procs, int g_rank){

	// Split all rows equally [base_size]
	int base_size = num_rows / g_num_procs;
	// and give rest to first n <= rest
	int num_big = num_rows % g_num_procs;
	/*
	 * number of processes that receive one row more = rest => (0 --> rest - 1) => rank < num_big
	 * calculate number of processes that receive one more, prior to current process
	 * eg: num_big = 4 & rank = 3 ==> rank 0, 1, 2 have one more also --> #rank prior to this one
	 * eg: num_big = 4 & rank = 4 ==> rank 0, 1, 2, 3 have one more --> num_big prior to this one
	 */
	int num_big_pres = (g_rank < num_big) ? g_rank : num_big;
	/*
	 * number of processes that do not receive one more, prior to current process
	 * max(rank - num_big_pres, 0)
	 * eg: num_big = 4 & rank = 3 ==> 3 bigger prior => 3 - 3 = 0; no small prior
	 * eg: num_big = 4 & rank = 5 ==> 4 bigger prior (0,1,2,3) => 5 - 4 = 1; 1 smaller prior (4)
	 */
	int num_sml_pres = (g_rank - num_big_pres > 0) ? g_rank - num_big_pres : 0;
	// number of bigger prior to this * bigger-size + number of smaller prior * base-size
	int g_minMat = (num_big_pres * (base_size + 1)) + num_sml_pres * base_size;
	// if bigger: minIndex + bigSie; else: minIndex + smallSize
	int my_size = (g_rank < num_big) ? base_size + 1 : base_size;
	int g_maxMat = g_minMat + my_size - 1;
	printf ("=================\n");
	printf ("rows: %d\n", num_rows);
	printf ("rank: %d\n", g_rank);
	printf ("size: %d\n", my_size);
	printf ("> min: %d\n", g_minMat);
	printf ("> max: %d\n", g_maxMat);

}

int
main (int argc, char** argv)
{
	setLowAndHigh(10, 3, 0);
	setLowAndHigh(10, 3, 1);
	setLowAndHigh(10, 3, 2);

	printf ("~~~~~~~~~~~~~~~~~\n");
	printf ("~~~~~~~~~~~~~~~~~\n");

	setLowAndHigh(10, 1, 0);

	printf ("~~~~~~~~~~~~~~~~~\n");
	printf ("~~~~~~~~~~~~~~~~~\n");

	setLowAndHigh(40, 4, 0);
	setLowAndHigh(40, 4, 1);
	setLowAndHigh(40, 4, 2);
	setLowAndHigh(40, 4, 3);
}


