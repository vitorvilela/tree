#ifndef RANDOM_H
#define RANDOM_H

#define IA 16807
#define IM 2147483647
#define AM_statistic_randomOne (1.0/IM)
#define IQ 127773
#define IR 2836
#define NDIV_statistic_randomOne (1+(IM-1)/NTAB)

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 64
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


//HEADERS OF STATISTICS' FUNCTIONS - BEGIN

float statistic_randomOne(int *);
float statistic_randomTwo(long int *);
float statistic_gasdev(int *);

//HEADERS OF STATISTICS' FUNCTIONS - END


#endif
