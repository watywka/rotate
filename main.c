#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "solve.h"

int debug;
int restr;

int main(int argc, char* argv[])
{
    int n;
    double  *a, *ta, *b, *tb, *x;
    int opt;
    struct timespec begin, end;
    FILE* fin = NULL;
    opterr = 0;
    restr = 5;
    debug = 0;
    if(argc == 1)
    {
        fprintf(stderr,"Usage: solve [options] \n   -f   file with the linear system of equations\n OR\n   -e   equation\n   -d   debug mode\n   -r   output restrictions\n");
        return 1;
    }
    while((opt = getopt(argc,argv,"df:r:")) != -1)
    {
        switch(opt) {
        case 'd':
            debug=1;
            break;
        case 'f':
            fin = fopen(optarg,"r");
            if(!fin) {
                fprintf(stderr,"File not found\n");
                return -1;
            }
            break;
        case 'r':
            if( sscanf(optarg,"%d",&restr) != 1) {
                fprintf(stderr,"Wrong option usage\n");
                return -1;
            }
            if(restr<1)
            {
                fprintf(stderr,"Restriction must be positive\n");
                return -1;
            }
            break;
        case '?':
            if((optopt == 'r') || (optopt == 'f')) {
                fprintf(stderr,"-%c requires an argument\n",optopt);
                return -1;
            }
            fprintf(stderr, "%c: unknown option\n", optopt);
            return -1;
            break;
        default:
            return -1;
        }
    }
    if(optind<argc)
    {
        fprintf(stderr,"Wrong options format\n");
		fclose(fin);
		return -1;
    } 
    if(fin)
    {
        if(fscanf(fin,"%d", &n)!=1)
        {
            fprintf(stderr,"Corrupted file\n");
			fclose(fin);
            return -1;
        }
        if(n<1)
        {
            fprintf(stderr,"Wrong matrix size\n");
			fclose(fin);
            return -1;
        }
        if(!(a=(double*) malloc(sizeof(double)*n*n)))
        {
            fprintf(stderr,"Can't allocate memory\n");
            return -1;
        }
        if(!(ta =(double*) malloc(sizeof(double)*n*n)))
        {
            fprintf(stderr,"Can't allocate memory\n");
        }
        if(!(b =(double*) malloc(sizeof(double)*n)))
        {
            fprintf(stderr,"Can't allocate memory\n");
        }
        if(!(tb =(double*) malloc(sizeof(double)*n)))
        {
            fprintf(stderr,"Can't allocate memory\n");
        }
        if(!(x =(double*) malloc(sizeof(double)*n)))
        {
            fprintf(stderr,"Can't allocate memory\n");
        }
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                if(fscanf(fin,"%lf", &a[i*n+j])!=1)
                {
                    fprintf(stderr,"Corrupted file\n");
					fclose(fin);
                    free(ta);
                    free(tb);
                    free(a);
                    free(b);
                    free(x);
                    return -1;
                }
                ta[i*n+j]=a[i*n+j];
            }
            if(fscanf(fin,"%lf", &b[i])!=1) {
                fprintf(stderr,"Corrupted file\n");
				fclose(fin);
                free(ta);
                free(tb);
                free(a);
                free(b);
                free(x);
                return -1;
            }
            tb[i]=b[i];
        }
	fclose(fin);
    }
    printf("Equation:\n");
    printm(stdout,n,a,b);
    clock_gettime(CLOCK_MONOTONIC,&begin);
    if(solve(n,a,b,x))
	{ 
        double eps=0;
        clock_gettime(CLOCK_MONOTONIC, &end);
        printf("\nTime spent:%fns\n",((end.tv_sec-begin.tv_sec)+(double)(end.tv_nsec-begin.tv_nsec)/100000000));
        printf("\n\nSolution:\n");
        if (restr > n-2)
        {
            for(int i=0;i<n;i++)
            {
                printf("%.2f\n",x[i]);
            }
        }
        else
        {
            for(int i=0;i<restr;i++)
            {
                printf("%.2f\n",x[i]);
            }
            printf("..\n");
            printf("%.2f\n",x[n-1]);
        }

        for(int i=0;i<n;i++)
        {
            double sq=0;
            for(int j=0;j<n;j++)
            {
                sq+=ta[i*n+j]*x[j];
            }
            sq-=tb[i];
            eps+=sq*sq;
        }
        printf("Residual : %f\n", sqrt(eps));
    }
    else
    {
        fprintf(stderr,"Unable to solve the system\n");
    }
    free(ta);
    free(tb);
    free(a);
    free(b);
    free(x);
    return 1;

}
