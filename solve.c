#include "solve.h"
#define A(i,j)  (i*n + j)
#define eps 1.0e-19
#define eq(a,b) ((a-b<eps) && (b-a<eps))

extern int debug;
extern int restr;

enum FUNC formula(char* str)
{
    if(strcmp(str, "symm") == 0)
        return  symm;    
    if(strcmp(str, "positive_symm") == 0)
        return positive_symm;
    if(strcmp(str, "hilbert") == 0)
        return hilbert;
	if(strcmp(str, "upper") == 0)
		return upper;
	if(strcmp(str,"disg")==0)
		return disg;
    return errf;
}

void fill(enum FUNC xin,int n, double* a)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            switch(xin) {
            case symm:
                a[A(i,j)] = abs(i-j);
                break;
            case positive_symm:
                a[A(i,j)] = 1+abs(i-j);
                break;
            case hilbert:
                a[A(i,j)] = 1./(double)(i+j+1);
                break;
			case upper:
				if(i==j) a[A(i,j)] = 1;
				if(i>j) a[A(i,j)] = 0;
				if(i<j) a[A(i,j)] = -1;
				break;
			case disg:
				if(i>j)
					a[A(i,j)] = n-i;
				else 
					a[A(i,j)] = n-j;
			case errf:
				break;
              }
        }
    }
}
void printm(FILE* fout, int n, double* a, double* b)
{
  if(restr > n-2)
    {
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
               fprintf(fout, "%.2f ",a[i*n+j]);
            }
           fprintf(fout, "   %.2f\n", b[i]);
        }
    }
    else
    {
        for(int i = 0;i<restr;i++)
        {
            for(int j = 0; j<restr;j++)
            {
               fprintf(fout,"%.2f ",a[i*n+j]);
            }
           fprintf(fout,".. %.2f    %.2f\n",a[i*n + n -1],b[i]);
        }
       fprintf(fout,"..\n");
        for(int j = 0; j<restr;j++)
        {
            fprintf(fout,"%.2f ",a[(n-1)*n+j]);
        }
       fprintf(fout,".. %.2f    %.2f\n",a[(n-1)*n + n -1],b[n-1]);
       
    }

}

int solve(int n, double* a, double* b, double* x)
{   
    double sq, sphi, cphi, olx, oly;
    for(int j = 0;j<n-1;j++)
    {
        for(int i = j+1; i<n;i++)
        {
            sq = sqrt(a[A(j,j)]*a[A(j,j)] + a[A(i,j)]*a[A(i,j)]);
            if(eq(sq,0)) return 0;
            cphi = a[A(j,j)]/sq;
            sphi = - a[A(i,j)]/sq;
            for(int k = j; k<n;k++)
            {
                olx = a[A(j,k)];
                oly = a[A(i,k)];
                a[A(j,k)] = olx*cphi - oly*sphi;
                a[A(i,k)] = olx*sphi + oly*cphi;
            }
            olx = b[j];
            oly = b[i];
            b[j] = olx*cphi - oly*sphi;
            b[i] = olx*sphi + oly*cphi;
            
        }
    }
    if(debug)
    {
        fprintf(stderr,"\n\nUpper triangular form:\n");
	printm(stderr,n,a,b);
    }
    for(int i = n-1;i>=0;i--)
    {
        if(eq(a[A(i,i)],0)) return 0;
        x[i]=b[i]/a[A(i,i)];
        for(int j = 0; j < i; j++)
        {
            b[j]-=x[i]*a[A(j,i)];
        }
    }
    return 1;
}
