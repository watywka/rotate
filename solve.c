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

void fill(enum FUNC xin,int n, long double* a)
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
                a[A(i,j)] = 1./(long double)(i+j+1);
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
void printm(FILE* fout, int n, long double* a, long double* b)
{
  if(restr > n-2)
    {
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
               fprintf(fout, "%.4Lf ",a[i*n+j]);
            }
           fprintf(fout, "   %.4Lf\n", b[i]);
        }
    }
    else
    {
        for(int i = 0;i<restr;i++)
        {
            for(int j = 0; j<restr;j++)
            {
               fprintf(fout,"%.2Lf ",a[i*n+j]);
            }
           fprintf(fout,".. %.2Lf    %.2Lf\n",a[i*n + n -1],b[i]);
        }
       fprintf(fout,"..\n");
        for(int j = 0; j<restr;j++)
        {
            fprintf(fout,"%.2Lf ",a[(n-1)*n+j]);
        }
       fprintf(fout,".. %.2Lf    %.2Lf\n",a[(n-1)*n + n -1],b[n-1]);
       
    }

}

int solve(int n, long double* a, long double* b, long double* x)
{   
    long double sq, sphi, cphi, olx, oly;
    for(int j = 0;j<n-1;j++)
    {
        for(int i = j+1; i<n;i++)
        {
            if(!eq(a[A(i,j)],0))
            {
                long double ratio = a[A(j,j)]/a[A(i,j)];
                //sq = a[A(i,j)]*sqrt(1+ratio*ratio);
                sphi = - 1 / sqrt(1+ratio*ratio);
                if(a[A(i,j)]<0) sphi = -sphi;
                cphi = sqrt(1-sphi*sphi);
                if (a[A(j,j)]<0) cphi = -cphi;
            }
            else if(!eq(a[A(j,j)],0))
            {
                long double ratio = a[A(i,j)]/a[A(j,j)];
                //sq = a[A(j,j)]*sqrt(1+ratio*ratio);
                cphi = 1/sqrt(1+ratio*ratio);
                if (a[A(j,j)]<0) cphi = -cphi;
                sphi = -sqrt(1-cphi*cphi);
                if(a[A(i,j)]<0) sphi = -sphi;
            }
            else return 0;
            //printf("%Lf %Lf\n",a[A(i,j)],a[A(i,j)]);
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
