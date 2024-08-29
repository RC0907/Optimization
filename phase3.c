#include <stdio.h>
#include <math.h>


	
double func(double x[],int n,double R) 
{
    
    double temp1=0;
    double temp2 = 0;
    //temp1=pow((pow(x[0],2)+x[1]-11),2)+pow((x[0]+pow(x[1],2)-7),2);
   // temp2=((pow((x[0]-5),2))+(pow(x[1],2))-26);
     /* double temp3=0;
     temp1=pow((x[0]-10),3)+pow((x[1]-20),3);
      temp2=pow((x[0]-5),2)+pow((x[1]-5),2)-100;
      temp3=pow((x[0]-6),2)+pow((x[1]-5),2)-82.81;
     return temp1+R*pow(temp2,2)+R*pow(temp3,2);
    */
    
    double temp3=0;
    /*
    double temp4=0;
    double temp5=0;
    double temp6=0;
    double temp7=0;
   
    temp1=x[0]+x[1]+x[2];
    temp2=-1+.0025*(x[3]+x[5]);
    temp3=-1+.0025*(-x[3]+x[4]+x[6])+100*x[0]-x[0]*x[5];
    temp4=-1+.01*(-x[5]+x[7]);
    temp5=x[1]*x[3]-x[1]*x[6]-1250*x[3]+1250*x[4];
    temp6=100*x[0]-x[0]*x[5]+833.33252*x[3]-83333.333;
    temp7=x[2]*x[4]-x[2]*x[6]-2500*x[4]+1250000;
    return temp1+R*pow(temp2,2)+R*pow(temp3,2)+R*pow(temp4,2)+R*pow(temp5,2)+R*pow(temp6,2)+R*pow(temp7,2);
    */
    temp1=-(pow(sin(2*3.14*x[0]),3)*sin(2*3.14*x[1])/(pow(x[0],3)*(x[0]+x[1])));
    temp2=pow(x[0],2)-x[1]+1;
    temp3=1-x[0]+pow((x[1]-4),2);
    return temp1+R*pow(temp2,2)+R*pow(temp3,2);
    
  // return temp1 + R*pow(temp2,2); 
}
double funct(double x[],int n){
double temp2=0;
             //temp2=pow((x[0]-5),2)+pow((x[1]-5),2)-100+pow((x[0]-6),2)+pow((x[1]-5),2)-82.81;  
		//temp2=((pow((x[0]-5),2))+(pow(x[1],2))-26);
		//return temp2;
		
	double temp3=0;
	/*
	double temp4=0;
	double temp5=0;
	double temp6=0;
	double temp7=0;
	 temp2=-1+.0025*(x[3]+x[5]);
    temp3=-1+.0025*(-x[3]+x[4]+x[6])+100*x[0]-x[0]*x[5];
    temp4=-1+.01*(-x[5]+x[7]);
    temp5=x[1]*x[3]-x[1]*x[6]-1250*x[3]+1250*x[4];
    temp6=100*x[0]-x[0]*x[5]+833.33252*x[3]-83333.333;
    temp7=x[2]*x[4]-x[2]*x[6]-2500*x[4]+1250000;
    return temp2+temp3+temp4+temp5+temp6+temp7;
    */
    
      temp2=pow(x[0],2)-x[1]+1;
    temp3=1-x[0]+pow((x[1]-4),2);
    return temp2+temp3;
    /* temp2=pow((x[0]-5),2)+pow(x[1]-5),2)-100;
      temp3=pow((x[0]-6),2)+pow((x[1]-5),2)-82.81;
      return temp2+temp3;
      */
}

double objective_function(double c,int n,double x[],double product[][1],double R) 
{
       double D[n];
    for (int i = 0; i < n; i++) 
    {
        D[i]=x[i];
        x[i] = x[i] - c* product[i][0];
       // printf("value of x taken %.20lf\n",x[i]);
    }
   
  
    double temp1 = 0.0;
    double temp2 = 0.0;
    
   /*
     double temp3=0;
     temp1=pow((x[0]-10),3)+pow((x[1]-20),3);
      temp2=pow((x[0]-5),2)+pow((x[1]-5),2)-100;
      temp3=pow((x[0]-6),2)+pow((x[1]-5),2)-82.81;
   
     */
    
    double temp3=0;
    /*
    double temp4=0;
    double temp5=0;
    double temp6=0;
    double temp7=0;
   
    temp1=x[0]+x[1]+x[2];
    temp2=-1+.0025*(x[3]+x[5]);
    temp3=-1+.0025*(-x[3]+x[4]+x[6])+100*x[0]-x[0]*x[5];
    temp4=-1+.01*(-x[5]+x[7]);
    temp5=x[1]*x[3]-x[1]*x[6]-1250*x[3]+1250*x[4];
    temp6=100*x[0]-x[0]*x[5]+833.33252*x[3]-83333.333;
    temp7=x[2]*x[4]-x[2]*x[6]-2500*x[4]+1250000;
   */
    
     //temp1=pow((pow(x[0],2)+x[1]-11),2)+pow((x[0]+pow(x[1],2)-7),2);
   // temp2=((pow((x[0]-5),2))+(pow(x[1],2))-26);
    
     temp1=-(pow(sin(2*3.14*x[0]),3)*sin(2*3.14*x[1])/(pow(x[0],3)*(x[0]+x[1])));
    temp2=pow(x[0],2)-x[1]+1;
    temp3=1-x[0]+pow((x[1]-4),2);
  
    
    for(int i=0;i<n;i++)
    {
        x[i]=D[i];
    }
  //  return temp1 + R*pow(temp2,2);
    //return temp1+R*pow(temp2,2)+R*pow(temp3,2);
 //  return temp1+R*pow(temp2,2)+R*pow(temp3,2)+R*pow(temp4,2)+R*pow(temp5,2)+R*pow(temp6,2)+R*pow(temp7,2);
     return temp1+R*pow(temp2,2)+R*pow(temp3,2);
}

double calculateDerivative(double c,int n, double x[],double product[][1],double R) 
{
    double f_x = objective_function(c, n, x, product,R);
    double h=0.0001;
    double f_x_plus_h = objective_function(c+ h, n, x, product,R);
    
    // Use the central difference formula for the derivative approximation
    double derivative = (f_x_plus_h - f_x) / h;
    
    return derivative;
}
double secant_Method (double a, double b, double e,int n,double x[],double product[][1],double R)
{
  double x1, x2, z, f3,c1,c2;
  int feval = 4;

  x1 = a;
  x2 = b;
  do
    {
      //z = (x1 + x2) / 2;
      z=x2-calculateDerivative (x2, n, x, product,R)*(x2-x1)/(calculateDerivative (x2, n, x, product,R)-calculateDerivative (x1, n, x, product,R));
      f3=calculateDerivative (z, n, x, product,R);
      
     
       if (f3 >0)
	{
	  x1 = z;
	}
      else  
	{
	  x2 = z;
	}
      feval++;
      c1=fabs (f3);
      c2=pow(10,-3);
    }
  while (c1>=c2);
 /*printf(" 167 %lf\n",fabs(f3));
  printf("168 %lf\n",pow(10,-3));
  
  printf (" 170 The minimum point lies between (%lf,%lf) \n", x1, x2);
  printf ("171 \nfunction evaluations(secant method): %d\n",feval);
  */
 //double ans= (x1+x2)/2;
  double ans=z;
  return ans;
}

double solution_Value(double c, double d, double k,double e,int n,double x[],double product[][1],double R) 
{
    double x_k = objective_function(c, n, x, product,R);
    double x_v = c-(pow(2, k-1)) * d;

    // Define a base case here to stop recursion
    
    c = c + (pow(2, k)) * d;

    double x_kplus_1 = objective_function(c, n, x, product,R);
    if (x_kplus_1 < x_k) 
    {
        k = k + 1;
        // Return the result of the recursive call
        return solution_Value(c, d, k,e, n, x, product,R);
    } 
    else 
    {
        printf(" \nInitial x: %.20lf\nFinal x: %.20lf\n", x_v, c);
        // Return the final value
         return secant_Method(x_v,c,e, n, x, product,R);
    }
}

double bounding_Phase_Method(double c, double d,double e, int n, double x[],double product[][1],double R) 
{
    double x_plus_Delta, x_neg_Delta, x_Value, delta;
    int fevel=3;
    d = fabs(d);
    x_plus_Delta = objective_function(c + d, n, x, product,R);
    x_neg_Delta = objective_function(c - d, n, x, product,R);
    x_Value = objective_function(c, n, x, product,R);
       /* printf("f+d= %.6lf\n",x_plus_Delta);
        printf("f-d= %.6lf\n",x_neg_Delta);
        printf("f= %6lf\n",x_Value);*/
        for(int i=0;i<n;i++)
        {
            printf(" the new value of x :%.5lf\n",x[i]);
            fevel++;
        }
        printf(" function evaluation(bounding phase method): %d",fevel);
    if (x_Value >= x_plus_Delta && x_Value <= x_neg_Delta) 
    {
        d = d;
        return solution_Value(c, d, 0,e, n, x, product,R);
    } 
    else if (x_Value <= x_plus_Delta && x_Value >= x_neg_Delta) 
    {
        d = -d;
        return solution_Value(c, d, 0, e, n, x, product,R);
    } 
    else 
    {
        printf("  Enter a new  Value.\n");
        return 0.0; // Modify this part to handle the case properly
    }
}

 
double guess_Value(double e,int n,double x[], double product[][1],int count,double R)
{
    double c;
    double d;
    if(count<1){
    
    printf("Enter initial guess c: ");
    scanf("%lf", &c);
    printf("Enter d: ");
    scanf("%lf", &d);
  
    bounding_Phase_Method(c, d, e, n, x, product,R);
  
   }
   else{
   bounding_Phase_Method(c,d,e,n,x,product,R);
   }
}
 double check_condition(int n,double product[][1],double arr1[][1],double R)
 {
     double sum=0;
     for(int i=0;i<n;i++)
     {
         //printf("%.2lf",arr1[i]);
         //printf("%.2lf",product[i]);
         
         //sum+=(arr1[i][0]*product[i][0]);
         //printf("%.20lf",sum);
     }
     return sum;
}
void product_matrix(int n, double identity[n][n], double arr1[n][1], double product[n][1],double R) 
{
    for (int i = 0; i < n; i++) 
    {
        double sum=0;
            for (int k = 0; k < n; k++) 
            {
                sum+=identity[i][k]* arr1[k][0];
            
            }
        product[i][0]=sum;
    }
    for(int i=0;i<n;i++)
    {
      //  printf("product matrix %.20lf",product[i][0]);
    }
}

void inverse(int n, double arr2[][n], double arr1[][1], double product[][1], double identity[][n],double R) 
{
    double pivot, factor;
    // Initialize the identity matrix
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            if (i == j)
                identity[i][j] = 1.0;
            else
                identity[i][j] = 0.0;
        }
    }

    // Perform Gaussian elimination with partial pivoting
    for (int i = 0; i < n; i++) 
    {
        pivot = arr2[i][i];

        // Divide both the current row of 'mat' and 'identity' by the pivot
        for (int j = 0; j < n; j++) 
        {
            arr2[i][j] /= pivot;
            identity[i][j] /= pivot;
        }

        for (int k = 0; k < n; k++) 
        {
            if (k != i) 
            {
                factor = arr2[k][i];

                for (int j = 0; j < n; j++) 
                {
                    arr2[k][j] -= factor * arr2[i][j];
                    identity[k][j] -= factor * identity[i][j];
                }
            }
        }
    }


   /* for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            printf("%.20lf",identity[i][j]);
        }
    }*/
}

double diff_second(double x[], int n, int i, int j,double R) 
{
    double h = 1e-5;
    double temp1 = x[i];
    double temp2 = x[j];

    x[i] += h;
    x[j] += h;

    double f1 = func(x,n,R);

    x[i] = temp1;
    x[j] = temp2;

    x[i] += h;
    double f2 = func(x,n,R);

    x[i] = temp1;

    x[j] += h;
    double f3 = func(x,n,R);

    x[j] = temp2;

    double f0 = func(x,n,R);

    return (f1 - f2 - f3 + f0) / (h * h);
}

void diff_first(int n,double x[], double arr1[n][1],double R) 
{
    // Define a small value for numerical differentiation
    double h = 1e-6;

    // Calculate partial derivatives using numerical differentiation
    for (int i = 0; i < n; i++) 
    {
        double temp1 = x[i];
        double temp = func(x,n,R);
        x[i] += h;
        arr1[i][0] = (func(x,n,R) - temp) / h;
        x[i] = temp1;
    }
}

int main() 
{
    int n;
    int m;
    int count=0;
    double R;
    double c1;
    printf("NO OF ITRATION :");
    scanf("%d", &m);
    printf("NO OF VARIABLE :");
    scanf("%d", &n);
    printf("VAlUE OF R :");
    scanf("%lf",&R);
    printf(" VALUE OF C1 :");
    scanf("%lf",&c1);
    double arr1[n][1];
    double arr2[n][n];
    double product[n][1];
    double identity[n][n];
    double x[n];
   
   
    double e = 1e-5;
   
    double ans5=100;
    
    double e3 = pow(10,-3);
    double c;
   FILE*fp;
fp = fopen("output.txt","w");
    int counting =0;
    printf("guess the intial value of all variable : ");
    for (int i = 0; i < n; i++) 
    {
        scanf("%lf", &x[i]);
    }

    /* diff_first(n,x,arr1);
       // int ans2=0;
       printf("first derivative\n");
        for(int i=0;i<n;i++){
            printf("%.20lf\n",arr1[i][0]);
        }
        //ans1=sqrt(ans2);
        printf("second derivative\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                arr2[i][j] = diff_second(x, n, i, j);
                printf("%.20lf\t",arr2[i][j]);
            }
            printf("\n");
        }
   printf("inverse\n");
        inverse(n, arr2, arr1, product, identity);
         for(int i=0;i<n;i++){
           for(int j=0;j<n;j++){
               printf("%.20lf\t",identity[i][j]);
           }
           printf("\n");
       }
        product_matrix(n, identity, arr1, product);
         for(int i=0;i<n;i++){
           printf("%.20lf\n",product[i][0]);
       }
       for(int i=0;i<n;i++){
           for(int j=0;j<n;j++){
               printf("%.20lf\t",identity[i][j]);
           }
           printf("\n");
       }*/
        double x2[n];
         double e1 = pow(10,-6);
    double e2 = pow(10,-6);
         double ans1=100;
         double ans2=100;
          int ittration=0;
        double x1[n];
        double R2=0;
   while(ans5>e3){
              
     
while (ans1 > e1 && ittration<=m && ans2 > e2) 
{
        diff_first(n,x, arr1,R);
        double ans3=0;
        ans1=100;
        for(int i=0;i<n;i++)
        {
        ans3+=(pow(arr1[i][0], 2));
        }
        ans1=sqrt(ans3);
       printf("ans 1 value %.6lf",ans1);
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                arr2[i][j] = diff_second(x, n, i, j,R);
            }
        }

        inverse(n, arr2, arr1, product, identity,R);
        product_matrix(n, identity, arr1, product,R);
         
       double  condition_check=check_condition(n,product,arr1,R);
     // printf("%.20lf",condition_check);
      if(condition_check<0)
      {
         printf("Please inter new guess value");
        break;
      }
       for(int i=0;i<n;i++)
       {
           x1[i]=x[i];
           printf("\n%.4lf\n",x1[i]); //output come after funcValue
       }
        double alpha = guess_Value(e, n, x, product,count,R);
       // printf("%.4lf\n",alpha);
       count++;
        for(int i=0;i<n;i++) //output come before funcValue
        {
            printf("%.4lf\n",x[i]);
            x[i]=x1[i]-alpha*product[i][0];
            printf("%.4lf\n",x[i]);
        }
        double ans_value=0;
        double ans1_value=0;
        ans2=100;
        for(int i=0;i<n;i++)
        {
          ans_value+=pow((x[i]-x1[i]),2);
        //  printf("ans_value %.2lf\n",ans_value);
          ans1_value+=pow((x1[i]),2);
        //  printf("ans1_value %.2lf\n",ans1_value);
        }
         printf("Value of voilation : %.8lf\n",funct(x,n));
         fprintf(fp,"Value of voilation : %.8lf\n",funct(x,n));
        ans_value=sqrt(ans_value);
        ans1_value=sqrt(ans1_value);
        double funcValue=func(x,n,R);
        printf("funcValue : %.20lf\n",funcValue);
        fprintf(fp,"funcValue : %.20lf\n",funcValue);
       //  printf("ans_value %.4lf\n",ans_value);
       // printf("ans1_value %.4lf\n",ans1_value);
        ans2=ans_value/ans1_value;
       printf("ans2 value %.4lf\n",ans2);
        ittration++;
         
    }
    for(int i=0;i<n;i++){
    printf("final value of x after ittration come : %.9lf\n",x[i]);
    fprintf(fp,"final value of x after ittration come : %.9lf\n",x[i]);
    }
 
    if(counting>0){
    double sum7=(func(x2,n,R2));
    printf("the valu eof sum7 %.20lf\n",sum7);
    fprintf(fp,"the valu eof sum7 %.20lf\n",sum7);
    double sum8=(func(x,n,R));
      printf("the valu eof sum8 %.20lf\n",sum8);
      fprintf(fp,"the valu eof sum8 %.20lf\n",sum8);
    ans5=fabs(sum8-sum7);
    }
    for(int i=0;i<n;i++){
    x2[i]=x[i];
    }
    printf("the valu eof the ans5 : %.20lf",ans5);
    fprintf(fp,"the valu eof the ans5 : %.20lf",ans5);
    counting++;
    R2=R;
     R=R*c1;
     printf("Value of r : %.20lf",R);
     fprintf(fp,"Value of r : %.20lf",R);
    
     
          ans1=100;
          ans2=100;
          ittration=0;
    }
   
    return 0;
}
