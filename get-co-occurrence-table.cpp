//usage: c++ get-co-occurrence-table.cpp -o co -fopenmp
// ./co yourdata.csv treatment start p_value> ct_result.tsv


#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <dirent.h>
#include <vector>
#include <sstream>
#include "ranker.h"
#include <cstdlib>
#include <iomanip>
#include <math.h>

#include <cmath>
#include <cstdio>
#include <algorithm>

#pragma GCC diagnostic ignored "-Wwrite-strings"
using namespace std;

int checkFile(ifstream &input)
// check if file is presence before open file
{
  if(input.fail()){                           //    Check open
    cerr << "Can't open file\n";
    exit(EXIT_FAILURE);
    //return 1;
  }else{return 0;}
}

int open_csv (string filenameDIR,vector <vector <string> > &data){
//file contents to data-vector 
  char const* fin = filenameDIR.c_str();
  ifstream input;
  string s;
  input.open(fin);
  checkFile(input);
  while(getline(input,s))
    {
      istringstream ss (s);
      //this add data into the column (second number)
      vector <string> record;
      while (ss)
	{
	  string s1;
	  if(!getline(ss,s1,',')) break;
	  if(s1!=""){
	    s1.erase(remove( s1.begin(), s1.end(), '\"' ),s1.end());
	    record.push_back(s1);
	  }
	}
      //this add data into the row (first number)
      data.push_back(record);
    }//while
  input.close();
  return 0;
}

void printMatrix(vector <vector <string> > &dad){
//print matrix
  for (int i=0;i<dad.size();i++){
    for (int j=0;j<dad[i].size();j++){
      cout<<dad[i][j]+" "<<flush;
    }
    cout<<endl<<flush;
  }
}

void printNumMatrix(vector <vector <double> > &dad){
//print matrix
  for (int i=0;i<dad.size();i++){
    for (int j=0;j<dad[i].size();j++){
      cout<<dad[i][j]<<" "<<flush;
    }
    cout<<endl<<flush;
  }
}


double sum(vector<double> a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {    
    s += a[i];
  }
  return s;
}

double mean(vector<double> a)
{
  return sum(a) / a.size();
}


double sqsum(vector<double> a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {
    s += pow(a[i], 2);
  }
  return s;
}

double stdev(vector<double> nums)
{
  double N = nums.size();
  return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

vector<double> operator-(vector<double> a, double b)
{
  vector<double> retvect;
  for (int i = 0; i < a.size(); i++)
  {
    retvect.push_back(a[i] - b);
  }
  return retvect;
}

vector<double> operator*(vector<double> a, vector<double> b)
{
  vector<double> retvect;
  for (int i = 0; i < a.size() ; i++)
  {
     retvect.push_back(a[i] * b[i]);
  }
  return retvect;
}

double pearsoncoeff(vector<double> X, vector<double> Y)
//calculate pearson coefficient
{
  return sum((X - mean(X))*(Y - mean(Y))) / (X.size()*stdev(X)* stdev(Y));
}

void nrerror( char error_text[]) {
// print error message
  fprintf( stderr, "Error: %s\n", error_text);
  exit(1);
}

float betai(float a, float b, float x)
// Returns the incomplete beta function Ix(a, b).
{
  float betacf(float a, float b, float x);
  float gammln(float xx);
  void nrerror(char error_text[]);
  float bt;
  if (x < 0.0 || x > 1.0) 
    nrerror("Bad x in routine betai");
  if (x == 0.0 || x == 1.0) 
    bt=0.0;
  else // Factors in front of the continued fraction.
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0)) // Use continued fraction directly.
    return bt*betacf(a,b,x)/a;
  else // Use continued fraction after making the sym
    return 1.0-bt*betacf(b,a,1.0-x)/b; // metry transformation.
}


#define MAXIT 1000
#define EPS 3.0e-7
#define FPMIN 1.0e-30
float betacf(float a, float b, float x)
// Used by betai: Evaluates continued fraction for incomplete beta function by modiﬁed Lentz’s method (§5.2).
{
  void nrerror(char error_text[]);
  int m,m2;
  float aa,c,d,del,h,qab,qam,qap;
  qab=a+b; // These q’s will be used in factors that occur
  qap=a+1.0; // in the coeﬃcients (6.4.6).
  qam=a-1.0;
  c=1.0; // First step of Lentz’s method.
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) 
    d=FPMIN;
  d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d; // One step (the even one) of the recurrence.
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d; // Next step of the recurrence (the odd one).
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break; // Are we done?
  }
  //if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
  return h;
}


float gammln(float xx)
// Returns the value ln[Γ(xx)] for xx > 0.
{
  // Internal arithmetic will be done in double precision, a nicety that you can omit if ﬁve-ﬁgure
  // accuracy is good enough.
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}


float pvalue( float t, float df ) 
// Compute p-value of t-statistic
{
  return betai(0.5*df,0.5,df/(df+t*t));
}


int main(int argc, char *argv[])
{

  int checkFile(ifstream &input);
  int open_csv (string filenameDIR,vector <vector <string> > &data);
  void printMatrix(vector <vector <string> > &dad);
  void printNumMatrix(vector <vector <double> > &dad);

  double sum(vector<double> a);
  double mean(vector<double> a);
  double sqsum(vector<double> a);
  double stdev(vector<double> nums);
  vector<double> operator-(vector<double> a, double b);
  vector<double> operator*(vector<double> a, vector<double> b);
  double pearsoncoeff(vector<double> X, vector<double> Y);

  int st = 0;
  int tr = 0;
  double pcut = 0.0;
  string filename = "";
  
  //Read argument
  if (argc < 5){
    cout << "Usage is ./co -f <infile> -t <treatment> -s <start> -p <pvalue>\n";
    exit(0);
  }else {
    for (int i=1;i < argc;i+=2){
      if (i + 1 != argc){
	if (string(argv[i]) == "-f") {
	    filename = argv[i + 1];
	} else if (string(argv[i]) == "-t"){
	  tr = atoi(argv[i + 1]);
	} else if (string(argv[i]) == "-s"){
	  st = atoi(argv[i + 1]);
	} else if (string(argv[i]) == "-p"){
	  pcut = atof(argv[i + 1]);
	} else {
	    std::cout << "Not enough or invalid arguments, please try again.\n";
	    exit(0);
	}
      }
    }
  }
  
  //Read data
  vector <vector <string> > data;
  open_csv(filename,data);

  vector <vector <string> > info;
  for (int i = 1; i< data.size();i++){
    vector <string> Dtemp;
    for (int j = 0;j<st-1;j++){
      Dtemp.push_back(data[i][j]);
    }
    info.push_back(Dtemp);
  }
  
  
  //get treatment
  vector <string> trt(info.size());
  for (int i = 0;i<info.size();i++){
    trt[i] = info[i][tr-1];
  }
  sort(trt.begin(), trt.end() );
  trt.erase (unique(trt.begin(), trt.end() ), trt.end());

  string tr1;
  for (int tri = 0;tri < trt.size();tri++){
    tr1 = trt[tri];
    vector <int> trnum;
    for (int i = 0;i<info.size();i++){
      if( info[i][tr-1] == tr1){
	trnum.push_back(i);
      }
    }

    vector <vector <double> > num;
    int row;
    for (int i = 0 ;i<trnum.size();i++){
      row = trnum[i]+1 ;
      vector <double> Dtemp;
      for (int j = st-1;j<data[row].size();j++){
		Dtemp.push_back(atof(data[row][j].c_str()));
      }
      num.push_back(Dtemp);
    }

    int nsize = num.size();
    string method = "average";
    
    vector <double> z(num[0].size());
    vector <vector <double> > rerho(num[0].size(),z);
    vector <vector <double> > rep(num[0].size(),z);
    //This is main for loop
    #pragma omp parallel
    #pragma omp for
    for (int nownum = 0;nownum < num[0].size()-1;nownum++){
       vector<double> a(nsize);
       for (int i = 0;i < nsize;i++){
         a[i] = num[i][nownum];
       }
       vector<double> aranks;
       rank(a, aranks, method);

       for (int nownumb = nownum+1;nownumb < num[0].size();nownumb++){
          vector<double> b(nsize);

	  for (int i = 0;i < nsize;i++){
            b[i] = num[i][nownumb];
          }
  
	  vector<double> branks;
	  rank(b, branks, method);
	  

	  int n = aranks.size();
	  double sumd = 0.0;
	  //get rho value
	  double rho = pearsoncoeff(aranks, branks);
	  rerho[nownum][nownumb] = rho;
	  
	  //get p value using t statistics
	  float t = rho * sqrt((n-2)/(1-rho*rho));
	  float df = n - 2;
	  float p = pvalue( t, df );
	  rep[nownum][nownumb] = p;

        }
     }
    
    //print output
    for (int i = 0;i< num[0].size()-1;i++){
      for (int j = i+1;j < num[0].size();j++){
	if (rep[i][j] < pcut){
	  cout <<tr1<<"\t"<<data[0][i+st-1]<<"\t"<<data[0][j+st-1]<<"\t"<< rerho[i][j] <<"\t"<<rep[i][j]<<endl;
	}
      }
    }
     
  }
  
}
