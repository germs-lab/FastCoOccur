//usage: c++ get-co-occurrence-table.cpp -o co -fopenmp
// ./co yourdata.csv treatment start> ct_result.tsv


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
using namespace std;

int main(int argc, char *argv[])
{

  int checkFile(ifstream &input);
  int open_csv (string filenameDIR,vector <vector <string> > &data);
  void printMatrix(vector <vector <string> > &dad);
  void printNumMatrix(vector <vector <double> > &dad);

  int open_csv_double (string filenameDIR,vector <vector <double> > &data);
  double calculate_p(vector <vector <double> > &oritable, int n, double t, int col);
  int find_near(vector <vector <double> > &oritable, int n, double t);

  vector <vector <double> > get_t_table();
  double get_p(int n, double t, vector <vector <double> > &oritable);

  vector <vector <double> > oritable = get_t_table();
  //printNumMatrix(oritable);

  int st = atoi (argv[3]);
  int tr = atoi (argv[2]);
  //open file
  string filename = argv[1];
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
  
	  for (uint i = 0; i < branks.size(); ++i){
             sumd = sumd + (aranks[i] - branks[i]) * (aranks[i] - branks[i]);
          }
	  double rho = 1 - (6*sumd)/(n*(n*n -1));

	  rerho[nownum][nownumb] = rho;
	  double p = get_p(n, rho, oritable);
	  rep[nownum][nownumb] = p;

        }
     }

    //output
    for (int i = 0;i< num[0].size()-1;i++){
      for (int j = i+1;j < num[0].size();j++){
	cout <<tr1<<"\t"<<data[0][i+st-1]<<"\t"<<data[0][j+st-1]<<"\t"<< rerho[i][j] <<"\t"<<rep[i][j]<<endl;
      }
    }
    
   


  }
    
}
// check openfile
int checkFile(ifstream &input)
{
  if(input.fail()){                           //    Check open
    cerr << "Can't open file\n";
    exit(EXIT_FAILURE);
    //return 1;
  }else{return 0;}
}

//file contents to data-vector
int open_csv (string filenameDIR,vector <vector <string> > &data){
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
	    record.push_back(s1);
	  }
	}
      //this add data into the row (first number)
      data.push_back(record);
    }//while
  input.close();
  return 0;
}

//print matrix
void printMatrix(vector <vector <string> > &dad){
  for (int i=0;i<dad.size();i++){
    for (int j=0;j<dad[i].size();j++){
      cout<<dad[i][j]+" "<<flush;
    }
    cout<<endl<<flush;
  }
}

//print matrix                                                                  
void printNumMatrix(vector <vector <double> > &dad){
  for (int i=0;i<dad.size();i++){
    for (int j=0;j<dad[i].size();j++){
      cout<<dad[i][j]<<" "<<flush;
    }
    cout<<endl<<flush;
  }
}


//file contents to data-vector
int open_csv_double (string filenameDIR,vector <vector <double> > &data){
  char const* fin = filenameDIR.c_str();
  ifstream input;
  string s;
  input.open(fin);
  checkFile(input);
  while(getline(input,s))
    {
      istringstream ss (s);
      //this add data into the column (second number)
      vector <double> record;
      while (ss)
	{
	  string s1;
	  if(!getline(ss,s1,',')) break;
	  if(s1!=""){
	    record.push_back(atof(s1.c_str()));
	  }
	}
      //this add data into the row (first number)
      data.push_back(record);
    }//while
  input.close();
  return 0;
}



// calculate p
double calculate_p(vector <vector <double> > &oritable, int n, double t, int col){

  double p = 0.0;
  double tempp;
  double t1 = oritable[n][col];
  double t2 = oritable[n][col+1];
  double p1 = oritable[0][col];
  double p2 = oritable[0][col+1];

  tempp = p2 - (p2-p1)*(t2-t)/(t2-t1);

  if (tempp <= 0){
    p = pow(10,tempp);
  }else{
    p =1-pow(10,-tempp);
  }

  return p;
}

// find near
int find_near(vector <vector <double> > &oritable, int n, double t){
  for (int i = 1; i< 11;i++){
    if( t > oritable[n][i] && t < oritable[n][i+1] ){
      return i;
    }
  } 
  return 0;
}


//get t table
vector <vector <double> > get_t_table(){ 
  vector <vector <double> > oritable;

  double p0[] = {0,-0.301029996,-0.397940009,-0.522878745,-0.698970004,-1,-1.301029996,-1.698970004,-2,-2.301029996,-2.698970004,-3};
  vector <double> v0(p0, p0+12);
  oritable.push_back(v0);

  double p1[] = {1,1,1.376,1.963,3.078,6.314,12.71,31.82,63.66,127.3,318.3,636.6};
  vector <double> v1(p1, p1+12);
  oritable.push_back(v1);

  double p2[] = {2,0.816,1.08,1.386,1.886,2.92,4.303,6.965,9.925,14.09,22.33,31.6};
  vector <double> v2(p2, p2+12);
  oritable.push_back(v2);

  double p3[] = {3,0.765,0.978,1.25,1.638,2.353,3.182,4.541,5.841,7.453,10.21,12.92};
  vector <double> v3(p3, p3+12);
  oritable.push_back(v3);

  double p4[] = {4,0.741,0.941,1.19,1.533,2.132,2.776,3.747,4.604,5.598,7.173,8.61};
  vector <double> v4(p4, p4+12);
  oritable.push_back(v4);

  double p5[] = {5,0.727,0.92,1.156,1.476,2.015,2.571,3.365,4.032,4.773,5.893,6.869};
  vector <double> v5(p5, p5+12);
  oritable.push_back(v5);

  double p6[] = {6,0.718,0.906,1.134,1.44,1.943,2.447,3.143,3.707,4.317,5.208,5.959};
  vector <double> v6(p6, p6+12);
  oritable.push_back(v6);

  double p7[] = {7,0.711,0.896,1.119,1.415,1.895,2.365,2.998,3.499,4.029,4.785,5.408};
  vector <double> v7(p7, p7+12);
  oritable.push_back(v7);

  double p8[] = {8,0.706,0.889,1.108,1.397,1.86,2.306,2.896,3.355,3.833,4.501,5.041};
  vector <double> v8(p8, p8+12);
  oritable.push_back(v8);

  double p9[] = {9,0.703,0.883,1.1,1.383,1.833,2.262,2.821,3.25,3.69,4.297,4.781};
  vector <double> v9(p9, p9+12);
  oritable.push_back(v9);

  double p10[] = {10,0.7,0.879,1.093,1.372,1.812,2.228,2.764,3.169,3.581,4.144,4.587};
  vector <double> v10(p10, p10+12);
  oritable.push_back(v10);

  double p11[] = {11,0.697,0.876,1.088,1.363,1.796,2.201,2.718,3.106,3.497,4.025,4.437};
  vector <double> v11(p11, p11+12);
  oritable.push_back(v11);

  double p12[] = {12,0.695,0.873,1.083,1.356,1.782,2.179,2.681,3.055,3.428,3.93,4.318};
  vector <double> v12(p12, p12+12);
  oritable.push_back(v12);

  double p13[] = {13,0.694,0.87,1.079,1.35,1.771,2.16,2.65,3.012,3.372,3.852,4.221};
  vector <double> v13(p13, p13+12);
  oritable.push_back(v13);

  double p14[] = {14,0.692,0.868,1.076,1.345,1.761,2.145,2.624,2.977,3.326,3.787,4.14};
  vector <double> v14(p14, p14+12);
  oritable.push_back(v14);

  double p15[] = {15,0.691,0.866,1.074,1.341,1.753,2.131,2.602,2.947,3.286,3.733,4.073};
  vector <double> v15(p15, p15+12);
  oritable.push_back(v15);

  double p16[] = {16,0.69,0.865,1.071,1.337,1.746,2.12,2.583,2.921,3.252,3.686,4.015};
  vector <double> v16(p16, p16+12);
  oritable.push_back(v16);

  double p17[] = {17,0.689,0.863,1.069,1.333,1.74,2.11,2.567,2.898,3.222,3.646,3.965};
  vector <double> v17(p17, p17+12);
  oritable.push_back(v17);

  double p18[] = {18,0.688,0.862,1.067,1.33,1.734,2.101,2.552,2.878,3.197,3.61,3.922};
  vector <double> v18(p18, p18+12);
  oritable.push_back(v18);

  double p19[] = {19,0.688,0.861,1.066,1.328,1.729,2.093,2.539,2.861,3.174,3.579,3.883};
  vector <double> v19(p19, p19+12);
  oritable.push_back(v19);

  double p20[] = {20,0.687,0.86,1.064,1.325,1.725,2.086,2.528,2.845,3.153,3.552,3.85};
  vector <double> v20(p20, p20+12);
  oritable.push_back(v20);

  double p21[] = {21,0.686,0.859,1.063,1.323,1.721,2.08,2.518,2.831,3.135,3.527,3.819};
  vector <double> v21(p21, p21+12);
  oritable.push_back(v21);

  double p22[] = {22,0.686,0.858,1.061,1.321,1.717,2.074,2.508,2.819,3.119,3.505,3.792};
  vector <double> v22(p22, p22+12);
  oritable.push_back(v22);

  double p23[] = {23,0.685,0.858,1.06,1.319,1.714,2.069,2.5,2.807,3.104,3.485,3.767};
  vector <double> v23(p23, p23+12);
  oritable.push_back(v23);

  double p24[] = {24,0.685,0.857,1.059,1.318,1.711,2.064,2.492,2.797,3.091,3.467,3.745};
  vector <double> v24(p24, p24+12);
  oritable.push_back(v24);

  double p25[] = {25,0.684,0.856,1.058,1.316,1.708,2.06,2.485,2.787,3.078,3.45,3.725};
  vector <double> v25(p25, p25+12);
  oritable.push_back(v25);

  double p26[] = {26,0.684,0.856,1.058,1.315,1.706,2.056,2.479,2.779,3.067,3.435,3.707};
  vector <double> v26(p26, p26+12);
  oritable.push_back(v26);

  double p27[] = {27,0.684,0.855,1.057,1.314,1.703,2.052,2.473,2.771,3.057,3.421,3.69};
  vector <double> v27(p27, p27+12);
  oritable.push_back(v27);

  double p28[] = {28,0.683,0.855,1.056,1.313,1.701,2.048,2.467,2.763,3.047,3.408,3.674};
  vector <double> v28(p28, p28+12);
  oritable.push_back(v28);

  double p29[] = {29,0.683,0.854,1.055,1.311,1.699,2.045,2.462,2.756,3.038,3.396,3.659};
  vector <double> v29(p29, p29+12);
  oritable.push_back(v29);

  double p30[] = {30,0.683,0.854,1.055,1.31,1.697,2.042,2.457,2.75,3.03,3.385,3.646};
  vector <double> v30(p30, p30+12);
  oritable.push_back(v30);

  double p40[] = {40,0.681,0.851,1.05,1.303,1.684,2.021,2.423,2.704,2.971,3.307,3.551};
  vector <double> v40(p40, p40+12);
  oritable.push_back(v40);

  double p50[] = {50,0.679,0.849,1.047,1.299,1.676,2.009,2.403,2.678,2.937,3.261,3.496};
  vector <double> v50(p50, p50+12);
  oritable.push_back(v50);

  double p60[] = {60,0.679,0.848,1.045,1.296,1.671,2,2.39,2.66,2.915,3.232,3.46};
  vector <double> v60(p60, p60+12);
  oritable.push_back(v60);

  double p80[] = {80,0.678,0.846,1.043,1.292,1.664,1.99,2.374,2.639,2.887,3.195,3.416};
  vector <double> v80(p80, p80+12);
  oritable.push_back(v80);

  double p100[] = {100,0.677,0.845,1.042,1.29,1.66,1.984,2.364,2.626,2.871,3.174,3.39};
  vector <double> v100(p100, p100+12);
  oritable.push_back(v100);

  double p120[] = {120,0.677,0.845,1.041,1.289,1.658,1.98,2.358,2.617,2.86,3.16,3.373};
  vector <double> v120(p120, p120+12);
  oritable.push_back(v120);

  double p1000[] = {1000,0.674,0.842,1.036,1.282,1.645,1.96,2.326,2.576,2.807,3.09,3.291};
  vector <double> v1000(p1000, p1000+12);
  oritable.push_back(v1000);

  return oritable;
}


double get_p(int n, double r,vector <vector <double> > &oritable){
  //get t
  double t = r * sqrt((n-2)/(1-(r*r)));
  
  if (n >= 1 && n <= 30) {
    n = n;
  }else if( n > 30 && n <=40){
    n = 31;
  }else if( n >40 && n<=50){
    n = 32;
  }else if( n >50 && n<=60){
    n = 33;
  }else if( n >60 && n<=80){
    n = 34;
  }else if( n >80 && n<=100){
    n = 35;
  }else if( n >100 && n<=120){
    n = 36;
  }else if( n >120 ){
    n = 37;
  }

  double p;
  int col;
  if (t > oritable[n][1] && t < oritable[n][11]){
    col = find_near(oritable, n, t);
    p = calculate_p(oritable, n, t, col);
  }else if(t <=  oritable[n][1]){
    col = 1;
    p = calculate_p(oritable, n, t, col);
  }else if(t >=  oritable[n][11]){
    col = 10;
    p = calculate_p(oritable, n, t, col);
  }
  return p;
}
