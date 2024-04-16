#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <chrono>
#define _USE_MATH_DEFINES
using namespace std;
std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

vector<vector<double>> DDx(vector<vector<double>> A, double* dx,int* gc,int* n,int order){
    vector<vector<double>> derivative(*n + 2 * *gc,vector<double>(*n + 2 * *gc));

		if (order == 2) {
			for (int j = 1; j < *n + 2 * *gc; j++){
				for (int i = 1; i < *n + 2 * *gc - 1; i++) {
				derivative[i][j] = (A[i][j+1] - A[i][j-1])/2.0/ *dx;
				} 
			}
		} else if (order == 4) {
		
		for (int j = 1; j < *n + 2 * *gc; j++){
			for (int i = 1; i < *n + 2 * *gc; i++) {

			derivative[i][j] = (-A[i][j+2] + 8*A[i][j+1] - 8*A[i][j-1] + A[i][j-2])/12.0/ *dx;
			}
		} 
			
		} else if (order == 6) {
		for (int j = 1; j < *n + 2 * *gc; j++){
			for (int i = 1; i < *n + 2 * *gc; i++) {

			derivative[i][j] = (A[i][j+3]/60 - 3*A[i][j+2]/20 + 3*A[i][j+1]/4 - 3*A[i][j-1]/4 + 3*A[i][j-2]/20 - A[i][j-3]/60)/ *dx;
			}
		} 
		}

		

	return derivative;
}

vector<vector<double>> DDy(vector<vector<double>> A, double* dy,int* gc,int* n,int order)
{
    vector<vector<double>> derivative(*n + 2 * *gc,vector<double>(*n + 2 * *gc));

		if (order == 2) {
			for (int i = 1; i < *n + 2 * *gc - 1; i++){
				for (int j = 1; j < *n + 2 * *gc - 1; j++) {
				derivative[i][j] = (A[i+1][j] - A[i-1][j])/2.0/ *dy;
				} 
			}
		} else if (order == 4) {
		
		for (int i = 1; i < *n + 2 * *gc; i++){
			for (int j = 1; j < *n + 2 * *gc; j++) {

			derivative[i][j] = (-A[i+2][j] + 8*A[i+1][j] - 8*A[i-1][j] + A[i-2][j])/12.0/ *dy;
			}
		} 
			
		} else if (order == 6) {
		for (int i = 1; i < *n + 2 * *gc; i++){
			for (int j = 1; j < *n + 2 * *gc; j++) {

			derivative[i][j] = (A[i+3][j]/60 - 3*A[i+2][j]/20 + 3*A[i+1][j]/4 - 3*A[i-1][j]/4 + 3*A[i-2][j]/20 - A[i-3][j]/60)/ *dy;
			}
		} 
		}

		

	return derivative;
}

void create(string file, vector<vector<double>> A,int n) 
{ 

	
    // file pointer 
    fstream fout; 

    // opens an existing csv file or creates a new file. 
    fout.open(file, ios::out | ios::app); 
  
    int i,j;
	// Read the input 
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
            if (j < (n - 1)) {
                fout << A[i][j] << ",";
            }
            else if (j == (n - 1)) {
                fout << A[i][j] << "\n";
            }
    }
}

vector<vector<double>> Poisson(vector<vector<double>> L,vector<vector<double>> B,double dx,double dy,int gc, int n, int order,vector<double> y) 
{
	double emax = 1e-6,epsilon_max = 1;
	vector<vector<double>> error(n + 2*gc,vector<double>(n + 2*gc));
    vector<vector<double>> L_new(n + 2*gc,vector<double>(n + 2*gc));
    struct location
	{
		int x,y;
	};
    location lL;
	double BC_L,BC_R,BC_T,BC_B;
	BC_L = 1.0;
	BC_R = 0.0;
	BC_T = 0.0;
	BC_B = 0.0;

	
	// for (int i = 0; i < n+gc; i++){
	// L[i][0]    = 2.0*BC_L;// - L[i][1];
	// L[i][n+gc] = 2.0*BC_R; //- L[i][n];
	// L[0][i]    = dy * BC_B + L[1][i];
	// L[n+gc][i] = dy * BC_T + L[n][i];
	// }
  while (emax < epsilon_max){
	
	for (int i = 0; i < n + gc; i++){
	L[i][0]    = y[i]*BC_L - L[i][1];
	L[i][n+gc] = 2.0*BC_R - L[i][n];
	// L[0][i]    = 2.0*BC_B;
	// L[n+gc][i] = 2.0*BC_T;
	L[0][i]    = dy * BC_B + L[1][i];
	L[n+gc][i] = dy * BC_T + L[n][i];
	}

	for (int i = 1; i < n+gc; i++) {
		for (int j = 1; j < n+gc; j++) {
			L_new[i][j] = (dy*dy*(L[i+1][j] + L[i-1][j]) + dx*dx*(L[i][j+1] + L[i][j-1]) - B[i][j]*dx*dx*dy*dy)/(2*dx*dx + 2*dy*dy);
		}
	}

	for (int i = gc; i < n + gc; i++){
		for (int j = gc; j < n + gc; j++){
			error[i][j] = abs(L_new[i][j] - L[i][j]);
		}
	}

	epsilon_max = 0;

   	for(int i = 1; i < n + gc; i++){
        for(int j = 1; j < n + gc; j++){
			if(epsilon_max < error[i][j])
            {
				epsilon_max = error[i][j];
                lL.x = i + 1;
                lL.y = j + 1;
            }
		}
	}
    for (int j = 1; j < n + gc; j++) {
		for (int i = 1; i < n + gc; i++) {
            L[i][j] = L_new[i][j];

		}
	}
	
}	
	return L;

}

vector<vector<double>> PoissonResidual(vector<vector<double>> u,vector<vector<double>> v,double rho,double dx,double dy,int gc, int n,int order)
{
	vector<vector<double>> dudx(n + 2*gc,vector<double>(n + 2*gc)),dvdx(n + 2*gc,vector<double>(n + 2*gc));
	vector<vector<double>> dudy(n + 2*gc,vector<double>(n + 2*gc)),dvdy(n + 2*gc,vector<double>(n + 2*gc));
	vector<vector<double>> B(n + 2*gc,vector<double>(n + 2*gc));
	dudx = DDx(u,&dx,&gc,&n,order);
	dvdx = DDx(v,&dx,&gc,&n,order);
	dudy = DDy(u,&dy,&gc,&n,order);
	dvdy = DDy(v,&dy,&gc,&n,order);
	for (int j = 1; j < n + 2 * gc; j++) {
		for (int i = 1; i < n + 2 * gc; i++) {	
			B[i][j] = -rho*(dudx[i][j]*dudx[i][j] + 2*dudy[i][j]*dvdx[i][j] + dvdy[i][j]*dvdy[i][j]);
		}
	}
	
	return B;
}

void createxy(string file, vector<double> x, vector<double> y,int n) 
{ 

	
    // file pointer 
    fstream fout; 

    // opens an existing csv file or creates a new file. 
    fout.open(file, ios::out | ios::app); 
  
    
  
    int i;
  
    // Read the input 
    for (i = 0; i < n; i++) { 
  
        fout << x[i] << ", "
            << y[i] << ", "
            <<"\n"; 
  
    } 
}

int main(){
	
	int nodes_x, nodes_y, order, gc, i,i_max,j,n,k;
	double x_low,x_high,cfl,dt,dx,nu,y_low,y_high,dy,time_check = 0,time = 0,rho;
	string type;
	ifstream input;
	// *** INPUTS ***
	input.open("Input.txt", ios::in);
	string line;

		while (1){
			input >> nodes_x;
			input >> nodes_y;
			input >> order;
			input >> x_low;
			input >> x_high;
			input >> y_low;
			input >> y_high;
			input >> nu;
            input >> rho;
			if (input.eof()){
				break;
			}
			
		}
	
	const int nx = nodes_x, ny = nodes_y;
	if (nx != ny) {
		if     (nx < ny){n = ny;}
		else if(ny < nx){n = nx;}	
	}else{ n = nx;}

	gc = order/2;		            // Find Number of GC

	// ***** Define Mesh *****
	dx = (x_high - x_low)/nx;
	dy = (y_high - y_low)/ny;

	vector<vector<double>> u(n + 2*gc,vector<double>(n + 2*gc));
	vector<vector<double>> v(n + 2*gc,vector<double>(n + 2*gc));
	vector<vector<double>> B(n + 2*gc,vector<double>(n + 2*gc)), P(n + 2*gc,vector<double>(n + 2*gc));
	vector<double> x(n + 2*gc),y(n + 2*gc);
	// ***** Initial Conditions *****	
		/***** Grid *****/
		for ( i = 0; i < n + 2*gc; i++) {
			for (j = 0; j < n + 2*gc; j++){
				x[i] = (x_low - gc*dx + i*dx + dx/2); 
				y[j] = (y_low - gc*dy + j*dy + dy/2);

			}
		}
		/***** Initial Velocity Field *****/
		for ( i = 0; i< n + 2*gc; i++) {	
			for (j = 0; j < n + 2*gc; j++){
				u[i][j] = 0.0;
				v[i][j] = 0.0;
				// u[i][j] = 4*pow(2.7182818,-(((x[i] - 3.14159)*(x[i] - 3.14159)) + ((y[j] - 3.14159)*(y[j] - 3.14159)))/2/0.5/0.5)/(2*3.14159*0.5*0.5);
				// v[i][j] = 4*pow(2.7182818,-(((x[i] - 3.14159)*(x[i] - 3.14159)) + ((y[j] - 3.14159)*(y[j] - 3.14159)))/2/0.5/0.5)/(2*3.14159*0.5*0.5);
				// u[i][j] = 2*sin(x[i]*x[i]*0.125 + y[j]*y[j]*0.125)+1;
				// v[i][j] = 4*sin(x[i])*cos(y[j]);
				// u[i][j] = x[j]*x[j];
				// v[i][j] = 0.0;
				
			}
		} 
		/***** Pressure Initialization and Boundary Conditions *****/

		std::filesystem::remove("xy.csv");
		std::filesystem::remove("init.csv");
		create("init.csv",u,n + 2*gc);
		createxy("xy.csv",x,y,n + 2*gc);
		string filename;
	// Solve

	B = PoissonResidual(u,v,rho,dx,dy,gc,n,order);
	P = Poisson(P,B,dx,dy,gc,n,order,y);		
    std::filesystem::remove("B.csv");
	std::filesystem::remove("Pressure.csv");
	create("B.csv",B,n + 2*gc);
    create("Pressure.csv",P,n + 2*gc);			
			
		std::filesystem::current_path("..");
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		double time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		cout << "Time Elapsed = " << time_elapsed << "[ms]" << endl;
		string file = "out.csv";
		fstream fout;
		fout.open(file, ios::out | ios::app);
		fout << time_elapsed << ",";
	
	return 0;

}