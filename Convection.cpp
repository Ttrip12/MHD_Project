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

vector<vector<double>> DDx(vector<vector<double>> A, double* dx,int* gc,int* n,int order)

{
    vector<vector<double>> derivative(*n + 2 * *gc,vector<double>(*n + 2 * *gc));

		if (order == 2) {
			for (int j = 1; j < *n + 2 * *gc; j++){
				for (int i = 1; i < *n + 2 * *gc - 1; i++) {
				derivative[i][j] = (A[i+1][j] - A[i-1][j])/2.0/ *dx;
				} 
			}
		} else if (order == 4) {
		
		for (int j = 1; j < *n + 2 * *gc; j++){
			for (int i = 1; i < *n + 2 * *gc; i++) {

			derivative[i][j] = (-A[i+2][j] + 8*A[i+1][j] - 8*A[i-1][j] + A[i-2][j])/12.0/ *dx;
			}
		} 
			
		} else if (order == 6) {
		for (int j = 1; j < *n + 2 * *gc; j++){
			for (int i = 1; i < *n + 2 * *gc; i++) {

			derivative[i][j] = (A[i+3][j]/60 - 3*A[i+2][j]/20 + 3*A[i+1][j]/4 - 3*A[i-1][j]/4 + 3*A[i-2][j]/20 - A[i-3][j]/60)/ *dx;
			}
		} 
		}

		

	return derivative;
}

vector<vector<double>> DDy(vector<vector<double>> A, double* dy,int* gc,int* n,int order)
{
    vector<vector<double>> derivative(*n + 2 * *gc,vector<double>(*n + 2 * *gc));

		if (order == 2) {
			for (int i = 1; i < *n + 2 * *gc; i++){
				for (int j = 1; j < *n + 2 * *gc; j++) {
				derivative[i][j] = (A[i][j+1] - A[i][j-1])/2.0/ *dy;
				} 
			}
		} else if (order == 4) {
		
		for (int i = 1; i < *n + 2 * *gc; i++){
			for (int j = 1; j < *n + 2 * *gc; j++) {

			derivative[i][j] = (-A[i][j+2] + 8*A[i][j+1] - 8*A[i][j-1] + A[i][j-2])/12.0/ *dy;
			}
		} 
			
		} else if (order == 6) {
		for (int i = 1; i < *n + 2 * *gc; i++){
			for (int j = 1; j < *n + 2 * *gc; j++) {

			derivative[i][j] = (A[i][j+3]/60 - 3*A[i][j+2]/20 + 3*A[i][j+1]/4 - 3*A[i][j-1]/4 + 3*A[i][j-2]/20 - A[i][j-3]/60)/ *dy;
			}
		} 
		}

		

	return derivative;
}

vector<vector<double>> DDxDDx(vector<vector<double>> A, double* dx,int* gc,int* n,int order)
{
    vector<vector<double>> derivative(*n + 2 * *gc,vector<double>(*n + 2 * *gc));

		if (order == 2) {
			for (int j = 1; j < *n + 2 * *gc; j++) {
				for (int i = 1; i < *n + 2 * *gc - 1; i++) {

				derivative[i][j] = (A[i+1][j] -2*A[i][j] + A[i-1][j])/ *dx/ *dx;
				
				}
			} 
		} else if (order == 4) {
			for (int j = 2; j < *n + 2 * *gc; j++) {
				for (int i = 2; i < *n + 2 * *gc; i++) {
				
				derivative[i][j] = (-A[i+2][j]/12.0 + 4.0*A[i+1][j]/3.0 - 5.0*A[i][j]/2.0 + 4*A[i-1][j]/3 - A[i-2][j]/12.0)/ *dx/ *dx;
				
				}
			} 
			
		} else if (order == 6) {
			for (int j = *gc; j < *n + 2 * *gc; j++) {
				for (int i = *gc; i < *n + 2 * *gc; i++) {

				derivative[i][j] = (A[i+3][j]/90 - 3*A[i+2][j]/20 + 3*A[i+1][j]/2 - 49*A[i][j]/18 + 3*A[i-1][j]/2 - 3*A[i-2][j]/20 + A[i-3][j]/90)/ *dx/ *dx;

				}
			} 
		}

	return derivative;
}

vector<vector<double>> DDyDDy(vector<vector<double>> A, double* dy,int* gc,int* n,int order)
{
  	vector<vector<double>> derivative(*n + 2 * *gc,vector<double>(*n + 2 * *gc));

		if (order == 2) {
			for (int i = *gc; i < *n + 2 * *gc; i++) {
				for (int j = *gc; j < *n + 2 * *gc; j++) {

				derivative[i][j] = (A[i][j+1] -2*A[i][j] + A[i][j-1])/ *dy/ *dy;
				
				}
			} 
		} else if (order == 4) {
			for (int i = *gc; i < *n + 2 * *gc; i++) {
				for (int j = *gc; j < *n + 2 * *gc; j++) {
				
				derivative[i][j] = (-A[i][j+2]/12.0 + 4.0*A[i][j+1]/3.0 - 5.0*A[i][j]/2.0 + 4*A[i][j-1]/3 - A[i][j-2]/12.0)/ *dy/ *dy;
				
				}
			} 
			
		} else if (order == 6) {
			for (int i = *gc; i < *n + 2 * *gc; i++) {
				for (int j = *gc; j < *n + 2 * *gc; j++) {

				derivative[i][j] = (A[i][j+3]/90 - 3*A[i][j+2]/20 + 3*A[i][j+1]/2 - 49*A[i][j]/18 + 3*A[i][j-1]/2 - 3*A[i][j-2]/20 + A[i][j-3]/90)/ *dy/ *dy;

				}
			} 
		}

	return derivative;
}

vector<vector<double>> Poisson(vector<vector<double>> L,vector<vector<double>> B,double dx,double dy,int gc, int n, int order) 
{
	double emax = 1;
	vector<vector<double>> error(n + 2*gc,vector<double>(n + 2*gc));
	
	for (int j = 1; j < n + 2 * gc; j++) {
		for (int i = 1; i < n + 2 * gc - 1; i++) {

			L[i][j] = (dy*dy*(L[i+1][j] + L[i-1][j]) + dx*dx*(L[i][j+1] + L[i][j-1]) - B[i][j]*dx*dx*dy*dy)/(2*dx*dx + 2*dy*dy);
				
		}
	}

	for (int i = 0; i < n+2*gc; i++){
		L[0][i] = 1 - L[1][i];
		L[n+gc][i] = 1 - L[n][i];
		L[i][0] = L[i][1];
		L[i][n+gc] = L[i][n];
	}


	// for (int i = 0; i < n+2*gc; i++){
	// 	for (int j = 0; j < n+2*gc; j++){
	// 		error[i][j] = 
	// 	}
	// }

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
		for (int i = 1; i < n + 2 * gc - 1; i++) {	
			B[i][j] = dudx[i][j]*dudx[i][j] + 2*dudy[i][j]*dvdx[i][j] + dvdy[i][j]*dvdy[i][j];
		}
	}
	return B;
}

vector<vector<double>> fill_gc(vector<vector<double>> A, int gc,int n){
	
	int sz = n + 2*gc;
	for(int i = 0; i < gc; i++){
		for(int j = 1; j < n; j++){
			A[i][j] = A[n + i][j];
			A[sz - i - 1][j] = A[2*gc - i - 1][j];
		}
	}	 
 	for(int j = 0; j < gc; j++){
		for(int i = 1; i < n; i++){
			A[i][j] = A[i][n + j];
			A[i][sz - j - 1] = A[i][2*gc - j - 1];
		}
	}	
	return A;
}

double find_dt(vector<vector<double>> u,vector<vector<double>> v,double dx,double cfl,double nu) {
	double dt,v_max,u_max;
	double temp = 0;
	struct location
	{
		int x,y;
	};
	
	
	location lu,lv;
	// double u_max = (*max_element(u.begin(),u.end()));
	// double v_max = (*max_element(v.begin(),v.end()));
	    
		for(int i = 0; i < u[0].size(); i++)
        for(int j = 0; j< v[0].size(); j++)
            if(u_max < u[i][j])
            {
                u_max = u[i][j];
                lu.x = i + 1;
                lu.y = j + 1;
            }
	
		for(int i = 0; i < 102; i++)
        for(int j = 0; j< 102; j++)
            if(v_max < v[i][j])
            {
                v_max = v[i][j];
                lv.x = i + 1;
                lv.y = j + 1;
            }
	
	
	double dtvelo;
	if (v_max > u_max){
		dtvelo = dx*cfl/v_max;
	}else{
		dtvelo = dx*cfl/u_max;
	}
	double dtvisc = (dx*dx*cfl)/nu;
	
	if (dtvisc > dtvelo){
		dt = dtvelo;
	}else{ 
		dt = dtvisc;
	}

	return dt;
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

		//while (getline(input,line)) {
		while (1){
			// getline(input,line);
			// if(line.find("#") == 0) {
			// 	continue;
			// }
			input >> nodes_x;
			input >> nodes_y;
			input >> order;
			input >> x_low;
			input >> x_high;
			input >> y_low;
			input >> y_high;
			input >> cfl;
			input >> type;
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
	dt = 1.0/n*cfl;                 // Timestep
	i_max = n*80;                   // Total Number of Iterations

	// ***** Define Mesh *****
	dx = (x_high - x_low)/nx;
	dy = (y_high - y_low)/ny;

	vector<vector<double>> u(n + 2*gc,vector<double>(n + 2*gc)), unew(n + 2*gc,vector<double>(n + 2*gc)), dudx(n + 2*gc,vector<double>(n + 2*gc)), d2u_d2x(n + 2*gc,vector<double>(n + 2*gc)), dudy(n + 2*gc,vector<double>(n + 2*gc)), d2u_d2y(n + 2*gc,vector<double>(n + 2*gc));
	vector<vector<double>> v(n + 2*gc,vector<double>(n + 2*gc)), vnew(n + 2*gc,vector<double>(n + 2*gc)), dvdx(n + 2*gc,vector<double>(n + 2*gc)), d2v_d2x(n + 2*gc,vector<double>(n + 2*gc)), dvdy(n + 2*gc,vector<double>(n + 2*gc)), d2v_d2y(n + 2*gc,vector<double>(n + 2*gc));
	vector<vector<double>> B(n + 2*gc,vector<double>(n + 2*gc)), P(n + 2*gc,vector<double>(n + 2*gc)),dudx(n + 2*gc,vector<double>(n + 2*gc)), dp_dx(n + 2*gc,vector<double>(n + 2*gc)),dudx(n + 2*gc,vector<double>(n + 2*gc)), dp_dy(n + 2*gc,vector<double>(n + 2*gc));
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
		for ( i = 0; i < n + 2*gc; i++) {
			
			for (j = 0; j < n + 2*gc; j++){
				// u[i][j] = sin(x[i])*sin(y[j]);
				// v[i][j] = 0.1;
				// u[i][j] = 4*pow(2.7182818,-(((x[i] - 3.14159)*(x[i] - 3.14159)) + ((y[j] - 3.14159)*(y[j] - 3.14159)))/2/0.5/0.5)/(2*3.14159*0.5*0.5);
				// v[i][j] = 4*pow(2.7182818,-(((x[i] - 3.14159)*(x[i] - 3.14159)) + ((y[j] - 3.14159)*(y[j] - 3.14159)))/2/0.5/0.5)/(2*3.14159*0.5*0.5);
				u[i][j] = sin(x[i]*x[i]*0.125 + y[j]*y[j]*0.125)+1;
				v[i][j] = sin(x[i])*cos(y[j]);
				
			}
		} 
		
		/***** Pressure Initialization and Boundary Conditions *****/

		std::filesystem::remove("xy.csv");
		std::filesystem::remove("init.csv");
		create("init.csv",u,n + 2*gc);
		createxy("xy.csv",x,y,n + 2*gc);
		// Create sub directory for csv
		std::filesystem::remove_all("DataU");
		std::filesystem::create_directory("DataU");
		std::filesystem::remove_all("DataV");
		std::filesystem::create_directory("DataV");
		std::filesystem::remove_all("DataMagnitude");
		std::filesystem::create_directory("DataMagnitude");
		string filename;

	// Solve
		k = 0;
		double snapshot_time = 4;
		//int num = i_max/1000;
		vector<vector<double>> V_mag(n + 2*gc,vector<double>(n + 2*gc));
		for (int iter = 0; iter < i_max; iter++){
			
			dudx = DDx(u,&dx,&gc,&n,order);
			d2u_d2x = DDxDDx(u,&dx,&gc,&n,order);
			dudy = DDy(u,&dy,&gc,&n,order);
			d2u_d2y = DDyDDy(u,&dy,&gc,&n,order);
			dvdx = DDx(v,&dx,&gc,&n,order);
			d2v_d2x = DDxDDx(v,&dx,&gc,&n,order);
			dvdy = DDy(v,&dy,&gc,&n,order);
			d2v_d2y = DDyDDy(v,&dy,&gc,&n,order);
			dp_dx = DDx(u,&dx,&gc,&n,order);
			dp_dy = DDy(v,&dy,&gc,&n,order);
				if(type == "burger"){

					dt =  find_dt(u,v,dx,cfl,nu);
					
					for ( j = 0; j < n + 2*gc; j++){
						for ( i = 0; i < n + 2*gc; i++){

						unew[i][j] = u[i][j] + dt*(-u[i][j]*dudx[i][j] - v[i][j]*dudy[i][j] + nu*(d2u_d2x[i][j] + d2u_d2y[i][j]) + dp_dx[i][j]/rho);

						vnew[i][j] = v[i][j] + dt*(-u[i][j]*dvdx[i][j] - v[i][j]*dvdy[i][j] + nu*(d2v_d2x[i][j] + d2v_d2y[i][j]) + dp_dy[i][j]/rho);

						u[i][j] = unew[i][j];

						v[i][j] = vnew[i][j];
						
						}
					}
				}					
			
				u = fill_gc(u,gc,n);
				v = fill_gc(v,gc,n);
				B = PoissonResidual(u,v,rho,dx,dy,gc,n,order);
				P = Poisson(P,B,dx,dy,gc,n,order);

				
				for ( j = 0; j < n + 2*gc; j++){
						for ( i = 0; i < n + 2*gc; i++){

						V_mag[i][j] = sqrt(u[i][j]*u[i][j] + v[i][j]*v[i][j]);
						
						}
					}
			
					time_check = time_check + dt;
					if (time_check >= time){
						cout << time_check << endl;
						string s = to_string(time);
						s.erase(s.find_last_not_of('0') + 1, s.length()-1);
   						s.erase(s.find_last_not_of('.') + 1, s.length()-1);
						if (time < 10) {
							s.insert(0,1,'0');
						}
						std::filesystem::current_path("DataU");
						filename = "TimeU_";
						filename.append(s);
						filename.append("s.csv");
						create(filename,u,n + 2*gc);
						std::filesystem::current_path("..");
						std::filesystem::current_path("DataV");
						filename = "TimeV_";
						filename.append(s);
						filename.append("s.csv");
						create(filename,v,n + 2*gc);
						std::filesystem::current_path("..");
						std::filesystem::current_path("DataMagnitude");
						filename = "TimeMag_";
						filename.append(s);
						filename.append("s.csv");
						create(filename,V_mag,n + 2*gc);
						std::filesystem::current_path("..");

						k = k+1;
						time = time + snapshot_time;
					}

			
			// if (iter%num == 0){         //Taking snapshots of runs
				
			// 	string s = to_string(k);   // Organizing files names so matlab can read them in order
			// 	if (k < 10) {
			// 		s.insert(0,1,'0');
			// 		s.insert(0,1,'0');
			// 		s.insert(0,1,'0');
			// 	}else if (k < 100){
			// 		s.insert(0,1,'0');
			// 		s.insert(0,1,'0');
			// 	}
			// 	else if (k < 1000){
			// 		s.insert(0,1,'0');
			// 	}

			// 	filename = "Iter_";
			// 	filename.append(s);
			// 	filename.append(".csv");
			// 	create(filename,u,n + 2*gc);
				
			// 	k = k + 1;
			// }
			
		}

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