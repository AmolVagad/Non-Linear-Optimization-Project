#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include<vector>
#include<time.h>
#include<sys/time.h>
#include<ctime>

#include"globals.h"
using namespace std;
using namespace dlib;


///////////////////////////// Function declarations///////////////////////////////////

//Function to calculate the point at shortest distance 
void cal_closest_points(const column_vector &rt);
// Function to calculate gradient of cost function
column_vector CalculateGradient(const column_vector rt, column_vector gradient);
// Function that uses Adam optimizer for optimization
column_vector AdamOptimizer(const column_vector rt);
//Function to calculate norm
double cal_norm(column_vector r);
//Function to calculate elementwise division 
dlib::matrix<double> pointwise_divide(dlib::matrix<double> m, dlib::matrix<double> v);
// Function to compute rotation of the data point
dlib::matrix<double> PerformRotation(dlib::matrix<double> R,dlib::matrix<double> t, dlib::matrix<double> point)

///////////////////////////////////////////////////////////////////////////////




//////////////////Initializations/////////////////////////////////

point_cloud_data measurement_data;
point_cloud_data model_data;
point_cloud_data transformed_data;




/////////////// Main function////////////////////////////////////////////////////////////
int main()
{
	

	// Load model data from .csv file that contains 3D point cloud data	
	
	ifstream infile1;
  	infile1.open ("icp_model.csv");
	
	char* pEnd;
	string x,y,z;

	 while(!infile1.eof()){
		getline(infile1,x, ',');
		getline(infile1,y, ',');
		getline(infile1,z);
		//getline(infile,index);
		model_data.x_coord.push_back(strtod(x.c_str(),&pEnd));
		model_data.y_coord.push_back(strtod(y.c_str(),&pEnd));
		model_data.z_coord.push_back(strtod(z.c_str(),&pEnd));
		measurement_data.index.push_back(-1);
	}

	model_data.size = model_data.x_coord.size();

	model_data.x_coord.pop_back();
	model_data.y_coord.pop_back();
	model_data.z_coord.pop_back();
	model_data.size = model_data.size - 1;

	//Rotational function test
	double theta = 0.03;
	double point_x = 0.003;
	double point_y = 0.005;
	double point_z = 0.0;
	dlib::matrix<double> R(3,3);
	dlib::matrix<double> t(3,1);

	R = cos(theta), -sin(theta), 0,
	    sin(theta), cos(theta), 0,
	    0, 0, 1;

	t = point_x, point_y, point_z;
	


	// Generate mesasurement data by rotating the model data
	PerformTransformationToAllPoints(R, t, &model_data, &measurement_data,1);

	cout<<"Measurement data size "<<measurement_data.size<<endl;

	//Calculating the closest distance for a point 
	column_vector rt(4), rt_lower(4), rt_upper(4), gradient(4);

	rt = -theta, -cos(theta)*point_x - sin(theta)*point_y, sin(theta)*point_x - cos(theta)*point_y, point_z;
	cout<<"rt: "<<rt<<endl;

	rt = -0.0,0.0,0.0,0.0;
	rt_lower = -1.0, -1.0,-1.0,-1.0;
	rt_upper = 1.0, 1.0, 1.0, 1.0;
	gradient = 0,0,0,0;
	double final_error = 0;
	// time measurement variables 

	int cpu_starttime , cpu_endtime;
	for(int i = 0; i<20; i++)
	{
		cout<<"iteration #: "<<i<<endl;
		cpu_starttime = clock();
		cal_closest_points(rt);    // Function call to find the shortest distance
		cpu_endtime = clock();
		cout<<"The time taken for calculation = "<<((cpu_endtime - cpu_starttime)/CLOCKS_PER_SEC)<<endl;
		rt = AdamOptimizer(rt);   // Function call to optimizer function

		cout<<"Rt parameters "<<rt<<endl;
		
		
	}
	cout<<"Error after optimization "<<final_error<<endl;

	
	


	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////







void cal_closest_points(const column_vector &rt)
{

	
	dlib::matrix<double> R(3,3);
	dlib::matrix<double> t(3,1);

	R = cos(rt(0)), -sin(rt(0)), 0,
	    sin(rt(0)), cos(rt(0)), 0,
	    0, 0, 1;

	t = rt(1), rt(2), rt(3);

	float distance = 0.0;
	int best_index ;
	float closest_distance;

	transformed_data.x_coord.clear();
	transformed_data.y_coord.clear();
	transformed_data.z_coord.clear();
	transformed_data.index.clear();		
	PerformTransformationToAllPoints(R, t, &measurement_data, &transformed_data,1);

	
	int numPointsUpdated = 0;

	
	for(int i = 0; i < transformed_data.size; i++)
	{
		best_index = 0;
		closest_distance = 65535;

		for(int j = 0; j < model_data.size; j++)
		{
			
			distance = sqrt(pow((transformed_data.x_coord.at(i) - model_data.x_coord.at(j)),2) + pow((transformed_data.y_coord.at(i) - model_data.y_coord.at(j)),2) + pow((transformed_data.z_coord.at(i) - model_data.z_coord.at(j)),2));
			
			if(distance < closest_distance)
			{
				closest_distance = distance;
				best_index = j;
			}
		
		}
		measurement_data.index[i] = best_index;
		
	}
	

}


column_vector CalculateGradient(column_vector rt, column_vector gradient)
{
	int curr_index;
	double p1_sq, p2_sq, p3_sq;
	//cout<<"Transformed data size "<<transformed_data.size<<endl;
	for(int i = 0; i < transformed_data.size; i++)
	{
		//cout<<"Test i "<<i<<endl;
		curr_index = measurement_data.index[i]; 
		//cout<<"Test index "<<curr_index<<endl;
		p1_sq = pow(transformed_data.x_coord[i] - model_data.x_coord[curr_index],2);
		p2_sq = pow(transformed_data.y_coord[i] - model_data.y_coord[curr_index],2);
		p3_sq = pow(transformed_data.z_coord[i] - model_data.z_coord[curr_index],2);
		
		gradient(0) += 2*sqrt(p1_sq)*(-measurement_data.x_coord[i]*sin(rt(0)) - measurement_data.y_coord[i]*cos(rt(0))) + 2*sqrt(p2_sq)*(measurement_data.x_coord[i]*cos(rt(0)) - measurement_data.y_coord[i]*sin(rt(0)));
		
		gradient(1) += 2*sqrt(p1_sq);
		gradient(2) += 2*sqrt(p2_sq);
		gradient(3) += 2*sqrt(p3_sq);		
	}
	return gradient;
}

column_vector AdamOptimizer(column_vector rt)
{
	double alpha = pow(10,-7);
	double beta1 = 0.9;
	double beta2 = 0.999;
	double epsilon = pow(10,-8);
	double alpha_curr = 0;
	int t = 0;	
	dlib::matrix<double> m_old(4,1);
	m_old = 0, 0, 0, 0;
	dlib::matrix<double> v_old(4,1);
	v_old = 0,0,0,0;
	dlib::matrix<double> v(4,1);
	v = 0,0,0,0;
	column_vector gradient(4,1);
	gradient = 0,0,0,0;
	dlib::matrix<double> m(4,1);
	m = 0,0,0,0;
	column_vector rt_old(4);
	rt_old = rt;
	double error = 0.0;
	dlib::matrix<double> s(4,1);
	dlib::matrix<double> R(3,3);
	dlib::matrix<double> tr(3,1);
	while(t < 100)
	{	
		R = cos(rt(0)), -sin(rt(0)), 0,
	    	    sin(rt(0)), cos(rt(0)), 0,
	    	    0, 0, 1;
		tr = rt(1), rt(2), rt(3);
		
		transformed_data.x_coord.clear();
		transformed_data.y_coord.clear();
		transformed_data.z_coord.clear();
		transformed_data.index.clear();
		PerformTransformationToAllPoints(R, tr, &measurement_data, &transformed_data,1);
	
		t = t+1;
		gradient = {0.0,0.0,0.0,0.0};
		
		gradient = CalculateGradient(rt, gradient);
		m = beta1*m_old + (1-beta1)*gradient;
		v = beta2*v_old + (1-beta2)*pow(gradient,2);
		
		alpha_curr = alpha*sqrt(1-pow(beta2,t))/(1 - pow(beta1,t));
		s = pointwise_divide(m, sqrt(v) + epsilon);		
		rt = rt_old - alpha_curr*pointwise_multiply(m,s);
		//rt = rt_old - alpha*gradient;
		if(cal_norm(rt - rt_old) < epsilon)
		{
			cout<<"Converged at iter "<<t<<endl;
			break;
		}		
		m_old = m;
		v_old = v;
		rt_old = rt;
		
		error = findTotalErrorInCloud(rt);
		cout<<"Error value is "<<error<<endl;
		
	}
	return rt;	

}

dlib::matrix<double> pointwise_divide(dlib::matrix<double> m, dlib::matrix<double> v)
{
	dlib::matrix<double> output(m.size(),1);
	for(int i = 0; i < m.size(); i++)
	{
		output(i) = m(i)/v(i);
	}
	return output;
}
double cal_norm(column_vector r)
{
	double norm_val = 0.0;
	for(int i = 0; i < r.size(); i++)
		norm_val += pow(r(i),2);

	return sqrt(norm_val);
}
	


