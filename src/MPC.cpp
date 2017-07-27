/**********************************************
* Self-Driving Car Nano-degree - Udacity
*  Created on: July 12, 2017
*      Author: Munir Jojo-Verge
**********************************************/
#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;
using namespace std;
using namespace Eigen;

// TODO: Set the timestep length and duration
const uint32_t N = 15; // prediction Horizon
const double dt = 0.08;
const double T = N*dt; // This is the Prediction Horizon in seconds. 


// Defining our references (setpoints). We want the cte = 0 the yaw error = 0 and...
const double ref_cte  = 0;
const double ref_epsi = 0;
const double ref_acc = 1;
AD<double> ref_v; // The Reference Velocity Magnitude will change depending on where we are on the track


// The solver takes all the state variables and actuator
// variables in a single vector. On the lectured we've seen how this structure is.
// Here, we establish when one variable starts and another ends to be able to address its indexes in an easy way.
uint32_t x_start = 0;
uint32_t y_start = x_start + N;
uint32_t psi_start = y_start + N;
uint32_t v_start = psi_start + N;
uint32_t cte_start = v_start + N;
uint32_t epsi_start = cte_start + N;
uint32_t delta_start = epsi_start + N;
uint32_t a_start = delta_start + N - 1;

// define the WEIGTHS that we will use to quantify how "costly" (bad) are each component of the COST function
// Basically HOW important is each element of the COST function: For instance, it's very important that cte remains close to 0
// but also it's veru important to make sure that the changes in commands (steering and throattle) are smooth.
// For more explanations look below on the COST function construntion

const double W_cte       = 0.75;
const double W_yaw_err   = 0.03;
const double W_vel_err   = 999999999;
const double W_steer_use = 0;
const double W_accel_use = 0;
const double W_cteSpeed  = 0;
const double W_dSteer    = 9999;
const double W_dAccel    = 0;
const double W_acc_err   = 0;

		

class FG_eval {
 public:
  // Fitted polynomial coefficients
  VectorXd coeffs;
  FG_eval(VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and the Solver function below.
	
	// The cost is stored is the first element of `fg`.
	// Any additions to the cost should be added to `fg[0]`.
	fg[0] = 0;

	// Cost function
	// TODO: Define the cost related the reference state and
	// any anything you think may be beneficial.

	/*
	What I'm trying to achieve is a Ref_Vel that changes depending on the radious of Curvature of the Road
	If I have a fitted Poly, I could calculate the RoC and with that estimate the MAX speed I should achieve 
	depending on the road conditions and the car/tire specs. In this case I'll tune it by hand and try to see 
	if I can improve the results
	*/
	
	AD<double> x = vars[x_start];
	AD<double> a = coeffs[3];
	AD<double> b = coeffs[2];
	AD<double> c = coeffs[1];

	const AD<double> ToMeters = 3; // from Pixels? to meters

	AD<double> RoC = ToMeters * fabs(CppAD::pow(1 + CppAD::pow(3 * a * CppAD::pow(x, 2) + 2 * b * x + c, 2), 1.5) / (4 * a * x + 2 * b));

	const AD<double> MaxAccel = 1.6; // m/s2: In theory, the max for confort driving is about 1.8m/s2 and beyond 3m/s2 is not even fun, but my 
	//conversion factors might be not correct and I'm adjusting for a "visual" fast response
	ref_v = CppAD::pow(MaxAccel * RoC, 0.5);
	ref_v = ref_v / 0.44704; //Converting speed back to mph

	// Saturation Control
	const AD<double> MaxTrackSpeed = 120.0;
	if (ref_v > MaxTrackSpeed) ref_v = MaxTrackSpeed;

	cout << " RoC          " << RoC << endl;
	cout << " Ref_v        " << ref_v << endl;
	

	// overriding 
	//ref_v = 70;

	// The part of the cost based on the reference state.
	for (uint32_t t = 0; t < N; t++) {
		fg[0] += W_cte * CppAD::pow(vars[cte_start + t], 2);
		fg[0] += W_yaw_err* CppAD::pow(vars[epsi_start + t], 2);
		fg[0] += W_vel_err * CppAD::pow(vars[v_start + t] - ref_v, 2);
		
	}

	// Minimize the use of actuators.	
	for (uint32_t t = 0; t < N - 1; t++) {
		fg[0] += W_steer_use * CppAD::pow(vars[delta_start + t], 2);
		fg[0] += W_accel_use * CppAD::pow(vars[a_start + t], 2);
		fg[0] += W_acc_err * CppAD::pow(vars[a_start + t] - ref_acc, 2);
		

		// As I did on the PID project, I'd like to penalize the ratio between the Cte and the speed.
		// I want to encourange HIGH speed while cte is close to 0
		AD<double> speed = vars[v_start + t];
		if (fabs(speed) > 1) {
			fg[0] += W_cteSpeed * CppAD::pow((vars[cte_start + t] / speed), 2);
		}				
	}

	// Minimize the value gap between sequential actuations. (This is actually to guarantee a min "snap" trajectory)
	// We could try to use even deeper derivatives (min Jerk trajectories), but for now we can see how this performs.
	for (uint32_t t = 0; t < N - 2; t++) {
		fg[0] += W_dSteer * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
		fg[0] += W_dAccel * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
	}

	/*
	By definition, fg has the following structure:
		fg[0] - is the accumulated COST
		1 to N*6 - The initial values

	*/

	// INITIAL VALUES
	fg[1 + x_start] = vars[x_start];
	fg[1 + y_start] = vars[y_start];
	fg[1 + psi_start] = vars[psi_start];
	fg[1 + v_start] = vars[v_start];
	fg[1 + cte_start] = vars[cte_start];
	fg[1 + epsi_start] = vars[epsi_start];

	// The rest of the constraints
	// Before we start defining the constraints, we need to recall that, as a discrete model, all the constraints
	// at time "t+1" depend on the values at "t" AND also that for simplicity we put all the constraints of the form
	// XX(t+1) = F(X(t))
	// as
	// F(X(t)) - XX(t+1) = 0 or XX(t+1) - F(X(t)) = 0
	// Therefore, we will start collecting all the actual (t+1) and previous (t) values
	for (uint32_t t = 1; t < N; t++) {
		// X
		AD<double> x_t_1 = vars[x_start + t - 1];
		AD<double> x_t = vars[x_start + t];

		// Y
		AD<double> y_t_1 = vars[y_start + t - 1];
		AD<double> y_t = vars[y_start + t];
		
		// YAW / HEADING
		AD<double> psi_t_1 = vars[psi_start + t - 1];
		AD<double> psi_t = vars[psi_start + t];
		
		// SPEED / VELOCITY MAGNITUDE
		AD<double> v_t_1 = vars[v_start + t - 1];
		AD<double> v_t = vars[v_start + t];
		
		// CTE
		AD<double> cte_t_1 = vars[cte_start + t - 1];
		AD<double> cte_t = vars[cte_start + t];
		
		// YAW ERROR
		AD<double> epsi_t_1 = vars[epsi_start + t - 1];
		AD<double> epsi_t = vars[epsi_start + t];
		
		// we are just interested in getting the previous accel (throttle) and steering
		AD<double> a_t_1 = vars[a_start + t - 1];
		AD<double> delta_t_1 = vars[delta_start + t - 1];
		
		/*
		if (t > 1) {   // use previous actuations (to account for latency)
			a = vars[a_start + t - 2];
			delta = vars[delta_start + t - 2];
		}
		*/

		
		AD<double> f_t_1 = coeffs[0] + c * x_t_1 + b * CppAD::pow(x_t_1, 2) + a * CppAD::pow(x_t_1, 3);
		AD<double> psides_t_1 = CppAD::atan(c + 2 * b * x_t_1 + 3 * a * CppAD::pow(x_t_1, 2));

		// Now we are ready to Setup the rest of the model constraints
		
		fg[1 + x_start + t]    = x_t - (x_t_1 + v_t_1 * CppAD::cos(psi_t_1) * dt);
		fg[1 + y_start + t]    = y_t - (y_t_1 + v_t_1 * CppAD::sin(psi_t_1) * dt);
		fg[1 + psi_start + t]  = psi_t - (psi_t_1 + (v_t_1 / Lf * delta_t_1 * dt));
		fg[1 + v_start + t]    = v_t - (v_t_1 + a_t_1 * dt);
		fg[1 + cte_start + t]  = cte_t - ((f_t_1 - y_t_1) + (v_t_1 * CppAD::sin(epsi_t_1) * dt));
		fg[1 + epsi_start + t] = epsi_t - ((psi_t_1 - psides_t_1) + (v_t_1/Lf * delta_t_1 * dt));

	}
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}

MPC::~MPC() {}

vector<double> MPC::Solve(VectorXd state, VectorXd coeffs) {
  bool ok = true;
 
  int StateSize = state.size();
  const int NumActuators = 2;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  // In "N" timesteps => "N - 1" actuations

  const uint32_t n_vars = N * StateSize + (N - 1) * NumActuators;

  // TODO: Set the number of constraints
  /* In theory, and since each state is 1 dimention, it needs to be bounded: Lower anf higher.
     We will assume that this variable is actually refering to the number of VARIBALES (from the State) that 
	 will have cosntrains.
	 Later we will apply the the corresponding Low and high boundary contrsiant.	 
  */
  const uint32_t n_constraints = N * StateSize;

  // For clarity
  double x = state[0]; // Always 0 since we moved to the Car Ref System
  double y = state[1]; // Always 0 since we moved to the Car Ref System
  double psi = state[2]; // Always 0 since we moved to the Car Ref System
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  // Initial State:
  // Set the initial variable values
  
  Dvector vars(n_vars);

  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  for (uint32_t i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  // Set all non-actuators (x,y,psi,v,cte,epsi) upper and lowerlimits to the max negative and positive values.
  // We can refine these limits but for simplicity we do this for now.
  for (uint32_t i = 0; i < delta_start; i++) {
	vars_lowerbound[i] = -1.0e19;
	vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits (for actuators) of delta (Steering) IS -25 to 25.
  for (uint32_t i = delta_start; i < a_start; i++) {
	  vars_lowerbound[i] = -0.436332; // 25 deg in radians
	  vars_upperbound[i] = 0.436332;
  }

  // The upper and lower limits (for actuators) of Accel (Throttle) ARE -1 to 1.
  for (uint32_t i = delta_start; i < n_vars; i++) {
	  vars_lowerbound[i] = -1.0;
	  vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (uint32_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  
  // Initial state should have smae upper and lower bounds since it will NOT change
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  cout << "Cost " << cost << endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x_t = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result;

  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);

  // This is done for visualization on the simulator. (to show a MPC trajectory)
  // We will add all the N step MPC {x,y} solutions.
  for (uint32_t i = 0; i < N - 1; i++) {
	  result.push_back(solution.x[x_start + i + 1]);
	  result.push_back(solution.x[y_start + i + 1]);
  }

  return result;
  
}
