#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;
using namespace std;
using namespace Eigen;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
VectorXd polyfit(VectorXd xvals, VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  MatrixXd A(xvals.size(), order + 1);
  
  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

MatrixXd transMap2CarRef(const double x, const double y, const double psi, const vector<double> & ptsx, const vector<double> & ptsy) {
	
	// Check for consistency
	assert(ptsx.size() == ptsy.size());
	int NumWayPoints = ptsx.size();

	//create the Matrix with 2 vectors for the waypoint x and y coordinates w.r.t. car 
	auto wp_vehRef = MatrixXd(2, NumWayPoints);

	double cos_psi = cos(-psi);
	double sin_psi = sin(-psi);
	for (int i = 0; i<NumWayPoints; ++i) {
		double dx = (ptsx[i] - x);
		double dy = (ptsy[i] - y);
		wp_vehRef(0, i) = cos_psi * dx - sin_psi * dy;
		wp_vehRef(1, i) = sin_psi * dx + cos_psi * dy;
	}
	
	return wp_vehRef;
	
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v_mph = j[1]["speed"];
		  double v = v_mph * 0.44704; // I'm converting into m/s because I'm caculating the Radious of curature in meters inside the MPC, 
		  // and generatung a ref_v using a MaxAccel in m/s2
		  double steering = j[1]["steering_angle"]; // in degrees
		  double throttle = j[1]["throttle"]; // -1 to 1


          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          double steer_cmd;
          double throttle_cmd;

		  // *************************************************
		  
		  // Since ptsx and ptsy are on MAP reference but the rest of the visualization (mpc_x_, mpc_y_, next_x,...)
		  // are in Car coordinates, it seems only reasonable to convert them all to the same system.
		  // The easiest way (only one conversion) is to put ptsx and ptsy on CAR coordinates.
		  // When we do so, px = py = psi = 0 (from the Car Ref Sys)
		  // I create this variables simply for clarity and teaching purposes. All could be simplified by a 0
		  double CarRef_px, CarRef_py, CarRef_psi;
		  CarRef_px = CarRef_py = CarRef_psi = 0.0;

		  MatrixXd wp_vehRef = transMap2CarRef(px, py, psi, ptsx, ptsy);
		  VectorXd wp_vehRef_x = wp_vehRef.row(0);
		  VectorXd wp_vehRef_y = wp_vehRef.row(1);

		  

		  // fit a 3rd order polynomial to the waypoints
		  auto coeffs = polyfit(wp_vehRef_x, wp_vehRef_y, 3);

		  // get cross-track error from fit
		  // In vehicle coordinates the cross-track error is the intercept at x = 0, That means that since we have a fit of the form:
		  // C0 + C1*X + C2*X2 + C3X3
		  // when we evaluate using x=0 we just get C0.
		  //But to understand where this is coming from I like to keep the whole evaluation, even though this is exactly C0
		  double cte = polyeval(coeffs, CarRef_px) - CarRef_py;


		  // get orientation error from fit ( Since we are trying a 3rd order poly, then, f' = a1 + 2*a2*x + 3*a3*x2)
		  // in this case and since we moved our reference sys to the Car, px = 0 and also psi = 0 
		  //double epsi = psi - atan(coeffs[1] + 2 * coeffs[2] * px + 3 * coeffs[3] * pow(px, 2));
		  double epsi = CarRef_psi - atan(coeffs[1]);

		  // state in vehicle coordinates: x,y and orientation are always zero
		  VectorXd state(6);
		  /*
		  I can send the ACTUAL state to the MPC or I can try to compensate for the latency by "predicting" what would 
		  be the state after the latency period.
		  */
		  const double latency = 0.1; // 100 ms

		  // Let's predict the state. Rembember that px, py and psi wrt car are all 0. 
		  const double pred_px = v * latency;
		  const double pred_py = 0;
		  const double pred_psi = -v * steering * latency / Lf;
		  const double pred_v = v + throttle * latency;
		  const double pred_cte = cte + v * sin(epsi) * latency;
		  const double pred_epsi = epsi + pred_psi;


		  //state << CarRef_px, CarRef_py, CarRef_psi, v, cte, epsi;
		  state << pred_px, pred_py, pred_psi, pred_v, pred_cte, pred_epsi;


		  // compute the optimal trajectory          
		  auto mpc_solution = mpc.Solve(state, coeffs);
		  steer_cmd = mpc_solution[0]; // This should be in dregrees since I used degrees before I sent it to the MPC
		  throttle_cmd = mpc_solution[1];

		  // Just for help...
		  //cout << " x            " << px << endl;
		  //cout << " y            " << py << endl;
		  //cout << " psi          " << psi << endl;
		  cout << " v            " << v_mph << endl;
		  cout << " cte          " << cte << endl;
		  //cout << " epsi         " << epsi << endl;
		  cout << " steer_cmd    " << steer_cmd << endl;
		  cout << " throttle_cmd " << throttle_cmd << endl;

		  //*****************************************************
          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
		  msgJson["steering_angle"] = -steer_cmd / (deg2rad(25));
          msgJson["throttle"] = throttle_cmd;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
		  // the points in the simulator are connected by a Green line

		  // For questions, take a look at the MPC return values.
		  for (uint32_t i = 2; i < mpc_solution.size(); i++) {
			  if (i % 2 == 0) {
				  mpc_x_vals.push_back(mpc_solution[i]);
			  }
			  else {
				  mpc_y_vals.push_back(mpc_solution[i]);
			  }
		  }
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line
		  for (unsigned i = 0; i < wp_vehRef_x.size(); ++i) {
			  next_x_vals.push_back(wp_vehRef_x(i));
			  next_y_vals.push_back(wp_vehRef_y(i));
		  }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //cout << msg << endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(int(latency*1000)));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
    const string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    cout << "Connected!!!" << endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    cout << "Disconnected" << endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    cout << "Listening to port " << port << endl;
  } else {
    cerr << "Failed to listen to port" << endl;
    return -1;
  }
  h.run();
}
