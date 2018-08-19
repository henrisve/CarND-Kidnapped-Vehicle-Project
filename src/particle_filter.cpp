/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	std::default_random_engine generator;
	normal_distribution<double> norm_x(x, std[0]);
	normal_distribution<double> norm_y(y, std[1]);
	normal_distribution<double> norm_theta(theta, std[2]);

	num_particles=100;
	for(int i=0;i<num_particles;i++){
		Particle new_particle;
		new_particle.id = i;
		new_particle.x = norm_x(generator);
		new_particle.y = norm_y(generator);
		new_particle.theta = norm_theta(generator);
		new_particle.weight =1;
		particles.push_back(new_particle);
		weights.push_back(1);//but why again?
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	std::default_random_engine generator;
	double delta_pos = velocity * delta_t;
	double delta_yaw = yaw_rate * delta_t;		
	normal_distribution<double> norm_x(0, std[0]);
	normal_distribution<double> norm_y(0, std[1]);
	normal_distribution<double> norm_theta(0, std[2]);
	for (auto &p : particles){
		if(yaw_rate==0){
			p.theta +=  delta_yaw + norm_theta(generator);
			p.x += cos(p.theta)*delta_pos + norm_x(generator);
			p.y += sin(p.theta)*delta_pos + norm_y(generator);
		}else{
			double t0 = p.theta ;
			p.theta += delta_yaw + norm_theta(generator);
			p.x += velocity/yaw_rate * (sin(p.theta) - sin(t0))+ norm_theta(generator); //Here p.theta already have the nosie, how does this affect it? better to add after?
			p.y += velocity/yaw_rate * (cos(t0) - cos(p.theta))+ norm_theta(generator);
		}	
	}
	
}

best_location ParticleFilter::dataAssociation(LandmarkObs observed_location, double sensor_range, double std_landmark[], const Map &map_landmarks) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

//	std::vector<best_location> best_locations;

//	for(int i=0; i < observations.size(); i++){
//		auto observed_location = observations[i];
//for (auto &observed_location : observations){
	double min_distance = sensor_range;//*sensor_range;//
	best_location best_location;
	best_location.id=-1;
	int i=0;
	for (auto &map_location : map_landmarks.landmark_list){
		double xdiff = map_location.x_f - observed_location.x;
		double ydiff = map_location.y_f - observed_location.y;
		double sq_diff = pow(xdiff, 2) / 2*pow(std_landmark[0], 2) +
		              pow(ydiff, 2) / 2*pow(std_landmark[1], 2);
		double diff = sqrt(xdiff*xdiff + ydiff*ydiff);
		if(diff < min_distance){
			min_distance=diff;
			best_location.id = map_location.id_i;
			best_location.distance = sq_diff; // Use square diff to reduce computations
			cout << "new best" << best_location.id << " . " << best_location.distance << " (" << diff << ")"<< endl;
		}
		
	}
	//	best_locations.push_back(best_location);
	//}
	return best_location;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	int i=0;
	long double total_weight=0;
	for(auto &p : particles){	
		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;
		p.weight = 1;
		cout << "init weight for p" << p.id << " = " << p.weight << endl;
		for(auto &o : observations){
			LandmarkObs trueObservation;
			trueObservation.x = p.x + cos(p.theta) * o.x - sin(p.theta) * o.y;
			trueObservation.y = p.y + sin(p.theta) * o.x + cos(p.theta) * o.y;
			best_location best_locations = dataAssociation(trueObservation,sensor_range,std_landmark,map_landmarks);
			if(best_locations.id != -1){

				long double vm = 1/(2 * M_1_PI * std_landmark[0]*std_landmark[1]) * exp(-best_locations.distance);
				cout << "vm " << vm << " __ " <<  1/(2 * M_1_PI * std_landmark[0]*std_landmark[1]) << " * " <<  exp(-best_locations.distance) << endl;
				p.weight *= vm;
			}
			associations.push_back(best_locations.id);
			sense_x.push_back(trueObservation.x);
			sense_y.push_back(trueObservation.y);
			cout << "weight for p" << p.id << " o" << o.id << " = " << p.weight << endl;		
		}	
		p.associations=associations;
		p.sense_x=sense_x;
		p.sense_y=sense_y;
		total_weight += p.weight;
	}
	for(auto &p : particles){	
		p.weight = p.weight/total_weight;
		weights[i++]=p.weight;
		cout << "after weight for p" << p.id << " = " << p.weight << endl;

		//SetAssociations(p,associations,sense_x,sense_y);
	}
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());

	vector<Particle> resamples;
	for(int i=0;i<num_particles;i++) 
		//Particle new_p = ;
		resamples.push_back(particles[distribution(gen)]);
	particles=resamples;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
