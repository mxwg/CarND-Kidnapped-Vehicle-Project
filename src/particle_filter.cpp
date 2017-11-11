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
    std::cout << "Initializing PF with: " << x << ", " << y << ", " << theta << "\n";
    std::cout << "Creating " << this->num_particles << " particles...\n";

    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    particles.clear();
    weights.clear();

    for (size_t i = 0; i < this->num_particles; ++i) {
        Particle particle;
        particle.id = i;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1.0;
        particles.push_back(particle);
        weights.push_back(1.0);
    }

    std::cout << "First particle: " << particles[0].x << ", " << particles[0].y << ", " << particles[0].theta << "\n";
    this->is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;

    for(size_t i = 0; i < num_particles; ++i)
    {
        double x = particles[i].x;
        double y = particles[i].y;
        double theta = particles[i].theta;
        double new_x, new_y, new_theta;

        // predict new state
        if (fabs(yaw_rate) > 0.0000001) {
            new_x = x + (velocity / yaw_rate) * (sin(theta + yaw_rate * delta_t) - sin(theta));
            new_y = y + (velocity / yaw_rate) * (cos(theta) - cos(theta + yaw_rate * delta_t));
            new_theta = (theta + yaw_rate * delta_t);
        } else{
            new_x = x + velocity * delta_t * cos(theta);
            new_y = y + velocity * delta_t * sin(theta);
            new_theta = theta;
        }

        // add noise
        normal_distribution<double> dist_x(new_x, std_pos[0]);
        normal_distribution<double> dist_y(new_y, std_pos[1]);
        normal_distribution<double> dist_t(new_theta, std_pos[2]);
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_t(gen);
        particles[i].id = i;
    }
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

    double sum_weights = 0;
    for (size_t i = 0; i < this->num_particles; ++i) {
        std::vector<LandmarkObs> predicted;
        double current_weight = 1.0;
        particles[i].associations.clear();
        particles[i].sense_x.clear();
        particles[i].sense_y.clear();
        for (auto obs : observations) {
            // transform observation into map frame (using orientation and position of the particle in the map)
            transform(obs.x, obs.y, particles[i].x, particles[i].y, particles[i].theta);
            particles[i].sense_x.push_back(obs.x);
            particles[i].sense_y.push_back(obs.y);

            // find nearest neighbor
            int id = -1;
            double lmx = 0;
            double lmy = 0;
            double nearest = 1e10;
            for (auto landmark: map_landmarks.landmark_list)
            {
                double d = dist(landmark.x_f, landmark.y_f, obs.x, obs.y);
                if( d > sensor_range)
                {
                    continue;   // ignore measurements out of sensor range
                }
                if (d < nearest)
                {
                    id = landmark.id_i;
                    nearest = d;
                    lmx = landmark.x_f;
                    lmy = landmark.y_f;
                }
            }
            // save the nearest landmark
            particles[i].associations.push_back(id);
            double p = mvg(obs.x, obs.y, lmx, lmy, std_landmark[0], std_landmark[1]);
            current_weight *= p;
        }

        // update weights
        particles[i].weight = current_weight;
        weights[i] = current_weight;
        sum_weights += current_weight;
    }

    // normalize weights
    for(size_t i = 0; i < num_particles; ++i)
    {
        particles[i].weight /= sum_weights;
        weights[i] /= sum_weights;
    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    default_random_engine gen;
    std::discrete_distribution<int> dist(weights.begin(), weights.end());

    std::vector<Particle> newParticles;
    for(size_t i = 0; i < num_particles; ++i)
    {
        newParticles.push_back(this->particles[dist(gen)]);
    }
    this->particles = newParticles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x,
                                         std::vector<double> sense_y) {
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best) {
    vector<int> v = best.associations;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseX(Particle best) {
    vector<double> v = best.sense_x;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseY(Particle best) {
    vector<double> v = best.sense_y;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

void ParticleFilter::transform(double &x, double &y, double tX, double tY, double theta) {
    double st = sin(theta);
    double ct = cos(theta);
    double new_x = x * ct - y * st + tX;
    double new_y = x * st + y * ct + tY;
    x = new_x;
    y = new_y;
}

double ParticleFilter::mvg(double x, double y, double mx, double my, double sx, double sy) {
    double xd = x-mx;
    double yd = y-my;
    double factor = 1.0 / (2 * M_PI * sx * sy);
    double exponent = - ((xd*xd)/(2*sx*sx) + (yd*yd)/(2*sy*sy) );

    return factor * exp(exponent);
}

void ParticleFilter::print(std::string where) {
    std::cout << "STATE: " << where << " ###########################\n";
    for(size_t i = 0; i < num_particles; ++i)
    {
        std::cout << particles[i].id << ": "
                                     << particles[i].x << ", "
                  << particles[i].y << " w(p)=" << particles[i].weight << ", w=" << weights[i]
                  << "\n";
    }
    std::cout << std::endl;
}
