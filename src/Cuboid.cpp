/*
 * Cub.cpp
 *
 *  Created on: 10.11.2013
 *      Author: son
 */

#include "Cuboid.h"
#include <vector>
#include <Particle.h>
#include "utils/Vector.h"
#include "MaxwellBoltzmannDistribution.h"
#include <cassert>

Cuboid::Cuboid() {

}

Cuboid::Cuboid(int height, int width, int depth, double distance, double mass,
		utils::Vector<double, 3> ori, utils::Vector<double, 3> startVelocity,
		double meanVelocity, int parType, double EPSILON, double SIGMA) {
	// Initialize first variables
	cHeight = height;
	cWidth = width;
	cDepth = depth;
	origin = ori;
	h = distance;
	m = mass;
	startV = startVelocity;
	meanV = meanVelocity;
	this->parType = parType;
	this->EPSILON = EPSILON;
	this->SIGMA = SIGMA;

	// Initialize cub
	cub.clear();

	// Initialize particles in cub
	for (int hei = 0; hei < height; hei++) {
		for (int w = 0; w < width; w++) {
			for (int d = 0; d < depth; d++) {
				// Must set each particle with its own coordinate
				// Ox along width, Oy along height, Oz along depth
				double addTemp[] = { w * h, hei * h, d * h };
				utils::Vector<double, 3> addVector(addTemp);
				utils::Vector<double, 3> vel(ori + addVector);

				int id = d * height * width + hei * width + w;
				Particle p(vel, startVelocity, mass, parType, id);

				// Movement of each particle superposed by Brownian Motion
//				MaxwellBoltzmannDistribution(p, meanVelocity, 2);

				cub.push_back(p);
			}
		}
	}
}

void Cuboid::initNeighbours() {
	//works only for 2D membranes
	Particle * pp = NULL;

	for (std::list<Particle>::iterator it = cub.begin(); it != cub.end();
			it++) {
		Particle& p = *it;
		int id = p.getID();
		bool isFirst = ((id % cWidth) == 0);
		bool isLast = (((id + 1) % cWidth) == 0);
		std::list<int> directNeighbours;
		std::list<int> diagNeighbours;

		//direct under
		int newID = id - cWidth;
		if (newID >= 0)
			directNeighbours.push_back(newID);

		//direct left
		newID = id - 1;
		if (id % cWidth != 0) {
			directNeighbours.push_back(newID);
		}

		//direct right
		newID = id + 1;
		if (newID % cWidth != 0) {
			directNeighbours.push_back(newID);
		}

		//direct above
		newID = id + cWidth;
		if (newID < cWidth * cHeight) {
			directNeighbours.push_back(newID);
		}

		//set the direct neighbours
		p.setDirectNeighbours(directNeighbours);

		//diagonal lower left
		newID = id - cWidth - 1;
		if (!isFirst && newID >= 0) {
			diagNeighbours.push_back(newID);
		}

		//diagonal lower right
		newID = id - cWidth + 1;
		if (!isLast && newID >= 0) {
			diagNeighbours.push_back(newID);
		}

		//diagonal upper left
		newID = id + cWidth - 1;
		if (!isFirst && newID < cWidth * cHeight) {
			diagNeighbours.push_back(newID);
		}

		//diagonal upper right
		newID = id + cWidth + 1;
		if (!isLast && newID < cWidth * cHeight) {
			diagNeighbours.push_back(newID);
		}

		//set the diagonal neighbours
		p.setDiagNeighbours(diagNeighbours);
	}

}

Particle * Cuboid::getParticleAtID(int id){
	//id starts from 0
	if ((id < 0) || (id >= cub.size()))
		return NULL;

	for (std::list<Particle>::iterator it = cub.begin();
			it != cub.end(); it++){
		if ((*it).getID() == id)
			return &(*it);
	}

	return NULL;
}

utils::Vector<double, 3>& Cuboid::getOrigin() {
	return origin;
}

utils::Vector<double, 3>& Cuboid::getStartV() {
	return startV;
}

std::list<Particle>& Cuboid::getCuboid() {
	return cub;
}

int Cuboid::getHeight() {
	return cHeight;
}

void Cuboid::setHeight(double newH) {
	cHeight = newH;
}

int Cuboid::getWidth() {
	return cWidth;
}

void Cuboid::setWidth(double newW) {
	cWidth = newW;
}

int Cuboid::getDepth() {
	return cDepth;
}

void Cuboid::setDepth(double newD) {
	cDepth = newD;
}

double Cuboid::getDistance() {
	return h;
}

void Cuboid::setDistance(double newD) {
	h = newD;
}

double Cuboid::getMass() {
	return m;
}

void Cuboid::setMass(double newM) {
	m = newM;
}

double Cuboid::getMeanV() {
	return meanV;
}

void Cuboid::setMeanV(double newV) {
	meanV = newV;
}

int Cuboid::getSize() {
	return cHeight * cWidth * cDepth;
}

int& Cuboid::getType() {
	return parType;
}

double& Cuboid::getEpsilon() {
	return EPSILON;
}

double& Cuboid::getSigma() {
	return SIGMA;
}

Cuboid::~Cuboid() {
}
