/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include <sstream>
#include <iostream>
#include <log4cxx/logger.h>
#include <log4cxx/propertyconfigurator.h>
#include <log4cxx/xml/domconfigurator.h>

log4cxx::LoggerPtr particlelogger(log4cxx::Logger::getLogger("particle"));

Particle::Particle(int type_arg) {
	type = type_arg;
	LOG4CXX_INFO(particlelogger, "Arrived @ Type Constructor.");
	f = 0.0;
	old_f = 0.0;
	parID = 0;
}

Particle::Particle(const Particle& other) {
	x = other.x;
	x0 = other.x0;
	v = other.v;
	f = other.f;
	old_f = other.old_f;
	temp_f = 0.0;
	m = other.m;
	type = other.type;
	parID = other.parID;
	directNeighbours = other.directNeighbours;
	diagNeighbours = other.diagNeighbours;
	LOG4CXX_TRACE(particlelogger, "Particle generated by copy!");
}

// Todo: maybe use initializater list instead of copy?
Particle::Particle(utils::Vector<double, 3> x_arg,
		utils::Vector<double, 3> v_arg, double m_arg, int type_arg, int id) {
	x = x_arg;
	x0 = x_arg;
	v = v_arg;
	m = m_arg;
	type = type_arg;
	f = 0.0;
	old_f = 0.0;
	parID = id;
	LOG4CXX_TRACE(particlelogger, "Particle generated!");
}

Particle::~Particle() {
	LOG4CXX_TRACE(particlelogger, "Particle destructed!");
}

utils::Vector<double, 3>& Particle::getX() {
	return x;
}

utils::Vector<double, 3>& Particle::getX0() {
	return x0;
}

utils::Vector<double, 3>& Particle::getV() {
	return v;
}

utils::Vector<double, 3>& Particle::getF() {
	return f;
}

utils::Vector<double, 3>& Particle::getOldF() {
	return old_f;
}

double& Particle::getM() {
	return m;
}

utils::Vector<double, 3> Particle::getTempF() {
	return temp_f;
}

void Particle::updateTempF(utils::Vector<double, 3> newF) {
#pragma omp critical
	{
		temp_f = temp_f + newF;
	}
}

int& Particle::getType() {
	return type;
}

void Particle::setF(utils::Vector<double, 3> newF) {
	old_f = f;
	f = newF;
}

void Particle::setDirectNeighbours(std::list<int> direct){
	directNeighbours = direct;
}

void Particle::setDiagNeighbours(std::list<int> diag){
	diagNeighbours = diag;
}

void Particle::deleteTempF() {
	temp_f = 0.0;
}

int& Particle::getID(){
        return parID;
}

bool Particle::isDirectNeighbour(Particle& p){
        //ID is unique for each particle
        for (std::list<int>::iterator it = directNeighbours.begin();
                        it != directNeighbours.end(); it++){
                if ((*it) == p.getID())
                        return true;
        }
        return false;
}

bool Particle::isDiagNeighbour(Particle& p){
        //ID is unique for each particle
        for (std::list<int>::iterator it = diagNeighbours.begin();
                        it != diagNeighbours.end(); it++){
                if ((*it) == p.getID())
                        return true;
        }
        return false;
}

std::list<int>& Particle::getDirectNeighbours(){
	return directNeighbours;
}

std::list<int>& Particle::getDiagNeighbours(){
	return diagNeighbours;
}

std::string Particle::toString() {
	std::stringstream stream;
	stream << "Particle: X:" << x << " v: " << v << " f: " << f << " old_f: "
			<< old_f << " type: " << type;
	return stream.str();
}

bool Particle::operator ==(Particle& other) {
	if ((x == other.x) && (v == other.v) && (f == other.f)
			&& (type == other.type) && (m == other.m)
			&& (old_f == other.old_f)) {
		return true;
	}

	return false;
}

std::ostream & operator<<(std::ostream & stream, Particle & p) {
	stream << p.toString();
	return stream;
}
