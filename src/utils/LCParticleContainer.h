/*
 * LCParticleContainer.h
 *
 *  Created on: Nov 27, 2013
 *      Author: andal
 */

#ifndef LCPARTICLECONTAINER_H_
#define LCPARTICLECONTAINER_H_

#include <list>
#include <vector>
#include <Particle.h>
#include "Vector.h"
#include "ParticleIterator.h"
#include "LCInnerParticleIterator.h"
#include "LCOuterParticleIterator.h"

namespace utils {

class LCParticleContainer {
public:
	LCParticleContainer();
	virtual ~LCParticleContainer();

	/**
	 * fills the particles list
	 * @param particles_arg the list containing the particles
	 * @param domain_size_arg the domain size
	 * @param cutoff_radius_arg the cutoff_radius
	 */
	void initialize(std::list<Particle>& particles_arg,
			Vector<double, 3> domain_size_arg, double cutoff_radius_arg);

	/**
	 * @return the first element of the container
	 */
	LCOuterParticleIterator beginOuter();

	/**
	 * @return the last element of the container
	 */
	LCOuterParticleIterator endOuter();

	////////////////Parallelization methods////////////////

	/**
	 * initializes the vector of subdomains with the cellnumbers
	 * they have to work on
	 */
	void initializeCellNumbers();

	/**
	 * @param i the subdomain to work with
	 * @return the first element of the i-th subdomain
	 */
	LCOuterParticleIterator beginOuter(int i);

	/**
	 * @param i the subdomain to work with
	 * @return the last element of the i-th subdomain
	 */
	LCOuterParticleIterator endOuter(int i);
	////////////////Parallelization methods////////////////

	/**
	 * updates the list of cells (also used to initialize)
	 */
	void updateCells();

	/**
	 * Initializes the cells
	 */
	void initializeCells();

	/**
	 * initializes the halo cells around the domain
	 */
	void initializeHaloCells();

	/**
	 * initializes the cells in each boundary
	 */
	void initializeBoundaryCells();

	/**
	 * initializes the first particle for the domain and all boundaries
	 */
	void initializeBeginOuter();

	/**
	 *initializes the last particle for the domain and all boundaries
	 */
	void initializeEndOuter();

	/**
	 * @return the first element of the neighboring particles
	 * @param the outer particle iterator
	 */
	LCInnerParticleIterator beginInner(LCOuterParticleIterator it);

	/**
	 * @return the endOfInners[i]
	 * @param the index of the outer iterator
	 */
	LCInnerParticleIterator& endInner(int i);

	LCOuterParticleIterator beginBoundary(int i);

	LCOuterParticleIterator endBoundary(int i);

	LCOuterParticleIterator beginHalo(int i);

	LCOuterParticleIterator endHalo(int i);

	std::vector<std::list<Particle *> *>& getBoundaryCells(int i);

	std::vector<std::list<Particle *> *>& getHaloCells(int i);

	/**
	 * @return the particles list
	 */
	std::list<Particle *>& getList();

	/**
	 * @return the domain sizes
	 */
	utils::Vector<double, 3>& getDomainSize();

	/**
	 * @return the width
	 */
	int& getWidth();

	/**
	 * @return the height
	 */
	int& getHeight();

	/**
	 * @return the depth
	 */
	int& getDepth();

	/**
	 * @return number of particles within the domain
	 */
	int size();

	/**
	 * initalizes the ends for the inner iterators
	 */
	void initializeEndInner();

private:
	/**
	 * contains the list of the given particles
	 */
	std::list<Particle *> particles;

	/**
	 * contains the vector of cells of the given particles
	 */
	std::vector<std::list<Particle *> *> cells;

	/**
	 * contains the size of the domain.
	 * domain_size[0] - valid region of the "x-axis"
	 * domain_size[1] - valid region of the "y-axis"
	 * domain_size[2] - valid region of the "z-axis"
	 */
	utils::Vector<double, 3> domain_size;

	/**
	 * length of the edges (3D or 2D) of the cells
	 */
	double cutoff_radius;

	/**
	 * number of cells in the "x-axis" = domain_size[0] / cutoff_radius (rounded down)
	 */
	int width;

	/**
	 * number of cells in the "y-axis" = domain_size[1] / cutoff_radius (rounded down)
	 */
	int height;

	/**
	 * number of cells in the "z-axis" = domain_size[2] / cutoff_radius (rounded down)
	 */
	int depth;

	/**
	 * number of cells = width * height * depth (in 3D) or width * height (in 2D)
	 */
	int num_of_cells;

	/**
	 * the first particle of the domain
	 */
	LCOuterParticleIterator beginDomain;
	/**
	 * the last particle of the domain
	 */
	LCOuterParticleIterator endDomain;

	/**
	 * the ends of the inner iterators
	 */
	std::vector<LCInnerParticleIterator> endOfInners;

	std::vector<LCOuterParticleIterator> beginOfBoundary;

	std::vector<LCOuterParticleIterator> beginOfHalo;

	std::vector<LCOuterParticleIterator> endOfBoundary;

	std::vector<LCOuterParticleIterator> endOfHalo;

	std::vector<std::vector<std::list<Particle *> *> > haloCells;

	std::vector<std::vector<std::list<Particle *> *> > boundaryCells;


	////////////////Parallelization variables////////////////

	int num_of_threads;

	std::vector<std::vector<int> > start;

	std::vector<std::vector<int> > end;

	std::vector<LCOuterParticleIterator> beginOfSubDomain;

	std::vector<LCOuterParticleIterator> endOfSubDomain;

	std::vector<std::vector<int> > cellnumber;

	////////////////Parallelization variables////////////////

	/**
	 * checks whether the iterator is in the right boundary
	 * @param the index
	 */
	bool checkRight(int i);

	/**
	 * checks whether the iterator is in the left boundary
	 * @param the index
	 */
	bool checkLeft(int i);

	/**
	 * checks whether the iterator is in the top boundary
	 * @param the index
	 */
	bool checkTop(int i);

	/**
	 * checks whether the iterator is on the boundary
	 * @param the index
	 */
	bool checkBottom(int i);

	/**
	 * checks whether the iterator is in the front boundary
	 * @param the index
	 */
	bool checkFront(int i);
};

} /* namespace utils */

#endif /* LCPARTICLECONTAINER_H_ */
