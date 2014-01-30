/*
 * LCParticleContainer.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: andal
 */

#include "LCParticleContainer.h"
#include <log4cxx/logger.h>
#include <math.h>
#include <omp.h>

log4cxx::LoggerPtr lcparticlecontainerlogger(
		log4cxx::Logger::getLogger("utils.lcparticlecontainer"));

namespace utils {

LCParticleContainer::LCParticleContainer() {
	// TODO Auto-generated constructor stub

}

LCParticleContainer::~LCParticleContainer() {
	// TODO Auto-generated destructor stub
	for (int i = 0; i < num_of_cells; i++) {
		delete cells[i];
	}

}

void LCParticleContainer::initialize(std::list<Particle>& particles_arg,
		Vector<double, 3> domain_size_arg, double cutoff_radius_arg) {

	//fill the particle list
	std::list<Particle>::iterator iterator = particles_arg.begin();
	while (iterator != particles_arg.end()) {
		particles.push_back(&(*iterator));
		++iterator;
	}

	//set the domain sizes and the cutoff radius
	domain_size = domain_size_arg;
	cutoff_radius = cutoff_radius_arg;

	//if the cutoff radius is no divisor of the domain sizes,
	//we have to adjust it(make it bigger)
	width = domain_size[0] / cutoff_radius;
	if (fmod(domain_size[0], cutoff_radius) != 0) {
		width++;
		domain_size[0] = width * cutoff_radius;
	}

	height = domain_size[1] / cutoff_radius;
	if (fmod(domain_size[1], cutoff_radius) != 0) {
		height++;
		domain_size[1] = height * cutoff_radius;
	}

	depth = domain_size[2] / cutoff_radius;
	if (fmod(domain_size[2], cutoff_radius) != 0) {
		depth++;
		domain_size[2] = depth * cutoff_radius;
	}

	//initialize the number of cells needed in our array
	if (depth > 0) {
		num_of_cells = width * height * depth;
	} else {
		num_of_cells = width * height;
	}
	cells.resize(num_of_cells);
	endOfInners.resize(num_of_cells);
	for (int i = 0; i < num_of_cells; i++) {
		cells[i] = new std::list<Particle *>();
	}

	//initialize the boundaries and halos
	int n = depth > 0 ? 6 : 4;

	boundaryCells.resize(n);
	haloCells.resize(n);
	beginOfBoundary.resize(n);
	endOfBoundary.resize(n);
	beginOfHalo.resize(n);
	endOfHalo.resize(n);

	//initialize the cells with the corresponding particles
	initializeCells();

	//if we are working parallel, initialize the parallelization variables
	cellnumber = std::vector<std::vector<int> >();
#ifdef _OPENMP
	num_of_threads = omp_get_max_threads();
	std::cout << num_of_threads << " threads exist " << std::endl;

	initializeCellNumbers();
#endif

	initializeBoundaryCells();

	initializeHaloCells();

	initializeBeginOuter();

	initializeEndOuter();

	initializeEndInner();

}

void LCParticleContainer::updateCells() {
	//initialize the cells with the corresponding particles
	initializeCells();

	initializeBeginOuter();

	initializeEndOuter();

	initializeEndInner();

}

void LCParticleContainer::initializeCells() {
	for (int i = 0; i < num_of_cells; i++) {
		cells[i]->clear();
	}
	std::list<Particle *>::iterator iterator = particles.begin();
	while (iterator != particles.end()) {
		//check if it's in the domain
		if (((*(*iterator)).getX())[0] <= domain_size[0]
				&& (*(*iterator)).getX()[1] <= domain_size[1]
				&& (*(*iterator)).getX()[2] <= domain_size[2]
				&& (*(*iterator)).getX()[0] >= 0
				&& ((*(*iterator)).getX())[1] >= 0
				&& (*(*iterator)).getX()[2] >= 0) {

			Particle& p = (*(*iterator));
			int i = (int) (p.getX()[0] / cutoff_radius)
					+ ((int) (p.getX()[1] / cutoff_radius)) * width
					+ ((int) (p.getX()[2] / cutoff_radius)) * width * height;
			cells[i]->push_back(*iterator);
		} else {
		}
		++iterator;
	}
	LOG4CXX_DEBUG(lcparticlecontainerlogger, "Initialized cells.");
}

void LCParticleContainer::initializeCellNumbers() {

	std::cout << num_of_threads << " threads exist 2" << std::endl;
	/*
	 * variables, to set the subCuboid of each subdomain
	 */

	//split along the depth
	if (depth >= width && depth >= height) {
		std::cout << "using depth" << num_of_threads << std::endl;
		if (num_of_threads >= depth)
			num_of_threads = depth - 1;

		std::cout << num_of_threads << " threads exist3 " << std::endl;
		cellnumber.resize(num_of_threads);
		start.resize(num_of_threads);
		end.resize(num_of_threads);
		beginOfSubDomain.resize(num_of_threads);
		endOfSubDomain.resize(num_of_threads);

		int j = 0;
		for (int i = 0; i < num_of_threads; i++) {
			(start[i]).resize(3);
			(start[i])[0] = 0;
			(start[i])[1] = 0;
			(start[i])[2] = j;
			std::cout << "Thread " << i << " : Start[" << (start[i])[0] << ", "
					<< (start[i])[1] << ", " << (start[i])[2] << " ]"
					<< std::endl;
			j = (depth * (i + 1)) / num_of_threads;
			(end[i]).resize(3);
			(end[i])[0] = width;
			(end[i])[1] = height;
			(end[i])[2] = j;
			++j;
		}
	}
	//split along the width
	else if (width >= height && width >= depth) {
		std::cout << "using width" << std::endl;
		if (num_of_threads >= width)
			num_of_threads = width - 1;

		cellnumber.resize(num_of_threads);
		start.resize(num_of_threads);
		end.resize(num_of_threads);
		beginOfSubDomain.resize(num_of_threads);
		endOfSubDomain.resize(num_of_threads);

		int j = 0;
		for (int i = 0; i < num_of_threads; i++) {
			(start[i]).resize(3);
			(start[i])[0] = j;
			(start[i])[1] = 0;
			(start[i])[2] = 0;
			std::cout << "Thread " << i << " : Start[" << (start[i])[0] << ", "
					<< (start[i])[1] << ", " << (start[i])[2] << " ]"
					<< std::endl;
			j = (width * (i + 1)) / num_of_threads;
			(end[i]).resize(3);
			(end[i])[0] = j;
			(end[i])[1] = height;
			(end[i])[2] = depth;
			++j;
		}
	}
	//split along the height
	else {
		std::cout << "using height" << std::endl;
		if (num_of_threads >= height)
			num_of_threads = height - 1;

		cellnumber.resize(num_of_threads);
		start.resize(num_of_threads);
		end.resize(num_of_threads);
		beginOfSubDomain.resize(num_of_threads);
		endOfSubDomain.resize(num_of_threads);

		int j = 0;
		for (int i = 0; i < num_of_threads; i++) {
			(start[i]).resize(3);
			(start[i])[0] = 0;
			(start[i])[1] = j;
			(start[i])[2] = 0;
			std::cout << "Thread " << i << " : Start[" << (start[i])[0] << ", "
					<< (start[i])[1] << ", " << (start[i])[2] << " ]"
					<< std::endl;
			j = (height * (i + 1)) / num_of_threads;
			(end[i]).resize(3);
			(end[i])[0] = width;
			(end[i])[1] = j;
			(end[i])[2] = depth;
			++j;
		}
	}

	//for each subdomain fill in the cell indexes
	for (int i = 0; i < num_of_threads; i++) {
		for (int index = 0; index < cells.size(); index++) {

			int width_index = index % width;			//the width of the index
			int height_index = (index / width) % height;//the height of the index
			int depth_index = index / (width * height);	//the depth of the index

			if (width_index >= (start[i])[0] && width_index <= (end[i])[0]
					&& height_index >= (start[i])[1]
					&& height_index <= (end[i])[1]
					&& depth_index >= (start[i])[2]
					&& depth_index <= (end[i])[2]) {
				(cellnumber[i]).push_back(index);
			}
		}
		std::cout << cellnumber[i].size() << std::endl;
	}

#ifndef NDEBUG
	for (int i = 0; i < num_of_cells; i++) {
		bool in = false;
		bool tooOftenIn = false;
		for (int j = 0; j < num_of_threads; j++) {
			for (int k = 0; k < cellnumber[j].size(); k++) {
				if ((cellnumber[j])[k] == i) {
					if (in == false) {
						in = true;
					} else {
						tooOftenIn = true;
					}
				}
			}
		}
		if (in == false) {
			std::cout << i << " ist nicht enthalten" << std::endl;
		}
		if (tooOftenIn == true) {
			std::cout << i << " ist zu oft enthalten" << std::endl;
		}
		assert(in == true);
		assert(tooOftenIn == false);
	}
#endif

}

void LCParticleContainer::initializeBoundaryCells() {

	//initialize vector of left bondary cells
	for (int i = 0; i < num_of_cells; i = i + width) {
		boundaryCells[0].push_back(cells[i]);
	}

	//initialize vector of right boundary cells
	for (int i = width - 1; i < num_of_cells; i = i + width) {
		boundaryCells[1].push_back(cells[i]);
	}

	//initialize vector of bottom boundary cells
	for (int i = 0; i < num_of_cells; i = i + width * height) {
		for (int j = 0; j < width; j++) {
			boundaryCells[2].push_back(cells[i + j]);
		}
	}

	//initialize vector of top boundary cells
	for (int i = width * (height - 1); i < num_of_cells;
			i = i + width * height) {
		for (int j = 0; j < width; j++) {
			boundaryCells[3].push_back(cells[i + j]);
		}
	}

	//only needed if it's a 3D simulation
	if (depth > 0) {
		//initialize vector of front boundary cells
		for (int i = num_of_cells - width * height; i < num_of_cells; i++) {
			boundaryCells[4].push_back(cells[i]);
		}

		//initialize vector of back boundary cells
		for (int i = 0; i < width * height; i++) {
			boundaryCells[5].push_back(cells[i]);
		}
	}
}

void LCParticleContainer::initializeHaloCells() {
	haloCells[0] = boundaryCells[1];
	haloCells[1] = boundaryCells[0];
	haloCells[2] = boundaryCells[3];
	haloCells[3] = boundaryCells[2];
	if (depth > 0) {
		haloCells[4] = boundaryCells[5];
		haloCells[5] = boundaryCells[4];
	}
}

void LCParticleContainer::initializeBeginOuter() {
	//initialize the beginning of the domain

	int i = 0;
	while (cells[i]->empty() == true) {
		++i;
	}
	assert(cells[i]->empty() == false);
	beginDomain = LCOuterParticleIterator(num_of_cells, &cells,
			cells[i]->begin(), i);

#ifdef _OPENMP
	for(int j = 0; j < num_of_threads; j++) {
		int i = 0;
		while(cells[(cellnumber[j])[i]]->empty() && i < cellnumber[j].size()) {
			++i;
		}
		if (i >= cellnumber[j].size()) {
			i = 0;
		}
		assert(&(cellnumber[j]) != 0);
		beginOfSubDomain[j] = LCOuterParticleIterator(num_of_cells, &cells, &(cellnumber[j]), cells[(cellnumber[j])[i]]->begin(), i);
	}
#endif

	//
	// Initialize boundary cells:
	//

	/*
	 * if depth > 0 we need 6 dimensions, otherwise we only need 4 dimensions
	 */
	int n = depth > 0 ? 6 : 4;

	/*
	 * i ~ 1 => leftBoundary
	 * i ~ 2 => rightBoundary
	 * i ~ 3 => bottomBoundary
	 * i ~ 4 => topBoundary
	 * i ~ 5 => frontBoundary
	 * i ~ 6 => backBoundary
	 *
	 */

	int j = 0;
	for (int i = 0; i < n; i++) {
		j = 0;
		while (j < (boundaryCells[i]).size()
				&& (boundaryCells[i])[j]->empty() == true) {
			++j;
		}
		if (j >= (boundaryCells[i]).size()) {
			j = 0;
		}
		beginOfBoundary[i] = LCOuterParticleIterator((boundaryCells[i]).size(),
				&(boundaryCells[i]), (boundaryCells[i])[j]->begin(), j);
	}

	beginOfHalo[0] = beginOfBoundary[1];
	beginOfHalo[1] = beginOfBoundary[0];
	beginOfHalo[2] = beginOfBoundary[3];
	beginOfHalo[3] = beginOfBoundary[2];
	if (depth > 0) {
		beginOfHalo[4] = beginOfBoundary[5];
		beginOfHalo[5] = beginOfBoundary[4];
	}
}

void LCParticleContainer::initializeEndInner() {
	for (int i = beginOuter().getCellNumber(); i <= endOuter().getCellNumber();
			i++) {
		if (!cells[i]->empty()) {
			int x;
			if (depth > 1) {
				if (num_of_cells > i + width * height + width + 1
						&& !cells[i + width * height + width + 1]->empty()
						&& checkRight(i) && checkFront(i) && checkTop(i)) {
					x = i + width * height + width + 1;
				} else if (num_of_cells > i + width * height + width
						&& !cells[i + width * height + width]->empty()
						&& checkFront(i) && checkTop(i)) {
					x = i + width * height + width;
				} else if (num_of_cells > i + width * height + width - 1
						&& !cells[i + width * height + width - 1]->empty()
						&& checkLeft(i) && checkFront(i) && checkTop(i)) {
					x = i + width * height + width - 1;
				} else if (num_of_cells > i + width * height + 1
						&& !cells[i + width * height + 1]->empty()
						&& checkRight(i) && checkFront(i)) {
					x = i + width * height + 1;
				} else if ((num_of_cells > i + width * height
						&& !cells[i + width * height]->empty())
						&& checkFront(i)) {
					x = i + width * height;
				} else if ((num_of_cells > i + width * height - 1
						&& !cells[i + width * height - 1]->empty())
						&& checkLeft(i) && checkFront(i)) {
					x = i + width * height - 1;
				} else if (num_of_cells > i + width * height - width + 1
						&& !cells[i + width * height - width + 1]->empty()
						&& checkRight(i) && checkFront(i) && checkBottom(i)) {
					x = i + width * height - width + 1;
				} else if ((num_of_cells > i + width * height - width
						&& !cells[i + width * height - width]->empty())
						&& checkFront(i) && checkBottom(i)) {
					x = i + width * height - width;
				} else if ((num_of_cells > i + width * height - width - 1
						&& !cells[i + width * height - width - 1]->empty()
						&& checkLeft(i) && checkFront(i) && checkBottom(i))) {
					x = i + width * height - width - 1;
				} else if (num_of_cells > i + width + 1
						&& !cells[i + width + 1]->empty() && checkRight(i)
						&& checkTop(i)) {
					x = i + width + 1;
				} else if (num_of_cells > i + width
						&& !cells[i + width]->empty() && checkTop(i)) {
					x = i + width;
				} else if (num_of_cells > i + width - 1
						&& !cells[i + width - 1]->empty() && checkLeft(i)
						&& checkTop(i)) {
					x = i + width - 1;
				} else if (num_of_cells > i + 1 && !cells[i + 1]->empty()
						&& checkRight(i)) {
					x = i + 1;
				} else {
					x = i;
				}
			} else {
				if (num_of_cells > i + width + 1
						&& !cells[i + width + 1]->empty() && checkTop(i)
						&& checkRight(i)) {
					x = i + width + 1;
				} else if (num_of_cells > i + width
						&& !cells[i + width]->empty() && checkTop(i)) {
					x = i + width;
				} else if (num_of_cells > i + width
						&& !cells[i + width - 1]->empty() && checkTop(i)
						&& checkLeft(i)) {
					x = i + width - 1;
				} else if (num_of_cells > i + 1 && !cells[i + 1]->empty()
						&& checkRight(i)) {
					x = i + 1;
				} else {
					x = i;
				}
			}

			assert(!cells[x]->empty());
			assert(cells[x]->size() > 0);
			assert(cells[i]->size() != 0);

			LCInnerParticleIterator ip(x, x, num_of_cells, width, height, depth,
					cells[x]->end(), &cells);

			endOfInners[i] = ip;
		}
	}
}

void LCParticleContainer::initializeEndOuter() {
	//
	// initialize the end of the domain
	//

	int i = num_of_cells - 1;
	while (i >= 0 && cells[i]->empty() == true) {
		assert(cells[i]->size() == 0);
		--i;
	}
	if (i < 0) {
		endDomain = LCOuterParticleIterator(num_of_cells, &cells,
				cells[0]->begin(), 0);
	} else {
		endDomain = LCOuterParticleIterator(num_of_cells, &cells,
				cells[i]->end(), i);
	}

#ifdef _OPENMP

	for (int j = 0; j < num_of_threads; j++) {
		int i = cellnumber[j].size() - 1;
		while (cells[(cellnumber[j])[i]]->empty() && i >= 0) {
			--i;
		}
		if (i < 0) {
			endOfSubDomain[j] = LCOuterParticleIterator(
					num_of_cells, &cells, &(cellnumber[j]),
					cells[(cellnumber[j])[0]]->begin(), 0);
		} else {
			assert(i > 0);
			assert(cells[(cellnumber[j])[i]]->empty() == false);
			endOfSubDomain[j] = LCOuterParticleIterator(
					num_of_cells, &cells, &(cellnumber[j]),
					cells[(cellnumber[j])[i]]->end(), i);
		}
	}
#endif

	/*
	 * if depth > 0 we need 6 dimensions, otherwise we only need 4 dimensions
	 */
	int n = depth > 0 ? 6 : 4;

	/*
	 * i ~ 1 => leftBoundary
	 * i ~ 2 => rightBoundary
	 * i ~ 3 => bottomBoundary
	 * i ~ 4 => topBoundary
	 * i ~ 5 => frontBoundary
	 * i ~ 6 => backBoundary
	 *
	 */
	for (int i = 0; i < n; i++) {
		int j = boundaryCells[i].size() - 1;
		while (j >= 0 && (boundaryCells[i])[j]->empty()) {
			assert((boundaryCells[i])[j]->size() == 0);
			--j;
		}
		if (j < 0) {
			endOfBoundary[i] = LCOuterParticleIterator(
					(boundaryCells[i]).size(), &(boundaryCells[i]),
					(boundaryCells[i])[0]->begin(), 0);
		} else {
			endOfBoundary[i] = LCOuterParticleIterator(
					(boundaryCells[i]).size(), &(boundaryCells[i]),
					(boundaryCells[i])[j]->end(), j);
		}
	}

	endOfHalo[0] = endOfBoundary[1];
	endOfHalo[1] = endOfBoundary[0];
	endOfHalo[2] = endOfBoundary[3];
	endOfHalo[3] = endOfBoundary[2];
	if (depth > 0) {
		endOfHalo[4] = endOfBoundary[5];
		endOfHalo[5] = endOfBoundary[4];
	}
}

std::vector<std::list<Particle *> *>& LCParticleContainer::getBoundaryCells(
		int i) {
	return boundaryCells[i];
}

std::vector<std::list<Particle *> *>& LCParticleContainer::getHaloCells(int i) {
	return haloCells[i];
}

LCOuterParticleIterator LCParticleContainer::beginBoundary(int i) {
	return beginOfBoundary[i];
}

LCOuterParticleIterator LCParticleContainer::endBoundary(int i) {
	return endOfBoundary[i];
}

LCOuterParticleIterator LCParticleContainer::beginHalo(int i) {
	return beginOfHalo[i];
}

LCOuterParticleIterator LCParticleContainer::endHalo(int i) {
	return endOfHalo[i];
}

LCOuterParticleIterator LCParticleContainer::beginOuter() {
	return beginDomain;
}

LCOuterParticleIterator LCParticleContainer::endOuter() {
	return endDomain;
}

LCOuterParticleIterator LCParticleContainer::beginOuter(int i) {
	return beginOfSubDomain[i];
}

LCOuterParticleIterator LCParticleContainer::endOuter(int i) {
	return endOfSubDomain[i];
}

LCInnerParticleIterator& LCParticleContainer::endInner(int i) {
	return endOfInners[i];
}

LCInnerParticleIterator LCParticleContainer::beginInner(
		LCOuterParticleIterator it) {
	int i = it.getCellNumber();
	assert(cells[i]->empty() == false);
	std::list<Particle *>::iterator iterator = it.getIterator();
	LCInnerParticleIterator inner(i, i, num_of_cells, width, height, depth,
			iterator, &cells);
	assert(inner != endOfInners[i]);
	return inner;
}

std::list<Particle *>& LCParticleContainer::getList() {
	return particles;
}

utils::Vector<double, 3>& LCParticleContainer::getDomainSize() {
	return domain_size;
}

int& LCParticleContainer::getWidth() {
	return width;
}

int& LCParticleContainer::getHeight() {
	return height;
}

int& LCParticleContainer::getDepth() {
	return depth;
}

int LCParticleContainer::size() {
	int size = 0;
	for (int i = 0; i < num_of_cells; i++) {
		size += cells[i]->size();
	}
	return size;
}

bool LCParticleContainer::checkLeft(int i) {
	return i % width > 0;
}
bool LCParticleContainer::checkRight(int i) {
	return (i + 1) % width > 0;
}
bool LCParticleContainer::checkBottom(int i) {
	return (width <= (i % (width * height)));
}
bool LCParticleContainer::checkFront(int i) {
	if (depth < 2) {
		return false;
	} else {
		return (i < (depth - 1) * (width * height));
	}
}
bool LCParticleContainer::checkTop(int i) {
	return i % (width * height) < width * height - width;
}

} /* namespace utils */
