/*
 * BoundaryHandler.cpp
 *
 *  Created on: Dec 8, 2013
 *      Author: andal
 */

#include "BoundaryHandler.h"

namespace utils {

BoundaryHandler::BoundaryHandler(std::vector<int> boundary_type_arg,
		LCParticleContainer& container_arg, double h_arg,
		void calculate_arg(Particle&, Particle&)) :
		container(container_arg), domain_size(container_arg.getDomainSize()), width(
				container_arg.getWidth()), height(container_arg.getHeight()), depth(
				container_arg.getDepth()) {
	boundary_type = boundary_type_arg;
	h = h_arg;
	calculate = calculate_arg;
	widhei = width * height;
	widdep = width * depth;
	heidep = height * depth;
	i = 0;
	j = 0;
	k = 0;
	l = 0;
}

BoundaryHandler::~BoundaryHandler() {
	// TODO Auto-generated destructor stub
}

void BoundaryHandler::applyPeriodicMoving() {
	//links, unten, vorne sind inklusiv Grenze,
	//rechts, oben, hinten sind exklusive Grenze
	if (boundary_type[0] == 3) {
		iterator = container.beginBoundary(0);
		while (iterator != container.endBoundary(0)) {
			if (((*iterator).getX())[0] < 0) {
				utils::Vector<double, 3> tempX = (*iterator).getX();
				tempX[0] = tempX[0] + domain_size[0];
				(*iterator).getX() = tempX;
			}
			assert(((*iterator).getX())[0] >= 0);
//			assert(((*iterator).getX())[0] < domain_size[0]);
//			assert(((*iterator).getX())[1] >= 0);
//			assert(((*iterator).getX())[1] < domain_size[1]);

			++iterator;
		}
	}
	if (boundary_type[1] == 3) {
		iterator = container.beginBoundary(1);
		while (iterator != container.endBoundary(1)) {
			if (((*iterator).getX())[0] >= domain_size[0]) {
				utils::Vector<double, 3> tempX = (*iterator).getX();
				tempX[0] = tempX[0] - domain_size[0];
				(*iterator).getX() = tempX;
			}

			assert(((*iterator).getX())[0] >= 0);
			assert(((*iterator).getX())[0] < domain_size[0]);
//			assert((*iterator).getX()[1] >= 0);
//			assert((*iterator).getX()[1] < domain_size[1]);
			++iterator;
		}
	}
	if (boundary_type[2] == 3) {
		iterator = container.beginBoundary(2);
		while (iterator != container.endBoundary(2)) {
			if (((*iterator).getX())[1] < 0) {
				utils::Vector<double, 3> tempX = (*iterator).getX();
				tempX[1] = tempX[1] + domain_size[1];
				(*iterator).getX() = tempX;
			}

			assert((*iterator).getX()[0] >= 0);
			assert((*iterator).getX()[0] < domain_size[0]);
			assert((*iterator).getX()[1] >= 0);
//			assert((*iterator).getX()[1] < domain_size[1]);
			++iterator;
		}
	}
	if (boundary_type[3] == 3) {
		iterator = container.beginBoundary(3);
		while (iterator != container.endBoundary(3)) {
			if (((*iterator).getX())[1] >= domain_size[1]) {
				utils::Vector<double, 3> tempX = (*iterator).getX();
				tempX[1] = tempX[1] - domain_size[1];
				(*iterator).getX() = tempX;
			}
			assert((*iterator).getX()[0] >= 0);
			assert((*iterator).getX()[0] < domain_size[0]);
			assert((*iterator).getX()[1] >= 0);
			assert((*iterator).getX()[1] < domain_size[1]);
			++iterator;
		}
	}
	if (depth > 0) {
		if (boundary_type[4] == 3) {
			iterator = container.beginBoundary(4);
			while (iterator != container.endBoundary(4)) {
				if (((*iterator).getX())[2] >= domain_size[2]) {
					utils::Vector<double, 3> tempX = (*iterator).getX();
					tempX[2] = tempX[2] - domain_size[2];
					(*iterator).getX() = tempX;
				}
				++iterator;
			}
		}
		if (boundary_type[5] == 3) {
			iterator = container.beginBoundary(5);
			while (iterator != container.endBoundary(5)) {
				if (((*iterator).getX())[2] < 0) {
					utils::Vector<double, 3> tempX = (*iterator).getX();
					tempX[2] = tempX[2] + domain_size[2];
					(*iterator).getX() = tempX;
				}
				++iterator;
			}
		}
	}
}

void BoundaryHandler::applyReflecting() {
	utils::Vector<double, 3> v(0.0);
	if (boundary_type[0] == 2) {
		iterator = container.beginBoundary(0);
		while (iterator != container.endBoundary(0)) {
			if ((*iterator).getX()[0] <= h) {
				double x_arg[3] = { 0, ((*iterator).getX())[1],
						((*iterator).getX())[2] };
				utils::Vector<double, 3> x(x_arg);
				Particle p(x, v, 0, (*iterator).getType(), -1);
				calculate((*iterator), p);
			}
			++iterator;
		}
	}

	if (boundary_type[1] == 2) {
		iterator = container.beginBoundary(1);
		while (iterator != container.endBoundary(1)) {
			if ((*iterator).getX()[0] >= domain_size[0] - h) {
				double x_arg[3] = { domain_size[0], ((*iterator).getX())[1],
						((*iterator).getX())[2] };
				utils::Vector<double, 3> x(x_arg);
				Particle p(x, v, 0, (*iterator).getType(), -1);
				calculate((*iterator), p);
			}
			++iterator;
		}
	}

	if (boundary_type[2] == 2) {
		iterator = container.beginBoundary(2);
		while (iterator != container.endBoundary(2)) {
			if ((*iterator).getX()[1] <= h) {
				double x_arg[3] = { ((*iterator).getX())[0], 0,
						((*iterator).getX())[2] };
				utils::Vector<double, 3> x(x_arg);
				Particle p(x, v, 0, (*iterator).getType(), -1);
				calculate((*iterator), p);
			}
			++iterator;
		}
	}

	if (boundary_type[3] == 2) {
		iterator = container.beginBoundary(3);
		while (iterator != container.endBoundary(3)) {
			if ((*iterator).getX()[1] >= domain_size[1] - h) {
				double x_arg[3] = { ((*iterator).getX())[0], domain_size[1],
						((*iterator).getX())[2] };
				utils::Vector<double, 3> x(x_arg);
				Particle p(x, v, 0, (*iterator).getType(), -1);
				calculate((*iterator), p);
			}
			++iterator;
		}
	}

	if (depth > 0) {

		if (boundary_type[4] == 2) {
			iterator = container.beginBoundary(4);
			while (iterator != container.endBoundary(4)) {
				if ((*iterator).getX()[2] >= domain_size[2] - h) {
					double x_arg[3] = { (*iterator).getX()[0],
							(*iterator).getX()[1], domain_size[2] };
					utils::Vector<double, 3> x(x_arg);
					Particle p(x, v, (*iterator).getM(), (*iterator).getType(),
							-1);
					calculate((*iterator), p);
				}
				++iterator;
			}
		}

		if (boundary_type[5] == 2) {
			iterator = container.beginBoundary(5);
			while (iterator != container.endBoundary(5)) {
				if ((*iterator).getX()[2] <= h) {
					double x_arg[3] = { (*iterator).getX()[0],
							(*iterator).getX()[1], 0 };
					utils::Vector<double, 3> x(x_arg);
					Particle p(x, v, (*iterator).getM(), (*iterator).getType(),
							-1);
					calculate((*iterator), p);
				}
				++iterator;
			}
		}
	}

}

void BoundaryHandler::applyOutflow() {
	if (boundary_type[0] == 1) {
		iterator = container.beginBoundary(0);
		while (iterator != container.endBoundary(0)) {
			if ((*iterator).getX()[0] <= 0) {
				(container.getList()).remove(&(*iterator));
			}
			++iterator;
		}
	}
	if (boundary_type[1] == 1) {
		iterator = container.beginBoundary(1);
		while (iterator != container.endBoundary(1)) {
			if ((*iterator).getX()[0] >= domain_size[0]) {
				(container.getList()).remove(&(*iterator));
			}
			++iterator;
		}
	}
	if (boundary_type[2] == 1) {
		iterator = container.beginBoundary(2);
		while (iterator != container.endBoundary(2)) {
			if ((*iterator).getX()[1] <= 0) {
				(container.getList()).remove(&(*iterator));
			}
			++iterator;
		}
	}
	if (boundary_type[3] == 1) {
		iterator = container.beginBoundary(3);
		while (iterator != container.endBoundary(3)) {
			if ((*iterator).getX()[1] >= domain_size[1]) {
				(container.getList()).remove(&(*iterator));
			}
			++iterator;
		}
	}
	if (depth > 0) {
		if (boundary_type[4] == 1) {
			iterator = container.beginBoundary(4);
			while (iterator != container.endBoundary(4)) {
				if ((*iterator).getX()[2] < 0) {
					(container.getList()).remove(&(*iterator));
				}
				++iterator;
			}
		}
		if (boundary_type[5] == 1) {
			iterator = container.beginBoundary(5);
			while (iterator != container.endBoundary(5)) {
				if ((*iterator).getX()[2] >= domain_size[2]) {
					(container.getList()).remove(&(*iterator));
				}
				++iterator;
			}
		}
	}
}

void BoundaryHandler::applyPeriodic() {
	if (boundary_type[0] == 3) {
		if (boundary_type[2] == 3) {
			if (boundary_type[4] == 3) {
				/* Periodicity for every boundary */
				iterator = container.beginBoundary(0);
				if (depth > 1) { /* 3D */
					while (iterator != container.endBoundary(0)) {
						i = iterator.getCellNumber();
						j = i - height - 1;
						/* Left wall */
						leftSide3D();
						/* Edge on the bottom left */
						if (i % height == 0) {
							j = i - 1;
							leftLowerEdge3D();
						}
						/* Edge on the top left */
						else if ((i + 1) % height == 0) {
							j = i + 1;
							leftUpperEdge3D();
						}

						/* Edge on the left back */
						if (i < height) {
							j = i + heidep - height - 1;
							leftBackEdge3D();
							/* Corner on the left bottom back */
							if (i == 0) {
								j = heidep - 1;
								leftLowerBackCorner3D();
							}
							/* Corner on the left top back */
							else if (i == height - 1) {
								j = heidep - height;
								leftUpperBackCorner3D();
							}
						}
						/* Edge on the left front */
						else if (i >= heidep - height) {
							j = (i % height) - 1;
							leftFrontEdge3D();
							/* Corner on the left bottom front */
							if (i == heidep - height) {
								j = height - 1;
								leftLowerFrontCorner3D();
							}
							/* Corner on the left top front */
							else if (i == heidep - 1) {
								j = 0;
								leftUpperFrontCorner3D();
							}
						}
						++iterator;
					}
				} else { /* 2D */
					while (iterator != container.endBoundary(0)) {
						i = iterator.getCellNumber();
						j = i - 1;
						/* Left wall */
						leftSide2D();
						/* Corner on the left bottom */
						if (i == 0) {
							j = height - 1;
							leftLowerCorner2D();
						}
						/* Corner on the left top */
						if (i == height - 1) {
							j = 0;
							leftUpperCorner2D();
						}
						++iterator;
					}
				}
				iterator = container.beginBoundary(2);
				if (depth > 1) { /* 3D */
					while (iterator != container.endBoundary(2)) {
						i = iterator.getCellNumber();
						j = i - width - 1;
						/* Lower wall */
						lowerSide3D();
						/* Edge on the lower back */
						if (i < width) {
							j = i + widdep - width - 1;
							lowerBackEdge3D();
						}
						/* Edge on the lower front */
						if (i >= widdep - width) {
							j = (i % width) - 1;
							lowerFrontEdge3D();
						}
						/* All corners already calculated */
						++iterator;
					}
				} else { /* 2D */
					while (iterator != container.endBoundary(2)) {
						i = iterator.getCellNumber();
						j = i - 1;
						/* Lower wall */
						lowerSide2D();
						/* All corners already calculated */
						++iterator;
					}
				}

				if (depth > 2) { /* 3D */
					iterator = container.beginBoundary(5);
					while (iterator != container.endBoundary(5)) {
						i = iterator.getCellNumber();
						j = i - width - 1;
						/* Rear wall */
						backSide3D();
						/* All edges and corners already calculated */
						++iterator;
					}
				}
			} else {
				/* Periodicity for the left, right, lower and upper boundaries */
				iterator = container.beginBoundary(0);
				if (depth > 1) { /* 3D */
					while (iterator != container.endBoundary(0)) {
						i = iterator.getCellNumber();
						j = i - height - 1;
						/* Left wall */
						leftSide3D();
						/* Edge on the left bottom */
						if (i % height == 0) {
							j = i - 1;
							leftLowerEdge3D();
						}
						/* Edge on the left top */
						if ((i + 1) % height == 0) {
							j = i + 1;
							leftUpperEdge3D();
						}
						++iterator;
					}
				} else { /* 2D */
					while (iterator != container.endBoundary(0)) {
						i = iterator.getCellNumber();
						j = i - 1;
						/* Left wall */
						leftSide2D();
						/* Corner on the left bottom */
						if (i == 0) {
							j = height - 1;
							leftLowerCorner2D();
						}
						/* Corner on the left top */
						else if (i == height - 1) {
							j = 0;
							leftUpperCorner2D();
						}
						++iterator;
					}
				}
				iterator = container.beginBoundary(2);
				if (depth > 1) { /* 3D */
					while (iterator != container.endBoundary(2)) {
						i = iterator.getCellNumber();
						j = i - width - 1;
						/* Lower wall */
						lowerSide3D();
						/* All edges already calculated */
						++iterator;
					}
				} else { /* 2D */
					while (iterator != container.endBoundary(2)) {
						i = iterator.getCellNumber();
						j = i - 1;
						/* Lower wall */
						lowerSide2D();
						/* All edges already calculated */
						++iterator;
					}
				}
			}
		} else if (boundary_type[4] == 3) {
			/* Periodicity for the left, right, front and rear boundaries */
			iterator = container.beginBoundary(0);
			if (depth > 1) { /* 3D */
				while (iterator != container.endBoundary(0)) {
					i = iterator.getCellNumber();
					j = i - height - 1;
					/* Left  wall */
					leftSide3D();
					/* Edge on the left back */
					if (i < height) {
						j = i + heidep - height - 1;
						leftBackEdge3D();
					}
					/* Edge on the left front */
					if (i >= heidep - height) {
						j = (i % height) - 1;
						leftFrontEdge3D();
					}
					++iterator;
				}
			} else { /* 2D */
				while (iterator != container.endBoundary(0)) {
					i = iterator.getCellNumber();
					j = i - 1;
					/* Left wall */
					leftSide2D();
					/* These edges do not exist */
					++iterator;
				}
			}
			if (depth > 2) { /* 3D */
				iterator = container.beginBoundary(5);
				while (iterator != container.endBoundary(5)) {
					i = iterator.getCellNumber();
					j = i - width - 1;
					/* Rear wall */
					backSide3D();
					/* All edges already calculated */
					++iterator;
				}
			}
		} else {
			/* Periodicity for the left and right boundaries */
			iterator = container.beginBoundary(0);
			if (depth > 1) { /* 3D */
				while (iterator != container.endBoundary(0)) {
					i = iterator.getCellNumber();
					j = i - height - 1;
					/* Left wall */
					leftSide3D();
					/* All edges already calculated */
					++iterator;
				}
			} else { /* 2D */
				while (iterator != container.endBoundary(0)) {
					i = iterator.getCellNumber();
					j = i - 1;
					/* Left wall */
					leftSide2D();
					/* All edges already calculated */
					++iterator;
				}
			}
		}
	} else if (boundary_type[2] == 3) {
		if (boundary_type[4] == 3) {
			/* Periodicity for the lower, upper, front and rear boundaries */
			iterator = container.beginBoundary(2);
			if (depth > 1) { /* 3D */
				while (iterator != container.endBoundary(2)) {
					i = iterator.getCellNumber();
					j = i - width - 1;
					/* Lower wall */
					lowerSide3D();
					/* Edge on the lower back */
					if (i < width) {
						j = i + widdep - width - 1;
						lowerBackEdge3D();
					}
					/* Edge on the lower front */
					if (i >= widdep - width) {
						j = (i % width) - 1;
						lowerFrontEdge3D();
					}
					++iterator;
				}
			} else { /* 2D */
				while (iterator != container.endBoundary(2)) {
					i = iterator.getCellNumber();
					j = i - 1;
					/* Lower wall */
					lowerSide2D();
					/* These edges do not exist */
					++iterator;
				}
			}
			if (depth > 2) { /* 3D */
				iterator = container.beginBoundary(5);
				while (iterator != container.endBoundary(5)) {
					i = iterator.getCellNumber();
					j = i - width - 1;
					backSide3D();
					/* All edges already calculated */
					++iterator;
				}
			}
		} else {
			/* Periodicity for the lower and upper boundaries */
			iterator = container.beginBoundary(2);
			if (depth > 1) { /* 3D */
				while (iterator != container.endBoundary(2)) {
					i = iterator.getCellNumber();
					j = i - width - 1;
					/* Lower wall */
					lowerSide3D();
					/* All edges already calculated */
					++iterator;
				}
			} else { /* 2D */
				while (iterator != container.endBoundary(2)) {
					i = iterator.getCellNumber();
					j = i - 1;
					/* Lower wall */
					lowerSide2D();
					/* All edges already calculated */
					++iterator;
				}
			}
		}
	} else if (boundary_type[4] == 3) {
		/* Periodicity for front and rear boundaries */
		if (depth > 2) { /* 3D */
			iterator = container.beginBoundary(5);
			while (iterator != container.endBoundary(5)) {
				i = iterator.getCellNumber();
				j = i - width - 1;
				/* Rear wall */
				backSide3D();

				++iterator;
			}
			/* All edges already calculated */
		}
	}
}

void BoundaryHandler::leftSide3D() {
	for (k = 0; k < 3; k++) {
		for (l = 0; l < 3; l++) {
			if (j >= 0 && j < heidep
					&& !(container.getHaloCells(0))[j]->empty()) {
				it2 = (container.getHaloCells(0))[j]->begin();
				while (it2 != (container.getHaloCells(0))[j]->end()) {

					utils::Vector<double, 3> tempX = (*it2)->getX();
					tempX[0] = tempX[0] - ((double) domain_size[0]);
					(*it2)->getX() = tempX;

					calculate((*(*it2)), *iterator);

					tempX[0] = tempX[0] + ((double) domain_size[0]);
					(*it2)->getX() = tempX;

					++it2;
				}
			}
			++j;
		}
		j = j + height - 3;
	}
}

void BoundaryHandler::lowerSide3D() {
	for (k = 0; k < 3; k++) {
		for (l = 0; l < 3; l++) {
			if (j >= 0 && j < widdep
					&& !(container.getHaloCells(2))[j]->empty()) {
				it2 = (container.getHaloCells(2))[j]->begin();
				while (it2 != (container.getHaloCells(2))[j]->end()) {
					utils::Vector<double, 3> tempX = (*it2)->getX();
					tempX[1] = tempX[1] - ((double) domain_size[1]);
					(*it2)->getX() = tempX;

					calculate((*(*it2)), *iterator);

					tempX[1] = tempX[1] + ((double) domain_size[1]);
					(*it2)->getX() = tempX;

					++it2;
				}
			}
			++j;
		}
		j = j + width - 3;
	}
}
void BoundaryHandler::backSide3D() {
	for (k = 0; k < 3; k++) {
		for (l = 0; l < 3; l++) {
			if (j >= 0 && j < widhei
					&& !(container.getHaloCells(5))[j]->empty()) {
				it2 = (container.getHaloCells(5))[j]->begin();
				while (it2 != (container.getHaloCells(5))[j]->end()) {

					utils::Vector<double, 3> tempX = (*it2)->getX();
					tempX[2] = tempX[2] - ((double) domain_size[2]);
					(*it2)->getX() = tempX;

					calculate((*(*it2)), *iterator);

					tempX[2] = tempX[2] + ((double) domain_size[2]);
					(*it2)->getX() = tempX;

					++it2;
				}
			}
			++j;
		}
		j = j + width - 3;
	}
}
void BoundaryHandler::leftBackEdge3D() {
	for (k = 0; k < 3; k++) {
		if (j >= heidep - height && j < heidep
				&& !((container.getHaloCells(0))[j]->empty())) {
			it2 = (container.getHaloCells(0))[j]->begin();
			while (it2 != (container.getHaloCells(0))[j]->end()) {

				utils::Vector<double, 3> tempX = (*it2)->getX();
				tempX[0] = tempX[0] - ((double) domain_size[0]);
				tempX[2] = tempX[2] - ((double) domain_size[2]);
				(*it2)->getX() = tempX;

				calculate((*(*it2)), *iterator);

				tempX[0] = tempX[0] + ((double) domain_size[0]);
				tempX[2] = tempX[2] + ((double) domain_size[2]);
				(*it2)->getX() = tempX;

				++it2;
			}
		}
		++j;
	}
}
void BoundaryHandler::leftFrontEdge3D() {
	//TODO
	for (k = 0; k < 3; k++) {
		if (j >= 0 && j < height
				&& !(container.getHaloCells(0))[j]->empty()) {
			it2 = (container.getHaloCells(0))[j]->begin();
			while (it2 != (container.getHaloCells(0))[j]->end()) {

				utils::Vector<double, 3> tempX = (*it2)->getX();
				tempX[0] = tempX[0] - ((double) domain_size[0]);
				tempX[2] = tempX[2] + ((double) domain_size[2]);
				(*it2)->getX() = tempX;

				calculate((*(*it2)), *iterator);

				tempX[0] = tempX[0] + ((double) domain_size[0]);
				tempX[2] = tempX[2] - ((double) domain_size[2]);
				(*it2)->getX() = tempX;

				++it2;
			}
		}
		++j;
	}
}
void BoundaryHandler::leftLowerEdge3D() {
	for (k = 0; k < 3; k++) {
		if (j >= 0 && j < heidep
				&& !(container.getHaloCells(0))[j]->empty()) {
			it2 = (container.getHaloCells(0))[j]->begin();
			while (it2 != (container.getHaloCells(0))[j]->end()) {

				utils::Vector<double, 3> tempX = (*it2)->getX();
				tempX[0] = tempX[0] - ((double) domain_size[0]);
				tempX[1] = tempX[1] - ((double) domain_size[1]);
				(*it2)->getX() = tempX;

				calculate((*(*it2)), *iterator);

				tempX[0] = tempX[0] + ((double) domain_size[0]);
				tempX[1] = tempX[1] + ((double) domain_size[1]);
				(*it2)->getX() = tempX;

				++it2;
			}
		}
		j = j + height;
	}
}
void BoundaryHandler::leftUpperEdge3D() {
	for (k = 0; k < 3; k++) {
		if (j >= 0 && j < heidep
				&& !(container.getHaloCells(0))[j]->empty()) {
			it2 = (container.getHaloCells(0))[j]->begin();
			while (it2 != (container.getHaloCells(0))[j]->end()) {

				utils::Vector<double, 3> tempX = (*it2)->getX();
				tempX[0] = tempX[0] - ((double) domain_size[0]);
				tempX[1] = tempX[1] + ((double) domain_size[1]);
				(*it2)->getX() = tempX;

				calculate((*(*it2)), *iterator);

				tempX[0] = tempX[0] + ((double) domain_size[0]);
				tempX[1] = tempX[1] - ((double) domain_size[1]);
				(*it2)->getX() = tempX;

				++it2;
			}
		}
		j = j - height;
	}
}
void BoundaryHandler::lowerBackEdge3D() {
	for (k = 0; k < 3; k++) {
		if (j >= widdep - width && j < widdep
				&& !(container.getHaloCells(2))[j]->empty()) {
			it2 = (container.getHaloCells(2))[j]->begin();
			while (it2 != (container.getHaloCells(2))[j]->end()) {
				utils::Vector<double, 3> tempX = (*it2)->getX();
				tempX[1] = tempX[1] - ((double) domain_size[1]);
				tempX[2] = tempX[2] - ((double) domain_size[2]);
				(*it2)->getX() = tempX;

				calculate((*(*it2)), *iterator);

				tempX[1] = tempX[1] + ((double) domain_size[1]);
				tempX[2] = tempX[2] + ((double) domain_size[2]);
				(*it2)->getX() = tempX;

				++it2;
			}
		}
		++j;
	}
}
void BoundaryHandler::lowerFrontEdge3D() {
	for (k = 0; k < 3; k++) {
		if (j >= 0 && j < width
				&& !(container.getHaloCells(2))[j]->empty()) {
			it2 = (container.getHaloCells(2))[j]->begin();
			while (it2 != (container.getHaloCells(2))[j]->end()) {
				utils::Vector<double, 3> tempX = (*it2)->getX();
				tempX[1] = tempX[1] - ((double) domain_size[1]);
				tempX[2] = tempX[2] + ((double) domain_size[2]);
				(*it2)->getX() = tempX;

				calculate((*(*it2)), *iterator);

				tempX[1] = tempX[1] + ((double) domain_size[1]);
				tempX[2] = tempX[2] - ((double) domain_size[2]);
				(*it2)->getX() = tempX;

				++it2;
			}
		}
		++j;
	}
}
void BoundaryHandler::leftLowerBackCorner3D() {
	if (!(container.getHaloCells(0))[j]->empty()) {
		it2 = (container.getHaloCells(0))[j]->begin();
		while (it2 != (container.getHaloCells(0))[j]->end()) {

			utils::Vector<double, 3> tempX = (*it2)->getX();
			tempX[0] = tempX[0] - ((double) domain_size[0]);
			tempX[1] = tempX[1] - ((double) domain_size[1]);
			tempX[2] = tempX[2] - ((double) domain_size[2]);
			(*it2)->getX() = tempX;

			calculate((*(*it2)), *iterator);

			tempX[0] = tempX[0] + ((double) domain_size[0]);
			tempX[1] = tempX[1] + ((double) domain_size[1]);
			tempX[2] = tempX[2] + ((double) domain_size[2]);
			(*it2)->getX() = tempX;

			++it2;
		}
	}
}
void BoundaryHandler::leftLowerFrontCorner3D() {
	if (!(container.getHaloCells(0))[j]->empty()) {
		it2 = (container.getHaloCells(0))[j]->begin();
		while (it2 != (container.getHaloCells(0))[j]->end()) {

			utils::Vector<double, 3> tempX = (*it2)->getX();
			tempX[0] = tempX[0] - ((double) domain_size[0]);
			tempX[1] = tempX[1] - ((double) domain_size[1]);
			tempX[2] = tempX[2] + ((double) domain_size[2]);
			(*it2)->getX() = tempX;

			calculate((*(*it2)), *iterator);

			tempX[0] = tempX[0] + ((double) domain_size[0]);
			tempX[1] = tempX[1] + ((double) domain_size[1]);
			tempX[2] = tempX[2] - ((double) domain_size[2]);
			(*it2)->getX() = tempX;

			++it2;
		}
	}
}
void BoundaryHandler::leftUpperBackCorner3D() {
	if (!(container.getHaloCells(0))[j]->empty()) {
		it2 = (container.getHaloCells(0))[j]->begin();
		while (it2 != (container.getHaloCells(0))[j]->end()) {

			utils::Vector<double, 3> tempX = (*it2)->getX();
			tempX[0] = tempX[0] - ((double) domain_size[0]);
			tempX[1] = tempX[1] + ((double) domain_size[1]);
			tempX[2] = tempX[2] - ((double) domain_size[2]);
			(*it2)->getX() = tempX;

			calculate((*(*it2)), *iterator);

			tempX[0] = tempX[0] + ((double) domain_size[0]);
			tempX[1] = tempX[1] - ((double) domain_size[1]);
			tempX[2] = tempX[2] + ((double) domain_size[2]);
			(*it2)->getX() = tempX;

			++it2;
		}
	}
}
void BoundaryHandler::leftUpperFrontCorner3D() {
	if (!(container.getHaloCells(0))[j]->empty()) {
		it2 = (container.getHaloCells(0))[j]->begin();
		while (it2 != (container.getHaloCells(0))[j]->end()) {

			utils::Vector<double, 3> tempX = (*it2)->getX();
			tempX[0] = tempX[0] - ((double) domain_size[0]);
			tempX[1] = tempX[1] + ((double) domain_size[1]);
			tempX[2] = tempX[2] + ((double) domain_size[2]);
			(*it2)->getX() = tempX;

			calculate((*(*it2)), *iterator);

			tempX[0] = tempX[0] + ((double) domain_size[0]);
			tempX[1] = tempX[1] - ((double) domain_size[1]);
			tempX[2] = tempX[2] - ((double) domain_size[2]);
			(*it2)->getX() = tempX;

			++it2;
		}
	}
}

void BoundaryHandler::leftSide2D() {
	for (k = 0; k < 3; k++) {
		if (j >= 0 && j < height
				&& !((container.getHaloCells(0))[j]->empty())) {

			it2 = (container.getHaloCells(0))[j]->begin();
			while (it2 != ((container.getHaloCells(0))[j])->end()) {

				utils::Vector<double, 3> tempX = (*it2)->getX();
				tempX[0] = tempX[0] - ((double) domain_size[0]);
				(*it2)->getX() = tempX;

				calculate((*(*it2)), *iterator);

				tempX[0] = tempX[0] + ((double) domain_size[0]);
				(*it2)->getX() = tempX;

				++it2;
			}
		}
		++j;
	}
}
void BoundaryHandler::lowerSide2D() {
	for (k = 0; k < 3; k++) {
		if (j >= 0 && j < width
				&& !(container.getHaloCells(2))[j]->empty()) {
			it2 = (container.getHaloCells(2))[j]->begin();
			while (it2 != (container.getHaloCells(2))[j]->end()) {

				utils::Vector<double, 3> tempX = (*it2)->getX();
				tempX[1] = tempX[1] - ((double) domain_size[1]);
				(*it2)->getX() = tempX;

				calculate((*(*it2)), *iterator);

				tempX[1] = tempX[1] + ((double) domain_size[1]);
				(*it2)->getX() = tempX;

				++it2;
			}
		}
		++j;
	}
}
void BoundaryHandler::leftLowerCorner2D() {
	if (!(container.getHaloCells(0))[j]->empty()) {
		it2 = (container.getHaloCells(0))[j]->begin();
		while (it2 != ((container.getHaloCells(0))[j])->end()) {

			utils::Vector<double, 3> tempX = (*it2)->getX();
			tempX[0] = tempX[0] - ((double) domain_size[0]);
			tempX[1] = tempX[1] - ((double) domain_size[1]);
			(*it2)->getX() = tempX;

			calculate((*(*it2)), *iterator);

			tempX[0] = tempX[0] + ((double) domain_size[0]);
			tempX[1] = tempX[1] + ((double) domain_size[1]);
			(*it2)->getX() = tempX;

			++it2;
		}
	}
}
void BoundaryHandler::leftUpperCorner2D() {
	if (!(container.getHaloCells(0))[j]->empty()) {
		it2 = (container.getHaloCells(0))[j]->begin();
		while (it2 != ((container.getHaloCells(0))[j])->end()) {

			utils::Vector<double, 3> tempX = (*it2)->getX();
			tempX[0] = tempX[0] - ((double) domain_size[0]);
			tempX[1] = tempX[1] + ((double) domain_size[1]);
			(*it2)->getX() = tempX;

			calculate((*(*it2)), *iterator);

			tempX[0] = tempX[0] + ((double) domain_size[0]);
			tempX[1] = tempX[1] - ((double) domain_size[1]);
			(*it2)->getX() = tempX;

			++it2;
		}
	}
}

} /* namespace utils */
