/*
 * LCParticleIterator.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: andal
 */

#include "LCInnerParticleIterator.h"

namespace utils {

LCInnerParticleIterator::LCInnerParticleIterator() {
	//Constructor only for resize
}

LCInnerParticleIterator::LCInnerParticleIterator(int index_arg,
		int original_index_arg, int cell_size_arg, int width_arg,
		int height_arg, int depth_arg,
		std::list<Particle *>::iterator iterator_arg,
		std::vector<std::list<Particle *> *>* cells_arg) :
		cells(cells_arg) {
	original_index = original_index_arg;
	index = index_arg;
	cell_size = cell_size_arg;
	width = width_arg;
	height = height_arg;
	depth = depth_arg;
	iterator = iterator_arg;

	depth > 1 ? neighbour.resize(13) : neighbour.resize(4);
	neighbour[0] = original_index + 1;			// right neighbour
	neighbour[1] = original_index + width - 1;	// left above neighbour
	neighbour[2] = neighbour[1] + 1;			// above neighbour
	neighbour[3] = neighbour[2] + 1;			// right above neighbour
	if (depth > 1) {

		/*      begin neighbours in the back: 		*/
		neighbour[4] = original_index + width * height - width - 1;	// lower left
		neighbour[5] = neighbour[4] + 1;						// lower middle
		neighbour[6] = neighbour[5] + 1;						// lower right
		neighbour[7] = original_index + width * height - 1;				// left
		neighbour[8] = neighbour[7] + 1;							// behind
		neighbour[9] = neighbour[8] + 1;								// right
		neighbour[10] = original_index + width * height + width - 1;// upper left
		neighbour[11] = neighbour[10] + 1;						// upper middle
		neighbour[12] = neighbour[11] + 1;						// upper right
		/*      end neighbours in the back: 		*/
	}
}

LCInnerParticleIterator::~LCInnerParticleIterator() {
	// TODO Auto-generated destructor stub
}

Particle& LCInnerParticleIterator::operator*() const {
	return *(*iterator);
}

void LCInnerParticleIterator::operator++() {

	++iterator;
	/**
	 * Checks whether the outer particle was already the last Particle in its cell
	 * which would mean iterator were now on the dummy end of the list.
	 */
	if (iterator != (*cells)[index]->end()) {
		assert(iterator != (*cells)[index]->end());
	} else {
		int old_index = index;
		bool done = false;
		/**
		 * Checks the neighoring cells of the outer particle in a predetermined order:
		 * right, top left, top, top right [for both dimensional cases],
		 * Back bottom left, back bottom, back bottom right, back left, back, back right,
		 * back top left, back top and back top right [for the 3-dimensional case]
		 */
		while (index < cell_size && (done == false || (*cells)[index]->empty())) {
			done = false;
			if (index == original_index) {
				index++;
				if (checkRight()) {
					done = true;
				}
			} else if (index == neighbour[0]) {
				index = neighbour[1];
				if (checkLeft() && checkTop()) {
					done = true;
				}
			} else if (index == neighbour[1]) {
				index++;
				if (checkTop() == true) {
					done = true;
				}
			} else if (index == neighbour[2]) {
				index++;
				if (checkTop() && checkRight()) {
					done = true;

				}
			} else if (index == neighbour[3]) {
				index = neighbour[4];
				if (depth > 1) {
					if (checkFront() && checkBottom() && checkLeft()) {
						done = true;
					}
				} else {
					index = cell_size + 1;
					done = true;
					//assert(true == false);
				}
			} else if (depth > 1 && index == neighbour[4]) {
				index++;
				if (checkFront() && checkBottom()) {
					done = true;
				}
			} else if (depth > 1 && index == neighbour[5]) {
				index++;
				if (checkFront() && checkBottom() && checkRight()) {
					done = true;
				}
			} else if (depth > 1 && index == neighbour[6]) {
				index = neighbour[7];
				if (checkFront() && checkLeft()) {
					done = true;
				}
			} else if (depth > 1 && index == neighbour[7]) {
				index++;
				if (checkFront()) {
					done = true;
				}
			} else if (depth > 1 && index == neighbour[8]) {
				index++;
				if (checkFront() && checkRight()) {
					done = true;
				}
			} else if (depth > 1 && index == neighbour[9]) {
				index = neighbour[10];
				if (checkFront() && checkTop() && checkLeft()) {
					done = true;
				}
			} else if (depth > 1 && index == neighbour[10]) {
				index++;
				if (checkFront() && checkTop()) {
					done = true;
				}
			} else if (depth > 1 && index == neighbour[11]) {
				index++;
				if (checkFront() && checkTop() && checkRight()) {
					done = true;
				}
			} else if (depth > 1 && index == neighbour[12]) {
				index++;
				done = true;
				break;
			} else {
				index = original_index + width * height + width + 2;
				//std::cout << index << std::endl;
				done = true;
			}

//			std::cout << index << std::endl;
		}
		if ((depth < 2 && index >= original_index + width + 2)) {
		} else if (index >= original_index + width * height + width + 2) {

		} else if (index >= cell_size) {

		} else if (index > old_index) {
			assert(index <= cell_size);
			assert((*cells)[index]->begin() != (*cells)[index]->end());
			assert((*cells)[index]->empty() == false);
			iterator = (*cells)[index]->begin();
			assert(index > old_index);
		} else {
			assert(index == old_index);
		}
	}
}

bool LCInnerParticleIterator::operator!=(const LCInnerParticleIterator b) {

	bool return_value = false;
	if (iterator != b.iterator) {
		return_value = true;
	}
	if (index > b.index) {
		return_value = false;
	}
	return return_value;
}

LCInnerParticleIterator& LCInnerParticleIterator::operator=(
		const LCInnerParticleIterator& cpy) {
	if (this == &cpy)
		return *this;
	cells = cpy.cells;
	cell_size = cpy.cell_size;
	iterator = cpy.iterator;
	original_index = cpy.original_index;
	index = cpy.index;
	width = cpy.width;
	height = cpy.height;
	depth = cpy.depth;
	neighbour = cpy.neighbour;
	return *this;
}

int LCInnerParticleIterator::getCellNumber() {
	return index;
}

bool LCInnerParticleIterator::checkLeft() {
	return ((original_index % width) > 0);
}
bool LCInnerParticleIterator::checkRight() {
	return ((original_index % width) != (width - 1));
}
bool LCInnerParticleIterator::checkBottom() {
	return (width <= (original_index % (width * height)));
}
bool LCInnerParticleIterator::checkFront() {
	if (depth  < 2) {
		return false;
	} else {
		return (original_index  < (depth - 1) * (width * height));
	}
}
bool LCInnerParticleIterator::checkTop() {
	return (width * height - width > (original_index % (width * height)));
}

}
/* namespace utils */
