/*
 * LCOuterParticleIterator.cpp
 *
 *  Created on: Nov 28, 2013
 *      Author: andal
 */

#include "LCOuterParticleIterator.h"
#include <log4cxx/logger.h>
#include <omp.h>

log4cxx::LoggerPtr lcouterparticleiteratorlogger(
		log4cxx::Logger::getLogger("utils.lcouterparticleiterator"));

namespace utils {

LCOuterParticleIterator::LCOuterParticleIterator() {

}

LCOuterParticleIterator::LCOuterParticleIterator(int cell_size_arg,
		std::vector<std::list<Particle *> *>* cells_arg,
		std::list<Particle *>::iterator iterator_arg, int index_arg) :
		cells(cells_arg) {
	cell_size = cell_size_arg;
	index = index_arg;
	iterator = iterator_arg;
	cellnumber = 0;
}

LCOuterParticleIterator::LCOuterParticleIterator(int cell_size_arg,
		std::vector<std::list<Particle *> *>* cells_arg,
		std::vector<int>* cellnumber_arg,
		std::list<Particle *>::iterator iterator_arg, int index_arg) :
		cells(cells_arg), cellnumber(cellnumber_arg) {
	cell_size = cell_size_arg;
	index = index_arg;
	iterator = iterator_arg;
	assert((*cellnumber).size() > 0);
}

LCOuterParticleIterator::~LCOuterParticleIterator() {
	// TODO Auto-generated destructor stub
}

Particle& LCOuterParticleIterator::operator*() const {
	return (*(*iterator));
}

bool LCOuterParticleIterator::operator!=(const LCOuterParticleIterator b) {
	bool return_value = false;
	if (iterator != b.iterator) {
		return_value = true;
	}
	if (index > b.index) {
		return_value = false;
	}
	if (cellnumber != 0) {
		if ((*cellnumber)[index] > (*b.cellnumber)[b.index]) {
			return_value = false;
		}

	}
//	if(cellnumber != 0 && return_value == true){
//#pragma omp critical
//		{
//			std::cout << "Thread: " << omp_get_thread_num() << std::endl;
//			std::cout << "Index: " << index << std::endl;
//			std::cout << "Cell Index: " << (*cellnumber)[index] << std::endl;
//		}
//	}
	return return_value;
}

void LCOuterParticleIterator::operator++() {
	if (cellnumber == 0) {
		assert((*cells)[index]->empty() == false);
		//std::cout << (*cells)[index].size() << " " << index << " " << (*(*iterator)).toString() << std::endl;

		++iterator;
		if (iterator != (*cells)[index]->end()) {
			assert(iterator != (*cells)[index]->end());
		} else {
			int old_index = index;
			index++;
			while (index < cell_size && (*cells)[index]->empty() == true) {
				index++;
			}
			LOG4CXX_TRACE(lcouterparticleiteratorlogger, index);
			assert(index > old_index);
			if (index < cell_size) {
				iterator = (*cells)[index]->begin();
			}
		}
	} else {

		assert((*cells)[(*cellnumber)[index]]->empty() == false);

		++iterator;
		if (iterator != (*cells)[(*cellnumber)[index]]->end()) {
			assert(iterator != (*cells)[(*cellnumber)[index]]->end());
		} else {
			int old_index = index;
			index++;
			while (index < cellnumber->size()
					&& (*cells)[(*cellnumber)[index]]->empty() == true) {
				index++;
			}
			LOG4CXX_TRACE(lcouterparticleiteratorlogger, (*cellnumber)[index]);

			assert(old_index < index);

			if (index < cellnumber->size()) {
				assert((*cells)[(*cellnumber)[index]]->empty() == false);
				iterator = (*cells)[(*cellnumber)[index]]->begin();
			}
		}
	}
}

int LCOuterParticleIterator::getCellNumber() {
	if (cellnumber == 0) {
		return index;
	} else {
		return (*cellnumber)[index];
	}
}

std::list<Particle *>::iterator LCOuterParticleIterator::getIterator() {
	std::list<Particle *>::iterator newIterator = iterator;
	return newIterator;
}

LCOuterParticleIterator& LCOuterParticleIterator::operator=(
		const LCOuterParticleIterator& cpy) {
	if (this == &cpy)
		return *this;

	cells = cpy.cells;
	cell_size = cpy.cell_size;
	index = cpy.index;
	iterator = cpy.iterator;
	cellnumber = cpy.cellnumber;
	return *this;
}

} /* namespace utils */
