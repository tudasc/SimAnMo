#include "stdafx.h"
#include "QualityLogger.h"

QualityLogger::QualityLogger() {
	this->max_cost = -1.0;
	tn = 1;
}

int QualityLogger::get_size(int tid)
{
	const LogVector& vec = logbookmap[tid];
	return vec.size();
}

std::pair<unsigned int, double> QualityLogger::get_entry_at(unsigned int pos, int tid)
{
	const LogVector& vec = logbookmap[tid];
	return vec[pos];
}

void QualityLogger::insertEntry(std::pair<unsigned int, double>& npair)
{
	insertEntry(npair, 0);
}

void QualityLogger::insertEntry(std::pair<unsigned int, double>& npair, int tid)
{
//#pragma omp critical 
	{
		LogVector& vec = logbookmap[tid];
		vec.push_back(npair);

		if (npair.second > this->max_cost)
			this->max_cost = npair.second;
	}
}

double QualityLogger::get_max_cost_to_print(int stepsize)
{
	double maxcost = -1.0;
	// Run over all vectors in the map to find the maximum
	for (LogVectorMap::const_iterator it = logbookmap.begin();
		it != logbookmap.end(); ++it) {
		// Search in all entries of the corresponding vector
		const LogVector& vec = it->second;
		for (int i = 0; i < vec.size(); i += stepsize) {
			std::pair<unsigned int, double> act_pair = get_entry_at(i, it->first);
			if (act_pair.second > maxcost)
				maxcost = act_pair.second;
		}
	}
	return maxcost;
}

double QualityLogger::get_max_cost_to_print(int stepsize, int tid)
{
	double maxcost = -1.0;
	// Run over all vectors in the map to find the maximum

		// Search in all entries of the corresponding vector
		const LogVector& vec = this->logbookmap[tid];

		for (int i = 0; i < vec.size(); i += stepsize) {
			std::pair<unsigned int, double> act_pair = vec[i];
			if (act_pair.second > maxcost)
				maxcost = act_pair.second;
		}
	return maxcost;
}
