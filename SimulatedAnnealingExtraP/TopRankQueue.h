#ifndef TOPRANKQUEUE_H
#define TOPRANKQUEUE_H

#include "AbstractSolution.h"
#include <vector>
#include<map>
#include <limits>
#include <iomanip>

using namespace std;

namespace DarmstadtSC {
	template <class SolutionType>
	class TopRankQueue {
	public:
		TopRankQueue(int size) : _size(size) {
			// Initialize map
			SolutionType sol = SolutionType();
			std::pair<double, SolutionType> solpair(numeric_limits<double>::max(), sol);

			//auto ret_it = map<double, SolutionType>::iterator it = _entries.begin();
			/*for (int i = 0; i < _size; i++) {
				_entries.insert(solpair);
				solpair.first -= 1e-11;
			}*/

			_no_entries = 0;
			_highest_cost = numeric_limits<double>::max();
		}

		~TopRankQueue() {
			_entries.clear();
		}

		int insert(SolutionType sol) {
			if (sol.get_costs() > _highest_cost)
				return -1;

			if (_no_entries < _size) {
				_no_entries++;
			}

			int prev_size = _entries.size();

			//SolutionType sol = in_sol;

			// Remove the highest entry in the map
			//map<double, SolutionType> it = _entries.rbegin();
			//_entries.erase(it);

			// Check if solution is already there
			typename map<double, SolutionType>::iterator search_it;
			search_it = _entries.find(sol.get_costs());

			if (search_it != _entries.end()) 
			{
				cout << "Already there" << endl;
				cout << setprecision (15) << "[" << sol.get_costs() << "] : " << endl;
				sol.printModelFunction();
				cout << "VS" << endl;
				_entries[sol.get_costs()].printModelFunction();
				cout << endl << endl;
			}

			// Add the new solution
			_entries[sol.get_costs()] = sol;
			int after_size = _entries.size();

			if (prev_size == after_size) {

			}

			else {
				if (_no_entries == _size) {
					_entries.erase(prev(_entries.end()));
					_highest_cost = sol.get_costs();
				}
			}

			return _no_entries;
		}

		bool contains(SolutionType& sol) {
			if (_entries.find (sol->sol.get_costs()) == _entries.end())
					return false;
			return true;
		}

		void printQueue() {
			int i = 0;
			for (auto it = _entries.begin(); it != _entries.end(); it++) {
				std::cout << "[" << i << "]:"  << it->first << " ";
				it->second.printModelFunction();
				cout << endl;
				i++;
			}
			std::cout << std::endl;
		}

	private:
		//std::vector<SolutionType> _entries;
		int _size;
		int _no_entries;
		int _nextpos;
		double _highest_cost;

		std::map<double, SolutionType> _entries;
	};
}

#endif