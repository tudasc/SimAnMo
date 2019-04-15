#ifndef RINGQUEUE_H
#define RINGQUEUE_H

#include <iostream>

namespace MichaelB {
	template <class T>
	class RingQueue {
	public:
		RingQueue(int size) : _size(size) {
			_entries = new T[_size];
			_no_entries = 0;
			_nextpos = 0;
		}

		~RingQueue() {
			delete[] _entries;
		}

		int insert(T val) {
			if (_no_entries < _size) {
				_no_entries++;
			}

			_entries[_nextpos] = val;
			_nextpos = (_nextpos + 1) % _size;

			return _nextpos;
		}

		bool contains(T val) {
			for (int i = 0; i < _no_entries; i++)
				if (_entries[i] == val)
					return true;
			return false;
		}

		void printQueue() {
			for (int i = 0; i < _no_entries; i++) {
				std::cout << "[" << i << "]:" << _entries[i] << " ";
			}
			std::cout << std::endl;
		}

	private:
		T * _entries;
		int _size;
		int _no_entries;
		int _nextpos;
	};
}

#endif // !RINGQUEUE_H



