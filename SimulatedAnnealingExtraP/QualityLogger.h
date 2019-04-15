#ifndef QUALITYLOGGER_H
#define QUALITYLOGGER_H

#include <vector>
#include <map>

typedef std::vector<std::pair<unsigned int, double>> LogVector;
typedef std::map<int, LogVector> LogVectorMap;

class QualityLogger {
public:
	static QualityLogger& getInstance()
	{
		static QualityLogger instance;
		return instance;
	}

	int get_size() { return get_size(0); }
	int get_size(int tid);

	std::pair<unsigned int, double> get_entry_at(unsigned int pos) { return get_entry_at(pos, 0); }
	std::pair<unsigned int, double> get_entry_at(unsigned int pos, int tid);

	void insertEntry(std::pair<unsigned int, double>&);
	void insertEntry(std::pair<unsigned int, double>&, int tid);

	double get_max_cost() { return this->max_cost; }

	static QualityLogger* _instance;

	int get_no_of_vectors() { return (int)this->logbookmap.size(); }

	double get_max_cost_to_print(int stepsize);
	double get_max_cost_to_print(int stepsize, int tid);

private:	
	QualityLogger();
	QualityLogger(const QualityLogger&);
	~QualityLogger() {}

	LogVectorMap logbookmap;
	//LogVector logbook;

	double max_cost;
	int tn;
};

#endif
