// ExpontialRuntimeTester.cpp : Diese Datei enthält die Funktion "main". Hier beginnt und endet die Ausführung des Programms.
//

#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <chrono>
#include <cmath>

/*#include < filesystem>
#include <fstream>
#include <boost/algorithm/string/classification.hpp>// Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp>
#include <map>
#include <vector>*/

#define RARRAY_SIZE 5000

using namespace std;
using namespace chrono;

/*namespace fs = std::filesystem;

map<int, vector<long long>> _miterations;

int readResFile(std::string filestr, long long& nodes, int& node_entries)
{
	std::string locpath = filestr;
	std::ifstream file(locpath.c_str());

	if (!file)
	{
		cerr << "No valid Resfile " << locpath << endl;
		return -1;
	}

	std::string   line;
	bool isnewdim = false;
	bool isnewpoint = false;
	int actdim = 0;
	int actdimcnt = 0;

	double timesperit[5];

	string str1 = "Running";
	string str2 = "The";

	string str3 = "[Thread:";
	string str4 = "ENUM:";

	while (std::getline(file, line))
	{

		std::vector<std::string> words;
		boost::split(words, line, boost::is_any_of(", "), boost::token_compress_on);

		if (words.size() == 0)	continue;

		if (str1.compare(words[0]) == 0) {
			int dim = (int)std::stoi(words[5]);
			long long its = (long long)std::stoi(words[3]);

			auto it = _miterations.find(dim);
			if (it != _miterations.end()) {
				_miterations[dim].push_back(its);
			}
			else {
				vector<long long> tvec;
				tvec.push_back(its);
				_miterations[dim] = tvec;
			}
		}
			

		if (str1.compare(words[0]) == 0 && actdim != (int)std::stoi(words[5])) {
			isnewdim = true;
			
			try {
				if (actdim > 0 && actdimcnt < 5) {
					cerr << "Inkonsistent input file. Ending." << endl;
					return -1;
				}

				if (actdim > 0 && actdimcnt == 5) {
					double sum = 0.0;
					for (int i = 0; i < 5; i++) {
						sum += timesperit[i];
					}
					double average = sum / 5.0;
					cout << "( " << actdim << " ; " << average << " )" << endl;
				}
				actdim = (int)std::stoi(words[5]);
				actdimcnt = 0;
				continue;

			}
			catch (const std::exception& e) {
				cout << words[5] << " is no number." << endl;
				cout << e.what() << endl;
			}
		}

		if (str2.compare(words[0]) == 0 && actdimcnt < 5) {
			double time = 0.0;
			try {
				time = (double)std::stof(words[5]);
			}
			catch (const std::exception& e) {
				cout << words[5] << " is no floating number." << endl;
				cout << e.what() << endl;
			}
			timesperit[actdimcnt] = time;
			actdimcnt++;
		}

	}

	// Print the iterations
	for (auto it = _miterations.begin(); it != _miterations.end(); ++it) {
		for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
			cout << "( " << it->first << " ; " << *it2 << " )" << endl;
		}
	}

	return 0;
}

int readRunningTimes(int argc, char** argv) {
	long long nodes = 0;
	int node_entries = 0;
	std::string filepath;

	bool do_result_reading = false;

	for (int i = 1; i < argc; i++) {
		std::string input = std::string(argv[i]);

		if (input == "--resultDir") {
			do_result_reading = true;

			if (argc <= i) {
				std::cerr << "Missing basisfile in argument. Terminating." << std::endl;
				return 1;
			}

			filepath = std::string(argv[i + 1]);
			break;
		}
	}

	if (!do_result_reading) {
		return 0;
	}

	cout << "Got filepath: " << filepath << endl;


	for (auto& p : fs::directory_iterator(filepath.c_str())) {
		std::cout << p.path() << '\n';
		readResFile(p.path().u8string(), nodes, node_entries);
	}

	//cout << "Nodes in average: " << nodes / node_entries << endl;
	exit(4);
	return 0;
	
}*/


int main(int argc, char** argv)
{
	//readRunningTimes(argc, argv);

    srand((unsigned int)time(NULL));
    long long probsize = 0;
    if (argc > 1) {
        probsize = atoi(argv[1]);
    }
    else {
        // Read a probsize from user
        cout << "Starting problem size?";
        cin >> probsize;
    }
    

    // Create and fill the Random numbers array
    double** rarray = new double* [RARRAY_SIZE];
    for (int i = 0; i < RARRAY_SIZE; i++) {
        rarray[i] = new double[RARRAY_SIZE];
    }
    
    for (int i = 0; i < RARRAY_SIZE-1; i++) {
        for (int j = 0; j < RARRAY_SIZE-1; j++) {
            rarray[i][j] = (double)rand() / 10034.346321;
            //cout << rarray[i][j] << endl;
        }
    }   

    while (probsize < 570)
    {
        for (int k = 0; k < 5; k++) {
            double res = 0.0;
            // Calculate the number of iterations that will be performed
            //RUN1 
			//long long its = (long long)pow(2.0, 0.32 * (double)probsize + 2.25);
			//RUN2		
			long long its = (long long)pow(2.0, 1.49 * sqrt((double)probsize) + 2.25);
			// RUN3			
			//long long its = (long long)pow(2.0, 0.11 * pow((double)probsize,1.71) + 1.12);
			// RUN4
			//long long its = log(probsize) * (long long)pow(2.0, 0.32 * (double)probsize);

            double jitter = (double)((int)rand() % 5001 - 2500) / 100000;
            
            its -= (long long)(its * jitter);
            cout << "Running with iterarations: " << its << " (Probsize: " << probsize << " )" << endl;
            high_resolution_clock::time_point t1 = high_resolution_clock::now();

			if (its < 1)
				exit(5);

            for (long long i = 0; i < its; i++) {
                for (long long j = 0; j < 1000000u; j++) {}
                int x = rand() % RARRAY_SIZE;
                int y = rand() % RARRAY_SIZE;
                int fac = (rand() % 3) - 1;
                res += fac * sqrt(fac * rarray[x][y] * fac * rarray[x][y] * 0.9999999999);
				
				x = rand() % RARRAY_SIZE;
                y = rand() % RARRAY_SIZE;
                fac = (rand() % 3) - 1;
				
				// Something < 1
				double betw = cos(pow(sin(rarray[x][y]),50));				
                res += fac * sqrt(fac * rarray[x][y] * fac * rarray[x][y] * 0.9999999999)
					+ fac * sqrt(fac * rarray[x][y] * fac * rarray[x][y] * 0.9999999999) * betw;
				
				if(its % 2 == 0)
				{
					x = rand() % RARRAY_SIZE;
					y = rand() % RARRAY_SIZE;
					fac = (rand() % 3) - 1;
					res += fac * sqrt(fac * rarray[x][y] * fac * rarray[x][y] * 0.9999999999);
				}
				
				else {
					x = rand() % (RARRAY_SIZE-1);
					y = rand() % (RARRAY_SIZE-1);
					fac = (rand() % 3) - 1;
					res += fac * sqrt(fac * rarray[x][y] * fac * rarray[x][y] * 0.99999);
				}
				
				x = rand() % RARRAY_SIZE;
                y = rand() % RARRAY_SIZE;
                fac = (rand() % 3) - 1;
                res += fac * sqrt(fac * rarray[x][y] * fac * rarray[x][y] * 0.9999999999);
				
				x = rand() % RARRAY_SIZE;
                y = rand() % RARRAY_SIZE;
                fac = (rand() % 3) - 1;
                res += fac * sqrt(fac * rarray[x][y] * fac * rarray[x][y] * 0.9999999999);
            }
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
            cout << "The result is: " << res << " after " << time_span.count() << " secs." << endl;
            
        }
        probsize++;
    }
    return 0;
}

// Programm ausführen: STRG+F5 oder Menüeintrag "Debuggen" > "Starten ohne Debuggen starten"
// Programm debuggen: F5 oder "Debuggen" > Menü "Debuggen starten"

// Tipps für den Einstieg: 
//   1. Verwenden Sie das Projektmappen-Explorer-Fenster zum Hinzufügen/Verwalten von Dateien.
//   2. Verwenden Sie das Team Explorer-Fenster zum Herstellen einer Verbindung mit der Quellcodeverwaltung.
//   3. Verwenden Sie das Ausgabefenster, um die Buildausgabe und andere Nachrichten anzuzeigen.
//   4. Verwenden Sie das Fenster "Fehlerliste", um Fehler anzuzeigen.
//   5. Wechseln Sie zu "Projekt" > "Neues Element hinzufügen", um neue Codedateien zu erstellen, bzw. zu "Projekt" > "Vorhandenes Element hinzufügen", um dem Projekt vorhandene Codedateien hinzuzufügen.
//   6. Um dieses Projekt später erneut zu öffnen, wechseln Sie zu "Datei" > "Öffnen" > "Projekt", und wählen Sie die SLN-Datei aus.
