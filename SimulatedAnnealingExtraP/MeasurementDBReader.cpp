#include "stdafx.h"
#include "MeasurementDBReader.h"
#include <fstream>
#include <sstream>
#include <iostream>

MeasurementDB * MeasurementDBReader::readInputFile(string file)
{
	MeasurementDB* mdb = new MeasurementDB();
	std::ifstream infile(file);
	std::string line;

	if (infile.fail()) {
		std::cout << file << " does not exist. Terminating." << std::endl;
		exit(-2);
	}

	double x, y;
	char c, d, e;

	// Mode 0 = Read the training data
	// Mode 1 = Read further measurements
	int mode = 0;

	while (std::getline(infile, line))
	{
		std::istringstream iss(line);

		iss >> c;
		// Ignore comments
		if (c == '!')
			break;

		if (line == "#ADDPOINTS")
			mode = 1;

		if (line == "#STOPREADING")
			break;

		if (mode == 0)
		{
			if (!(iss >> x >> d >> y >> e)) { continue; } // ignore error
			std::pair<double, double> pair(x, y);
			mdb->addTrainingPoint(pair);
		}

		else if (mode == 1) {
			if (!(iss >> x >> d >> y >> e)) { continue; } // ignore error
			std::pair<double, double> pair(x, y);
			mdb->addMeasurementPoint(pair);
		}

	}
	return mdb;
}

MeasurementDB * MeasurementDBReader::giveExampleMeasurementDB01()
{
	MeasurementDB* mdb = new MeasurementDB();
	std::pair<double, double> pair1(0, 1);
	std::pair<double, double> pair2(1, 2);
	std::pair<double, double> pair3(2, 4);
	std::pair<double, double> pair4(3, 8);
	std::pair<double, double> pair5(4, 16);

	mdb->addTrainingPoint(pair1);
	mdb->addTrainingPoint(pair2);
	mdb->addTrainingPoint(pair3);
	mdb->addTrainingPoint(pair4);
	mdb->addTrainingPoint(pair5);

	return mdb;
}

MeasurementDB * MeasurementDBReader::giveExampleMeasurementDB02()
{
	MeasurementDB* mdb = new MeasurementDB();
	std::pair<double, double> pair1(0, 3.0);
	std::pair<double, double> pair2(1, 3.5);
	std::pair<double, double> pair3(2, 4.5);
	std::pair<double, double> pair4(3, 6.5);
	std::pair<double, double> pair5(4, 10.5);

	mdb->addTrainingPoint(pair1);
	mdb->addTrainingPoint(pair2);
	mdb->addTrainingPoint(pair3);
	mdb->addTrainingPoint(pair4);
	mdb->addTrainingPoint(pair5);

	return mdb;
}

MeasurementDB * MeasurementDBReader::giveExampleMeasurementDB03()
{
	MeasurementDB* mdb = new MeasurementDB();
	std::pair<double, double> pair1(0, 3.0);
	std::pair<double, double> pair2(1, 3.401250463);
	std::pair<double, double> pair3(2, 4.124504793);
	std::pair<double, double> pair4(3, 5.428171392);
	std::pair<double, double> pair5(4, 7.778031643);
	std::pair<double, double> pair6(4, 12.01365692);
	std::pair<double, double> pair7(4, 19.6483754);

	mdb->addTrainingPoint(pair1);
	mdb->addTrainingPoint(pair2);
	mdb->addTrainingPoint(pair3);
	mdb->addTrainingPoint(pair4);
	mdb->addTrainingPoint(pair5);

	return mdb;
}

MeasurementDB * MeasurementDBReader::giveExampleMeasurementDB04()
{
	MeasurementDB* mdb = new MeasurementDB();
	std::pair<double, double> pair1(0, 5.25);
	std::pair<double, double> pair2(1, 8.095840979);
	std::pair<double, double> pair3(2, 64.2004833);
	std::pair<double, double> pair4(3, 5612.740253);
	std::pair<double, double> pair5(4, 9128384.525);
	std::pair<double, double> pair6(3, 4.33738E+11);
	std::pair<double, double> pair7(4, 8.98308E+17);

	mdb->addTrainingPoint(pair1);
	mdb->addTrainingPoint(pair2);
	mdb->addTrainingPoint(pair3);
	mdb->addTrainingPoint(pair4);
	mdb->addTrainingPoint(pair5);
	mdb->addTrainingPoint(pair6);
	mdb->addTrainingPoint(pair7);

	return mdb;
}

MeasurementDB * MeasurementDBReader::giveExampleMeasurementDB05()
{
	MeasurementDB* mdb = new MeasurementDB();
	std::pair<double, double> pair1(40, 0.024149247310092378);
	std::pair<double, double> pair2(42, 0.025475941939138307);
	std::pair<double, double> pair3(44, 0.026775241285610468);
	std::pair<double, double> pair4(46, 0.032408996281100916);
	std::pair<double, double> pair5(48, 0.052947304775672076);
	std::pair<double, double> pair6(50, 0.08474154964755179);
	std::pair<double, double> pair7(52, 0.14185554886103377);
	std::pair<double, double> pair8(54, 0.2113750428841837);
	std::pair<double, double> pair9(56, 0.39196785122755506);
	std::pair<double, double> pair10(58, 0.6354884127322109);
	std::pair<double, double> pair11(60, 1.143539019691136);
	std::pair<double, double> pair12(62, 2.185860654746768);
	std::pair<double, double> pair13(64, 4.0220848467324055);
	std::pair<double, double> pair14(66, 6.864727794107658);
	std::pair<double, double> pair15(68, 12.094192494629402);
	std::pair<double, double> pair16(70, 20.391819644117003);
	std::pair<double, double> pair17(72, 34.94810368046828);
	std::pair<double, double> pair18(74, 68.81139729593603);
	std::pair<double, double> pair19(76, 128.55942520311189);
	std::pair<double, double> pair20(78, 249.09854780212135);
	std::pair<double, double> pair21(80, 289.3069826699564);
	std::pair<double, double> pair22(82, 589.8312887869633);
	std::pair<double, double> pair23(84, 1074.5646286076976);
	std::pair<double, double> pair24(86, 1971.1372180479016);
	std::pair<double, double> pair25(88, 3279.4092065538193);
	std::pair<double, double> pair26(90, 5901.874941082804);

	mdb->addTrainingPoint(pair1);
	mdb->addTrainingPoint(pair2);
	mdb->addTrainingPoint(pair3);
	mdb->addTrainingPoint(pair4);
	mdb->addTrainingPoint(pair5);
	mdb->addTrainingPoint(pair6);
	mdb->addTrainingPoint(pair7);
	mdb->addTrainingPoint(pair8);
	mdb->addTrainingPoint(pair9);
	mdb->addTrainingPoint(pair10);
	mdb->addTrainingPoint(pair11);
	mdb->addTrainingPoint(pair12);
	mdb->addTrainingPoint(pair13);
	mdb->addTrainingPoint(pair14);
	mdb->addTrainingPoint(pair15);
	mdb->addTrainingPoint(pair16);
	mdb->addTrainingPoint(pair17);
	mdb->addTrainingPoint(pair18);
	mdb->addTrainingPoint(pair19);
	mdb->addTrainingPoint(pair20);
	mdb->addTrainingPoint(pair21);
	mdb->addTrainingPoint(pair22);
	mdb->addTrainingPoint(pair23);
	mdb->addTrainingPoint(pair24);
	mdb->addTrainingPoint(pair25);
	mdb->addTrainingPoint(pair26);

	return mdb;
}
