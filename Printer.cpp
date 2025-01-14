#include "Printer.h"
#include <fstream>

void savetxt(std::vector<double>& data, std::string filename) {
	std::ofstream ofs(filename);
	for (int i = 0; i < data.size(); i++) {
		ofs << data[i] << '\n';
	}
}



void savetxt(std::vector<double>& x, std::vector<double>& y, std::string filename) {
	std::ofstream ofs(filename);
	for (int i = 0; i < x.size(); i++) {
		ofs << x[i] <<' ' << y[i] << '\n';
	}
}



