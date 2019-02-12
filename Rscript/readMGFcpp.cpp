#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>

// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace Rcpp;

const double TOLERANCE = 1e-5;
const double TAGMASS = 229.162932 + 1.0072765;
const double TAGMASSTOL = TAGMASS * TOLERANCE;
const double TMT6_126 = 126.127726;
const double TMT6_127 = 127.124761;
const double TMT6_128 = 128.134436;
const double TMT6_129 = 129.131471;
const double TMT6_130 = 130.141145;
const double TMT6_131 = 131.138180;
const double TMT6_126TOL = 126.127726 * TOLERANCE;
const double TMT6_127TOL = 127.124761 * TOLERANCE;
const double TMT6_128TOL = 128.134436 * TOLERANCE;
const double TMT6_129TOL = 129.131471 * TOLERANCE;
const double TMT6_130TOL = 130.141145 * TOLERANCE;
const double TMT6_131TOL = 131.138180 * TOLERANCE;


bool scanLine(const string line);
bool rtLine(const string line);
bool massIntLine(const string line);
string getScanNum(const string line);
string getRT(const string line);
string getMZ(const string line);
string getIntensity(const string line);
bool chargeEqPlusOne(const string line);
bool isTag(const double mass);
double getMS2mass(const string line);


// [[Rcpp::export]]
Rcpp::DataFrame readMGF(const string pathfile, const string tag) {

  string line;
  ifstream mgf (pathfile.c_str());
  Rcpp::IntegerVector scanNums;
  // Rcpp::LogicalVector tagP;
  std::array<Rcpp::LogicalVector, 7> tagP;
  Rcpp::NumericVector retTime, mz, intensity;
  // for MS2 masses to find reporter ions
  double mass = 0;
  
  if (mgf.is_open()) {
    while ( getline (mgf, line) ) {
      if (scanLine(line)) {
        int snT = std::stoi(getScanNum(line));
        // Get next line and the retention time
        getline (mgf, line);
        double rtT = std::stod(getRT(line));
        // Get next line and m/z and intensity
        getline (mgf, line);
        double mzT = std::stod(getMZ(line));
        double iT = std::stod(getIntensity(line));
        // only add charge=+1 scans
        getline (mgf, line);
        if (chargeEqPlusOne(line)) {
          scanNums.push_back(snT);
          retTime.push_back(rtT);
          mz.push_back(mzT);
          intensity.push_back(iT);
          // Check MS2 for tags
          if (!tag.compare(0, 4, "tmt6")) {
            array<bool, 7> tagsFound;
            tagsFound.fill(false);
            while (getline (mgf, line) &&
                   line.compare(0, 3, "END")) {
              mass = getMS2mass(line);
              if (abs(mass - TMT6_126) < TMT6_126TOL) {
                tagsFound[0] = true;
              } else if (abs(mass - TMT6_127) < TMT6_127TOL) {
                tagsFound[1] = true;
              } else if (abs(mass - TMT6_128) < TMT6_128TOL) {
                tagsFound[2] = true;
              } else if (abs(mass - TMT6_129) < TMT6_129TOL) {
                tagsFound[3] = true;
              } else if (abs(mass - TMT6_130) < TMT6_130TOL) {
                tagsFound[4] = true;
              } else if (abs(mass - TMT6_131) < TMT6_131TOL) {
                tagsFound[5] = true;
              } else if (abs(mass - TAGMASS) < TAGMASSTOL) {
                tagsFound[6] = true;
                break;
              } else if (mass > TAGMASS)
                break;
            } // while (getline (mgf, line) && not "END")
            for (unsigned int i = 0; i < tagsFound.size(); ++i)
              tagP[i].push_back(tagsFound[i]);
          } // if (!tag.compare(0, 4, "tmt6"))
        } // if (chargeEqPlusOne(line))
      } // if (scanLine(line))
    } // while ( getline (mgf, line) )
    mgf.close();
  } // if (mgf.is_open())

  else std::cout << "Unable to open file"; 

  if (!tag.compare(0, 4, "tmt6"))    
    return Rcpp::DataFrame::create(Rcpp::Named("ScanNum")=scanNums,
                                   Rcpp::Named("RT")=retTime,
                                   Rcpp::Named("MZ")=mz,
                                   Rcpp::Named("intensity")=intensity,
                                   Rcpp::Named("tmt1")=tagP[0],
                                   Rcpp::Named("tmt2")=tagP[1],
                                   Rcpp::Named("tmt3")=tagP[2],
                                   Rcpp::Named("tmt4")=tagP[3],
                                   Rcpp::Named("tmt5")=tagP[4],
                                   Rcpp::Named("tmt6")=tagP[5],
                                   Rcpp::Named("tag")=tagP[6]);
  else
    return Rcpp::DataFrame::create(Rcpp::Named("ScanNum")=scanNums,
                                   Rcpp::Named("RT")=retTime,
                                   Rcpp::Named("MZ")=mz,
                                   Rcpp::Named("intensity")=intensity);
  //  return 0;
}

bool scanLine(const string line) {
  std::string prefix("TITLE");
  return !line.compare(0, prefix.size(), prefix);
}
bool rtLine(const string line) {
  std::string prefix("RT");
  return !line.compare(0, prefix.size(), prefix);
}
bool massIntLine(const string line)
{
  std::string prefix("PEPMASS");
  return !line.compare(0, prefix.size(), prefix);
}
string getScanNum(const string line) {
  int startOfScanNum = line.find_last_of("=")+1;
  return line.substr(startOfScanNum,
                     line.size()-(startOfScanNum+1));
}
string getRT(const string line) {
  int startOfRT = line.find_last_of("=")+1;
  return line.substr(startOfRT, line.size());
}
string getMZ(const string line) {
  int startOfMZ = line.find_first_of("=")+1;
  int endOfMZ = line.find_first_of(" ");
  return line.substr(startOfMZ, endOfMZ-startOfMZ);
}
string getIntensity(const string line) {
  auto startOfIntens = line.find_first_of(" ")+1;
  if (line.find_first_of(" ") > line.size())
    return "0";
  else 
    return line.substr(startOfIntens, line.size());
}
bool chargeEqPlusOne(const string line) {
  auto startOfCharge = line.find_last_of("=")+1;
  // either +1 or 1+
  if (startOfCharge+1 >= line.size())
    return false;
  
  if ((line.at(startOfCharge) == '1' &&
      line.at(startOfCharge+1) == '+') ||
      (line.at(startOfCharge) == '+' &&
       line.at(startOfCharge+1) == '1'))
    return true;
  else 
    return false;
}
bool isTag(const double mass) {
  return abs(mass - TAGMASS) < TAGMASSTOL;
}
double getMS2mass(const string line) {
  auto endMass = line.find_first_of(" ");
  return std::stod(line.substr(0, endMass));
}

