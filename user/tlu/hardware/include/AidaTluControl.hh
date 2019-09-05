#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <csignal>
#include <memory>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <ctime>
#include <thread>
#include <map>
#include <math.h>
#include <numeric>
#include "AidaTluController.hh"
#include "AidaTluHardware.hh"
#include "AidaTluPowerModule.hh"
#include "eudaq/OptionParser.hh"


#include <TROOT.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TBrowser.h>
#include <TFrame.h>
#include <TFile.h>
#include <TApplication.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TGraph2D.h>
#include <TMath.h>

class AidaTluControl {
public:
    AidaTluControl();
    void DoStartUp();
    void SetPMTVoltage(double voltage);
    void SetPMTVoltage(std::vector<double> voltage);
    void SetTLUThreshold(double threshold);
    void SetTLUThreshold(std::vector<double> threshold, std::vector<bool> connection, std::string mode);
    std::vector<std::vector<std::vector<double>>> readFiles(std::string filename);
    std::vector<std::vector<double> > GetOptimalThreshold(std::string filename);
    void PlotTrigger(std::string filename);
    std::vector<double> MeasureRate(std::vector<bool> connection);
    void WriteOutputFile(std::string filename, double voltage, std::vector<std::vector<double>> rates, std::vector<double> thresholds);
    void WriteOutputFile(std::string filename, std::vector<double> voltage, std::vector<std::vector<double>> rates, std::vector<double> thresholds);
    void WriteOutputFileTrigger(int channel, std::string filename, std::vector<double> voltage, std::vector<std::vector<double>> thresholdsTrigger, std::vector<std::vector<std::vector<double>>> ratesTrigger, std::vector<double> optimalThresholds);
    void WriteParameters(int in_numThresholdValues, int in_numTriggerInputs, int in_time, int in_numStepsTrigger);
    void WriteParameters(int in_numThresholdValues, int in_numTriggerInputs, int in_time);
    bool flagPlateauError;

private:
    std::unique_ptr<tlu::AidaTluController> m_tlu;
    uint8_t m_verbose;
    uint64_t m_starttime;
    uint64_t m_lasttime;
    int m_nTrgIn;
    int numThresholdValues;
    int numTriggerInputs;
    int time;
    int numStepsTrigger;

    double m_duration;
    //    bool m_exit_of_run;
};

