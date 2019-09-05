#include "AidaTluControl.hh"
#include "AidaTluController.hh"
#include "AidaTluHardware.hh"
#include "AidaTluPowerModule.hh"
#include "eudaq/OptionParser.hh"

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

// ROOT includes
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

int main(int /*argc*/, char **argv) {
    eudaq::OptionParser op("EUDAQ Command Line FileReader modified for TLU", "2.1", "EUDAQ FileReader (TLU)");
    eudaq::Option<double> thrMin(op, "tl", "thresholdlow", -0.2, "double", "threshold value low [V]");
    eudaq::Option<double> thrMax(op, "th", "thresholdhigh", -0.01, "double", "threshold value high [V]");
    eudaq::Option<int> thrNum1(op, "tn1", "thresholdsteps 1", 10, "int", "number of threshold steps");
    eudaq::Option<int> thrNum2(op, "tn2", "thresholdsteps 2", 30, "int", "number of threshold steps");
    eudaq::Option<double> volt1(op, "v1", "pmtvoltage 1", 0.85, "double", "PMT voltage [V]");
    eudaq::Option<double> volt2(op, "v2", "pmtvoltage 2", 0.85, "double", "PMT voltage [V]");
    eudaq::Option<double> volt3(op, "v3", "pmtvoltage 3", 0.85, "double", "PMT voltage [V]");
    eudaq::Option<double> volt4(op, "v4", "pmtvoltage 4", 0.85, "double", "PMT voltage [V]");
    eudaq::Option<int> acqtime1(op, "t1", "acquisitiontime 1", 10, "int", "acquisition time");
    eudaq::Option<int> acqtime2(op, "t2", "acquisitiontime 2", 10, "int", "acquisition time");
    eudaq::Option<std::string> name(op, "f", "filename", "output", "string", "filename");
    eudaq::Option<std::string> con(op, "c", "connectionmap", "111100", "string", "connection map");

    try{
        op.Parse(argv);
    }
    catch (...) {
        return op.HandleMainException();
    }

    // pass all values:
    double thresholdMin = thrMin.Value();
    double thresholdMax = thrMax.Value();
    int numThresholdValues1 = thrNum1.Value();
    int numThresholdValues2 = thrNum2.Value();
    std::vector<double> voltage = {volt1.Value(), volt2.Value(), volt3.Value(), volt4.Value()};
    int time1 = acqtime1.Value(); //time in seconds
    int time2 = acqtime2.Value(); //time in seconds
    std::string filename = name.Value();
    std::string connection = con.Value();

    double standardThreshold = -0.07;
    double standardVoltage = 0.85;

    if (filename == "output") {
        std::cout << "---------------CAUTION: FILENAME IS SET TO DEFAULT. DANGER OF DATA LOSS!---------------" <<std::endl;
        std::cout << "Enter 'c' to continue" << std::endl;
        int k;
        std::cin >> k;
        std::cin.clear();
        std::cin.ignore(512, '\n');
    }

    double thresholdDifference = thresholdMax - thresholdMin;

    int numTriggerInputs = 0;
    std::vector<bool> connectionBool(6, false);
    for (int i = 0; i < 6; i++){
        if(connection[i] == '1'){
            numTriggerInputs++;
            connectionBool[i] = true;
        }
    }

    if(numTriggerInputs < 2){
        std::cout << "CAUTION: number of TLU inputs is smaller than 2! Correction for beam fluctuation is not possible."<< std::endl;
        std::cout << "Enter 'c' to continue" << std::endl;
        int k;
        std::cin >> k;
        std::cin.clear();
        std::cin.ignore(512, '\n');
    }

    int numThresholdValues;
    int time;
    AidaTluControl myTlu;

    ///////////////////////////////////////////////
    // Pre Run:
    ///////////////////////////////////////////////
    int preRun = 1;
    int iterations = 1;
    while (preRun == 1){
        numThresholdValues = numThresholdValues1;
        time = time1;        
        // create array of threshold
        std::vector<double> thresholds(numThresholdValues);

        if (numThresholdValues < 2) thresholds[0] = thresholdMin;
        else{
            for (int i = 0; i < numThresholdValues; i++){
                thresholds[i] = thresholdMin + i * thresholdDifference / (numThresholdValues-1);
            }
        }

        myTlu.WriteParameters(numThresholdValues, numTriggerInputs, time);

        std::vector<std::vector<double>> rates(numThresholdValues, std::vector<double>(numTriggerInputs + 2));



        myTlu.DoStartUp();
        for (int i = 0; i < numThresholdValues; i++){
            myTlu.SetPMTVoltage({standardVoltage, voltage[1], voltage[2], voltage[3]});
            myTlu.SetTLUThreshold(thresholds[i]);
            myTlu.SetTLUThreshold({standardThreshold}, connectionBool, "first");

            rates[i] = myTlu.MeasureRate(connectionBool);
        }
        std::string filenameFirst = filename + "_" + std::to_string(iterations) + std::string("_first");
        myTlu.WriteOutputFile(filenameFirst, voltage, rates, thresholds);

        // Repeat Measurement for first input, now the second input is constant
        for (int i = 0; i < numThresholdValues; i++){
            myTlu.SetPMTVoltage({voltage[0], standardVoltage,voltage[2], voltage[3]});
            myTlu.SetTLUThreshold(thresholds[i]);
            myTlu.SetTLUThreshold({standardThreshold}, connectionBool, "second");

            rates[i] = myTlu.MeasureRate(connectionBool);
        }
        std::string filenameSecond = filename + "_" + std::to_string(iterations) + std::string("_second");
        myTlu.WriteOutputFile(filenameSecond, voltage, rates, thresholds);



        std::vector<std::vector<double>> optimalReturn = myTlu.GetOptimalThreshold(filename+ "_" + std::to_string(iterations));

        std::cout << "Do you want to change the voltages and rerun the measurement?" << std::endl;
        int decision;
        std::cin >> decision;
        std::cin.clear();
        std::cin.ignore(512, '\n');
        if (decision > 0){
            double v0,v1,v2,v3;
            std::cout << "Enter v1" << std::endl;
            std::cin >> v0;
            std::cin.clear();
            std::cin.ignore(512, '\n');
            voltage[0] = v0;
            std::cout << "Enter v2" << std::endl;
            std::cin >> v1;
            std::cin.clear();
            std::cin.ignore(512, '\n');
            voltage[1] = v1;
            std::cout << "Enter v3" << std::endl;
            std::cin >> v2;
            std::cin.clear();
            std::cin.ignore(512, '\n');
            voltage[2] = v2;
            std::cout << "Enter v4" << std::endl;
            std::cin >> v3;
            std::cin.clear();
            std::cin.ignore(512, '\n');
            voltage[3] = v3;
            iterations++;
        }
        else preRun = 0;
    }

    ///////////////////////////////////////////////
    // Main Run:
    ///////////////////////////////////////////////

    numThresholdValues = numThresholdValues2;
    time = time2;
    // create array of threshold
    std::vector<double> thresholds(numThresholdValues);

    if (numThresholdValues < 2) thresholds[0] = thresholdMin;
    else{
        for (int i = 0; i < numThresholdValues; i++){
            thresholds[i] = thresholdMin + i * thresholdDifference / (numThresholdValues-1);
        }
    }

    myTlu.WriteParameters(numThresholdValues, numTriggerInputs, time);

    std::vector<std::vector<double>> rates(numThresholdValues, std::vector<double>(numTriggerInputs + 2));


    myTlu.DoStartUp();
    for (int i = 0; i < numThresholdValues; i++){
        myTlu.SetPMTVoltage({standardVoltage, voltage[1], voltage[2], voltage[3]});
        myTlu.SetTLUThreshold(thresholds[i]);
        myTlu.SetTLUThreshold({standardThreshold}, connectionBool, "first");

        rates[i] = myTlu.MeasureRate(connectionBool);
    }
    std::string filenameFirst = filename + std::string("_first");
    myTlu.WriteOutputFile(filenameFirst, voltage, rates, thresholds);

    // Repeat Measurement for first input, now the second input is constant
    for (int i = 0; i < numThresholdValues; i++){
        myTlu.SetPMTVoltage({voltage[0], standardVoltage,voltage[2], voltage[3]});
        myTlu.SetTLUThreshold(thresholds[i]);
        myTlu.SetTLUThreshold({standardThreshold}, connectionBool, "second");

        rates[i] = myTlu.MeasureRate(connectionBool);
    }
    std::string filenameSecond = filename + std::string("_second");
    myTlu.WriteOutputFile(filenameSecond, voltage, rates, thresholds);

    std::vector<std::vector<double>> optimalReturn = myTlu.GetOptimalThreshold(filename);
    std::vector<double> optimalThresholds;
    std::vector<double> thresholdMinOpt;
    std::vector<double> thresholdMaxOpt;

    optimalThresholds = optimalReturn[0];
    thresholdMinOpt = optimalReturn[1];
    thresholdMaxOpt = optimalReturn[2];

    std::cout << "__________________________" << std::endl;
    std::cout << "Optimal Threshold Values:" << std::endl;
    for (int i = 0; i < numTriggerInputs; i++){
        std::cout << "PMT " << i+1 << ":   " << optimalThresholds[i]<< " V" << "   (" << "Plateau:  "<< thresholdMinOpt[i] << " - " << thresholdMaxOpt[i] << ")"<< std::endl;
    }
    std::cout << "__________________________" << std::endl;

    return 1;

}

