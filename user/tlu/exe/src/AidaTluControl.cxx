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
    eudaq::Option<double> thrMin(op, "tl", "thresholdlow", -0.21, "double", "threshold value low [V]");
    eudaq::Option<double> thrMax(op, "th", "thresholdhigh", -0.01, "double", "threshold value high [V]");
    eudaq::Option<int> thrNum(op, "tn", "thresholdsteps", 30, "int", "number of threshold steps");
    eudaq::Option<double> volt(op, "v", "pmtvoltage", 0.8, "double", "PMT voltage [V]");
    eudaq::Option<int> acqtime(op, "t", "acquisitiontime", 10, "int", "acquisition time");
    eudaq::Option<std::string> name(op, "f", "filename", "output", "string", "filename");
    eudaq::Option<std::string> con(op, "c", "connectionmap", "111100", "string", "connection map");
    eudaq::Option<bool> tluOn(op, "tlu", "tluconnected", true, "bool", "tlu connected");
    eudaq::Option<bool> optTrig(op, "opt", "optimizetrig", false, "bool", "optimize for trigger rates");
    eudaq::Option<bool> optPlot(op, "optplot", "optplot", false, "bool", "only plot optimization stuff");
    eudaq::Option<int> trigNum(op, "otn", "optimizetrignum", 11, "int", "number of steps for trigger optimization");
    eudaq::Option<int> deb(op, "deb", "debug", 0, "int", "debug number");


    try{
        op.Parse(argv);
    }
    catch (...) {
        return op.HandleMainException();
    }

    // Threshold in [-1.3V,=1.3V] with 40e-6V presision
    int numStepsTrigger = trigNum.Value();

    bool optimizationPlot = optPlot.Value();

    int debugNumber = deb.Value();
    double thresholdMin = thrMin.Value();
    double thresholdMax = thrMax.Value();
    double thresholdDifference = thresholdMax - thresholdMin;
    int numThresholdValues = thrNum.Value();
    int time = acqtime.Value(); //time in seconds
    double voltage = volt.Value();
    std::string filename = name.Value();
    if (filename == "output") {
        std::cout << "---------------CAUTION: FILENAME IS SET TO DEFAULT. DANGER OF DATA LOSS!---------------" <<std::endl;
        std::cout << "Enter 'c' to continue" << std::endl;
        int k;
        std::cin >> k;
    }
    std::string connection = con.Value();

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
    }


    // create array of thresholds
    std::vector<double> thresholds(numThresholdValues);

    if (numThresholdValues < 2) thresholds[0] = thresholdMin;
    else{
        for (int i = 0; i < numThresholdValues; i++){
            thresholds[i] = thresholdMin + i * thresholdDifference / (numThresholdValues-1);
        }
    }

    AidaTluControl myTlu;
    myTlu.WriteParameters(numThresholdValues, numTriggerInputs, time, numStepsTrigger);

    std::vector<std::vector<double>> rates(numThresholdValues, std::vector<double>(numTriggerInputs + 2));
    bool tluConnected = tluOn.Value();
    bool optimizeTriggerRates = optTrig.Value();


    // Define standard parameters for reference channel:
    //double standardVoltage = 0.8;
    //std::vector<double> standardThreshold = {-0.08};


    /////////////////////////////////////////////////////////////////
    // Determine optimal Threshold values
    /////////////////////////////////////////////////////////////////
//    std::vector<double> optimalVoltages;
//    std::vector<double> optimalThresholds;
//    optimalVoltages = {0.9, 0.95, 0.9, 0.85};
//    optimalThresholds = {-0.064, -0.078, -0.078, -0.068};
    double standardThreshold = -0.07;
    double standardVoltage = 0.85;

    if(tluConnected){
        //std::vector<double> voltages = {0.8, 0.85, 0.9, 0.95};
        //std::vector<std::string> filenames = {};
        //std::cout << "CAUTION!! CHANGE FILENAMES!!"<< std::endl;


//        for (int k = 0; k<4; k++){
//            voltage = voltages[k];
//            filename = filenames[k];
            // background measurements
            if (debugNumber == 2){
                std::cout << "Measuring Background"<< std::endl;
                myTlu.DoStartUp();
                for (int i = 0; i < numThresholdValues; i++){
                    myTlu.SetPMTVoltage(voltage);
                    myTlu.SetTLUThreshold(thresholds[i]);
                    rates[i] = myTlu.MeasureRate(connectionBool);
                }
                std::string filenameFirst = filename + std::string("_background");
                myTlu.WriteOutputFile(filenameFirst, voltage, rates, thresholds);
            }
            else{
                myTlu.DoStartUp();
                for (int i = 0; i < numThresholdValues; i++){
                    myTlu.SetPMTVoltage({standardVoltage, voltage, voltage, voltage});

                    myTlu.SetTLUThreshold(thresholds[i]);
                    myTlu.SetTLUThreshold({standardThreshold}, connectionBool, "first");

                    rates[i] = myTlu.MeasureRate(connectionBool);
                }
                std::string filenameFirst = filename + std::string("_first");
                myTlu.WriteOutputFile(filenameFirst, voltage, rates, thresholds);

                // Repeat Measurement for first input, now the second input is constant
                for (int i = 0; i < numThresholdValues; i++){
                    myTlu.SetPMTVoltage({voltage, standardVoltage,voltage, voltage});
                    myTlu.SetTLUThreshold(thresholds[i]);
                    myTlu.SetTLUThreshold({standardThreshold}, connectionBool, "second");

                    rates[i] = myTlu.MeasureRate(connectionBool);
                }
                std::string filenameSecond = filename + std::string("_second");
                myTlu.WriteOutputFile(filenameSecond, voltage, rates, thresholds);
            }
//        }
    }


    std::vector<double> optimalVoltages;
    std::vector<double> optimalThresholds;
    std::vector<double> thresholdMinOpt;
    std::vector<double> thresholdMaxOpt;

    if((! tluConnected) && (!optimizeTriggerRates) && (!optimizationPlot)){
        std::vector<std::vector<double>> optimalReturn = myTlu.GetOptimalThreshold(filename);

        optimalThresholds = optimalReturn[0];
        thresholdMinOpt = optimalReturn[1];
        thresholdMaxOpt = optimalReturn[2];

        std::cout << "__________________________" << std::endl;
        std::cout << "Optimal Threshold Values:" << std::endl;
        for (int i = 0; i < numTriggerInputs; i++){
            std::cout << "PMT " << i+1 << ":   " << optimalThresholds[i]<< " V" << "   (" << "Plateau:  "<< thresholdMinOpt[i] << " - " << thresholdMaxOpt[i] << ")"<< std::endl;
        }
        std::cout << "__________________________" << std::endl;

    }

    //////////////////////////////////////////////////////////////////////////
    // Vary threshold of single channel and see how TriggerRate changes
    //////////////////////////////////////////////////////////////////////////


    optimalVoltages = {0.9, 0.95, 0.9, 0.85};
    optimalThresholds = {-0.064, -0.078, -0.078, -0.068};

    if(optimizeTriggerRates){


        std::vector<double> voltages = {0.8, 0.85, 0.9, 0.95};
        std::vector<std::string> filenames = {"200819_1", "200819_2", "200819_3", "200819_4"};
        std::cout << "CAUTION!! CHANGE FILENAMES!!"<< std::endl;


        for (int k = 0; k<4; k++){
            voltage = voltages[k];
            filename = filenames[k];
            myTlu.DoStartUp();


            //optimalVoltages = {0.9, 0.95, 0.9, 0.85};
            //optimalThresholds = {-0.064, -0.078, -0.078, -0.068};
            //thresholdMinOpt = {-0.078, -0.107, -0.103, -0.098};
            //thresholdMaxOpt = {-0.049, -0.053, -0.054, -0.039};
            thresholdMinOpt = {-0.2, -0.2, -0.2, -0.2};
            thresholdMaxOpt = {-0.0001, -0.0001, -0.0001, -0.0001};



            //three dimensional vector [numTriggerInputs] [numStepsTrigger] [numTriggerInputs+2]
            std::vector<std::vector<std::vector<double>>> ratesTrigger(numTriggerInputs, std::vector<std::vector<double>>(numStepsTrigger, std::vector<double>(numTriggerInputs + 2)));
            std::vector<std::vector<double>> thresholdsTrigger(numTriggerInputs, std::vector<double>(numStepsTrigger));

            for (int channelNo = 0; channelNo < numTriggerInputs; channelNo++){
                for (int step = 0; step < numStepsTrigger; step++){
                    thresholdsTrigger[channelNo][step] = thresholdMinOpt[channelNo] + step * (thresholdMaxOpt[channelNo] - thresholdMinOpt[channelNo]) / (numStepsTrigger-1);
                }
            }

            for (int channelNo = 0; channelNo < numTriggerInputs; channelNo++){
                std::vector<double> varThresholds = optimalThresholds;

                for (int step = 0; step < numStepsTrigger; step++){
                    varThresholds[channelNo] = thresholdsTrigger[channelNo][step];
                    //HARDCODE
                    if(channelNo==0) myTlu.SetPMTVoltage({voltage,optimalVoltages[1], optimalVoltages[2],optimalVoltages[3]});
                    else if(channelNo==1) myTlu.SetPMTVoltage({optimalVoltages[0], voltage, optimalVoltages[2],optimalVoltages[3]});
                    myTlu.SetTLUThreshold(varThresholds, connectionBool, "normal");
                    ratesTrigger[channelNo][step] = myTlu.MeasureRate(connectionBool);
                }
            }

            //for (int channel = 0; channel < numTriggerInputs; channel++){
            for (int channel = 0; channel < numTriggerInputs; channel++){
                myTlu.WriteOutputFileTrigger(channel, filename, optimalVoltages, thresholdsTrigger, ratesTrigger, optimalThresholds);
            }



            myTlu.PlotTrigger(filename);
        }
    }

    return 1;

}
