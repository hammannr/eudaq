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


//#include "gnuplot-iostream.h"

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
    void WriteOutputFileTrigger(int channel, std::string filename, std::vector<double> voltage, std::vector<std::vector<double>> thresholdsTrigger, std::vector<std::vector<std::vector<double>>> ratesTrigger, std::vector<double> optimalThresholds);
    void WriteParameters(int in_numThresholdValues, int in_numTriggerInputs, int in_time, int in_numStepsTrigger);
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

AidaTluControl::AidaTluControl(){
    m_verbose = 0x0;
    m_duration = 0;
    m_starttime = 0;
    m_lasttime = 0;
    flagPlateauError = false;
}

void AidaTluControl::WriteParameters(int in_numThresholdValues, int in_numTriggerInputs, int in_time, int in_numStepsTrigger){
    numThresholdValues = in_numThresholdValues;
    numTriggerInputs = in_numTriggerInputs;
    time = in_time;
    numStepsTrigger = in_numStepsTrigger;
}

// Initialize TLU
void AidaTluControl::DoStartUp(){
    /* Establish a connection with the TLU using IPBus.
       Define the main hardware parameters.
    */

    std::string uhal_conn = "file:///opt/eudaq2/user/eudet/misc/hw_conf/aida_tlu/aida_tlu_connection.xml";

    std::string uhal_node = "aida_tlu.controlhub";
    //std::string uhal_conn;
    //    std::string uhal_node;
    //    uhal_conn = ini->Get("ConnectionFile", uhal_conn);
    //    uhal_node = ini->Get("DeviceName",uhal_node);
    m_tlu = std::unique_ptr<tlu::AidaTluController>(new tlu::AidaTluController(uhal_conn, uhal_node));

    m_verbose = 0x1;
    m_tlu->DefineConst(0, 6); // 0 DUTs, 6 Triggers



    // Populate address list for I2C elements
    m_tlu->SetI2C_core_addr(0x21);
    m_tlu->SetI2C_clockChip_addr(0x68);
    m_tlu->SetI2C_DAC1_addr(0x13);
    m_tlu->SetI2C_DAC2_addr(0x1f);
    m_tlu->SetI2C_EEPROM_addr(0x50);
    m_tlu->SetI2C_expander1_addr(0x74);
    m_tlu->SetI2C_expander2_addr(0);
    m_tlu->SetI2C_pwrmdl_addr(0x1C, 0x76, 0x77, 0x51);
    m_tlu->SetI2C_disp_addr(0x3A);

    // Initialize TLU hardware
    m_tlu->InitializeI2C(m_verbose);
    m_tlu->InitializeIOexp(m_verbose);
    // 1.3V = reference voltage
    m_tlu->InitializeDAC(false, 1.3, m_verbose);


    // Initialize the Si5345 clock chip using pre-generated file
    //std::string defaultCfgFile= "file:///opt/eudaq2/user/eudet/misc/hw_conf/aida_tlu/fmctlu_clock_config.txt";
    std::string defaultCfgFile= "/opt/eudaq2/user/eudet/misc/hw_conf/aida_tlu/aida_tlu_clk_config.txt";
    int clkres;
    clkres= m_tlu->InitializeClkChip( defaultCfgFile, m_verbose  );
    if (clkres == -1){
        std::cout << "TLU: clock configuration failed." << std::endl;
    }

    // Set trigger stretch

    std::vector<unsigned int> stretcVec = {(unsigned int) 1,
                                           (unsigned int) 1,
                                           (unsigned int) 1,
                                           (unsigned int) 1,
                                           (unsigned int) 1,
                                           (unsigned int) 1};
    m_tlu->SetPulseStretchPack(stretcVec, m_verbose);

    // Set trigger mask (high,low)
    // MIGHT BE ADJUSTED!!
    m_tlu->SetTriggerMask((uint32_t)0x00000000,  (uint32_t)0x0008000);

    // Reset IPBus registers
    m_tlu->ResetSerdes();
    m_tlu->ResetCounters();

    m_tlu->SetTriggerVeto(1, m_verbose); // no triggers
    m_tlu->ResetFIFO();
    m_tlu->ResetEventsBuffer();
    //m_tlu->ResetBoard();

    m_tlu->ResetTimestamp();

    m_tlu->SetDUTMask(0x1, m_verbose);
    m_tlu->SetDUTIgnoreBusy(0xF,m_verbose);
    m_tlu->SetDUTMaskMode(0xFF,m_verbose);
    m_tlu->enableClkLEMO(true, m_verbose);

    m_tlu->SetEnableRecordData(0x0);
    m_tlu->GetEventFifoCSR();
    m_tlu->GetEventFifoFillLevel();

}


// Set PMT Voltage
void AidaTluControl::SetPMTVoltage(double val){
    m_tlu->pwrled_setVoltages(float(val), float(val), float(val), float(val), m_verbose);
}

// even if it outputs something else, voltage[0] corresponds to PMT 1, voltage[1] to PMT2 ,... (checked with multimeter)
void AidaTluControl::SetPMTVoltage(std::vector<double> voltage){
    m_tlu->pwrled_setVoltages(float(voltage[0]), float(voltage[1]), float(voltage[2]), float(voltage[3]), m_verbose);
}

// Set TLU threshold
void AidaTluControl::SetTLUThreshold(double val){
    m_tlu->SetThresholdValue(7, float(val) , m_verbose); //all channels
}

// Set TLU threshold (individual values and only connected channels)
void AidaTluControl::SetTLUThreshold(std::vector<double> val, std::vector<bool> connection, std::string mode){
    // mode first: Writes voltage only to first connected
    int i = 0;
    bool done = false;

    if (connection[0] && (!done)) {
        if (mode.compare("second") != 0) {
            m_tlu->SetThresholdValue(0, float(val[i]) , m_verbose);
            i++;
        }
        if (mode.compare("first") == 0){
            done = true;
        }
        else if (mode.compare("second") == 0){
            mode = "first";
        }
    }
    if (connection[1] && (!done)) {
        if (mode.compare("second") != 0) {
            m_tlu->SetThresholdValue(1, float(val[i]) , m_verbose);
            i++;
        }
        if (mode.compare("first") == 0){
            done = true;
        }
        else if (mode.compare("second") == 0){
            mode = "first";
        }
    }
    if (connection[2] && (!done)) {
        if (mode.compare("second") != 0) {
            m_tlu->SetThresholdValue(2, float(val[i]) , m_verbose);
            i++;
        }
        if (mode.compare("first") == 0){
            done = true;
        }
        else if (mode.compare("second") == 0){
            mode = "first";
        }
    }
    if (connection[3] && (!done)) {
        if (mode.compare("second") != 0) {
            m_tlu->SetThresholdValue(3, float(val[i]) , m_verbose);
            i++;
        }
        if (mode.compare("first") == 0){
            done = true;
        }
        else if (mode.compare("second") == 0){
            mode = "first";
        }
    }
    if (connection[4] && (!done)) {
        if (mode.compare("second") != 0) {
            m_tlu->SetThresholdValue(4, float(val[i]) , m_verbose);
            i++;
        }
        if (mode.compare("first") == 0){
            done = true;
        }
        else if (mode.compare("second") == 0){
            mode = "first";
        }
    }
    if (connection[5] && (!done)) {
        if (mode.compare("second") != 0) {
            m_tlu->SetThresholdValue(5, float(val[i]) , m_verbose);
            i++;
        }
        if (mode.compare("first") == 0){
            done = true;
        }
        else if (mode.compare("second") == 0){
            mode = "first";
        }
    }

}

// Measure rate
std::vector<double> AidaTluControl::MeasureRate(std::vector<bool> connectionBool){

    std::vector<uint32_t> sl={0,0,0,0,0,0};
    // Output: First for: TLU, last: Trigger Rate (pre & post veto)
    // Convert String input into bool input & count trigger inputs

    std::vector<double> output(numTriggerInputs + 2, 0);

    std::this_thread::sleep_for (std::chrono::seconds(1));
    m_tlu->ResetCounters();
    m_tlu->ResetSerdes();

    m_tlu->SetRunActive(1, 1); // reset internal counters
    m_starttime = m_tlu->GetCurrentTimestamp()*25;
    m_tlu->SetTriggerVeto(0, m_verbose); //enable trigger
    m_tlu->ReceiveEvents(m_verbose);

    std::this_thread::sleep_for (std::chrono::milliseconds(time*1000));

    m_tlu->SetTriggerVeto(1, m_verbose); //disable trigger
    // Set TLU internal logic to stop.
    m_tlu->SetRunActive(0, 1);
    m_lasttime = m_tlu->GetCurrentTimestamp()*25;

    m_tlu->GetScaler(sl[0], sl[1], sl[2], sl[3], sl[4], sl[5]);
    output[numTriggerInputs] = m_tlu->GetPreVetoTriggers() / double(time);
    output[numTriggerInputs + 1] = m_tlu->GetPostVetoTriggers() / double(time);


    std::cout << "Scalers:  "<<std::dec << sl[0] << "  " << sl[1]<< "  " << sl[2]<< "  " << sl[3]<< "  " << sl[4]<< "  " << sl[5] << std::endl;
    std::cout << "Triggers:  " <<  m_tlu->GetPreVetoTriggers() <<std::endl;
    m_tlu->ResetCounters();
    m_tlu->ResetSerdes();



    m_duration = double(m_lasttime - m_starttime) / 1000000000; // in seconds
    std::cout << "Run duration [s]" << m_duration << std::endl;

    int i = 0;
    for (int k = 0; k < 6; k++){
        if (connectionBool[k]){
            output[i] = sl[k] / double(time);
            i++;
        }
    }
    return {output};
}

std::vector<std::vector<std::vector<double>>> AidaTluControl::readFiles(std::string filename){
//    Double_t threshold[numThresholdValues];
//    Double_t rateFirst[numTriggerInputs][numThresholdValues];
//    Double_t rateSecond[numTriggerInputs][numThresholdValues];


    std::vector<std::vector<std::vector<double>>> returnValue (3,std::vector<std::vector<double>>(numTriggerInputs,std::vector<double>(numThresholdValues))) ;//[3][numTriggerInputs][numThresholdValues];
    std::string line;
    int skiplines = 6;
    int lineCounter = 0;

    // read first file
    std::ifstream infileFirst;
    infileFirst.open(filename + "_first" + ".txt");

    if (infileFirst.is_open())
    {
        while (getline(infileFirst,line))
        {
            if (lineCounter >= skiplines){
                std::istringstream lineS;
                lineS.str(line);

                for (int i = 0; i < numTriggerInputs + 1; i++){
                    std::string val;
                    lineS >> val;
                    if(val=="") continue;
                    if (i == 0){
                        returnValue[0][0][lineCounter - skiplines] = std::stod(val);
                    }
                    else{
                        returnValue[1][i-1][lineCounter - skiplines] = std::stod(val);
                    }
                }
            }
            lineCounter++;
        }
        infileFirst.close();
    }
    else std::cout << "Unable to open first file";

    // read second file
    lineCounter = 0;
    std::ifstream infileSecond;
    infileSecond.open(filename + "_second" + ".txt");

    if (infileSecond.is_open())
    {
        while (getline(infileSecond,line))
        {
            if (lineCounter >= skiplines){
                std::istringstream lineS;
                lineS.str(line);

                for (int i = 0; i < numTriggerInputs + 1; i++){
                    std::string val;
                    lineS >> val;
                    if(val=="") continue;
                    if (i > 0){
                        returnValue[2][i-1][lineCounter - skiplines] = std::stod(val);
                    }
                }
            }
            lineCounter++;
        }
        infileSecond.close();
    }
    else std::cout << "Unable to open second file";
    return returnValue;
}

/*
// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
   return par[3] + par[0]/sqrt(2*TMath::Pi()*pow(par[2],2)) * exp(-pow((x - par[1]),2) / (2*pow(par[2],2)));
}
*/
std::vector<std::vector<double>> AidaTluControl::GetOptimalThreshold(std::string filename){
    //std::vector<std::string> appendices = {"_first", "_second"};
    // Open File with data and write them into arrays

    //    Double_t threshold[numThresholdValues];
    Double_t rate[numTriggerInputs][numThresholdValues];
    Double_t threshold[numThresholdValues];
    Double_t rateFirst[numTriggerInputs][numThresholdValues];
    Double_t rateSecond[numTriggerInputs][numThresholdValues];
    std::vector<std::vector<std::vector<double>>> files (3,std::vector<std::vector<double>>(numTriggerInputs,std::vector<double>(numThresholdValues))) ;
    //    Double_t rateFirst[numTriggerInputs][numThresholdValues];
    //    Double_t rateSecond[numTriggerInputs][numThresholdValues];

    files = readFiles(filename);
    for (int i = 0; i < numThresholdValues; i++){
        threshold[i] = files[0][0][i];
        for (int j = 0; j< numTriggerInputs; j++){
            rateFirst[j][i] = files[1][j][i];
            rateSecond[j][i] = files[2][j][i];
        }
    }

    // Correct rates:
    for(int i = 0; i < numTriggerInputs; i++){
        for(int j = 0; j < numThresholdValues; j++){
            if(i == 0){
                rate[i][j] = rateSecond[i][j] / rateSecond[1][j];
            }
            if(i==1){
                rate[i][j] = rateFirst[i][j] / rateFirst[0][j];
            }
            // for non-calibration channels use both measurements
            else{
                rate[i][j] = (rateFirst[i][j] / rateFirst[0][j] + rateSecond[i][j]/ rateSecond[1][j]) / 2 ;
            }
        }
    }



    // Calculate linearized derivative of rate
    Double_t derivative[numTriggerInputs][numThresholdValues];
    for (int i = 0; i < numTriggerInputs; i++){
        for (int j = 0; j < numThresholdValues - 1; j++){
            derivative[i][j] = (rate[i][j+1] - rate[i][j]) / (threshold[1] - threshold[0]);
        }
        derivative[i][numThresholdValues - 1] = derivative[i][numThresholdValues-2]; //extrapolate last point
    }



    //Plot Data and determine optimum

    std::vector<double> optimalThreshold(numTriggerInputs);
    std::vector<double> thresholdMinOpt(numTriggerInputs);
    std::vector<double> thresholdMaxOpt(numTriggerInputs);

    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 0,10,940,1100);
    TCanvas *c2 = new TCanvas("c2", "Graph Draw Options", 980,10,940,1100);

    //    c1->SetRightMargin(0.09);
    //    c1->SetLeftMargin(0.1);
    //    c1->SetBottomMargin(0.15);
    //c1->Divide(numTriggerInputs%3+1,numTriggerInputs/3+1);
    if (numTriggerInputs == 1) c1->Divide(1,1);
    else if (numTriggerInputs == 2) c1->Divide(2,1);
    else if (numTriggerInputs == 3) c1->Divide(3,1);
    else if (numTriggerInputs == 4) c1->Divide(2,2);
    else if (numTriggerInputs == 5) c1->Divide(3,2);
    else if (numTriggerInputs == 6) c1->Divide(3,2);
    else if (numTriggerInputs > 6) c1->Divide(3,3);
    else c1->Divide(5,3);

    if (numTriggerInputs == 1) c2->Divide(1,1);
    else if (numTriggerInputs == 2) c2->Divide(2,1);
    else if (numTriggerInputs == 3) c2->Divide(3,1);
    else if (numTriggerInputs == 4) c2->Divide(2,2);
    else if (numTriggerInputs == 5) c2->Divide(3,2);
    else if (numTriggerInputs == 6) c2->Divide(3,2);
    else if (numTriggerInputs > 6) c2->Divide(3,3);
    else c2->Divide(5,3);

    TGraph *gr[numTriggerInputs];
    TGraph *gr2[numTriggerInputs];
    TGraph *grPlateau[numTriggerInputs];
    TGraph *grMidpoint[numTriggerInputs];
    TGraph *grPlateauDeriv[numTriggerInputs];





    for (int i = 0; i < numTriggerInputs; i++){\

        c1->cd(i+1);

        gr[i] = new TGraph (numThresholdValues,threshold,rate[i]);
        gr[i]->Draw("AL*");

        gr[i]->SetMarkerStyle(20);
        gr[i]->SetMarkerSize(0.5);
        gr[i]->SetMarkerColor(kBlue - 2);
        gr[i]->SetLineColor(kBlue - 2);
        std::string title = std::string("PMT ") + std::to_string(i+1) + std::string("; Threshold / V; Rate / Hz");
        gr[i]->SetTitle(title.c_str());

        c2->cd(i+1);
        gr2[i] = new TGraph (numThresholdValues,threshold,derivative[i]);

        // Fit Gaussian to first derivative
        TF1 *f;        

        //f = new TF1("f",  "[3]+[0]/sqrt(2*pi*[2]^2)*exp(-1*(x-[1])^2/(2*[2]^2))", -0.2, -0.05);
        f = new TF1("f",  "gaus", -0.2, -0.05);
        //f = new TF1("f",  "[0] * exp(-1*(x - [1])^2)", -0.2, -0.05);
        f->SetParNames("Constant","Mean","Sigma","Offset");
        //f = new TF1("f",  "[0]*x^2*sqrt([1])*pi+exp(-x)", -0.2, -0.05);

        gr2[i]->Fit(f, "R"); //R: fit only in predefined range Q:Quiet
        double meanGaus = f->GetParameter(1);
        double constGaus = f->GetParameter(0);
        gr2[i]->Draw("AL*");
        gr2[i]->SetMarkerStyle(20);
        gr2[i]->SetMarkerSize(0.5);
        gr2[i]->SetMarkerColor(kBlue - 2);
        gr2[i]->SetLineColor(kBlue - 2);
        std::string titleDerivative = std::string("PMT ") + std::to_string(i+1) + std::string("; Threshold / V; a.U.");
        gr2[i]->SetTitle(titleDerivative.c_str());

        // Find Plateau
        grPlateau[i] = new TGraph();
        grPlateauDeriv[i] = new TGraph();
        double coefficient = 0.5;
        //constGaus = 10;

        //Make sure at least one Plateau point is found, otherwise: lower condition
        while ((grPlateau[i]->GetN() == 0) & (coefficient <= 1)){
            int k = 0;
            // shift derivative one to the right
            for (int j = 1; j < numThresholdValues; j++){
                if ((threshold[j] > meanGaus) & (derivative[i][j - 1] <= coefficient * constGaus)){
                    grPlateau[i]->SetPoint(k, threshold[j], rate[i][j]);
                    grPlateauDeriv[i]->SetPoint(k, threshold[j - 1], derivative[i][j - 1]);
                    k++;
                }
                // also add point if both next neighbours fulfill plateau conditions
                else if ((j > 0) & (j < numThresholdValues-2) & (threshold[j-2] > meanGaus) & (derivative[i][j-2] <= coefficient * constGaus) & (threshold[j] > meanGaus) & (derivative[i][j] <= coefficient * constGaus)){
                    grPlateau[i]->SetPoint(k, threshold[j], rate[i][j]);
                    grPlateauDeriv[i]->SetPoint(k, threshold[j - 1], derivative[i][j - 1]);
                    k++;
                }

            }
            coefficient +=0.05;
        }

        if(grPlateau[i]->GetN() == 0){
            std::cout << "No Plateau could be found for PMT " << i+1 << " (no points found)"<< std::endl;
            flagPlateauError = true;
        }
        else if(grPlateau[i]->GetN() == numThresholdValues){
            std::cout << "No Plateau could be found for PMT "<< i+1 << " (all points identified as plateau)" << std::endl;
            grPlateau[i] = new TGraph();
            flagPlateauError = true;
        }

        else{
            // Find Midpoint of Plateau
            grMidpoint[i] = new TGraph();
            int lenPlateau = grPlateau[i]->GetN();

            int indexMidpoint = (lenPlateau + (lenPlateau % 2))/2 - (lenPlateau % 2);
            grMidpoint[i]->SetPoint(0, grPlateau[i]->GetX()[indexMidpoint], grPlateau[i]->GetY()[indexMidpoint]);

            if (lenPlateau < 3){
                double diff = gr[i]->GetX()[1] - gr[i]->GetX()[0];
                thresholdMinOpt[i] = grPlateau[i]->GetX()[indexMidpoint] - diff;
                thresholdMaxOpt[i] = grPlateau[i]->GetX()[indexMidpoint] + diff;
            }
            else{
                thresholdMinOpt[i] = grPlateau[i]->GetX()[0];
                thresholdMaxOpt[i] = grPlateau[i]->GetX()[lenPlateau - 1];
            }



            // Plot Plateau and Midpoint
            c1->cd(i+1);
            grPlateau[i]->Draw("*");
            grPlateau[i]->SetMarkerStyle(20);
            grPlateau[i]->SetMarkerSize(1);
            grPlateau[i]->SetMarkerColor(kRed + 1);
            //grMidpoint[i]->SetTitle("Optimal Threshold");
            grMidpoint[i]->Draw("*");
            grMidpoint[i]->SetMarkerStyle(20);
            grMidpoint[i]->SetMarkerSize(1);
            grMidpoint[i]->SetMarkerColor(kGreen + 1);

            //Plot legend
            auto legend = new TLegend(0.1,0.8,0.58,0.9);
            std::string midPointString = std::to_string(grMidpoint[i]->GetX()[0]);
            midPointString.erase(midPointString.begin()+6, midPointString.end()); //limit number of values after comma
            std::string labelMidpoint = std::string("Optimal Threshold:  ") + midPointString + std::string(" V");
            legend->AddEntry(grMidpoint[i],labelMidpoint.c_str(),"p");
            legend->Draw();

            // Plot Plateau and Midpoint
            c2->cd(i+1);
            grPlateauDeriv[i]->Draw("*");
            grPlateauDeriv[i]->SetMarkerStyle(20);
            grPlateauDeriv[i]->SetMarkerSize(1);
            grPlateauDeriv[i]->SetMarkerColor(kRed + 1);



            optimalThreshold[i] = grMidpoint[i]->GetX()[0];
        }

    }


    TGaxis::SetMaxDigits(3);
    c1->SetLogy();
    c1->Update();
    c1->Modified();
    c2->Update();
    c2->Modified();

    std::string exportFile = filename + (std::string)"_rates.pdf";
    c1->SaveAs(exportFile.c_str());
    exportFile = filename + (std::string)"_rates.root";
    c1->SaveAs(exportFile.c_str());
    std::string exportFile2 = filename + (std::string)"_derivative.pdf";
    c2->SaveAs(exportFile2.c_str());
    exportFile2 = filename + (std::string)"_derivative.root";
    c2->SaveAs(exportFile2.c_str());
    std::cout << "Enter 'q' to continue" << std::endl;
    int k;
    std::cin >> k;

    //myApp->Run();
    //myApp->Terminate();


    return {optimalThreshold, thresholdMinOpt, thresholdMaxOpt};
}



void AidaTluControl::PlotTrigger(std::string filename){

    // Open File with data and write them into arrays
    std::ifstream infile;
    Double_t threshold[numTriggerInputs][numStepsTrigger];
    Double_t rate[numTriggerInputs][numStepsTrigger];
    Double_t correctionFactor;

    for (int channel = 0; channel < numTriggerInputs; channel++){

        infile.open(filename + "_OptimizeTrigger_" + std::to_string(channel+1) + ".txt");
        std::string line;
        int skiplines = 7;
        int lineCounter = 0;

        // read file
        if (infile.is_open())
        {
            while ( getline (infile,line) )
            {
                if (lineCounter >= skiplines){

                    std::istringstream lineS;
                    lineS.str(line);

                    for (int i = 0; i < numTriggerInputs + 2; i++){
                        std::string val;
                        lineS >> val;
                        //if(val=="") continue;
                        // get rate of first channel as normalization for each point (except for first channel where you get 2nd channel as reference)
                        if (channel == 0){
                            if (i == 2){
                                correctionFactor = std::stod(val);
                            }
                        }
                        else{
                            if (i == 1){
                                correctionFactor = std::stod(val);
                            }
                        }

                        if (i == 0){
                            threshold[channel][lineCounter - skiplines] = std::stod(val);
                        }

                        // Take PostVetoTrigger:
                        else if (i == numTriggerInputs + 1){
                            rate[channel][lineCounter - skiplines] = std::stod(val) / correctionFactor;
                        }
                    }
                }
                lineCounter++;
            }
            infile.close();
        }
        else std::cout << "Unable to open file";
    }



    //Plot Data


    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 0,10,940,1100);


    if (numTriggerInputs == 1) c1->Divide(1,1);
    else if (numTriggerInputs == 2) c1->Divide(2,1);
    else if (numTriggerInputs == 3) c1->Divide(3,1);
    else if (numTriggerInputs == 4) c1->Divide(2,2);
    else if (numTriggerInputs == 5) c1->Divide(3,2);
    else if (numTriggerInputs == 6) c1->Divide(3,2);
    else if (numTriggerInputs > 6) c1->Divide(3,3);
    else c1->Divide(5,3);

    TGraph *gr[numTriggerInputs];

    for (int i = 0; i < numTriggerInputs; i++){
        c1->cd(i+1);

        gr[i] = new TGraph (numStepsTrigger,threshold[i],rate[i]);
        gr[i]->Draw("AL*");

        gr[i]->SetMarkerStyle(20);
        gr[i]->SetMarkerSize(0.5);
        gr[i]->SetMarkerColor(kBlue - 2);
        gr[i]->SetLineColor(kBlue - 2);
        std::string title = std::string("Vary PMT ") + std::to_string(i+1) + std::string("; Threshold ") + std::to_string(i+1) + std::string("/ V; Trigger Rate / a.U.");
        gr[i]->SetTitle(title.c_str());
    }


    TGaxis::SetMaxDigits(3);
    c1->Update();
    c1->Modified();

    std::string exportFile = filename + (std::string)"_OptimizeTrigger.pdf";
    c1->SaveAs(exportFile.c_str());

    exportFile = filename + (std::string)"_OptimizeTrigger.root";
    c1->SaveAs(exportFile.c_str());

    std::cout << "Enter 'q' to continue" << std::endl;
    int k;
    std::cin >> k;

}





void AidaTluControl::WriteOutputFile(std::string filename, double voltage, std::vector<std::vector<double>> rates, std::vector<double> thresholds){
    std::ofstream outFile;
    outFile.open (filename + ".txt");

    auto now = std::chrono::system_clock::now();
    std::time_t timeNow = std::chrono::system_clock::to_time_t(now);
    outFile << "Date:\t" << std::ctime(&timeNow) << "\n";
    outFile << "PMT Voltage [V]:\t   " << voltage << "\n";
    outFile << "Acquisition Time [s]\t:   " << time << "\n" << "\n";

    outFile << "Thr [V]\t\t";
    for (int i = 0; i < numTriggerInputs; i++){
        outFile << "PMT " << i+1 << " [Hz]\t";
    }
    outFile << "PreVeto [Hz]\t";
    outFile << "PostVeto [Hz]\t";
    outFile << "\n";


    for (int i = 0; i < numThresholdValues; i++){
        outFile << thresholds[i] << "\t\t";
        for (auto r:rates[i]) outFile << r << "\t\t";
        outFile << "\n";
    }

    outFile.close();
}


void AidaTluControl::WriteOutputFileTrigger(int channel, std::string filename, std::vector<double> voltage, std::vector<std::vector<double>> thresholdsTrigger, std::vector<std::vector<std::vector<double>>> ratesTrigger, std::vector<double> optimalThresholds){
    std::ofstream outFile;
    outFile.open (filename + "_OptimizeTrigger_" + std::to_string(channel+1) + ".txt");

    auto now = std::chrono::system_clock::now();
    std::time_t timeNow = std::chrono::system_clock::to_time_t(now);
    outFile << "Date:\t" << std::ctime(&timeNow) << "\n";
    outFile << "Acquisition Time [s]:\t   " << time << "\n";

    outFile << "Optimal PMT Voltages [V]:\t   ";
    for (int channel = 0; channel < numTriggerInputs; channel++){
        outFile <<"Channel " << channel + 1<< ")  " << voltage[channel] << "\t";
    }
    outFile << "\n";


    outFile << "Optimal Thresholds [V]:\t   ";
    for (int channel = 0; channel < numTriggerInputs; channel++){
        outFile <<"Channel " << channel + 1<< ")  " << optimalThresholds[channel] << "\t";
    }
    outFile << "\n\n";

    outFile << "Thr " << channel+1 << " [V]\t";
    for (int channel_ = 0; channel_ < numTriggerInputs; channel_++){
        outFile << "PMT " << channel_+1 << " [Hz]\t";
    }
    outFile << "PreVeto [Hz]\t";
    outFile << "PostVeto [Hz]\t";
    outFile << "\n";

    for (int step = 0; step < thresholdsTrigger[0].size(); step++){
        outFile << thresholdsTrigger[channel][step] << "\t\t";
        for (auto r:ratesTrigger[channel][step]) outFile << r << "\t\t";
        outFile << "\n";
    }


    outFile.close();
}


int main(int /*argc*/, char **argv) {
    eudaq::OptionParser op("EUDAQ Command Line FileReader modified for TLU", "2.1", "EUDAQ FileReader (TLU)");
    eudaq::Option<double> thrMin(op, "tl", "thresholdlow", -0.21, "double", "threshold value low [V]");
    eudaq::Option<double> thrMax(op, "th", "thresholdhigh", -0.01, "double", "threshold value high [V]");
    eudaq::Option<int> thrNum(op, "tn", "thresholdsteps", 9, "int", "number of threshold steps");
    eudaq::Option<double> volt(op, "v", "pmtvoltage", 0.8, "double", "PMT voltage [V]");
    eudaq::Option<int> acqtime(op, "t", "acquisitiontime", 10, "int", "acquisition time");
    eudaq::Option<std::string> name(op, "f", "filename", "output", "string", "filename");
    eudaq::Option<std::string> con(op, "c", "connectionmap", "111100", "string", "connection map");
    eudaq::Option<bool> tluOn(op, "tlu", "tluconnected", true, "bool", "tlu connected");
    eudaq::Option<bool> optTrig(op, "opt", "optimizetrig", false, "bool", "optimize for trigger rates");
    eudaq::Option<bool> optPlot(op, "optplot", "optplot", false, "bool", "only plot optimization stuff");
    eudaq::Option<int> trigNum(op, "otn", "optimizetrignum", 11, "int", "number of steps for trigger optimization");
    eudaq::Option<int> deb(op, "deb", "debug", 0, "int", "debug number");
    //eudaq::Option<bool> setOpt(op, "setopt", "setoptparams", false, "bool", "set optimal parameters");
    //eudaq::Option<std::vector<double>> voltOpt(op, "ov", "optpmtvoltages", (0.8,0.8,0.8,0.8), "vector double", "optimal PMT voltages [V]");
    //eudaq::Option<std::vector<double>> thrOpt(op, "ot", "optthresholds", (-0.08,-0.08,-0.08,-0.08), "vector double", "optimal threshold values [V]");
    //eudaq::Option<std::vector<double>> thrMinOpt(op, "otl", "optthresholdslow", (-0.1,-0.1,-0.1,-0.1), "vector double", "optimal threshold values low [V]");
    //eudaq::Option<std::vector<double>> thrMaxOpt(op, "oth", "optthresholdshigh", (-0.06,-0.06,-0.06,-0.06), "vector double", "optimal threshold values high [V]");


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
        std::cout << "Enter 'q' to continue" << std::endl;
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
        std::cout << "Enter 'q' to continue" << std::endl;
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
    double standardVoltage = 0.9;

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
            std::cout << "TLU " << i+1 << ":   " << optimalThresholds[i]<< " V" << "   (" << "Plateau:  "<< thresholdMinOpt[i] << " - " << thresholdMaxOpt[i] << ")"<< std::endl;
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
/*
    if (optimizationPlot){
        optimalVoltages = {0.9, 0.95, 0.9, 0.85};
        optimalThresholds = {-0.064, -0.078, -0.078, -0.068};
        thresholdMinOpt = {-0.078, -0.107, -0.103, -0.098};
        thresholdMaxOpt = {-0.049, -0.053, -0.054, -0.039};

        myTlu.PlotTrigger(filename);
    }*/


    /*
    if (debugNumber == 1){
        double_t voltages = {0.85, 0.9, 0.95};
        TFile *f;
        TCanvas *c;

        f = new TFile("060819_7_rates.root");
        c = (TCanvas *)f ->Get("c1");
        c->cd(0);
        c->GetListOfPrimitives()->Print();
        //TGraph *g = (TGraph *) c->Get("Graph");

    }*/

    return 1;

}
