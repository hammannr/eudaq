std::string line;
int skiplines = 6;

std::vector<std::string> inFileVector;
std::string inFileList;

Double_t threshold[numThresholdValues];
Double_t rate[numTriggerInputs][numThresholdValues];
Double_t rateFirst[numTriggerInputs][numThresholdValues];
Double_t rateSecond[numTriggerInputs][numThresholdValues];
std::ifstream infile;

//// Read multiple files
infile.open(inFileList + ".txt");
if (infile.is_open())
{
    while (getline(infile,line))
    {
      std::istringstream lineS;
      lineS.str(line);

      for (int i = 0; i < numTriggerInputs + 1; i++){
        std::string val;
        lineS >> val;
        inFileVector.push_back();s
            }

        lineCounter++;
    }
    infile.close();
}
else std::cout << "Unable to open infile list";



for ()

infile.open(filename + ".txt");




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
                //std::cout << i <<"\t" <<val<< std::endl;
                if(val=="") continue;
                if (i == 0){
                    threshold[lineCounter - skiplines] = std::stod(val);
                }
                else{
                    rateFirst[i-1][lineCounter - skiplines] = std::stod(val);
                }
            }
        }
        lineCounter++;
    }
    infileFirst.close();
}
else std::cout << "Unable to open first file";

// read first file
int lineCounter = 0;
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
                //std::cout << i <<"\t" <<val<< std::endl;
                if(val=="") continue;
                if (i == 0){
                    threshold[lineCounter - skiplines] = std::stod(val);
                }
                else{
                    rateFirst[i-1][lineCounter - skiplines] = std::stod(val);
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
                //std::cout << i <<"\t" <<val<< std::endl;
                if(val=="") continue;
                if (i > 0){
                    rateSecond[i-1][lineCounter - skiplines] = std::stod(val);
                }
            }
        }
        lineCounter++;
    }
    infileSecond.close();
}
else std::cout << "Unable to open second file";



// Correct rates:
for(int i = 0; i < numTriggerInputs; i++){
    for(int j = 0; j < numThresholdValues; j++){
        if(i == 0){
            rate[i][j] = rateSecond[i][j] * rateSecond[1][0] / rateSecond[1][j];
        }
        if(i==1){
            rate[i][j] = rateFirst[i][j] * rateFirst[0][0] / rateFirst[0][j];
        }
        // for non-calibration channels use both measurements
        else{
            rate[i][j] = (rateFirst[i][j] * rateFirst[0][0] / rateFirst[0][j] + rateSecond[i][j] * rateSecond[1][0] / rateSecond[1][j]) / 2 ;
        }
    }
}



// Calculate linearized derivative of rate
Double_t derivative[numTriggerInputs][numThresholdValues];
for (int i = 0; i < numTriggerInputs; i++){
    for (int j = 0; j < numThresholdValues - 1; j++){
        derivative[i][j] =  rate[i][j+1] - rate[i][j];
    }
    derivative[i][numThresholdValues-1] = derivative[i][numThresholdValues-2]; //extrapolate last point
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





for (int i = 0; i < numTriggerInputs; i++){
    c1->cd(i+1);

    gr[i] = new TGraph (numThresholdValues,threshold,rate[i]);
    gr[i]->Draw("AL*");

    gr[i]->SetMarkerStyle(20);
    gr[i]->SetMarkerSize(0.5);
    gr[i]->SetMarkerColor(kBlue + 1);
    std::string title = std::string("PMT ") + std::to_string(i+1) + std::string("; Threshold / V; Rate / Hz");
    gr[i]->SetTitle(title.c_str());

    c2->cd(i+1);
    gr2[i] = new TGraph (numThresholdValues,threshold,derivative[i]);

    // Fit Gaussian to first derivative
    TF1 *f;
    if((i==0) || (i==1)){
        f = new TF1("f", "gaus", -0.2, -0.02);
    }
    else{
        f = new TF1("f", "gaus", -0.2, -0.05);
    }
    gr2[i]->Fit(f, "RQ"); //R: fit only in predefined range Q:Quiet
    double meanGaus = f->GetParameter(1);
    double constGaus = f->GetParameter(0);
    gr2[i]->Draw("AL*");
    gr2[i]->SetMarkerStyle(20);
    gr2[i]->SetMarkerSize(0.5);
    gr2[i]->SetMarkerColor(kBlue + 1);
    std::string titleDerivative = std::string("PMT ") + std::to_string(i+1) + std::string("; Threshold / V; a.U.");
    gr2[i]->SetTitle(titleDerivative.c_str());

    // Find Plateau
    grPlateau[i] = new TGraph();
    double coefficient = 0.5;

    //Make sure at least one Plateau point is found, otherwise: lower condition
    while ((grPlateau[i]->GetN() == 0) & (coefficient <= 1)){
        int k = 0;
        for (int j = 0; j < numThresholdValues; j++){
            if ((threshold[j] > meanGaus) & (derivative[i][j] <= coefficient * constGaus)){
                grPlateau[i]->SetPoint(k, threshold[j], rate[i][j]);
                k++;
            }

        }
        coefficient +=0.05;
    }

    if(grPlateau[i]->GetN() == 0){
        std::cout << "No Plateau could be found for PMT " << i+1 << std::endl;
        flagPlateauError = true;
    }

    else{
        // Find Midpoint of Plateau
        grMidpoint[i] = new TGraph();
        int lenPlateau = grPlateau[i]->GetN();
        int indexMidpoint = (lenPlateau + 1*lenPlateau%2)/2 - 1;
        std::cout << "index mid   " << indexMidpoint << std::endl;
        std::cout << "len plateau   " << lenPlateau << std::endl;
        grMidpoint[i]->SetPoint(0, grPlateau[i]->GetX()[indexMidpoint], grPlateau[i]->GetY()[indexMidpoint]);

        std::cout << ":::::::::::::::::::::::::::::::::::::" << std::endl;
        for (int j = 0; j < lenPlateau; j++){
            std::cout << grPlateau[i]->GetX()[j] << std::endl;
        }
        std::cout << ":::::::::::::::::::::::::::::::::::::" << std::endl;


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
        grMidpoint[i]->Draw("*");
        grMidpoint[i]->SetMarkerStyle(20);
        grMidpoint[i]->SetMarkerSize(1);
        grMidpoint[i]->SetMarkerColor(kGreen + 1);
        grMidpoint[i]->SetTitle("Optimal Threshold");


        //Plot legend
        auto legend = new TLegend(0.1,0.8,0.58,0.9);
        std::string midPointString = std::to_string(grMidpoint[i]->GetX()[0]);
        midPointString.erase(midPointString.begin()+6, midPointString.end()); //limit number of values after comma
        std::string labelMidpoint = std::string("Optimal Threshold:  ") + midPointString + std::string(" V");
        legend->AddEntry(grMidpoint[i],labelMidpoint.c_str(),"p");
        legend->Draw();

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
