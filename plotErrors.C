#include "../inc/std_headers.h"
#include "../inc/root.h"

using namespace std;

void plotErrors(string inFile) { 

  ifstream dataFile;

  int numEvents = 0;
  int fileCount = 0;
  int duplicateEvent[5] = {0, 0, 0, 0, 0};
  int skippedEvent[5] = {0, 0, 0, 0, 0};
  int runNum = 0;


  std::vector<float> duplicateFraction[4];
  std::vector<float> skippedFraction[4];
  std::vector<float> runNumber;
  std::vector<float> tmpvect;

  dataFile.open(inFile.c_str());
  if (dataFile.bad()) { 
    std::cout << "Error opening data file: " << inFile << std::endl;
    return;
  } 

  std::string line;
  std::string delimiter =":";
  std::string token;

  while (!dataFile.eof() && !dataFile.bad()) { 

    // numEvents = 0;

    fileCount++;

    for (int i = 0; i < 15; i++) { 
      getline(dataFile, line);
      if (dataFile.eof()) break;
      if (i == 2) { 
	cout << line << endl;
	token = line.substr(38, 8);
	runNum = atoi(token.c_str());
      } else if (i == 5) { 
	token = line.substr(28, line.length());
	numEvents = atoi(token.c_str());
      } else if (i == 9 ) { 
	int start = 0;
	int end = line.length();
	int pos = 0;
	
	for (int j = 0; j < 8; j++) {
	  pos = line.find(delimiter, start);
	  token = line.substr(start, pos-start);
	  start = pos+1;	 
	  if (j > 3) 
	    duplicateEvent[j-4] = atoi(token.c_str());
	}
      } else if (i == 10) { 
	int start = 0;
	int end = line.length();
	int pos = 0;
	
	for (int j = 0; j < 8; j++) {
	  pos = line.find(delimiter, start);
	  token = line.substr(start, pos-start);
	  start = pos+1;	 
	  if (j > 3)
	    skippedEvent[j-4] = atoi(token.c_str());
	}
      }
      
    }
    
    float frac; 
    if (numEvents > 0) {
      runNumber.push_back((float)runNum);
      std::cout << "Run Number: " << runNum << std::endl;
      for (int i = 0; i < 4; i++) {
	frac = (float) duplicateEvent[i]/(float) numEvents * 100.;
	duplicateFraction[i].push_back(frac);
	frac = (float) skippedEvent[i]/(float) numEvents * 100.;
	skippedFraction[i].push_back(frac);
      } 
    }

  }

  TGraph *duplicateGraph[4];
  TGraph *skippedGraph[4];

  TString title;

  for (int i = 0; i < 4; i++) {  
    tmpvect = std::vector<float>();
    tmpvect = duplicateFraction[i];
    duplicateGraph[i] = new TGraph((int)runNumber.size(), &runNumber[0], &tmpvect[0]);

    if ( i == 0)
      title = "XMIT % Duplicate Events";
    else { 
      title = "FEM-";
      title += i-1;
      title += " % Duplicate Events";
    }
    duplicateGraph[i]->SetTitle(title);

    duplicateGraph[i]->GetXaxis()->SetTitle("Run Number");
    duplicateGraph[i]->GetXaxis()->SetTitleOffset(0.3);
    duplicateGraph[i]->GetXaxis()->SetLabelSize(0.025);
    duplicateGraph[i]->GetXaxis()->SetLabelOffset(0.05);
    duplicateGraph[i]->GetYaxis()->SetTitle("% Duplicated");
    duplicateGraph[i]->GetYaxis()->SetTitleOffset(0.3);
    duplicateGraph[i]->GetYaxis()->SetLabelSize(0.025);
    duplicateGraph[i]->GetYaxis()->SetLabelOffset(0.05);


    tmpvect = std::vector<float>();
    tmpvect = skippedFraction[i];
    skippedGraph[i] = new TGraph((int)runNumber.size(), &runNumber[0], &tmpvect[0]);
    if ( i == 0)
      title = "XMIT % Skipped Events";
    else { 
      title = "FEM-";
      title += i-1;
      title += " % Skipped Events";
    }
    skippedGraph[i]->SetTitle(title);

    skippedGraph[i]->GetXaxis()->SetTitle("Run Number");
    skippedGraph[i]->GetXaxis()->SetTitleOffset(0.3);
    skippedGraph[i]->GetXaxis()->SetLabelSize(0.025);
    skippedGraph[i]->GetXaxis()->SetLabelOffset(0.05);
    skippedGraph[i]->GetYaxis()->SetTitle("% Skippedd");
    skippedGraph[i]->GetYaxis()->SetTitleOffset(0.3);
    skippedGraph[i]->GetYaxis()->SetLabelSize(0.025);
    skippedGraph[i]->GetYaxis()->SetLabelOffset(0.05);
  }

  TCanvas *C_duplicate, *C_skipped;
  TString CId, CTitle;
  
  CId = "c_duplicate";
  CTitle = "Duplicate Event Fraction";
  C_duplicate = new TCanvas(CId, CTitle, 0, 0, 800, 600);
  C_duplicate->Divide(2,2);
  for ( int i = 0; i < 4; i++) {
    C_duplicate->cd(i+1);
    duplicateGraph[i]->Draw("AL");
  }

  CId = "c_skipped";
  CTitle = "Skipped Event Fraction";
  C_skipped = new TCanvas(CId, CTitle, 0, 0, 800, 600);
  C_skipped->Divide(2,2);
  for ( int i = 0; i < 4; i++) {
    C_skipped->cd(i+1);
    skippedGraph[i]->Draw("AL");
  }



  return;

}

