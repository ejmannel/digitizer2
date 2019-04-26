#include <iostream>
#include <cmath>
#include <cstdlib>
#include <map>

#include "Event.h"
#include "Eventiterator.h"
#include "fileEventiterator.h"
#include <packet_A.h>
#include <packetHeaders.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TDatime.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TObjString.h"

#include "digitizerrawdata_2.h"

using namespace std;

// program to look at digitzer raw data
// E.J. Mannel, BNL
// 16-Sept-2013

// Prototypes

void bookHist(void);
void plotHist(void);
void saveHist(void);
int digitizer(char *);

const int MASK = 0x0000FFFF;
const int MASK2 = 0x00003FFF;

// Global variables

TH1 *h_ADC[NMOD][MaxChannels];
TH1 *h_Mean[NMOD][MaxChannels];
TH1 *h_RMS[NMOD][MaxChannels];

int digitizer(int pkt,  char *file, int nsamples, int nmodules) {

  int buffer[32768];
  int nw;
  int iword;
  int n_events = 0;
  int n_recs = 0;
  int status;
  int bco_err = 0;
  int evtNumDelta[NMOD];;

  Packet *p;
  Event *evt;
  Eventiterator *it;

  TString filename;

  adcdata *adc;
  adc = new adcdata();
  
  int verbose = 0;
  FILE *errorFile;


  // a bunch of error counters, etc.
  int parity_err[2] = {0,0};
  int event_err[5] = {0, 0, 0, 0, 0};
  int event_err2[5] = {0, 0, 0, 0, 0};
  int event_err3[5] = {0, 0, 0, 0, 0};
  int evtNum_last[4] = {-1, -1, -1, -1};
  int FEMBCO_Last[4] = {0, 0, 0, 0};
  int BCO_Last = 0;
  int EventNum_last = 0;
  int FEMBCO_Delta[4] = {0, 0, 0, 0};
  int BCO_Delta = 0;
  int BCO_err[4] = {0, 0, 0, 0};
  int BCO_err2[4] = {0, 0, 0, 0};
  int parityLast[NMOD][2][2];

  errorFile = fopen("digitizerError.log", "a+");

  // first make a nice filename by making use of sthe string class

  TString rootfilename = file;
  TString delimiter = "/";
  TObjArray *tempString = rootfilename.Tokenize(delimiter);
  TIter iString(tempString);
  TObjString *os = 0;

  while ((os = (TObjString *)iString())){
    rootfilename = os->GetString().Data();
  }

  rootfilename.ReplaceAll(".PRDFF",".root");
  rootfilename.ReplaceAll(".prdf",".root");
  rootfilename.ReplaceAll(".evt",".root");

  rootfilename = "root/" + rootfilename;

  cout << "creating root file: " << rootfilename << endl;

  // now make the root file
  TFile f(rootfilename,"recreate");

  bookHist();

  TString branchString;

  TTree ADC("ADC", "Digitizer Raw Data");
  ADC.SetMaxTreeSize(50000000000LL);

  ADC.Branch("eventNum", &adc->eventNumber, "eventNumber/I");

  ADC.Branch("flag", &adc->flag, "flag/I");

  branchString = "detId/I";
  ADC.Branch("detId", &adc->detId, branchString);

  branchString = "module/I";
  ADC.Branch("module", &adc->module, branchString);

  branchString = "bco/I";
  ADC.Branch("bco", &adc->bco, branchString);

  branchString = "FEMHeader[";
  branchString += NMOD;
  branchString += "]/I";
  ADC.Branch("FEMHeader", adc->FEMHeader, branchString);

  branchString = "FEMSlot[";
  branchString += NMOD;
  branchString += "]/I";
  ADC.Branch("FEMSlot", adc->FEMSlot, branchString);

  branchString = "FEMEvtNum[";
  branchString += NMOD;
  branchString += "]/I";
  ADC.Branch("FEMEvtNum", adc->FEMEvtNum, branchString);

  branchString = "FEMBCO[";
  branchString += NMOD;
  branchString += "]/I";
  ADC.Branch("FEMBCO", adc->FEMBCO, branchString);

  ADC.Branch("nsamples", &adc->nsamples, "nsamples/I");

  branchString = "parity[";
  branchString += NMOD;
  branchString += "][2][2]/I";
  ADC.Branch("parity", &adc->parity, branchString);

  branchString = "chanADC[";
  branchString += NMOD*MaxChannels;
  branchString += "][";
  branchString += MaxSamples;
  branchString += "]/I";
  ADC.Branch("chanADC", adc->chanADC, branchString);

  branchString = "mean[";
  branchString += NMOD*MaxChannels;
  branchString += "]/F";
  ADC.Branch("mean", adc->mean, branchString);

  branchString = "rms[";
  branchString += NMOD*MaxChannels;
  branchString += "]/F";
  ADC.Branch("rms", adc->rms, branchString);
  
  filename = file;

  std::cout << "file: " << filename << std::endl;
  std::cout << "nmodules = " << nmodules << "  nsamples = " << nsamples << std::endl;

  fprintf(errorFile, "\n--------------------------------------------\n");
  fprintf(errorFile, "Processing file: %s\n", file);
  fprintf(errorFile, "nmodules = %d nsamples = %d\n\n", nmodules,nsamples);

  it =  new fileEventiterator (filename, status);

  if (status) {
    std::cout  << "Couldn't open input file " << filename << std::endl;
    delete it;
    exit(1);
  }
  
  std::cout << "Processing file" << std::endl;

  int imodule = 0;
  int errorInPkt = 0;
  int firstPacket = 1;

  for (int i = 0; i < NMOD; i++) {
    evtNumDelta[i] = 1;
  }
  
  while ( (evt = it->getNextEvent())) {

    n_events++;

    if ((n_events < 1000 && n_events%100 == 0) ||
	(n_events%1000 == 0 )) {
      std::cout<< "processing event " << n_events << std::endl;
    }
      
    for (Int_t raw_pkt = pkt; raw_pkt < pkt+1; raw_pkt++) {
      
      p = evt->getPacket( raw_pkt );
      
      if ( p ) {

      	n_recs++;

      	p->fillIntArray(buffer, 32767, &nw, "DATA");

	if (verbose) 
	  std::cout << "Packet word count: " << nw << std::endl; 
	
	for (int i = 0; i < NMOD*MaxChannels; i++) { 
	  adc->mean[i] = 0.;
	  adc->rms[i] = 0.;
	  for (int j = 0; j < MaxSamples; j++) { 
	    adc->chanADC[i][j] = 0;
	  }
	}

	// there are a lot of kludges in this since there
      	// are still issues wit the packet format 9-Dec-2016
	
      	// this extracts the number of samples based on packet format 
      	// as of 9-Dec-2016

	iword = 0;
	  
	adc->eventNumber = buffer[iword] & MASK;
	iword++;

	// iword == 1
	adc->flag = buffer[iword] & MASK;
	iword++;
	
	// iword == 2
	adc->detId = buffer[iword] & MASK;
	iword++;
	
	// iword == 3
	adc->module = buffer[iword] & MASK;
	iword++;
	
	// iword == 4
	adc->bco = buffer[iword] & MASK;

	iword++;
	
	if (verbose) { 
	  std::cout << std::endl << "New Event: " << adc->eventNumber 
		    << " : " << adc->bco << std::endl;
	}

	for (imodule = 0; imodule < nmodules; imodule++) {

	  if (verbose) { 
	    std::cout << "   Module " << imodule <<": ";
	  }

	  adc->nsamples = nsamples;
	  
	  adc->parity[imodule][0][0] = 0;
	  adc->parity[imodule][0][1] = 0;
	  adc->parity[imodule][1][0] = 0;
	  adc->parity[imodule][1][1] = 0;
	  	  
	  // iword == 5
	  adc->FEMHeader[imodule] = buffer[iword] & MASK;
	  adc->parity[imodule][1][iword%2] ^= buffer[iword] & MASK;

	  iword++;

	  // iword == 6
	  adc->FEMSlot[imodule] = buffer[iword] & MASK;

	  adc->parity[imodule][1][iword%2] ^= buffer[iword] & MASK;
	  iword++;

	  // iword == 7
	  adc->FEMEvtNum[imodule] = buffer[iword] & MASK;
	  adc->parity[imodule][1][iword%2] ^= buffer[iword] & MASK;
	  iword++;

	  // iword == 8
	  adc->FEMBCO[imodule] = buffer[iword] & MASK;
	  adc->parity[imodule][1][iword%2] ^= buffer[iword] & MASK;
	  iword++;
	
	  int ichan;
	  int isample;

	  if (verbose) { 
	    
	    std::cout << "Module event number: " << adc->FEMEvtNum[imodule] << std::endl;
	    std::cout << "             Module slot number: " 
		      << adc->FEMSlot[imodule] << std::endl;
	    std::cout << "             Module FEM BCO: " 
		      << adc->FEMBCO[imodule] << std::endl;
	  }	    

	  // iword == 9...
	  // int sequence, lastSequence;
	  float fadc;

	  for(ichan = 0; ichan < 64; ichan+=2) { 
	    isample = 0;

	    while (isample < nsamples) {

	      adc->chanADC[imodule*MaxChannels+ichan][isample] = 
		buffer[iword] & MASK;
	      adc->parity[imodule][1][iword%2] ^= buffer[iword] & MASK;
	    
	      adc->chanADC[imodule*MaxChannels+ichan+1][isample] = 
		buffer[iword+1] & MASK;
	      adc->parity[imodule][1][(iword+1)%2] ^= buffer[iword+1] & MASK;
	 
	      fadc = (float)(adc->chanADC[imodule*MaxChannels+ichan][isample] & MASK2);
	      adc->mean[imodule*MaxChannels+ichan] += fadc;
	      adc->rms[imodule*MaxChannels+ichan] += (fadc * fadc);
	      h_ADC[imodule][ichan]->Fill(fadc, 1.);

	      fadc = (float)(adc->chanADC[imodule*MaxChannels+ichan+1][isample] & MASK2);
	      adc->mean[imodule*MaxChannels+ichan+1] += fadc;
	      adc->rms[imodule*MaxChannels+ichan+1] += (fadc *fadc);
	      h_ADC[imodule][ichan+1]->Fill(fadc, 1.);	    
	    
	      isample++;
	      iword+=2;
	    }

	    // done looping over 2 channels of data.

	    adc->mean[imodule*MaxChannels+ichan] /= (float)isample;
	    adc->mean[imodule*MaxChannels+ichan+1] /= (float)isample;
	  
	    h_Mean[imodule][ichan]->Fill(adc->mean[imodule*MaxChannels+ichan], 1.);
	    h_Mean[imodule][ichan+1]->Fill(adc->mean[imodule*MaxChannels+ichan+1], 1.);
	  
	    adc->rms[imodule*MaxChannels+ichan] = 
	      sqrt( (adc->rms[imodule*MaxChannels+ichan]/ (float)isample) -
		    (adc->mean[imodule*MaxChannels+ichan]*
		     adc->mean[imodule*MaxChannels+ichan]));

	    adc->rms[imodule*MaxChannels+ichan+1] = 
	      sqrt( (adc->rms[imodule*MaxChannels+ichan+1]/(float)isample) -
		    (adc->mean[imodule*MaxChannels+ichan+1]*
		     adc->mean[imodule*MaxChannels+ichan+1]));
	  
	    h_RMS[imodule][ichan]->Fill(adc->rms[imodule*MaxChannels+ichan], 1.);
	    h_RMS[imodule][ichan+1]->Fill(adc->rms[imodule*MaxChannels+ichan+1], 1.);
	  
	  }

	  adc->parity[imodule][0][1] = buffer[iword] & MASK;
	  iword++;
	  adc->parity[imodule][0][0] = buffer[iword] & MASK;
	  iword++;

	  if (adc->parity[imodule][0][0] != adc->parity[imodule][1][0]) {
	    parity_err[0]++;
	    errorInPkt |= 0x0001;
	  }
	  if (adc->parity[imodule][0][1] != adc->parity[imodule][1][1]) {
	    parity_err[1]++;
	    errorInPkt |= 0x0001;
	  }

	  if ( iword >= nw && imodule+1 != nmodules) { 
	    std::cout << "Module count error: " << imodule << " - " 
	  	      << nmodules << std::endl;
	    break;
	  }

	} //end of module loop.

 	// Global packet analysis
	if (firstPacket == 1) { 
	  EventNum_last = adc->eventNumber;
	  BCO_Last = adc->bco;

	  for (imodule = 0; imodule < nmodules; imodule ++) { 
	    evtNum_last[imodule] = adc->FEMEvtNum[imodule];
	    FEMBCO_Last[imodule] = adc->FEMBCO[imodule];	  
	    parityLast[imodule][0][0]  = adc->parity[imodule][0][0];
	    parityLast[imodule][0][1]  = adc->parity[imodule][0][1];
	    parityLast[imodule][1][0]  = adc->parity[imodule][1][0];
	    parityLast[imodule][1][1]  = adc->parity[imodule][1][1];
	  }

	  firstPacket = 0;
	} else { 
	  if ((adc->eventNumber == 0 && EventNum_last != 0x0000FFFF) ||
	      (adc->eventNumber != 0 && adc->eventNumber - EventNum_last != 1)) {
	    event_err[0]++;
	    errorInPkt |= 0x0004;
	  }

	  if (adc->bco < BCO_Last)
	    BCO_Delta = adc->bco - BCO_Last + 65536;
	  else
	    BCO_Delta = adc->bco - BCO_Last;

 	  bco_err = 0;
	  for (imodule = 0; imodule < nmodules; imodule ++) {

	    if ((adc->eventNumber+evtNumDelta[imodule])%65536 != 
		adc->FEMEvtNum[imodule]) { 
	      event_err2[imodule]++;
	      errorInPkt |= 0x0008;
	      evtNumDelta[imodule] = (adc->FEMEvtNum[imodule] - 
				      adc->eventNumber)%65536;
	      std::cout << "evtNumDelta = "  << evtNumDelta[imodule] << std::endl;
	    }
	    if ((adc->FEMEvtNum[imodule] == 0 && evtNum_last[imodule]
		 != 0x0000FFFF) ||
		(adc->FEMEvtNum[imodule] != 0 && 
		 adc->FEMEvtNum[imodule] - evtNum_last[imodule] != 1)) {
	      event_err[imodule+1]++;
	      errorInPkt |= 0x0010;
	    }

	    if (adc->FEMBCO[imodule] >= FEMBCO_Last[imodule]) 
		FEMBCO_Delta[imodule] = 
		  adc->FEMBCO[imodule] - FEMBCO_Last[imodule];
	    else
		FEMBCO_Delta[imodule] = 
		  adc->FEMBCO[imodule] - FEMBCO_Last[imodule] + 65536; 

	    if (BCO_Delta != FEMBCO_Delta[imodule]) {
	      BCO_err[imodule]++;
	      bco_err++;
	    }
	  }	  

	  if ((nmodules > 1) && (adc->FEMEvtNum[0] != adc->FEMEvtNum[1])) {
	    event_err3[0]++;
	    errorInPkt |= 0x0040;
	  }
	  if ((nmodules > 2) && (adc->FEMEvtNum[0] != adc->FEMEvtNum[2])) { 
	    event_err3[1]++;
	    errorInPkt |= 0x0040;
	  }

	  if ((nmodules > 2) && adc->FEMEvtNum[1] != adc->FEMEvtNum[2]) {
	    event_err3[2]++;
	    errorInPkt |= 0x0040;
	  }

	  if ((nmodules > 1) && FEMBCO_Delta[0] != FEMBCO_Delta[1]) {
	    BCO_err2[0]++;
	    bco_err++;
	  }
	  if ((nmodules > 2) && FEMBCO_Delta[0] != FEMBCO_Delta[2]){
	    BCO_err2[1]++;
	    bco_err++;
	  }
	  if ((nmodules > 2) && FEMBCO_Delta[1] != FEMBCO_Delta[2]) {  
	    BCO_err2[2]++;
	    bco_err++;
	  }

	  if (bco_err) { 
	    errorInPkt |= 0x0020;
	  }

	  if (errorInPkt) {
	    std::cout << "File Event umber: " << n_events 
		      << " Events w/ PKT = 21351: " << n_recs 
		      << " XMIT Current: " << adc->eventNumber
		      << " XMIT Previous: " <<  EventNum_last << std::endl;
	    fprintf(errorFile, 
		    "File Event Number: %d  Events w/ PKT 21351 %d Current: %d Previous: %d \n",
		    n_events, n_recs, 
adc->eventNumber, EventNum_last);

	    printf("Error Code: %4.4X\n", errorInPkt);
	    fprintf(errorFile, "Error Code: %4.4X\n", errorInPkt);
	    if ((errorInPkt & 0x0001) == 0x0001) { 
	      std::cout << "Parity Error, ";
	      fprintf(errorFile, "Parity Error, ");
	    }
	    if ((errorInPkt & 0x0002) == 0x0002) {
	      std::cout << "Event Number Error, ";
	      fprintf(errorFile, "Event Number Error, ");
	    } 
	    if ((errorInPkt & 0x0004) == 0x0004) {
	      std::cout << "XMIT Number Increment Error, ";
	      fprintf(errorFile, "XMIT Number Error, ");
	    } 
	    if ((errorInPkt & 0x0008) == 0x0008) {
	      std::cout << "XMIT-FEM Event Number Error, ";
	      fprintf(errorFile, "Global FEM Number Error, ");
	    } 
	    if ((errorInPkt & 0x0010) == 0x0010) {
	      std::cout << "FEM Number Increment Error, ";
	      fprintf(errorFile, "FEM Event Number Increment Error, ");
	    } 
	    if ((errorInPkt & 0x0020) == 0x0020) {
	      std::cout << "BCO Error cnt =  " << bco_err << ", ";
	      fprintf(errorFile, "BCO Error Cnt = %d, ", bco_err);
	    }
	    if ((errorInPkt & 0x0040) == 0x0040) {
	      std::cout << "FEM Number Mismatch Error, ";
	      fprintf(errorFile, "FEM Event Number Mismatch Error, ");
	    } 

	    std::cout << std::endl;
	    fprintf(errorFile, "\n");
	    std::cout << "XMIT Counter Number: current= " << adc->eventNumber
		      << " previous= " <<  EventNum_last << std::endl;
	    fprintf(errorFile, "XMIT Counter Number current= %d previous= %d\n",
		    adc->eventNumber,EventNum_last);

	    for (imodule = 0; imodule < nmodules; imodule++) { 

	      std::cout << "FEM " << imodule << " Counter Number: current: " 
			<< adc->FEMEvtNum[imodule] 
			<< " previous: " <<evtNum_last[imodule] << std::endl;
	      fprintf(errorFile, "FEM %d Counter: current = %d previous = %d\n", 	
	      imodule, adc->FEMEvtNum[imodule], evtNum_last[imodule]);

	    }

	    std::cout << "BCO Error Count: " << bco_err << " Deltas: " 
		      << BCO_Delta << " ";
	    for (int i = 0; i < nmodules; i++) { 
	      std::cout << FEMBCO_Delta[i] << "  ";
	    }
	    std::cout << std::endl;

	    std::cout << "            Current BCO   BCO_Last  BCO_Delta  Difference" << std::endl; 
	    std::cout << "Global BCO: " << adc->bco << " " << BCO_Last 
		      << " " << BCO_Delta << " " << adc->bco - BCO_Last << std::endl;
	    for (int i = 0; i < nmodules; i++) { 
	      std::cout << "FEM " << i << ": " << adc->FEMBCO[i] << " BCO  "
			<< FEMBCO_Last[i] << "  " << FEMBCO_Delta[i] << " " 
			<< adc->FEMBCO[i] - FEMBCO_Last[i] << std::endl;
	    }

	    fprintf(errorFile, "\n      Current BCO_Last Delta Difference\n");
	    fprintf(errorFile, "Global BCO: %d %d %d %d\n", BCO_Last, adc->bco,
		    BCO_Delta, adc->bco - BCO_Last);
	    for (int i = 0; i < nmodules; i++) { 
	      fprintf(errorFile, "FEM %d BCO: %d %d %d %d\n",i, adc->FEMBCO[i],
		      FEMBCO_Last[i],FEMBCO_Delta[i], 
		      adc->FEMBCO[i] - FEMBCO_Last[i]);
	    }

	    std::cout << " Packet size: " << nw << std::endl;
	    fprintf(errorFile, "\nPacket size: %d\n", nw);

	    for (imodule = 0; imodule < nmodules; imodule++) { 

	      printf("Parity Module %d  slot: %d \n", 
		      imodule,  adc->FEMSlot[imodule]);
	      printf ("  Even Parity Words: %4.4X %4.4X %4.4X %4.4X %4.4X\n", 
		      adc->parity[imodule][0][0], adc->parity[imodule][1][0], 
		      adc->parity[imodule][0][0] ^ adc->parity[imodule][1][0],
		      parityLast[imodule][0][0], parityLast[imodule][1][0]);
	      printf ("   Odd Parity Words: %4.4X %4.4X %4.4X %4.4X %4.4X\n", 
		      adc->parity[imodule][0][1], adc->parity[imodule][1][1], 
		      adc->parity[imodule][0][1] ^ adc->parity[imodule][1][1],
		      parityLast[imodule][0][1], parityLast[imodule][1][1]);


	      fprintf(errorFile, "Parity Module %d  slot: %d \n", 
		      imodule,  adc->FEMSlot[imodule]);
	      fprintf(errorFile, "  Even Parity Words: %4.4X %4.4X %4.4X %4.4X %4.4X\n", 
		      adc->parity[imodule][0][0], adc->parity[imodule][1][0], 
		      adc->parity[imodule][0][0] ^ adc->parity[imodule][1][0],
		      parityLast[imodule][0][0], parityLast[imodule][1][0]);

	      fprintf(errorFile, "   Odd Parity Words: %4.4X %4.4X %4.4X %4.4X %4.4X\n", 
		      adc->parity[imodule][0][1], adc->parity[imodule][1][1], 
		      adc->parity[imodule][0][1] ^ adc->parity[imodule][1][1], 
		      parityLast[imodule][0][1], parityLast[imodule][1][1]);
   
	    }

	    fprintf(errorFile,"\n");
	    std::cout << std::endl;
	  }

	  // Update "previous" event information
	  EventNum_last = adc->eventNumber;
	  BCO_Last = adc->bco;
	  for (int imodule = 0; imodule < nmodules; imodule ++) {
	    evtNum_last[imodule] = adc->FEMEvtNum[imodule];	    
	    FEMBCO_Last[imodule] = adc->FEMBCO[imodule];
	    parityLast[imodule][0][0]  = adc->parity[imodule][0][0];
	    parityLast[imodule][0][1]  = adc->parity[imodule][0][1];
	    parityLast[imodule][1][0]  = adc->parity[imodule][1][0];
	    parityLast[imodule][1][1]  = adc->parity[imodule][1][1];
	  }

	} // End of global packet analysis 


	errorInPkt = 0;
	
      	ADC.Fill();

      } else { 
	std::cout<< "Packet: " << raw_pkt << " not found" << std::endl;
      } 
      
      delete p;
      
    } //end of packet loop

    delete evt;

  } // end of event loop
  
  if (n_events > 0) { 

    std::cout << "Building Tree" << std::endl;

    ADC.BuildIndex("eventNumber");
    ADC.Write();
  
    saveHist();
    f.Close();

  } 

  delete adc;
  delete it;
  
  cout << "Digitizer raw: Final n_recs: " << n_recs << endl;
  cout << "Digitizer raw: Total events: " << n_events << endl;

  std::cout << "Parity Error Count even/odd: " << parity_err[0] << " : " << parity_err[1]
	    << std::endl;

  std::cout << "Counter Errors: Number did not increment properly (Global: FEM 0: FEM 1: FEM2) " ;
  for (int i = 0; i < nmodules+1; i++) 
    std::cout<< ": " << event_err[i] << " ";
  std::cout << std::endl;
  std::cout << "Counter Errors: XMIT Counter != FEM Counter" ;
  for (int i = 0; i < nmodules; i++) 
    std::cout<< ": " << event_err2[i] << " ";
  std::cout << std::endl;
  std::cout << "Event Number Errors: FEM Counter Mismatch" ;
  for (int i = 0; i < nmodules; i++) 
    std::cout<< ": " << event_err3[i] << " ";
  std::cout << std::endl;

  std::cout << "BCO Errors: Global-FEM" ;
  
  for (int i = 0; i < nmodules; i++) 
    std::cout<< ": " << BCO_err[i] << " ";
  std::cout << std::endl;

  std::cout << "BCO Errors: FEM-FEM" ;
  for (int i = 0; i < nmodules; i++) 
    std::cout<< ": " << BCO_err2[i] << " ";
  std::cout << std::endl;

  fprintf(errorFile, "Digitizer raw: Final n_recs: %d\n", n_recs);
  fprintf(errorFile, "Digitizer raw: Total events: %d\n", n_events);

  fprintf(errorFile, "Parity Error Count even/odd: %d : %d\n", 
	  parity_err[0],parity_err[1]);
  
  fprintf(errorFile, "Event Number Errors: Run Number did not increment properly (Global: FEM 0: FEM 1: FEM2)");
  for (int i = 0; i < nmodules+1; i++) 
    fprintf(errorFile, ": %d ", event_err[i]);
  fprintf(errorFile,"\n");

  fprintf(errorFile, "Event Number Errors: Global Event Number != FEM Event Number");
  for (int i = 0; i < nmodules; i++) 
    fprintf(errorFile, ": %d ",event_err2[i]);
  fprintf(errorFile,"\n");

  fprintf(errorFile, "Event Number Errors: FEM Run Number Mismatch");
  for (int i = 0; i < nmodules; i++) 
    fprintf(errorFile, ": %d ",event_err3[i]);
  fprintf(errorFile,"\n");

  fprintf(errorFile, "BCO Errors: Global-FEM");
  
  for (int i = 0; i < nmodules; i++) 
    fprintf(errorFile,  ": %d ", BCO_err[i]);
  fprintf(errorFile,"\n");

  fprintf(errorFile, "BCO Errors: FEM-FEM");
  for (int i = 0; i < nmodules; i++) 
    fprintf(errorFile, ": %d ", BCO_err2[i]);
  fprintf(errorFile,"\n");


  return 0;  
}



void saveHist(void) { 
  
  for (int imodule = 0; imodule < NMOD; imodule++){ 
    for ( int i = 0; i < MaxChannels; i++) { 
      h_ADC[imodule][i]->Write();
      h_Mean[imodule][i]->Write();
      h_RMS[imodule][i]->Write();
    }
  }
  return;
}


void bookHist(void) {

  TString name, title;

  for (int imodule = 0; imodule < NMOD; imodule++) {
    for (int i = 0; i < MaxChannels; i++) { 
      name = "adc_";
      name += imodule; name += "-"; name += i;
      title = "ADC ";
      title += imodule; title += "-"; title += i;
      h_ADC[imodule][i] = new TH1D(name, title, 300, 1300., 1600.);
      
      name = "mean_";
      name += imodule; name += "-"; name += i;
      title = "mean ";
      title += imodule; title += "-"; title += i;
      h_Mean[imodule][i] = new TH1D(name, title, 16384, 0., 1684.);
      
      name = "rms_";
      name += imodule; name += "-"; name += i;
      title = "rms ";
      title += imodule; title += "-"; title += i;
      h_RMS[imodule][i] = new TH1D(name, title, 100, 0., 10.);
    }
  }

  return;

}

