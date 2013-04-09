#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include <iostream>
#include <sstream>

#include <sys/time.h>
#include <time.h>
#include <inttypes.h>
#include <fstream>
#include <unistd.h>
#include <getopt.h>    /* for getopt_long; standard getopt is in unistd.h */
#include <vector>
#include <algorithm>
#include <ctime>
#include <climits>
#include <cmath>
#include <iomanip>

#include "globalConstants.h"


#include "TFile.h"
#include "TNtuple.h"
#include "TObject.h"

using namespace std;

int deleteFile(const char *fileName){
  cout << yellow;
  cout << "Will overwrite: " << fileName << endl << endl;
  cout << normal;
  return unlink(fileName);
}

bool fileExist(const char *fileName){
  ifstream in(fileName,ios::in);
  
  if(in.fail()){
    //cout <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
    in.close();
    return false;
  }
  
  in.close();
  return true;
}

/*========================================================
  ASCII progress bar
==========================================================*/
void showProgress(unsigned int currEvent, unsigned int nEvent) {

  const int nProgWidth=50;

  if ( currEvent != 0 ) {
    for ( int i=0;i<nProgWidth+8;i++)
      cout << "\b";
  }

  double percent = (double) currEvent/ (double) nEvent;
  int nBars = (int) ( percent*nProgWidth );

  cout << " |";
  for ( int i=0;i<nBars-1;i++)
    cout << "=";
  if ( nBars>0 )
    cout << ">";
  for ( int i=nBars;i<nProgWidth;i++)
    cout << " ";
  cout << "| " << setw(3) << (int) (percent*100.) << "%";
  cout << flush;

}

void printCopyHelp(const char *exeName, bool printFullHelp=false){
  
  if(printFullHelp){
    cout << bold;
    cout << endl;
    cout << "This program computes overscan mean and subtracts it line by line.\n";
    cout << "It handles all the available HDUs. The HDUs in the output fit file\n";
    cout << "will be double (64bits) for 64bits images and float (32bis) in all\n";
    cout << "the other cases.\n";
    cout << "The card \"TRIMSEC\" must be present in the header to use this program.\n";
    cout << normal;
  }
  cout << "==========================================================================\n";
  cout << yellow;
  cout << "\nUsage:\n";
  cout << "  "   << exeName << " <input file> -o <output filename> \n\n";
  cout << "\nOptions:\n";
  cout << "  -v for verbosity\n";
  cout << "  -s <HDU number> for processing a single HDU \n\n";
  cout << normal;
  cout << blue;
  cout << "For any problems or bugs contact Javier Tiffenberg <javiert@fnal.gov>\n\n";
  cout << normal;
  cout << "==========================================================================\n\n";
}

string bitpix2TypeName(int bitpix){
  
  string typeName;
  switch(bitpix) {
      case BYTE_IMG:
          typeName = "BYTE(8 bits)";
          break;
      case SHORT_IMG:
          typeName = "SHORT(16 bits)";
          break;
      case LONG_IMG:
          typeName = "INT(32 bits)";
          break;
      case FLOAT_IMG:
          typeName = "FLOAT(32 bits)";
          break;
      case DOUBLE_IMG:
          typeName = "DOUBLE(64 bits)";
          break;
      default:
          typeName = "UNKNOWN";
  }
  return typeName;
}


const double kSigma = 18.1;
const double kSeedThr = kSigma*4;
const double kAddThr = kSigma*3;
const int kSkirtSize = 3;
const int kExtractedMask = -1000;
const int kEdgeFlag = 2;
const int kSatFlag = 4;
const double kSat  = 5e9;

struct track_t{
  TNtuple &nt;
  vector<int> xPix;
  vector<int> yPix;
  
  double eCore;
  int nCore;
  int flag;
  int nSat;
  
  int xMin;
  int xMax;
  int yMin;
  int yMax;
  
  track_t(TNtuple &n) : nt(n), eCore(0), nCore(0), flag(0), nSat(0), xMin(0), xMax(0), yMin(0), yMax(0) {};
};

void extractTrack(double* outArray, const int i, const int nX, const int nY, track_t &hit, const char* mask){
  
  int hitX = i%nX;
  int hitY = i/nX;
  const double &Ei = outArray[i];
  hit.flag = hit.flag|mask[i];
  
  if(Ei>kSat){
    hit.flag = hit.flag|kSatFlag;
    hit.nSat += 1;
  }
  hit.nt.Fill(hitX,hitY, Ei,0);
  hit.eCore += Ei;
  hit.xPix.push_back(hitX);
  hit.yPix.push_back(hitY);
  
  outArray[i] = kExtractedMask;
  
  //West
  if(hitX>0){
    const double &En = outArray[i-1];
    if(En>kAddThr){
      extractTrack(outArray, i-1, nX, nY, hit, mask);
    }
  }
  
  //South
  if(hitY>0){
    const double &En = outArray[i-nX];
    if(En>kAddThr){
      extractTrack(outArray, i-nX, nX, nY, hit, mask);
    }
  }
  
  //East
  if(hitX<nX-1){
    const double &En = outArray[i+1];
    if(En>kAddThr){
      extractTrack(outArray, i+1, nX, nY, hit, mask);
    }
  }

  //North
  if(hitY<nY-1){
    const double &En = outArray[i+nX];
    if(En>kAddThr){
      extractTrack(outArray, i+nX, nX, nY, hit, mask);
    }
  }
  
  return;
}


void addSkirt(const double* outArray, const int nX, const int nY, track_t &hit, const char* mask){
  const int xMin = hit.nt.GetMinimum("x");
  const int xMax = hit.nt.GetMaximum("x");
  const int yMin = hit.nt.GetMinimum("y");
  const int yMax = hit.nt.GetMaximum("y");
  
  const int xScanMin = (xMin-kSkirtSize >  0) ? xMin-kSkirtSize : 0;
  const int xScanMax = (xMax+kSkirtSize < nX) ? xMax+kSkirtSize : nX-1;
  const int yScanMin = (yMin-kSkirtSize >  0) ? yMin-kSkirtSize : 0;
  const int yScanMax = (yMax+kSkirtSize < nY) ? yMax+kSkirtSize : nY-1;
  
  if(xScanMin==0 || xScanMax==nX-1 || yScanMin==0 || yScanMax==nY-1) hit.flag = hit.flag|kEdgeFlag;
  
  const unsigned int nPix = hit.xPix.size();
  
  hit.xMin = xMin;
  hit.xMax = xMax;
  hit.yMin = yMin;
  hit.yMax = yMax;
  hit.nCore = nPix;
  
//   cout << hit.nt.GetName() << "\t" << xScanMin << "\t" << xScanMax << "\t" << yScanMin << "\t" << yScanMax << endl;
  
  for(int y=yScanMin;y<=yScanMax;++y){
    for(int x=xScanMin;x<=xScanMax;++x){
      int dMin = kSkirtSize+1;
      const double &En = outArray[x+y*nX];
      hit.flag = hit.flag|mask[x+y*nX];
      for(unsigned int i=0;i<nPix;++i){
	const int dx = abs(x-hit.xPix[i]);
	if( dx   >kSkirtSize) continue;
	const int dTot = abs(y-hit.yPix[i]) + dx;
	if( dTot >kSkirtSize) continue;
	
	if(En == kExtractedMask) continue;
	if(dMin>dTot){
	  dMin = dTot;  
	}
      }
      if(dMin >kSkirtSize) continue;
      else hit.nt.Fill(x,y,En,dMin);
    }
  }
  
  return;
}

int searchForTracks(TFile *outF, TNtuple &hitSumm, double* outArray, const int runID, const int ext, const int totpix, const int nX, const int nY, char* mask){
  
  outF->cd();
  
  unsigned int hitN = hitSumm.GetEntries();
  outF->cd("hits");
  if(gVerbosity){
    cout << "\nProcessing runID " << runID << " ext " << ext << ":\n";
  }
  for(unsigned int i=0;i<totpix;++i){
    
    if(outArray[i]>kSeedThr){
      
      ostringstream hitName;
      hitName << "hit_" << hitN;
      TNtuple nt(hitName.str().c_str(),hitName.str().c_str(),"x:y:E:level");
      track_t hit(nt);
      ++hitN;
      
      vector<int> xPix;
      vector<int> yPix;
      
      extractTrack(outArray, i, nX, nY, hit, mask);
      addSkirt(outArray, nX, nY, hit, mask);
      
      //nt.Write();
      hitSumm.Fill(runID, ext, hit.eCore*5.3363e-4, hit.nCore, hit.nSat, hit.flag, hit.xMin, hit.xMax, hit.yMin, hit.yMax);
    }
    if(gVerbosity){
      if(i%1000 == 0) showProgress(i,totpix);
    }
  }
  outF->cd();
  hitSumm.Write(hitSumm.GetName(),TObject::kOverwrite);
  if(gVerbosity){
    showProgress(1,1);
  }
}


bool readCardValue(fitsfile  *infptr, const char *keyName, double &value){
  
  int status = 0;
  char record[1024] = "";
  fits_read_card(infptr, keyName, record, &status);
  if(status==KEY_NO_EXIST){
    status=0;
    return false;
  }
  else{
    string sRec(record);
    size_t tPosEq = sRec.find("=");
    size_t tPosSl = sRec.find("/");
    istringstream recISS( sRec.substr(tPosEq+1, tPosSl-tPosEq-1) );
    recISS >> value;
    return true;
  }

}

int readMask(const char* maskName, vector <char*> &masks){
  int status = 0;
  int nhdu = 0;
  double nulval = 0.;
  int anynul = 0;

  
  fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
  fits_open_file(&infptr, maskName, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  fits_get_num_hdus(infptr, &nhdu, &status);
  if (status != 0) return(status);
    
  for (int n=1; n<=nhdu; ++n)  /* Main loop through each extension */
  {
     /* get input image dimensions and total number of pixels in image */
      int hdutype, bitpix, naxis = 0;
      long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      fits_movabs_hdu(infptr, n, &hdutype, &status);
      for (int i = 0; i < 9; ++i) naxes[i] = 1;
      fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
      long totpix = naxes[0] * naxes[1];
//       double bzero;
//       ffgky(infptr, TBYTE, "BZERO", &bzero, NULL, &status);
//       if (status){
// 	status = 0;
// 	bzero = 0.0;
//       }
      
      /* Don't try to process data if the hdu is empty */    
//       cout << (hdutype != IMAGE_HDU) << (naxis == 0) << (totpix == 0) << endl;
      if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
	continue;
      }
      
      char* maskArray = new char[totpix];
//       for(int i=0;i<totpix;++i) outArray[i] = 0;
      
      /* Open the input file */
      fits_movabs_hdu(infptr, n, &hdutype, &status);
      if (status != 0) return(status);
      int xMin=1;
      int xMax=naxes[0];
      int yMin=1;
      int yMax=naxes[1];
      
      /* Read the images as doubles, regardless of actual datatype. */
      long fpixel[2]={xMin,yMin};
      long lpixel[2]={xMax,yMax};
      long inc[2]={1,1};
      fits_read_subset(infptr, TBYTE, fpixel, lpixel, inc, &nulval, maskArray, &anynul, &status);
      if (status != 0){
	fits_report_error(stderr, status);
	return(status);
      }
      masks.push_back(maskArray);
  }
  fits_close_file(infptr,  &status);
  return status;
}


int computeImage(const vector<string> &inFileList,const char *maskName, const char *outFile, const vector<int> &singleHdu){
  int status = 0;
  double nulval = 0.;
  int anynul = 0;
  int single = 0;
  
  
  // Read mask
  vector <char*> masks;
  //const char maskName[]="mask6_sel_R.fits";
  readMask(maskName, masks);
  
  int nhdu = 0;
  const unsigned int nFiles  = inFileList.size();
  unsigned int nUseHdu = singleHdu.size();
  
  if (nUseHdu>0){
    single = 1; /* Copy only a single HDU if a specific extension was given */
  }
  
  TFile outRootFile(outFile, "RECREATE");
  outRootFile.mkdir("hits");
  TNtuple hitSumm("hitSumm","hitSumm","runID:ext:eCore:nCore:nSat:flag:xMin:xMax:yMin:yMax");
  for(unsigned int fn=0; fn < nFiles; ++fn){
    
    fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
    fits_open_file(&infptr, inFileList[fn].c_str(), READONLY, &status); /* Open the input file */
    if (status != 0) return(status);
    fits_get_num_hdus(infptr, &nhdu, &status);
    if (status != 0) return(status);
      
    if(nUseHdu==0) nUseHdu=nhdu;
    for (int eN=0; eN<nUseHdu; ++eN)  /* Main loop through each extension */
    {
      
      const unsigned int n = single? singleHdu[eN] : eN+1;
      
      if(n>nhdu){
	cerr << red << "\nError, requested hdu does not exist!\n\n" << normal;
	continue;
      }

      /* get input image dimensions and total number of pixels in image */
      int hdutype, bitpix, naxis = 0;
      long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      fits_movabs_hdu(infptr, n, &hdutype, &status);
      for (int i = 0; i < 9; ++i) naxes[i] = 1;
      fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
      long totpix = naxes[0] * naxes[1];
      double bzero;
      ffgky(infptr, TDOUBLE, "BZERO", &bzero, NULL, &status);
      if (status){
	status = 0;
	bzero = 0.0;
      }
      
      /* Don't try to process data if the hdu is empty */    
//       cout << (hdutype != IMAGE_HDU) << (naxis == 0) << (totpix == 0) << endl;
      if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
	if(single) break;
	continue;
      }
      
      double* outArray = new double[totpix];
      for(int i=0;i<totpix;++i) outArray[i] = 0;
      
      
      /* Open the input file */
      fits_movabs_hdu(infptr, n, &hdutype, &status);
      if (status != 0) return(status);
      int xMin=1;
      int xMax=naxes[0];
      int yMin=1;
      int yMax=naxes[1];
      
      /* Read the images as doubles, regardless of actual datatype. */
      long fpixel[2]={xMin,yMin};
      long lpixel[2]={xMax,yMax};
      long inc[2]={1,1};
      
      fits_read_subset(infptr, TDOUBLE, fpixel, lpixel, inc, &nulval, outArray, &anynul, &status);
      if (status != 0) return(status);
            
      double runID = 0;
      readCardValue(infptr, "RUNID", runID);
      double ext = n;
      readCardValue(infptr, "OHDU", ext);
      searchForTracks(&outRootFile, hitSumm, outArray, runID, ext,totpix, naxes[0], naxes[1], masks[n-1]);
      
      /* clean up */
      delete[] outArray;
    }

    /* Close the input file */
    fits_close_file(infptr,  &status);   
    
  }
  outRootFile.Close();
  
  return status;
}


void checkArch(){
  if(sizeof(float)*CHAR_BIT!=32 || sizeof(double)*CHAR_BIT!=64){
    cout << red;
    cout << "\n ========================================================================================\n";
    cout << "   WARNING: the size of the float and double variables is non-standard in this computer.\n";
    cout << "   The program may malfunction or produce incorrect results\n";
    cout << " ========================================================================================\n";
    cout << normal;
  }
}

int processCommandLineArgs(const int argc, char *argv[], 
                           vector<int> &singleHdu, vector<string> &inFileList, string &maskFile, string &outFile){
  
  if(argc == 1) return 1;
  inFileList.clear();
  singleHdu.clear();
  bool outFileFlag = false;
  bool maskFileFlag = false;
  int opt=0;
  while ( (opt = getopt(argc, argv, "i:m:o:s:vVhH?")) != -1) {
    switch (opt) {
      case 'm':
        if(!maskFileFlag){
          maskFile = optarg;
          maskFileFlag = true;
        }
        else{
          cerr << red << "\nError, can not set more than one mask file!\n\n" << normal;
          return 2;
        }
      break;
      case 'o':
        if(!outFileFlag){
          outFile = optarg;
          outFileFlag = true;
        }
        else{
          cerr << red << "\nError, can not set more than one output file!\n\n" << normal;
          return 2;
        }
        break;
      case 's':
        singleHdu.push_back(atoi(optarg));
        break;
      case 'V':
      case 'v':
        gVerbosity = 1;
        break;
      case 'h':
      case 'H':
      default: /* '?' */
        return 1;
    }
  }
  
  if(!outFileFlag){
    cerr << red << "\nError: output filename missing.\n" << normal;
    return 2;
  }
  
  if(!maskFileFlag){
    cerr << red << "\nError: mask filename missing.\n" << normal;
    return 2;
  }

  for(int i=optind; i<argc; ++i){
    inFileList.push_back(argv[i]);
    if(!fileExist(argv[i])){
      cout << red << "\nError reading input file: " << argv[i] <<"\nThe file doesn't exist!\n\n" << normal;
      return 1;
    }
  }
  
  if(inFileList.size()==0){
    cerr << red << "Error: no input file(s) provided!\n\n" << normal;
    return 1;
  }
  
  return 0;
}

int main(int argc, char *argv[])
{
  
  checkArch(); //Check the size of the double and float variables.
  
  time_t start,end;
  double dif;
  time (&start);
  
  string maskFile;
  string outFile;
  vector<string> inFileList;
  vector<int> singleHdu;

  int returnCode = processCommandLineArgs( argc, argv, singleHdu, inFileList, maskFile, outFile);
  if(returnCode!=0){
    if(returnCode == 1) printCopyHelp(argv[0],true);
    if(returnCode == 2) printCopyHelp(argv[0]);
    return returnCode;
  }
  
  if(gVerbosity){
    cout << bold << "\nWill read the following files:\n" << normal;
    for(unsigned int i=0; i<inFileList.size();++i) cout << "\t" << inFileList[i] << endl;
    if(singleHdu.size()>0){
      cout << bold << "\nAnd the following extension:" << normal << endl << "\t";
      for(unsigned int i=0; i<singleHdu.size();++i) cout << singleHdu[i] << ",";
      cout << "\b " << endl;
    }
    cout << bold << "\nThe output will be saved in the file:\n\t" << normal << outFile << endl;
  }
  
  int status = computeImage( inFileList, maskFile.c_str(),  outFile.c_str(), singleHdu);
  if (status != 0){ 
    fits_report_error(stderr, status);
    return status;
  }
  
  /* Report */
  time (&end);
  dif = difftime (end,start);
  if(gVerbosity) cout << green << "\nAll done!\n" << bold << "-> It took me " << dif << " seconds to do it!\n\n" << normal;

  return status;
}
