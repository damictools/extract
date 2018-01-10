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
#include <sys/resource.h>

#include "globalConstants.h"


#include "TFile.h"
#include "TNtuple.h"
#include "TObject.h"
#include "tinyxml2.h"
#include "gConfig.h"


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
    cout << "This program extracts hits and tracks, computes their relevant parameters\n";
    cout << "and saves them in a root file.\n";
    cout << normal;
  }
  cout << "==========================================================================\n";
  cout << yellow;
  cout << "\nUsage:\n";
  cout << "  "   << exeName << " <input file> <mask file> -o <output filename> \n\n";
  cout << "\nOptions:\n";
  cout << "  -q for quiet (no screen output)\n";
  cout << "  -b for computing and subtracting baseline from the image \n";
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

struct track_t{
//   TNtuple &nt;
  vector<Int_t>    xPix;
  vector<Int_t>    yPix;
  vector<Float_t>  adc;
  vector<Int_t>    level;
  
  Float_t eCore;
  Int_t nCore;
  Int_t flag;
  Int_t nSat;
  Int_t id;
  
  Int_t xMin;
  Int_t xMax;
  Int_t yMin;
  Int_t yMax;
  
  track_t() : eCore(0), nCore(0), flag(0), nSat(0), xMin(0), xMax(0), yMin(0), yMax(0) {};
  void fill(const Int_t &x , const Int_t &y, const Float_t &adcVal, const Int_t &l){ xPix.push_back(x); yPix.push_back(y); adc.push_back(adcVal); level.push_back(l); };
  
  void reset(){ xPix.clear(); yPix.clear(); adc.clear(); level.clear(); 
                eCore=0; nCore=0; flag=0; nSat=0; xMin=0; xMax=0; yMin=0; yMax=0; };
};

bool testForBugHit(const double* outArray, const int &i, const int &nX, const int &nY, track_t &hit, const char* mask, const double &kAddThr){
  
  int hitX = i%nX;
  int hitY = i/nX;
  
  //West
  if(hitX>0)
  {
    
    const double En = outArray[i-1];
    if (En<kBugEnergy && En>kExtractedMask)
    {
      hit.flag = hit.flag|kBugFlag;
      //cout << "W" << hitX << "\t" << hitY << "\t" << Ei << "\t" << En << endl;
      return true;
    }
  }

  //North-West
  if(hitX>0 && hitY<nY-1 /* && has kBugFlag */)
  {
    const double En = outArray[i+nX-1];
    if (En<kBugEnergy && En>kExtractedMask)
    {
      hit.flag = hit.flag|kBugFlag;
      //cout << "NW" << hitX << "\t" << hitY << "\t" << Ei << "\t" << En << endl;
      return true;
    }
  }
  
  //North
  if(hitY<nY-1 /* && has kBugFlag */)
  {
    const double En = outArray[i+nX];
    if (En<kBugEnergy && En>kExtractedMask)
    {
      hit.flag = hit.flag|kBugFlag;
      //cout << "N" << hitX << "\t" << hitY << "\t" << Ei << "\t" << En << endl;
      return true;
    }
  }
  
  //North-East
  if(hitX<nX-1 && hitY<nY-1 /* && has kBugFlag */)
  {
    const double En = outArray[i+nX+1];
    if (En<kBugEnergy && En>kExtractedMask)
    {
      hit.flag = hit.flag|kBugFlag;
      //cout << "NE" << hitX << "\t" << hitY << "\t" << Ei << "\t" << En << endl;
      return true;
    }
  }
  
  //East
  if(hitX<nX-1 /* && has kBugFlag */)
  {
    const double En = outArray[i+1];
    if (En<kBugEnergy && En>kExtractedMask)
    {
      hit.flag = hit.flag|kBugFlag;
      //cout << "E" << hitX << "\t" << hitY << "\t" << Ei << "\t" << En << endl;
      return true;
    }
  }

  //South-East
  if(hitX<nX-1 && hitY>0 /* && has kBugFlag */)
  {
    const double En = outArray[i-nX+1];
    if (En<kBugEnergy && En>kExtractedMask)
    {
      hit.flag = hit.flag|kBugFlag;
      //cout << "SE" << hitX << "\t" << hitY << "\t" << Ei << "\t" << En << endl;
      return true;
    }
  }
  
  //South
  if(hitY>0 /* && has kBugFlag */)
  {
    const double En = outArray[i-nX];
    if (En<kBugEnergy && En>kExtractedMask)
    {
      hit.flag = hit.flag|kBugFlag;
      //cout << "S" << hitX << "\t" << hitY << "\t" << Ei << "\t" << En << endl;
      return true;
    }
  }
  
  //South-West
  if(hitY>0 && hitX>0 /* && has kBugFlag */)
  {
    const double En = outArray[i-nX-1];
    if (En<kBugEnergy && En>kExtractedMask)
    {
      hit.flag = hit.flag|kBugFlag;
      //cout << "W" << hitX << "\t" << hitY << "\t" << Ei << "\t" << En << endl;
      return true;
    }
  }
  
  return false;
}

void extractTrack(double* outArray, const int &i, const int &nX, const int &nY, track_t &hit, const char* mask, const double &kAddThr, int recLevel = 0){
  
  if(recLevel>gHitMaxSize)
    return;
  
  int hitX = i%nX;
  int hitY = i/nX;
  const double &Ei = outArray[i];
  hit.flag = hit.flag|mask[i];
  if (Ei>kSat){
    if( (hit.flag & kBugFlag) == 0){
        testForBugHit(outArray, i, nX, nY, hit, mask, kAddThr);
    }
  }

  if(Ei>kSat){
    hit.flag = hit.flag|kSatFlag;
    hit.nSat += 1;

    //Check if it is a bugHit
    if( (hit.flag & kBugFlag) == 0){
      testForBugHit(outArray, i, nX, nY, hit, mask, kAddThr);
    }
  }
  
  if(Ei<kSat) hit.eCore += Ei;
  hit.fill(hitX, hitY, Ei, 0);
  
  outArray[i] = kExtractedMask-hit.id;
  
  //West
  if(hitX>0){
    const double &En = outArray[i-1];
    if(En>kAddThr){
      extractTrack(outArray, i-1, nX, nY, hit, mask, kAddThr, recLevel+1);
    }
  }
  
  //South
  if(hitY>0){
    const double &En = outArray[i-nX];
    if(En>kAddThr){
      extractTrack(outArray, i-nX, nX, nY, hit, mask, kAddThr, recLevel+1);
    }
  }
  
  //East
  if(hitX<nX-1){
    const double &En = outArray[i+1];
    if(En>kAddThr){
      extractTrack(outArray, i+1, nX, nY, hit, mask, kAddThr, recLevel+1);
    }
  }
  
  //North
  if(hitY<nY-1){
    const double &En = outArray[i+nX];
    if(En>kAddThr){
      extractTrack(outArray, i+nX, nX, nY, hit, mask, kAddThr, recLevel+1);
    }
  }
  
  return;
}


void addSkirt(const double* outArray, const int nX, const int nY, track_t &hit, const char* mask){
  
  gConfig &gc = gConfig::getInstance();
  const int kSkirtSize = gc.getSkirtSize();
  
  const int xMin = *( min_element(hit.xPix.begin(),hit.xPix.end()) );
  const int xMax = *( max_element(hit.xPix.begin(),hit.xPix.end()) );
  const int yMin = *( min_element(hit.yPix.begin(),hit.yPix.end()) );
  const int yMax = *( max_element(hit.yPix.begin(),hit.yPix.end()) );
  
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

  const int thisTrackMask = kExtractedMask-hit.id;
  for(int y=yScanMin;y<=yScanMax;++y){
    for(int x=xScanMin;x<=xScanMax;++x){
      
      const double &En = outArray[x+y*nX];
      if(En == thisTrackMask) continue;
      float dMin = kSkirtSize+1;
      bool doneSearch = false;
      
      for(int dxi=-kSkirtSize;dxi<=kSkirtSize;++dxi){
        int xi = x+dxi;
        if(xi<xScanMin || xi>xScanMax) continue;
        for(int dyi=-kSkirtSize;dyi<=kSkirtSize;++dyi){
          int yi = y+dyi;
          if(yi<yScanMin || yi>yScanMax) continue;
          
          const double &Ei = outArray[xi+yi*nX];
          if(Ei!=thisTrackMask) continue;
          
          const float dTot = sqrt(dxi*dxi+dyi*dyi);
          if(dMin>dTot){
            dMin = dTot;
            if(dMin==1){
              doneSearch=true;
              break;
            }
          }
        }
        if(doneSearch) break;  
      }
      
      if(dMin<kSkirtSize+1){
        hit.flag = hit.flag|mask[x+y*nX];
	if(En<kExtractedMask){
	  hit.flag = hit.flag|kMerFlag;
	}
	if(En>kSat){
	  hit.flag = hit.flag|kSatFlag;
	}
	hit.fill(x,y,En,(int)dMin);
      }

    }
  }
  
  return;
}

void computeHitParameters(const track_t &hit, const int &l, const double &kCal, Float_t *hitParam){
  
  double xSum   = 0;
  double ySum   = 0;
  double wSum   = 0;
  double adcSum = 0;
  int    nPix   = 0;
  for(unsigned int i=0;i<hit.xPix.size();++i){
    const double &adc = hit.adc[i]; 
    if(hit.level[i]<=l && adc>kExtractedMask){
      if(hit.adc[i]>0){
	xSum += hit.xPix[i]*adc;
	ySum += hit.yPix[i]*adc;
	wSum += adc;
      }
      if(adc<kSat) adcSum += adc;
      ++nPix;
    }
  }
  float xBary = xSum/wSum;
  float yBary = ySum/wSum;
  
  double x2Sum = 0;
  double y2Sum = 0;
  for(unsigned int i=0;i<hit.xPix.size();++i){
    const double &adc = hit.adc[i];
    if(hit.level[i]<=l && adc>kExtractedMask){
      if(hit.adc[i]>0){
	double dx = (hit.xPix[i] - xBary);
	double dy = (hit.yPix[i] - yBary);
	x2Sum   += dx*dx*adc;
	y2Sum   += dy*dy*adc;
      }
    }
  }
  float xVar = x2Sum/wSum;
  float yVar = y2Sum/wSum;
  
  hitParam[0] = adcSum*kCal;
  hitParam[1] = nPix;
  hitParam[2] = xBary;
  hitParam[3] = yBary;
  hitParam[4] = xVar;
  hitParam[5] = yVar;
}

void writeHit(const int &hitN, const track_t &hit, const double &kCal){
  ostringstream hitName;
  hitName << "hit_" << hitN;
  TNtuple nt(hitName.str().c_str(),hitName.str().c_str(),"x:y:E:level");
  
  for(unsigned int i=0;i<hit.xPix.size();++i)
    nt.Fill(hit.xPix[i], hit.yPix[i], hit.adc[i]*kCal, hit.level[i]);
  
  nt.Write();
}

struct hitTreeEntry_t{
  track_t hit;
  Int_t   runID;
  Int_t   ohdu;
  Int_t   expoStart;
  Int_t   nSavedPix;
  
  Float_t *hitParam;
  
  hitTreeEntry_t(const Int_t nPar): runID(-1), ohdu(-1), expoStart(-1),nSavedPix(0){ hitParam = new Float_t[nPar]; };
  ~hitTreeEntry_t(){ delete[] hitParam; };
};

void initHitTree(TTree &hitSumm, hitTreeEntry_t &evt ){
  hitSumm.Branch("runID",    &(evt.runID),    "runID/I");
  hitSumm.Branch("ohdu",     &(evt.ohdu),     "ohdu/I");
  hitSumm.Branch("expoStart",     &(evt.expoStart),     "expoStart/I");
  
  hitSumm.Branch("nSat", &(evt.hit.nSat), "nSat/I");
  hitSumm.Branch("flag", &(evt.hit.flag), "flag/I");
  hitSumm.Branch("xMin", &(evt.hit.xMin), "xMin/I");
  hitSumm.Branch("xMax", &(evt.hit.xMax), "xMax/I");
  hitSumm.Branch("yMin", &(evt.hit.yMin), "yMin/I");
  hitSumm.Branch("yMax", &(evt.hit.yMax), "yMax/I");
  
  gConfig &gc = gConfig::getInstance();
  const int kSkirtSize  = gc.getSkirtSize();
  const int kHitMaxSize = gc.getHitMaxSize();
  
  for(int l=0;l<=kSkirtSize;++l){
    for(int n=0;n<gNExtraTNtupleVars;++n){
      ostringstream varName;
      varName << gExtraTNtupleVars[n] << l;
      hitSumm.Branch(varName.str().c_str(),  &(evt.hitParam[n + gNExtraTNtupleVars*l]),  (varName.str()+"/F").c_str() );
    }
  }
  
  hitSumm.Branch("nSavedPix", &(evt.nSavedPix), "nSavedPix/I");
  evt.hit.xPix.reserve(kHitMaxSize+1);
  hitSumm.Branch("xPix",  &(evt.hit.xPix[0]), "xPix[nSavedPix]/I");
  evt.hit.yPix.reserve(kHitMaxSize+1);
  hitSumm.Branch("yPix",  &(evt.hit.yPix[0]), "yPix[nSavedPix]/I");
  evt.hit.level.reserve(kHitMaxSize+1);
  hitSumm.Branch("level", &(evt.hit.level[0]), "level[nSavedPix]/I");
  evt.hit.adc.reserve(kHitMaxSize+1);
  hitSumm.Branch("ePix", &(evt.hit.adc[0]), "ePix[nSavedPix]/F");
  
}

void refreshTreeAddresses(TTree &hitSumm, hitTreeEntry_t &evt)
{
  hitSumm.SetBranchAddress("xPix",  &(evt.hit.xPix[0]));
  hitSumm.SetBranchAddress("yPix",  &(evt.hit.yPix[0]));
  hitSumm.SetBranchAddress("level", &(evt.hit.level[0]));
  hitSumm.SetBranchAddress("ePix", &(evt.hit.adc[0]));
}


double stableMean(const double *v, const int &N){
  
  std::vector<double> temparray(v, v+N);
  std::sort(temparray.begin(), temparray.end());
  
  const int nMin = N/3;
  const int nMax = 2*N/3;
  
  double sum=0;
  for(int i=nMin;i<nMax;++i){
    sum+=temparray[i];
  }
  return sum/(nMax-nMin);
}

double subtractImageBaseline(double* outArray, const long totpix){
  double mean = stableMean(outArray, totpix);
  for (int i = 0; i < totpix; ++i)
  {
    outArray[i]-=mean;
  }
  return mean;
}


int searchForTracks(TFile *outF, TTree &hitSumm, hitTreeEntry_t &evt, double* outArray, const int runID, const int ohdu, const int expoStart, const long totpix, const int nX, const int nY, char* mask){

  gConfig &gc = gConfig::getInstance();
  const double kSeedThr   = gc.getExtSigma(ohdu) * gc.getSeedThr();
  const double kAddThr    = gc.getExtSigma(ohdu) * gc.getAddThr();
  const double kCal       = gc.getExtCal(ohdu);
  const int    kSkirtSize = gc.getSkirtSize();
  const bool  kSaveTracks = gc.getSaveTracks();
  const char *kTrackCuts  = gc.getTracksCuts().c_str();
  
  const int kNVars = gNBaseTNtupleVars + gNExtraTNtupleVars*(kSkirtSize+1);
  Float_t *ntVars = new Float_t[kNVars];
  if(mask == 0){
  	mask = new char[nX*nY]();
  } 
  
  outF->cd();
  
  unsigned int hitN = hitSumm.GetEntries();
  
  evt.runID = runID;
  evt.ohdu = ohdu;
  evt.expoStart = expoStart;
  track_t &hit = evt.hit;
  
  TTree hitSummAux("hitSummAux","hitSummAux");
  if(kSaveTracks){
    initHitTree(hitSummAux, evt);
    hitSummAux.SetCircular(1);
  }
  
  
  if(gVerbosity){
    cout << "\nProcessing runID " << runID << " ohdu " << ohdu << " -> sigma: " << gc.getExtSigma(ohdu) << ":\n";
  }
  for(long i=0;i<totpix;++i){
    
    if(outArray[i]>kSeedThr){
      
      evt.nSavedPix = 0;
      hit.reset();
      hit.id = hitN;
      ++hitN;
      
      extractTrack(outArray, i, nX, nY, hit, mask, kAddThr);
      addSkirt(outArray, nX, nY, hit, mask);
      
      for(int l=0;l<=kSkirtSize;++l){
        computeHitParameters( hit, l, kCal, &(evt.hitParam[gNExtraTNtupleVars*l]) );
      }
      
      if(kSaveTracks){
        hitSummAux.Fill();
        if( hitSummAux.GetEntries(kTrackCuts) == 1 ){
          evt.nSavedPix = hit.xPix.size();
          refreshTreeAddresses(hitSumm, evt);
          //writeHit(hitN, hit, kCal);
        }
      }
      
      hitSumm.Fill();
    }
    if(gVerbosity){
      if(i%1000 == 0) showProgress(i,totpix);
    }
  }
  delete[] ntVars;
  outF->cd();
  hitSumm.Write(hitSumm.GetName(),TObject::kOverwrite);
  if(gVerbosity){
    showProgress(1,1);
  }
  
  return 0;
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

int readMask(const char* maskName, vector <char*> &masks, const vector<int> &singleHdu){
  int status = 0;
  int nhdu = 0;
  double nulval = 0.;
  int anynul = 0;

  
  fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
  fits_open_file(&infptr, maskName, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  fits_get_num_hdus(infptr, &nhdu, &status);
  if (status != 0) return(status);
  
  /* check the extensions to process*/
  for(unsigned int i=0;i<singleHdu.size();++i){
    if(singleHdu[i] > nhdu){
      fits_close_file(infptr,  &status);
      cerr << red << "\nError: the file does not have the required HDU!\n\n" << normal;
      return -1000;
    }
  }
  
  vector<int> useHdu(singleHdu);
  if(singleHdu.size() == 0){
    for(int i=0;i<nhdu;++i){
      useHdu.push_back(i+1);
    }
  }
  const unsigned int nUseHdu=useHdu.size();
  
  
    
  for (unsigned int eN=0; eN<nUseHdu; ++eN)  /* Main loop through each extension */
  {
    const int n = useHdu[eN];
    
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
      masks.push_back(0);
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

void writeConfigTree(TFile *outF){
  
  gConfig &gc = gConfig::getInstance();
  
  outF->cd();
  TTree configTree("config","config");
  
  Float_t kSigma[101];
  for(int i=0;i<=100;++i){
    kSigma[i] = gc.getExtSigma(i);
  }
  configTree.Branch("sigma", &kSigma, "sigma[101]/F");
  
  Float_t kCal[101];
  for(int i=0;i<=100;++i){
    kCal[i] = gc.getExtCal(i);
  }
  configTree.Branch("sigma", &kCal, "cal[101]/F");
  
  Float_t kSeedThr = gc.getSeedThr();
  configTree.Branch("seedThr", &kSeedThr, "seedThr/F");
  
  Float_t kAddThr = gc.getAddThr();
  configTree.Branch("addThr", &kAddThr, "addThr/F");
  
  Int_t kSkirtSize  = gc.getSkirtSize();
  configTree.Branch("skirtSize", &kSkirtSize, "skirtSize/I");
  
  Int_t kStackSize  = gc.getStackSize();
  configTree.Branch("stackSize", &kStackSize, "stackSize/I");
  
  Int_t kHitMaxSize  = gc.getHitMaxSize();
  configTree.Branch("hitMaxSize", &kHitMaxSize, "hitMaxSize/I");
 
  Bool_t kSaveTracks  = gc.getSaveTracks();
  configTree.Branch("saveTracks", &kSaveTracks, "saveTracks/B");
  
  TString kTrackCuts  = gc.getTracksCuts();
  configTree.Branch("trackCuts", &kTrackCuts);

  configTree.Fill();
  configTree.Write();
}

int computeImage(const vector<string> &inFileList,const char *maskName, const char *outFile, const vector<int> &singleHdu, const int optFlag){
  int status = 0;
  double nulval = 0.;
  int anynul = 0;
  
  // Read mask
  bool noMask = false;
  vector <char*> masks;
  if(strlen(maskName) != 0) readMask(maskName, masks, singleHdu);
  else{
  	noMask = true;
  } 
  
  int nhdu = 0;
  const unsigned int nFiles  = inFileList.size();
  
  gConfig &gc = gConfig::getInstance();
  const int   kSkirtSize  = gc.getSkirtSize();
  
  TFile outRootFile(outFile, "RECREATE");
  
  writeConfigTree(&outRootFile);
  
  TTree hitSumm("hitSumm","hitSumm");
  hitTreeEntry_t evt((kSkirtSize+1)*gNExtraTNtupleVars);
  initHitTree(hitSumm, evt);
  
  for(unsigned int fn=0; fn < nFiles; ++fn){
    
    fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
    fits_open_file(&infptr, inFileList[fn].c_str(), READONLY, &status); /* Open the input file */
    if (status != 0) return(status);
    fits_get_num_hdus(infptr, &nhdu, &status);
    if (status != 0) return(status);
      
    /* check the extensions to process*/
    for(unsigned int i=0;i<singleHdu.size();++i){
      if(singleHdu[i] > nhdu){
	fits_close_file(infptr,  &status);
	cerr << red << "\nError: the file does not have the required HDU!\n\n" << normal;
	return -1000;
      }
    }
    
    vector<int> useHdu(singleHdu);
    if(singleHdu.size() == 0){
      for(int i=0;i<nhdu;++i){
	useHdu.push_back(i+1);
      }
    }
    const unsigned int nUseHdu=useHdu.size();
    
    
    for (unsigned int eN=0; eN<nUseHdu; ++eN)  /* Main loop through each extension */
    {
      
      const int n = useHdu[eN];

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
      double expoStart = 0;
      readCardValue(infptr, "EXPSTART", expoStart);

      if(optFlag&kCompBL){
        cout << "\n====================================\n";
        cout << "WARNING: subtracting image baseline.\n";
        double imBL = subtractImageBaseline(outArray, totpix);
        cout << "  Im BL = " << imBL << endl;
        cout << "Finished computing image baseline.\n";
        cout << "====================================\n\n";
      }

      if(noMask) searchForTracks(&outRootFile, hitSumm, evt, outArray, (int)runID, (int)ext, (int)expoStart, totpix, naxes[0], naxes[1], 0);
      else       searchForTracks(&outRootFile, hitSumm, evt, outArray, (int)runID, (int)ext, (int)expoStart, totpix, naxes[0], naxes[1], masks[eN]);
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
                           vector<int> &singleHdu, vector<string> &inFileList, string &maskFile, string &outFile, string &confFile, int &optFlag){
  
  if(argc == 1) return 1;
  optFlag = 0;
  inFileList.clear();
  singleHdu.clear();
  bool outFileFlag = false;
  bool maskFileFlag = false;
  int opt=0;
  while ( (opt = getopt(argc, argv, "bi:m:o:s:c:qhH?")) != -1) {
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
      case 'c':
        confFile = optarg;
        break;
      case 'b':
        optFlag |= kCompBL;
        break;
      case 'q':
        gVerbosity = 0;
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
    cerr << yellow << "\nMask filename missing. Will use empty mask\n" << normal;
    maskFile = "";
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
  
  int optFlag;
  string maskFile;
  string outFile;
  string confFile = "extractConfig.xml";
  vector<string> inFileList;
  vector<int> singleHdu;

  int returnCode = processCommandLineArgs( argc, argv, singleHdu, inFileList, maskFile, outFile, confFile, optFlag);
  if(returnCode!=0){
    if(returnCode == 1) printCopyHelp(argv[0],true);
    if(returnCode == 2) printCopyHelp(argv[0]);
    return returnCode;
  }

  /* Create configuration singleton and read configuration file */
  gConfig &gc = gConfig::getInstance();
  if(gc.readConfFile(confFile.c_str()) == false){
    return 1;
  }
  if(gVerbosity){
    cout << "\nConfig file: " << confFile << endl;
    gc.printVariables();
  }
  gHitMaxSize = gc.getHitMaxSize();

  
  /* Increase the stack size to be able to use a deeply
   * nested recursive function */
  const rlim_t kStackSize = gc.getStackSize() * 1024 * 1024;   // min stack size = gc.getStackSize() MB
  struct rlimit rl;
  int result;
  result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0){
      if (rl.rlim_cur < kStackSize){
          rl.rlim_cur = kStackSize;
          result = setrlimit(RLIMIT_STACK, &rl);
          if (result != 0){
	    cerr << "Error increasing the stack size: setrlimit returned result = " << result << endl;
          }
      }
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
  
  int status = computeImage( inFileList, maskFile.c_str(),  outFile.c_str(), singleHdu, optFlag);
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
