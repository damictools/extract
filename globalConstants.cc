#include "globalConstants.h"

int gVerbosity = true;

const int    kMaxLine = 10000;

const int kExtractedMask = -1000;
const int kEdgeFlag = 2;
const int kSatFlag = 4;
const double kSat  = 5e9;

const double kPrescan = 6;
const float kSatValue = 1e10;
const float kSatMargin = 0.9;

//colors
const char cyan[] = { 0x1b, '[', '1', ';', '3', '6', 'm', 0 };
const char magenta[] = { 0x1b, '[', '1', ';', '3', '5', 'm', 0 };
const char red[] = { 0x1b, '[', '1', ';', '3', '1', 'm', 0 };
const char green[] = { 0x1b, '[', '1', ';', '3', '2', 'm', 0 };
const char yellow[] = { 0x1b, '[', '1', ';', '3', '3', 'm', 0 };
const char blue[] = "\x1b[1;34m";

const char bold[] = "\x1b[1;39m";

const char whiteOnRed[]    = "\x1b[1;41m";
const char whiteOnGreen[]  = "\x1b[1;42m";
const char whiteOnPurple[] = "\x1b[1;45m";
const char whiteOnViolet[] = "\x1b[1;44m";
const char whiteOnBrown[]  = "\x1b[1;43m";
const char whiteOnGray[]   = "\x1b[1;47m";

const char normal[] = { 0x1b, '[', '0', ';', '3', '9', 'm', 0 };
