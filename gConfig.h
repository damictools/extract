#include <map>
#include <string>

class gConfig
{
    public:
        static gConfig& getInstance()
        {
            static gConfig instance; // Guaranteed to be destroyed.
                                     // Instantiated on first use.
            return instance;
        }
        bool readConfFile(const char* confFileName);
        
	int getStackSize(){return fStackSize;};
        
	float getExtSigma(const int ext);
	float getExtCal(const int ext);
	
	float getSeedThr()   {return fSeedThr;};
	float getAddThr()    {return fAddThr;};
	float getSkirtSize() {return fSkirtSize;};
	
	bool  getSaveTracks(){return fSaveTracks;};
	std::string getTracksCuts(){return fTracksCuts;};
	
        void printVariables();
        
    private:
        gConfig();
        
        int fStackSize;
	
        std::map<int, float> fExt2Sigma;
        float fDefaultSigma;
        
	std::map<int, float> fExt2Cal;
        float fDefaultCal;
	
        float fSeedThr;
        float fAddThr;
        int fSkirtSize;
        
        bool fSaveTracks;
        std::string fTracksCuts;
        
        // Dont forget to declare these two. You want to make sure they
        // are unaccessable otherwise you may accidently get copies of
        // your singleton appearing.
        gConfig(gConfig const&);              // Don't Implement
        void operator=(gConfig const&); // Don't implement
};