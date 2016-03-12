
#include "points.h"

//DESCRIPTOR TYPE
typedef unsigned char   DTYPE;
//FEATURE LOCATION TYPE
typedef float LTYPE;

enum
{
//      READ_BUFFER_SIZE = 0x100000,
    SIFT_NAME= ('S'+ ('I'<<8)+('F'<<16)+('T'<<24)),
    MSER_NAME= ('M'+ ('S'<<8)+('E'<<16)+('R'<<24)),
    RECT_NAME= ('R'+ ('E'<<8)+('C'<<16)+('T'<<24)), 
    //SIFT_VERSION_2=('V'+('2'<<8)+('.'<<16)+('0'<<24)),
    //SIFT_VERSION_3=('V'+('3'<<8)+('.'<<16)+('0'<<24)),
    SIFT_VERSION_4=('V'+('4'<<8)+('.'<<16)+('0'<<24)),
    SIFT_VERSION_5=('V'+('5'<<8)+('.'<<16)+('0'<<24)),
    SIFT_EOF = (0xff+('E'<<8)+('O'<<16)+('F'<<24)),
};

static inline int IsValidVersionName(int value)
{
    return value == SIFT_VERSION_4 || value == SIFT_VERSION_5;
}

static inline int IsValidFeatureName(int value)
{
    return value == SIFT_NAME || value == MSER_NAME;
}

class LocationData: public Points<LTYPE>
{
public:
    int              _file_version;
public:
    //each feature has x,y, z, size, orientation
    //for rect feature, there is x, y, width, height
    //for eclips feature, there is u,v,a,b,c +z
public:
    void glPaint2D(int style);
    LocationData(int d, int n):Points<LTYPE>(d,n),  _file_version(0){   };
    LocationData( LocationData& sup, int index[], int n): Points<LTYPE>(sup,index,n), _file_version(0){};
    static void glPaintTexSIFT(const float * loc, float td[3]);
    static void glPaintTexFrontalVIP(const float * loc, float td[3]);
    static void glPaintSIFT(const LTYPE* loc);
    static void glPaintSIFTSQ(LTYPE*  loc);
    static void glPaintELIPS(LTYPE*  loc);
    static void glPaintRECT(LTYPE* loc);
    static void SetPaintColor(LTYPE sz);
};

class DescriptorData: public Points<DTYPE>
{
    //generally d is 128
public:
    DescriptorData(DescriptorData& sup, int index[], int n): Points<DTYPE>(sup,index,n){};
    DescriptorData(int d, int n):Points<DTYPE>(d,n){};
};
static float gSiftDisplayScale;
static int   gSiftVisualStyle;

LocationData *   _locData;
DescriptorData * _desData;
int              _npoint;
int              _updated;
void SetUpdated(){_updated = 1;}

void ResizeFeatureData(int npoint, int locDim = 5, int desDim = 128)
{
    if(npoint ==0)
    {
        if(_locData) delete _locData;
        if(_desData) delete _desData;
        _locData = NULL;
        _desData = NULL;
    }else
    {
        if(_locData)
            _locData->resize(locDim, npoint);
        else
            _locData = new LocationData(locDim, npoint);
        if(_desData)
            _desData->resize(desDim, npoint);
        else
            _desData = new DescriptorData(desDim, npoint);
        _locData->_file_version  = SIFT_VERSION_4;
    }
    _npoint = npoint;

}