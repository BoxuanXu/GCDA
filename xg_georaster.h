/*
 * $Id$
 */
#ifndef XG_GEORASTER_BASE_H
#  define XG_GEORASTER_BASE_H

#include <vector>
#include <iostream>
#include <map>
//#include "xg_debug.h"
#include "gdal_priv.h"
#include "assert.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"

using namespace std;
GDALDataset* OpenRasterData(string filename);
int GetXSize (GDALDataset*dataset);// get raster's X size
int GetYSize (GDALDataset*dataset);// get raster's Y size
void ReadLine(GDALDataset*dataset, void * p,int lineID, int bandID = 1);// read all the pixel value in line lineID, of band bandID. malloc by else
bool isNodata(GDALDataset* dataset, double v, int bandID=1);
double GetNodata(GDALDataset* dataset, int bandID=1);
int BandCount(GDALDataset*dataset);

class GeoRaster{
    public:
        GeoRaster(const string& filename, GDALAccess opentype = GA_ReadOnly);
        ~GeoRaster();
        void init(const string& filename, GDALAccess opentype);
        void close()const;
        int GetXSize ()const;// get raster's X size
        int GetYSize ()const;// get raster's Y size
        void * readLine(int lineID, int bandID = 0);// read all the pixel value in line lineID, of band bandID. malloc by GeoRaster
        void readLine(void * p,int lineID, int bandID = 0)const;// read all the pixel value in line lineID, of band bandID. malloc by else
        void EmptyTrash();
        bool isNodata(double v, int bandID=0)const;
        double GetNodata(int bandID=0)const;
        template<typename eT>
        map<eT, long long> GetHist()const;
        void setCoord(double upleftx, double uplefty, double cellwidth, double cellheght);
        GDALRasterBand * GetBand(int bandID)const;/* bandID is counted from 1 */
        int GetRasterCount(void) const;
    private:
        GDALDataset *dataset;
        string Filename;
        int bandcount;
        vector<GDALRasterBand *> bands;
        vector<void *> memtrash;
};

template <typename eT>
bool equal_assert(const GeoRaster & a, const GeoRaster & b, double eps = 1e-9){
    int xsizea = a.GetXSize(), ysizea = a.GetYSize();
    int xsizeb = b.GetXSize(), ysizeb = b.GetYSize();
    assert(xsizea == xsizeb);
    assert(ysizea == ysizeb);
    eT* linebufa = (eT*)malloc(sizeof(eT) * xsizea);
    eT* linebufb = (eT*)malloc(sizeof(eT) * xsizeb);
    cout<<" xsizea = " << xsizea;
    cout<<" ysizea = " << ysizea;
    for(int i = 0; i < ysizea; ++i){
        a.readLine(linebufa, i);
        b.readLine(linebufb, i);
        for(int j = 0; j < xsizea; ++j){
            assert(fabs(linebufa[j] - linebufb[j] ) < eps);
        }
    }
    free(linebufa);
    free(linebufb);
    return true;
}

template<typename eT>
map<eT, long long> GeoRaster::GetHist()const{
    map<eT, long long> ans;
    int xsize = GetXSize();
    int ysize = GetYSize();
    eT* linebuf = (eT*)malloc(sizeof(eT) * xsize);
    for(int i = 0; i < ysize; ++i){
        readLine(linebuf, i);
        for(int j = 0; j < xsize; ++j){
            if(!isNodata(linebuf[j]))
                ans[linebuf[j]] += 1;
        }
    }
    free(linebuf);
    return ans;
}

class GeoTiffWriter{
    public:
        GeoTiffWriter(const string& filename, const string& tempfile, GDALDataType eDataType, double NaN = -9999);
        ~GeoTiffWriter();
        void writeLine(void * p,int lineID);// write all the pixel value in line lineID
        GDALRasterBand *band;
    private:
        int XSize,YSize;
        GDALDataset *dataset;
};

GDALDataset* GeoTiffWriterOpen(const string& filename, const string& tempfile, double NaN = -9999);
void writeLine(GDALDataset* dataset, void * p,int lineID);// write all the pixel value in line lineID

bool isNodata(GDALRasterBand * band, double v);
double GetNodata(GDALRasterBand * band);
int GetXSize (GDALRasterBand * band);// get raster's X size
int GetYSize (GDALRasterBand * band);// get raster's Y size
void ReadLine(GDALRasterBand * band, void * p,int lineID);// read all the pixel value in line lineID, of band bandID. malloc by user
void WriteLine(GDALRasterBand * band, void * p,int lineID);// write all the pixel value in line lineID

#endif /* ifndef XG_GEORASTER_BASE_H */

