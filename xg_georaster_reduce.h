/*
 * $Id$
 */
#ifndef XG_GEORASTER_REDUCE_H
#  define XG_GEORASTER_REDUCE_H

#include "xg_georaster_base.h"
#include "omp.h"
template <typename T1, typename T2>
void rasterreduce(const vector<string>& GRnames,T2 (*rf)(const vector<T1>&) , const string & outfilename, GDALDataType eDataType, bool skipNodata = true, T2 dstNodata = -10000){
#ifdef DEBUG
    printf("start rasterreduce(string)...\n");
#endif
    vector < GeoRaster > GRs;
    for(size_t i = 0; i < GRnames.size(); ++i){
#ifdef DEBUG
        printf("%lu/%lu::Creating input Tif(%s)...\n",i,GRnames.size(), GRnames[i].c_str());
#endif
        GRs.push_back(GeoRaster(GRnames[i]));
    }
#ifdef DEBUG
    printf("Create input Tif finished...\n");
#endif
    GeoTiffWriter outRaster(outfilename, GRnames[0], eDataType, dstNodata);
#ifdef DEBUG
    printf("Create output Tif finished...\n");
#endif

    rasterreduce(GRs, rf, outRaster, skipNodata);

    for(size_t k = 0; k < GRnames.size(); ++k)
        GRs[k].close();
    //GDALClose(outRaster);
}

template <typename T1, typename T2>
void rasterreduce(const vector<GeoRaster>& GRs,T2 (*rf)(const vector<T1> & ) , GeoTiffWriter& outRaster, bool skipNodata = true){

    vector<GDALRasterBand * > bands;

    for(size_t k= 0; k < GRs.size(); ++k){
        int bandcount = GRs[k].GetRasterCount();
        for(int i = 0; i < bandcount; ++i){
            bands.push_back(GRs[k].GetBand(i+1));
        }
    }
    rasterreduce(bands, rf, outRaster.band, skipNodata);
}

template <typename T1, typename T2>
void rasterreduce(const vector<GDALRasterBand * >& bands,T2 (*rf)(const vector<T1> & ) , GDALRasterBand* outBand, bool skipNodata = true){
    T1** lines = (T1**)malloc(sizeof(T1*)*bands.size());
    int linecount = bands[0]->GetYSize();
    int xsize = bands[0]->GetXSize();

#ifdef DEBUG
    printf("getting sizes(%d)...\n",bands[0]->GetXSize());
    printf("xsize = %d\n", xsize);
    printf("malloc memory...\n");
#endif

    for(size_t k= 0; k < bands.size(); ++k){
        lines[k] = (T1*)malloc(sizeof(T1) * xsize);
        assert(lines[k] != NULL);
    }
    T2* ans = (T2*)malloc(sizeof(T2) * xsize);


    for(int line = 0; line < linecount; ++line){

        fprintf(stderr, "\rprocessing %d of %d lines...", line , linecount);

        for(size_t k = 0; k < bands.size(); ++k)
            ReadLine(bands[k], lines[k], line);

#pragma omp parallel for schedule (dynamic)
        for(int i = 0; i < xsize; ++i){
            vector<T1> inp;
            bool hasNaN = false;
            for(size_t k = 0; k < bands.size(); ++k){
                hasNaN = hasNaN || ( isNodata(bands[k], lines[k][i]) );
                inp.push_back(lines[k][i]);
            }
            /*
             * 这里对空值处理办法是：
             *      一个像元的序列中只要有一个空值，就将整个序列视为无效
             */
            if(skipNodata && hasNaN)
                ans[i] = GetNodata(bands[0]);
            else 
                ans[i] = rf(inp);
        }
        WriteLine(outBand, ans, line);
    }

#ifdef DEBUG
        printf("finished computing...\n");
        printf("memory freeing...\n");
#endif
    free(ans);
    for(size_t k= 0; k < bands.size(); ++k){
        free(lines[k]);
    }
    free(lines);
#ifdef DEBUG
        printf("memory freed...\n");
#endif
}

template <typename T1, typename T2, typename T>
void rasterreduce(const vector<string>& GRnames1, const vector<string>& GRnames2, 
        T (*rf)(const vector<T1>&, const vector<T2>&) ,
        const string & outfilename, GDALDataType eDataType, bool skipNodata = true){
    assert(GRnames1.size() == GRnames2.size());
    vector < GeoRaster > GRs1, GRs2;
    for(size_t i = 0; i < GRnames1.size(); ++i){
        GRs1.push_back(GeoRaster(GRnames1[i]));
        GRs2.push_back(GeoRaster(GRnames2[i]));
    }
    GeoTiffWriter outRaster(outfilename, GRnames1[0], eDataType, GRs1[0].GetNodata());

    rasterreduce(GRs1, GRs2, rf, outRaster, skipNodata);

    for(size_t k = 0; k < GRnames1.size(); ++k){
        GRs1[k].close();
        GRs2[k].close();
    }
    //GDALClose(outRaster);
}

template <typename T1, typename T2, typename T>
void rasterreduce(const vector<GeoRaster>& GRs1,const vector<GeoRaster>& GRs2,
        T (*rf)(const vector<T1> &, const vector<T2> & ) ,
        GeoTiffWriter& outRaster, bool skipNodata = true){
    assert(GRs1.size() == GRs2.size());

    vector<GDALRasterBand * > bands1, bands2;

    for(size_t k= 0; k < GRs1.size(); ++k){
        int band1count = GRs1[k].GetRasterCount();
        int band2count = GRs2[k].GetRasterCount();
        for(int i = 0; i < band1count; ++i){
            bands1.push_back(GRs1[k].GetBand(i+1));
        }
        for(int i = 0; i < band2count; ++i){
            bands2.push_back(GRs2[k].GetBand(i+1));
        }
    }
    rasterreduce(bands1, bands2, rf, outRaster.band, skipNodata);
}

template <typename T1, typename T2, typename T>
void rasterreduce(const vector<GDALRasterBand * >& bands1, const vector<GDALRasterBand * >& bands2,
        T (*rf)(const vector<T1> &, const vector<T2> &) , 
        GDALRasterBand* outBand, bool skipNodata = true){
    assert(bands1.size() == bands2.size());
    T1** lines1 = (T1**)malloc(sizeof(T1*)*bands1.size());
    T2** lines2 = (T2**)malloc(sizeof(T2*)*bands2.size());
    assert(bands1[0]->GetXSize() == bands2[0]->GetXSize());
    assert(bands1[0]->GetYSize() == bands2[0]->GetYSize());
    
    int linecount = bands1[0]->GetYSize();
    int xsize = bands1[0]->GetXSize();


    for(size_t k= 0; k < bands1.size(); ++k){
        lines1[k] = (T1*)malloc(sizeof(T1) * xsize);
        lines2[k] = (T2*)malloc(sizeof(T2) * xsize);
        assert(lines1[k] != NULL);
        assert(lines2[k] != NULL);
    }
    T* ans = (T*)malloc(sizeof(T) * xsize);


    for(int line = 0; line < linecount; ++line){

        fprintf(stderr, "\rprocessing %d of %d lines...", line , linecount);

        for(size_t k = 0; k < bands1.size(); ++k){
            ReadLine(bands1[k], lines1[k], line);
            ReadLine(bands2[k], lines2[k], line);
        }

#pragma omp parallel for schedule (dynamic)
        for(int i = 0; i < xsize; ++i){
            vector<T1> inp1;
            vector<T2> inp2;
            bool hasNaN = false;
            for(size_t k = 0; k < bands1.size(); ++k){
                hasNaN = hasNaN || ( isNodata(bands1[k], lines1[k][i]) );
                hasNaN = hasNaN || ( isNodata(bands2[k], lines2[k][i]) );
                inp1.push_back(lines1[k][i]);
                inp2.push_back(lines2[k][i]);
            }
            /*
             * 这里对空值处理办法是：
             *      一个像元的序列中只要有一个空值，就将整个序列视为无效
             */
            if(skipNodata && hasNaN)
                ans[i] = GetNodata(bands1[0]);
            else 
                ans[i] = rf(inp1, inp2);
        }
        WriteLine(outBand, ans, line);
    }

    free(ans);
    for(size_t k= 0; k < bands1.size(); ++k){
        free(lines1[k]);
        free(lines2[k]);
    }
    free(lines1);
    free(lines2);
}

#endif /* ifndef XG_GEORASTER_REDUCE_H */

