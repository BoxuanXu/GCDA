/*
 * $Id$
 */
#ifndef XG_GEORASTER_TRANSFORM_H
#  define XG_GEORASTER_TRANSFORM_H

#include "xg_georaster.h"
#include "omp.h"

template <typename T1, typename T2>
void raster_transform(const string & inputfile, T2 (*rf)(const T1 & ) , const string & outfilename, GDALDataType outDataType,  bool skipNodata = true){
    GeoRaster inp(inputfile);
    GeoTiffWriter outRaster(outfilename, inputfile, outDataType, inp.GetNodata());
    raster_transform(inp, rf, outRaster,skipNodata);
    inp.close();
}

template <typename T1, typename T2>
void raster_transform(const GeoRaster& GRs,T2 (*rf)(const T1 & ) , GeoTiffWriter& outRaster, bool skipNodata = true){

    int linecount = GRs.GetYSize();
    int xsize = GRs.GetXSize();

    T1* linebuf = (T1*)malloc(sizeof(T1) * xsize);
    assert(linebuf != NULL);
    T2* ans = (T2*)malloc(sizeof(T2) * xsize);
    assert(ans != NULL);


    for(int line = 0; line < linecount; ++line){

        fprintf(stderr, "\rprocessing %d of %d lines...", line , linecount);

        GRs.readLine(linebuf, line);

#pragma omp parallel for schedule (dynamic)
        for(int i = 0; i < xsize; ++i){
#ifdef DEBUG
            //printf(" processing i = %d by thread-%d\n", i, omp_get_thread_num());
#endif
            if(skipNodata && GRs.isNodata(linebuf[i]))
                ans[i] = GRs.GetNodata();
            else 
                ans[i] = rf(linebuf[i]);
        }
        outRaster. writeLine(ans, line);
    }

    free(ans);
    free(linebuf);
}

template <typename T1, typename T2, typename T3>
void raster_transform(const string & inputfile1, const string & inputfile2,
        T3 (*rf)(const T1 & , const T2&) , const string & outfilename, GDALDataType outDataType, bool skipNodata = true){
    GeoRaster inp1(inputfile1), inp2(inputfile2);
    GeoTiffWriter outRaster(outfilename, inputfile1, outDataType, inp1.GetNodata());
    raster_transform(inp1, inp2, rf, outRaster,skipNodata);
    inp1.close();
    inp2.close();
}

template <typename T1, typename T2, typename T3>
void raster_transform(const GeoRaster& GR1, const GeoRaster& GR2, T3 (*rf)(const T1 &, const T2 & ) , GeoTiffWriter& outRaster, bool skipNodata = true){

    assert(GR1.GetXSize() == GR2.GetXSize());
    assert(GR1.GetYSize() == GR2.GetYSize());
    int linecount = GR1.GetYSize();
    int xsize = GR1.GetXSize();

    T1* linebuf1 = (T1*)malloc(sizeof(T1) * xsize);
    T2* linebuf2 = (T2*)malloc(sizeof(T2) * xsize);
    assert(linebuf1 != NULL);
    assert(linebuf2 != NULL);
    T3* ans = (T3*)malloc(sizeof(T3) * xsize);
    assert(ans != NULL);


    for(int line = 0; line < linecount; ++line){

        fprintf(stderr, "\rprocessing %d of %d lines...", line , linecount);

        GR1.readLine(linebuf1, line);
        GR2.readLine(linebuf2, line);

#pragma omp parallel for schedule (dynamic)
        for(int i = 0; i < xsize; ++i){
#ifdef DEBUG
            printf(" processing i = %d by thread-%d\n", i, omp_get_thread_num());
#endif
            if(skipNodata && GR1.isNodata(linebuf1[i]))
                ans[i] = GR1.GetNodata();
            else if(skipNodata && GR2.isNodata(linebuf2[i]))
                ans[i] = GR2.GetNodata();
            else 
                ans[i] = rf(linebuf1[i], linebuf2[i]);
        }
        outRaster. writeLine(ans, line);
    }

    free(ans);
    free(linebuf1);
    free(linebuf2);
}

#endif /* ifndef XG_GEORASTER_TRANSFORM_H */

