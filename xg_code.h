/*
 * $Id$
 */
#ifndef XG_CODE_H
#  define XG_CODE_H


#include "omp.h"
#include "gdal_priv.h"
#include "assert.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "xg_math_vector.h"

vector<vector<double> > ReadAsArray(GDALRasterBand* band, int startX, int startY, int XSize, int YSize){
    void * value = malloc(GDALGetDataTypeSize ( band -> GetRasterDataType() ) / 8  * XSize * YSize);
    assert(band -> RasterIO( GF_Read, startX, startY, XSize, YSize, value, XSize , YSize,  band -> GetRasterDataType(), 0, 0 ) == CE_None);
    vector<vector<double> > ans(YSize);
    switch(band -> GetRasterDataType() ){
        case GDT_Byte://Eight bit unsigned integer 
            for(int i = 0; i < YSize; ++i){
                ans[i].resize(XSize);
                for(int j = 0; j < XSize; ++j){
                    ans[i][j] = ((unsigned char*)value)[i*XSize + j];
                }
            }break;
        case GDT_UInt16://Sixteen bit unsigned integer
            for(int i = 0; i < YSize; ++i){
                ans[i].resize(XSize);
                for(int j = 0; j < XSize; ++j){
                    ans[i][j] = ((unsigned short int *)value)[i*XSize + j];
                }
            }break;
        case GDT_Int16://Sixteen bit signed integer
            for(int i = 0; i < YSize; ++i){
                ans[i].resize(XSize);
                for(int j = 0; j < XSize; ++j){
                    ans[i][j] = ((signed short int*)value)[i*XSize + j];
                }
            }break;
        case GDT_UInt32://Thirty two bit unsigned integer
            for(int i = 0; i < YSize; ++i){
                ans[i].resize(XSize);
                for(int j = 0; j < XSize; ++j){
                    ans[i][j] = ((unsigned int*)value)[i*XSize + j];
                }
            }break;
        case GDT_Int32://Thirty two bit signed integer
            for(int i = 0; i < YSize; ++i){
                ans[i].resize(XSize);
                for(int j = 0; j < XSize; ++j){
                    ans[i][j] = ((signed int*)value)[i*XSize + j] ;
                }
            }break;
        case GDT_Float32://Thirty two bit floating point
            for(int i = 0; i < YSize; ++i){
                ans[i].resize(XSize);
                for(int j = 0; j < XSize; ++j){
                    ans[i][j] = ((float*)value)[i*XSize + j];
                }
            }break;
        case GDT_Float64://Sixty four bit floating point
            for(int i = 0; i < YSize; ++i){
                ans[i].resize(XSize);
                for(int j = 0; j < XSize; ++j){
                    ans[i][j] = ((double*)value)[i*XSize + j];
                }
            }break;
        case GDT_CInt16://Complex Int16
            fprintf(stderr, "Don't support Complex Int16 now. \n");
            exit(1);
        case GDT_CInt32://Complex Int32
            fprintf(stderr, "Don't support Complex Int32 now. \n");
            exit(1);
        case GDT_CFloat32://Complex Float32
            fprintf(stderr, "Don't support Complex Float32 now. \n");
            exit(1);
        case GDT_CFloat64://Complex Float64 
            fprintf(stderr, "Don't support Complex Float64 now. \n");
            exit(1);
        default:
            fprintf(stderr, "Unknown or unspecified type. \n");
            exit(1);
    }
    free(value);
    return ans;
}

vector<double> ReadLine(GDALRasterBand* band, int lineid){
    int XSize = band -> GetXSize();
    void * value = malloc(GDALGetDataTypeSize ( band -> GetRasterDataType() ) / 8  * XSize);
    band -> RasterIO( GF_Read, 0, lineid, XSize, 1, value, XSize,1,  band -> GetRasterDataType(), 0, 0 );
    vector<double>ans;
    switch(band -> GetRasterDataType() ){
        case GDT_Byte://Eight bit unsigned integer 
            for(int i = 0; i < XSize; ++i){
                ans.push_back(((unsigned char*)value)[i]);
            }break;
        case GDT_UInt16://Sixteen bit unsigned integer
            for(int i = 0; i < XSize; ++i){
                ans.push_back(((unsigned short int *)value)[i]);
            }break;
        case GDT_Int16://Sixteen bit signed integer
            for(int i = 0; i < XSize; ++i){
                ans.push_back(((signed short int *)value)[i]);
            }break;
        case GDT_UInt32://Thirty two bit unsigned integer
            for(int i = 0; i < XSize; ++i){
                ans.push_back(((unsigned int *)value)[i]);
            }break;
        case GDT_Int32://Thirty two bit signed integer
            for(int i = 0; i < XSize; ++i){
                ans.push_back(((signed int*)value)[i]);
            }break;
        case GDT_Float32://Thirty two bit floating point
            for(int i = 0; i < XSize; ++i){
                ans.push_back(((float*)value)[i]);
            }break;
        case GDT_Float64://Sixty four bit floating point
            for(int i = 0; i < XSize; ++i){
                ans.push_back(((double*)value)[i]);
            }break;
        case GDT_CInt16://Complex Int16
            fprintf(stderr, "Don't support Complex Int16 now. \n");
            exit(1);
        case GDT_CInt32://Complex Int32
            fprintf(stderr, "Don't support Complex Int32 now. \n");
            exit(1);
        case GDT_CFloat32://Complex Float32
            fprintf(stderr, "Don't support Complex Float32 now. \n");
            exit(1);
        case GDT_CFloat64://Complex Float64 
            fprintf(stderr, "Don't support Complex Float64 now. \n");
            exit(1);
        default:
            fprintf(stderr, "Unknown or unspecified type. \n");
            exit(1);
    }
    free(value);
    return ans;
}

void WriteLine(GDALRasterBand * band, void * p,int lineID){
   assert(band -> RasterIO( GF_Write, 0, lineID, band->GetXSize(), 1, 
                      p, band->GetXSize(),1,  band -> GetRasterDataType(), 
                      0, 0 ) == CE_None);
}// write all the pixel value in line lineID

/*
 * 多波段到多波段
 */
template <typename T>
void Reduce(vector<T> (*rf)(const vector<double>&, const vector<double>&), const vector<string>& infs, const string & outfs, int outBandCount, GDALDataType outType, double DstNodata){
    /* read input GDALDataset */
    vector<GDALDataset* > SrcDats;
    //for(size_t i = 0; i < infs.size(); ++i){
#if __cplusplus >= 201103L
    for(auto inf : infs){
        SrcDats.push_back( (GDALDataset*) GDALOpenShared(inf.c_str(), GA_ReadOnly));
    }
    
    

#else
    for(size_t i = 0; i < infs.size(); ++i){
        SrcDats.push_back( (GDALDataset*) GDALOpenShared(infs[i].c_str(), GA_ReadOnly));
    }
#endif
    /* read input GDALRasterBand */
    vector<GDALRasterBand * > bands;
    //for(size_t i = 0; i < SrcDats.size(); ++i){
#if __cplusplus >= 201103L
    for(auto ds: SrcDats){
        int bandcount = ds ->GetRasterCount();
        for(int j = 0; j < bandcount; ++j){
            fprintf(stderr, "\rOpening file[%s] ...%d/%d", ds -> GetFileList() [0], j, bandcount);
            bands.push_back( ds ->GetRasterBand( j + 1) );
        }
        printf("\n");
    }
#else
    for(size_t i = 0; i < SrcDats.size(); ++i){
        int bandcount = SrcDats[i] ->GetRasterCount();
        for(int j = 0; j < bandcount; ++j){
            fprintf(stderr, "\rOpening file[%s] ...%d/%d", SrcDats[i] -> GetFileList() [0], j, bandcount);
            bands.push_back( SrcDats[i] ->GetRasterBand( j + 1) );
        }
        printf("\n");
    }
#endif
    fprintf(stderr, "%lu input bands ready...\n", bands.size());
    /* read input Nodata value of each bands*/
    vector<double> NaNs(bands.size());
    for(size_t i = 0; i < bands.size(); ++i){
        fprintf(stderr, "\rreading nodata value of band[%lu/%lu]...", i, bands.size());
        NaNs[i] = bands[i] -> GetNoDataValue();
    }
    assert(NaNs.size() == bands.size());
    int YSize = bands[0]->GetYSize(), XSize = bands[0]->GetXSize();
    /* 
     * prepare output raster bands 
     *  take input raster[0] as a tempelate
     */
    fprintf(stderr, "preparing output raster...\n");
    GDALDataset *tempdata;
    tempdata = (GDALDataset *) GDALOpenShared(infs[0].c_str(), GA_ReadOnly);
    assert(tempdata != NULL);
    GDALRasterBand *tempband = tempdata->GetRasterBand(1);
    assert(tempband != NULL);
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    assert(poDriver != NULL);
    char **	papszOptions = 0; //String list must be init as 0, or you'll get segementation fault
    /*压缩后arcgis打不开*/
    papszOptions = CSLSetNameValue( papszOptions, "COMPRESS", "DEFLATE" );
    papszOptions = CSLSetNameValue( papszOptions, "BIGTIFF", "YES" );
    GDALDataset *dataset = poDriver->Create( outfs.c_str(), XSize, YSize, outBandCount, outType , papszOptions);  
    dataset -> SetProjection(tempdata -> GetProjectionRef());
    double adfGeoTransform[6];
    tempdata -> GetGeoTransform( adfGeoTransform ) ;
    dataset -> SetGeoTransform( adfGeoTransform );
    vector<GDALRasterBand *>outBand(outBandCount);
    for(int i = 0; i < outBandCount; ++i){
        outBand[i] = (dataset -> GetRasterBand(i + 1));
        outBand[i] -> SetNoDataValue (DstNodata);
    }
    CPLFree(papszOptions);
    GDALClose(tempdata);
    
    /*
     * Start processing each line
     */
    fprintf(stderr, "start to process each line...");
#pragma omp parallel for schedule (dynamic)
    for(int line = 0; line < YSize; ++line){
        vector< vector<double> > lines(bands.size());
        //T* ans = (T*)malloc(sizeof(T) * XSize);
        vector<vector<T> > ans(outBandCount, vector<T>(XSize));

        fprintf(stderr, "\rprocessing %d of %d lines...", line , YSize);

        for(size_t k = 0; k < bands.size(); ++k){
#pragma omp critical 
            lines[k] = ReadLine(bands[k], line);
        }

        for(int i = 0; i < XSize; ++i){
            vector<double> inp;
            for(size_t k = 0; k < bands.size(); ++k){
                inp.push_back(lines[k][i]);
            }
            /*
             *  NaN is given to `rf` to handle
             */
            vector<T> tmp_mbands = rf(inp, NaNs);
            for(int s = 0; s < outBandCount; ++s){
                ans[s][i] = tmp_mbands[s];
            }
        }
        for(int s = 0; s < outBandCount; ++s){
#pragma omp critical 
            WriteLine(outBand[s], ans[s].data(), line);
        }
    }
#if __cplusplus >= 201103L
    for(auto ds: SrcDats){
        GDALClose(ds);
    }
#else
    for(size_t i = 0 ; i < SrcDats.size(); ++i){
        GDALClose(SrcDats[i]);
    }
#endif
    GDALClose(dataset);
}

/*
 * reduce multi-raster -> 1-band
 */
template <typename T>
void Reduce(T (*rf)(const vector<double>&, const vector<double>&), const vector<string>& infs, const string & outfs, GDALDataType outType, double DstNodata){
    /* read input GDALDataset */
    vector<GDALDataset* > SrcDats;
    //for(size_t i = 0; i < infs.size(); ++i){
#if __cplusplus >= 201103L
    for(auto inf : infs){
        SrcDats.push_back( (GDALDataset*) GDALOpenShared(inf.c_str(), GA_ReadOnly));
    }
#else
    for(size_t i = 0; i < infs.size(); ++i){
        SrcDats.push_back( (GDALDataset*) GDALOpenShared(infs[i].c_str(), GA_ReadOnly));
    }
#endif
    /* read input GDALRasterBand */
    vector<GDALRasterBand * > bands;
    //for(size_t i = 0; i < SrcDats.size(); ++i){
#if __cplusplus >= 201103L
    for(auto ds: SrcDats){
        int bandcount = ds ->GetRasterCount();
        for(int j = 0; j < bandcount; ++j){
            fprintf(stderr, "\rOpening file[%s] ...%d/%d", ds -> GetFileList() [0], j, bandcount);
            bands.push_back( ds ->GetRasterBand( j + 1) );
        }
        printf("\n");
    }
#else
    for(size_t i = 0; i < SrcDats.size(); ++i){
        int bandcount = SrcDats[i] ->GetRasterCount();
        for(int j = 0; j < bandcount; ++j){
            fprintf(stderr, "\rOpening file[%s] ...%d/%d", SrcDats[i] -> GetFileList() [0], j, bandcount);
            bands.push_back( SrcDats[i] ->GetRasterBand( j + 1) );
        }
        printf("\n");
    }
#endif
    fprintf(stderr, "%lu input bands ready...\n", bands.size());
    /* read input Nodata value of each bands*/
    vector<double> NaNs(bands.size());
    for(size_t i = 0; i < bands.size(); ++i){
        fprintf(stderr, "\rreading nodata value of band[%lu/%lu]...", i, bands.size());
        NaNs[i] = bands[i] -> GetNoDataValue();
    }
    assert(NaNs.size() == bands.size());
    int YSize = bands[0]->GetYSize(), XSize = bands[0]->GetXSize();
    /* 
     * prepare output raster bands 
     *  take input raster[0] as a tempelate
     */
    fprintf(stderr, "preparing output raster...\n");
    GDALDataset *tempdata;
    tempdata = (GDALDataset *) GDALOpenShared(infs[0].c_str(), GA_ReadOnly);
    assert(tempdata != NULL);
    GDALRasterBand *tempband = tempdata->GetRasterBand(1);
    assert(tempband != NULL);
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    assert(poDriver != NULL);
    char **	papszOptions = 0; //String list must be init as 0, or you'll get segementation fault
    //papszOptions = CSLSetNameValue( papszOptions, "COMPRESS", "DEFLATE" );
    //papszOptions = CSLSetNameValue( papszOptions, "BIGTIFF", "YES" );
    GDALDataset *dataset = poDriver->Create( outfs.c_str(), XSize, YSize, 1, outType , papszOptions);  
    dataset -> SetProjection(tempdata -> GetProjectionRef());
    double adfGeoTransform[6];
    tempdata -> GetGeoTransform( adfGeoTransform ) ;
    dataset -> SetGeoTransform( adfGeoTransform );
    GDALRasterBand *outBand = dataset -> GetRasterBand(1);
    outBand -> SetNoDataValue (DstNodata);
    CPLFree(papszOptions);
    GDALClose(tempdata);
    
    /*
     * Start processing each line
     */
    fprintf(stderr, "start to process each line...");
#pragma omp parallel for schedule (dynamic)
    for(int line = 0; line < YSize; ++line){
        vector< vector<double> > lines(bands.size());
        //T* ans = (T*)malloc(sizeof(T) * XSize);
        vector<T> ans(XSize);

        fprintf(stderr, "\rprocessing %d of %d lines...", line , YSize);

        for(size_t k = 0; k < bands.size(); ++k){
            #pragma omp critical 
            lines[k] = ReadLine(bands[k], line);
        }

        for(int i = 0; i < XSize; ++i){
            vector<double> inp;
            for(size_t k = 0; k < bands.size(); ++k){
                inp.push_back(lines[k][i]);
            }
            /*
             *  NaN is given to `rf` to handle
             */
            ans[i] = rf(inp, NaNs);
        }
        #pragma omp critical 
        WriteLine(outBand, ans.data(), line);
    }
#if __cplusplus >= 201103L
    for(auto ds: SrcDats){
        GDALClose(ds);
    }
#else
    for(size_t i = 0 ; i < SrcDats.size(); ++i){
        GDALClose(SrcDats[i]);
    }
#endif
    GDALClose(dataset);
}


/*
 * reduce multi-band -> value
 * considering position of pixels when calculating
 */
template <typename T>
T Reduce(
        void (*rf)(
            const vector<double>&,
            const vector<double>&,
            T&,
            double, //x_pos
            double, // y_pos
            double, // x_wd
            double // y_wd
            
        ),
        const vector<string>& infs,
        T init_value = T())
{
    /* read input GDALDataset */
    vector<GDALDataset* > SrcDats;
    //for(size_t i = 0; i < infs.size(); ++i){
#if __cplusplus >= 201103L
    for(auto inf : infs){
        SrcDats.push_back( (GDALDataset*) GDALOpenShared(inf.c_str(), GA_ReadOnly));
    }
#else
    for(size_t i = 0; i < infs.size(); ++i){
        SrcDats.push_back( (GDALDataset*) GDALOpenShared(infs[i].c_str(), GA_ReadOnly));
    }
#endif
    /* read input GDALRasterBand */
    vector<GDALRasterBand * > bands;
    //for(size_t i = 0; i < SrcDats.size(); ++i){
#if __cplusplus >= 201103L
    for(auto ds: SrcDats){
        int bandcount = ds ->GetRasterCount();
        for(int j = 0; j < bandcount; ++j){
            fprintf(stderr, "\rOpening file[%s] ...%d/%d", ds -> GetFileList() [0], j, bandcount);
            bands.push_back( ds ->GetRasterBand( j + 1) );
        }
        printf("\n");
    }
#else
    for(size_t i = 0; i < SrcDats.size(); ++i){
        int bandcount = SrcDats[i] ->GetRasterCount();
        for(int j = 0; j < bandcount; ++j){
            fprintf(stderr, "\rOpening file[%s] ...%d/%d", SrcDats[i] -> GetFileList() [0], j, bandcount);
            bands.push_back( SrcDats[i] ->GetRasterBand( j + 1) );
        }
        printf("\n");
    }
#endif
    fprintf(stderr, "%lu input bands ready...\n", bands.size());
    /* read input Nodata value of each bands*/
    vector<double> NaNs(bands.size());
    for(size_t i = 0; i < bands.size(); ++i){
        fprintf(stderr, "\rreading nodata value of band[%lu/%lu]...", i, bands.size());
        NaNs[i] = bands[i] -> GetNoDataValue();
    }
    assert(NaNs.size() == bands.size());
    int YSize = bands[0]->GetYSize(), XSize = bands[0]->GetXSize();
    
    /*
     * read projection info 
     */
    double adfGeoTransform[6];
    assert( SrcDats[0] -> GetGeoTransform( adfGeoTransform ) == CE_None );


    /*
     * Start processing each line
     */
    fprintf(stderr, "start to process each line...");
    for(int line = 0; line < YSize; ++line){
        vector< vector<double> > lines(bands.size());
        //T* ans = (T*)malloc(sizeof(T) * XSize);
        vector<T> ans(XSize);

        fprintf(stderr, "\rprocessing %d of %d lines...", line , YSize);

        for(size_t k = 0; k < bands.size(); ++k){
            lines[k] = ReadLine(bands[k], line);
        }

        for(int i = 0; i < XSize; ++i){

            double xpos = adfGeoTransform[0] + i*adfGeoTransform[1] + line*adfGeoTransform[2];
            double ypos = adfGeoTransform[3] + i*adfGeoTransform[4] + line*adfGeoTransform[5];


            vector<double> inp;
            for(size_t k = 0; k < bands.size(); ++k){
                inp.push_back(lines[k][i]);
            }
            /*
             *  NaN is given to `rf` to handle
             */
            rf(inp, NaNs, init_value, xpos, ypos, adfGeoTransform[1], adfGeoTransform[5]);
        }
    }
#if __cplusplus >= 201103L
    for(auto ds: SrcDats){
        GDALClose(ds);
    }
#else
    for(size_t i = 0 ; i < SrcDats.size(); ++i){
        GDALClose(SrcDats[i]);
    }
#endif
    return init_value;
}

/*
 * reduce multi-band -> value
 */
template <typename T>
T Reduce(void (*rf)(const vector<double>&, const vector<double>&, T&,int), const vector<string>& infs, T init_value = T()){
    /* read input GDALDataset */
    vector<GDALDataset* > SrcDats;
    //for(size_t i = 0; i < infs.size(); ++i){
#if __cplusplus >= 201103L
    for(auto inf : infs){
        SrcDats.push_back( (GDALDataset*) GDALOpenShared(inf.c_str(), GA_ReadOnly));
    }
#else
    for(size_t i = 0; i < infs.size(); ++i){
        SrcDats.push_back( (GDALDataset*) GDALOpenShared(infs[i].c_str(), GA_ReadOnly));
    }
#endif
    /* read input GDALRasterBand */
    vector<GDALRasterBand * > bands;
    //for(size_t i = 0; i < SrcDats.size(); ++i){
#if __cplusplus >= 201103L
    for(auto ds: SrcDats){
        int bandcount = ds ->GetRasterCount();
        for(int j = 0; j < bandcount; ++j){
            fprintf(stderr, "\rOpening file[%s] ...%d/%d", ds -> GetFileList() [0], j, bandcount);
            bands.push_back( ds ->GetRasterBand( j + 1) );
        }
        printf("\n");
    }
#else
    for(size_t i = 0; i < SrcDats.size(); ++i){
        int bandcount = SrcDats[i] ->GetRasterCount();
        for(int j = 0; j < bandcount; ++j){
            fprintf(stderr, "\rOpening file[%s] ...%d/%d", SrcDats[i] -> GetFileList() [0], j, bandcount);
            bands.push_back( SrcDats[i] ->GetRasterBand( j + 1) );
        }
        printf("\n");
    }
#endif
    fprintf(stderr, "%lu input bands ready...\n", bands.size());
    /* read input Nodata value of each bands*/
    vector<double> NaNs(bands.size());
    for(size_t i = 0; i < bands.size(); ++i){
        fprintf(stderr, "\rreading nodata value of band[%lu/%lu]...", i, bands.size());
        NaNs[i] = bands[i] -> GetNoDataValue();
    }
    assert(NaNs.size() == bands.size());
    int YSize = bands[0]->GetYSize(), XSize = bands[0]->GetXSize();
    
    /*
     * Start processing each line
     */
    fprintf(stderr, "start to process each line...");
    for(int line = 0; line < YSize; ++line){
        vector< vector<double> > lines(bands.size());
        //T* ans = (T*)malloc(sizeof(T) * XSize);
        vector<T> ans(XSize);

        fprintf(stderr, "\rprocessing %d of %d lines...", line , YSize);

        for(size_t k = 0; k < bands.size(); ++k){
            lines[k] = ReadLine(bands[k], line);
        }

        for(int i = 0; i < XSize; ++i){
            vector<double> inp;
            for(size_t k = 0; k < bands.size(); ++k){
                inp.push_back(lines[k][i]);
            }
            /*
             *  NaN is given to `rf` to handle
             */
            rf(inp, NaNs, init_value,line);
        }
    }
#if __cplusplus >= 201103L
    for(auto ds: SrcDats){
        GDALClose(ds);
    }
#else
    for(size_t i = 0 ; i < SrcDats.size(); ++i){
        GDALClose(SrcDats[i]);
    }
#endif
    return init_value;
}

// 降低分辨率重采样，而且输出图层是多波段的
template <typename T>
void Filter(vector<T> (*rf)(const vector<double>&,const vector<double>&, const double&), 
        //函数返回值是一个数组
        const string & infile, // only one input file
        const string & outfile,//only one output file
        int dstXSize, int dstYSize, //the output raster size
        GDALDataType outType, // Data type of output file
        int outBandCount, //输出影像波段数
        double DstNodata // output nodata
        ){
    /*
     * read input raster
     */
    GDALDataset* srcD =  (GDALDataset*) GDALOpenShared(infile.c_str(), GA_ReadOnly);
    assert(srcD != NULL);
    GDALRasterBand *srcB = srcD ->GetRasterBand(1);
    assert(srcB != NULL);
    double nodata = srcB -> GetNoDataValue();
    /*
     * read src/dst SIZE
     */
    long long srcXSize = srcB -> GetXSize();
    long long srcYSize = srcB -> GetYSize();
    //int dstXSize = atoi(argv[3]);
    //int dstYSize = atoi(argv[4]);
    double srcADF[6], dstADF[6];
    srcD -> GetGeoTransform( srcADF ) ;
    srcD -> GetGeoTransform( dstADF ) ;
    /*  
     *  0: top left x 
     *  1: w-e pixel resolution 
     *  2: rotation, 0 if image is "north up" 
     *  3: top left y 
     *  4: rotation, 0 if image is "north up" 
     *  5: n-s pixel resolution 

     *  Xp = [0] + P*[1] + L*[2];
     *  Yp = [3] + P*[4] + L*[5];
     */
    dstADF[1] = srcADF[1] * srcXSize / dstXSize;
    dstADF[2] = srcADF[2] * srcXSize / dstXSize;
    dstADF[4] = srcADF[4] * srcYSize / dstYSize;
    dstADF[5] = srcADF[5] * srcYSize / dstYSize;
    /* 
     * prepare output raster bands 
     *  take input raster[0] as a tempelate
     */
    fprintf(stderr, "preparing output raster...\n");
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    assert(poDriver != NULL);
    char **	papszOptions = 0; //String list must be init as 0, or you'll get segementation fault
    papszOptions = CSLSetNameValue( papszOptions, "COMPRESS", "DEFLATE" );
    papszOptions = CSLSetNameValue( papszOptions, "BIGTIFF", "YES" );
    GDALDataset *dstD = poDriver->Create( outfile.c_str(), dstXSize, dstYSize, outBandCount, outType, papszOptions);  
    dstD -> SetProjection(srcD -> GetProjectionRef());
    dstD -> SetGeoTransform( dstADF);
    //准备输出波段
    vector<GDALRasterBand *>dstB(outBandCount);
    for(int i = 0; i < outBandCount; ++i){
        dstB[i] = (dstD -> GetRasterBand(i + 1));
        dstB[i] -> SetNoDataValue (DstNodata);
    }
    CPLFree(papszOptions);

#pragma omp parallel for schedule (dynamic) 
    for (int dstI = 0; dstI < dstYSize; ++dstI){
        fprintf(stderr, "\rprocessing line %d/%d...", dstI, dstYSize);
        vector<vector<T> > dstLine(outBandCount);
        for(int i = 0; i < outBandCount; ++i){
            dstLine[i].resize(dstXSize);
        }
        for(int dstJ = 0; dstJ < dstXSize; ++dstJ){
            vector<double> pixels, weights;

            long long minI = dstI * srcYSize / dstYSize;
            long long  minJ = dstJ * srcXSize / dstXSize;
            long long  maxI = min(srcYSize, (dstI+1ll) * srcYSize / dstYSize+ ( ((dstI+1ll) * srcYSize % dstYSize) != 0));
            long long  maxJ = min(srcXSize, (dstJ+1ll) * srcXSize / dstXSize+ ( ((dstJ+1ll) * srcXSize % dstXSize) != 0));

            //if(maxJ < minJ){
                //out(srcXSize);
                //out(srcYSize);
                //out(dstXSize);
                //out(dstYSize);
                //out((dstJ+1) );
                //out((dstJ+1) * srcXSize );
                //out((dstJ+1) * srcXSize / dstXSize);
                //out(maxJ);
                //out(minJ);
                //out(dstI);
                //out(dstJ);
                //out((dstJ+1) * srcXSize / dstXSize);
            //}
            //fprintf(stderr, "start(%lld,%lld), size(%lld,%lld)\n", minJ, minI, maxJ - minJ , maxI-minI);
            vector<vector<double> >srcValue ;
            #pragma omp critical 
            srcValue = ReadAsArray(srcB, minJ, minI, maxJ - minJ , maxI-minI);

            for(int srcI = dstI * srcYSize / dstYSize; srcI < maxI; ++srcI){
                for(int srcJ = dstJ * srcXSize / dstXSize; srcJ < maxJ; ++srcJ){
                    double x = 1.0 * srcJ / srcXSize;
                    double y = 1.0 * srcI / srcYSize;
                    pixels.push_back(srcValue[srcI - minI][srcJ - minJ]);
                    weights.push_back( ( min( (dstJ+1.0) / dstXSize , x + 1.0 / srcXSize ) - max(1.0*dstJ / dstXSize , x ) ) *
                            ( min( (dstI+1.0) / dstYSize,  y + 1.0 / srcYSize ) - max(1.0*dstI / dstYSize,  y ) ) *
                            dstXSize * dstYSize);
                }
            }
            //printf("start(%d,%d), size(%d,%d)\n", minJ, minI, maxJ - minJ , maxI-minI);
            //out(pixels);
            //out(weights);
            //if(abs(sum(weights) - 1.0 )>1e-9){
            //out(pixels);
            //out(weights);
            //out(1.0 / srcXSize / srcYSize);
            //out(1.0 / dstXSize / dstYSize);
            //out(sum(weights));
            //}
            assert(abs(sum(weights) - 1.0 )<1e-9);
            //dstLine[dstJ] = max(pixels);
            vector<T> ans_mBands = rf(pixels, weights, nodata);
            for(int i = 0; i < outBandCount; ++i){
                dstLine[i][dstJ] = ans_mBands[i];
            }

            //if(max(pixels) < 1e-5 && max(pixels) > 0){
            //printf("start(%d,%d), size(%d,%d)\n", minJ, minI, maxJ - minJ , maxI-minI);
            //for(int srcI = dstI * srcYSize / dstYSize; srcI < maxI; ++srcI){
            //out(srcValue[srcI]);
            //}
            //out(pixels);
            //out(max(pixels));
            //out(dstLine[dstJ]);
            //}
            //assert(!(max(pixels) < 1e-5 && max(pixels) > 0));
            //assert(srcValue[0].size() > 0);
        }
        #pragma omp critical 
        for(int i = 0; i < outBandCount; ++i){
            WriteLine(dstB[i], dstLine[i].data(), dstI);
        }
    }
    GDALClose(srcD );
    GDALClose(dstD );
}

// 一对一重采样
template <typename T>
void Filter(T (*rf)(const vector<double>&,const vector<double>&, const double&), 
        // filter function take 3 arguments: input pixels, area/weights of each pixel, nodata value
        const string & infile, // only one input file
        const string & outfile,//only one output file
        int dstXSize, int dstYSize, //the output raster size
        GDALDataType outType, // Data type of output file
        double DstNodata // output nodata
        ){
    /*
     * read input raster
     */
    GDALDataset* srcD =  (GDALDataset*) GDALOpenShared(infile.c_str(), GA_ReadOnly);
    assert(srcD != NULL);
    GDALRasterBand *srcB = srcD ->GetRasterBand(1);
    assert(srcB != NULL);
    double nodata = srcB -> GetNoDataValue();
    /*
     * read src/dst SIZE
     */
    long long srcXSize = srcB -> GetXSize();
    long long srcYSize = srcB -> GetYSize();
    //int dstXSize = atoi(argv[3]);
    //int dstYSize = atoi(argv[4]);
    double srcADF[6], dstADF[6];
    srcD -> GetGeoTransform( srcADF ) ;
    srcD -> GetGeoTransform( dstADF ) ;
    /*  
     *  0: top left x 
     *  1: w-e pixel resolution 
     *  2: rotation, 0 if image is "north up" 
     *  3: top left y 
     *  4: rotation, 0 if image is "north up" 
     *  5: n-s pixel resolution 

     *  Xp = [0] + P*[1] + L*[2];
     *  Yp = [3] + P*[4] + L*[5];
     */
    dstADF[1] = srcADF[1] * srcXSize / dstXSize;
    dstADF[2] = srcADF[2] * srcXSize / dstXSize;
    dstADF[4] = srcADF[4] * srcYSize / dstYSize;
    dstADF[5] = srcADF[5] * srcYSize / dstYSize;
    /* 
     * prepare output raster bands 
     *  take input raster[0] as a tempelate
     */
    fprintf(stderr, "preparing output raster...\n");
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    assert(poDriver != NULL);
    char **	papszOptions = 0; //String list must be init as 0, or you'll get segementation fault
    papszOptions = CSLSetNameValue( papszOptions, "COMPRESS", "DEFLATE" );
    papszOptions = CSLSetNameValue( papszOptions, "BIGTIFF", "YES" );
    GDALDataset *dstD = poDriver->Create( outfile.c_str(), dstXSize, dstYSize, 1, outType, papszOptions);  
    dstD -> SetProjection(srcD -> GetProjectionRef());
    dstD -> SetGeoTransform( dstADF);
    GDALRasterBand *dstB = dstD -> GetRasterBand(1);
    dstB -> SetNoDataValue (DstNodata);
    CPLFree(papszOptions);

#pragma omp parallel for schedule (dynamic) 
    for (int dstI = 0; dstI < dstYSize; ++dstI){
        fprintf(stderr, "\rprocessing line %d/%d...", dstI, dstYSize);
        vector<T> dstLine(dstXSize);
        for(int dstJ = 0; dstJ < dstXSize; ++dstJ){
            vector<double> pixels, weights;

            long long minI = dstI * srcYSize / dstYSize;
            long long  minJ = dstJ * srcXSize / dstXSize;
            long long  maxI = min(srcYSize, (dstI+1ll) * srcYSize / dstYSize+ ( ((dstI+1ll) * srcYSize % dstYSize) != 0));
            long long  maxJ = min(srcXSize, (dstJ+1ll) * srcXSize / dstXSize+ ( ((dstJ+1ll) * srcXSize % dstXSize) != 0));

            //if(maxJ < minJ){
                //out(srcXSize);
                //out(srcYSize);
                //out(dstXSize);
                //out(dstYSize);
                //out((dstJ+1) );
                //out((dstJ+1) * srcXSize );
                //out((dstJ+1) * srcXSize / dstXSize);
                //out(maxJ);
                //out(minJ);
                //out(dstI);
                //out(dstJ);
                //out((dstJ+1) * srcXSize / dstXSize);
            //}
            //fprintf(stderr, "start(%lld,%lld), size(%lld,%lld)\n", minJ, minI, maxJ - minJ , maxI-minI);
            vector<vector<double> >srcValue ;
            #pragma omp critical 
            srcValue = ReadAsArray(srcB, minJ, minI, maxJ - minJ , maxI-minI);

            for(int srcI = dstI * srcYSize / dstYSize; srcI < maxI; ++srcI){
                for(int srcJ = dstJ * srcXSize / dstXSize; srcJ < maxJ; ++srcJ){
                    double x = 1.0 * srcJ / srcXSize;
                    double y = 1.0 * srcI / srcYSize;
                    pixels.push_back(srcValue[srcI - minI][srcJ - minJ]);
                    weights.push_back( ( min( (dstJ+1.0) / dstXSize , x + 1.0 / srcXSize ) - max(1.0*dstJ / dstXSize , x ) ) *
                            ( min( (dstI+1.0) / dstYSize,  y + 1.0 / srcYSize ) - max(1.0*dstI / dstYSize,  y ) ) *
                            dstXSize * dstYSize);
                }
            }
            //printf("start(%d,%d), size(%d,%d)\n", minJ, minI, maxJ - minJ , maxI-minI);
            //out(pixels);
            //out(weights);
            //if(abs(sum(weights) - 1.0 )>1e-9){
            //out(pixels);
            //out(weights);
            //out(1.0 / srcXSize / srcYSize);
            //out(1.0 / dstXSize / dstYSize);
            //out(sum(weights));
            //}
            assert(abs(sum(weights) - 1.0 )<1e-9);
            //dstLine[dstJ] = max(pixels);
            dstLine[dstJ] = rf(pixels, weights, nodata);

            //if(max(pixels) < 1e-5 && max(pixels) > 0){
            //printf("start(%d,%d), size(%d,%d)\n", minJ, minI, maxJ - minJ , maxI-minI);
            //for(int srcI = dstI * srcYSize / dstYSize; srcI < maxI; ++srcI){
            //out(srcValue[srcI]);
            //}
            //out(pixels);
            //out(max(pixels));
            //out(dstLine[dstJ]);
            //}
            //assert(!(max(pixels) < 1e-5 && max(pixels) > 0));
            //assert(srcValue[0].size() > 0);
        }
        #pragma omp critical 
        WriteLine(dstB, dstLine.data(), dstI);
    }
    GDALClose(srcD );
    GDALClose(dstD );
}



#endif /* ifndef XG_CODE_H */

