/*
 * $Id$
 */
#ifndef BPCH_H
#  define BPCH_H

#include <fstream>
#include "xg_math.h"
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "ogrsf_frmts.h"

/*
 * swap the bit of variable v
 * convert between Big-endian and Small-endian
 */
template<class T>
T bitswap(T v){
    T ans;
    unsigned char *dst = (unsigned char *)&ans;
    unsigned char *src = (unsigned char *)&v;
    size_t s = sizeof(T);
    for(size_t i = 0; i < s; ++i){
        dst[i] = src[s - i - 1];
    }
    return ans;
}
void show_str_list(char * s){
    int l = strlen(s);
    printf("[");
    for(int i = 0; i < l; ++i)
        printf("%c(%x),", s[i], s[i]);
    printf("]");
}
void str_cpy(char * s, const char *p){
    //printf("\ns = ");
    //show_str_list(s);
    int l = strlen(p);
    for(int i = 0; i < l; ++i) s[i] = p[i];
    //printf("\ns = ");
    //show_str_list(s);
}

class datablock{
    public:
        datablock():
            halfpolar(1), center180(1), tracer(1)
        {
            clear_str();
        }
        datablock(int d1, int d2, int d3, int d4, int d5, int d6):
            halfpolar(1), center180(1), tracer(1)
        {
            setdim(d1,d2,d3,d4,d5,d6);
            clear_str();
        }
        char modelname[21] , category[41] , unit[41] , reserved[41] ;
        float modelres_lon, modelres_lat;
        int halfpolar, center180, tracer, dim[6], skip;
        double tau0, tau1;
        /*
         * [i][j][k] --> [x][y][z]
         * x[0-dim1] --> [-180 --> 180]
         * y[0-dim0] --> [-90 --> 90]
         */
        vector< vector< vector<float> > > data;  // 
        vector<vector<float> >& operator[](size_t i){
            return data[i];
        }

        void toGeoTIFF(const char *s){
            float * dat_buf = (float*)malloc(sizeof(float) * dim[0] * dim[1]);

            /* set up TIF output */
            GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
            assert(poDriver != NULL);
            char **	papszOptions = 0; //String list must be init as 0, or you'll get segementation fault
            papszOptions = CSLSetNameValue( papszOptions, "BIGTIFF", "YES" );
            GDALDataset *dstD = poDriver->Create( s, dim[0], dim[1], dim[2], GDT_Float32, papszOptions);  

            /* set projection as EPSG 4326 */
            OGRSpatialReference oSRS;
            char *pszSRS_WKT = NULL;
            oSRS.importFromEPSG(4326);
            oSRS.exportToWkt( &pszSRS_WKT );
            dstD -> SetProjection(pszSRS_WKT);

            /* set resolution , start point*/
            double adfGeoTransform[6] = { modelres_lon * (dim[3] - 1) - 180, modelres_lon, 0, modelres_lat * (dim[4] - 1) + 90, 0, -modelres_lat  };
            dstD -> SetGeoTransform( adfGeoTransform);

            //output
            for(int b = 0; b < dim[2]; ++b){
                /* prepare data*/
                for(int x = 0; x < dim[0]; ++x){
                    for (int y = 0; y < dim[1]; ++y){
                        //if(abs(data[x][dim[1] - 1 - y][b] - 1 ) > 1){
                            //printf("data[%d][%d][%d] = %f\n", x, dim[1] - 1 - y, b, data[x][dim[1] - 1 - y][b]);
                        //}
                        dat_buf[y*dim[0] + x] = data[x][dim[1] - 1 - y][b]; //
                    }
                }

                assert(dstD -> GetRasterBand(b + 1) -> 
                       RasterIO( GF_Write, 0, 0, dim[0], dim[1], 
                            dat_buf, dim[0], dim[1],  GDT_Float32, 
                            0, 0 ) == CE_None);
            }


            GDALClose(dstD);
            CPLFree(papszOptions);
            free(pszSRS_WKT);
            free(dat_buf);
        }
        void showheader(){
            cout << "modelname=" << modelname<< endl;
            cout << "modelres(lon,lat)=" << modelres_lon << "," << modelres_lat<< endl;
            cout << "halfpolar=" << halfpolar<< endl;
            cout << "center180=" << center180<< endl;
            cout << "category=" << category<< endl;
            cout << "tracer=" << tracer<< endl;
            cout << "unit=" << unit<< endl;
            printf("tau0 = %f, tau1 = %f\n", tau0, tau1);
            cout << "reserved=" << reserved<< endl;
            cout << "dim=" ;
            for(int i = 0; i < 6; ++i)
                cout << dim[i] << " ";
            cout << endl;
            cout << "skip=" << skip<< endl;
        }
        void showdata(){
            printf("[\n");
            for(int i = 0; i < dim[0]; ++i){
                printf("[");
                for(int j = 0; j < dim[1]; ++j){
                    printf("[");
                    for(int k = 0; k < dim[2]; ++k){
                        cout << data[i][j][k] << ",";
                    }
                    printf("]");
                }
                printf("]\n");
            }
            printf("]\n");
        }
        void setdim(int d1, int d2, int d3, int d4, int d5, int d6){
            data.resize(d1, vector<vector<float> >(d2, vector<float>(d3) ) );
            dim[3] = d4;
            dim[4] = d5;
            dim[5] = d6;
            uniform();
        }
        void setvalue(float v){
            for(size_t i = 0; i < data.size(); ++i){
                for(size_t j = 0; j < data[i].size(); ++j){
                    for(size_t k = 0; k < data[i][j].size(); ++k){
                        data[i][j][k] = v;
			cout<<"i:"<<i<<"j:"<<j<<"k"<<k<<"v:"<<v<<endl;
                    }
                }
            }
        }
        void uniform(){
            dim[0] = dim[1] = dim[2] = data.size();
            if(dim[0] > 0){
                dim[1] = data[0].size();
                if(dim[1] > 0)
                    dim[2] = data[0][0].size();
            }
            skip = dim[0] * dim[1] * dim[2] * 4;
        }
    private:
        void clear_str(){
            for(int i = 0; i < 21; ++i) modelname[i] = ' '; modelname[20] = 0;
            for(int i = 0; i < 41; ++i) category[i] = ' ';  category[40]  = 0;
            for(int i = 0; i < 41; ++i) unit[i] = ' ';      unit[40]      = 0;
            for(int i = 0; i < 41; ++i) reserved[i] = ' ';  reserved[40]  = 0;
        }
};

class bpch{
    public:
        bpch(const char* _f, const char* _t)
        {
            clear_str();
            str_cpy(ftype, _f);
            str_cpy(toptitle, _t);
        }
        bpch(const string& f ){
            readF(f);
        }
        void showheader(){
            cout << "ftype=" << ftype << endl;
            cout << "toptitle=" << toptitle<< endl;
        }
        datablock& operator[](size_t index){
            return datablocks[index];
        }
        void push(const datablock& db){
            datablocks.push_back(db);
        }
        size_t size(){
            return datablocks.size();
        }
        void readF(const string& f ){
            //cout << "start Read " << f << endl;
            clear_str();
            filename = f;
            bpch_f.open (filename.c_str(), ios::in | ios::binary);
            if(!bpch_f){
                //printf("error!");
            }
            else{
                float f_buf;
                int i_buf;
                bpch_f.read((char*)&i_buf, 4);
                //cout << "ibuf = " << bitswap (i_buf) << endl;

                bpch_f.read(ftype, 40);
                //out(ftype);

                bpch_f.read((char*)&i_buf, 4);
                //cout << "ibuf = " << bitswap (i_buf) << endl;

                bpch_f.read((char*)&i_buf, 4);
                //cout << "ibuf = " << bitswap (i_buf) << endl;

                bpch_f.read(toptitle, 80);
                //out(toptitle);

                bpch_f.read((char*)&i_buf, 4);
                //cout << "ibuf = " << bitswap (i_buf) << endl;

                datablock _db;
                while(!bpch_f.eof()){
                    bpch_f.read((char*)&i_buf, 4);
                    //cout << "ibuf = " << bitswap (i_buf) << endl;

                    if(bpch_f.eof()){
                        break;
                    }

                    bpch_f.read(_db.modelname, 20);
                    //out(_db.modelname);

                    bpch_f.read((char*)&(_db.modelres_lon), 4);
                    //out(bitswap(_db.modelres_lon));
                    bpch_f.read((char*)&(_db.modelres_lat), 4);
                    //out(bitswap(_db.modelres_lat));
                    bpch_f.read((char*)&(_db.halfpolar), 4);
                    //out(bitswap(_db.halfpolar));
                    bpch_f.read((char*)&(_db.center180), 4);
                    //out(bitswap(_db.center180));

                    bpch_f.read((char*)&i_buf, 4);
                    //out(bitswap(i_buf));
                    bpch_f.read((char*)&i_buf, 4);
                    //out(bitswap(i_buf));

                    bpch_f.read(_db.category, 40);
                    //out(_db.category);

                    bpch_f.read((char*)&(_db.tracer), 4);
                    //out(bitswap(_db.tracer));

                    bpch_f.read(_db.unit, 40);
                    //out(_db.unit);

                    _db.tau0 = -1;
                    _db.tau1 = -1;
                    bpch_f.read((char*)&(_db.tau0), 8);
                    //out(bitswap(_db.tau0));
                    bpch_f.read((char*)&(_db.tau1), 8);
                    //out(bitswap(_db.tau1));

                    bpch_f.read(_db.reserved, 40);
                    //out(_db.reserved);

                    //cout << "dim=[" ;
                    for(int i = 0; i < 6; ++i){
                        _db.dim[i] = -1;
                        bpch_f.read((char*)&(_db.dim[i]), 4);
                        _db.dim[i] = bitswap(_db.dim[i]);
                        //cout << _db.dim[i] << ",";
                    }
                    //cout << endl;

                    bpch_f.read((char*)&(i_buf), 4);
                    //cout << "i_buf= " << bitswap(i_buf) << endl;
                    bpch_f.read((char*)&(i_buf), 4);
                    //cout << "i_buf= " << bitswap(i_buf) << endl;

                    bpch_f.read((char*)&(_db.skip), 4);
                    //out(bitswap(_db.skip));


                    _db.modelres_lon = bitswap(_db.modelres_lon);
                    _db.modelres_lat = bitswap(_db.modelres_lat);
                    _db.halfpolar = bitswap(_db.halfpolar);
                    _db.center180 = bitswap(_db.center180);
                    _db.tracer = bitswap(_db.tracer);
                    _db.skip = bitswap(_db.skip);
                    _db.tau0 = bitswap(_db.tau0);
                    _db.tau1 = bitswap(_db.tau1);

                    /* resize first */
                    _db.data.resize(_db.dim[0]);
                    for(int i = 0; i < _db.dim[0]; ++i){
                        _db.data[i].resize(_db.dim[1]);
                        for(int j = 0; j < _db.dim[1]; ++j){
                            _db.data[i][j].resize(_db.dim[2]);
                        }
                    }

                    /* read into data*/
                    for(int k = 0; k < _db.dim[2]; ++k){    // z direction
                        for(int j = 0; j < _db.dim[1]; ++j){ // y direction
                            for(int i = 0; i < _db.dim[0]; ++i){ /// x direction
                                if(bpch_f.eof()){
                                    printf("ERROR to END!");
                                    exit(1);
                                }
                                bpch_f.read((char*)&f_buf, 4);
                                _db.data[i][j][k] = bitswap(f_buf);
                            }
                        }
                    }
                    bpch_f.read((char*)&i_buf, 4);
                    //cout << "i_buf= " << bitswap(i_buf) << endl;
                    datablocks.push_back(_db);
                }
                bpch_f.close();
            }
        }
        void writeF(const string& f ){
            filename = f;
            bpch_f.open (filename.c_str(), ios::out | ios::binary);
            if(!bpch_f){
                printf("error!");
            }
            else{
                float f_buf;
                double d_buf;
                int i_buf;

                /***********************************/
                i_buf = bitswap(40); bpch_f.write((char*)&i_buf, 4);
                bpch_f.write(ftype, 40);
                i_buf = bitswap(40); bpch_f.write((char*)&i_buf, 4);
                /***********************************/

                /***********************************/
                i_buf = bitswap(80); bpch_f.write((char*)&i_buf, 4);
                bpch_f.write(toptitle, 80);
                i_buf = bitswap(80); bpch_f.write((char*)&i_buf, 4);
                /***********************************/

                for(vector<datablock>::iterator it = datablocks.begin(); it != datablocks.end(); ++it){
                    it->uniform();

                    /***********************************/
                    i_buf = bitswap(36); bpch_f.write((char*)&i_buf, 4);
                    bpch_f.write(it->modelname, 20);
                    f_buf = bitswap(it->modelres_lon);bpch_f.write((char*)&f_buf, 4);
                    f_buf = bitswap(it->modelres_lat);bpch_f.write((char*)&f_buf, 4);
                    i_buf = bitswap(it->halfpolar);   bpch_f.write((char*)&i_buf, 4);
                    i_buf = bitswap(it->center180);   bpch_f.write((char*)&i_buf, 4);
                    i_buf = bitswap(36); bpch_f.write((char*)&i_buf, 4);
                    /***********************************/


                    /***********************************/
                    i_buf = bitswap(168); bpch_f.write((char*)&i_buf, 4);
                    bpch_f.write(it->category, 40);

                    i_buf = bitswap(it->tracer); bpch_f.write((char*)&i_buf, 4);

                    bpch_f.write(it->unit, 40);

                    d_buf = bitswap(it->tau0); bpch_f.write((char*)&d_buf, 8);
                    d_buf = bitswap(it->tau1); bpch_f.write((char*)&d_buf, 8);

                    bpch_f.write(it->reserved, 40);

                    for(int i = 0; i < 6; ++i){
                        i_buf = bitswap(it->dim[i]); bpch_f.write((char*)&i_buf, 4);
                    }
                    i_buf = bitswap(it->skip + 8); bpch_f.write((char*)&(i_buf), 4);
                    i_buf = bitswap(168); bpch_f.write((char*)&(i_buf), 4);
                    /***********************************/

                    /***********************************/
                    i_buf = bitswap(it->skip); bpch_f.write((char*)&i_buf, 4);
                    for(int k = 0; k < it->dim[2]; ++k){
                        for(int j = 0; j < it->dim[1]; ++j){
                            for(int i = 0; i < it->dim[0]; ++i){
                                f_buf = bitswap(it->data[i][j][k]); bpch_f.write((char*)&f_buf, 4);
                            }
                        }
                    }
                    i_buf = bitswap(it->skip); bpch_f.write((char*)&i_buf, 4);
                    /***********************************/
                }
                bpch_f.close();
            }
        }
    private:
        string filename;
        fstream bpch_f;
        vector<datablock> datablocks;
        char ftype[41] , toptitle[81] ;
        void clear_str(){
            for(int i = 0; i < 41; ++i) ftype[i] = ' ';
            for(int i = 0; i < 81; ++i) toptitle[i] = ' ';
        }

};


#endif /* ifndef BPCH_H */

