//
// Created by uno on 2019/9/20.
//

#ifndef LATITUDELONGITUDE_2_ENU_COORDCOVERT_H
#define LATITUDELONGITUDE_2_ENU_COORDCOVERT_H
#include "stdlib.h"
bool xyz_center_flag = 0;
bool lbh_center_flag = 0;

//WGS-84椭球体参数
const double a=6378137.0;//长半轴，单位m
const double b=6356752.314245179;//短半轴，单位m
const double flattening=1/298.257223563;//扁率,0.003352811
const double delta=0.0000001;
typedef struct tagCRDCARTESIAN{
    double x;
    double y;
    double z;
}CRDCARTESIAN;

typedef CRDCARTESIAN *PCRDCARTESIAN;
//笛卡尔坐标系

typedef struct tagCRDGEODETIC{
    double longitude;
    double latitude;
    double height;
}CRDGEODETIC;

typedef CRDGEODETIC *PCRDGEODETIC;
//大地坐标系

typedef struct tagCRDTOPOCENTRIC{
    double northing;
    double easting;
    double upping;
}CRDTOPOCENTRIC;

typedef CRDTOPOCENTRIC *PCRDTOPOCENTRIC;
//站心地平坐标系（线坐标形式）

typedef struct tagCRDTOPOCENTRICPOLAR{
    double range;
    double azimuth;
    double elevation;
}CRDTOPOCENTRICPOLAR;

typedef CRDTOPOCENTRICPOLAR *PCRDTOPOCENTRICPOLAR;
//站心地平坐标系（极坐标形式）

//由笛卡尔坐标转换为大地坐标
void CartesianToGeodetic (PCRDGEODETIC pcg, PCRDCARTESIAN pcc,
                          double dSemiMajorAxis, double dFlattening);
//pcg：指向所转换出的大地坐标的指针；
//pcc：指向待转换的笛卡尔坐标的指针；
//dSemiMajorAxis：参考椭球的长半轴；
//dFlattening：参考椭球的扁率。

//由大地坐标转换为笛卡尔坐标
void GeodeticToCartesian (PCRDCARTESIAN pcc, PCRDGEODETIC pcg,
                          double dSemiMajorAxis, double dFlattening);
//pcc：指向所转换出的笛卡尔坐标的指针；
//pcg：指向待转换的大地坐标的指针；
//dSemiMajorAxis：参考椭球的长半轴；
//dFlattening：参考椭球的扁率。

//由笛卡尔坐标转换为站心地平坐标
void topocentric(PCRDTOPOCENTRIC pct, PCRDCARTESIAN pcc, PCRDCARTESIAN pccCenter, PCRDGEODETIC pd);
//pct：指向所转换出的站心地平坐标的指针；
//pcc：指向待转换的笛卡尔坐标的指针；
//pccCenter：指向站心的笛卡尔坐标的指针；
//dSemiMajorAxis：参考椭球的长半轴；
//dFlattening：参考椭球的扁率。

//由站心地平直角坐标转换为站心地平极坐标
void TopocentricToTopocentricPolar (PCRDTOPOCENTRICPOLAR pctp,
                                    PCRDTOPOCENTRIC pct);
//pctp：指向所转换出的站心地平极坐标的指针；
//pct：指向待转换的站心地平坐标的指针；


//由站心地平极坐标转换为站心地平直角坐标
void TopocentricPolarToTopocentric (PCRDTOPOCENTRIC pct,PCRDTOPOCENTRICPOLAR pctp);
//pct：指向所转换的站心地平坐标的指针；
//pctp：指向待转换的站心地平极坐标的指针；


#endif //LATITUDELONGITUDE_2_ENU_COORDCOVERT_H
