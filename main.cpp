//
// Created by uno on 2019/9/20.
//
#include <iostream>
#include "CoordCovert.h"
#include "math.h"
#include <fstream>
#include <cassert>
#include <string>
#include <sstream>
#include <iomanip>
//#include <backward/strstream>

void CartesianToGeodetic (PCRDGEODETIC pcg, PCRDCARTESIAN pcc,
                          double dSemiMajorAxis, double dFlattening)
{
    double e2;//第一偏心率的平方
    e2=2*dFlattening-dFlattening*dFlattening;

    pcg->longitude=atan(pcc->y/pcc->x);
    double W,N,N1=0,B,B1;
    B1=atan(pcc->z/sqrt(pcc->x*pcc->x+pcc->y*pcc->y));
    while(1)
    {
        W=sqrt(1-e2*sin(B1)*sin(B1));
        N1=dSemiMajorAxis/W;
        B=atan((pcc->z+N1*e2*sin(B1))/sqrt(pcc->x*pcc->x+pcc->y*pcc->y));

        if(fabs(B-B1)<delta)
            break;
        else
            B1=B;
    }

    pcg->latitude=B;
    N=dSemiMajorAxis/sqrt(1-e2*sin(pcg->latitude)*sin(pcg->latitude));
    pcg->height=sqrt(pcc->x*pcc->x+pcc->y*pcc->y)/cos(B)-N;


}

//由大地坐标转换为笛卡尔坐标
void GeodeticToCartesian (PCRDCARTESIAN pcc, PCRDGEODETIC pcg,
                          double dSemiMajorAxis, double dFlattening)
{
    double e2;//第一偏心率的平方
    double N;//卯酉圈半径
    e2=2*dFlattening-dFlattening*dFlattening;
    N=dSemiMajorAxis/sqrt(1-e2*sin(pcg->latitude/180*M_PI)*sin(pcg->latitude/180*M_PI));

    pcc->x = (N+pcg->height)*cos(pcg->latitude/180*M_PI)*cos(pcg->longitude/180*M_PI);
    pcc->y=(N+pcg->height)*cos(pcg->latitude/180*M_PI)*sin(pcg->longitude/180*M_PI);
    pcc->z=(N*(1-e2)+pcg->height)*sin(pcg->latitude/180*M_PI);
}

//ecef坐标系转换为enu
void topocentric(PCRDTOPOCENTRIC pct, PCRDCARTESIAN pcc, PCRDCARTESIAN pccCenter, PCRDGEODETIC pd)
{
    double dx,dy,dz;
    dx=pcc->x-pccCenter->x;
    dy=pcc->y-pccCenter->y;
    dz=pcc->z-pccCenter->z;

    //PCRDGEODETIC pd;
    //pd=(PCRDGEODETIC)malloc(sizeof(CRDGEODETIC));

    //CartesianToGeodetic (pd,pccCenter,dSemiMajorAxis,dFlattening);

    pct->northing=-sin(pd->latitude/180*M_PI)*cos(pd->longitude/180*M_PI)*dx
                  -sin(pd->latitude/180*M_PI)*sin(pd->longitude/180*M_PI)*dy
                  +cos(pd->latitude/180*M_PI)*dz;
    pct->easting=-sin(pd->longitude/180*M_PI)*dx
                 +cos(pd->longitude/180*M_PI)*dy;
    pct->upping=cos(pd->latitude/180*M_PI)*cos(pd->longitude/180*M_PI)*dx
                +cos(pd->latitude/180*M_PI)*sin(pd->longitude/180*M_PI)*dy
                +sin(pd->latitude/180*M_PI)*dz;
    //free(pd);

}

//由站心地平直角坐标转换为站心地平极坐标
void TopocentricToTopocentricPolar (PCRDTOPOCENTRICPOLAR pctp,
                                    PCRDTOPOCENTRIC pct)
{

    pctp->range=sqrt(pct->northing*pct->northing+pct->easting*pct->easting+pct->upping*pct->upping);
    pctp->azimuth=atan(pct->easting/pct->northing);
    pctp->elevation=asin(pct->upping/pctp->range);


}

//由站心地平极坐标转换为站心地平直角坐标
void TopocentricPolarToTopocentric (PCRDTOPOCENTRIC pct,
                                    PCRDTOPOCENTRICPOLAR pctp)
{
    pct->northing=pctp->range*cos(pctp->elevation)*cos(pctp->azimuth);
    pct->easting=pctp->range*cos(pctp->elevation)*sin(pctp->azimuth);
    pct->upping=pctp->range*sin(pctp->elevation);

}


int main(int argc, char *argv[]) {
    std::string fileIn = (std::string)argv[1];
    std::string fileOutXYZ = (std::string)argv[2];
    std::string fileOutENU = (std::string)argv[3];
    //std::string file_base = "/home/wayne/my_study/BIT_UAV/LatitudeLongitude_2_ENU";
    //std::string file = file_base+"/Data_ori/for_enu/data1.txt";
    std::ifstream infile;
    infile.open(fileIn.data());
    assert(infile.is_open());

    time_t t = time(0);
    char ch[64];
    strftime(ch, sizeof(ch), "%Y%m%d_%H%M%S", localtime(&t));
    //年-月-日 时-分-秒
    std::string s_ch = ch;

    //std::string file2 = file_base + "/Data_out/" + s_ch + "_xyz.txt";
    //std::string file2 = file_base + "/Data_out/data1_xyz.txt";
    //std::cout << file2 << std::endl;
    std::ofstream outfile_xyz;
    outfile_xyz.open(fileOutXYZ.data(),std::ios::out | std::ios::app);

    assert(outfile_xyz.is_open());


    //std::string file3 = file_base + "/Data_out/"+s_ch+"_enu.txt";
    //std::string file3 = file_base + "/Data_out/data1_enu.txt";
    //std::cout << file3 << std::endl;
    std::ofstream outfile_enu;
    outfile_enu.open(fileOutENU.data(),std::ios::out | std::ios::app);

    assert(outfile_enu.is_open());


    std::string s;

    CRDCARTESIAN xyz_center;
    PCRDCARTESIAN p_xyz_center;
    p_xyz_center = &xyz_center;

    while(getline(infile,s)) {

        size_t pos1 = s.find(",");
        size_t pos2 = s.find(",", pos1 + 1);
        size_t pos3 = s.find(",", pos2 + 1);
        size_t pos4 = s.find(",", pos3 + 1);
        size_t pos5 = s.find(",", pos4 + 1);
        size_t pos6 = s.find(",", pos5 + 1);
        size_t pos7 = s.find(",", pos6 + 1);

        if (s.substr(0, pos1 - 0 ) == "GPGGA") {

            std::string time_stamp = s.substr(pos2 + 1, pos3 - pos2 - 1);
            std::string B = s.substr(pos5 + 1, pos6 - pos5 - 1);
            std::string L = s.substr(pos6 + 1, pos7 - pos6 - 1);
            std::string H = s.substr(pos7 + 1);

            std::cout << "LBH:" << L << "," << B << "," << H <<std::endl;

            std::stringstream ss,ss1,ss2,ss3,ss_time_stamp;

            CRDGEODETIC lbh;
            PCRDGEODETIC p_lbh;
            p_lbh = &lbh;


            ss << L;
            ss >> p_lbh->longitude;

            ss1 << B;
            ss1 >> p_lbh->latitude;

            ss2 << H;
            ss2 >> p_lbh->height;

            CRDCARTESIAN xyz;
            PCRDCARTESIAN p_xyz;
            p_xyz = &xyz;

            CRDGEODETIC lbh_center;
            PCRDGEODETIC p_lbh_center;
            p_lbh_center = &lbh_center;

            if (lbh_center_flag == 0) {
                p_lbh_center->longitude = p_lbh->longitude;
                p_lbh_center->latitude = p_lbh->latitude;
                p_lbh_center->height = p_lbh->height;
                lbh_center_flag = 1;
            }

            GeodeticToCartesian(p_xyz, p_lbh, a, flattening);


            ss.str("");
            ss1.str("");
            ss2.str("");
            ss_time_stamp.str("");
            ss.clear();
            ss1.clear();
            ss2.clear();
            ss_time_stamp.clear();
            ss << std::fixed << p_xyz->x;
            ss1 << std::fixed << p_xyz->y;
            ss2 << std::fixed << p_xyz->z;
            std::cout << "ecef:" << time_stamp << "," << ss.str() << "," << ss1.str() << "," << ss2.str()<<std::endl;

            outfile_xyz << time_stamp << ","
                        <<ss.str() << ","
                        << ss1.str() << ","
                        << ss2.str()
                        << std::endl;



            CRDTOPOCENTRIC enu;
            PCRDTOPOCENTRIC p_enu;
            p_enu = &enu;

            if (xyz_center_flag == 0) {
                p_xyz_center->x = p_xyz->x;
                p_xyz_center->y = p_xyz->y;
                p_xyz_center->z = p_xyz->z;
                xyz_center_flag = 1;
            }


            topocentric(p_enu, p_xyz, p_xyz_center, p_lbh_center);



            ss.str("");
            ss1.str("");
            ss2.str("");
            ss.clear();
            ss1.clear();
            ss2.clear();
            ss << std::fixed << p_enu->easting;
            ss1 << std::fixed << p_enu->northing;
            ss2 << std::fixed << p_enu->upping;

            std::cout << "enu:" << time_stamp << "," << ss.str() << "," << ss1.str() << "," << ss2.str() <<std::endl<<std::endl;

            outfile_enu << time_stamp << ","
                        << ss.str() << ","
                        << ss1.str() << ","
                        << ss2.str()
                        << std::endl;

        }
    }
    infile.close();

    return 0;
}