//This code can be run using the following steps
//.L data_crossing_3D.C+
//data_crossing_3D k
//k.Loop()
//This code uses the anode-cathode-anode crossing tracks to find the Spatial distortion in 3D bins across the TPC
//https://lss.fnal.gov/archive/thesis/2000/fermilab-thesis-2021-02.pdf [For details of the method, chapter 4]

#define data_crossing_3D_cxx
#include "data_crossing_3D.h"
#include "driftvel.h"
#include "Efield_input.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <iostream>
#include <TGraphErrors.h>

int countold=0;
int countnew=0;
double x0=-358.6;//cm
double x1=358.6;//cm
double v_nominal=0.1565;//cm/us

double dotp(double A[3],double B[3],double dotAB){
  return A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
}

double crossp(double vect_A[3],double vect_B[3],double cross_P[3]){
  //std::cout<<"vectors are "<<vect_A[0]<<" "<<vect_A[1]<<" "<<vect_A[2]<<" "<<vect_B[0]<<" "<<vect_B[1]<<" "<<vect_B[2]<<std::endl;
  cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]; 
  cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]; 
  cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]; 
  //std::cout<<"cross products "<<cross_P[0]<<" "<<cross_P[1]<<"  "<<cross_P[2]<<std::endl;
  return *cross_P;
}

double shortest_int(double A1[3],double B1[3],double A2[3],double B2[3],double result[7]){//C includes point of intersections, x1,y1,z1,x2,y2,z2,dist
  double p1[3]={A1[0],A1[1],A1[2]};
  double d1[3]={-A1[0]+B1[0],-A1[1]+B1[1],-A1[2]+B1[2]};
  double p2[3]={A2[0],A2[1],A2[2]};
  double d2[3]={-A2[0]+B2[0],-A2[1]+B2[1],-A2[2]+B2[2]};
  double p2p1[3]={p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]};
  double n[3];
  crossp(d1,d2,n);
  double n2[3];
  crossp(d2,n,n2);
  double pvalue=(p2p1[0]*n2[0]+p2p1[1]*n2[1]+p2p1[2]*n2[2])/(d1[0]*n2[0]+d1[1]*n2[1]+d1[2]*n2[2]);
  double n1[3];
  crossp(d1,n,n1);
  double pvalue1=(-p2p1[0]*n1[0]-p2p1[1]*n1[1]-p2p1[2]*n1[2])/(d2[0]*n1[0]+d2[1]*n1[1]+d2[2]*n1[2]);
  double k[6]={p1[0]+pvalue*d1[0],p1[1]+pvalue*d1[1],p1[2]+pvalue*d1[2],p2[0]+pvalue1*d2[0],p2[1]+pvalue1*d2[1],p2[2]+pvalue1*d2[2]};
  double dist=sqrt(pow(k[0]-k[3],2)+pow(k[1]-k[4],2)+pow(k[2]-k[5],2));
  result[0]=k[0];
  result[1]=k[1];
  result[2]=k[2];
  result[3]=k[3];
  result[4]=k[4];
  result[5]=k[5];
  result[6]=dist;
  //result[7]={p1[0]+pvalue*d1[0],p1[1]+pvalue*d1[1],p1[2]+pvalue*d1[2],p2[0]+pvalue1*d2[0],p2[1]+pvalue1*d2[1],p2[2]+pvalue1*d2[2],dist};
  return *result;
}


void data_crossing_3D::Loop()
{
  ////just a test to see the function gives shortest distance
  double A1[3]={1,2,3};
  double B1[3]={4,5,6};
  double A2[3]={1,4,3};
  double B2[3]={3,3,-4};
  
  //std::cout<<"distance "<<smallest_dist(A1,B1,A2,B2)<<std::endl;

  double result0[7];
  shortest_int(A1,B1,A2,B2,result0);
  // std::cout<<"dist old formula "<<result0[6]<<std::endl;
  //*******************************************//
 

  double Ymin=0;
  double Ymax=600;
  
  double Zmin=0;
  double Zmax=700;

  int nbinX=18;//40cm bins along all direction
  int nbinY=30;
  int nbinZ=35;
  int nbinYy=30;
  int nbinZy=35;

  TH3F *XYZ_Zdist_bkwd=new TH3F("XYZ_Zdist_bkwd","3D Z dist_bkwd;X[cm];Y[cm];Z[cm];deltaZ[cm]",nbinX,-360,360,nbinY,Ymin,Ymax,nbinZ,Zmin,Zmax);
  TH3F *XYZ_Ydist_bkwd=new TH3F("XYZ_Ydist_bkwd","3D Y dist_bkwd;X[cm];Y[cm];Z[cm];deltaZ[cm]",nbinX,-360,360,nbinY,Ymin,Ymax,nbinZ,Zmin,Zmax);
  TH3F *XYZ_Xdist_bkwd=new TH3F("XYZ_Xdist_bkwd","3D X dist_bkwd;X[cm];Y[cm];Z[cm];deltaZ[cm]",nbinX,-360,360,nbinY,Ymin,Ymax,nbinZ,Zmin,Zmax);

  TH3F *XYZ_Zdist_fwd=new TH3F("XYZ_Zdist_fwd","3D Z dist_fwd;X[cm];Y[cm];Z[cm];deltaZ[cm]",nbinX,-360,360,nbinY,Ymin,Ymax,nbinZ,Zmin,Zmax);
  TH3F *XYZ_Ydist_fwd=new TH3F("XYZ_Ydist_fwd","3D Y dist_fwd;X[cm];Y[cm];Z[cm];deltaZ[cm]",nbinX,-360,360,nbinY,Ymin,Ymax,nbinZ,Zmin,Zmax);
  TH3F *XYZ_Xdist_fwd=new TH3F("XYZ_Xdist_fwd","3D X dist_fwd;X[cm];Y[cm];Z[cm];deltaZ[cm]",nbinX,-360,360,nbinY,Ymin,Ymax,nbinZ,Zmin,Zmax);

  TH3F *XYZ_entries_fwd=new TH3F("XYZ_entries_fwd","entries in each bin;X[cm];Y[cm];Z[cm];deltaZ[cm]",nbinX,-360,360,nbinY,Ymin,Ymax,nbinZ,Zmin,Zmax);
  TH3F *XYZ_entries_bkwd=new TH3F("XYZ_entries_bkwd","entries in each bin;X[cm];Y[cm];Z[cm];deltaZ[cm]",nbinX,-360,360,nbinY,Ymin,Ymax,nbinZ,Zmin,Zmax);

  //**************************************distortion histograms***************************************************//
  //fd or foreward distortions are distortions w.r.t to true X, Y and Z positions
  TH2F *Zdist_neg_fd=new TH2F("Zdist_neg_fd","Z distortion at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH2F *Zdist_pos_fd=new TH2F("Zdist_pos_fd","Z distortion at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH2F *Ydist_neg_fd=new TH2F("Ydist_neg_fd","Y distortion at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Ydist_pos_fd=new TH2F("Ydist_pos_fd","Y distortion at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Xdist_neg_fd=new TH2F("Xdist_neg_fd","Y distortion at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Xdist_pos_fd=new TH2F("Xdist_pos_fd","Y distortion at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);

  //bd or backward distortions are distortions w.r.t to measured X, Y and Z positions
  TH2F *Zdist_neg_bd=new TH2F("Zdist_neg_bd","Z distortion at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH2F *Zdist_pos_bd=new TH2F("Zdist_pos_bd","Z distortion at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH2F *Ydist_neg_bd=new TH2F("Ydist_neg_bd","Y distortion at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Ydist_pos_bd=new TH2F("Ydist_pos_bd","Y distortion at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Xdist_neg_bd=new TH2F("Xdist_neg_bd","Y distortion at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Xdist_pos_bd=new TH2F("Xdist_pos_bd","Y distortion at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);

  //*******************True distortion-measured distortion plots*************************************//
  //fd_delta or foreward distortions are distortions w.r.t to true X, Y and Z positions
  TH2F *Zdist_neg_fd_delta=new TH2F("Zdist_neg_fd_delta","Z distortion (input-measured) at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH2F *Zdist_pos_fd_delta=new TH2F("Zdist_pos_fd_delta","Z distortion (input-measured) at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH2F *Ydist_neg_fd_delta=new TH2F("Ydist_neg_fd_delta","Y distortion (input-measured) at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Ydist_pos_fd_delta=new TH2F("Ydist_pos_fd_delta","Y distortion (input-measured) at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Xdist_neg_fd_delta=new TH2F("Xdist_neg_fd_delta","Y distortion (input-measured) at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Xdist_pos_fd_delta=new TH2F("Xdist_pos_fd_delta","Y distortion (input-measured) at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);

  //bd_delta or backward distortion (input-measured)s are distortion (input-measured)s w.r.t to measured X, Y and Z positions
  TH2F *Zdist_neg_bd_delta=new TH2F("Zdist_neg_bd_delta","Z distortion (input-measured) at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH2F *Zdist_pos_bd_delta=new TH2F("Zdist_pos_bd_delta","Z distortion (input-measured) at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH2F *Ydist_neg_bd_delta=new TH2F("Ydist_neg_bd_delta","Y distortion (input-measured) at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Ydist_pos_bd_delta=new TH2F("Ydist_pos_bd_delta","Y distortion (input-measured) at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Xdist_neg_bd_delta=new TH2F("Xdist_neg_bd_delta","Y distortion (input-measured) at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);
  TH2F *Xdist_pos_bd_delta=new TH2F("Xdist_pos_bd_delta","Y distortion (input-measured) at cathode beam left;Z coordinate[cm];Y coordinate[cm]",nbinZy,Zmin,Zmax,nbinYy,Ymin,Ymax);

  //*************************************************************************************************************//
  TH2F *Zdist_neg_cov=new TH2F("Zdist_neg_cov","Z distortion at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH2F *Zdist_pos_cov=new TH2F("Zdist_pos_cov","Z distortion at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH2F *Zdist_neg_cov_bkwd=new TH2F("Zdist_neg_cov_bkwd","Z distortion at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH2F *Zdist_pos_cov_bkwd=new TH2F("Zdist_pos_cov_bkwd","Z distortion at cathode beam right;Z coordinate[cm];Y coordinate[cm]",nbinZ,Zmin,Zmax,nbinY,Ymin,Ymax);
  TH1F *xbin_plot=new TH1F("xbin_plot","",nbinX,-360,360);
  TH1D *dZ_pos[35];
  TH1D *dZ_neg[35];
  for(int i=0;i<35;i++){
    dZ_pos[i]=new TH1D(Form("dZ_pos_%d",i),"",200,-100,100);
    dZ_neg[i]=new TH1D(Form("dZ_neg_%d",i),"",200,-100,100);
  }

  TH1F *deltaT=new TH1F("deltaT","PeakT_max-PeakT_min;deltaT;entries",1000,0,10000);
  TH1F *deltaT_final=new TH1F("deltaT_final","PeakT_max-PeakT_min;deltaT;entries",1000,0,10000);
  TH1F *deltaT_tanode=new TH1F("deltaT_tanode","PeakT_max-PeakT_min using anode time;deltaT;entries",1000,0,10000);
  TH1F *deltaT_pos=new TH1F("deltaT_pos","PeakT_max-PeakT_min;deltaT;entries",3000,0,6000);
  TH1F *deltaT_neg=new TH1F("deltaT_neg","PeakT_max-PeakT_min;deltaT;entries",3000,0,6000);

  TH1F *dr_all=new TH1F("dr_all","3D distance;dr[cm];entries",200,0,200);
  TH1F *dr_crossing=new TH1F("dr_crossing","3D distance;dr[cm];entries",200,0,200);
  TH1F *dr_value_plot=new TH1F("dr_value_plot","3D distance;dr[cm];entries",200,0,200);
 

  TH2F *zcov_old=new TH2F("zcov_old","Z vs Y for T0 and beam triggered;Z[cm];Y[cm]",139,0,695,120,0,600);
  TH2F *zcov_new=new TH2F("zcov_new","Z vs Y T0 for new stitched;Z[cm];Y[cm]",139,0,695,120,0,600);
  TH1F *shortestdist_true=new TH1F("shortestdist_true","Shortest distance true;dist[cm];entries",200,0,100);
  TH1F *dist_measured=new TH1F("dist_measured","Measured shortest distance;dist[cm];entries",200,0,100);

  TH2F *dz_vs_x=new TH2F("dz_vs_x","dZ vs X at Y=280-320 and Z=280-320cm;X[cm];dZ",nbinX,-360,360,200,-50,50);
  TH2F *t_vs_x=new TH2F("t_vs_x","measured T vs true X;X[cm];time[ticks]",720,-360,360,360,0,5400); 

  // TFile *f=new TFile("data_anode_3D_march28_nbinX18_june15.root","RECREATE");
  // TFile *f=new TFile("data_anode_3D_oct5.root","RECREATE");
  TFile *f=new TFile("data_anode_3D_oct25_bcwd_fwd.root","RECREATE");

  //**********************************defining vectors to store values****************//
  std::vector<std::vector<std::vector<double> > > Zdistneg_fd, Zdistpos_fd,Ydistneg_fd,Ydistpos_fd,Xdistneg_fd,Xdistpos_fd,Zdistneg_bd, Zdistpos_bd,Ydistneg_bd,Ydistpos_bd,Xdistneg_bd,Xdistpos_bd;
  std::vector<std::vector<std::vector<std::vector<double> > > >  XYZdist_fwd, XYZYdist_fwd, XYZXdist_fwd, XYZdist_bkwd, XYZYdist_bkwd, XYZXdist_bkwd;

  Zdistpos_fd.resize(nbinZ);
  Zdistneg_fd.resize(nbinZ);
  for(int i=0;i<nbinZ;i++){
    Zdistpos_fd[i].resize(nbinY);
    Zdistneg_fd[i].resize(nbinY);
  }


  Ydistpos_fd.resize(nbinZy);
  Ydistneg_fd.resize(nbinZy);
  for(int i=0;i<nbinZy;i++){
    Ydistpos_fd[i].resize(nbinYy);
    Ydistneg_fd[i].resize(nbinYy);
  }

  Xdistpos_fd.resize(nbinZy);
  Xdistneg_fd.resize(nbinZy);
  for(int i=0;i<nbinZy;i++){
    Xdistpos_fd[i].resize(nbinYy);
    Xdistneg_fd[i].resize(nbinYy);
  }



  Zdistpos_bd.resize(nbinZ);
  Zdistneg_bd.resize(nbinZ);
  for(int i=0;i<nbinZ;i++){
    Zdistpos_bd[i].resize(nbinY);
    Zdistneg_bd[i].resize(nbinY);
  }

  Ydistpos_bd.resize(nbinZy);
  Ydistneg_bd.resize(nbinZy);
  for(int i=0;i<nbinZy;i++){
    Ydistpos_bd[i].resize(nbinYy);
    Ydistneg_bd[i].resize(nbinYy);
  }

  Xdistpos_bd.resize(nbinZy);
  Xdistneg_bd.resize(nbinZy);
  for(int i=0;i<nbinZy;i++){
    Xdistpos_bd[i].resize(nbinYy);
    Xdistneg_bd[i].resize(nbinYy);
  }


  XYZdist_fwd.resize(nbinX);
  XYZYdist_fwd.resize(nbinX);
  XYZXdist_fwd.resize(nbinX);
  XYZdist_bkwd.resize(nbinX);
  XYZYdist_bkwd.resize(nbinX);
  XYZXdist_bkwd.resize(nbinX);


  for(int i=0;i<nbinX;i++){
    XYZdist_fwd[i].resize(nbinY);
    XYZYdist_fwd[i].resize(nbinY);
    XYZXdist_fwd[i].resize(nbinY);
    XYZdist_bkwd[i].resize(nbinY);
    XYZYdist_bkwd[i].resize(nbinY);
    XYZXdist_bkwd[i].resize(nbinY);


  }
  for(int i=0;i<nbinX;i++){
    for(int j=0;j<nbinY;j++){
      XYZdist_fwd[i][j].resize(nbinZ);
      XYZYdist_fwd[i][j].resize(nbinZ);
      XYZXdist_fwd[i][j].resize(nbinZ);
      XYZdist_bkwd[i][j].resize(nbinZ);
      XYZYdist_bkwd[i][j].resize(nbinZ);
      XYZXdist_bkwd[i][j].resize(nbinZ);
    }
  }


  
  ///*****************************************************************///
  std::vector<std::vector<double> > t_neg_tracks,z_neg_tracks,y_neg_tracks,tpc_neg_tracks,t_pos_tracks,z_pos_tracks,y_pos_tracks,tpc_pos_tracks,t_all_tracks,z_all_tracks,y_all_tracks,tpc_all_tracks;

 

  if (fChain == 0) return;
  std::vector<double> t_neg,t_pos,z_pos,z_neg,y_pos,y_neg, t_all,z_all,y_all,tpc_pos,tpc_neg,tpc_all;
  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  //  for(Long64_t jentry=0; jentry<100000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if(jentry%1000==0) std::cout<<nentries<<"/"<<jentry<<std::endl;
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

      
     
    if(!trkhitz_wire2->size()) continue;
    std::vector<std::vector<double> > t_pos_buff,z_pos_buff,y_pos_buff,tpc_pos_buff,t_neg_buff,z_neg_buff,y_neg_buff,tpc_neg_buff,t_all_buff,z_all_buff,y_all_buff,tpc_all_buff;
    t_neg_buff.clear();t_pos_buff.clear();z_pos_buff.clear();z_neg_buff.clear();y_pos_buff.clear();y_neg_buff.clear();z_all_buff.clear();t_all_buff.clear();y_all_buff.clear();tpc_pos_buff.clear();tpc_neg_buff.clear();tpc_all_buff.clear();
    for(size_t i=0;i<trkhitz_wire2->size();i++){ //This is the loop for each track, it's going to be a big loop
      if(trkhitz_wire2->at(i).size()<5) continue;
      t_neg.clear();t_pos.clear();z_pos.clear();z_neg.clear();y_pos.clear();y_neg.clear();z_all.clear();t_all.clear();y_all.clear();tpc_pos.clear();tpc_neg.clear();tpc_all.clear();
      for(size_t j=0;j<trkhitz_wire2->at(i).size();j++){
   	if(hit_tpc2->at(i)[j]==2||hit_tpc2->at(i)[j]==6||hit_tpc2->at(i)[j]==10){
   	  t_pos.push_back(hit_peakT2->at(i)[j]);
   	  z_pos.push_back(trkhitz_wire2->at(i)[j]);
   	  y_pos.push_back(trkhity2->at(i)[j]);
   	  tpc_pos.push_back(hit_tpc2->at(i)[j]);
   	}

   	if(hit_tpc2->at(i)[j]==1||hit_tpc2->at(i)[j]==5||hit_tpc2->at(i)[j]==9){
   	  t_neg.push_back(hit_peakT2->at(i)[j]);
   	  z_neg.push_back(trkhitz_wire2->at(i)[j]);
   	  y_neg.push_back(trkhity2->at(i)[j]);
   	  tpc_neg.push_back(hit_tpc2->at(i)[j]);
   	}
      }//j loop
      if(t_pos.size()){
   	deltaT_pos->Fill(*std::max_element(t_pos.begin(),t_pos.end())-*std::min_element(t_pos.begin(),t_pos.end()));
      }
      if(t_neg.size()){
   	deltaT_neg->Fill(*std::max_element(t_neg.begin(),t_neg.end())-*std::min_element(t_neg.begin(),t_neg.end()));
      }
      if(t_pos.size()==0){
   	t_neg_buff.push_back(t_neg);
   	z_neg_buff.push_back(z_neg);
   	y_neg_buff.push_back(y_neg);
   	tpc_neg_buff.push_back(tpc_neg);
      }
      if(t_neg.size()==0){
   	t_pos_buff.push_back(t_pos);
   	z_pos_buff.push_back(z_pos);
   	y_pos_buff.push_back(y_pos);
   	tpc_pos_buff.push_back(tpc_pos);
      }
      if(t_pos.size()==0||t_neg.size()==0) continue;
      if(t_neg[0]>t_neg[t_neg.size()-1]){
   	std::reverse(t_neg.begin(),t_neg.end());
   	std::reverse(z_neg.begin(),z_neg.end());
   	std::reverse(y_neg.begin(),y_neg.end());
   	std::reverse(tpc_neg.begin(),tpc_neg.end());
      }
      if(t_pos[0]<t_pos[t_pos.size()-1]){
   	std::reverse(t_pos.begin(),t_pos.end());
   	std::reverse(z_pos.begin(),z_pos.end());
   	std::reverse(y_pos.begin(),y_pos.end());
   	std::reverse(tpc_pos.begin(),tpc_pos.end());
      }

      double t_max_pos=*std::max_element(t_pos.begin(),t_pos.end());
      double t_min_pos=*std::min_element(t_pos.begin(),t_pos.end());
      double t_max_neg=*std::max_element(t_neg.begin(),t_neg.end());
      double t_min_neg=*std::min_element(t_neg.begin(),t_neg.end());
      double dT=t_max_pos-t_min_pos+t_max_neg-t_min_neg;
      double dr_value=sqrt(pow(z_neg[z_neg.size()-1]-z_pos[0],2)+pow(y_neg[y_neg.size()-1]-y_pos[0],2)+pow(2*x1-(t_neg[t_neg.size()-1]+t_pos[0]-t_neg[0]-t_pos[t_pos.size()-1])*0.5*v_nominal,2));
      dr_all->Fill(dr_value);
      if(!(abs(t_min_pos-t_min_neg)<10 && dT>9000 && dT<9220 && t_max_pos<5985 && t_max_neg<5985)) continue;
      //if(!(dT>9120 && dT<9160)) continue;
      //std::cout<<"event min max later dT "<<event<<"  "<<t_min_pos<<"  "<<t_min_neg<<"  "<<t_max_pos<<"  "<<t_max_neg<<"  "<<dT<<" dr "<<dr_value<<std::endl;
      dr_crossing->Fill(dr_value);
      deltaT->Fill(dT);
      deltaT_final->Fill(dT);
      for(size_t t=0;t<t_neg.size();t++){
   	t_all.push_back(t_neg[t]);
   	z_all.push_back(z_neg[t]);
   	y_all.push_back(y_neg[t]);
   	tpc_all.push_back(tpc_neg[t]);
      }
      for(size_t t=0;t<t_pos.size();t++){
   	t_all.push_back(t_pos[t]);
   	z_all.push_back(z_pos[t]);
   	y_all.push_back(y_pos[t]);
   	tpc_all.push_back(tpc_pos[t]);
      }
      t_neg_tracks.push_back(t_neg);
      z_neg_tracks.push_back(z_neg);
      y_neg_tracks.push_back(y_neg);
      tpc_neg_tracks.push_back(tpc_neg);
       
      t_pos_tracks.push_back(t_pos);
      z_pos_tracks.push_back(z_pos);
      y_pos_tracks.push_back(y_pos);
      tpc_pos_tracks.push_back(tpc_pos);

      t_all_tracks.push_back(t_all);
      z_all_tracks.push_back(z_all);
      y_all_tracks.push_back(y_all);
      tpc_all_tracks.push_back(tpc_all);
      countold++;
      zcov_old->Fill(z_pos[0],y_pos[0]);

    }//i loop

    if(!(t_neg_buff.size() && t_pos_buff.size())) continue;
    for(size_t t=0;t<t_neg_buff.size();t++){
      for(size_t t1=0;t1<t_pos_buff.size();t1++){
   	double t_max_pos=*std::max_element(t_pos_buff[t1].begin(),t_pos_buff[t1].end());
   	double t_min_pos=*std::min_element(t_pos_buff[t1].begin(),t_pos_buff[t1].end());
   	double t_max_neg=*std::max_element(t_neg_buff[t].begin(),t_neg_buff[t].end());
   	double t_min_neg=*std::min_element(t_neg_buff[t].begin(),t_neg_buff[t].end());
   	double dT=t_max_pos-t_min_pos+t_max_neg-t_min_neg;
   	if(!(abs(t_min_pos-t_min_neg)<10 && dT>9000 && dT<9220 && t_max_pos<5985 && t_max_neg<5985)) continue;
   	//if(!(dT>9120 && dT<9160)) continue;
   	//std::cout<<"event min max later dT "<<event<<"  "<<t_min_pos<<"  "<<t_min_neg<<"  "<<t_max_pos<<"  "<<t_max_neg<<"  "<<dT<<std::endl;
   	if(t_neg_buff[t][0]>t_neg_buff[t][t_neg_buff[t].size()-1]){
   	  std::reverse(t_neg_buff[t].begin(),t_neg_buff[t].end());
   	  std::reverse(z_neg_buff[t].begin(),z_neg_buff[t].end());
   	  std::reverse(y_neg_buff[t].begin(),y_neg_buff[t].end());
   	  std::reverse(tpc_neg_buff[t].begin(),tpc_neg_buff[t].end());
   	}
   	if(t_pos_buff[t1][0]<t_pos_buff[t1][t_pos_buff[t1].size()-1]){
   	  std::reverse(t_pos_buff[t1].begin(),t_pos_buff[t1].end());
   	  std::reverse(z_pos_buff[t1].begin(),z_pos_buff[t1].end());
   	  std::reverse(y_pos_buff[t1].begin(),y_pos_buff[t1].end());
   	  std::reverse(tpc_pos_buff[t1].begin(),tpc_pos_buff[t1].end());
   	}
   	double dr_value=sqrt(pow(z_neg_buff[t][z_neg_buff[t].size()-1]-z_pos_buff[t1][0],2)+pow(y_neg_buff[t][y_neg_buff[t].size()-1]-y_pos_buff[t1][0],2)+pow(2*x1-(t_neg_buff[t][t_neg_buff[t].size()-1]+t_pos_buff[t1][0]-t_neg_buff[t][0]-t_pos_buff[t1][t_pos_buff[t1].size()-1])*0.5*0.1565,2));
   	dr_value_plot->Fill(dr_value);
   	if(dr_value>30.0) continue;
   	deltaT_final->Fill(dT);
   	std::vector<double> t_all_b,z_all_b,y_all_b,tpc_all_b;
   	t_all_b.clear(); z_all_b.clear(); y_all_b.clear(); tpc_all_b.clear();
   	for(size_t k1=0;k1<t_neg_buff[t].size();k1++){
   	  t_all_b.push_back(t_neg_buff[t][k1]);
   	  z_all_b.push_back(z_neg_buff[t][k1]);
   	  y_all_b.push_back(y_neg_buff[t][k1]);
   	  tpc_all_b.push_back(tpc_neg_buff[t][k1]);
   	}

   	for(size_t k2=0;k2<t_pos_buff[t1].size();k2++){
   	  t_all_b.push_back(t_pos_buff[t1][k2]);
   	  z_all_b.push_back(z_pos_buff[t1][k2]);
   	  y_all_b.push_back(y_pos_buff[t1][k2]);
   	  tpc_all_b.push_back(tpc_pos_buff[t1][k2]);
   	}

   	t_neg_tracks.push_back(t_neg_buff[t]);
   	z_neg_tracks.push_back(z_neg_buff[t]);
   	y_neg_tracks.push_back(y_neg_buff[t]);
   	tpc_neg_tracks.push_back(tpc_neg_buff[t]);
       
   	t_pos_tracks.push_back(t_pos_buff[t1]);
   	z_pos_tracks.push_back(z_pos_buff[t1]);
   	y_pos_tracks.push_back(y_pos_buff[t1]);
   	tpc_pos_tracks.push_back(tpc_pos_buff[t1]);

   	t_all_tracks.push_back(t_all_b);
   	z_all_tracks.push_back(z_all_b);
   	y_all_tracks.push_back(y_all_b);
   	tpc_all_tracks.push_back(tpc_all_b);
   	zcov_new->Fill(z_pos_buff[t1][0],y_pos_buff[t1][0]);
   	countnew++;
      }
    }
      
  }//jentry loop
  ///////////measuring the distortion values here////////////////////
  for(size_t t=0;t<t_all_tracks.size();t++){
    //calculate the Z distortion here
    double z0=z_all_tracks[t][0];
    double z1=z_all_tracks[t][z_all_tracks[t].size()-1];
    double y0=y_all_tracks[t][0];
    double y1=y_all_tracks[t][y_all_tracks[t].size()-1];

    double x_cathode_neg=-0.15;//cm
    double x_cathode_pos=0.15;//cm
    double z_true_pos=z0+(z1-z0)*(x_cathode_pos-x0)/(x1-x0);
    double z_true_neg=z0+(z1-z0)*(x_cathode_neg-x0)/(x1-x0);

    double y_true_pos=y0+(y1-y0)*(x_cathode_pos-x0)/(x1-x0);
    double y_true_neg=y0+(y1-y0)*(x_cathode_neg-x0)/(x1-x0);

    ////***********forward binning************////
    int y_bin_pos=y_true_pos/20;
    int z_bin_pos=z_true_pos/20;
    int y_bin_neg=y_true_neg/20;
    int z_bin_neg=z_true_neg/20;

    ////************backward binning***************////
    int y_bin_pos_bd=y_pos_tracks[t][0]/20;
    int z_bin_pos_bd=z_pos_tracks[t][0]/20;
    int y_bin_neg_bd=y_neg_tracks[t][y_neg_tracks[t].size()-1]/20;
    int z_bin_neg_bd=z_neg_tracks[t][z_neg_tracks[t].size()-1]/20;
        
    if(!((y_bin_pos<0||y_bin_pos>=30)||(y_bin_neg<0||y_bin_neg>=30)||(z_bin_pos<0||z_bin_pos>=35)||(z_bin_neg<0||z_bin_neg>=35))){
      Ydistpos_fd[z_bin_pos][y_bin_pos].push_back(y_pos_tracks[t][0]-y_true_pos);
      Ydistneg_fd[z_bin_neg][y_bin_neg].push_back(y_neg_tracks[t][y_neg_tracks[t].size()-1]-y_true_neg);
      Zdist_pos_cov->Fill(z_true_pos,y_true_pos);
      Zdist_neg_cov->Fill(z_true_neg,y_true_neg);
      Zdist_pos_cov_bkwd->Fill(z_pos_tracks[t][0],y_pos_tracks[t][0]);
      Zdist_neg_cov_bkwd->Fill(z_neg_tracks[t][z_neg_tracks[t].size()-1],z_neg_tracks[t][z_neg_tracks[t].size()-1]);
      Zdistpos_fd[z_bin_pos][y_bin_pos].push_back(z_pos_tracks[t][0]-z_true_pos);
      Zdistneg_fd[z_bin_neg][y_bin_neg].push_back(z_neg_tracks[t][z_neg_tracks[t].size()-1]-z_true_neg);
      Xdistpos_fd[z_bin_pos][y_bin_pos].push_back(x1-(t_pos_tracks[t][0]-t_pos_tracks[t][t_pos_tracks[t].size()-1])*0.5*v_nominal-x_cathode_pos);
      Xdistneg_fd[z_bin_neg][y_bin_neg].push_back(x1-(t_neg_tracks[t][t_neg_tracks[t].size()-1]-t_neg_tracks[t][0])*0.5*v_nominal-x_cathode_neg);
    }
    if(!((y_bin_pos_bd<0||y_bin_pos_bd>=30)||(y_bin_neg_bd<0||y_bin_neg_bd>=30)||(z_bin_pos_bd<0||z_bin_pos_bd>=35)||(z_bin_neg_bd<0||z_bin_neg_bd>=35))){
      Ydistpos_bd[z_bin_pos_bd][y_bin_pos_bd].push_back(-y_pos_tracks[t][0]+y_true_pos);
      Ydistneg_bd[z_bin_neg_bd][y_bin_neg_bd].push_back(-y_neg_tracks[t][y_neg_tracks[t].size()-1]+y_true_neg);
      Zdistpos_bd[z_bin_pos_bd][y_bin_pos_bd].push_back(-z_pos_tracks[t][0]+z_true_pos);
      Zdistneg_bd[z_bin_neg_bd][y_bin_neg_bd].push_back(-z_neg_tracks[t][z_neg_tracks[t].size()-1]+z_true_neg);
      Xdistpos_bd[z_bin_pos_bd][y_bin_pos_bd].push_back(-(x1-(t_pos_tracks[t][0]-t_pos_tracks[t][t_pos_tracks[t].size()-1])*0.5*v_nominal-x_cathode_pos));
      Xdistneg_bd[z_bin_neg_bd][y_bin_neg_bd].push_back(-(x1-(t_neg_tracks[t][t_neg_tracks[t].size()-1]-t_neg_tracks[t][0])*0.5*v_nominal-x_cathode_neg));
    }
  }
  ///***************filling the XYZ histograms here**************************/


  for(int i=0;i<nbinZ;i++){
    for(int j=0;j<nbinY;j++){
      //Z pos fd


      if(Zdistpos_fd[i][j].size()){
	Zdist_pos_fd->SetBinContent(i+1,j+1,TMath::Median(Zdistpos_fd[i][j].size(),&Zdistpos_fd[i][j][0]));
	Zdist_pos_fd_delta->SetBinContent(i+1,j+1,zoffsetfd(0.15,j*20.0+10,i*20+10)-TMath::Median(Zdistpos_fd[i][j].size(),&Zdistpos_fd[i][j][0]));
      }
      if(Zdistpos_fd[i][j].size()==0){
	Zdist_pos_fd->SetBinContent(i+1,j+1,-999999);
	Zdist_pos_fd_delta->SetBinContent(i+1,j+1,-999999);
      }
      ///Z dist neg fd
      if(Zdistneg_fd[i][j].size()){
	Zdist_neg_fd->SetBinContent(i+1,j+1,TMath::Median(Zdistneg_fd[i][j].size(),&Zdistneg_fd[i][j][0]));
	Zdist_neg_fd_delta->SetBinContent(i+1,j+1,zoffsetfd(-0.15,j*20.0+10,i*20+10)-TMath::Median(Zdistneg_fd[i][j].size(),&Zdistneg_fd[i][j][0]));
      }
      if(Zdistneg_fd[i][j].size()==0){
	Zdist_neg_fd->SetBinContent(i+1,j+1,-999999);
	Zdist_neg_fd_delta->SetBinContent(i+1,j+1,-999999);
      }


      ///////////////////////////////////
      if(Ydistpos_fd[i][j].size()){
	Ydist_pos_fd->SetBinContent(i+1,j+1,TMath::Median(Ydistpos_fd[i][j].size(),&Ydistpos_fd[i][j][0]));
	Ydist_pos_fd_delta->SetBinContent(i+1,j+1,yoffsetfd(0.15,j*20.0+10,i*20+10)-TMath::Median(Ydistpos_fd[i][j].size(),&Ydistpos_fd[i][j][0]));
      }
      if(Ydistpos_fd[i][j].size()==0){
	Ydist_pos_fd->SetBinContent(i+1,j+1,-999999);
	Ydist_pos_fd_delta->SetBinContent(i+1,j+1,-999999);
      }
      ///////////////////Y neg fd///////////////////////
      if(Ydistneg_fd[i][j].size()){
	Ydist_neg_fd->SetBinContent(i+1,j+1,TMath::Median(Ydistneg_fd[i][j].size(),&Ydistneg_fd[i][j][0]));
	Ydist_neg_fd_delta->SetBinContent(i+1,j+1,yoffsetfd(-0.15,j*20.0+10,i*20+10)-TMath::Median(Ydistneg_fd[i][j].size(),&Ydistneg_fd[i][j][0]));
      }
      if(Ydistneg_fd[i][j].size()==0){
	Ydist_neg_fd->SetBinContent(i+1,j+1,-999999);
	Ydist_neg_fd_delta->SetBinContent(i+1,j+1,-999999);
      }


      //////////////////////////////////////
      if(Xdistpos_fd[i][j].size()){
	Xdist_pos_fd->SetBinContent(i+1,j+1,TMath::Median(Xdistpos_fd[i][j].size(),&Xdistpos_fd[i][j][0]));
	Xdist_pos_fd_delta->SetBinContent(i+1,j+1,xoffsetfd(0.15,j*20.0+10,i*20+10)-TMath::Median(Xdistpos_fd[i][j].size(),&Xdistpos_fd[i][j][0]));
      }
      if(Xdistpos_fd[i][j].size()==0){
	Xdist_pos_fd->SetBinContent(i+1,j+1,-999999);
	Xdist_pos_fd_delta->SetBinContent(i+1,j+1,-999999);
      }
      ////X neg fd
      if(Xdistneg_fd[i][j].size()){
	Xdist_neg_fd->SetBinContent(i+1,j+1,TMath::Median(Xdistneg_fd[i][j].size(),&Xdistneg_fd[i][j][0]));
	Xdist_neg_fd_delta->SetBinContent(i+1,j+1,xoffsetfd(-0.15,j*20.0+10,i*20+10)-TMath::Median(Xdistneg_fd[i][j].size(),&Xdistneg_fd[i][j][0]));
      }
      if(Xdistneg_fd[i][j].size()==0){
	Xdist_neg_fd->SetBinContent(i+1,j+1,-999999);
	Xdist_neg_fd_delta->SetBinContent(i+1,j+1,-999999);
      }



      ///////backward displacement values here////////

      //Z pos bd
      if(Zdistpos_bd[i][j].size()){
	Zdist_pos_bd->SetBinContent(i+1,j+1,TMath::Median(Zdistpos_bd[i][j].size(),&Zdistpos_bd[i][j][0]));
	Zdist_pos_bd_delta->SetBinContent(i+1,j+1,zoffsetbd(0.15,j*20.0+10,i*20+10)-TMath::Median(Zdistpos_bd[i][j].size(),&Zdistpos_bd[i][j][0]));
	 

      }
      if(Zdistpos_bd[i][j].size()==0){
	Zdist_pos_bd->SetBinContent(i+1,j+1,-999999);
	Zdist_pos_bd_delta->SetBinContent(i+1,j+1,-999999);
      }
      ///Z dist neg bd
      if(Zdistneg_bd[i][j].size()){
	Zdist_neg_bd->SetBinContent(i+1,j+1,TMath::Median(Zdistneg_bd[i][j].size(),&Zdistneg_bd[i][j][0]));
	Zdist_neg_bd_delta->SetBinContent(i+1,j+1,zoffsetbd(-0.15,j*20.0+10,i*20+10)-TMath::Median(Zdistneg_bd[i][j].size(),&Zdistneg_bd[i][j][0]));
      }
      if(Zdistneg_bd[i][j].size()==0){
	Zdist_neg_bd->SetBinContent(i+1,j+1,-999999);
	Zdist_neg_bd_delta->SetBinContent(i+1,j+1,-999999);
      }


      ///////////////////////////////////
      if(Ydistpos_bd[i][j].size()){
	Ydist_pos_bd->SetBinContent(i+1,j+1,TMath::Median(Ydistpos_bd[i][j].size(),&Ydistpos_bd[i][j][0]));
	Ydist_pos_bd_delta->SetBinContent(i+1,j+1,yoffsetbd(0.15,j*20.0+10,i*20+10)-TMath::Median(Ydistpos_bd[i][j].size(),&Ydistpos_bd[i][j][0]));
      }
      if(Ydistpos_bd[i][j].size()==0){
	Ydist_pos_bd->SetBinContent(i+1,j+1,-999999);
	Ydist_pos_bd_delta->SetBinContent(i+1,j+1,-999999);
      }
      ///////////////////Y neg bd///////////////////////
      if(Ydistneg_bd[i][j].size()){
	Ydist_neg_bd->SetBinContent(i+1,j+1,TMath::Median(Ydistneg_bd[i][j].size(),&Ydistneg_bd[i][j][0]));
	Ydist_neg_bd_delta->SetBinContent(i+1,j+1,yoffsetbd(-0.15,j*20.0+10,i*20+10)-TMath::Median(Ydistneg_bd[i][j].size(),&Ydistneg_bd[i][j][0]));
      }
      if(Ydistneg_bd[i][j].size()==0){
	Ydist_neg_bd->SetBinContent(i+1,j+1,-999999);
	Ydist_neg_bd_delta->SetBinContent(i+1,j+1,-999999);
      }


      //////////////////////////////////////
      if(Xdistpos_bd[i][j].size()){
	Xdist_pos_bd->SetBinContent(i+1,j+1,TMath::Median(Xdistpos_bd[i][j].size(),&Xdistpos_bd[i][j][0]));
	Xdist_pos_bd_delta->SetBinContent(i+1,j+1,xoffsetbd(0.15,j*20.0+10,i*20+10)-TMath::Median(Xdistpos_bd[i][j].size(),&Xdistpos_bd[i][j][0]));
      }
      if(Xdistpos_bd[i][j].size()==0){
	Xdist_pos_bd->SetBinContent(i+1,j+1,-999999);
	Xdist_pos_bd_delta->SetBinContent(i+1,j+1,-999999);
      }
      ////X neg bd
      if(Xdistneg_bd[i][j].size()){
	Xdist_neg_bd->SetBinContent(i+1,j+1,TMath::Median(Xdistneg_bd[i][j].size(),&Xdistneg_bd[i][j][0]));
	Xdist_neg_bd_delta->SetBinContent(i+1,j+1,xoffsetbd(-0.15,j*20.0+10,i*20+10)-TMath::Median(Xdistneg_bd[i][j].size(),&Xdistneg_bd[i][j][0]));
      }
      if(Xdistneg_bd[i][j].size()==0){
	Xdist_neg_bd->SetBinContent(i+1,j+1,-999999);
	Xdist_neg_bd_delta->SetBinContent(i+1,j+1,-999999);
      }
      ///////////////////////////////
    }
  }//nbinZ



 

  ///////////////////////////////**********************************//////
  int isize=t_all_tracks.size();
  for(size_t t=0;t<t_all_tracks.size();t++){
    if(t%100==0) std::cout<<t<<"/"<<isize<<std::endl;
    for(size_t t1=0;t1<t_all_tracks.size();t1++){
      if(t==t1) continue;
      double xyzdist[7];
      double A1[3]={x0,y_all_tracks[t][0],z_all_tracks[t][0]};
      double B1[3]={x1,y_all_tracks[t][y_all_tracks[t].size()-1],z_all_tracks[t][z_all_tracks[t].size()-1]};
      double A2[3]={x0,y_all_tracks[t1][0],z_all_tracks[t1][0]};
      double B2[3]={x1,y_all_tracks[t1][y_all_tracks[t1].size()-1],z_all_tracks[t1][z_all_tracks[t1].size()-1]};
      shortest_int(A1,B1,A2,B2,xyzdist);
      shortestdist_true->Fill(xyzdist[6]);
      if(xyzdist[6]>1.0 || abs(xyzdist[0])<-360 ||abs(xyzdist[3])>360) continue; 
      double xyzdist_measured[7]={99999,99999,99999,99999,99999,99999,99999};
      double dist_0=99999;
      if(xyzdist[0]<0){
   	for(int ik=0;ik<t_neg_tracks[t].size();ik++){
   	  double xt=x0+(t_neg_tracks[t][ik]-t_neg_tracks[t][0])*0.5*v_nominal;
   	  if(abs(xt-xyzdist[0])>100) continue;
   	  double yt=y_neg_tracks[t][ik];
   	  double zt=z_neg_tracks[t][ik];
   	  for(int im=0;im<t_neg_tracks[t1].size();im++){
   	    double xt1=x0+(t_neg_tracks[t1][im]-t_neg_tracks[t1][0])*0.5*v_nominal;
   	    if(abs(xt1-xyzdist[0])>100 || abs(xt1-xt)>5.0) continue;
   	    double yt1=y_neg_tracks[t1][im];
   	    double zt1=z_neg_tracks[t1][im];
   	    double dist_1=sqrt(pow(xt-xt1,2)+pow(yt-yt1,2)+pow(zt-zt1,2));
   	    if(dist_1<dist_0){
   	      xyzdist_measured[0]=xt; xyzdist_measured[1]=yt; xyzdist_measured[2]=zt;
   	      xyzdist_measured[3]=xt1; xyzdist_measured[4]=yt1; xyzdist_measured[5]=zt1;
   	      xyzdist_measured[6]=dist_1;
   	      dist_0=dist_1;
   	    }//dist_1
   	  }//im
   	}//ik
      }//if
       
      if(xyzdist[0]>0){
   	for(int ik=0;ik<t_pos_tracks[t].size();ik++){
   	  double xt=x1-(t_pos_tracks[t][ik]-t_pos_tracks[t][t_pos_tracks[t].size()-1])*0.5*v_nominal;
   	  if(abs(xt-xyzdist[0])>100) continue;
   	  double yt=y_pos_tracks[t][ik];
   	  double zt=z_pos_tracks[t][ik];
   	  for(int im=0;im<t_pos_tracks[t1].size();im++){
   	    double xt1=x1-(t_pos_tracks[t1][im]-t_pos_tracks[t1][t_pos_tracks[t1].size()-1])*0.5*v_nominal;
   	    if(abs(xt1-xyzdist[0])>100 || abs(xt1-xt)>5.0) continue;
   	    double yt1=y_pos_tracks[t1][im];
   	    double zt1=z_pos_tracks[t1][im];
   	    double dist_1=sqrt(pow(xt-xt1,2)+pow(yt-yt1,2)+pow(zt-zt1,2));
   	    if(dist_1<dist_0){
   	      xyzdist_measured[0]=xt; xyzdist_measured[1]=yt; xyzdist_measured[2]=zt;
   	      xyzdist_measured[3]=xt1; xyzdist_measured[4]=yt1; xyzdist_measured[5]=zt1;
   	      xyzdist_measured[6]=dist_1;
   	      dist_0=dist_1;
   	    }//dist_1
   	  }//im
   	}//ik
      }//if


       //  dist_measured->Fill(xyzdist_measured[6]);
      if(xyzdist_measured[6]<2.0){
	//Backward displacement
	int xbin=xbin_plot->FindBin((xyzdist_measured[0]+xyzdist_measured[3])/2.0);
	int ybin=(xyzdist_measured[1]+xyzdist_measured[4])/(2*20);
	int zbin=(xyzdist_measured[2]+xyzdist_measured[5])/(2*20);
	if(xbin>0 && xbin<=18 && ybin>=0 && ybin<30 && zbin>=0 && zbin<35) XYZdist_bkwd[xbin-1][ybin][zbin].push_back(-(xyzdist_measured[2]+xyzdist_measured[5])/2.0+(xyzdist[2]+xyzdist[5])/2.0);
	if(xbin>0 && xbin<=18 && ybin>=0 && ybin<30 && zbin>=0 && zbin<35) XYZYdist_bkwd[xbin-1][ybin][zbin].push_back(-(xyzdist_measured[1]+xyzdist_measured[4])/2.0+(xyzdist[1]+xyzdist[4])/2.0);
	if(xbin>0 && xbin<=18 && ybin>=0 && ybin<30 && zbin>=0 && zbin<35) XYZXdist_bkwd[xbin-1][ybin][zbin].push_back(-(xyzdist_measured[0]+xyzdist_measured[3])/2.0+(xyzdist[0]+xyzdist[3])/2.0);
	//Forward displacement
	int xbinf=xbin_plot->FindBin((xyzdist[0]+xyzdist[3])/2.0);
	int ybinf=(xyzdist[1]+xyzdist[4])/(2*20);
	int zbinf=(xyzdist[2]+xyzdist[5])/(2*20);
	if(xbinf>0 && xbinf<=18 && ybinf>=0 && ybinf<30 && zbinf>=0 && zbinf<35) XYZdist_fwd[xbinf-1][ybinf][zbinf].push_back((xyzdist_measured[2]+xyzdist_measured[5])/2.0-(xyzdist[2]+xyzdist[5])/2.0);
	if(xbinf>0 && xbinf<=18 && ybinf>=0 && ybinf<30 && zbinf>=0 && zbinf<35) XYZYdist_fwd[xbinf-1][ybinf][zbinf].push_back((xyzdist_measured[1]+xyzdist_measured[4])/2.0-(xyzdist[1]+xyzdist[4])/2.0);
	if(xbinf>0 && xbinf<=18 && ybinf>=0 && ybinf<30 && zbinf>=0 && zbinf<35) XYZXdist_fwd[xbinf-1][ybinf][zbinf].push_back((xyzdist_measured[0]+xyzdist_measured[3])/2.0-(xyzdist[0]+xyzdist[3])/2.0);

      }
    }
  }

 

  /////////////////bkward displacement filling/////////////////////
  for(int i=0;i<nbinX;i++){
    for(int j=0;j<nbinY;j++){
      for(int k=0;k<nbinZ;k++){
	if(XYZdist_bkwd[i][j][k].size()){
	  XYZ_Zdist_bkwd->SetBinContent(i+1,j+1,k+1,TMath::Median(XYZdist_bkwd[i][j][k].size(),&XYZdist_bkwd[i][j][k][0]));
	  //  std::cout<<"i, j, k and content "<<i<<" "<<j<<" k "<<k<<TMath::Median(XYZdist_bkwd[i][j][k].size(),&XYZdist_bkwd[i][j][k][0])<<std::endl;
	  if(j==20 && i==8){
	    float zvarb=20*k+10;
	    for(int l=0;l<XYZdist_bkwd[i][j][k].size();l++){
	      dZ_neg[k]->Fill(XYZdist_bkwd[i][j][k][l]);
	    }
	    dZ_neg[k]->Write(Form("bincenter_X%d_Y%d_Z%0.1f",-20,410,zvarb));
	  }
	  if(j==20 && i==9){
	    float zvarb=20*k+10;
	    for(int l=0;l<XYZdist_bkwd[i][j][k].size();l++){
	      
	      dZ_pos[k]->Fill(XYZdist_bkwd[i][j][k][l]);
	    }
	    dZ_pos[k]->Write(Form("bincenter_X%d_Y%d_Z%0.1f",20,410,zvarb));
	  }
	     
	}
	if(XYZdist_bkwd[i][j][k].size()==0) XYZ_Zdist_bkwd->SetBinContent(i+1,j+1,k+1,-99999999);
      }
    }
  }

  for(int i=0;i<nbinX;i++){
    for(int j=0;j<nbinY;j++){
      for(int k=0;k<nbinZ;k++){
	if(XYZYdist_bkwd[i][j][k].size()){
	  XYZ_Ydist_bkwd->SetBinContent(i+1,j+1,k+1,TMath::Median(XYZYdist_bkwd[i][j][k].size(),&XYZYdist_bkwd[i][j][k][0]));
	  //  std::cout<<"i, j, k and content "<<i<<" "<<j<<" k "<<k<<TMath::Median(XYZdist_bkwd[i][j][k].size(),&XYZdist_bkwd[i][j][k][0])<<std::endl;
	}
	if(XYZYdist_bkwd[i][j][k].size()==0) XYZ_Ydist_bkwd->SetBinContent(i+1,j+1,k+1,-99999999);
      }
    }
  } 
  for(int i=0;i<nbinX;i++){
    for(int j=0;j<nbinY;j++){
      for(int k=0;k<nbinZ;k++){
	double sizebkwd=XYZXdist_bkwd[i][j][k].size();
	XYZ_entries_bkwd->SetBinContent(i+1,j+1,k+1,sizebkwd);
	if(XYZXdist_bkwd[i][j][k].size()){
	  XYZ_Xdist_bkwd->SetBinContent(i+1,j+1,k+1,TMath::Median(XYZXdist_bkwd[i][j][k].size(),&XYZXdist_bkwd[i][j][k][0]));
	  //  std::cout<<"i, j, k and content "<<i<<" "<<j<<" k "<<k<<TMath::Median(XYZdist_bkwd[i][j][k].size(),&XYZdist_bkwd[i][j][k][0])<<std::endl;
	}
	if(XYZXdist_bkwd[i][j][k].size()==0) XYZ_Xdist_bkwd->SetBinContent(i+1,j+1,k+1,-99999999);
      }
    }
  } 

  //////////////////foward displacement filling/////////////////

  /////////////////bkward displacement filling/////////////////////
  for(int i=0;i<nbinX;i++){
    for(int j=0;j<nbinY;j++){
      for(int k=0;k<nbinZ;k++){
	if(XYZdist_fwd[i][j][k].size()){
	  XYZ_Zdist_fwd->SetBinContent(i+1,j+1,k+1,TMath::Median(XYZdist_fwd[i][j][k].size(),&XYZdist_fwd[i][j][k][0]));
	}
	if(XYZdist_fwd[i][j][k].size()==0) XYZ_Zdist_fwd->SetBinContent(i+1,j+1,k+1,-99999999);
      }
    }
  }

  for(int i=0;i<nbinX;i++){
    for(int j=0;j<nbinY;j++){
      for(int k=0;k<nbinZ;k++){
	if(XYZYdist_fwd[i][j][k].size()){
	  XYZ_Ydist_fwd->SetBinContent(i+1,j+1,k+1,TMath::Median(XYZYdist_fwd[i][j][k].size(),&XYZYdist_fwd[i][j][k][0]));
	  //  std::cout<<"i, j, k and content "<<i<<" "<<j<<" k "<<k<<TMath::Median(XYZdist_fwd[i][j][k].size(),&XYZdist_fwd[i][j][k][0])<<std::endl;
	}
	if(XYZYdist_fwd[i][j][k].size()==0) XYZ_Ydist_fwd->SetBinContent(i+1,j+1,k+1,-99999999);
      }
    }
  } 
  for(int i=0;i<nbinX;i++){
    for(int j=0;j<nbinY;j++){
      for(int k=0;k<nbinZ;k++){
	double sizefwd=XYZXdist_fwd[i][j][k].size();
	XYZ_entries_fwd->SetBinContent(i+1,j+1,k+1,sizefwd);
	if(XYZXdist_fwd[i][j][k].size()){
	  XYZ_Xdist_fwd->SetBinContent(i+1,j+1,k+1,TMath::Median(XYZXdist_fwd[i][j][k].size(),&XYZXdist_fwd[i][j][k][0]));
	  //  std::cout<<"i, j, k and content "<<i<<" "<<j<<" k "<<k<<TMath::Median(XYZdist_fwd[i][j][k].size(),&XYZdist_fwd[i][j][k][0])<<std::endl;
	}
	if(XYZXdist_fwd[i][j][k].size()==0) XYZ_Xdist_fwd->SetBinContent(i+1,j+1,k+1,-99999999);
      }
    }
  } 


  ////////////////////////////////////////////////////////


  //  std::cout<<"vector size "<<t_all_tracks.size()<<std::endl;
  deltaT->Write();
  deltaT_pos->Write();
  deltaT_neg->Write();
  deltaT_tanode->Write();
 
  
  dr_all->Write();
  dr_crossing->Write();
  zcov_old->Write();
  zcov_new->Write();
  dr_value_plot->Write();
  deltaT_final->Write();
  shortestdist_true->Write();
  dist_measured->Write();
  dz_vs_x->Draw("colz");
  dz_vs_x->Write();
  t_vs_x->Write();
  std::cout<<"number of stitched earlier, new stitched "<<countold<<"  "<<countnew<<std::endl;

  //////writing distortion histograms here//////
  Xdist_pos_fd->Write();
  Ydist_pos_fd->Write();
  Zdist_pos_fd->Write();
  Xdist_pos_bd->Write();
  Ydist_pos_bd->Write();
  Zdist_pos_bd->Write();
  Xdist_neg_fd->Write();
  Ydist_neg_fd->Write();
  Zdist_neg_fd->Write();
  Xdist_neg_bd->Write();
  Ydist_neg_bd->Write();
  Zdist_neg_bd->Write();

  Xdist_pos_fd_delta->Write();
  Ydist_pos_fd_delta->Write();
  Zdist_pos_fd_delta->Write();
  Xdist_pos_bd_delta->Write();
  Ydist_pos_bd_delta->Write();
  Zdist_pos_bd_delta->Write();
  Xdist_neg_fd_delta->Write();
  Ydist_neg_fd_delta->Write();
  Zdist_neg_fd_delta->Write();
  Xdist_neg_bd_delta->Write();
  Ydist_neg_bd_delta->Write();
  Zdist_neg_bd_delta->Write();


  Zdist_neg_cov->Write();
  Zdist_pos_cov->Write();
  Zdist_neg_cov_bkwd->Write();
  Zdist_pos_cov_bkwd->Write();


  XYZ_Zdist_bkwd->Write();
  XYZ_Ydist_bkwd->Write(); 
  XYZ_Xdist_bkwd->Write();
  XYZ_Zdist_fwd->Write();
  XYZ_Ydist_fwd->Write(); 
  XYZ_Xdist_fwd->Write();
  XYZ_entries_fwd->Write();
  XYZ_entries_bkwd->Write();
  std::cout<<"code finished successfully "<<std::endl;
}//Loop()
