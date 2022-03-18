#define velocityonly_cxx
#include "velocityonly.h"
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
#include <TGraphAsymmErrors.h>
#include <set>
#include <algorithm>
#include <TLatex.h>

int countold=0;
int countnew=0;
double x0=-358.6;//cm
double x1=358.6;//cm
double v_nominal=0.1565;//cm/us
 
double bin_size=20.0;
double y_min=200;
double y_max=400;
int count_1=0;
double z_min=260;
double z_max=440;
//double z_min=125;
//double z_max=225;
double dT_minvalue=9000;//setting lower limit on ACA tracks
double dT_maxvalue=9300;//setting upper limit on ACA tracks


int nbinX=36;
void velocityonly::Loop()
{
 

  double Ymin=0;
  double Ymax=600;
  
  double Zmin=0;
  double Zmax=700;

  int nbinX=36;//40cm bins along all direction
  int nbinY=30;
  int nbinZ=35;
  int nbinYy=30;
  int nbinZy=35;

  TLatex t1[50];

  ///drift velocity plot
  TH1F *X_plot=new TH1F("X_plot","",nbinX,-360,360);
  TH2F *median_vel=new TH2F("median_vel","velocity vs X;X [cm];velocity[mm/us]",nbinX,-360,360,200,1.0,3.0);
  TH1F *med1d_vel=new TH1F("med1d_vel","median drift velocity;X[cm]",nbinX,-360,360);
  TH1F *med_Ef=new TH1F("med_Ef","median Electric field;X[cm]",nbinX,-360,360);
  TH1F *med_ef_measured=new TH1F("med_ef_measured","median Electric field;X[cm]",nbinX,-360,360);
  TH1F *input_Ef=new TH1F("input_Ef","input Electric field;X[cm]",nbinX,-360,360);
  TH1F *input_Ef_same=new TH1F("input_Ef_same","input Electric field with similar coverage;X[cm]",nbinX,-360,360);

  TH1F *deltaX=new TH1F("deltaX","true X-reco X;deltaX[cm];entries",2000,-100,100);
  TH1F *deltaX_cal=new TH1F("deltaX_cal","true X-reco X;deltaX[cm];entries",2000,-100,100);
  TH2F *deltaX_vs_X=new TH2F("deltaX_vs_X","deltaX vs true X",720,-360,360,2000,-100,100);
  TH2F *deltaX_vs_X_cal=new TH2F("deltaX_vs_X_cal","deltaX vs true X",720,-360,360,2000,-100,100);

  TH1F *deltaZ_Y=new TH1F("deltaZ_Y","trueZ-measured Z;deltaZ;entries",2000,-100,100);
  TH1F *X_pos=new TH1F("X_pos","X closest to cathode beam left;X[values];entries",100,-50,50);
  TH1F *X_neg=new TH1F("X_neg","X closest to cathode beam right;X[values];entries",100,-50,50);
  
  TH1F *xmin_plot=new TH1F("xmin_plot","minimum true X;min_X [cm];entries",800,-400,-320);
  TH1F *xmax_plot=new TH1F("xmax_plot","maximum true X;max_X [cm];entries",800,320,400);

  TH2F *deltaZ_X=new TH2F("deltaZ_X","deltaZ vs calculated X;X[cm];deltaZ[cm]",720,-360,360,1000,-50,50);
  TH2F *Zoffset_X=new TH2F("Zoffset_X","Zoffset vs X;X[cm];Zoffset[cm]",720,-360,360,1000,-50,50);

  TH1F *delta_z0=new TH1F("delta_z0","delta z0;deltaZ[cm];entries",1000,-50,50);
  TH1F *delta_z1=new TH1F("delta_z1","delta z1;deltaZ[cm];entries",1000,-50,50);

  TH2F *YZdist_pos=new TH2F("YZdist_pos","Z vs Y coordinate beam left;Z[cm];Y[cm]",1460,0,699.632,1050,0,605);
  TH2F *YZdist_neg=new TH2F("YZdist_neg","Z vs Y coordinate beam right;Z[cm];Y[cm]",1460,0,699.632,1050,0,605);

  TH2F *YZdist_pos_after=new TH2F("YZdist_pos_after","Z vs Y coordinate beam left;Z[cm];Y[cm]",1460,0,699.632,1050,0,605);
  TH2F *YZdist_neg_after=new TH2F("YZdist_neg_after","Z vs Y coordinate beam right;Z[cm];Y[cm]",1460,0,699.632,1050,0,605);


  TH2F *deltaZvsX_pos=new TH2F("deltaZvsX_pos","Z dist vs X value at cathode beam left;deltaZ[cm];X[cm]",80,0,40,400,-200,200);
  TH2F *deltaZvsX_neg=new TH2F("deltaZvsX_neg","Z dist vs X value at cathode beam right;deltaZ[cm];X[cm]",80,0,40,400,-200,200);

  TH1F *delta_tmin=new TH1F("delta_tmin","Time difference between hits closeset to anode;deltaT_anode[ticks];entries",3000,-1500,1500);

  //comment out the 3 lines below if you do not have enough entries to make SCE correction map (1,000,000 entries or more)

  /*  TFile *crossing_dist=new TFile("Zdist_3D_histograms_test1.root");
      TH3F *pos_X=(TH3F*)crossing_dist->Get("XYZ_Zdist_pos");
      TH3F *neg_X=(TH3F*)crossing_dist->Get("XYZ_Zdist_neg");*/


  TH1F *gr_velocity[nbinX];
  TH1F *gr_efieldx[nbinX];
  for(int i=0;i<nbinX;i++){
    double xval=-360+i*20+10.0;
    gr_velocity[i]=new TH1F(Form("gr_%d",i),Form("Graph for X = %f",xval),120,1.0,2.2);
    gr_velocity[i]->GetXaxis()->SetTitle("drift velocity[mm/us]");
    gr_efieldx[i]=new TH1F(Form("gr_%d",i),Form("Graph for X = %f",xval),120,0.22,1.14);
    gr_efieldx[i]->GetXaxis()->SetTitle("Efield_X[kV/cm]");
  }
  TFile *f=new TFile("velocity_100kevents.root","RECREATE");
 
 
  /////////////

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

  
  ///*****************************************************************///
  std::vector<std::vector<double> > velocities,xvalue,tvalue,efields,efieldxs;
  velocities.resize(nbinX);
  efieldxs.resize(nbinX);
  efields.resize(nbinX);



  std::vector<double> t_neg,t_pos,z_pos,z_neg,y_pos,y_neg, t_all,z_all,y_all,tpc_pos,tpc_neg,tpc_all, x_pos, x_neg, x_all, x_pos_true, x_neg_true, x_all_true, z_pos_true, z_neg_true, z_all_true, y_pos_true, y_neg_true, y_all_true;
  std::vector<int> event_number, run_number, entry_number;

  std::vector<std::vector<double> > t_neg_tracks,z_neg_tracks,y_neg_tracks,tpc_neg_tracks,t_pos_tracks,z_pos_tracks,y_pos_tracks,tpc_pos_tracks,t_all_tracks,z_all_tracks,y_all_tracks,tpc_all_tracks;

  // std::vector<double> testX, testT; 
  TCanvas *canv[50];
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  //for (Long64_t jentry=0;jentry<nentries;jentry++) {
    for(Long64_t jentry=0;jentry<=200000;jentry++){
    Long64_t ientry = LoadTree(jentry);
    if(jentry%10000==0) std::cout<<nentries<<"/"<<jentry<<std::endl;
    //std::cout<<"Test 0 "<<std::endl;
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
   
    if(!trkhitz_wire2->size()) continue;

    //////////////
    t_neg_tracks.clear();
    z_neg_tracks.clear();
    y_neg_tracks.clear();
    tpc_neg_tracks.clear();
       
    t_pos_tracks.clear();
    z_pos_tracks.clear();
    y_pos_tracks.clear();
    tpc_pos_tracks.clear();

    t_all_tracks.clear();
    z_all_tracks.clear();
    y_all_tracks.clear();
    tpc_all_tracks.clear();


    /////////////

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
      // std::cout<<"tpos size, tneg size "<<t_pos_buff.size()<<" "<<t_neg_buff.size()<<std::endl;
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
      if(!(abs(t_min_pos-t_min_neg)<10 && dT>dT_minvalue && dT<dT_maxvalue && t_max_pos<5985 && t_max_neg<5985)) continue;
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
      event_number.push_back(event);
      run_number.push_back(run);
      t_all_tracks.push_back(t_all);
      entry_number.push_back(jentry);
      z_all_tracks.push_back(z_all);
      y_all_tracks.push_back(y_all);
      tpc_all_tracks.push_back(tpc_all);
      countold++;
      zcov_old->Fill(z_pos[0],y_pos[0]);

    }//i loop

    if(t_neg_buff.size() && t_pos_buff.size()){
      for(size_t t=0;t<t_neg_buff.size();t++){
	for(size_t t1=0;t1<t_pos_buff.size();t1++){
	  double t_max_pos=*std::max_element(t_pos_buff[t1].begin(),t_pos_buff[t1].end());
	  double t_min_pos=*std::min_element(t_pos_buff[t1].begin(),t_pos_buff[t1].end());
	  double t_max_neg=*std::max_element(t_neg_buff[t].begin(),t_neg_buff[t].end());
	  double t_min_neg=*std::min_element(t_neg_buff[t].begin(),t_neg_buff[t].end());
	  double dT=t_max_pos-t_min_pos+t_max_neg-t_min_neg;
	  delta_tmin->Fill(t_min_pos-t_min_neg);
	  if(!(abs(t_min_pos-t_min_neg)<10 && dT>dT_minvalue && dT<dT_maxvalue && t_max_pos<5985 && t_max_neg<5985)) continue;
	 
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
	  event_number.push_back(event);
	  run_number.push_back(run);
	  entry_number.push_back(jentry);
	  tpc_all_tracks.push_back(tpc_all_b);
	  zcov_new->Fill(z_pos_buff[t1][0],y_pos_buff[t1][0]);
	  countnew++;
	}
      }
    }
    //  std::cout<<"tracks number is "<<t_all_tracks.size()<<std::endl;
    ///////////////////////calculating velocity for each of the anode-anode tracks
    ///////////measuring the drift velocity here////////////////////
    std::vector<double> t_all_new,z_all_new,y_all_new,x_all_new,tpc_all_new;
    std::vector<double> xvalues_all_new;
    // if(t%1000==0) std::cout<<"tracks "<<t<<"/"<<t_all_tracks.size()<<std::endl;
    if(t_all_tracks.size()==0) continue;
    for(int t=0;t<t_all_tracks.size();t++){
      t_all_new.clear();z_all_new.clear();x_all_new.clear();y_all_new.clear();tpc_all_new.clear();
      if(t_all_tracks[t].size()<10) continue;
      if(t_pos_tracks[t].size()<5) continue;
      if(t_neg_tracks[t].size()<5) continue;
      //std::cout<<"tracks, event number, run number "<<t<<"/"<<t_all_tracks.size()<<" "<<event_number[t]<<"  "<<run_number[t]<<" jentry "<<entry_number[t]<<std::endl;
      // std::cout<<"test1 "<<std::endl;
      xvalue.clear();
      tvalue.clear();
      xvalue.resize(nbinX);
      tvalue.resize(nbinX);
   
      // TestT.clear();TestX.clear();

      // std::cout<<"test2 "<<std::endl;
      ////velocity code
      t_all_new.push_back(t_all_tracks[t][0]);
      z_all_new.push_back(z_all_tracks[t][0]);
      y_all_new.push_back(y_all_tracks[t][0]);
      tpc_all_new.push_back(tpc_all_tracks[t][0]);
      //  std::cout<<"test3 "<<std::endl;
 
      for(size_t s1=1;s1<t_all_tracks[t].size()-1;s1++){
	bool test=true;
	//	std::cout<<"test4 "<<std::endl;
	for(size_t s2=0;s2<t_all_tracks[t].size();s2++){
	  if(z_all_tracks[t][s2]==z_all_tracks[t][s1] && s2!=s1){
	    test=false;
	    break;
	  }
	}
	if(test){
	  // std::cout<<"test5 "<<std::endl;
	  t_all_new.push_back(t_all_tracks[t][s1]);
	  z_all_new.push_back(z_all_tracks[t][s1]);
	  y_all_new.push_back(y_all_tracks[t][s1]);
	  tpc_all_new.push_back(tpc_all_tracks[t][s1]);

	}
      }
      // std::cout<<"test6 "<<std::endl;
      t_all_new.push_back(t_all_tracks[t][t_all_tracks[t].size()-1]);
      z_all_new.push_back(z_all_tracks[t][t_all_tracks[t].size()-1]);
      y_all_new.push_back(y_all_tracks[t][t_all_tracks[t].size()-1]);
      tpc_all_new.push_back(tpc_all_tracks[t][t_all_tracks[t].size()-1]);
      // std::cout<<"test7 "<<std::endl;
    
      double z0=z_all_new[0];
      double z1=z_all_new[z_all_new.size()-1];
      if(abs(z0-z1)<30) continue;     
      // std::cout<<"test8 "<<std::endl;
      //Calculate distortion at cathode here
      double y0=y_all_new[0];
      double y1=y_all_new[y_all_new.size()-1];
      double x_cathode=0;
      xvalues_all_new.clear();
      // std::cout<<"step 9"<<std::endl;

      for(size_t it=0;it<t_all_new.size();it++){
	// 	std::cout<<"test10 "<<std::endl;
	//  std::cout<<"n hits in a track "<<t_all_new.size()<<std::endl;
	if(z_all_new[it]>z_min && z_all_new[it]<z_max && y_all_new[it]>y_min && y_all_new[it]<y_max){
	  if(it<5||it>t_all_new.size()-5) continue;
	  double xcal=0;
	  double zcorrection=0;
	  // std::cout<<"test11 "<<std::endl;
	  if((tpc_all_new[it]==2 ||tpc_all_new[it]==6||tpc_all_new[it]==10)) xcal=x1-(t_all_new[it]-t_pos_tracks[t][t_pos_tracks[t].size()-1])*0.5*v_nominal;
	  if((tpc_all_new[it]==1 ||tpc_all_new[it]==5||tpc_all_new[it]==9)) xcal=x0+(t_all_new[it]-t_neg_tracks[t][0])*0.5*v_nominal;
	  if(xcal<-360||xcal>360) continue;
	 
	  //comment the two lines below if you do not have SCE maps
	  //	  if(tpc_all_new[it]==1 ||tpc_all_new[it]==5||tpc_all_new[it]==9) zcorrection=neg_X->Interpolate(xcal,y_all_new[it],z_all_new[it]);
	  //	  if(tpc_all_new[it]==2 ||tpc_all_new[it]==6||tpc_all_new[it]==10) zcorrection=pos_X->Interpolate(xcal,y_all_new[it],z_all_new[it]);
        
	  double zcorrected=z_all_new[it]+zcorrection;
	  double xvalues=x1+(x1-x0)*(zcorrected-z1)/(z1-z0);
	  double zvalues=z1+(z1-z0)*(y_all_new[it]-y1)/(y1-y0);
	  xvalues_all_new.push_back(xvalues);
	  int xbin=X_plot->FindBin(xvalues);
	  // std::cout<<"test13 "<<std::endl;
	  if(((tpc_all_new[it]==2 ||tpc_all_new[it]==6||tpc_all_new[it]==10) && xbin<19)/*||abs(zdistpos)>30.0*/) continue;
	  if(((tpc_all_new[it]==1 ||tpc_all_new[it]==5||tpc_all_new[it]==9) && xbin>18)/*|| abs(zdistneg)>30.0*/) continue;
	  int counter=0;
	  if(xbin<1||xbin>nbinX) continue;
	  efields[xbin-1].push_back(Efield_input(xvalues,y_all_new[it],z_all_new[it]));
	  xvalue[xbin-1].push_back(xvalues);
	  tvalue[xbin-1].push_back(t_all_new[it]);
	  
	  // std::cout<<"test14 "<<std::endl;
	}//if z<zmin.... loop
      }//it loop
      // std::cout<<"test15 "<<std::endl;
      std::vector<double> xbuff, xbuff1, tbuff, tbuff1;
      for(int k=0;k<nbinX;k++){
	xbuff.clear();xbuff1.clear();tbuff.clear();tbuff1.clear();
	//	std::cout<<"test16 "<<std::endl;
	if(xvalue[k].size()>5){
	  //Try some data cleaning, look remove part of tracks with repeating Z and also time and Z not in order
	  for(int i1=0;i1<xvalue[k].size();i1++){
	    xbuff.push_back(xvalue[k][i1]);
	    tbuff.push_back(tvalue[k][i1]);
	  }
	  //  std::cout<<"test17 "<<std::endl;
	  set<double> s(xbuff.begin(),xbuff.end()); //checking if the entries are repeated
	  if(s.size()!=xbuff.size()) continue;
	  xbuff1=xbuff;
	  tbuff1=tbuff;
	  std::sort(xbuff.begin(),xbuff.end());
	  if(xbuff1!=xbuff) continue;
	  if(tbuff[0]<tbuff[tbuff.size()-1]){
	    tbuff1=tbuff;
	    std::sort(tbuff.begin(),tbuff.end());
	    if(tbuff1!=tbuff) continue;
	  }
	  if(tbuff[0]>tbuff[tbuff.size()-1]){
	    tbuff1=tbuff;
	    std::sort(tbuff.rbegin(),tbuff.rend());
	    if(tbuff1!=tbuff) continue;
	  }
	  if(tvalue[k].size()<2) continue;
  
	  TGraph *gr=new TGraph(tvalue[k].size(),&tvalue[k][0],&xvalue[k][0]);
	  gr->Fit("pol1","Q");
	  TF1 *f1=gr->GetFunction("pol1");
	  double slope=f1->GetParameter(1);
	  delete gr;
	  velocities[k].push_back(abs(20.0*slope));
	  gr_velocity[k]->Fill(abs(20.0*slope));
	  if(abs(20.0*slope)>1.0 && abs(20.0*slope)<2.2){
	    gr_efieldx[k]->Fill(driftvel(abs(20.0*slope)));
	    efieldxs[k].push_back(driftvel(abs(20.0*slope)));
	  }
	  median_vel->Fill(-360+k*20+10.0,abs(20.0*slope));
	}
      }//Nbinx
    }//t loop tracks
    ///////////////////////////////////////////////////////////////////////////
  }//jentry loop
 
  std::vector<double> xerr, fitted_vel, std_err, xposition, measured_efield,std_err_ef; 
  for(int i=0;i<nbinX;i++){
    // std::cout<<"velocities size , nbin "<<velocities[i].size()<<"  "<<i<<std::endl;
    if(velocities[i].size()){
      med1d_vel->SetBinContent(i+1,TMath::Median(velocities[i].size(),&velocities[i][0]));
      med_Ef->SetBinContent(i+1,driftvel(TMath::Median(velocities[i].size(),&velocities[i][0])));
      med_ef_measured->SetBinContent(i+1,TMath::Median(efieldxs[i].size(),&efieldxs[i][0]));
      input_Ef->SetBinContent(i+1,Efield_input(-360+i*20+10.0,300,347.5));
      input_Ef_same->SetBinContent(i+1,TMath::Median(efields[i].size(),&efields[i][0]));
      // std::cout<<"bin vs median "<<i<<"  "<<TMath::Median(velocities[i].size(),&velocities[i][0])<<std::endl;
      gr_efieldx[i]->Write();
      gr_velocity[i]->Write();
      gr_velocity[i]->Fit("gaus","Q");
      TF1 *fun=gr_velocity[i]->GetFunction("gaus");
      xposition.push_back(-360+i*20+10.0);
      xerr.push_back(0.0);
      fitted_vel.push_back(driftvel(fun->GetParameter(1)));
      std_err.push_back(0.31*fun->GetParError(1));
      gr_efieldx[i]->Write();
      gr_efieldx[i]->Fit("gaus","Q");
      TF1 *fune=gr_efieldx[i]->GetFunction("gaus");
      measured_efield.push_back(fune->GetParameter(1));
      std_err_ef.push_back(fune->GetParError(1));
    }
  }//nbinX
  TGraphErrors *fittedefield=new TGraphErrors(fitted_vel.size(),&xposition[0],&fitted_vel[0],&xerr[0],&std_err[0]);
  TGraphErrors *fitted_meas_ef=new TGraphErrors(measured_efield.size(),&xposition[0],&measured_efield[0],&xerr[0],&std_err_ef[0]);


  ////New TGraphAsymmErrors including error on median values
  double vlow[nbinX];double vhigh[nbinX];double vmed[nbinX];double xpos[nbinX];double errx[nbinX];
  double Elow[nbinX];double Ehigh[nbinX];double Emed[nbinX];
  for(int i=0;i<nbinX;i++){
    if(velocities[i].size()){
      double qout[3];
      double p[3];
      p[0]=0.5-0.5/sqrt(velocities[i].size()+2);
      p[1]=0.5;
      p[2]=0.5+0.5/sqrt(velocities[i].size()+2);
      TMath::Quantiles(velocities[i].size(),3,&velocities[i][0],qout,p,0,0,7);
      vlow[i]=qout[1]-qout[0];
      vhigh[i]=qout[2]-qout[1];
      vmed[i]=qout[1];
      //Efields
      double qout1[3];
      double p1[3];
      p1[0]=0.5-0.5/sqrt(efieldxs[i].size()+2);
      p1[1]=0.5;
      p1[2]=0.5+0.5/sqrt(efieldxs[i].size()+2);
      TMath::Quantiles(efieldxs[i].size(),3,&efieldxs[i][0],qout1,p1,0,0,7);
      Elow[i]=qout1[1]-qout1[0];
      Ehigh[i]=qout1[2]-qout1[1];
      Emed[i]=qout1[1];
      xpos[i]=xposition[i];
      errx[i]=10;
    }
  }
  ////////////////////////////////////////////////////
  TGraphAsymmErrors *vel_med_error=new TGraphAsymmErrors(nbinX,xpos,vmed,errx,errx,vlow,vhigh);
  TGraphAsymmErrors *Ef_med_error=new TGraphAsymmErrors(nbinX,xpos,Emed,errx,errx,Elow,Ehigh);


  deltaT->Write();
  deltaT->Draw();
  median_vel->Write();
  median_vel->Draw("colz");
  med1d_vel->Write();
  input_Ef->Write();
  med_Ef->Write();
  input_Ef_same->Write();
  fittedefield->Write("fitted_efield");
  deltaX->Write();
  deltaX_cal->Write();
  deltaX_vs_X->Write();
  deltaX_vs_X_cal->Write();
  deltaZ_Y->Write();
  deltaZ_X->Write();
  X_pos->Write();
  X_neg->Write();
  xmin_plot->Write();
  xmax_plot->Write();
  Zoffset_X->Write();
  delta_z0->Write();
  delta_z1->Write();
  deltaZvsX_pos->Write();
  deltaZvsX_neg->Write();
  YZdist_pos->Write();
  YZdist_neg->Write();
  YZdist_pos_after->Write();
  YZdist_neg_after->Write();
  fitted_meas_ef->Write();
  med_ef_measured->Write();
  vel_med_error->Write("vel_med_error");
  Ef_med_error->Write("Ef_med_error");
  delta_tmin->Write();
  std::cout<<"code finished successfully "<<std::endl;
  f->Close();
}//Loop()
