//This code can be run using the following steps
//.L deltaT_plots.C+
//deltaT_plots k
//k.Loop()
//This code uses the anode-cathode-anode crossing tracks to find the Spatial distortion in 3D bins across the TPC
//https://lss.fnal.gov/archive/thesis/2000/fermilab-thesis-2021-02.pdf [For details of the method, chapter 4]

#define deltaT_plots_cxx
#include "deltaT_plots.h"
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

void deltaT_plots::Loop()
{
  
  double Ymin=0;
  double Ymax=600;
  
  double Zmin=0;
  double Zmax=700;

  int nbinX=18;//40cm bins along all direction
  int nbinY=30;
  int nbinZ=35;
  int nbinYy=30;
  int nbinZy=35;

 
  TH1D *dZ_pos[35];
  TH1D *dZ_neg[35];
  for(int i=0;i<35;i++){
    dZ_pos[i]=new TH1D(Form("dZ_pos_%d",i),"",200,-100,100);
    dZ_neg[i]=new TH1D(Form("dZ_neg_%d",i),"",200,-100,100);
  }

  TH1F *deltaT=new TH1F("deltaT","PeakT_max-PeakT_min;deltaT;entries",5000,0,10000);
  TH1F *deltaT_after_stitching=new TH1F("deltaT_after_stitching","PeakT_max-PeakT_min;deltaT;entries",5000,0,10000);
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

  TFile *f=new TFile("plots_for_making_deltaT_cut.root","RECREATE");

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
      deltaT->Fill(dT);
      deltaT_after_stitching->Fill(dT);
      if(!(abs(t_min_pos-t_min_neg)<10 && dT>9000 && dT<9250 && t_max_pos<5985 && t_max_neg<5985)) continue;
      //if(!(dT>9120 && dT<9160)) continue;
      //std::cout<<"event min max later dT "<<event<<"  "<<t_min_pos<<"  "<<t_min_neg<<"  "<<t_max_pos<<"  "<<t_max_neg<<"  "<<dT<<" dr "<<dr_value<<std::endl;
      dr_crossing->Fill(dr_value);
    
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
   	if(!(abs(t_min_pos-t_min_neg)<10 && dT>9000 && dT<9250 && t_max_pos<5985 && t_max_neg<5985)) continue;
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
	deltaT_after_stitching->Fill(dT);
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
        
  ////////////////////////////////////////////////////////

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
  deltaT_after_stitching->Write();
  shortestdist_true->Write();
  dist_measured->Write();
  // dz_vs_x->Draw("colz");
  //  dz_vs_x->Write();
  // t_vs_x->Write();
  std::cout<<"number of stitched earlier, new stitched "<<countold<<"  "<<countnew<<std::endl;
 
  f->Close();
  std::cout<<"code finished successfully "<<std::endl;
}//Loop()
