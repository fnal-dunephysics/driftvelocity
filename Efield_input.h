#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h> 
#include <string>
#include <TImage.h>
#include <iomanip>
#include <algorithm>
//TFile *ef=new TFile("SCE_DataDriven_180kV_v3.root");
//TFile *ef=new TFile("/dune/data/users/apaudel/R5770_DataDriven_180kV_v2.root");
TFile *ef=new TFile("/cvmfs/dune.opensciencegrid.org/products/dune/dune_pardata/v01_65_01/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");
TH3F *xneg=(TH3F*)ef->Get("Reco_ElecField_X_Neg");
TH3F *yneg=(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
TH3F *zneg=(TH3F*)ef->Get("Reco_ElecField_Z_Neg");
TH3F *xpos=(TH3F*)ef->Get("Reco_ElecField_X_Pos");
TH3F *ypos=(TH3F*)ef->Get("Reco_ElecField_Y_Pos");
TH3F *zpos=(TH3F*)ef->Get("Reco_ElecField_Z_Pos");


TH3F *zpos_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Pos");
TH3F *zneg_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Neg");
TH3F *zpos_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Pos");
TH3F *zneg_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Neg");
TH3F *ypos_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Y_Pos");
TH3F *yneg_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Y_Neg");
TH3F *ypos_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Y_Pos");
TH3F *yneg_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Y_Neg");

TH3F *xpos_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_X_Pos");
TH3F *xneg_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_X_Neg");

TH3F *xpos_fd=(TH3F*)ef->Get("RecoFwd_Displacement_X_Pos");
TH3F *xneg_fd=(TH3F*)ef->Get("RecoFwd_Displacement_X_Neg");


float zoffsetfd(float xval,float yval,float zval){
  if(xval>=0){
    return zpos_fd->Interpolate(xval,yval,zval);
  }
if(xval<0){
  return zneg_fd->Interpolate(xval,yval,zval);
  }
 return 0;
}

float zoffsetbd(float xval,float yval,float zval){
   if(xval>=0){
     //return zpos_bd->GetBinContent(zpos_bd->FindBin(xval,yval,zval));
     return zpos_bd->Interpolate(xval,yval,zval);
  }
if(xval<0){
  // return zneg_bd->GetBinContent(zneg_bd->FindBin(xval,yval,zval));
  return zneg_bd->Interpolate(xval,yval,zval);
  }
 return 0;
}



float xoffsetfd(float xval,float yval,float zval){
  if(xval>=0){
    return xpos_fd->Interpolate(xval,yval,zval);
  }
if(xval<0){
  return xneg_fd->Interpolate(xval,yval,zval);
  }
 return 0;
}


float xoffsetbd(float xval,float yval,float zval){
   if(xval>=0){
     //return zpos_bd->GetBinContent(zpos_bd->FindBin(xval,yval,zval));
     return xpos_bd->Interpolate(xval,yval,zval);
  }
if(xval<0){
  // return zneg_bd->GetBinContent(zneg_bd->FindBin(xval,yval,zval));
  return xneg_bd->Interpolate(xval,yval,zval);
  }
 return 0;
}


float yoffsetfd(float xval,float yval,float zval){
  if(xval>=0){
    return ypos_fd->Interpolate(xval,yval,zval);
  }
if(xval<0){
  return yneg_fd->Interpolate(xval,yval,zval);
  }
 return 0;
}




float yoffsetbd(float xval,float yval,float zval){
   if(xval>=0){
     //return zpos_bd->GetBinContent(zpos_bd->FindBin(xval,yval,zval));
     return ypos_bd->Interpolate(xval,yval,zval);
  }
if(xval<0){
  // return zneg_bd->GetBinContent(zneg_bd->FindBin(xval,yval,zval));
  return yneg_bd->Interpolate(xval,yval,zval);
  }
 return 0;
}





float Efield_input(float xval,float yval,float zval){
  float E0value=0.4867;
  if(xval>=0){
    //float ex=E0value+E0value*xpos->GetBinContent(xpos->FindBin(xval,yval,zval));
    float ex=E0value+E0value*xpos->Interpolate(xval,yval,zval);
    float ey=0.0+E0value*ypos->Interpolate(xval,yval,zval);
    float ez=0.0+E0value*zpos->GetBinContent(zpos->FindBin(xval,yval,zval));
    //return sqrt(ex*ex+ey*ey+ez*ez);
     return ex;
  }
  if(xval<0){
    float ex=E0value+E0value*xneg->Interpolate(xval,yval,zval);
    float ey=0.0+E0value*yneg->GetBinContent(yneg->FindBin(xval,yval,zval));
    float ez=0.0+E0value*zneg->GetBinContent(zneg->FindBin(xval,yval,zval));
    // return sqrt(ex*ex+ey*ey+ez*ez);
    return ex;
  }
  return E0value;
}
