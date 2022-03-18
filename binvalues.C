void binvalues(){
  //This script can be run using root -l binvalues.C
  //This uses the 3D spatial distortion map made usign anode-cathode-anode crossing tracks and and combines it with the distortion at the cathode; so as to make it continuous across the cathode; and 3D interpolation can be used 
  //The map is good for central region of the TPC;boundary (Z=240-450 cm; Y=200-400 cm); it is also useful for regions outside but close the boundary.




  TFile *f=new TFile("/dune/app/users/apaudel/code/velocity_code_tests/data_anode_3D_oct25_bcwd_fwd.root");
  TH3F *xyzdist=(TH3F*)f->Get("XYZ_Zdist_bkwd"); //3D spatial distortion using anode-cathode-anode tracks

  TH2F *negZ=(TH2F*)f->Get("Zdist_neg_bd");//distortion at cathod X<0 using using cathode-anode tracks
  TH2F *posZ=(TH2F*)f->Get("Zdist_pos_bd");//distortion at cathod X>0 using using cathode-anode tracks

  TFile *newfile=new TFile("Zdist_3D_histograms_test.root","RECREATE");//output root file with the final 3D Z distortion map
  TH3F *XYZ_Zdist_pos=new TH3F("XYZ_Zdist_pos","3D Z dist pos;X[cm];Y[cm];Z[cm];deltaZ[cm]",11,-40,400,30,0,600,35,0,700);
  TH3F *XYZ_Zdist_neg=new TH3F("XYZ_Zdist_neg","3D Z dist neg;X[cm];Y[cm];Z[cm];deltaZ[cm]",11,-400,40,30,0,600,35,0,700);
  TH1F *zvsx=new TH1F("zvsx","",36,-360,360);

  for(int i=0;i<11;i++){
    for(int j=0;j<30;j++){
      for(int k=0;k<35;k++){
	if(i==0){
	  XYZ_Zdist_pos->SetBinContent(i+1,j+1,k+1,2*posZ->GetBinContent(k+1,j+1)-xyzdist->GetBinContent(10,j+1,k+1));//bin near anode, at the anode distortion is assumed to be 0
	  // std::cout<<"neg bin1 value "<<-(xyzdist->GetBinContent(1,j+1,k+1));
	  XYZ_Zdist_neg->SetBinContent(i+1,j+1,k+1,-(xyzdist->GetBinContent(1,j+1,k+1)));
	  continue;
	}
	if(i==10){
	  XYZ_Zdist_neg->SetBinContent(i+1,j+1,k+1,2*negZ->GetBinContent(k+1,j+1)-xyzdist->GetBinContent(9,j+1,k+1));//bin near cathode 
	  XYZ_Zdist_pos->SetBinContent(i+1,j+1,k+1,-xyzdist->GetBinContent(18,j+1,k+1));
	  continue;
	}
	  XYZ_Zdist_pos->SetBinContent(i+1,j+1,k+1,xyzdist->GetBinContent(9+i,j+1,k+1));
	  XYZ_Zdist_neg->SetBinContent(i+1,j+1,k+1,xyzdist->GetBinContent(i,j+1,k+1));

      }
    }
  }



  //just an exampe how to extract the results
  for(int i=0;i<36;i++){
    //if(xyzdist->GetBinContent(i,7,13)==-9999) continue;
    if(i<18)  zvsx->SetBinContent(i+1,XYZ_Zdist_neg->Interpolate(-360+10+i*20,300,400));
    if(i>=18) zvsx->SetBinContent(i+1,XYZ_Zdist_pos->Interpolate(10+(i-18)*20,300,400));
     //std::cout<<"i "<<i<<" "<<xyzdist->GetBinContent(i,12,20)<<std::endl;
  }
  zvsx->Draw();
  //end of example
  XYZ_Zdist_neg->Write();
  XYZ_Zdist_pos->Write();
  newfile->Close();	
}
