BEGIN{
# This program is used to extract the same distribution 
# of P and S wave data
  num=0;
  ip=0;
  is=0;
}
{
  # read absolute data or differential data
  #determine if num==0
  if($1=="#"){
    if(num==0){
      event=$0;
    }
    else{ #num!=0
      nump=ip;
      nums=is;
      if(nump+nums!=num){
	print "Something is wrong in reading data";
      }
      else{
	if(nums==0){
	  # No S data. skip this event or event pair
	  num=0;
	  is=0;
	  ip=0;
	  nump=0;
	  nums=0;
          event=$0;
	}
	else{
	  # do have some S data. then choose P-wave data according to S data
	  print event > "dt_sp.ct";
	  for(i=1;i<=nums;i++){
	    #print data_S[i] > "dtct_repicked_7644_sp.dat";
	    for(j=1;j<=nump;j++){
	      if(staP[j]==staS[i]){
		print staP[j],S1[i]-S2[i]-P1[j]+P2[j],(S_wgh[i]+P_wgh[j])/2>"dt_sp.ct";		
		#print data_P[j] > "dtct_repicked_7644_sp.dat";
	      }
	    }
	  }
	  num=0;
	  is=0;
	  ip=0;
	  nump=0;
	  nums=0;
	  event=$0;
	}
      }
    }
  }
  else{ #read station data
    if($NF=="S"){
      # read S data
      is=is+1;
      data_S[is]=$0;
      staS[is]=$1;
      S1[is]=$2;
      S2[is]=$3;
      S_wgh[is]=$4;
      num=num+1;
    }
    else{
      # read P data
      ip=ip+1;
      data_P[ip]=$0;
      staP[ip]=$1;
      P1[ip]=$2;
      P2[ip]=$3;
      P_wgh[ip]=$4;
      num=num+1;
    }
  }
}
END{
  if(is>0) {    
    print event > "dt_sp.ct";
    for(i=1;i<=is;i++){
      for(j=1;j<=ip;j++){
	if(staP[j]==staS[i]){
	  print staP[j],S1[i]-S2[i]-P1[j]+P2[j],(S_wgh[i]+P_wgh[j])/2>"dt_sp.ct";		
	}
      }
    }
  }
}
