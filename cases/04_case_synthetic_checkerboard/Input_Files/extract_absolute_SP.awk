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
      #print event;
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
	  print event > "absolute_sp.dat";
	  for(i=1;i<=nums;i++){
	    #print data_S[i] > "absolute_repicked_7644_sp.dat";
	    for(j=1;j<=nump;j++){
	      if(staP[j]==staS[i]){
		print staP[j],S_trav[i]-P_trav[j],(S_wgh[i]+P_wgh[j])/2>"absolute_sp.dat";		
		#print data_P[j] > "absolute_repicked_7644_sp.dat";
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
       #phase=$NF;
       #print phase;
    if($NF=="S"){
      # read S data
      is=is+1;
      data_S[is]=$0;
      staS[is]=$1;
      S_trav[is]=$2;
      S_wgh[is]=$3;
      num=num+1;
      #print is,$0;
    }
    else{
      # read P data
      ip=ip+1;
      data_P[ip]=$0;
      staP[ip]=$1;
      P_trav[ip]=$2;
      P_wgh[ip]=$3;
      num=num+1;
      #print num,phase,$4;
    }
  }
}
END{
  if(is>0) {    
    print event > "absolute_sp.dat";
    for(i=1;i<=is;i++){
      for(j=1;j<=ip;j++){
	if(staP[j]==staS[i]){
	  print staP[j],S_trav[i]-P_trav[j],(S_wgh[i]+P_wgh[j])/2>"absolute_sp.dat";		
	}
      }
    }
  }
}
