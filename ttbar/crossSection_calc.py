import ROOT

#Define function to calc signal strength: XSBR(data)/XSBR(expected)
def mu_calc( N_data, N_bkg, N_sig ): return (N_data-N_bkg)/N_sig

#Function to calculate uncertainty on mu
# N_sig = XSBR(theory)*eff*Acc*Lumi
# Here u_sig is a vector of systematic uncertainties related to MC signal (eleID,muID,QCDscale,PU,lumi etc)
def mu_unc_calc( N_data, u_data, N_bkg, u_bkg, N_sig, u_sig_list ): #
  
  #Sum up the uncertainties in quadrature
  u_sig_combined = 0
  for unc in u_sig_list: u_sig_combined+=unc*unc
  u_sig_combined = u_sig_combined**0.5

  u_mu = mu_calc(N_data,N_bkg,N_sig)*(((u_data**2+u_bkg**2)/(N_data-N_bkg)**2)+(u_sig_combined)**2)**0.5
  return u_mu
  
# Function to read input file from systematics_calc
def readSystematics( channel="temp" ):
  
  #Initiate objects
  N_data_ = 0
  u_data_ = 0
  N_bkg_ = 0
  u_bkg_ = 0
  N_sig_ = 0
  uStat_sig_ = 0
  u_sig_ = []

  f_input = open("systematics_%s.txt"%channel,"r")
  for line in f_input:
    #split line by spaces
    line_info = line.split(" ")
    
    if( line_info[0] == 'channel' ): continue
    elif( line_info[0] == 'data' ): 
      N_data_ = float( line_info[1] )
      u_data_ = float( line_info[2] )
    elif( line_info[0] == 'bkg' ): 
      N_bkg_ = float( line_info[1] )
      u_bkg_ = float( line_info[2] )
    elif( line_info[0] == "sig" ):
      N_sig_ = float( line_info[1] )
      uStat_sig_ = float( line_info[2] )
    elif( line_info[0][0] == "-" ): continue
    # Extract systematics and add to list
    else:
      u_sig_.append( float( line_info[1] ) )

  #Add fractional MC stat uncertainty to list
  u_sig_.append( uStat_sig_/N_sig_ )

  #Run mu calc and get uncertainty
  mu = mu_calc( N_data_, N_bkg_, N_sig_ )
  u_mu = mu_unc_calc( N_data_, u_data_, N_bkg_, u_bkg_, N_sig_, u_sig_ )

  print "Channel %s"%channel
  print "  * mu = %5.4f +- %5.4f"%(mu,u_mu) 

readSystematics( channel="temp" )

    

