import ROOT

# Global variables
intLumi = 296.08 # pb^-1, uncertainty = 3.5 %

#Channels: ee, mumu, emu
channelString = "me"

#Open file to read yield histograms
f_input = ROOT.TFile("/home/users/lguzzi/LongEx/CMSSW_10_2_10/src/ttCMSDAS/ttbar/temp/nominal/TT.root")

#Define systematic dictionary
syst_dict = {"TT%s_pdf" %channelString:'Y%s_Pdf'%channelString, "TT%s_alpha" %channelString:'Y%s_Alphas'%channelString, "TT%s_muonID" %channelString:"Y%s_SFmuon"%channelString, "TT%s_electronID"%channelString:"Y%s_SFelec"%channelString, "TT%s_QCDscale"%channelString:"Y%s_MatrixEl"%channelString, "TT%s_pileup"%channelString:"Y%s_SFPU"%channelString}

#Dictionary to hold central uncertainty
syst_values = {"TT_lumi":0.035}

#For nicer output to group D
syst_strings = {"TT_lumi":"1.035 0.965"}

#Extract nominal events from first histogram
h_initialize = f_input.Get( syst_dict["TT%s_muonID"%channelString] )
n_nominal = h_initialize.GetBinContent(1)
n_nominal_statUnc = h_initialize.GetBinError(1)

print "This is for channel: %s"%channelString
print "N_signal = %5.4f +- %5.4f"%(n_nominal*intLumi,n_nominal_statUnc*intLumi)
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#Iterate over keys, values in dictionary and extract the yield variations
for syst, hSyst_string in syst_dict.iteritems():

  #Initiate string to hold yield variations for systematic
  syst_string = "%s "%syst

  #Get the corresponding histogram
  h_syst = f_input.Get("%s"%hSyst_string)

  #For "specific" systematics: treat differently
  if( syst == "TT_QCDscale" ):
    #Loop over bins and extract n_up as maximum and n_down as minimum (not equal to zero)
    #Initialise n_up and n_down to nominal
    n_up = n_nominal
    n_down = n_nominal
    for bin_idx in range(1,h_syst.GetNbinsX()+1):
      if h_syst.GetBinContent(bin_idx)>n_up: n_up = h_syst.GetBinContent(bin_idx)
      if h_syst.GetBinContent(bin_idx)<n_down and h_syst.GetBinContent(bin_idx)!= 0: n_down = h_syst.GetBinContent(bin_idx)
      
  #For statdard (up/down) systematics
  else:
    #Get number of events with up and down variation
    n_up = h_syst.GetBinContent(2)
    n_down = h_syst.GetBinContent(3)
  
  #Yield variations (percentage)
  yieldVar_up = n_up/n_nominal
  yieldVar_down = n_down/n_nominal
  yieldVar_central = (0.5*abs(n_up-n_down))/n_nominal

  #Add value to dictionary
  syst_values['%s'%syst] = yieldVar_central
  syst_strings['%s'%syst] = "%5.4f %5.4f"%(yieldVar_up,yieldVar_down)

  #Add information to string
  syst_string += "%5.4f/%5.4f (%4.2f%%/%4.2f%%) (central = %4.2f%%)"%(yieldVar_up,yieldVar_down,abs(yieldVar_up-1)*100,abs(yieldVar_down-1)*100,(yieldVar_central*100))

  #Print the yielf variations for the systematic
  print syst_string
    
print ""
print "############################################################"
print " NICE OUTPUT FOR THE GROUP D GUYS                           "
print ""
print " CHANNEL: %s"%channelString
print ""
print " Nominal Signal Yield = %5.4f +- %5.4f"%(n_nominal*intLumi,n_nominal_statUnc*intLumi)
print ""
print " SYSTEMATCS UNCERTAINTY YIELD VARIATIONS:"
print ""
for key, value in syst_strings.iteritems():
  print "  * %s %s"%(key,value)
print ""
print "############################################################"
print ""

#Output information to file

# INPUT FROM GROUP A!!!
N_data = 5
u_data = 1.2
N_bkg = 3.5
u_bkg = 0.9

f_output = open( "systematics_%s.txt"%channelString, "w" )
f_output.write( "channel %s\n"%channelString )
f_output.write( "data %5.4f %5.4f\n"%(N_data,u_data)) ### UPDATE THIS WHEN NUMBERS ARE AVAILABLE
f_output.write( "bkg %5.4f %5.4f\n"%(N_bkg,u_bkg)) ### UPDATE THIS WHEN NUMBERS ARE AVAILABLE
f_output.write( "sig %5.4f %5.4f\n"%(n_nominal*intLumi,n_nominal_statUnc))
f_output.write( "-----------------\n")
for syst, unc in syst_values.iteritems():
  f_output.write( "%s %5.4f\n"%(syst,unc) )
f_output.write( "-----------------\n")
f_output.close()
