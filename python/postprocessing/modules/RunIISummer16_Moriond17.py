def getXsec(sample):
	# ttbar
	if sample.find("TTJets_DiLept"                                ) !=-1 : return 831.76*((3*0.108)*(3*0.108));
	elif sample.find("TTJets_SingleLeptFromT"                     ) !=-1 : return 831.76*(3*0.108)*(1-3*0.108);
	elif sample.find("TTJets_SingleLeptFromTbar"                  ) !=-1 : return 831.76*(3*0.108)*(1-3*0.108);
	elif sample.find("TTJets_HT-600to800"                         ) !=-1 : return 1.610*831.76/502.2;
	elif sample.find("TTJets_HT-800to1200"                        ) !=-1 : return 0.663*831.76/502.2;
	elif sample.find("TTJets_HT-1200to2500"                       ) !=-1 : return 0.12*831.76/502.2;
	elif sample.find("TTJets_HT-2500toInf"                        ) !=-1 : return 0.001430*831.76/502.2;
	#W+Jets                                                       
	elif sample.find("WJetsToLNu_HT-1200To2500"                   ) !=-1 : return 1.329*1.21;
	elif sample.find("WJetsToLNu_HT-2500ToInf"                    ) !=-1 : return 0.03216*1.21;
	elif sample.find("WJetsToLNu_HT-400To600"                     ) !=-1 : return 48.91*1.21;
	elif sample.find("WJetsToLNu_HT-600To800"                     ) !=-1 : return 12.05;
	elif sample.find("WJetsToLNu_HT-800To1200"                    ) !=-1 : return 5.501;
	#QCD                                                          
	elif sample.find("QCD_HT500to700"                             ) !=-1 : return 31630.0;
	elif sample.find("QCD_HT700to1000"                            ) !=-1 : return 6802.0;
	elif sample.find("QCD_HT1000to1500"                           ) !=-1 : return 1206.0;
	elif sample.find("QCD_HT1500to2000"                           ) !=-1 : return 120.4;
	elif sample.find("QCD_HT2000toInf"                            ) !=-1 : return 25.25;
	#Single Top                                                   
	elif sample.find("TTTo2L2Nu_"                                 ) !=-1 : return 831.76*((3*0.108)**2);
	elif sample.find("ST_s-channel_4f_leptonDecays"               ) !=-1 : return 3.365;
	elif sample.find("ST_t-channel_antitop_4f_inclusiveDecays"    ) !=-1 : return 80.95;
	elif sample.find("ST_t-channel_top_4f_inclusiveDecays"        ) !=-1 : return 80.95;
	elif sample.find("ST_tW_antitop_5f_NoFullyHadronicDecays"     ) !=-1 : return 19.55;
	elif sample.find("ST_tW_top_5f_NoFullyHadronicDecays"         ) !=-1 : return 19.55;
	#Drell-Yan
	elif sample.find("DYJetsToLL_M-50_HT-100to200_"               ) !=-1 : return 147.4*1.23;
	elif sample.find("DYJetsToLL_M-50_HT-1200to2500_"             ) !=-1 : return 0.1514*1.23;
	elif sample.find("DYJetsToLL_M-50_HT-200to400_"	              ) !=-1 : return 40.99*1.23;
	elif sample.find("DYJetsToLL_M-50_HT-2500toInf_"              ) !=-1 : return 0.003565*1.23;
	elif sample.find("DYJetsToLL_M-50_HT-400to600_"               ) !=-1 : return 5.678*1.23;
	elif sample.find("DYJetsToLL_M-50_HT-600to800_"               ) !=-1 : return 1.367*1.23;
	elif sample.find("DYJetsToLL_M-50_HT-800to1200_"              ) !=-1 : return 0.6304*1.23;
	#TTV
	elif sample.find("TTWJetsToLNu_"                              ) !=-1 : return 0.2043;
	elif sample.find("TTWJetsToQQ_"                               ) !=-1 : return 0.40620;
	elif sample.find("TTZToLLNuNu_M-10_"                          ) !=-1 : return 0.2529;
	elif sample.find("TTZToQQ_"                                   ) !=-1 : return 0.5297;
	#DiBoson
	elif sample.find("WWTo2L2Nu_"                                 ) !=-1 : return 10.481;
	elif sample.find("WWToLNuQQ_"                                 ) !=-1 : return 43.53;
	elif sample.find("WZTo1L1Nu2Q_"                               ) !=-1 : return 10.71;
	elif sample.find("WZTo2L2Q_"                                  ) !=-1 : return 5.60;
	elif sample.find("ZZTo2L2Nu_"                                 ) !=-1 : return 0.564;
	elif sample.find("ZZTo2L2Q_ "                                 ) !=-1 : return 3.28;
	#Signal Samples 
	elif sample.find( "SMS-T1tttt_mGluino-1200_mLSP-800"            ) !=-1 : return 0.04129;
	elif sample.find( "SMS-T1tttt_mGluino-1500_mLSP-100"            ) !=-1 : return 0.006889;
	elif sample.find( "SMS-T1tttt_mGluino-2000_mLSP-100"            ) !=-1 : return 0.0004488;    
	# put 1 for data 
	elif sample.find("SingleMuon")!=-1  or sample.find("SingleElectron") !=-1 or sample.find("JetHT") !=-1 or sample.find("MET") !=-1 or sample.find("MuonEG") !=-1 or sample.find("DoubleMuon") !=-1 or sample.find("DoubleEG") !=-1: return 1.
	else:
		print "Cross section not defined for this sample!"
		return 0.0


