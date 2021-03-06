import socket


  
MYCC='gcc'
MYCXX='g++'
OPT = "-O3"
STD_FLAGS      = OPT + " -g --std=c++0x -Wall " 
STD_LINK_FLAGS = OPT + " -g --std=c++0x "


	       
BOOST_THREAD = "boost_thread" if socket.gethostname() == "moose" else "boost_thread-mt"

CPLEX_FLAGS   = OPT + " -g --std=c++0x -pipe -fexceptions -DNDEBUG -DIL_STD  -mfpmath=sse "
CPLEX_HOME    = "/opt/ibm/ILOG/CPLEX_Studio124/"
CPLEX_INCLUDE = [ CPLEX_HOME + 'cplex/include',
                  CPLEX_HOME + 'concert/include' ]                  
BELOV_INCLUDE = CPLEX_INCLUDE + ['./src/pack/belov_inc']               
CPLEX_LIBPATH = [ CPLEX_HOME + 'cplex/lib/x86-64_sles10_4.1/static_pic',
                 CPLEX_HOME + 'concert/lib/x86-64_sles10_4.1/static_pic']
CPLEX_LIBS    = ['ilocplex', 'cplex', 'concert', 'm', 'pthread',  BOOST_THREAD, 'boost_program_options']

env = Environment(CC=MYCC, CXX=MYCXX,
                  CCFLAGS=STD_FLAGS,  
                  #LIBS=['boost_program_options', 'boost_system','boost_filesystem'],
		  LIBS=['boost_program_options'],	
                  LINKFLAGS=STD_LINK_FLAGS,                      
                  LIBPATH = ['/usr/lib64/'])

boost_env = Environment(CC=MYCC, CXX=MYCXX,
                          CCFLAGS=STD_FLAGS,  
                          LIBS=['boost_program_options',  BOOST_THREAD],
                          LINKFLAGS=STD_LINK_FLAGS,                      
                          LIBPATH = ['/usr/lib64/'])
              
env_belov = Environment(CC=MYCC, CXX=MYCXX,
                        CCFLAGS=CPLEX_FLAGS,
                        CPPPATH=BELOV_INCLUDE,                        
                        LIBS=CPLEX_LIBS,                        
                        LIBPATH=CPLEX_LIBPATH,
                        LINKFLAGS=STD_LINK_FLAGS)
                        

vd        			= "./bin/"
vd_pack       	= vd + "pack/"
vd_utils        = vd + "utils/"
vd_ss           = vd + "ss/"
vd_ss_extended  = vd_ss + "extended/"
vd_completion 	= vd_pack + "completions/"

vd_belov      	= vd_pack + "belov/"
vd_test   			= vd + "test/"
vd_partition 		= vd + "partition/" 
vd_main         = vd + "main/"

VariantDir(vd,"./src",duplicate=0)
VariantDir(vd_belov,"./pack/belov/cpp",duplicate=0)


###################################################################
# For Cached Iterative Weakening, both normal and low cardinality #
###################################################################
ciw_source = [  					    		       
  vd             + "Utils.cpp",
  vd_main 			 + "CLI.cpp",
	vd_main        + "MainUtils.cpp",	
	vd_pack        + "PackingUtils.cpp", 	   
  vd_partition   + "CachedIETrees.cpp",
  vd_partition   + "CachedSetsBuffered.cpp",
  vd_partition   + "GenerateCachedIECardinality.cpp",
  vd_partition   + "GenerateIEBitset.cpp",
  vd_partition   + "Partition.cpp",
  vd_partition   + "PartitionIterativeWeakening.cpp",  
  vd_partition   + "PartitionIterativeWeakeningLowCardinality.cpp",  
  vd_partition   + "PartitionUtils.cpp",
  vd_ss_extended + "HorowitzSahniLowCardinalityMSubsets.cpp", 
  vd_ss_extended + "Schroeppel_Shamir_M_Subsets.cpp",  
  vd_ss          + "KK.cpp",      
  vd_ss          + "SchroeppelShamir.cpp",		   
  vd_ss          + "SubsetSum.cpp",
  vd_utils       + "MemoryUsage.cpp", 
]

###############
# For Moffitt #
###############
moffitt_source = [  					    		       
  vd             + "Utils.cpp",
	vd_main 			 + "CLI.cpp",
	vd_main        + "MainUtils.cpp",	
	vd_pack        + "PackingUtils.cpp", 	     
  vd_partition   + "Moffitt.cpp",  
  vd_partition   + "PartitionUtils.cpp",
  vd_ss          + "KK.cpp",
]

###############
# For CGA_MW #
###############
cga_mw_source = [  					    		       
  vd             + "Utils.cpp",
	vd_main 	 + "CLI.cpp",
	vd_main      + "MainUtils.cpp",	
	vd_pack      + "PackingUtils.cpp", 	     
  vd_partition   + "CGA_MW.cpp",  
  vd_partition   + "PartitionUtils.cpp",  
]

###########
# For RNP #
###########
rnp_source = [                           
  vd             + "Utils.cpp",
  vd_main        + "CLI.cpp",
  vd_main        + "MainUtils.cpp", 
  vd_pack        + "PackingUtils.cpp",       
  vd_partition   + "RNP.cpp",  
  vd_partition   + "PartitionUtils.cpp",
  vd_ss          + "KK.cpp",
]

###############
# For RNP2009 #
###############
rnp2009_source = [                           
  vd             + "Utils.cpp",
  vd_main        + "CLI.cpp",
  vd_main        + "MainUtils.cpp", 
  vd_pack        + "PackingUtils.cpp",       
  vd_partition   + "RNP_2009.cpp",  
  vd_partition   + "PartitionUtils.cpp",
  vd_ss          + "KK.cpp",
]

###########
# For SNP #
###########
snp_source = [  					    		       
  vd             + "Utils.cpp",
	vd_main 			 + "CLI.cpp",
	vd_main        + "MainUtils.cpp",	
	vd_pack        + "PackingUtils.cpp", 	     
  vd_partition   + "SNP.cpp",  
  vd_partition   + "PartitionUtils.cpp",
  vd_ss          + "KK.cpp",
]

############
# For BSBC #
############
bsbc_source = [  					    		       
  vd             + "Utils.cpp",
  vd_completion  + "VectorCompletionGenerator.cpp",
  vd_completion  + "IECompletionGenerator.cpp",
  vd_completion  + "BufferedCompletionGenerator.cpp",
  vd_completion  + "OPSCompletionGenerator.cpp",
	vd_completion  + "SSCompletionGenerator.cpp",
	vd_main 			 + "CLI.cpp",
	vd_main        + "MainUtils.cpp",	
	vd_pack        + "BinCompletion.cpp",
	vd_pack        + "BinCompletionUtils.cpp",
	vd_pack        + "PackingUtils.cpp", 	       
  vd_partition   + "BinarySearch.cpp",
  vd_partition   + "BinarySearchBC.cpp",  
  vd_partition   + "PartitionUtils.cpp",
  vd_ss          + "KK.cpp",  
	vd_ss          + "SchroeppelShamir.cpp",
	vd_ss          + "SubsetSum.cpp",
	vd_ss_extended + "Extended_Schroeppel_Shamir.cpp",
	vd_utils 		   + "MemoryUsage.cpp",	  
]

#############
# For BSBCP #
#############
belov_source = [
  vd_belov + "bb.cpp",
  vd_belov + "bb_mos.cpp",
  vd_belov + "cutgen.cpp",
  vd_belov + "main.cpp",
  vd_belov + "probl_csp1.cpp",
  vd_belov + "probl_cp22.cpp",
  vd_belov + "probl_pmp1.cpp",
  vd_belov + "timer.cpp",
  vd_belov + "bbcuts.cpp",
  vd_belov + "mydebug.cpp",
  vd_belov + "stdafx.cpp",
  vd_belov + "bcp.cpp",
  vd_belov + "bcp_branch.cpp",
  vd_belov + "bcp_cuts.cpp",
  vd_belov + "bcp_lp.cpp",
  vd_belov + "bcpstuff.cpp",
  vd_belov + "bcp2.cpp",
  vd_belov + "wrpcpx.cpp",
  vd_belov + "mytools.cpp",
  vd_belov + "subgr.cpp",
  vd_belov + "raster.cpp",
  vd_belov + "spprc_dp.cpp",
  vd_belov + "spprc_dp1.cpp",
  vd_belov + "bkp_dp.cpp",   
  vd_belov + "ethan.cpp",
  ]

bsbcp_source = [  					    		       
  env.Object(vd             + "Utils.cpp"),
  env.Object(vd_completion  + "VectorCompletionGenerator.cpp"),
  env.Object(vd_completion  + "IECompletionGenerator.cpp"),
  env.Object(vd_completion  + "BufferedCompletionGenerator.cpp"),
  env.Object(vd_completion  + "OPSCompletionGenerator.cpp"),
	env.Object(vd_completion  + "SSCompletionGenerator.cpp"),
	env.Object(vd_main 			 + "CLI.cpp"),
	env.Object(vd_main        + "MainUtils.cpp"),	
	env.Object(vd_pack        + "BinCompletion.cpp"),
	env.Object(vd_pack        + "BinCompletionUtils.cpp"),
	env.Object(vd_pack        + "PackingUtils.cpp"), 	       
  env.Object(vd_partition   + "BinarySearch.cpp"),
  env.Object(vd_partition   + "BinarySearchBCP.cpp"),  
  env.Object(vd_partition   + "PartitionUtils.cpp"),
  env.Object(vd_ss          + "KK.cpp"),  
	env.Object(vd_ss          + "SchroeppelShamir.cpp"),
	env.Object(vd_ss          + "SubsetSum.cpp"),
	env.Object(vd_ss_extended + "Extended_Schroeppel_Shamir.cpp"),
	env.Object(vd_utils 		   + "MemoryUsage.cpp"),	  
] + belov_source

############################
# For two-way partitioning #
############################
partition_source = [
 	vd + "PartitioningMain.cpp",  					    
 	env.Object(vd_partition + "PartitionUtils.cpp"),		   
 	vd_partition + "Partition.cpp",
	vd_partition + "GenerateCachedIECardinality.cpp",
	vd_partition + "CachedIETrees.cpp",
  vd_partition + "CachedSetsBuffered.cpp",			
  env.Object(vd_ss_extended + "Horowitz_Sahni_Cardinality.cpp"),
  env.Object(vd_ss_extended + "HorowitzSahniLowCardinalityMSubsets.cpp"),
	env.Object(vd_ss_extended + "Schroeppel_Shamir_M_Subsets.cpp"),
 	env.Object(vd_ss + "CGA.cpp"),
	env.Object(vd_ss + "CKK.cpp"),
	env.Object(vd_ss + "DPPartition.cpp"),
 	env.Object(vd_ss + "HorowitzSahni.cpp"),
 	env.Object(vd_ss + "SchroeppelShamir.cpp"),
 	env.Object(vd_utils + "MemoryUsage.cpp"),
	env.Object(vd_ss + "SubsetSum.cpp"),
	env.Object(vd_ss + "KK.cpp"),
	env.Object(vd + "Utils.cpp"),	
	env.Object(vd_main + "MainUtils.cpp"),
	env.Object(vd_pack + "PackingUtils.cpp"),		
]    

env.Program(vd + 'ciw'    , [  vd_main + "MainCIW.cpp"] + ciw_source)
env.Program(vd + 'mof'    , [  vd_main + "MainMoffitt.cpp"] + moffitt_source)
env.Program(vd + 'cga_mw' , [  vd_main + "MainCGAMW.cpp"] + cga_mw_source)

env.Program(vd + 'rnp'    , [  vd_main + "MainRNP.cpp"] + rnp_source)
env.Program(vd + 'rnp2009', [  vd_main + "MainRNP2009.cpp"] + rnp2009_source)
env.Program(vd + 'snp'    , [  vd_main + "MainSNP.cpp"] + snp_source)
env.Program(vd + 'bsbc'   , [  vd_main + "MainBSBC.cpp"] + bsbc_source)
env.Program(vd + 'Partitioning'             , partition_source)
env_belov.Program(vd + 'bsbcp', [  vd_main + "MainBSBCP.cpp"] + bsbcp_source)




#env.Program(vd + 'CreatePackingProblems'          , [vd_utils + "CreatePackingProblems.cpp"         , vd_pack + "PackingUtils.cpp", vd + "Utils.cpp"])
#env.Program(vd + 'CreatePartitioningProblems'     , [vd_utils + "CreatePartitioningProblems.cpp"    , vd_pack + "PackingUtils.cpp", vd + "Utils.cpp"])
#env.Program(vd + 'ProcessPackingExperiments'      , [vd_utils + "ProcessPackingExperiments.cpp"     , vd_pack + "PackingUtils.cpp", vd + "Utils.cpp"])
#env.Program(vd + 'ProcessSolution'                , [vd_utils + "ProcessSolution.cpp"])





#################################### old STUFF ####################################

subset_sum_source = [  
  env.Object(vd_ss + "SubsetSum.cpp"),
  env.Object(vd_ss + "SchroeppelShamir.cpp"),	
	env.Object(vd_utils + "MemoryUsage.cpp"),
  env.Object(vd_ss + "CGA.cpp"),
  env.Object(vd_ss + "CKK.cpp"),
  env.Object(vd_ss + "KK.cpp"),  
  env.Object(vd_ss + "OrderedPowerSet.cpp"),
  env.Object(vd_ss_extended + "Extended_Schroeppel_Shamir.cpp"),
  env.Object(vd_ss_extended + "Schroeppel_Shamir_M_Subsets.cpp"),
  env.Object(vd_ss_extended + "Extended_Horowitz_Sahni.cpp"),
  env.Object(vd_ss_extended + "Horowitz_Sahni_Cardinality.cpp"),
  env.Object(vd_ss_extended + "HorowitzSahniLowCardinalityMSubsets.cpp"),
  env.Object(vd_ss_extended + "InclusionExclusion.cpp"),
  env.Object(vd_pack + "BinCompletionUtils.cpp"),
  env.Object(vd_completion + "SSCompletionGenerator.cpp"),  
  env.Object(vd_completion + "VectorCompletionGenerator.cpp"),
  env.Object(vd_completion + "IECompletionGenerator.cpp"),
  env.Object(vd_completion + "BufferedCompletionGenerator.cpp"),
  env.Object(vd_completion + "OPSCompletionGenerator.cpp")
  
]

bin_packing_source = [
  env.Object(vd_pack + "BinCompletion.cpp"),
  env.Object(vd_pack + "BinPackingTest.cpp"),
  env.Object(vd_partition + "RNP.cpp"),   
  #env.Object(vd_ss + "CGA.cpp"),     
  env.Object(vd_pack + "PackingUtils.cpp"), 	
  env.Object(vd + "Utils.cpp"), 	   
]   + subset_sum_source	# + belov_source  # For belov


#multiway_source = [
#  #vd_pack + "PackingUtils.cpp",					   
#  env.Object(vd_partition + "PartitionUtils.cpp"),		   
#  vd_partition + "Moffitt.cpp",
#  vd_partition + "SNP.cpp",
#  vd_partition + "CachedSetsBuffered.cpp",
#  vd_partition + "Partition.cpp",
#  vd_partition + "PartitionIterativeWeakening.cpp",
#  vd_partition + "PartitionIterativeWeakeningLowCardinality.cpp",
#  vd_partition + "CachedIETrees.cpp",
#  vd_partition + "GenerateIEBitset.cpp",
#  vd_partition + "GenerateCachedIECardinality.cpp",
##	vd_partition + "BinarySearch.cpp",
#	env.Object(vd_main + "MainUtils.cpp"),
#]  + bin_packing_source  




extended_partition_source = [
  vd + "ExtendedPartitioningMain.cpp",
  env.Object(vd_ss + "KK.cpp"),  					    
  vd_partition + "PartitionUtils.cpp",		   
  vd_partition + "Partition.cpp",
	env.Object(vd + "Utils.cpp"),
	vd_partition + "GenerateCachedIECardinality.cpp",
	vd_partition + "CachedIETrees.cpp",
  vd_partition + "CachedSetsBuffered.cpp",			
	env.Object(vd_main + "MainUtils.cpp"),
	env.Object(vd_pack + "PackingUtils.cpp"),
	env.Object(vd_ss + "SubsetSum.cpp"),
  env.Object(vd_ss + "SchroeppelShamir.cpp"),
  env.Object(vd_utils + "MemoryUsage.cpp"),	
  env.Object(vd_ss_extended + "Extended_Schroeppel_Shamir.cpp"),
	env.Object(vd_ss_extended + "InclusionExclusion.cpp"),	
]    


#env.Program(vd + 'SubsetSum'                      , [vd + "SubsetSumMain.cpp", vd_pack + "PackingUtils.cpp", vd + "Utils.cpp"] + subset_sum_source)
#env_belov.Program(vd + 'SubsetSelector'            , [vd_partition + "SubsetSelector.cpp"] +  multiway_source)
#env.Program(vd + 'TestDistribution'               , [vd + "TestDistribution.cpp", vd_pack + "PackingUtils.cpp", vd + "Utils.cpp"] + subset_sum_source)
#env_belov.Program(vd + 'kPartitioning'            , [  vd + "kPartitioningMain.cpp"] + multiway_source)

#env_belov.Program(vd + 'ExtendedPartitioning'     , extended_partition_source)

#env.Program(vd + 'BinPacking'               , [vd + "BinPackingMain.cpp"] + bin_packing_source)

