from mpi4py import MPI
from subprocess import Popen, PIPE, STDOUT, call, check_output, CalledProcessError
import random

comm = MPI.COMM_WORLD

size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

sim_path = '/sciclone/home00/geogdan/SimTests/demo/pySims.R'

iterations = 1000
c = rank
nonfit_cnt = 0
while c < iterations:
    print "Worker - rank %d on %s."%(rank, name) 
    out_path = '/sciclone/home00/geogdan/BetaThree/sim_'+str(c)+'.csv'

    version = "1"
    nrandom = "10000"#str(100 + random.random()*10000)
    xvar_psill = "1.0"#str(0.1 + random.random()*0.9)
    minx = "-45.0"
    maxx = "45.0"
    miny = "-22.5"
    maxy = "22.5"

    var1_vrange = str(0.5 + random.random()*100000)
    var1_error = "0.1"#str(0.1 + random.random()*2)
    prop_acc = str(0.1 + random.random()*.89)
    var1_error_vrange = "0.1"#str(0.1 + random.random()*2)
    mod_error_magnitude = "0"#str(0.1 + random.random()*2)
    trt_prc = str(0.1 + random.random() * 0.3)
    theta = "1.0"
    beta = str(0.2 + random.random()*5)
    spill_vrange = str(0.5 + random.random()*100000)
    spill_magnitude= str(0.01 + random.random()*2)
    cal= str(0.25 + random.random()*1.75)
    sample_size = "0.5"#str(0.1 + random.random()*0.9)
    tree_split_lim= str(0.05 + random.random()*0.1)
    mod_error_vrange= "0.1"#str(0.1 + random.random()*2)
    xvar_error_psill ="0.1" #str(0.1 + random.random()*0.9)
    mod_error_psill ="0.1" #str(0.1 + random.random()*0.9)
    trt_spill_sill = "1.0"#str(0.1 + random.random()*0.9)


    try:
      R_ret = check_output("Rscript " +
                        sim_path + " " +
                        version + " " +
                        nrandom+ " " +
                        xvar_psill+ " " +
                        minx+ " "+
                        maxx+ " "+
                        miny+ " "+
                        maxy+ " "+
                        var1_vrange+ " "+
                        var1_error+ " "+
                        prop_acc+ " "+
                        var1_error_vrange+ " "+
                        mod_error_magnitude+ " "+
                        trt_prc+ " "+
                        theta+ " "+
                        beta+ " "+
                        spill_vrange+ " "+
                        spill_magnitude+ " "+
                        cal+ " "+
                        sample_size+ " "+
                        tree_split_lim+ " "+
                        mod_error_vrange+ " "+
                        out_path+ " " +
                        xvar_error_psill+ " "+
                        mod_error_psill+ " "+
                        trt_spill_sill,
		                  	stderr=STDOUT,
                        shell=True)

      print c
      print R_ret
      print "========================"

    except CalledProcessError as sts_err:                                                                                                   
	    print ">> subprocess error code:", sts_err.returncode, '\n', sts_err.output
	    print "ERROR: Insufficient Matches or another error occured in the R Script. Specifics unknown.  Good luck."
	    print version
	    print nrandom 
	    print xvar_psill 
	    print minx
	    print maxx 
	    print miny
	    print maxy
	    print var1_vrange
	    print var1_error
	    print prop_acc 
	    print var1_error_vrange 
	    print mod_error_magnitude 
	    print trt_prc 
	    print theta 
	    print beta 
	    print spill_vrange 
	    print spill_magnitude
	    print cal
	    print sample_size
	    print tree_split_lim
	    print mod_error_vrange
	    print xvar_error_psill 
	    print mod_error_psill 
	    print trt_spill_sill
      
    c += size

