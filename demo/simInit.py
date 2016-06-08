from mpi4py import MPI
from subprocess import Popen, PIPE, STDOUT, call, check_output, CalledProcessError
import random
import time

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

sim_path = '/sciclone/home00/geogdan/SimTests/demo/pySims_h.R'

iterations = 120

for c in range(0,size):
	if rank == c:
		for i in range((iterations*rank/size + 1), (iterations*(rank+1)/size+1)):
			print "Worker - rank %d on %s."%(rank, name) 
			print "i:%d"%i
			out_path = '/sciclone/home00/geogdan/jTest/sim_'+str(i)+'.csv'
			version = "1"
			""" Init set up for the simulation 
			Parameters
			----------
			nrandom : total numbers of points in the map
			
			xvar_psill : maximum amount of spatial covariation for covariates (spatial auotcorrealtion) 
			
			minx, maxx, miny, maxy : latitude and longitude bound
			
			var1_vrange : distance threshold for covariates
			
			var1_error : pecentage between 0 and 1, total error in the covarite 
			
			prop_acc : (1-prop_acc) is how much error added to the covarite 
			
			var1_error_vrange : the range of error for the covariate
			
			mod_error_magnitude : error coefficient 
			
			trt_prc : the pecentage of treated points
			
			theta : treatment coefficient
			
			beta : covariate coefficient
			 
			spill_vrange : distance of treatment effect   
			
			spill_magnitude : coefficient of spill_vrange
			
			cal : ecnon model 
			
			sample_size : data for prediction
			
			tree_split_lim : the minimum size of node to split
			
			mod_error_vrange : distance of cov error
			
			xvar_error_psill : maximum amount of spatial covariation for covariates for errors
			
			mod_error_psill : maximum amount of saptial cov error 
			
			trt_spill_sill : magnitude of treatment spillover

			tree_thresh: CT drops top and bottom percent of P-scores equal to this value

			thresh_est: For GWR and Spatial PSM, the additive error for the spatial threshold.

			trtcon_overlap: The amount of overlap between treatment and control units in covariate space (0-1)

			
			
			Returns
			-------
			
			"""
			nrandom = "5000"#str(500 + random.random()*4000)
			xvar_psill = "1.0"#str(0.1 + random.random()*0.9)
			minx = "-45.0"
			maxx = "45.0"
			miny = "-22.5"
			maxy = "22.5"

			var1_vrange = "3000"#str(250 + random.random()*2500)
			var1_error = "0.1"#str(0.1 + random.random()*0.9)
			prop_acc = "1.0"#str(0.1 + random.random()*.75)
			var1_error_vrange = "1000"#str(250 + random.random()*2500)
			mod_error_magnitude = "0"#str(0.1 + random.random()*2)
			trt_prc = str(0.05 + random.random() * 0.05)
			theta = "1.0"
			beta = "1.0"#str(0.2 + random.random()*5)
			spill_vrange = "1500"#str(250 + random.random()*2500)
			spill_magnitude= str(0.05 + random.random()*1.0)
			cal= str(0.25 + random.random()*1.75)
			sample_size = "0.5"#str(0.2 + random.random()*0.8)
			tree_split_lim= "5"#str(5 + random.random()*5)
			mod_error_vrange= "1000"#str(250 + random.random()*2500)
			xvar_error_psill = "1"#str(0.1 + random.random()*0.9)
			mod_error_psill = "1"#str(0.1 + random.random()*0.9)
			trt_spill_sill = "1"#str(0.1 + random.random()*0.9)
			tree_thresh = str(0.10 + random.random() * 0.30)
			thresh_est = str(0.05 + random.random() * 0.95)
			trtcon_overlap = str(0.2 + random.random() * 0.75)
			
			outputList = [("version: ",version),("nrandom: ",nrandom),("xvar_psill: ",xvar_psill),("minx: ",minx),
			("maxx: ",maxx),("miny: ", miny),("maxy: ",maxy),("var1_vrange: ",var1_vrange),("var1_error: ",var1_error),
			("prop_acc: ",prop_acc),("var1_error_vrange: ",var1_error_vrange),("mod_error_magnitude: ",mod_error_magnitude),
			("trt_prc: ", trt_prc), ("theta: ", theta),("beta: ", beta),("spill_vrange: ", spill_vrange), 
			("spill_magnitude: ",spill_magnitude),("cal: ",cal), ("sample_size: ",sample_size),
			("tree_split_lim: ",tree_split_lim),("mod_error_vrange: ", mod_error_vrange),
			("xvar_error_psill: ",xvar_error_psill), ("mod_error_psill: ", mod_error_psill),
			("trt_spill_sill: ", trt_spill_sill)]
			
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
				trt_spill_sill+" "+
				tree_thresh+ " "+
				thresh_est+ " "+
				trtcon_overlap,
				stderr=STDOUT,
				shell=True)

				print i
				print R_ret
				print "========================"
			
			except CalledProcessError as sts_err:
				print outputList
				print sts_err
