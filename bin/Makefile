# This submit a job in bigred2

JOBSCR = iu_cluster_job_script.pbs
TESTDiskSCR = test_disk.pbs
TESTSphereSCR = test_sphere.pbs

submit: 
	@echo "Launching the job";
	qsub $(JOBSCR)

test_sphere: 
	@echo "Launching the test job for Sphere shaped NP";
	qsub $(TESTSphereSCR)

test_disk: 
	@echo "Launching the test job for disk shaped NP";
	qsub $(TESTDiskSCR)	
