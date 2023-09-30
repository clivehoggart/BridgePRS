#!/usr/bin/env python3 

import sys  
#from ..Util.BridgeProgress import BridgeProgress
#from ..Util.BridgeIO import BridgeIO





def bridge_error(eString):
    if type(eString) in [list,tuple]:
        sys.stderr.write('\nBridgeJobError: '+eString[0]+'\n')                     
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeJobError: '+eString+'\n')                           
    sys.exit()





class BridgeJobs:
        def __init__(self,bridge):
            self.args, self.progress = bridge.args, bridge.io.progress 
            


        def run(self, job_target, job_args, job_comment = None, job_marker = 5): 
            self.progress.mark(1) 
            if len(job_args) == 0:                                   
                bridge_error(['No jobs found for '+self.args.module+': '+self.args.cmd+'\n']) 
            elif len(job_args) == 1 or self.args.platform == 'mac':  
                self.progress.mark(5) 
                self.progress.mark(5) 
                self.q_serial(job_target, job_args, job_comment, job_marker) 
            elif self.args.cores > 1: 
                
                self.progress.write('....(Parallelizing Across '+str(self.args.cores)+' Cores)........') 
                
                self.q_parralel(job_target, job_args, job_comment, job_marker) 
            else: 
                self.progress.write('...(Serialized On One Core)') 
                self.q_serial(job_target, job_args, job_comment, job_marker) 
            return 


        def q_serial(self,job_target,job_args,job_comment=None,job_marker=5): 
            
            self.progress.mark(5) 
            for j,(job_key,job_vals) in enumerate(job_args):
                if j % 4 == 0: self.progress.mark(1) 
                job_target(**{a: b for (a,b) in zip(job_key,job_vals)})
            return 


        def q_parralel(self,job_target,job_args,job_comment=None,job_marker=5): 
            import multiprocessing
            for i in range(0,len(job_args),self.args.cores): 
                p_args = [job_arg[1] for job_arg in job_args[i:i+self.args.cores]]
                p_jobs = [multiprocessing.Process(target=job_target,args=p_arg) for p_arg in p_args]
                self.progress.mark(2)  
                for j in p_jobs: j.start() 
                for j in p_jobs: j.join() 
                for j in p_jobs: j.terminate()         
            return


