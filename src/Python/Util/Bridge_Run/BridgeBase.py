#!/usr/bin/env python3

import sys, os
from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################

# check 


class BridgeBase:
    def __init__(self, bridge):
        self.module, self.args, self.io, self.progress = bridge.io.module, bridge.args, bridge.io, bridge.io.progress 
        self.pop, self.cores = self.io.pop, str(self.args.cores) 
        self.bd, self.ss, self.pt = self.pop.bdata, self.pop.sumstats, self.pop.genopheno
        self.P = self.io.pop.gen 
        self.model, self.base_paths = {}, dd(bool) 
        if self.args.platform == 'mac':     self.make_job = self.make_mac_job
        else:                               self.make_job = self.make_linux_job         
       


        ### GOTTA ADD EVERYTHINGTO USER SETTINGS --- GET RID OF MOST OF THE X-FIELDS ###
        ### DO AN EXAMPLE WITH --thinned snp list
        ### DO AN EXAMPLE WITH COVARIATES 
        
        ### MAKE DIFFS JUST MINOR OPTIONS IN THE JOBS DIFFS HERE CAN JUST BE MINOR OPTIONS ####
        ### MAYBE ADD USER SETTING ### combine mark 

        if self.module.split('-')[0] != 'prs': 
            self.clump_args = ['--clump-p1','1e-1','--clump-p2','1e-1','--clump-kb','1000','--clump-r2','0.01','--maf','0.001'] 
            self.lambda_val = '0.1,0.2,0.5,1,2,5,10,20'
        else:                                  
            self.clump_args = ['--clump-p1','1e-2','--clump-p2','1e-1','--clump-kb','1000','--clump-r2','0.01','--maf','0.001'] 
            if self.pop.clump_value is not None: 
                try: cval = float(self.pop.clump_value)
                except: cval = 0.1 
                if cval > 0.5: cval = 0.5 
                cval2 = cval*2 
                if cval2 > 0.6: cval2 = 0.6
                self.clump_args = ['--clump-p1',str(cval),'--clump-p2',str(cval2),'--clump-kb','1000','--clump-r2','0.01','--maf','0.001'] 
            
            self.lambda_val = '0.2,0.5,1,2,5,10,20,50'   
            if self.module.split('-')[-1] != 'single': 
                fh = open(self.io.pop.gen['model']) 
                for line in fh: 
                    line = line.strip().split('=')
                    if len(line) != 2: continue 
                    elif line[0] == 'SUMSTATS_SIZE': self.model['pop_size'] = str(int(line[1])) 
                    elif len(line[0].split('_')) == 3 and line[0].split('_')[-1] == 'PREFIX': self.model[line[0].split('_')[1].lower()] = line[1]
                    else: continue 
                fh.close()  
                try:        self.model['pop_size'] = str(int(self.model['pop_size'])) 
                except:     self.model['pop_size'] = '50000'    
                if self.module.split('-')[-1] == 'prior': self.run_beta  = self.prior_beta
                 
            



    
    def make_linux_job(self, run_job, SILENT = False):
        from subprocess import call 
        if not SILENT: self.progress.start_rJob(run_job, self.name) 
        #self.progress.start_rJob(run_job, self.name) mark  
        #if self.name != 'clump': self.progress.start_rJob(run_job) 
        #self.progress.marq(1)  
        if not self.base_paths[self.name]:
            for k, paths in self.base_paths.items(): 
                if paths: 
                    paths[0].close()
                    paths[1].close() 
                    self.base_paths[k] = False 
            self.base_paths[self.name] = [open(self.io.paths[self.name]+'/'+self.name+'.stdout','w'), open(self.io.paths[self.name]+'/'+self.name+'.stderr','w')] 
        call(run_job,stdout=self.base_paths[self.name][0], stderr=self.base_paths[self.name][1]) 
        return        


    def make_mac_job(self, run_job, SILENT = False): 
        #self.progress.marq(1)  
        if not self.base_paths[self.name]:
            for k, paths in self.base_paths.items(): 
                if paths: self.base_paths[k] = False 
            self.base_paths[self.name] = [self.io.paths[self.name]+'/'+self.name+'.stdout',self.io.paths[self.name]+'/'+self.name+'.stderr'] 
        if not SILENT: self.progress.start_rJob(run_job, self.name) 
        if self.name == 'clump': my_job = " ".join(run_job + ['>>',self.base_paths[self.name][0], '2>>', self.base_paths[self.name][1]]) 
        else:                    my_job = " ".join(run_job + ['>',self.base_paths[self.name][0], '2>', self.base_paths[self.name][1]]) 
        os.system(my_job) 
        return  



    def close_all(self): 
        if self.args.platform != 'mac': 
            for k, paths in self.base_paths.items(): 
                if paths: 
                    paths[0].close()
                    paths[1].close() 
                    self.base_paths[k] = False 
        return 



    def run_clump(self,chromosome, name = 'clump'):
        self.name, pp =  name, self.io.programs['plink'] 
        if self.ss.TESTS['NOSNPS']: X = ['--bfile',self.bd.map[chromosome],'--keep', self.bd.id_file] + self.bd.X_fields + self.clump_args
        else:                       X = ['--bfile',self.bd.map[chromosome],'--extract', self.ss.snp_file, '--keep', self.bd.id_file] + self.bd.X_fields + self.clump_args
        rJOB = [pp,'--clump',self.ss.map[chromosome],'--out',self.io.paths['clump']+'/'+self.pop.name+'_clump_'+str(chromosome)] + X 
        self.make_job(rJOB, SILENT = (chromosome != '1')) 
        fname = self.io.paths['clump']+'/'+self.pop.name+'_clump_'+str(chromosome)+'.clumped' 
        if os.path.exists(fname): 
            rJOB = ['gzip','-f',fname] 
            self.make_job(rJOB, SILENT = True) 

        return 

        

    def run_beta(self, name = 'beta'): 
        self.name, pp = name, self.io.programs['est_beta_bychr'] 
        

        #print(self.ss.X_fields) 
        
        #self.ss.X_fields = ['--sumstats.allele0ID', 'REF', '--sumstats.allele1ID', 'A1', '--sumstats.betaID', 'BETA', '--sumstats.frqID', 'A1_FREQ', '--sumstats.nID', 'OBS_CT', '--sumstats.seID', 'SE', '--sumstats.snpID', 'ID']
        #self.ss.X_fields = ['--sumstats.allele0ID', 'REF', '--sumstats.allele1ID', 'A1', '--sumstats.betaID', 'BETA', '--sumstats.snpID', 'ID']
        
        # REMOVE A1_FREQ, OBS_CT, SE 

        X = ['--bfile',self.bd.prefix,'--ld.ids',self.bd.id_file,'--sumstats',self.ss.prefix, '--clump.stem',self.P['clump']] + self.ss.X_fields 

        X.extend(['--beta.stem',self.io.paths['beta']+'/'+self.pop.name+'_beta','--by.chr.sumstats',self.ss.suffix,'--n.cores',str(self.args.cores)]) 
        X.extend(['--S','0,0.25,0.5,0.75,1','--n.max.locus',self.ss.max_clump_size,'--thinned.snplist',self.ss.thin_snps])
        X.extend(['--by.chr', str(int(self.bd.BYCHR)), '--strand.check', '1'])
        X.extend(['--lambda',self.lambda_val]) 
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.make_job(rJOB) 
        return
        

    def run_predict(self, name = 'predict'): 
        self.name, pp = name, self.io.programs['predict_bychr'] 
        X = ['--bfile',self.pt.genotype_prefix,'--out.file',self.io.paths['predict']+'/'+self.pop.name+'_predict','--n.cores',str(self.args.cores)] + self.pt.X_fields
        if self.args.module == 'prs-port': X.extend(['--beta.stem',self.model['beta'],'--ranking','pv','--by.chr', str(int(self.pt.BYCHR)), '--strand.check', '1'])
        else:                              X.extend(['--beta.stem',self.P['beta'],'--ranking','pv','--by.chr', str(int(self.pt.BYCHR)), '--strand.check', '1'])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.make_job(rJOB) 
        return 
    
    def run_quantify(self, name = 'quantify'): 
        self.name = name 
        X = ['--pop2',self.pop.name,'--outfile',self.io.paths['quantify']+'/'+self.pop.name+'_quantify','--n.cores',self.cores]  + self.pt.X_fields
        if self.args.module == 'prs-port':  X.extend(['--pred1',self.P['predict'],'--models1',self.model['beta']]) 
        else:                               X.extend(['--pred1',self.P['predict'],'--models1',self.P['beta']]) 
        if self.args.module == 'prs-prior': pp = self.io.programs['pred_combine_prior_only'] 
        else:                               pp = self.io.programs['pred_combine_en'] 
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X  
        self.make_job(rJOB) 
        return  

    # PRS PRIOR # 
    def prior_beta(self, name = 'beta'): 
        self.name, pp = name, self.io.programs['est_beta_InformPrior_bychr'] 
        


        #X = ['--bfile',self.bd.prefix,'--ld.ids',self.bd.id_file,'--sumstats',self.ss.prefix, '--clump.stem',self.P['clump'],'--n.cores',str(self.args.cores)] + self.ss.X_fields 
        X = ['--bfile',self.bd.prefix,'--ld.ids',self.bd.id_file,'--sumstats',self.ss.prefix, '--n.cores',str(self.args.cores)] + self.ss.X_fields 
        X.extend(['--beta.stem',self.io.paths['beta']+'/'+self.pop.name+'_beta','--by.chr.sumstats',self.ss.suffix,'--ranking','f.stat']) 
        X.extend(['--param.file',self.model['predict']+'_best_model_params.dat','--prior',self.model['prior'],'--fst',str(self.args.fst)])
        
        X.extend(['--sumstats.P',self.ss.fields['P']]) 
        X.extend(['--by.chr', str(int(self.bd.BYCHR)), '--strand.check', '1'])
        


        #X.extend(['--N.pop1',str(self.ss.model_size), '--N.pop2', str(self.ss.pop_size)]) 
        
        X.extend(['--N.pop1',str(self.model['pop_size']), '--N.pop2', str(self.ss.pop_size)]) 
        
        # HMMMM ??? # 
        #alt_opts = ['S', 'n.max.locus', 'thinned.snplist']
        #X.extend(['--S','0,0.25,0.5,0.75,1','--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',str(self.args.thinned_snplist)])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        
        #for x in rJOB: 
        #    print(x, type(x)) 
        #
        #print(" ".join(rJOB)) 


        self.make_job(rJOB) 
        return 


    
    def run_prior(self, name = 'prior'): 
        self.name, pp = name, self.io.programs['est_beta_bychr'] # PRECISION = TRUE 
        X = ['--bfile',self.bd.prefix,'--ld.ids',self.bd.id_file,'--sumstats',self.ss.prefix,'--clump.stem',self.P['clump'],'--beta.stem',self.io.paths['prior']+'/'+self.pop.name+'_prior','--n.cores',str(self.args.cores)] 
        X += self.ss.X_fields
        opt_params = self.P['predict']+'_best_model_params.dat'
        X.extend(['--precision','TRUE','--param.file',opt_params,'--by.chr.sumstats',self.ss.suffix,'--S','1','--lambda','1']) 
        #X.extend(['--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',str(self.args.thinned_snplist)])
        X.extend(['--n.max.locus',str(self.ss.max_clump_size),'--thinned.snplist',self.ss.thin_snps]) 
        X.extend(['--by.chr', str(int(self.bd.BYCHR)), '--strand.check', '1'])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.make_job(rJOB) 
        return 
        
       
        # make prior (for build) prior requires clump and predict 
        # priors beta requires just clump and model stuff !!! 
        #  predict requires beta [ but not for port - that only needs the model ] 
        # quantify requires -> [predict and beta ] and for port it only needs predict 

        # frqid 














