#!/usr/bin/env python3

import os
import sys
from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################


class BridgeBase:
    def __init__(self, bridge):
        self.module, self.args, self.io, self.settings, self.progress = bridge.io.module, bridge.args, bridge.io, bridge.io.settings, bridge.io.progress 
        #self.popname = popname
       

        self.pop = self.settings.pop  
       
        self.bd, self.ss, self.pt = self.pop.bdata, self.pop.sumstats, self.pop.phenotypes 
        
        self.P, self.F = self.settings.prefixes, self.settings.files 


        
        self.model, self.base_paths = {}, dd(bool) 
        if self.args.platform == 'mac':     self.make_job = self.make_mac_job
        else:                               self.make_job = self.make_linux_job         
        
        
        
        ### MAYBE ADD USER SETTING ### combine mark 
        if self.module.split('-')[0] != 'prs': 
            self.clump_args = ['--clump-p1','1e-1','--clump-p2','1e-1','--clump-kb','1000','--clump-r2','0.01','--maf','0.001'] 
            self.lambda_val = '0.1,0.2,0.5,1,2,5,10,20'
        else:                                  
            self.clump_args = ['--clump-p1','1e-2','--clump-p2','1e-1','--clump-kb','1000','--clump-r2','0.01','--maf','0.001'] 
            self.lambda_val = '0.2,0.5,1,2,5,10,20,50' 
            
            if self.module.split('-')[-1] != 'single': 
                fh = open(self.settings.files['model']) 
                for line in fh: 
                    line = line.strip().split('=')
                    if len(line[0].split('_')) != 3: continue 
                    x1,x2,x3 = line[0].split('_') 
                    y = line[1] 
                    if x3 != 'PREFIX': continue 
                    self.model[x2.lower()] = y 
                fh.close()  
                if self.module.split('-')[-1] == 'port': self.run_predict, self.run_quantify   = self.port_predict, self.port_quantify
                else:                                    self.run_beta, self.run_quantify      = self.prior_beta, self.prior_quantify 





    
    def make_linux_job(self, run_job, SILENT = False):
        from subprocess import call 
        
        #print('yooooo', SILENT) 

        
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
        BOOL1 = chromosome!='1'
        self.name, pp =  name, self.io.programs['plink'] 
        X = ['--bfile',self.bd.map[chromosome],'--extract', self.ss.snp_file, '--keep', self.bd.id_file] + self.bd.X_fields + self.clump_args
        rJOB = [pp,'--clump',self.ss.map[chromosome],'--out',self.io.paths['clump']+'/'+self.pop.name+'_clump_'+str(chromosome)] + X 
        self.make_job(rJOB, SILENT = (chromosome != '1')) 
        rJOB = ['gzip','-f',self.io.paths['clump']+'/'+self.pop.name+'_clump_'+str(chromosome)+'.clumped'] 
        self.make_job(rJOB, SILENT = True) 
        return 

        

    def run_beta(self, name = 'beta'): 
        self.name, pp = name, self.io.programs['est_beta_bychr'] 
        X = ['--bfile',self.bd.prefix,'--ld.ids',self.bd.id_file,'--sumstats',self.ss.prefix, '--clump.stem',self.P['clump']] + self.ss.X_fields 
        X.extend(['--beta.stem',self.io.paths['beta']+'/'+self.pop.name+'_beta','--by.chr.sumstats',self.ss.suffix,'--n.cores',str(self.args.cores)]) 
        X.extend(['--S','0,0.25,0.5,0.75,1','--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',self.ss.thin_snps])
        X.extend(['--by.chr', str(int(self.bd.BYCHR)), '--strand.check', '1'])
        
        X.extend(['--lambda',self.lambda_val]) 
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        
        
        

        
        self.make_job(rJOB) 
        return
        


    
    def run_predict(self, name = 'predict'): 
        self.name, pp = name, self.io.programs['predict_bychr'] 
        

        #print('bro bro bro') 
        #print('bro bro bro') 
        #print('bro bro bro')
        #print(self.P['beta']) 
        X = ['--bfile',self.pt.genotype_prefix,'--beta.stem',self.P['beta'],'--out.file',self.io.paths['predict']+'/'+self.pop.name+'_predict','--n.cores',str(self.args.cores)] + self.pt.X_fields
        X.extend(['--ranking','pv','--by.chr', str(int(self.pt.BYCHR)), '--strand.check', '1'])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.make_job(rJOB) 
        return 
    
    
    


    def run_quantify(self, name = 'quantify'): 
        self.name, pp = name, self.io.programs['pred_combine_en']  
        X = ['--pop2',self.pop.name,'--outfile',self.io.paths['quantify']+'/'+self.pop.name+'_quantify','--n.cores',str(self.args.cores)] + self.pt.X_fields  
        X.extend(['--pred1',self.P['predict'],'--models1',self.P['beta']]) 
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X  
        self.make_job(rJOB)     
        return 








    #### PRS PORT ### 
    def port_predict(self, name='predict'): 
        self.name, pp = name, self.io.programs['predict_bychr'] 
        X = ['--bfile',self.pt.genotype_prefix,'--beta.stem',self.model['beta'],'--out.file',self.io.paths['predict']+'/'+self.pop.name+'_predict','--n.cores',str(self.args.cores)]  + self.pt.X_fields  
        X.extend(['--ranking','pv','--by.chr', str(int(self.pt.BYCHR)), '--strand.check', '1'])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.make_job(rJOB) 
        return 
    
    
    def port_quantify(self, name = 'quantify'): 
        self.name, pp = name, self.io.programs['pred_combine_en']  
        X = ['--pop2',self.pop.name,'--outfile',self.io.paths['quantify']+'/'+self.pop.name+'_quantify','--n.cores',str(self.args.cores)] + self.pt.X_fields 
        X.extend(['--pred1',self.P['predict'],'--models1',self.model['beta']]) 
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X  
        self.make_job(rJOB) 
        return  






    # PRS PRIOR # 
    def prior_beta(self, name = 'beta'): 
        
        self.name, pp = name, self.io.programs['est_beta_InformPrior_bychr'] 
        X = ['--bfile',self.bd.prefix,'--ld.ids',self.bd.id_file,'--sumstats',self.ss.prefix, '--clump.stem',self.P['clump'],'--n.cores',str(self.args.cores)] + self.ss.X_fields 
        X.extend(['--beta.stem',self.io.paths['beta']+'/'+self.pop.name+'_beta','--by.chr.sumstats',self.ss.suffix,'--ranking','f.stat']) 
        X.extend(['--param.file',self.model['optimize']+'_best_model_params.dat','--prior',self.model['prior'],'--fst',str(self.args.fst)])
        
        X.extend(['--sumstats.P',self.ss.fields['P']]) 
        X.extend(['--by.chr', str(int(self.bd.BYCHR)), '--strand.check', '1'])
        # HMMMM ??? # 
        #alt_opts = ['S', 'n.max.locus', 'thinned.snplist']
        #X.extend(['--S','0,0.25,0.5,0.75,1','--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',str(self.args.thinned_snplist)])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.make_job(rJOB) 
        return 


    def prior_quantify(self, name = 'quantify'): 
        self.name, pp = name, self.io.programs['pred_combine_prior_only'] 
        X = ['--pop2',self.pop.name,'--outfile',self.io.paths['quantify']+'/'+self.pop.name+'_quantify','--n.cores',str(self.args.cores)] + self.pt.X_fields 
        X.extend(['--pred1',self.P['predict'],'--models1',self.P['beta']]) 
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X  
        self.make_job(rJOB) 
        return  




    # PRS MODEL # 
    def run_optimize(self, name = 'optimize'): 
        self.name, pp = name, self.io.programs['predict_bychr'] 

        X = ['--bfile',self.pt.genotype_prefix,'--beta.stem',self.P['beta'],'--out.file',self.io.paths['optimize']+'/'+self.pop.name+'_optimize','--n.cores',str(self.args.cores)]  + self.pt.X_fields  
        X.extend(['--ranking','pv','--by.chr', str(int(self.pt.BYCHR)), '--strand.check', '1'])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.make_job(rJOB) 
        return 


    def run_prior(self, name = 'prior'): 
        self.name, pp = name, self.io.programs['est_beta_bychr'] # PRECISION = TRUE 
        X = ['--bfile',self.bd.prefix,'--ld.ids',self.bd.id_file,'--sumstats',self.ss.prefix,'--clump.stem',self.P['clump'],'--beta.stem',self.io.paths['prior']+'/'+self.pop.name+'_prior','--n.cores',str(self.args.cores)] 
        X += self.ss.X_fields
        opt_params = self.P['optimize']+'_best_model_params.dat'
        X.extend(['--precision','TRUE','--param.file',opt_params,'--by.chr.sumstats',self.ss.suffix,'--S','1','--lambda','1']) 
        #X.extend(['--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',str(self.args.thinned_snplist)])
        X.extend(['--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',self.ss.thin_snps]) 
        X.extend(['--by.chr', str(int(self.bd.BYCHR)), '--strand.check', '1'])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.make_job(rJOB) 
        return 
        
        
        
##########################################################################################################################################
