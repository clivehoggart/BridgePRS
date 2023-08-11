#!/usr/bin/env python3

import os
import sys
from collections import defaultdict as dd
##########################################################################################################################################
##########################################################################################################################################


class BridgeBase:
    def __init__(self, bridge):
        self.module, self.args, self.io, self.settings, self.progress = bridge.io.module, bridge.args, bridge.io, bridge.io.settings, bridge.io.progress 
        self.F, self.P, self.M = self.settings.files, self.settings.prefixes, self.settings.maps 
        
        self.model, self.base_paths = {}, dd(bool) 
        
        
        if self.args.platform == 'mac':     self.make_job = self.make_mac_job
        else:                               self.make_job = self.make_linux_job         

        if 'sumstats' in self.settings.fields: self.prep_sumstats() 

        ### MAYBE ADD USER SETTING ### combine  
        if self.module.split('-')[0] != 'prs': 
            self.clump_args = ['--clump-p1','1e-1','--clump-p2','1e-1','--clump-kb','1000','--clump-r2','0.01'] 
            self.lambda_val = '0.1,0.2,0.5,1,2,5,10,20'
        else:                                  
            self.clump_args = ['--clump-p1','1e-2','--clump-p2','1e-1','--clump-kb','1000','--clump-r2','0.01'] 
            self.lambda_val = '0.2,0.5,1,2,5,10,20,50' 
            #--clump-p1 1e-2 --clump-p2 1e-1 --clump-kb 1000 --clump-r2 0.01 \
            #--clump-p1 1e-1 --clump-p2 1e-1 --clump-kb 1000 --clump-r2 0.01 \


        if self.module.split('-')[0] == 'prs' and self.module.split('-')[-1] != 'single': 
            fh = open(self.F['model']) 
            for line in fh: 
                line = line.strip().split('=')
                if len(line[0].split('_')) != 3: continue 
                x1,x2,x3 = line[0].split('_') 
                y = line[1] 
                if x3 != 'PREFIX': continue 
                self.model[x2.lower()] = y 
            fh.close()  
            if self.module.split('-')[-1] == 'port': self.run_predict, self.run_quantify = self.port_predict, self.port_quantify 
            else:                                    self.run_eval = self.prior_eval 



    def prep_sumstats(self): 
        SF = self.settings.fields['sumstats'] 
        self.c_stat = ['--clump-field',SF['P'],'--clump-snp-field',SF['SNPID']]
        self.e_stat = ['--sumstats.allele0ID',SF['REF'],'--sumstats.allele1ID',SF['ALT'],'--sumstats.betaID',SF['BETA'],'--sumstats.frqID',SF['MAF'],'--sumstats.nID',SF['SS'],'--sumstats.seID',SF['SE'], '--sumstats.snpID',SF['SNPID']]
        return 


    
    def make_linux_job(self, run_job):
        from subprocess import call 
        if self.name != 'clump': self.progress.start_rJob(run_job) 
        if not self.base_paths[self.name]:
            for k, paths in self.base_paths.items(): 
                if paths: 
                    paths[0].close()
                    paths[1].close() 
                    self.base_paths[k] = False 
            self.base_paths[self.name] = [open(self.io.paths[self.name]+'/'+self.name+'.stdout','w'), open(self.io.paths[self.name]+'/'+self.name+'.stderr','w')] 
        call(run_job,stdout=self.base_paths[self.name][0], stderr=self.base_paths[self.name][1]) 
        return        


    def make_mac_job(self, run_job): 
        if not self.base_paths[self.name]:
            for k, paths in self.base_paths.items(): 
                if paths: self.base_paths[k] = False 
            self.base_paths[self.name] = [self.io.paths[self.name]+'/'+self.name+'.stdout',self.io.paths[self.name]+'/'+self.name+'.stderr'] 
        if self.name == 'clump': 
            my_job = " ".join(run_job + ['>>',self.base_paths[self.name][0], '2>>', self.base_paths[self.name][1]]) 
        else: 
            self.progress.start_rJob(run_job) 
            my_job = " ".join(run_job + ['>',self.base_paths[self.name][0], '2>', self.base_paths[self.name][1]]) 
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
        X = ['--bfile',self.P['bfile'],'--extract',self.F['snp'],'--keep',self.F['id']] + self.c_stat 
        X.extend(self.clump_args) 
        rJOB = [pp,'--clump',self.M['sumstats'][chromosome],'--out',self.io.paths['clump']+'/'+self.args.popname+'_clump_'+str(chromosome)] + X 
        self.make_job(rJOB) 
        self.progress.mark()  
        rJOB = ['gzip','-f',self.io.paths['clump']+'/'+self.args.popname+'_clump_'+str(chromosome)+'.clumped'] 
        self.make_job(rJOB) 
        return 


    def run_eval(self, name = 'eval'): 
        self.name, pp = name, self.io.programs['est_beta_bychr'] 
        X = ['--bfile',self.P['bfile'],'--ld.ids',self.F['id'],'--sumstats',self.P['sumstats'],'--clump.stem',self.P['clump']]  
        X.extend(['--beta.stem',self.io.paths['eval']+'/'+self.args.popname+'_eval','--by.chr.sumstats',self.args.sumstats_suffix,'--n.cores',str(self.args.cores)]) 
        X.extend(self.e_stat) 
        # HARDCODE # 
        X.extend(['--S','0,0.25,0.5,0.75,1','--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',str(self.args.thinned_snplist)])
        X.extend(['--by.chr', '0', '--strand.check', '1'])
        X.extend(['--lambda',self.lambda_val]) 
        
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.make_job(rJOB) 


        
        #r_log, r_err = open(self.io.paths['eval']+'/eval.stdout','w'), open(self.io.paths['eval']+'/eval.stderr','w')
        #call(rJOB,stdout=r_log, stderr=r_err)
        #r_log.close() 
        #r_err.close() 
        return 
    


    
    def run_predict(self, name = 'predict'): 
        self.name, pp = name, self.io.programs['predict_bychr'] 
        X = ['--bfile',self.P['bfile'],'--beta.stem',self.P['eval'],'--out.file',self.io.paths['predict']+'/'+self.args.popname+'_predict','--n.cores',str(self.args.cores)] 
        X.extend(['--test.data',self.F['pheno'], '--valid.data',self.F['validation']]) 
        
        X.extend(['--ranking','pv','--pheno.name',self.settings.fields['pheno']['NAME'], '--cov.names',self.settings.fields['pheno']['COVARIATES']])
        # HARDCODE # 
        X.extend(['--by.chr', '0', '--strand.check', '1'])
        
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        
        self.make_job(rJOB) 
        
        
        self.progress.start_rJob(rJOB) 
        #r_log, r_err = open(self.io.paths['predict']+'/predict.stdout','w'), open(self.io.paths['predict']+'/predict.stderr','w')
        #call(rJOB,stdout=r_log, stderr=r_err)
        #r_log.close() 
        #r_err.close() 
        return 
    
    
    
    def run_quantify(self, name = 'quantify'): 
        self.name, pp = name, self.io.programs['quantify']  
        X = ['--pop2',self.args.popname,'--outdir',self.io.paths['quantify'],'--n.cores',str(self.args.cores)] 
        X.extend(['--test.data',self.F['pheno'], '--valid.data',self.F['validation']]) 
        X.extend(['--pred1',self.P['predict'],'--eval1',self.P['eval']]) 
        X.extend(['--pheno.name',self.settings.fields['pheno']['NAME'], '--cov.names',self.settings.fields['pheno']['COVARIATES']])
        
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X  
        
        self.make_job(rJOB) 
        #self.progress.start_rJob(rJOB) 
        #r_log, r_err = open(self.io.paths['quantify']+'/quantify.stdout','w'), open(self.io.paths['quantify']+'/quantify.stderr','w')
        #call(rJOB,stdout=r_log, stderr=r_err)
        #r_log.close() 
        #r_err.close() 
        return 




    #### PRS PORT ### 

    def port_predict(self, name='predict'): 
        
        self.name, pp = name, self.io.programs['predict_bychr'] 
        X = ['--bfile',self.P['bfile'],'--beta.stem',self.model['eval'],'--out.file',self.io.paths['predict']+'/'+self.args.popname+'_predict','--n.cores',str(self.args.cores)] 
        X.extend(['--test.data',self.F['pheno'], '--valid.data',self.F['validation']]) 
        
        #run_opts = ['by.chr','strand.check','ranking','pheno.name','cov.names','non.overlapping','p.thresh'] 
        #run_opts = ['by.chr','strand.check','ranking','pheno.name','cov.names','p.thresh'] 
        #X.extend(['--pheno.name',self.settings.fields['pheno']['NAME'], '--cov.names',self.settings.fields['pheno']['COVARIATES']])
        
        X.extend(['--ranking','pv','--pheno.name',self.settings.fields['pheno']['NAME'], '--cov.names',self.settings.fields['pheno']['COVARIATES']])
        # HARDCODE # 
        X.extend(['--by.chr', '0', '--strand.check', '1'])
        
        
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        
        self.make_job(rJOB) 
        #self.progress.start_rJob(rJOB) 
        #r_log, r_err = open(self.io.paths['predict']+'/predict.stdout','w'), open(self.io.paths['predict']+'/predict.stderr','w')
        #call(rJOB,stdout=r_log, stderr=r_err)
        #r_log.close() 
        #r_err.close() 
        return 



    def port_quantify(self, name = 'quantify'): 
        self.name, pp = name, self.io.programs['quantify']  
        #F, P = self.settings.files, self.settings.prefixes 
        X = ['--pop2',self.args.popname,'--outdir',self.io.paths['quantify'],'--n.cores',str(self.args.cores)] 
        X.extend(['--test.data',self.F['pheno'], '--valid.data',self.F['validation']]) 
        X.extend(['--pred1',self.P['predict'],'--eval1',self.model['eval']]) 
        X.extend(['--pheno.name',self.settings.fields['pheno']['NAME'], '--cov.names',self.settings.fields['pheno']['COVARIATES']])
        
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X  
        self.make_job(rJOB) 
        #self.progress.start_rJob(rJOB) 
        #r_log, r_err = open(self.io.paths['quantify']+'/quantify.stdout','w'), open(self.io.paths['quantify']+'/quantify.stderr','w')
        #call(rJOB,stdout=r_log, stderr=r_err)
        #r_log.close() 
        #r_err.close() 
        return 



    #### PRS BRIDGE ### 

    def prior_eval(self, name = 'eval'): 
        
        self.name, pp = name, self.io.programs['est_beta_InformPrior_bychr'] 
        X = ['--bfile',self.P['bfile'],'--ld.ids',self.F['id'],'--sumstats',self.P['sumstats'],'--clump.stem',self.P['clump'], '--n.cores',str(self.args.cores)]  
        
        X.extend(['--beta.stem',self.io.paths['eval']+'/'+self.args.popname+'_eval','--by.chr.sumstats','.Phenotype.glm.linear.gz']) 
        X.extend(['--param.file',self.model['optimize']+'_best_model_params.dat','--prior',self.model['prior']])
        X.extend(['--ranking','f.stat','--fst',str(self.args.fst)])   
        X.extend(self.e_stat) 
        X.extend(['--sumstats.P',self.settings.fields['sumstats']['P']]) 
        # HARDCODE #
        X.extend(['--by.chr', '0', '--strand.check', '1'])
        
        # HMMMM # 
        #alt_opts = ['S', 'n.max.locus', 'thinned.snplist']
        #X.extend(['--S','0,0.25,0.5,0.75,1','--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',str(self.args.thinned_snplist)])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.make_job(rJOB) 
        #self.progress.start_rJob(rJOB) 
        #r_log, r_err = open(self.io.paths['eval']+'/eval.stdout','w'), open(self.io.paths['eval']+'/eval.stderr','w')
        #call(rJOB,stdout=r_log, stderr=r_err)
        #r_log.close() 
        #r_err.close() 
        return 



    def run_optimize(self, name = 'optimize'): 
        
        self.name, pp = name, self.io.programs['predict_bychr'] 
        X = ['--bfile',self.P['bfile'],'--beta.stem',self.P['eval'],'--out.file',self.io.paths['optimize']+'/'+self.args.popname+'_optimize','--n.cores',str(self.args.cores)] 
        X.extend(['--test.data',self.F['pheno'], '--valid.data','0']) 
        X.extend(['--ranking','pv','--pheno.name',self.settings.fields['pheno']['NAME'], '--cov.names',self.settings.fields['pheno']['COVARIATES']])
        
        # HARDCODE # 
        X.extend(['--by.chr', '0', '--strand.check', '1'])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        
        self.make_job(rJOB) 
        
        #self.progress.start_rJob(rJOB) 
        #r_log, r_err = open(self.io.paths['optimize']+'/optimize.stdout','w'), open(self.io.paths['optimize']+'/optimize.stderr','w')
        #call(rJOB,stdout=r_log, stderr=r_err)
        #r_log.close() 
        #r_err.close() 
        return 


    def run_prior(self, name = 'prior'): 
        self.name, pp = name, self.io.programs['est_beta_bychr'] # PRECISION = TRUE 
        F, P, M = self.settings.files, self.settings.prefixes, self.settings.maps 

        X = ['--bfile',self.P['bfile'],'--ld.ids',self.F['id'],'--sumstats',self.P['sumstats'],'--clump.stem',self.P['clump'],'--beta.stem',self.io.paths['prior']+'/'+self.args.popname+'_prior','--n.cores',str(self.args.cores)] 
        opt_params = self.P['optimize']+'_best_model_params.dat' 
        X.extend(['--precision','TRUE','--param.file',opt_params,'--by.chr.sumstats','.Phenotype.glm.linear.gz']) 
        X.extend(['--S','1','--lambda','1','--help','FALSE']) 
        X.extend(self.e_stat) 
        X.extend(['--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',str(self.args.thinned_snplist)])
        
        X.extend(['--by.chr', '0', '--strand.check', '1'])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        
        self.make_job(rJOB) 

        #self.progress.start_rJob(rJOB) 
        #r_log, r_err = open(self.io.paths['prior']+'/prior.stdout','w'), open(self.io.paths['prior']+'/prior.stderr','w')
        #call(rJOB,stdout=r_log, stderr=r_err)
        #r_log.close() 
        #r_err.close() 
        return 
        
        
        
##########################################################################################################################################
