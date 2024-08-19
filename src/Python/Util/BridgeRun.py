#!/usr/bin/env python3

import os, sys
from collections import defaultdict as dd
from .Bridge_Run.BridgeResult import BridgeResult 
from .Bridge_Run.BridgePlot   import BridgePlot 


#from .BridgeResult import BridgeResult
#from .BridgePlot import BridgePlot


def combine_error(eString):
    if type(eString) in [list,tuple]:
        sys.stderr.write('\nBridgeCombineError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeCombineError: '+eString+'\n')
    sys.exit(2)
##########################################################################################################################################
##########################################################################################################################################

# bridgeSummary analyze_snp_dists  pdf  print() settings current

class BridgeRun: 
    def __init__(self, bridge):
        self.module, self.args, self.io, self.progress = bridge.io.module, bridge.args, bridge.io, bridge.io.progress 
        
    
    def analyze(self, cmd, results_files, path): 

        self.results_files = results_files 
        self.results = [BridgeResult().read_file(r) for r in results_files] 
        pop_names = list(set([r.pop for r in self.results])) 
        if len(pop_names) > 1: combine_error('Cannot combine across different populations: '+",".join(pop_names)+'\n') 
        else:                  self.pop = pop_names[0]  

        if cmd == 'result': 
            if len(self.results) == 1: 
                self.progress.start_minor('Analyzing/Plotting Single PRS Results (analyze result)', self.io) 
                self.single_plot(self.results, path) 
            
            else: 



        if cmd == 'combine':  
            try:    self.source = results_files[-1].split('/')[-1].split('.')[1].lower()  
            except: self.source = self.pop.lower() 

            self.comboPath = path + '/prs-combined_'+self.source.upper() 
            self.io.progress.start_minor('Combining The Results (analyze combine)', self.io)     
            self.results.extend(self.pred_combine(self.results, path)) 
            if self.args.noplots: sys.stderr.write('\nSkipping Plotting Step!\n') 
            else: 
                self.io.progress.start_minor('Plotting Results', self.io) 
                self.multi_plot(self.results, self.comboPath, path) 
            

    def pred_combine(self, res, path): 
        prsR =      [r for r in res if r.name == 'prs-single'] 
        prsBridge = [r for r in res if r.name == 'prs-prior'] 
        if len(prsR) != 1 or len(prsBridge) != 1: combine_error('Exactly one run for prs-single and prs-bridge required') 
        combine = BridgeCombine(self, self.comboPath) 
        if combine.FIN:     self.io.progress.write('SKIPPING-JOB') 
        else:               combine.run(prsR[0], prsBridge[0]) 
        self.io.progress.end() 
        return BridgeResult().read_combo(combine)

    

    def multi_plot(self, BR, path, sPath):  
        plot_names = [path+'/bridgePRS-combo.pdf', sPath+'/bridgeSummary.'+self.source.lower()+'.pdf'] 
        self.bPlot = BridgePlot(self.args, BR, self.pop, plot_names).setup(TYPE='MEGA') 
        self.bPlot.full_var_bars()
        self.bPlot.add_pred_scatter('weighted')  
        self.bPlot.add_model(BR[1].modelpath) # SNPS = True) 
        self.bPlot.analyze_snp_dists('single','weighted') #,MODEL=True)  
        self.bPlot.add_summary_table('single')  
        self.bPlot.add_logo(2) 
        self.bPlot.finish() 
        return 

         
    def single_plot(self, BR, path): 
        br, method = BR[0], BR[0].name.split('-')[-1] 
        self.bPlot = BridgePlot(self.args, BR, self.pop, [path+'/bridgePRS-'+method+'.pdf']).setup(BR[0].name.split('-')[-1]) 
        self.bPlot.full_var_bars()
        self.bPlot.add_pred_scatter(method) 
        if method != 'single': self.bPlot.add_model(br.modelpath)
        self.bPlot.analyze_snp_dists(method, method) #,MODEL=True)  
        self.bPlot.add_summary_table(method)  
        self.bPlot.add_logo(2) 
        self.bPlot.finish() 
        return 













class BridgeMore:
    def __init__(self, bridge):
        self.module, self.args, self.io, self.progress = bridge.io.module, bridge.args, bridge.io, bridge.io.progress 
        
    def run(self, cmd, prs_results, path): 
        results = [BridgeResult().read_file(r) for r in prs_results] 
        pop_names = list(set([r.pop for r in results])) 
        if len(pop_names) > 1: combine_error('Cannot combine across different populations: '+",".join(pop_names)+'\n') 
        else:                  self.pop = pop_names[0] 
        if cmd == 'result':    self.single_plot(results, path) 
        else: 
            try: 
                self.source_lower = prs_results[-1].split('/')[-1].split('.')[1] 
                self.source_upper = self.source_lower.upper() 
            except: 
                self.source_lower = self.pop.lower() 
                self.source_upper = self.pop.upper() 

            comboPath = path + '/prs-combined_'+self.source_upper  
            results.extend(self.run_combine(results, comboPath)) 
            if self.args.noplots: sys.stderr.write('\nSkipping Plotting Step!\n') 
            else: 
                self.io.progress.start_minor('Plotting Results', self.io) 
                self.mega_plot(results, comboPath, path) 
                
    









































class BridgeCombine: 
    def __init__(self, bridge, path):  
        self.pop, self.module, self.args, self.io, self.progress = bridge.pop, bridge.io.module, bridge.args, bridge.io, bridge.io.progress 
        self.path, self.FIN = path, False 
        if not os.path.exists(self.path): os.makedirs(self.path) 
        else:                             self.finish() 

    def run(self, p1, p2): 
        pp = self.io.programs['pred_combine_en']
        X = ['--pop2',self.pop,'--n.cores',str(self.args.cores)]
        
        p_files, p_name = p1.PT['FILES'].split(','), p1.PT['FIELD-NAME'] 
        


        X.extend(['--test.data',p_files[0],'--valid.data',p_files[1]])
        X.extend(['--outfile',self.path+'/'+self.pop]) 
        X.extend(['--models1',p1.paths['BETA'],'--models2',p2.paths['BETA']])
        X.extend(['--pred1',p1.paths['PREDICT'],'--pred2',p2.paths['PREDICT']]) 
        X.extend(['--ids.col','TRUE']) 
        X.extend(['--pheno.name',p1.PT['FIELD-NAME']]) 
        if 'FIELD-COVARIATES' in p1.PT:  X.extend(['--cov.names',p1.PT['FIELD-COVARIATES']]) 
        #if 'COVARIATES' in p1.phenotype: X.extend(['--cov.names',p1.phenotype['COVARIATES']]) 
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.progress.start_rJob(rJOB, 'combine')
        out_file, err_file = self.path+'/combine.stdout', self.path+'/combine.stderr' 
        my_job = " ".join(rJOB + ['>',out_file,'2>', err_file]) 
        os.system(my_job) 
        self.finish() 
        return 

    
    def finish(self): 
        self.key = {} 
        for f in os.listdir(self.path): 
            if self.pop in f and f.split('.')[-1] not in ['log','png','pdf']: 
                pn = f.split('.')[0].split(self.pop+"_")[-1].strip('_') 
                self.key[pn]  = self.path+'/'+f 
        
        if len(self.key) > 5: self.FIN = True 
        return 
        




