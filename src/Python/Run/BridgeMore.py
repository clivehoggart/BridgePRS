#!/usr/bin/env python3

import os, sys
from collections import defaultdict as dd
from .BridgeResult import BridgeResult
from .BridgePlot import BridgePlot


def combine_error(eString):
    if type(eString) in [list,tuple]:
        sys.stderr.write('\nBridgeCombineError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeCombineError: '+eString+'\n')
    sys.exit(2)
##########################################################################################################################################
##########################################################################################################################################

# bridgeSummary analyze_snp_dists  pdf  print() 


class BridgeMore:
    def __init__(self, bridge):
        self.module, self.args, self.io, self.settings, self.progress = bridge.io.module, bridge.args, bridge.io, bridge.io.settings, bridge.io.progress 
        
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
            if self.args.noPlots: sys.stderr.write('\nSkipping Plotting Step!\n') 
            else: 
                self.io.progress.start_minor('Plotting Results', self.io.settings) 
                self.mega_plot(results, comboPath, path) 
                
    def run_combine(self, res, comboPath):
        prsR =      [r for r in res if r.name == 'prs-single'] 
        prsBridge = [r for r in res if r.name == 'prs-prior'] 
        if len(prsR) != 1 or len(prsBridge) != 1: combine_error('Exactly one run for prs-single and prs-bridge required') 
        combine = BridgeCombine(self, comboPath) 
        if combine.FIN:     self.io.progress.write('SKIPPING-JOB') 
        else:               combine.run(prsR[0], prsBridge[0]) 
        self.io.progress.end() 
        return BridgeResult().read_combo(combine)

    

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


    
    def mega_plot(self, BR, path, sPath): 
        plot_names = [path+'/bridgePRS-combo.pdf', sPath+'/bridgeSummary.'+self.source_lower+'.pdf'] 
        self.bPlot = BridgePlot(self.args, BR, self.pop, plot_names).setup(TYPE='MEGA') 
        self.bPlot.full_var_bars()
        self.bPlot.add_pred_scatter('weighted')  
        self.bPlot.add_model(BR[1].modelpath) # SNPS = True) 
        self.bPlot.analyze_snp_dists('single','weighted') #,MODEL=True)  
        self.bPlot.add_summary_table('single')  
        self.bPlot.add_logo(2) 
        self.bPlot.finish() 
        return 

         


    
    def run_plotter(self, BR, path, summaryPath = None, cols =1): 
        self.bPlot = BridgePlot(self.args, self.pop, prefix = path, summaryPath = summaryPath).setup(cols) 
        
        for i,br in enumerate(BR): 
            self.bPlot.set_result(br,i) 
            if i < 3: 
                self.bPlot.add_prs() 
                self.bPlot.add_table(my_key = 'Ridge') 
            else: 
                z=3
                #JOB2: Plotting Results...1000 0 dict_keys(['names', 'pheno', 'prs.Stage1', 'prs.Stage2', 'prs.Stages1+2', 'prs.weighted', 'prs']) prs
                self.bPlot.add_prs('prs.Stages1+2', my_title = 'prs-combined') 
                self.bPlot.add_table(my_key = 'Stage1+2') 
                
                self.bPlot.add_prs('prs.weighted', my_title = 'prs-weighted') 
                self.bPlot.add_table(my_key = 'Weighted') 
        self.bPlot.finish_plot() 
        return 

         




class BridgeCombine: 
    def __init__(self, bridge, path):  
        self.pop, self.module, self.args, self.io, self.settings, self.progress = bridge.pop, bridge.io.module, bridge.args, bridge.io, bridge.io.settings, bridge.io.progress 
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
        




