#!/usr/bin/env python3

import os, sys
from collections import defaultdict as dd
from .Bridge_Run.BridgeResult import BridgeResult 
from .Bridge_Run.BridgePlot   import BridgePlot 
from .Bridge_Run import RunTools as rtools  



def combine_error(eString):
    if type(eString) in [list,tuple]:
        sys.stderr.write('\nBridgeCombineError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeCombineError: '+eString+'\n')
    sys.exit(2)
##########################################################################################################################################
##########################################################################################################################################


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
                self.progress.start_minor('Analyzing/Plotting PRS Results (analyze result)', self.io) 
                self.single_plot(self.results, path) 
        if cmd == 'combine':  
            try:    self.source = results_files[-1].split('/')[-1].split('.')[1].lower()  
            except: self.source = self.pop.lower() 
            r1, r2 = [r for r in self.results if r.name == 'prs-single'], [r for r in self.results if r.name == 'prs-prior'] 
            if len(r1) != 1 or len(r2) != 1: self.progress.fail('Improper Input, Analyze Combine Requires Exactly one stage-1 one stage-2 result') 
            self.comboPath = path + '/prs-combined_'+self.source.upper() 
            progress_file = self.comboPath+'/'+'bridge.'+self.source.lower()+'.prs-combined.result'
            combine = BridgeCombine(self, self.comboPath, progress_file) 
            self.io.progress.start_minor('Combining The Results (analyze combine)', self.io)     
            if combine.FIN: self.io.progress.write('SKIPPING-JOB') 
            else:           combine.run(r1[0], r2[0])     
            self.io.progress.end()
            self.results = [r1[0], r2[0]] + BridgeResult().read_combo(combine) 

            if self.args.noplots: sys.stderr.write('\nSkipping Plotting Step!\n') 
            else: 
                self.io.progress.start_minor('Plotting Results', self.io) 
                self.multi_plot(self.results, self.comboPath, path) 
            

    def pred_combine2(self, prs1, prs2, path): 
        combine = BridgeCombine(self, self.comboPath) 
        if combine.FIN:     self.io.progress.write('SKIPPING-JOB') 
        else:               combine.run(prs1, prs2) 
        self.io.progress.end()
        return BridgeResult().read_combo(combine)

    

    def multi_plot(self, BR, path, sPath):  
        plot_names = [path+'/bridgePRS-combo.pdf', sPath+'/bridgeSummary.'+self.source.lower()+'.pdf'] 
        self.bPlot = BridgePlot(self.args, BR, self.pop, plot_names).setup(TYPE='MEGA') 
        self.bPlot.fill_in('single','weighted') 
        
        return 

         
    def single_plot(self, BR, path): 
        br, method = BR[0], BR[0].name.split('-')[-1] 
        self.bPlot = BridgePlot(self.args, BR, self.pop, [path+'/bridgePRS-'+method+'.pdf']).setup(BR[0].name.split('-')[-1]) 
        self.bPlot.fill_in(method, method) 

        #self.bPlot.full_var_bars()
        #self.bPlot.add_pred_scatter(method) 
        #if method != 'single': self.bPlot.add_model(br.modelpath)
        #self.bPlot.analyze_snp_dists(method, method) #,MODEL=True)  
        #self.bPlot.add_summary_table(method)  
        #self.bPlot.add_logo(2) 
        #self.bPlot.finish() 
        return 







class BridgeCombine: 
    def __init__(self, bridge, path, progress_file):  
        self.pop, self.module, self.args, self.io, self.progress = bridge.pop, bridge.io.module, bridge.args, bridge.io, bridge.io.progress 
        self.path, self.progress_file, self.FIN = path, progress_file, False 
        if not os.path.exists(self.path): os.makedirs(self.path) 
        self.key = self.read_progress(progress_file) 


    def read_progress(self, pf): 
        key = {}
        if os.path.isfile(pf): 
            with open(pf) as my_file: 
                for a,b in [rl.strip().split('=') for rl in my_file.readlines()]: key[a] = b  
        if 'COMBINE_FIN' in key and key.pop('COMBINE_FIN') == 'TRUE': self.FIN = True 
        return key 



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
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.progress.start_rJob(rJOB, 'combine')
        out_file, err_file = self.path+'/combine.stdout', self.path+'/combine.stderr' 
        my_job = " ".join(rJOB + ['>',out_file,'2>', err_file]) 
        os.system(my_job) 
        self.finish() 
        return 

    
    def finish(self): 
        w = open(self.progress_file, 'w') 
        for f in os.listdir(self.path): 
            if self.pop in f and f.split('.')[-1] not in ['log','png','pdf']: 
                pn = f.split('.')[0].split(self.pop+"_")[-1].strip('_') 
                self.key[pn]  = self.path+'/'+f 
                w.write(pn.upper()+'='+self.path+'/'+f+'\n') 
        if len(self.key) > 5: w.write('COMBINE_FIN=TRUE\n') 
        w.close() 
        return 
        




