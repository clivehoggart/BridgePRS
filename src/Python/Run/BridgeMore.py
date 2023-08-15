#!/usr/bin/env python3

import os, sys
from collections import defaultdict as dd

try: 
    import matplotlib
    import matplotlib.pyplot as plt
    MAKEPLOT=True 
except: 
    MAKEPLOT=False 


def combine_error(eString):
    if type(eString) in [list,tuple]:
        sys.stderr.write('\nBridgeCombineError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeCombineError: '+eString+'\n')
    sys.exit(2)
##########################################################################################################################################
##########################################################################################################################################




class BridgeMore:
    def __init__(self, bridge):
        self.module, self.args, self.io, self.settings, self.progress = bridge.io.module, bridge.args, bridge.io, bridge.io.settings, bridge.io.progress 
        
    def run(self, cmd, prs_results, path): 
        results = [BridgeResult().read_file(r) for r in prs_results] 
        pop_names = list(set([r.pop for r in results])) 
        if len(pop_names) > 1: combine_error('Cannot combine across different populations: '+",".join(pop_names)+'\n') 
        else:                  self.pop = pop_names[0] 
        if cmd == 'result':  self.run_plotter(results, path, cols = len(results)) 
        else:
            comboPath = path + '/prs-combined_'+self.pop 
            results.append(self.run_combine(results, comboPath)) 
            
            if not self.args.skipAnalysis: 
                if MAKEPLOT: 
                    self.io.progress.start_minor('Plotting Results') 
                    self.run_plotter(results, comboPath, cols = 5) 
                else: sys.stderr.write('\nWarning: Module matplotlib not found (plots not created)!\n') 
            return  
                
        
    def run_combine(self, res, comboPath): 
        prsR =      [r for r in res if r.method == 'prs-single'] 
        prsBridge = [r for r in res if r.method == 'prs-prior'] 
        if len(prsR) != 1 or len(prsBridge) != 1: combine_error('Exactly one run for prs-single and prs-bridge required') 
        
        combine = BridgeCombine(self, comboPath) 
        if combine.FIN and not self.args.repeatSteps: self.io.progress.update_minor('SKIPPING-JOB')
        else:                                         combine.run(prsR[0], prsBridge[0]) 
        self.io.progress.end() 
        br = BridgeResult().read_combo(combine) 
        return br 

    
    def run_plotter(self, BR, path, cols =1): 
        self.bPlot = BridgePlot(self.args, self.pop, prefix = path).setup(cols) 
        for i,br in enumerate(BR): 
            self.bPlot.set_result(br,i) 
            self.bPlot.add_prs() 
            self.bPlot.add_table(my_key = 'Ridge') 
            if i == 3: 
                self.bPlot.add_prs('weighted', my_title = 'prs-weighted') 
                self.bPlot.add_table(my_key = 'Ridge.w') 
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
        X = ['--pop2',self.pop, '--test.data', p1.files['PHENO'], '--valid.data', p1.files['VALIDATION'], '--n.cores',str(self.args.cores)] 
        X.extend(['--outfile',self.path+'/'+self.pop]) 
        X.extend(['--pred1',p1.paths['PREDICT']]) 
        X.extend(['--pred2',p2.paths['PREDICT']]) 
        #X.extend(['--eval1',p1.paths['EVAL'],'--eval2',p2.paths['EVAL']])
        X.extend(['--models1',p1.paths['EVAL'],'--models2',p2.paths['EVAL']])
        X.extend(['--ids.col','TRUE']) 
        X.extend(['--pheno.name',p1.fields['NAME'], '--cov.names',p1.fields['COVARIATES']])
        rJOB = ['Rscript','--vanilla',pp,'--fpath',self.io.programs['functions']] + X
        self.progress.start_rJob(rJOB)
        out_file, err_file = self.path+'/combine.stdout', self.path+'/combine.stderr' 
        my_job = " ".join(rJOB + ['>',out_file,'2>', err_file]) 
        os.system(my_job) 
        self.finish() 
        return 

    
    def finish(self): 
        self.key = {} 
        for f in os.listdir(self.path): 
            if self.pop in f and f.split('.')[-1] != 'log': 
                pn = f.split('.')[0].split(self.pop+"_")[-1].strip('_') 
                self.key[pn]  = self.path+'/'+f 
        if len(self.key) > 5: self.FIN = True 
        return 
        




class BridgeResult:
    def __init__(self, name = None):
        self.name   = name
        self.fields = dd(bool) 



    def read_combo(self, combo): 
        self.pop = combo.pop 
        self.method = 'prs-combine' 
        for k,f in combo.key.items():
            if k in ['weighted_combined_preds']: self.read_prs_predictions(f) 
            elif k in ['weighted_combined_var_explained']: self.read_var_explained(f) 
            else: continue 
        return self 



    def read_file(self, res): 
        self.key, self.paths, self.files = self.read_result(res)
        for k,f in self.key.items():
            if k in ['preds','prs_predictions']: self.read_prs_predictions(f) 
            elif k in ['var_explained']: self.read_var_explained(f) 
            else: continue  
        return self 
    
    
    def read_result(self, res): 
        K, P, F = {} , dd(bool), dd(bool) 
        self.model = None 
        f = open(res) 
        for line in f: 
            a,b = line.strip().split('=') 
            if a == 'POP_NAME':      self.pop = b 
            elif a == 'MODULE_NAME': self.method = b
            elif a.split('_')[-1] == 'FILE': F[a.split('_')[0]] = b 
            elif a.split('_')[-1] == 'PREFIX': P[a.split('_')[1]] = b 
            elif a.split('_')[0] == 'FIELD' and a.split('_')[1].split('-')[0] == 'PHENO':  self.fields[a.split('-')[-1]] = b 
            else: continue 
        f.close()
        q_path, qn =  "/".join(P['QUANTIFY'].split('/')[0:-1]), P['QUANTIFY'].split('/')[-1] 
        for f in os.listdir(q_path): 
            if qn in f and f.split('.')[-1] != 'log': 
                pn = f.split('.')[0].split(qn)[-1].strip('_') 
                K[pn] = q_path+'/'+f 
        return K, P, F  



    def read_var_explained(self,res): 
        f = open(res)
        self.var_key = dd(lambda: {}) 
        header = [x.strip("\"") for x in f.readline().strip().split(',')][1::] 
        for i,line in enumerate(f): 
            line = line.strip().split(',') 
            name, vals = line[0].strip('\"'), [float(x) for x in line[1::]] 
            for h,v in zip(header,vals): self.var_key[name][h] = v 
        f.close() 
        return 



    def read_prs_predictions(self, res): 
        f = open(res) 
        self.preds = dd(list) 
        header = [x.strip('\"') for x in f.readline().split()]

        if header[0] == '---': header[0] = 'names' 
        for i,line in enumerate(f): 
            line = line.split() 
            line = [float(x) if j > 0 else x for j,x in enumerate(line)]
            for h,v in zip(header,line): self.preds[h].append(v) 
        return 
        





class BridgePlot: 
    def __init__(self, args, pop, prefix):  
        self.args, self.pop, self.prefix, self.methods  = args, pop, prefix+'/plot_'+pop, []  
        
        
        


    # SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP #
    def setup(self, cols): 
        self.axes, self.ax_index, self.rows, self.cols = [] , 0, 3, cols  
        self.fig, self.HT, self.WD, self.fs1, self.fs2 = matplotlib.pyplot.gcf(), 8 + (0.25*cols), 10 + (8*cols), 25, 22 
        for j in range(0,self.cols): 
            self.axes.append(plt.subplot2grid((self.rows,self.cols), (0,j), rowspan = 2, colspan =1))
            self.axes.append(plt.subplot2grid((self.rows,self.cols), (2,j), rowspan = 1, colspan =1))
        self.fig.set_size_inches(self.WD, self.HT)  
        return self  


    def set_result(self, R, index = 0): 
        self.R, self.R_idx = R, index 
        self.pop, self.method = R.pop, R.method 
        self.methods.append(R.method) 
        return

    def finish_plot(self, suffix = 'NA'): 
        if len(self.methods) == 1: self.p_name = self.prefix+'_'+self.methods[0]+'.png' 
        elif 'prs-combine' not in self.methods: self.p_name = self.prefix+'_prs_independent_runs.png'
        else:                                   self.p_name = self.prefix+'_prs_combined_runs.png'
        plt.subplots_adjust(left=0.04, bottom=0.01, right=0.96, top=0.93,wspace=0.20,hspace=0.20) 
        plt.savefig(self.p_name,dpi=200) 
        plt.clf() 
        plt.close() 
        return


    def add_prs(self, k = 'prs',my_title=None): 
        ax = self.axes[self.ax_index] 
        self.ax_index += 1
        X = self.R.preds['pheno'] 
        Y = self.R.preds[k]
        ax.scatter(X,Y)
        ax.set_xticks([]) 
        ax.set_yticks([]) 
        ax.set_xlabel('Phenotypes', fontsize = self.fs2) 
        ax.set_ylabel('PRS', fontsize = self.fs2) 
        if my_title is None: ax.set_title(self.pop+': '+self.method, fontsize = self.fs1)
        else:                ax.set_title(self.pop+': '+my_title, fontsize = self.fs1) 
        return

    def add_table(self, my_key = None): 

        ax = self.axes[self.ax_index] 
        self.ax_index += 1 
        dt = DataTable(self, ax) 
        dt.add_prs_result(self.R, key = my_key, fs1 = self.fs1, fs2 = self.fs2)  
        return






class DataTable:
    def __init__(self, plot, ax,  clr='white'): 
        self.plot, self.args = plot, plot.args 
        self.ax, self.rows, self.cols, self.clr  = ax, [], [], clr 
        self.widths = [22,12,10,11,11,17,17,10,10,10]
        self.vwid = [0] + [sum(self.widths[0:i+1]) for i in range(len(self.widths))]   

    def get_loc(self,X,Y):
        x1,x2 = X[0]/100.0 , X[1]/100.0
        y1,y2 = Y[0]/100.0 , Y[1]/100.0
        return (x1,y1,x2-x1,y2-y1)

    def add_row(self,row_data,X=None,Y=None,COLORS=[],WIDTHS=[],FS=13,ALPHA=0,TITLE=False, CL = 'center',CLEAR=False): 
        if X == None: X = self.X_SPAN
        if Y == None: Y = self.Y_SPAN
        cL,rL,rD,xL = CL,None,[row_data],len(row_data)  
        bl = self.get_loc(X,Y) 
        while len(WIDTHS) < len(row_data): WIDTHS.append(10) 
        while len(COLORS) < len(row_data): COLORS.append('white') 
        if CL != 'center': row = self.ax.table(cellText=rD,rowLabels=rL,cellColours = [COLORS[0:xL]],colWidths=WIDTHS[0:xL],bbox=bl,loc = cL, cellLoc=cL, alpha=ALPHA, clip_on=False) 
        else:              row = self.ax.table(cellText=rD,rowLabels=rL,cellColours = [COLORS[0:xL]],colWidths=WIDTHS[0:xL],bbox=bl,cellLoc=cL, alpha=ALPHA, clip_on=False) 
        row.auto_set_font_size(False)
        row.set_fontsize(FS)
        table_props = row.properties()
        if ALPHA > 0: 
            for cell in row._cells: row._cells[cell].set_alpha(ALPHA)
        self.rows.append(row) 



    def add_prs_result(self, R, fs1=18, fs2=14, fs3 = 12, key = None): 
        vk = ['Est', '2.5%', '97.5%', 'cv.dev', 'cv.dev.sd']
        data_key = {rk: [round(R.var_key[rk][v],4) for v in vk] for rk in R.var_key.keys()} 
        header = ['variance\nexplained','95\%\nConfidence\nInterval', 'CV\ndeviations', 'CV\ndeviation\nstd'] 
        self.add_row(header, COLORS=['whitesmoke' for h in header], X = (0,100), Y = (50,100), FS = fs2, WIDTHS = [22,34,22,22])         
        my_data = data_key[key] 
        self.add_row(data_key[key], COLORS=['snow' for k in data_key[key]], X = (0,100), Y = (20,50), FS = fs2, WIDTHS = [22,17,17,22,22])         
        self.ax.axis('off') 
        return 






