import os, sys, gzip
from collections import defaultdict as dd
from math import log 
from . import RunTools as rtools  
try: 
    import matplotlib 
    import matplotlib.pyplot as plt 
except: pass 

# bar Result 
# print 


##########################################################################################################################################




        



class BridgePlot: 
    def __init__(self, args, BR, pop, fig_names):  

        self.args, self.pop, self.fig_names  = args, pop, fig_names 
        self.fig, self.axes, self.ax_index, self.fs1, self.fs2 = matplotlib.pyplot.gcf(), [], 0, 25, 22 
        self.results = BR 
        self.names = [bR.name.split('-')[-1] for bR in BR] 
        self.data  = {bR.name.split('-')[-1]: bR for bR in BR}                                                                                                                                                                                                                           
        #self.col_headers = BR[0].SL  





    # SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP #
    def setup(self, TYPE, cols = 1): 
        self.TYPE, self.MODEL, self.DATA_KEY = TYPE, False, dd(lambda: dd(lambda: 'NA')) 
        

        c1, c2 = 6,7 
        r1, r2,r3 = 3,3,3
        self.rows, self.cols, self.ax_index, self.HT, self.WD = r1+r2+r3+r3+r3,c1+c2, 1, 15, 16 
    

        for j,i in enumerate([0, r1, r1+r2, r1+r2+1, r1+r2+r3+1]): 
            if j == 0: 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (0,0), rowspan = 1, colspan = 2)) 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (0,2), rowspan = 1, colspan = 4)) 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (1,0), rowspan = r1+2, colspan = 6)) 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,c1), rowspan = r1, colspan = c2)) 
            elif j == 1: 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,c1), rowspan = r2, colspan = c2)) 
            elif j == 2: 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,0), rowspan = 1, colspan = self.cols)) 

            else: 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,0), rowspan = r3, colspan = self.cols)) 
        self.fig.set_size_inches(self.WD, self.HT)  
        return self  


    
    def fill_in(self,method1, method2): 
        
        
        for i,ax in enumerate(self.axes): 
            
            if i > 6: 
                ax.set_xticks([]) 
                ax.set_yticks([]) 
                ax.set_title(str(i+1)) 
        self.full_var_bars(3) 
        self.add_pred_scatter(method2,4)
        
        self.axes[5].axis('off') 

        if self.TYPE == 'MEGA' and len(self.results) > 1:  self.add_model(self.results[1].modelpath) 
        elif method1 != 'single':                          self.add_model(self.results[0].modelpath) 
        
        self.ax_index = 6
        self.analyze_snp_dists(method1, method2) 
         
       # dt = rtools.DataTable(self, self.axes[2]).add_bridge_summary(method1) 
        self.add_summary_table(method1, 1,2)
        self.add_logo(0)                                                                                                                                                                                             
        self.finish()     

    
    
    def finish(self):  
        if self.TYPE.upper() == 'MEGA': tString = 'bridgePRS: Weighted-Combined Analysis' 
        else:                   tString = 'bridgePRS: PRS-'+self.TYPE.upper()+' Analysis'
        plt.suptitle(tString, fontweight='bold', fontsize=22) 
        
        plt.subplots_adjust(left=0.04, bottom=-0.10, right=0.95, top=0.93,wspace=0.0001,hspace=0.001) 
        
        for fn in self.fig_names: 
            if fn.split('.')[-1] == 'pdf': plt.savefig(".".join(fn.split('.')[0:-1])+'.png',dpi=200) 
            plt.savefig(fn,dpi=200) 
        plt.clf() 
        plt.close()
        return


    def add_model(self, mp): 
        if not(mp): return 
        self.MODEL = True 
        with open(mp, 'rt') as my_file: self.model_key = {a: b for a,b in [line.strip().split('=') for line in my_file.readlines()]} 
        return 
        


    def full_var_bars(self, idx = None):
        if idx is None: 
            idx = self.ax_index
            ax = self.axes[self.ax_index] 
            self.ax_index += 1
        else: ax = self.axes[idx] 

        cands       = ['port','single', 'prior', 'combine', 'weighted'] 
        name_key    = ['port','stage-1','stage-2','stage-1+2','weighted'] 
        self.colors = ['tab:gray','yellow','cyan','lime','tab:green']  
        color_key   = {m: c for m,c in zip(cands, self.colors)}
        b_vals = [self.data[c].varexp.interval[-1] for c in self.names] 
        yMax, ys = round(max(b_vals)+0.01,2), 0 - (max(b_vals) * 0.04) 
        for i,n in enumerate(cands): 
            if n not in self.names: continue 
            ht, [a,b] = self.data[n].varexp.val, self.data[n].varexp.interval
            if a < 0: a = 0 
            for xp,yp in [[[i,i],[a,b]],[[i-0.1,i+0.1],[a,a]],[[i-0.1,i+0.1],[b,b]]]: ax.plot([xp[0],xp[1]],[yp[0],yp[1]], color='k',zorder=2, lw=1.5) 
            ax.plot([i,i],[ys*0.3,0], color='k', lw = 1, zorder=0)  
            ax.text(i,ys,name_key[i], ha='center',va='top',fontsize=13) 
            if n != 'weighted':  ax.bar(i,ht, color = self.colors[i], edgecolor='k', alpha = 0.65, zorder = 2) 
            else:
                ax.bar(i,ht, color = 'white', edgecolor='k', alpha = 0.5, zorder = 1) 
                cuts, btm, FK = [], 0, [[c,self.data[n].varexp.frac[c], color_key[c]] for c in cands]
                for cand, frac, clr in FK: 
                    if cand == n or frac == 0: continue 
                    cuts.append([ht*frac, clr, cand]) 
                for piece, clr, cand in cuts: 
                    ax.bar(i, piece, bottom = btm , color = clr, alpha = 0.7, zorder = 2) 
                    btm += piece 

       
        for xp,yp in [[[-0.6,4.6],[0,0]],[[4.6,4.6],[0,yMax*1.1]]]: ax.plot([xp[0],xp[1]],[yp[0],yp[1]], color='k', clip_on=False, zorder=2) 
        if yMax < 0.07: yCuts = [x/100.0 for x in range(1,10,1) if x/100.0 < yMax] 
        elif yMax < 0.15:  yCuts = [x/100.0 for x in range(2,20,2) if x/100.0 < yMax+0.01] 
        else:  yCuts = [x/100.0 for x in range(5,100,2) if x/100.0 < yMax] 
        for yc in yCuts: ax.text(4.60, yc, '-'+str(yc)) 
        ax.set_title('Model Performance', fontweight='bold',fontsize=20, y = 0.9)  
        ax.text(5.15,yMax/1.5, 'Variance Explained',rotation=-90, va='center',ha='center',fontsize=15)  
        ax.axis('off') 
        yMin, yMax = ax.get_ylim()
        ax.set_xlim(-1.4,5.1) 
        ax.set_ylim(0 - (yMax/10.0), yMax) 
        return 
         





    def add_pred_scatter(self, method, idx = None):
        if idx is None: 
            ax = self.axes[self.ax_index] 
            self.ax_index += 1 
        else:
            ax = self.axes[idx] 
        X = self.data[method].preds['prs'] 
        Y = self.data[method].preds['pheno'] 
        
        if len(Y) > 20 and len(list(set(Y[0:20]))) == 2: 
            
            cases, controls = [], [] 
            colors = ['gray', 'red', 'tomato']
            for x,y in zip(X,Y): 
                if y == 0.0: controls.append(x) 
                elif y == 1.0: cases.append(x)  
                else: continue 

            bplot = ax.boxplot([controls,cases], 0, 'rs', 0, patch_artist=True) 

            for patch, color in zip(bplot['boxes'] , colors): 
                patch.set_facecolor(color)
            
            self.get_lims(ax, BORDER=4, xLab = 'PRS\n('+method+')', yLab = None)
            ax.set_xlim(self.xMin - self.xSpan/6.0, self.xMax + self.xSpan/13.0) 
            ax.set_ylim(self.yMin - self.ySpan/15.0, self.yMax + self.ySpan/3.0) 
            ax.set_title('PRS Predictions', fontweight='bold',fontsize=20, ha='center',x=0.52, y = 0.79)  
            
            ax.text(self.xMax,2,'-Cases') 
            ax.text(self.xMax,1,'-Controls') 
            

        else: 
            try:    ax.scatter(X,Y, alpha=0.5) 
            except TypeError: 
                Xm, Ym, Xp, Yp = [], [], [], []     
                for x,y in zip(X,Y): 
                    if x != "NA" and y != "NA": 
                        Xp.append(x) 
                        Yp.append(y) 
                    elif y == 'NA': Xm.append(x) 
                    elif x == 'NA': Ym.append(y) 
                    else:           continue 
                ax.scatter(Xp, Yp, alpha= 0.5) 
                xMin, xMax = ax.get_xlim() 
                yMin, yMax = ax.get_ylim() 
                if len(Xm) > 0: ax.scatter(Xm, [yMin for x in Xm], color = 'orange', alpha=0.6) 
                if len(Ym) > 0: ax.scatter([xMin for y in Ym], Ym,  color = 'orange', alpha=0.6) 
            self.get_lims(ax, BORDER=4, xLab = 'PRS\n('+method+')', yLab = 'Phenotype')
            ax.set_xlim(self.xMin - self.xSpan/6.0, self.xMax + self.xSpan/13.0) 
            ax.set_ylim(self.yMin - self.ySpan/15.0, self.yMax + self.ySpan/3.0) 
            ax.set_title('PRS Predictions', fontweight='bold',fontsize=20, ha='center',x=0.51, y = 0.79)  



    def analyze_snp_dists(self, method, method2):

        base_scores, base_key, self.DATA_KEY['LEN']['target'] = self.load_base_scores(self.data[method].SS['PREFIX'],self.data[method].SS['SUFFIX'],[self.data[method].SL[k] for k in ['ID','P','BETA']], idx=0) 
        self.draw_manhattan(base_scores, 'Target GWAS')
        model_scores, model_key = None, 'NA' 
        if self.MODEL: model_scores, model_key, self.DATA_KEY['LEN']['model'] = self.load_base_scores(self.model_key['SUMSTATS_PREFIX'],self.model_key['SUMSTATS_SUFFIX'],[self.model_key[k] for k in ['SSFIELD_ID','SSFIELD_P','SSFIELD_BETA']],idx=1)
        else:          model_scores, model_key = None, 'NA' 
        self.draw_manhattan(model_scores, 'Model GWAS') 
        try: 
            snp_weights = self.data[method2].snp_weights 
            w_bb =  self.merge_snp_scores(base_key, snp_weights, model_key) 
            w, b1 = [x[0] for x in w_bb],[x[1] for x in w_bb]
            self.DATA_KEY['LEN']['weight'] = len(w) 
            self.DATA_KEY['CORR'] = {'wTarget': rtools.R_CORR(w, b1), 'wModel': 'NA', 'MT': 'NA'} 
            if self.MODEL: 
                b2 = [x[2] for x in w_bb] 
                self.DATA_KEY['CORR'] = {'wTarget': rtools.R_CORR(w, b1), 'wModel': rtools.R_CORR(w, b2), 'MT': rtools.R_CORR(b1,b2)} 
        except: 
            pass  
        return 
        

       
            
 

    def merge_snp_scores(self, base, weights, model_key = 'NA'): 
        ssw_flat = [] 
        chromosomes = sorted(base.keys())
        for ci,c in enumerate(sorted(base.keys())):
            if model_key == 'NA': ssw_flat.extend([[weights[k], base[c][k][1]] for k in base[c] if k in weights]) 
            else:                 ssw_flat.extend([[weights[k], base[c][k][1], model_key[c][k][1]] for k in base[c] if k in weights and k in model_key[c]]) 

        return ssw_flat








        
    def get_lims(self, ax, BORDER=0, xLab = None, yLab = None, xLims = [], yLims = [] ):  
        self.ax = ax 
        if len(xLims) == 2: self.xMin,self.xMax = xLims
        else:               self.xMin,self.xMax = ax.get_xlim() 
        if len(yLims) == 2: self.yMin,self.yMax = yLims  
        else:               self.yMin,self.yMax = ax.get_ylim() 
        
        
        
        self.xSpan, self.ySpan = self.xMax - self.xMin, self.yMax - self.yMin
        if self.xMin != self.xMax: ax.set_xlim(self.xMin, self.xMax) 
        if self.yMin != self.yMax: ax.set_ylim(self.yMin, self.yMax) 
        if  BORDER > 0: 
            ax.axis('off') 
            ax.plot([self.xMin, self.xMax],[self.yMin, self.yMin], color='k') 
            ax.plot([self.xMax, self.xMax],[self.yMin, self.yMax], color='k', linewidth=3) 
            if BORDER > 2: 
                ax.plot([self.xMin, self.xMin],[self.yMin, self.yMax], color='k') 
                ax.plot([self.xMin, self.xMax],[self.yMax, self.yMax], color='k') 
        
        if xLab is not None:  ax.text(self.xMin + (self.xMax-self.xMin)/2.0, self.yMin - self.ySpan/25.0, xLab,ha='center', va = 'top', fontsize=15, clip_on=False)  
        if yLab is not None:  
            if yLab == 'Phenotype': ax.text(self.xMax + self.xSpan/33.0, self.yMin + self.ySpan/2.0, yLab,ha='center', rotation = -90, va = 'center', fontsize=15, clip_on=False)  
            else: ax.text(self.xMin - self.xSpan/25.0, self.yMin + self.ySpan/2.0, yLab,ha='center', rotation = 90, va = 'center', fontsize=15, clip_on=False)  
            

    def draw_manhattan(self, scores, Tx): 
        ax = self.axes[self.ax_index] 
        self.ax_index += 1 
        if scores is None or len(scores) == 0:  
            ax.axis('off') 
            return
        
        chromosomes = sorted(scores.keys())
        maxY, minY, xticks, xlabs, p_offset = 0, 100,  [] , [], 0 
        for ci,c in enumerate(chromosomes):
            bs = scores[c] 
            
            xlocs = [b[0] for b in bs] 
            if xlocs[0] == xlocs[-1]: xlocs = [i for i in range(len(xlocs))] 
            X,Y = [p_offset+xp for xp in xlocs], [-log(b[2],10) for b in bs]
                



            if ci == 0: minX = X[0] 
            if max(Y) > maxY: maxY = max(Y) 
            if min(Y) < minY: minY = min(Y) 
            if ci % 2 == 0: colors = ['tab:pink', 'tab:cyan'] 
            else:           colors = ['red', 'blue']
            c1, c2 = colors          
            ax.scatter(X,Y, clip_on=False) 
            xm = (p_offset + X[-1]) / 2
            xticks.append(xm) 
            xlabs.append(c) 
            p_offset = X[-1]+1 
        
        p_range = p_offset - minX 
        self.get_lims(ax, BORDER=2, xLims = [0, p_offset+p_range/75.0])  
        for x,p in zip(xticks, xlabs): ax.text(x,0-self.ySpan/20.0,p, va = 'top', ha = 'center', fontsize=11) 
        if self.yMax < 5:     step = 2 
        elif self.yMax < 20:  step  = 5 
        elif self.yMax < 50:  step = 10 
        elif self.yMax < 100: step = 25 
        else:            step = 50 
        pt = step 
        while pt < self.yMax: 
            ax.text(self.xMax+(0*self.xSpan/10.0), pt, int(pt), ha='left', clip_on = False) 
            pt += step 
        ax.text(self.xMax+(self.xSpan*0.011), self.yMin + self.ySpan/18.0, 'logP', ha='left',  rotation = 90, fontsize=15,clip_on = False) 
        ax.axis('off') 
        ax.set_title('Manhattan Plot: '+Tx, x = 0.5, y = 0.82, fontsize=20, fontweight='bold')  
        ax.set_ylim(0 - (self.yMax/8.0), self.yMax*1.2) 
        return 












    def load_base_scores(self, f_prefix, suffix, h_triple, idx=0): 
        path   = "/".join(f_prefix.split('/')[0:-1])
        prefix = f_prefix.split('/')[-1] 
        full_len, full_scores, full_key    = 0, {}, {} 
        for f in os.listdir(path): 
            if prefix in f and suffix in f:  
                chr_key, chr_snps = {}, [] 
                chr_name = int(f.split(prefix)[1].split(suffix)[0]) 
                H, gf = rtools.zip_open(path+'/'+f)
                HK = {h: i for i,h in enumerate(H)} 
                SPB = [HK[cn] for cn in h_triple]
                if len(SPB) == 3: 
                    for line in gf: 
                        line = line.split() 
                        snp, pv, beta = [float(line[k]) if i > 0 else line[k] for i,k in enumerate(SPB)]
                        chr_snps.append([0, snp, pv, beta]) 
                        chr_key[snp] = [pv, beta, 0] 
                
                    full_scores[chr_name] = sorted(chr_snps)  
                    full_key[chr_name] = chr_key 
                    full_len += len(chr_snps) 
        return full_scores,  full_key, full_len 

    
    
    def add_summary_table(self, method, idx1, idx2): 
        #ax = self.axes[idx] 

        dt1 = rtools.DataTable(self, self.axes[idx1]).add_summary_header(method) 
        dt2 = rtools.DataTable(self, self.axes[idx2]).add_smart_summary(method) 



    def add_logo(self, ax_index): 
        ax = self.axes[ax_index]
        bp = rtools.BridgePic(ax).draw(self.TYPE.upper()) 


        #if self.TYPE == 'MEGA': val = 0 
        #elif self.TYPE == 'single': val = 1
        #elif self.TYPE == 'port': val =  2 
        #elif self.TYPE == 'prior': val = 3 
        #else: val = random.randrange(5) 











