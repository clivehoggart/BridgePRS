#!/usr/bin/env python3

import os, sys, gzip
from collections import defaultdict as dd
import matplotlib 
import matplotlib.pyplot as plt 
from math import log 
import numpy as np 
from scipy import stats
import math
import random 
from matplotlib.patches import Rectangle as Rect


def combine_error(eString):
    if type(eString) in [list,tuple]:
        sys.stderr.write('\nBridgeCombineError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeCombineError: '+eString+'\n')
    sys.exit(2)
##########################################################################################################################################
##########################################################################################################################################



def zip_open(fp, HEADER = True):                                                                                                                                                                                                                                                       
    if fp.split('.')[-1] == 'gz': gf = gzip.open(fp, 'rt')                                                                                                                                                                                                                              
    else:                         gf = open(fp, 'rt')                                                                                                                                                                                                                                   
    if HEADER:
        HL = gf.readline().split()                                                                                                                                                                                                                                                      
        return HL, gf 
    return gf 


def S_CORR(x,y): 
    Rs, pv = stats.spearmanr(x, y)
    Rs = str(round(Rs,2)) 
    pv = '%5.1es ' % pv
    return Rs 



class BridgePlot: 
    def __init__(self, args, BR, pop, fig_names):  
        self.args, self.pop, self.fig_names  = args, pop, fig_names 
        self.fig, self.axes, self.ax_index, self.fs1, self.fs2 = matplotlib.pyplot.gcf(), [], 0, 25, 22 
        self.names = [bR.name.split('-')[-1] for bR in BR] 
        self.data  = {bR.name.split('-')[-1]: bR for bR in BR}                                                                                                                                                                                                                           
        


    # SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP #
    def setup(self, TYPE, cols = 1): 
        
        self.TYPE, self.MODEL, self.DATA_KEY = TYPE, False, dd(lambda: dd(lambda: 'NA')) 
        
        self.rows, self.cols, self.ax_index, self.HT, self.WD = 4,2, 1, 16, 15 
        for i in range(0, self.rows): 
            if i == 0: 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,0), rowspan = 1, colspan = 1)) 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,1), rowspan = 1, colspan = 1)) 
            elif i == 1: 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,0), rowspan = 1, colspan = 1)) 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,1), rowspan = 1, colspan = 1)) 
            else: self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,0), rowspan = 1, colspan = 2)) 
    
        self.fig.set_size_inches(self.WD, self.HT)  
        return self  






    
    
    def finish(self): 
        
        if self.TYPE == 'MEGA': tString = 'bridgePRS: Weighted-Combined Analysis' 
        else:                   tString = 'bridgePRS: PRS-'+self.TYPE.upper()+' Result'
        plt.suptitle(tString, fontweight='bold', fontsize=22) 
        if self.TYPE == 'MEGA': plt.subplots_adjust(left=0.04, bottom=0.05, right=0.95, top=0.93,wspace=0.02,hspace=0.001) 
        else:                   plt.subplots_adjust(left=0.04, bottom=0.05, right=0.95, top=0.93,wspace=0.02,hspace=0.001) 
        for fn in self.fig_names: plt.savefig(fn,dpi=200) 
        plt.clf() 
        plt.close()
        return


    def add_model(self, mp): 
        self.MODEL, f = True, open(mp, 'rt') 
        self.model_key = {a: b for a,b in [line.strip().split('=') for line in f.readlines()]} 
        f.close() 
        return 
        


    def full_var_bars(self): 
        ax = self.axes[self.ax_index] 
        self.ax_index += 1

        #self.colors = ['tab:cyan','tab:gray','tab:green','tab:pink','tab:orange']  
        #self.colors = ['y','tab:gray','tab:blue','tab:green','tab:green']  
        self.colors = ['y','tab:gray','cornflowerblue','tab:green','tab:green']  
        self.colors = ['yellow','tab:gray','cyan','lime','tab:green']  
        self.colors = ['yellow','tab:pink','cyan','lime','tab:green']  
        cands       = ['single', 'port', 'prior', 'combine', 'weighted'] 
        color_key   = {m: c for m,c in zip(cands, self.colors)}
       
        b_vals = [self.data[c].varexp.interval[-1] for c in self.names] 
        
        for i,n in enumerate(cands): 
            if n not in self.names: continue 
            br = self.data[n] 
            ht, [a,b] = br.varexp.val, br.varexp.interval

            if a < 0: a = 0 
            ax.plot([i,i],[a,b], color='k', zorder=1, linewidth=0.85)  
            ax.plot([i-0.1,i+0.1],[a,a], color='k', zorder=1, linewidth=0.85)  
            ax.plot([i-0.1,i+0.1],[b,b], color='k', zorder=1,linewidth=0.85)  
            if len(self.names) > 1: ax.text(i-0.2,b+0.001,'prs\n'+n, ha='center',rotation=20, fontsize=14) 
            

            if n != 'weighted':  
                ax.bar(i,ht, color = self.colors[i], edgecolor='k', alpha = 0.85, zorder = 2) 
            else:
                ax.bar(i,ht, color = 'white', edgecolor='k', alpha = 0.5, zorder = 1) 
                cuts, btm, FK = [], 0, [[c,br.varexp.frac[c], color_key[c]] for c in cands]
                for cand, frac, clr in FK: 
                    if cand == n or frac == 0: continue 
                    cuts.append([ht*frac, clr, cand]) 
                for piece, clr, cand in cuts: 
                    ax.bar(i, piece, bottom = btm , color = clr, alpha = 0.9, zorder = 2) 
                    btm += piece 

                    
        ax.plot([-0.6,4.6],[0,0], color='k', clip_on=False, zorder=0)   
        ax.plot([4.6,4.6], [0, b*1.5], color='k') 
        bm = round((round(max(b_vals) + (min(b_vals) / 5.0) ,2))/4.0,2) 
        for z in [bm, bm*2, bm*3, bm*4]: 
            ax.text(4.60, z, '-', ha='left') 
            ax.text(4.60, z, '-', ha='right') 
            ax.text(4.63, z, round(z,2), ha='left') 
        ax.axis('off') 
        
        ax.set_title('Variance Explained', fontweight='bold',fontsize=20, y = 0.9)  
        yMin, yMax = ax.get_ylim()
        
        ax.set_ylim(0 - (yMax/10.0), yMax) 
        
         





    def add_pred_scatter(self, method): 
        self.ax_index += 1 
        X = self.data[method].preds['prs'] 
        Y = self.data[method].preds['pheno'] 
        ax = self.axes[self.ax_index] 
        self.ax_index += 1 
        
        ax.scatter(X,Y, alpha=0.5) 
        self.get_lims(ax, BORDER=4, xLab = 'PRS\n('+method+')', yLab = 'Phenotype')
        ax.set_xlim(self.xMin - self.xSpan/20.0, self.xMax + self.xSpan/20.0) 
        ax.set_ylim(self.yMin - self.ySpan/10.0, self.yMax + self.ySpan/10.0) 
        ax.set_title('PRS Predictions', fontweight='bold',fontsize=20, y = 0.92)  



    def analyze_snp_dists(self, method, method2): 
        base_scores, base_key, self.DATA_KEY['LEN']['target'] = self.load_base_scores(self.data[method].SS['PREFIX'], self.data[method].SS['SUFFIX']) 
        snp_weights, model_scores  = self.data[method2].snp_weights, None
        self.draw_manhattan(base_scores, 'Target GWAS')
        
        if self.MODEL: 
            model_scores, model_key, self.DATA_KEY['LEN']['model'] = self.load_base_scores(self.model_key['SUMSTATS_PREFIX'],self.model_key['SUMSTATS_SUFFIX']) 
            w_bb =  self.merge_snp_scores(base_key, snp_weights, model_key) 
            w, b1, b2 = [x[0] for x in w_bb],[x[1] for x in w_bb],  [x[2] for x in w_bb] 
            self.DATA_KEY['LEN']['weight'] = len(w) 
            self.DATA_KEY['CORR'] = {'wTarget': S_CORR(w, b1), 'wModel': S_CORR(w, b2), 'MT': S_CORR(b1,b2)} 
        
        else: 

            w_bb =  self.merge_snp_scores(base_key, snp_weights, 'NA') 
            w, b1 = [x[0] for x in w_bb],[x[1] for x in w_bb]
            self.DATA_KEY['LEN']['weight'] = len(w) 
            self.DATA_KEY['CORR'] = {'wTarget': S_CORR(w, b1), 'wModel': 'NA', 'MT': 'NA'} 
            


        self.draw_manhattan(model_scores, 'Model GWAS') 
       

        #self.merge_snp_scores(base_key, snp_weights) 
        
        
        
            #w_bb =  self.merge_snp_scores(base_key, snp_weights, model_key) 
            
 

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
        ax.set_xlim(self.xMin, self.xMax) 
        ax.set_ylim(self.yMin, self.yMax) 
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
        
        if scores is None: 
            ax.axis('off') 
            return
        chromosomes = sorted(scores.keys())
        maxY, minY, xticks, xlabs, p_offset = 0, 100,  [] , [], 0 
        for ci,c in enumerate(chromosomes):
            bs = scores[c] 
            X,Y = [p_offset+b[0] for b in bs], [-log(b[2],10) for b in bs]
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
            p_offset = X[-1]
        
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
        
        ax.set_ylim(0 - (self.yMax/10.0), self.yMax*1.15) 
        return 




    def load_base_scores(self, f_prefix, suffix): 
        path   = "/".join(f_prefix.split('/')[0:-1])
        prefix = f_prefix.split('/')[-1] 
        c_names = [vars(self.args)[x] for x in ['ssf-snpid','ssf-p','ssf-beta']]
        full_len, full_scores, full_key    = 0, {}, {} 
        for f in os.listdir(path): 
            if prefix in f and suffix in f:  
                chr_key, chr_snps = {}, [] 
                chr_name = int(f.split(prefix)[1].split(suffix)[0]) 
                H, gf = zip_open(path+'/'+f)
                HK = {h: i for i,h in enumerate(H)} 
                if 'POS' in HK:   pk = HK['POS'] 
                elif 'pos' in HK: pk = HK['POS'] 
                else:             pk = None 
                SPB = [HK[cn] for cn in c_names] 
                for line in gf: 
                    line = line.split() 
                    snp, pv, beta = [float(line[k]) if i > 0 else line[k] for i,k in enumerate(SPB)]
                    if pk is None: pos = 0 
                    else:          pos = int(line[pk])  
                    chr_snps.append([pos, snp, pv, beta]) 
                    chr_key[snp] = [pv, beta, pos] 
            
            full_scores[chr_name] = sorted(chr_snps)  
            full_key[chr_name] = chr_key 
            full_len += len(chr_snps) 

        return full_scores,  full_key, full_len 

    
    
    def add_summary_table(self, method): 
        ax = self.axes[0] 
        dt = DataTable(self, ax).add_bridge_summary(method) 


    def add_logo(self, ax_index): 
        
        ax = self.axes[ax_index]
        if self.TYPE == 'MEGA': val = 0 
        elif self.TYPE == 'single': val = 1
        elif self.TYPE == 'port': val =  2 
        elif self.TYPE == 'prior': val = 3 
        else: val = random.randrange(5) 
        
        bp = BridgePic(ax).draw(val) 




class DataTable:
    def __init__(self, plot, ax,  clr='white'): 
        self.plot, self.args = plot, plot.args 
        self.ax, self.rows, self.cols, self.clr  = ax, [], [], clr 
        self.widths = [22,12,10,11,11,17,17,10,10,10]
        self.vwid = [0] + [sum(self.widths[0:i+1]) for i in range(len(self.widths))]   
        self.names, self.data  = self.plot.names, self.plot.data 
        

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



    def add_bridge_summary(self, method, fs1 = 20, fs2 = 14): 
        
        
        pd, K_CORR, K_LEN = self.plot, self.plot.DATA_KEY['CORR'], self.plot.DATA_KEY['LEN'] 
        
        
        tp = self.data[method] 
        mk = dd(lambda: 'NA') 
        if pd.MODEL: mk = self.plot.model_key
        
        p_lens = [len(open(pf).readlines())-1 for pf in tp.PT['FILES'].split(',')]
        p_name, p_type = tp.PT['FIELD-NAME'], tp.PT['TYPE'] 

        WD, yL, yJ = [36,32,32], 82, 23 
        self.add_row(['','Target','Model'], COLORS=['whitesmoke' for i in range(4)], X = (0,100), Y = (yL,100), FS = fs1, WIDTHS = WD)
        d1 = ['Populations:\n(GWAS, Geno/Pheno)', tp.ldpop+', '+tp.pop, mk['LDPOP']+', '+mk['POP']]
        d2 = ['GWAS Size\n(Rank Correlation: '+K_CORR['MT']+')', K_LEN['target'], K_LEN['model']]
        d3 = ['Result: SNP Weights\n('+str(K_LEN['weight'])+')','GWAS Corr= '+K_CORR['wTarget'], 'GWAS Corr='+K_CORR['wModel']]
        d4 = ['Geno/Pheno Samples \n('+p_name+', '+p_type+')', str(p_lens[0])+', '+str(p_lens[-1])+'\n (Test/Validate)', 'NA']
        
        
        for d in [d1,d2,d3,d4]: 
            self.add_row(d, COLORS=['white' for h in range(4)], X = (0,100), Y = (yL-yJ,yL), FS = fs2, WIDTHS = WD) 
            yL-=yJ
        self.ax.axis('off')
        return self 











class BridgePic: 
    def __init__(self, ax): 
        self.ax = ax


    

    def make_circle(self, h,k,r): 
        coord_list, H1, H2 = [], 1+h-r, h+r-1
        for x in range(H1,H2,1):
            y1 = 140 - math.sqrt(r**2 - (x-h)**2)
            coord_list.append([x, y1])
        X,Y = np.array([coord_list]).T
        X = [x[0] for x in X] 
        Y = [y[0] for y in Y] 
        return X,Y 

        

    
    def draw_posts(self,XY = [[0,0]], ym = 40): 
        for x,y in XY:
            R1 = Rect((x-2.5,y-ym),5,ym+107,facecolor='darkslategrey',edgecolor='k',zorder=6) 
            self.ax.add_patch(R1)
            self.ax.plot([x,x],[y,y+105],color='lightgray',zorder=7,linestyle='--')
            R2 = Rect((x,y-ym+2),5,2+ym+107,facecolor='darkslategrey',edgecolor='k',zorder=4) 
            self.ax.add_patch(R2)
            self.ax.plot([x+2.5,x+2.5],[y+2,2+y+105],color='lightgray',zorder=5,linestyle='--')
    




    def draw(self, VAL):
        #ax, LX,LY,RX,RY = self.ax, [],[],[],[]  
        LX,LY,RX,RY = [],[],[],[]  
        X,Y = self.make_circle(0,20,100)         
        for i,(x,y) in enumerate(zip(X,Y)): 
            rx,lx,ry,ly= x+200,x-200, y,y 
            SEND = i%10 == 0 
            if SEND: 
                self.ax.plot([x,x],[20,y],color='gray')
            if i < 100: 
                if i > 0:
                    sl = (100-i)/100.0
                    sp = sl + (1 - sl)/2.0
                    ry = y * sp 
                RX.append(rx) 
                RY.append(ry) 
                if SEND: self.ax.plot([rx,rx],[20,ry],color='gray') 
            else:
                sl = (i - 100) / 100.0 
                sp = sl + (1 - sl)/2.0
                ly = y * sp 
                LX.append(x-200) 
                LY.append(ly) 
                if SEND: self.ax.plot([lx,lx],[20,ly],color='gray') 

        self.ax.plot(RX,RY,color='k') 
        self.ax.plot(LX,LY,color='k') 
        self.ax.plot(X, Y,color='k')
        self.ax.plot([-200,200],[20,20],color='k') 
        self.draw_posts([[-100,20],[100,20]])
        R = Rect((-200,15),400,5,facecolor='gray',edgecolor='k',zorder=5) 
        self.ax.add_patch(R) 
        self.ALPHA=0.9 
        self.CMAP = self.select_map(VAL) 
        self.gbar1([-190],[20], bottom=-100) 
        self.CMAP = self.select_map(VAL) 
        self.gbar1([-190],[150], bottom=-100) 
        
        self.ax.set_ylim(-30,140)
        xMin, xMax = -179, 180 
        self.ax.plot([xMin,xMin],[-19,135],color='k',zorder=200,linewidth=2) 
        self.ax.plot([xMax,xMax],[-19,135],color='k',zorder=200,linewidth=2) 
        self.ax.plot([xMin,xMax],[135,135],color='k',linewidth=2) 
        self.ax.text(0,110,'BridgePRS',fontsize=26,ha='center',va='center',fontweight='bold') 
        self.ax.set_ylim(-21,135)
        self.ax.set_xlim(-179,180)
        self.draw_letters()
        
        xMin, xMax = self.ax.get_xlim() 
        yMin, yMax = self.ax.get_ylim() 
        
        self.ax.plot([xMin,xMax],[yMin, yMin], color='k', linewidth=3) 
        self.ax.plot([xMin,xMax],[yMax-1, yMax-1], color='k', linewidth=4,zorder=300) 
        
        self.ax.axis('off') 



    def select_map(self, VAL): 
        
        #if VAL == 0: return plt.cm.YlOrRd
        #if VAL == 1: return plt.cm.Blues
        MAP_CANDS = [plt.cm.copper, plt.cm.spring, plt.cm.autumn, plt.cm.cool, plt.cm.twilight, plt.cm.Wistia]
        MAP_CANDS = [plt.cm.Blues, plt.cm.Blues, plt.cm.YlOrRd, plt.cm.spring, plt.cm.autumn, plt.cm.cool, plt.cm.twilight, plt.cm.Wistia]
        MAP_CANDS = [plt.cm.Blues, plt.cm.Blues,  plt.cm.spring, plt.cm.autumn, plt.cm.cool, plt.cm.twilight, plt.cm.Wistia]
        return random.choice(MAP_CANDS) 
        
        
        

        
        
        


    # SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP #
    def gbar1(self,  x, y, width=400, bottom=20, VAL = 0):
        #X = [[.7, .7], [.6, .6]]
        X = [[.6, .6], [.8, .8]]
        #X = [[.3, .3], [.7, .7]]
        for left, top in zip(x, y):
            right = left + width
            self.ax.imshow(X, interpolation='bicubic', cmap=self.CMAP,extent=(left, right, bottom, top), alpha=self.ALPHA)
   








    def draw_letters(self): 
        xpt, yp = -174, 23
        k = 0
        letters = ['A','C','G','T','G','G','C','T'] 
        ypts    = [23,26,31,38,46,57,72,90] 
        while xpt < -100: 
            if k == 7: xpt -= 1.35 
            self.ax.text(xpt, ypts[k], letters[k], ha= 'center', va= 'center') 
            k+=1 
            xpt += 10
        letters = ['A','T','T','C','G','T','G','C','C','T'] 
        ypts    = [90,71,56,46,37,29,26,22.5,38,46,57,72,90] 
        k, xpt = 0, 107.75 
        while xpt < 180: 
            self.ax.text(xpt, ypts[k], letters[k], ha= 'center', va= 'center') 
            if k == 0: xpt += 8 
            else:      xpt += 10
            k+=1 










    def finish(self, suffix = 'NA'): 
        if self.ICO: 
            self.p_name = 'bridgePRS_ico.png' 
            plt.subplots_adjust(left=0.01, bottom=0, right=0.99, top=0.99,wspace=0.0,hspace=0.0) 
        else: 
            self.p_name = 'bridgePRS_LOGO.png' 
            plt.subplots_adjust(left=0, bottom=0, right=1.0, top=1.0,wspace=0.0,hspace=0.0) 
        plt.savefig(self.p_name,dpi=400) 
        plt.clf() 
        plt.close() 
        return









    
    
    
#bp = BridgePlot(sys.argv[-1]) 
    


#bp.draw()
#bp.finish() 

