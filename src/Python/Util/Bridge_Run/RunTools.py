import sys, os 
import gzip
import math
import random
from collections import defaultdict as dd
try:                                                                                                                                                                                                                        
    import matplotlib                                                                                                                                                                                                       
    import matplotlib.pyplot as plt                                                                                                                                                                                         
    from matplotlib.patches import Rectangle as Rect                                                                                                                                                                        
except: pass        




def bridge_error(eString):
    if type(eString) in [list,tuple]:
        sys.stderr.write('\n     BridgeToolBoxError: '+eString[0]+'\n')                     
        for es in eString[1::]: sys.stderr.write('                         '+es+'\n')
    else: sys.stderr.write('\n    BridgeToolBoxError: '+eString+'\n')                           
    sys.exit()

                                                                                                                
def bridge_warning(eString):                                                                           
    if type(eString) in [list,tuple]:                                                                           
        sys.stderr.write('\n     BridgeToolBoxWarning: '+eString[0]+'\n')                                           
        for es in eString[1::]: sys.stderr.write('                           '+es+'\n')                                   
    else: sys.stderr.write('\n     BridgeToolBoxWarning: '+eString+'\n')                                            


def zip_open(fp, HEADER = False):                                                                               
    if fp.split('.')[-1] == 'gz': gf = gzip.open(fp, 'rt')                                                      
    else:                         gf = open(fp, 'rt')                                                           
    if HEADER:                                                                                                  
        HL = gf.readline().split()                                                                              
        gf.close()                                                                                              
        return HL                                                                                               
    else: return gf  



def zip_open(fp, HEADER = True):                                                                                                                                                                                                                                                       
    if fp.split('.')[-1] == 'gz': gf = gzip.open(fp, 'rt')                                                                                                                                                                                                                              
    else:                         gf = open(fp, 'rt')                                                                                                                                                                                                                                   
    if HEADER:
        HL = gf.readline().split()                                                                                                                                                                                                                                                      
        return HL, gf 
    return gf 


def R_CORR(x,y): 
    Rc = round(my_pearson_corr(x,y),2) 
    return str(Rc) 
        
def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

def my_pearson_corr(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / math.sqrt(xdiff2 * ydiff2)







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



    def add_summary_header(self, method, fs1 = 20, fs2 = 14): 

        pd, K_CORR, K_LEN = self.plot, self.plot.DATA_KEY['CORR'], self.plot.DATA_KEY['LEN'] 
        tp = self.data[method] 
        mk = dd(lambda: 'NA') 

        WD, yL, yJ = [36,32,32], 82, 23 
        self.add_row(['Target','Model'], COLORS=['whitesmoke','whitesmoke','whitesmoke'], X = (0,100), Y = (0,100), FS = fs1, WIDTHS = [50,50]) 
        self.ax.axis('off')
        self.ax.set_ylim(0,100)  
        self.ax.plot([-50,100],[99.2,99.2],clip_on=False, color='k', zorder=1, lw=2) 
        R = Rect((-1.07,-2),1.1,3,facecolor='k',edgecolor='k',zorder=5, clip_on=False) 
        self.ax.add_patch(R) 
        self.ax.set_xlim(0,100) 

        return self    







    def add_smart_summary(self, method, fs1 = 20, fs2 = 14): 
        
        
        pd, K_CORR, K_LEN = self.plot, self.plot.DATA_KEY['CORR'], self.plot.DATA_KEY['LEN'] 
        
        

        tp = self.data[method] 
        p_data, p_res = pd.data[method], pd.results[0] 
        mk = dd(lambda: 'NA') 
        if pd.MODEL: mk = self.plot.model_key
        
        my_phenos =  p_data.preds['pheno'] 
        my_prs =  p_data.preds['prs'] 
        my_weights = p_data.snp_weights.values() 
        big_data = {'pheno': p_data.preds['pheno'], 'prs': p_data.preds['prs'], 'weights': p_data.snp_weights.values()} 
        big_stats = {} 
        for b,V in big_data.items(): 
            my_avg = average(V) 
            my_var = sum([(v-my_avg)*(v-my_avg) for v in V])/len(V) 
            my_std = my_var ** 0.5 
            big_stats[b] = [round(my_avg,2), round(my_std,2), round(min(V),1), round(max(V),1)] 
            
        p_lens = [len(open(pf).readlines())-1 for pf in tp.PT['FILES'].split(',')]
        d1 = ['Populations:\n(Reference, GWAS)', tp.ldpop+', '+tp.pop, mk['LDPOP']+', '+mk['POP']]
        target_gwas_sample_size, base_gwas_sample_size = p_data.SS['SIZE'] , mk['SUMSTATS_SIZE'] 
        d2 = ['GWAS Sample Size:', target_gwas_sample_size, base_gwas_sample_size] 
        d3 = ['GWAS Total SNPs:\n(Correlation: '+K_CORR['MT']+')', K_LEN['target'], K_LEN['model']]
        d4 = ['Phenotype:\n(type)', tp.PT['FIELD-NAME']+'\n('+tp.PT['TYPE']+')','NA']
        if tp.PT['TYPE'] == 'continuous': 
            my_res = ', '.join([str(xx) for xx in big_stats['pheno']]) 
            d5 = ['Phenotype Dist:\n(Mean,Std,Min,Max)', my_res, 'NA'] 
        else: 
            treatment = int(sum(big_data['pheno']))
            controls = int(len(big_data['pheno']) - treatment) 
            d5 = ['Phenotype Count:\n(Case vs Controls)', str(treatment)+', '+str(controls), 'NA'] 

        my_res = ', '.join([str(xx) for xx in big_stats['prs']]) 
        d6 = ['PRS Dist:\n(Mean,Std,Min,Max)', my_res, 'NA'] 
        d7 = ['Result: SNP Weights\n('+str(K_LEN['weight'])+')','GWAS Corr= '+K_CORR['wTarget'], 'GWAS Corr='+K_CORR['wModel']]
        WD, yL, yJ = [33.4,33.3,33.3], 100, 14.3
        for d in [d1,d2,d3,d4,d5,d6,d7]: 
            self.add_row(d, COLORS=['white' for h in range(4)], X = (0,100), Y = (yL-yJ,yL), FS = fs2, WIDTHS = WD) 
            yL-=yJ
        self.ax.axis('off')
        self.ax.set_xlim(0,100) 
        self.ax.set_ylim(0,100) 
        clr = 'k' 
        self.ax.plot([0,100],[100.4,100.4],clip_on=False, color=clr, zorder=1, lw=2) 
        self.ax.plot([0,100],[0,0],clip_on=False, color=clr, zorder=1, lw=2) 
        self.ax.plot([0,0],[0,119.5],clip_on=False, color=clr, zorder=1, lw=2) 
        self.ax.plot([100,100],[0,119.5],clip_on=False, color=clr, zorder=1, lw=2) 
        return self 























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
        d2 = ['GWAS Size\n(Correlation: '+K_CORR['MT']+')', K_LEN['target'], K_LEN['model']]
        d3 = ['Result: SNP Weights\n('+str(K_LEN['weight'])+')','GWAS Corr= '+K_CORR['wTarget'], 'GWAS Corr='+K_CORR['wModel']]
        d4 = ['Geno/Pheno Samples \n('+p_name+', '+p_type+')', str(p_lens[0])+', '+str(p_lens[-1])+'\n (Test/Validate)', 'NA']
        
        
        for d in [d1,d2,d3,d4]: 
            self.add_row(d, COLORS=['white' for h in range(4)], X = (0,100), Y = (yL-yJ,yL), FS = fs2, WIDTHS = WD) 
            yL-=yJ
        self.ax.axis('off')
       
        self.ax.set_ylim(100,-100) 

        return self 
























class BridgePic: 
    def __init__(self, ax): 
        self.ax = ax


    def make_circle(self, h,k,r): 
        coord_list, H1, H2 = [], 1+h-r, h+r-1
        for x in range(H1,H2,1):
            y1 = 140 - math.sqrt(r**2 - (x-h)**2)
            coord_list.append([x, y1])
       
        X = [float(c[0]) for c in coord_list] 
        Y = [c[1] for c in coord_list] 
        return X,Y 

    
    def draw_posts(self,XY = [[0,0]], ym = 40): 
        for x,y in XY:
            R1 = Rect((x-2.5,y-ym),5,ym+107,facecolor='darkslategrey',edgecolor='k',zorder=6) 
            self.ax.add_patch(R1)
            self.ax.plot([x,x],[y,y+105],color='lightgray',zorder=7,linestyle='--')
            R2 = Rect((x,y-ym+2),5,2+ym+107,facecolor='darkslategrey',edgecolor='k',zorder=4) 
            self.ax.add_patch(R2)
            self.ax.plot([x+2.5,x+2.5],[y+2,2+y+105],color='lightgray',zorder=5,linestyle='--')
    




    def draw(self, TYPE): 
        
        if TYPE == 'MEGA': VAL = 0
        elif TYPE == 'single': VAL = 1 
        else:                  VAL = random.randrange(5) 


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
        #self.ax.plot([xMin,xMin],[-19,135],color='k',zorder=200,linewidth=2) 
        #self.ax.plot([xMax,xMax],[-19,135],color='k',zorder=200,linewidth=2) 
        #self.ax.plot([xMin,xMax],[135,135],color='k',linewidth=2) 
        self.ax.text(0,110,'BridgePRS',fontsize=13,ha='center',va='center',fontweight='bold') 
        self.ax.set_ylim(-21,135)
        self.ax.set_xlim(-179,180)
        
        self.draw_letters()
        xMin, xMax = self.ax.get_xlim() 
        yMin, yMax = self.ax.get_ylim() 
        #self.ax.plot([xMin,xMax],[yMin, yMin], color='k', linewidth=3) 
        #self.ax.plot([xMin,xMax],[yMax-1, yMax-1], color='k', linewidth=4,zorder=300) 
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
   








    def draw_letters(self, fs=5): 
        xpt, yp = -174, 23
        k = 0
        letters = ['A','C','G','T','G','G','C','T'] 
        ypts    = [23,26,31,38,46,57,72,90] 
        while xpt < -100: 
            if k == 7: xpt -= 1.35 
            self.ax.text(xpt, ypts[k], letters[k], ha= 'center', va= 'center',fontsize=fs) 
            k+=1 
            xpt += 10
        letters = ['A','T','T','C','G','T','G','C','C','T'] 
        ypts    = [90,71,56,46,37,29,26,22.5,38,46,57,72,90] 
        k, xpt = 0, 107.75 
        while xpt < 180: 
            self.ax.text(xpt, ypts[k], letters[k], ha= 'center', va= 'center', fontsize=fs) 
            if k == 0: xpt += 8 
            else:      xpt += 10
            k+=1 













    
    
