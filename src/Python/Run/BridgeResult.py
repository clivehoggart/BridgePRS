#!/usr/bin/env python3
import os, sys
from collections import defaultdict as dd
##########################################################################################################################################
class BridgeResult:
    def __init__(self, name = None, pop= None, SYN = False):
        self.name, self.pop, self.SYNTHETIC   = name, pop, SYN

    def read_file(self, res): 
        f, X, self.output, self.paths, self.SS, self.PT = open(res), dd(bool), dd(bool), dd(bool), dd(bool), dd(bool) 
        for line in f: 
            a,b = line.strip().split('=') 
            a_init, a_tail = a.split('_')[0], a.split('_')[-1] 
            if len(a.split('_')) == 3 and a_tail == 'PREFIX': self.paths[a.split('_')[1]] = b 
            elif a_init in ['GENOTYPE','PHENOTYPE']:          self.PT[a_tail] = b 
            elif a_init in ['SUMSTATS','SNP']:                self.SS[a_tail] = b 
            else:                                             X[a] = b  
        f.close() 
        self.pop, self.ldpop, self.name, self.modelpath, self.ldpath = X['POP'], X['LDPOP'], X['MODULE_NAME'], X['MODEL_FILE'], X['LDPATH']
        q_path, qn =  "/".join(self.paths['QUANTIFY'].split('/')[0:-1]), self.paths['QUANTIFY'].split('/')[-1] 
        for f in os.listdir(q_path): 
            if qn in f and f.split('.')[-1] not in ['log','png','pdf'] and '_lambda_' not in f and '_alpha_' not in f and '_tau_' not in f:  
                if 'var_explained' in f:    self.varexp     = BridgeOutput('VAR').read(q_path+'/'+f).process() 
                elif 'preds' in f:          self.preds      = BridgeOutput('PRED').read(q_path+'/'+f).process()            
                elif 'snp_weights' in f:    self.snp_weights    = BridgeOutput('WEIGHT').read(q_path+'/'+f).process()                    
                continue 
        return self
        
    def read_combo(self, combo): 
        #cands = ['prs.Stages1+2','prs.weighted'] 
        cands = ['Stage1+2','Weighted'] 
        for k,f in combo.key.items(): 
            if 'var_explained' in f:    varexp         =   BridgeOutput('VAR').read(f).process(cands) 
            elif 'preds' in f:          preds          =   BridgeOutput('PRED').read(f).process(cands)            
            elif 'snp_weights' in f:    snp_weights    =   BridgeOutput('WEIGHT').read(f).process(cands)                    
            else: continue  
        return [BridgeResult('prs-combine', combo.pop, True).add_triple([varexp[0], preds[0], snp_weights[0]]), BridgeResult('prs-weighted', combo.pop, True).add_triple([varexp[1], preds[1], snp_weights[1]])]  

    def add_triple(self, triple_data): 
        self.varexp, self.preds, self.snp_weights = triple_data 
        return self 






class BridgeOutput:
    def __init__(self, RULE):
        self.RULE = RULE
    
    def read(self, F):  
        f = open(F,'rt') 
        if self.RULE == 'VAR':
            rows = [[x.strip('\"').strip('\%') for x in lp.strip().split(',')] for lp in f.readlines()] 
            labels = rows[0][1::] 
            self.K  = {row[0]: {labels[i]: float(x) for i,x in enumerate(row[1::])} for row in rows[1::]} 
        else:
            self.K, header = dd(list), ['sample_id' if x == '---' else x for x in f.readline().split()]
            hk = [True if i > 0 and h.split('.')[-1] != 'allele' else False for i,h  in enumerate(header)]
            for line in f: 
                line = [float(x) if hk[i] else x for i,x in enumerate(line.split())]  
                if self.RULE == 'WEIGHT': 
                    self.K[line[0]] = line[-1] 
                elif self.RULE == 'PRED':
                    for i,x in enumerate(line): self.K[header[i]].append(x) 
        f.close() 
        return self 
    
    def create_var_key(self, V):
        
        try: self.frac, self.val, self.interval, self.devs =  V['Frac'], V['Est'], [V['2.5'],V['97.5']], [V['cv.dev'], V['cv.dev.sd']]
        except KeyError:  self.frac, self.val, self.interval, self.devs =  1, V['Est'], [V['2.5'],V['97.5']], [V['cv.dev'], V['cv.dev.sd']]
        return self 

    def process(self, cands =[]): 
        if len(cands) == 0: 
            if self.RULE != 'VAR': return self.K 
            else:                  return self.create_var_key(self.K[[z for z in self.K.keys()][0]])
        elif self.RULE == 'PRED':
            cand_dicts = [dd(list) for c in cands]  
            for k,V in self.K.items(): 
                for i in range(len(cand_dicts)): 
                    if k == 'pheno':  cand_dicts[i][k] = V 
                    elif k == cands[i]:  cand_dicts[i]['prs'] = V 
                    
                    #if k.split('.')[0] != 'prs': cand_dicts[i][k] = V 
                    #elif k == cands[i]:          cand_dicts[i]['prs'] = V 
            return cand_dicts 
        elif self.RULE == 'VAR': 
            MK, FK = {'Stage1': 'single', 'Stage2': 'prior',  'Stage1+2': 'combine', 'Weighted': 'weighted'}, dd(float)  
            for k in self.K.keys(): FK[MK[k]] = self.K[k]['Prob']
            VK = [self.K[c] for c in cands] 
            VK[-1]['Frac'] = FK 
    
                
        
            # ' combined' 

            #MK = {'prs.Stages1+2': 'Stage1+2', 'prs.weighted': 'Weighted'} 
            #VK = [self.K[MK[c]] for c in cands] 
            #VK = [self.K[c] for c in cands]
            #print(cands)
            #print(self.K.keys()) 

            #sys.exit() 


            
            
            return [BridgeOutput(self.RULE).create_var_key(vk) for vk in VK] 
        else: return([None, self.K]) 






















