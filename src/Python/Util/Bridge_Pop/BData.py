import sys, os 
from collections import defaultdict as dd
from collections import Counter as cc 
from . import PopTools as ptools



class BData: 
    def __init__(self, pop):
        self.parent, self.mps, self.progress = pop, pop.mps, pop.progress 
        self.ld_pop, self.ld_path = self.parent.ref_pop, self.parent.ld_path
        self.X_fields = ['--clump-field',self.parent.field_key['P'], '--clump-snp-field',self.parent.field_key['ID']] 
        self.map = dd(list) 
        self.VALID, self.BYCHR = False, False 


        
    def load_panel(self): 
        self.VALID, self.BYCHR = True, True
        self.id_file, self.prefix, self.map = self.get_ld_data(self.ld_pop, self.ld_path) 
        return self 


    def get_ld_data(self, ldpop, ld_path): 
        ref_panel = ld_path.split('/')[-1]                                                                                                                                                                                  
        key, cands, id_files = {}, [], []                                                                                                                                                                                   
        for f in os.listdir(ld_path):                                                                                                                                                                                       
            fp = f.split('.')[0].upper().split('_')                                                                                                                                                                         
            if f[0] in ['.','_']: continue                                                                                                                                                                                  
            if f.split('.')[-1] in ['bed','bim','fam']: cands.append(".".join(f.split('.')[0:-1]))                                                                                                                          
            elif 'IDS' in fp and ldpop.upper() in fp: id_files.append(f)                                                                                                                                                    
            else: continue                                                                                                                                                                                                  
        cands = [cn for cn in cc(cands) if cc(cands)[cn] == 3]                                                                                                                                                              
        if    len(cands) == 0: self.mps.error('LDError: No plink triples (bed, bim, fam) found in LD-Reference Path: '+ld_path)                                                                                                     
        elif  len(cands) == 1: self.mps('LD-Reference files must be split by chromosome: '+ld_path)                                                                                                                  
        elif  len(id_files) == 0: self.mps('No Corresponding ID file found for population '+ldpop+' in: '+ld_path)                                                                                                    
        elif  len(id_files) > 1: self.mps('Multiple ID files found for pop '+ldpop+' in: '+ld_path+'  ('+",".join(id_files)+')')                                                                                      
        BYCHR = True          
        my_prefix, my_suffix = ptools.get_prefix_suffix(cands)                                                                                                                                                                
        for cand in cands:                                                                                                                                                                                                  
            chr_cand = cand.split(my_prefix)[-1].split(my_suffix)[0]                                                                                                                                                        
            key[chr_cand] = ld_path+'/'+cand                                                                                                                                                                                
        return ld_path+'/'+id_files[0], ld_path+'/'+my_prefix, key 









