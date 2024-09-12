import sys, os, gzip, shutil
from collections import defaultdict as dd
from collections import Counter as cc 
#from .BridgeProgress  import BridgeProgress
#import BridgeTools 





class GenoPheno: 
    def __init__(self, pop, phenotype): 

        self.parent, self.mps, self.progress, self.paths, self.pop_name, self.pop_type = pop, pop.mps, pop.progress, pop.paths, pop.name, pop.type 
        self.genotype_prefix, self.phenotype_file, self.validation_file = self.parent.genotype_prefix, self.parent.phenotype_file, self.parent.validation_file
        self.phenotype, self.covariates = phenotype, self.parent.covariates 
        self.TESTS, self.VALID, self.BYCHR = dd(bool), False, False 
        self.fields, self.map = {}, dd(list) 
        


    def load(self): 
        self.ph_get_genotypes() 
        if self.genotype_prefix and self.phenotype_file: self.ph_fill(self.phenotype_file, self.validation_file)
        return self 

    
    def ph_get_genotypes(self): 
        CK, g_path, g_prefix = dd(int), "/".join(self.genotype_prefix.split('/')[0:-1]), self.genotype_prefix.split('/')[-1] 
        for f in [x for x in os.listdir(g_path) if x[0:len(g_prefix)] == g_prefix]:
            if f.split('.')[-1] in ['bed','bim','fam']: CK[".".join(f.split('.')[0:-1])]+=1                                                                                                                                                                                   
        cands = [k for k,v in CK.items() if v == 3]
        if len(cands) == 0: self.mps.error('Invalid genotype data prefix '+self.genotype_prefix) 
        if len(cands) == 1:
            self.genotype_prefix = "/".join(self.genotype_prefix.split('/')[0:-1])+'/'+cands[0]
            self.BYCHR = False 
        else: 
            self.BYCHR, k  = True, 0 
            cands = list(set(cands))                                                                                                                                                                                                                                                        
            while True:
                my_prefix = list(set([c[0:k] for c in cands]))[0] 
                ck1 = list(set([c[0:k+1] for c in cands]))                                                                                                                                                                                                                                  
                if len(ck1) == 1: k+=1                                                                                                                                                                                                                                                      
                else: break
        return


    def split_phenotype_file(self, savepath, fname):
        self.TESTS['NOVALID'] = True 
        with open(fname, "r") as f: 
            header = f.readline().strip() 
            data   = [ln.strip() for ln in f]
            h_len = int(len(data)/2)
        d1,d2 = data[0:h_len], data[h_len::]
        test_file, valid_file = savepath+'/'+self.pop_name+'.test_phenos.dat', savepath+'/'+self.pop_name+'.valid_phenos.dat'
        w1, w2 = open(test_file,'w'), open(valid_file,'w') 
        w1.write(header+'\n') 
        w1.write("\n".join(d1)+'\n') 
        w2.write(header+'\n') 
        w2.write("\n".join(d2)+'\n') 
        w1.close() 
        w2.close() 
        return [test_file, valid_file] 


    def ph_fill(self, phenotype_file, validation_file = None):
        self.VALID, self.type, col_data = True, 'binary', dd(list) 
        if self.pop_type.upper() == 'TARGET': 
            if validation_file is not None:     self.files = [phenotype_file, validation_file]
            else:                               self.files = self.split_phenotype_file(self.parent.paths['save'],phenotype_file) 
            self.X_fields = ['--test.data',self.files[0],'--valid.data',self.files[1]]
        else:
            self.files = [phenotype_file] 
            self.X_fields = ['--test.data',self.files[0],'--valid.data','0'] 
        for i,fn in enumerate(self.files): 
            with open(fn, 'rt') as f: 
                lp = f.readline().split() 
                if i == 0: self.header, lz = [x for x in lp] , ",".join(lp) 
                elif ",".join(lp) != lz: self.mps.error('Phenotype File Headers Do Not Match: '+lz+' AND '+",".join(lp))
                for k,line in enumerate(f):
                    line = line.split() 
                    for j,c in enumerate(self.header): col_data[c].append(line[j]) 
                    if k > 100: break 

        pheno_cands = [h for h in self.header[1::] if 'ID' not in h]
        if self.covariates is not None: 
            self.fields['COVARIATES'] = self.covariates
            self.X_fields.extend(['--cov.names',self.covariates]) 
            for c in self.covariates.split(','): 
                if c not in self.header:  self.mps.error('Invalid covariate field name(s) supplied '+self.covariates+', Available Fields: '+','.join(self.header)) 
                else:                     pheno_cands.pop(pheno_cands.index(c)) 

        self.availible_fields = ','.join(self.header) 
        #if self.phenotype is None: self.mps.error('Phenotype field required (--phenotype)\n Available Fields: '+','.join(self.header)) 
        if self.phenotype is not None: 
            self.fields['NAME'] = self.phenotype 
            self.X_fields.extend(['--pheno.name',self.phenotype]) 
            if self.phenotype not in self.header: self.mps.error('Invalid phenotype field name(s) supplied '+self.phenotype+', Available Fields: '+','.join(self.header)) 
            if len((list(set(col_data[self.phenotype])))) > 2:   self.type = 'continuous' 
            else:                                                self.X_fields.extend(['--binary','1']) 
        return self 
        
        

