import sys, os, gzip, shutil
from collections import defaultdict as dd
from collections import Counter as cc 
#from .BridgeProgress  import BridgeProgress
#import BridgeTools 



#from .BridgePop      import BridgePop

# model_file
REV_COMP =  {'A': 'T', 'C': 'G','G': 'C', 'T': 'A'}
NUM_STRS =  ['1','2','3','4','5','6','7','8','9','0'] 


def bridge_debug_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('BridgeDebugError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('BridgeDebugError: '+eString+'\n')
    sys.exit(2) 


# Ambiguous debug_level Ref/Alt Missing snps.txt check printout

def bridge_pop_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgePopError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('\nBridgePopError: '+eString+'\n')
    sys.exit(2) 


def bridge_sumstats_error(eString): 
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeSumstatsError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('\nBridgeSumstatsError: '+eString+'\n')
    sys.exit(2) 


def bridge_sumstats_warning(eString): 
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeSumstatsWarning: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('\nBridgeSumstatsWarning: '+eString+'\n')




            
def zip_open(fp, HEADER = False): 
    if fp.split('.')[-1] == 'gz': gf = gzip.open(fp, 'rt') 
    else:                         gf = open(fp, 'rt')  
    if HEADER: 
        HL = gf.readline().split() 
        gf.close() 
        return HL 
    else: return gf 
        

















def bridge_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeSettingsError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                     '+es+'\n')
    else: sys.stderr.write('\nBridgeSettingsError: '+eString+'\n')
    sys.exit(2) 

def bridge_ld_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeLDPanelError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                     '+es+'\n')
    else: sys.stderr.write('\nBridgeLDPanelError: '+eString+'\n')
    sys.exit(2) 


def get_prefix_suffix(cands):                                                                                                                                                                                     
    my_prefix,my_suffix, k = ' ', ' ', 0                                                                                                                                                                                
    while True:                                                                                                                                                                                                         
        my_prefixes = list(set([c[0:k] for c in cands]))                                                                                                                                                                
        if len(my_prefixes) == 1:                                                                                                                                                                                       
            my_prefix = my_prefixes[0]                                                                                                                                                                                  
            k+=1                                                                                                                                                                                                        
        else:                                                                                                                                                                                                           
            break                                                                                                                                                                                                       
    k = 0                                                                                                                                                                                                               
    while True:                                                                                                                                                                                                         
        my_suffixes = list(set([c[len(c)-k::] for c in cands]))                                                                                                                                                         
        if len(my_suffixes) == 1:                                                                                                                                                                                       
            my_suffix = my_suffixes[0]                                                                                                                                                                                  
            k+=1                                                                                                                                                                                                        
        else:                                                                                                                                                                                                           
            break                                                                                                                                                                                                       
    
    if len(my_prefix) == 0: my_prefix = ' '                                                                                                                                                                             
    if len(my_suffix) == 0: my_suffix = ' '  
    return my_prefix, my_suffix    







#POP_ARGS = ['pop','ldpop','sumstats_file','sumstats_prefix','sumstats_suffix','genotype_prefix','phenotype_file']
#POP_ARGS = ['pop','ldpop','ld_path','ldref','sumstats_file','sumstats_prefix','sumstats_suffix','genotype_prefix','phenotype_file']
#POP_ARGS = ['pop','ldpop','ld_path','sumstats_prefix','sumstats_suffix','genotype_prefix','phenotype_file']
#DEFAULT = dd(lambda: None) 
#for a,b in [['fst',0.15]]: DEFAULT[a] = b 
#for a,b in [["p","P"],["snpid","ID"],["se","SE"],["n","OBS_CT"],["beta","BETA"],["ref","REF"],["alt","A1"],["maf","A1_FREQ"]]: DEFAULT['ssf-'+a] = b 


#DEFAULT = dd(lambda: None)                                                                                                                                                                                                  
#    BridgePop.py:        self.X_fields = ['--clump-field',vars(self.args)['ssf-p'], '--clump-snp-field',vars(self.args)['ssf-snpid']] 
#    BridgePop.py:        sum_stats_fields    =  ['ssf-alt', 'ssf-beta', 'ssf-maf', 'ssf-p', 'ssf-ref', 'ssf-se', 'ssf-snpid', 'ssf-n']
#    BridgeSettings.py:for a,b in [["p","P"],["snpid","ID"],["se","SE"],["n","OBS_CT"],["beta","BETA"],["ref","REF"],["alt","A1"],["maf","A1_FREQ"]]: DEFAULT['ssf-'+a] = b 
#BridgePop.py:            self.X_fields.extend(['--pheno.name',self.args.phenotype]) 
#BridgePop.py:            else:                                                     self.X_fields.extend(['--binary','1']) 
#BridgePop.py:            self.X_fields.extend(['--cov.names',self.args.covariates]) 
#  $cov.names
# [1] "000"
#$pheno.name
#[1] "y.binary"
#$ranking
#[1] "pv"
#$by.chr
#[1] 1
#$strand.check
#[1] 1
#$binary
#[1] 1
        #BridgePop.py:        if self.args.thinned_snp_file is None: self.thin_snps = '0'
        #BridgePop.py:        else:                                  self.thin_snps = self.args.thinned_snp_file
        # BridgeBase.py:        ### DO AN EXAMPLE WITH --thinned snp list
        # BridgeBase.py:        X.extend(['--S','0,0.25,0.5,0.75,1','--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',self.ss.thin_snps])
        #BridgeBase.py:        #alt_opts = ['S', 'n.max.locus', 'thinned.snplist']
        #    BridgeBase.py:        #X.extend(['--S','0,0.25,0.5,0.75,1','--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',str(self.args.thinned_snplist)])
        #        BridgeBase.py:        #X.extend(['--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',str(self.args.thinned_snplist)])
        #            BridgeBase.py:        X.extend(['--n.max.locus',str(self.args.max_clump_size),'--thinned.snplist',self.ss.thin_snps]) 
        #
        # 

#$help
#[1] FALSE
#
#       bridge_error('Universal Argument '+k+' Cannot Be Included In Configuration File (must be passed on command line --'+k+')')
#
#
#POP_ARGS = ['pop','ldpop','ld_path','sumstats_prefix','sumstats_suffix','genotype_prefix','phenotype_file']
#DEFAULT = dd(lambda: None) 
#for a,b in [['fst',0.15]]: DEFAULT[a] = b 







#DEFAULT = dd(lambda: None)                                                                                                                                                                                                  
#for a,b in [["p","P"],["snpid","ID"],["se","SE"],["n","OBS_CT"],["beta","BETA"],["ref","REF"],["alt","A1"],["maf","A1_FREQ"]]:                                                                                              


#POP_DOUBLE   =    ['pop','ldpop','ld_path','sumstats_prefix','sumstats_suffix','genotype_prefix','phenotype_file','validation_file','snp_file']                                                                                                                                             
#POP_SINGLE   =    ['genotype_prefix','phenotype_file','validation_file','snp_file','model_file','phenotype','covariates','thinned_snp_file','max_clump_size']                                                         





#for a,b in [["p","P"],["snpid","ID"],["se","SE"],["n","OBS_CT"],["beta","BETA"],["ref","REF"],["alt","A1"],["maf","A1_FREQ"]]: DEFAULT['ssf-'+a] = [b] 



#POP_REQ = ['pop','ldpop','ld_path','sumstats_prefix','sumstats_suffix','genotype_prefix','phenotype_file','validation_file','snp_file','covariates','thinned_snp_file','max_clump_size']

#if self.args.phenotype is not None:                                                                                                                                                                                 
#    self.fields['NAME'] = self.args.phenotype                                                                                                                                                                       
#    self.X_fields.extend(['--pheno.name',self.args.phenotype])                                                                                                                                                      
#    if self.args.phenotype not in self.header: bridge_pop_error('Invalid phenotype field name(s) supplied '+self.args.phenotype+', Available Fields: '+','.join(self.header))                                       
#        if len((list(set(col_data[self.args.phenotype])))) > 2:   self.type = 'continuous'                                                                                                                              
#        else:                                                     self.X_fields.extend(['--binary','1'])                                                                                                                

#if self.args.covariates is not None:                                                                                                                                                                                
#    self.fields['COVARIATES'] = self.args.covariates                                                                                                                                                                
#    self.X_fields.extend(['--cov.names',self.args.covariates])                                                                                                                                                      
#    for c in self.args.covariates.split(','):                                                                                                                                                                       
#        if c not in self.header:  bridge_pop_error('Invalid covariate field name(s) supplied'+self.args.covariates+', Available Fields: '+','.join(self.header))                                                   
#    return self                                                                                                                                                                                                         
                                                                                                                                                                                                                                                                                           

POP_MUST     = ['pop','ldpop','ld_path','sumstats_prefix','sumstats_suffix', 'genotype_prefix','phenotype_file'] 
POP_OPTS     = ['validation_file','snp_file','covariates','thinned_snp_file','max_clump_size'] 
#SS_FIELDS    = [["p","P"],["snpid","ID"],["se","SE"],["n","OBS_CT"],["beta","BETA"],["ref","REF"],["alt","A1"],["maf","A1_FREQ"]] 
SS_FIELDS    = [["p","P"],["snpid","ID"],["beta","BETA"],["ref","REF"],["alt","A1"]] 
POP_DEFAULTS = {'fst': [0.1], 'phenotype': [None]} 
for a,b in SS_FIELDS: POP_DEFAULTS['ssf-'+a] = [b] 
for k in POP_MUST + POP_OPTS:    POP_DEFAULTS[k] = [None]



#X = POP_MUST + POP_OPTS + [k for k in POP_DEFAULTS.keys()]
#POP_OPTS = ['pop','ldpop','ld_path','sumstats_prefix','sumstats_suffix','genotype_prefix','phenotype_file','validation_file','snp_file','covariates','thinned_snp_file','max_clump_size']







    
class BridgePop: 
    def __init__(self, P, args, paths, pop_type, prevPop = None, popDict=None): 
        
        print(P) 

        self.name, self.ref_pop, self.file, self.type  = P['pop'], P['ldpop'], P['config_name'], pop_type
        
        
        print('whhho') 
        
        self.bdata = BData(self, P)
        
        print(self.bdata.VALID) 

        self.genopheno = GenoPheno(self, P) 
        self.sumstats = SumStats(self, P, self.genopheno, prevPop)  


        sys.exit() 

        #for p in ['pop','ld_path','sumstats_prefix']:
        # sumstats_suffix','sumstats_size','sumstats_fields']'

        #LEVEL_1 = ['config_name', 'pop', 'ldpop','ld_path','sumstats_prefix','sumstats_suffix','sumstats_size','sumstats_fields']                                                                                           
        #LEVEL_2 = ['genotype_prefix','phenotype_file','validation_file']                                                                                                                                                    
        #LEVEL_3 = ['validation_file','snp_file','max_clump_size','thinned_snp_file']                                                                                                                                        
        #LEVEL_4 = ['phenotype'] 


        #sys.exit() 




class BridgePop2: 
    def __init__(self, P, args, paths, pop_type, prevPop = None, popDict=None): 
        
        self.name = P['pop'] 
        self.ref_pop = P['ldpop']
        self.args, self.paths, self.type, self.prevPop = args, paths, pop_type, prevPop 
        


        self.bdata = BData(self, P)
        self.genopheno = GenoPheno(self, P) 
        self.sumstats = SumStats(self, P, self.genopheno, prevPop)  
        
        self.valids = [self.bdata.VALID, self.sumstats.VALID, self.genopheno.VALID]
        self.valid_str  = ",".join([str(v) for v in self.valids]) 
        self.verify_chromosomes() 

        #if self.args.debug_level == 2:                                                                                                                                                                                      
        #    if prevPop:                                              self.verify_ss_data(prevPop, self.sumstats)                                                                                                                    
        #    elif self.args.module not in ['pipeline','easyrun']:     self.verify_ss_data(self.sumstats, None)    



    def validate(self, F, P, L):
        mInputs, mGen = [], [] 
        if self.args.module == 'check': return True
        if self.args.cmd in ['go','run','clump','beta'] and not self.sumstats.VALID: mInputs.append('Sumstats Data')  
        if self.args.cmd in ['go','run','predict','quantify'] and not self.genopheno.VALID: mInputs.append('Genotype/Phenotype Data')  
        if self.args.cmd in ['beta']: 
            if 'clump' not in P: mGen.append('Clump Data (Hint: Run '+self.args.module+' clump)')
        if self.args.cmd in ['predict','quantify']: 
            if 'beta' not in P and self.args.module != 'prs-port': mGen.append('Weight Data (Hint: Run '+self.args.module+' beta)')
            if self.args.cmd == 'quantify' and 'predict' not in P: mGen.append('Pred Data (Hint: Run '+self.args.module+' predict') 
        if self.args.module in ['prs-port','prs-prior'] and 'model' not in F: mInputs.append('Base Model --model_file (Hint: Run build-model)') 
        if len(mInputs) + len(mGen) == 0: return True 
        if len(mInputs) > 0: bridge_pop_error(['Missing Input Data:']+mInputs) 
        if len(mGen) > 0:    bridge_pop_error(['Missing Run Data:']+mGen) 
        return True  
        

    def get_chr_strs(self, itbl): 
        c_int, c_str = [], [] 
        for c in itbl: 
            try:               c_int.append(int(c))  
            except ValueError: c_str.append(c) 
        return [str(c) for c in sorted(c_int) + sorted(c_str)]



    def verify_chromosomes(self): 
        chromosomes       = list(set([k for k in self.sumstats.map.keys()]+[k for k in self.bdata.map.keys()]))
        ld_chrs = self.get_chr_strs(self.bdata.map.keys()) 
        ss_chrs = self.get_chr_strs(self.sumstats.map.keys())
        valid_chromosomes = [k for k in self.sumstats.map.keys() if k in self.bdata.map.keys()] 
        self.chromosomes       = self.get_chr_strs(list(set([k for k in self.sumstats.map.keys()]+[k for k in self.bdata.map.keys()])))
        if len(self.chromosomes) == len(ld_chrs) and len(self.chromosomes) == len(ss_chrs) and len(self.chromosomes) == len(valid_chromosomes): return
        my_error = ['Missing Chromosomes'] 
        my_error.append('The Following Chromosomes Are Found in the ld-panel: '+",".join(ld_chrs)) 
        my_error.append('The Following Chromosomes Are Found in the sumstats: '+",".join(ss_chrs))
        if len(valid_chromosomes) == 0: my_error.append('The Following Chromosomes Are Found in both sources: None') 
        else:                           my_error.append('The Following Chromosomes Are Found in both sources: '+",".join(valid_chromosomes)) 
        my_error.append('Note: A one-to-one mapping is required, consider changing sumstats or ld-panel chromosome names')  
        bridge_pop_error(my_error) 
        return 


    def verify_ss_data(self, ss1, ss2): 
        self.progress.quikout('\nBridgeDebugData(Lvl'+str(self.args.debug_level)+'):')
        pTypes = [] 
        for si,s in enumerate([ss1, ss2]): 
            if s is not None: 
                pk, pn, pt = s.pop_type, s.pop_name, s.phenoTYPE, 
                eR = min(s.WT),max(s.WT)
                self.progress.quikout('\nPOP'+str(si+1)+' '+pn+' ('+pk+')')
                try: 
                    self.progress.quikout('Phenotype: '+self.args.phenotype+' ('+pt+')') 
                    eStr = 'Effect Size Range: '+str(eR[0])+', '+str(eR[1]) 
                    self.progress.quikout(eStr) 
                    if pt == 'binary' and eR[0] > 0:  bridge_debug_error('Truncated Effect Size Range, log(odds) is required.') 
                except: 
                    self.progress.quikout('Phenotype: '+str("Not Declared")) 
                
                g1, g2 = [], [] 
                for nk,n in [['MATCH','Matching'],['SWAPREF','Ref/Alt Swap'],['REVCOMP','Reverse Complement'],['INVALID','Invalid Bases']]: 
                    if s.CK[nk] > 0 :        g1.append(str(s.CK[nk])+' ('+n+')') 
                    if s.CK['GENO_'+nk] > 0: g2.append(str(s.CK['GENO_'+nk])+' ('+n+')') 
                if si == 0: 
                    self.progress.quikout('Initial Genotype Variants: '+str(s.snp_list_len)+' (SNP LIST), '+str(s.genome_snps)+' (VALID GENOTYPED SNPS)') 
                    self.progress.quikout('Shared GWAS Variants:      '+', '.join(g1))  
                else: 
                    self.progress.quikout('Shared Genotype Variants: '+','.join(g2)) 
                    self.progress.quikout('Shared GWAS Variants:     '+', '.join(g1)) 
        return 


























class BData: 
    def __init__(self, pop, P):
        self.map = dd(list) 
        self.VALID, self.BYCHR = False, False 
        if P['ld_path'] != None: self.load_bdata(P) 


    def load_bdata(self, P): 
        self.VALID, self.BYCHR = True, True
        self.ld_pop, self.ld_path = P['ldpop'], P['ld_path']
        self.id_file, self.prefix, self.map = self.get_ld_data(P['ldpop'], P['ld_path']) 
        self.X_fields = ['--clump-field',P['ssf-p'], '--clump-snp-field',P['ssf-snpid']] 
        return         


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
        if    len(cands) == 0: bridge_ld_error('No plink triples (bed, bim, fam) found in LD-Reference Path: '+ld_path)                                                                                                     
        elif  len(cands) == 1: bridge_ld_error('LD-Reference files must be split by chromosome: '+ld_path)                                                                                                                  
        elif len(id_files) == 0: bridge_ld_error('No Corresponding ID file found for population '+ldpop+' in: '+ld_path)                                                                                                    
        elif len(id_files) > 1: bridge_ld_error('Multiple ID files found for pop '+ldpop+' in: '+ld_path+'  ('+",".join(id_files)+')')                                                                                      
        BYCHR = True                                                                                                                                                                                                        
        my_prefix, my_suffix = get_prefix_suffix(cands)                                                                                                                                                                
        for cand in cands:                                                                                                                                                                                                  
            chr_cand = cand.split(my_prefix)[-1].split(my_suffix)[0]                                                                                                                                                        
            key[chr_cand] = ld_path+'/'+cand                                                                                                                                                                                
        return ld_path+'/'+id_files[0], ld_path+'/'+my_prefix, key 




class GenoPheno: 
    def __init__(self, pop, P): 
        self.TESTS = dd(bool) 
        self.VALID, self.BYCHR = False, False 
        self.fields, self.map = {}, dd(list) 
        self.pop_name, self.pop_type = pop.name, pop.type 
        self.phenotype, self.covariates = P['phenotype'], P['covariates']
        if P['genotype_prefix'] and P['phenotype_file']: self.ph_fill(pop, P, P['validation_file']) 


    def ph_fill(self, pop, P, validation_file = None):
        self.VALID, self.type, self.file, col_data = True, 'binary', P['phenotype_file'], dd(list) 
        self.ph_get_genotypes(P['genotype_prefix']) 
        if pop.type.upper() == 'TARGET': 
            if validation_file is not None: self.files = [P['phenotype_file'], validation_file]
            else:                           self.files = self.split_phenotype_file(pop.paths['save'],P['phenotype_file']) 
            self.X_fields = ['--test.data',self.files[0],'--valid.data',self.files[1]]
        else: 
            self.files = [P['phenotype_file']]
            self.X_fields = ['--test.data',self.files[0],'--valid.data','0'] 
        for i,fn in enumerate(self.files): 
            with open(fn, 'rt') as f: 
                lp = f.readline().split() 
                if i == 0: self.header, lz = [x for x in lp] , ",".join(lp) 
                elif ",".join(lp) != lz: bridge_pop_error('Phenotype File Headers Do Not Match: '+lz+' AND '+",".join(lp))
                for k,line in enumerate(f):
                    line = line.split() 
                    for j,c in enumerate(self.header): col_data[c].append(line[j]) 
                    if k > 100: break 

        pheno_cands = [h for h in self.header[1::] if 'ID' not in h]
        if self.covariates is not None: 
            self.fields['COVARIATES'] = self.covariates
            self.X_fields.extend(['--cov.names',self.covariates]) 
            for c in self.covariates.split(','): 
                if c not in self.header:  bridge_pop_error('Invalid covariate field name(s) supplied '+self.covariates+', Available Fields: '+','.join(self.header)) 
                else:                     pheno_cands.pop(pheno_cands.index(c)) 



        if self.phenotype is None: bridge_pop_error(['Phenotype field required (--phenotype)',', Available Fields: '+','.join(self.header)]) 


        if self.phenotype is not None: 
            self.fields['NAME'] = self.phenotype 
            self.X_fields.extend(['--pheno.name',self.phenotype]) 
            


            if self.phenotype not in self.header: bridge_pop_error('Invalid phenotype field name(s) supplied '+self.phenotype+', Available Fields: '+','.join(self.header)) 
            if len((list(set(col_data[self.phenotype])))) > 2:   self.type = 'continuous' 
            else:                                                self.X_fields.extend(['--binary','1']) 
       




        return self 
        
        

    def ph_get_genotypes(self, genotype_prefix): 
        self.genotype_prefix = genotype_prefix
        CK, g_path, g_prefix = dd(int), "/".join(self.genotype_prefix.split('/')[0:-1]), self.genotype_prefix.split('/')[-1] 
        for f in [x for x in os.listdir(g_path) if x[0:len(g_prefix)] == g_prefix]:
            if f.split('.')[-1] in ['bed','bim','fam']: CK[".".join(f.split('.')[0:-1])]+=1                                                                                                                                                                                   
        cands = [k for k,v in CK.items() if v == 3]
        if len(cands) == 0: bridge_pop_error('Invalid genotype data prefix '+self.genotype_prefix) 
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



















class SumStats: 
    def __init__(self, pop, P, genoPheno, prevPop): 
            
        #self.args, self.progress, self.datatype, self.paths, self.pop_name, self.pop_type = args, progress, datatype, paths, pop_name, pop_type
        
        self.debug_level = pop.args.debug_level 
        
        # X_fields  

        self.pop_name, self.pop_type, self.savepath = pop.name, pop.type, pop.paths['save'] 
        self.VALID, self.TESTS, self.map = False, dd(bool), dd(list) 
        self.fields, self.map, self.total = {}, dd(list), 0 
        if P['sumstats_prefix']: self.add_sumstats(P, genoPheno, prevPop)  

        self.max_clump_size = P['max_clump_size'] 
        if self.max_clump_size == None: self.max_clump_size = '0' 





    ####################################   SUMSTATS   #######################################
    
    def base_comp(self, g1, g2): 
        try: 
            if g1[0] == g2[0] and g1[1] == g2[1]: return 'MATCH' 
            if g1[0] == g2[1] and g1[1] == g2[0]: return 'SWAPREF' 
            if REV_COMP[g1[0]] == g2[0] and REV_COMP[g1[1]] == g2[1]: return 'REVCOMP' 
        except: pass 
        return 'INVALID'  


    def begin_sumstats_file(self, pf, CHR=None):
        required, locs, names = ['CHR','SNPID','REF','ALT','MAF','N','BETA','SE','P'], [], [] 
        
        required = ['CHR','SNPID','REF','ALT','P','BETA']      





        FK, NK, IK, ZK = {b: a for a,b in self.fields.items()}, {'CHR': 'CHR'}, {}, []  
        p_handle = zip_open(pf) 
        p_init   = p_handle.readline()  
        p_header = p_init.split()   
        
        for i,h in enumerate(p_header): 
            if   'CHR' in h.upper(): lk = 'CHR' 
            elif  h in FK:    lk = FK[h] 
            else: 
                ZK.append(h)     
                continue 
            if lk in IK: bridge_sumstats_error(['Repeated Column']) 
            IK[lk], NK[lk] = i, h     
        for r in required: 
            if r not in IK: 
                if r != 'CHR':                  
                    my_error = ['Invalid Sumstats File: '+pf,'    Missing Header Field: --SSF-'+r+' '+self.fields[r]]
                    if len(ZK) > 0: my_error.append('   See Available Choices:  '+",".join(ZK))  
                    else:           my_error.append('    Please Add Column to Sumstats File')  
                    bridge_sumstats_error(my_error) 
                elif self.debug_level < 2: 
                    bridge_sumstats_error(['Cannot Split Sumstats, CHR field not found','    Fields Found: '+",".join(p_header)]) 
                else:                           
                    IK['CHR'] = 'NA'  
            locs.append(IK[r])
            names.append(NK[r]) 
        if IK['CHR'] == 'NA': self.TESTS['INFER_CHR'] = True 
        return p_handle, "\t".join(names), locs, IK, IK['CHR'], IK['SNPID'] 


    def add_sumstats_line(self,fC,LD,snpID): 
        self.total += 1 
        if fC not in self.s_key: 
            self.s_key[fC] = [self.prefix+str(fC)+'.out', open(self.prefix+str(fC)+'.out','w')]
            self.s_key[fC][1].write(self.header+'\n')
        self.s_key[fC][1].write("\t".join(LD)+'\n') 
        
        if self.TESTS['NOSNPS']: 
            self.snp_handle.write(snpID+'\n')
        return




    def split_sumstats(self, p_file, CHR = 'NA'):
        self.s_key  = {} 
        p_handle, self.header, locs, IK, iC, iS = self.begin_sumstats_file(p_file)

        if self.debug_level < 1: 
            for line in p_handle: 
                lp = line.split()     
                try: fC, LD = int(lp[iC]), [lp[j] for j in locs] 
                except: bridge_sumstats_error('Increase Debug Level (--debug_level 2) or supply numerical chromsome information in the sumstats files') 
                self.add_sumstats_line(fC, LD, lp[iS]) 
        else:
            for li,line in enumerate(p_handle): 
                lp = line.split() 
                if lp[iS] not in self.rs_key: 
                    self.CK['NO_GENOTYPE'] += 1 
                    continue 
                rsData = self.rs_key[lp[iS]] 
                if len(rsData) == 4:   rsChr, rsLoc, rsRef, rsAlt = self.rs_key[lp[iS]] 
                elif len(rsData) == 5: rsChr, rsLoc, rsRef, rsAlt, rsBool = self.rs_key[lp[iS]] 
                else:                  bridge_debug_error('Invalid RS Data:'+','.join([str(xxx) for xxx in rsData])) 
                
                LD = [lp[j] if j != 'NA' else j for j in locs]

                #print(p_file) 
                #print(self.header)  
                #print(locs) 
                    
                #print(lp) 

                #print(LD) 
                #sys.exit() 

                #['9', '286593', 'rs529208', 'A', 'C', 'C', '0.2944', 'ADD', '20000', '-0.00971996', '0.0109591', '-0.886932', '0.375126', '.']
                #['9', 'rs529208', 'A', 'C', '0.375126', '-0.00971996']



                #[0, 2, 3, 5, 12, 9]

                #CHROM POS ID REF ALT A1 A1_FREQ TEST OBS_CT BETA SE T_STAT P ERRCODE
                #9 286593 rs529208 A C C 0.2944 ADD 20000 -0.00971996 0.0109591 -0.886932 0.375126 .
                #9 407227 rs4741867 G A A 0.48745 ADD 20000 0.00266047 0.0100489 0.264751 0.791204 .




                if 'NA' in LD: 
                    self.CK['NA_LINE'] += 1 
                    continue 
                if rsRef == rsRef.upper(): 
                    LD[2] = LD[2].upper() 
                    LD[3] = LD[3].upper() 
                chr_cands = list(set([c for c in [str(LD[0]), str(CHR), str(rsChr)] if c != 'NA'])) 
                if len(chr_cands) > 1: bridge_sumstats_error('Ambiguous Chromosomes For '+lp[iS]+': '+str(LD[0])+','+str(rsChr)+','+str(CHR)+' (Sumstats, Genotype, Filename(s))') 
                try:               
                    fC, LD[0] = int(chr_cands[0]), chr_cands[0]
                except ValueError: 
                    bridge_sumstats_error(['Nonnumerical Chromosome ('+chr_cands[0]+')','    Sumstats File: '+p_file])  
                try:                
                    
                    #lpRef, lpAlt, lpMaf, lpN, lpWt, lpSE, lpP = LD[2], LD[3], float(LD[4]), float(LD[5]), float(LD[6]), float(LD[7]), float(LD[8])
                    lpRef, lpAlt,  lpP, lpWt = LD[2], LD[3], float(LD[4]), float(LD[5]) 
                except ValueError:  
                    bridge_sumstats_error(['Sumstats File Error(s), Incorrect DataType:',self.header,'\t'.join(LD)]) 
               
                relationship = self.base_comp([lpRef, lpAlt], [rsRef.upper(), rsAlt.upper()])
                self.CK[relationship] += 1 
                if relationship != 'INVALID': 
                    if len(self.WT) < 2 or lpWt < min(self.WT) or lpWt > max(self.WT):  self.WT.append(lpWt) 
                    self.add_sumstats_line(fC, LD, lp[iS]) 
                    self.rs_key[lp[iS]].append(True) 
        
            cMatch, cInvalid, cSwap, cRev, cTotal = self.CK['MATCH'], self.CK['INVALID'], self.CK['SWAPREF'], self.CK['REVCOMP'], sum([kk for kk in self.CK.values()]) 
            if cInvalid == cTotal: bridge_sumstats_error(['Sumstats File Error(s), No Matching Ref/Alt Bases']) 
            elif cInvalid > cTotal/2.0: bridge_sumstats_error(['Sumstats File Error(s), Majority MisMatching Ref/Alt Bases - Please Check Builds']) 

     
        p_handle.close() 
        for sc,sd in self.s_key.items(): 
            sd[1].close() 
            os.system('gzip -f '+sd[0]) 
            self.map[str(sc)] = sd[0]+'.gz' 
        return 



    def initialize_sumstats(self,pk,genoPheno):

        if pk['thinned_snp_file'] is None: self.thin_snps = '0' 
        else:                                  self.thin_snps = pk['thinned_snp_file']
        self.rs_key, self.snp_file = {}, pk['snp_file'] 
        if self.snp_file is None:
            self.TESTS['NEWSNPS']  = True 
            self.TESTS['NOSNPS']   = True
            self.snp_file          = self.savepath+'/snps.'+self.pop_name.lower()+'_valid.txt'
            self.snp_handle        = open(self.snp_file,'w') 
        
        if self.debug_level >= 1: 
            if not self.TESTS['NOSNPS']: 
                with open(self.snp_file) as f: self.rs_key = {line.strip(): [] for line in f} 
            else:                              self.rs_key = {} 
            geno_path, geno_name = '/'.join(genoPheno.genotype_prefix.split('/')[0:-1]), genoPheno.genotype_prefix.split('/')[-1] 
            geno_files = [f for f in os.listdir(geno_path) if f[0:len(geno_name)] == geno_name and f.split('.')[-1]=='bim'] 
            for gf in geno_files: 
                with open(geno_path+'/'+gf) as f: 
                    for line in f: 
                        line = line.split() 
                        try: chr_name, rs, loc, ref, alt = int(line[0]), line[1], int(line[3]), line[4], line[5] 
                        except: bridge_pop_error('Invalid genotype file: '+gf) 
                        if self.TESTS['NOSNPS']: self.rs_key[rs] = [chr_name, loc, ref, alt] 
                        elif rs in self.rs_key:  self.rs_key[rs] = [chr_name, loc, ref, alt]  
                        else:                    continue  
        
            self.snp_list_len = len(self.rs_key) 
            self.rs_key = {a:b for a,b in self.rs_key.items() if len(b) > 0} 
            self.genome_snps = len(self.rs_key)


    def continue_sumstats(self,pk,genoPheno, prevPop): 


        self.TESTS['NOSNPS'] = prevPop.TESTS['NOSNPS'] 
        self.rs_key, self.snp_file, self.thin_snps = {}, prevPop.snp_file, prevPop.thin_snps 
        


        if self.debug_level >= 1:
            self.snp_list_len = prevPop.snp_list_len 
            geno_path, geno_name = '/'.join(genoPheno.genotype_prefix.split('/')[0:-1]), genoPheno.genotype_prefix.split('/')[-1] 
            geno_files = [f for f in os.listdir(geno_path) if f[0:len(geno_name)] == geno_name and f.split('.')[-1]=='bim'] 
            for gf in geno_files: 
                with open(geno_path+'/'+gf) as f: 
                    for line in f: 
                        line = line.split() 
                        try: chr_name, rs, loc, ref, alt = int(line[0]), line[1], int(line[3]), line[4], line[5] 
                        except: bridge_pop_error('Invalid genotype file: '+gf) 
                        
                        if rs not in prevPop.rs_key or len(prevPop.rs_key[rs]) < 5: self.CK['GENO_EXTRA'] += 1 
                        else:
                            rC, rL, rRef, rAlt, rBool = prevPop.rs_key[rs] 
                            if chr_name != rC or loc != rL: 
                                azd = str(rC)+':'+str(rL)+', '+str(chr_name)+':'+str(loc) 
                                if   chr_name != rC: bridge_sumstats_error('Invalid Genome Builds For Base/Target Genotype, '+rs+' On Multiple Chromosomes: '+azd) 
                                else:                bridge_sumstats_error('Invalid Genome Builds For Base/Target Genotype, '+rs+' At Multiple Locations: '+azd) 
                            else:
                                self.CK['GENO_'+self.base_comp([ref,alt],[rRef, rAlt])] += 1 
                                self.rs_key[rs] = [rC, rL, rRef, rAlt] 
            self.genome_snps = len(self.rs_key)+self.CK['GENO_EXTRA']

    
    
    
    def add_sumstats(self, pk, genoPheno, prevPop): 
        NN_ERROR = 'Nonnumerical Chromosomes or Ambiguous Sumstats Names' 
        self.phenoTYPE      = genoPheno.type
        sum_stats_fields    =  ['ssf-alt', 'ssf-beta', 'ssf-maf', 'ssf-p', 'ssf-ref', 'ssf-se', 'ssf-snpid', 'ssf-n']
        sum_stats_args      =  [pk[ks] for ks in sum_stats_fields] 
        self.fields         =  {a.split('-')[-1].upper(): b for a,b in zip(sum_stats_fields, sum_stats_args)}
        self.CK, self.WT, self.VALID, self.BYCHR = dd(int), [], True, True
        if not prevPop: self.initialize_sumstats(pk, genoPheno) 
        else:           self.continue_sumstats(pk, genoPheno, prevPop) 
        self.source_prefix, self.source_suffix = pk['sumstats_prefix'], pk['sumstats_suffix'] 
        prefix_path, prefix_name = '/'.join(self.source_prefix.split('/')[0:-1]), self.source_prefix.split('/')[-1] 
        prefix_files = [f for f in os.listdir(prefix_path) if f[0:len(prefix_name)] == prefix_name] 

        if len(prefix_files) < 5: d_files = ",".join(sorted(prefix_files)) 
        else:                     
            d_files = "\n                                      ".join(sorted(prefix_files[0:3]))+'\n                                            ' 
            d_files +="\n                                            ".join(['.','.','.'])+'\n                                      ' 
            d_files += "\n                                      ".join(sorted(prefix_files[-3::]))+'\n' 

        if   len(prefix_files) == 0: bridge_sumstats_error('Invalid Sumstats Prefix: '+prefix) 
        elif self.debug_level == 0 and len(prefix_files) > 1: self.prefix, self.suffix = self.source_prefix, self.source_suffix 
        else: 
            new_prefix_path = self.savepath+'/sumstats'

            if not os.path.exists(new_prefix_path): os.makedirs(new_prefix_path) 
            self.prefix, self.suffix = new_prefix_path +'/ss.'+self.pop_name+'.', '.out.gz'
        
        if self.source_suffix == 'FILE': prefix_files = [prefix_name] 
        elif self.source_suffix:
            suffix_files = [pf for pf in prefix_files if pf[-1*len(self.source_suffix)::] == self.source_suffix] 
            if len(suffix_files) == 0: bridge_sumstats_error(['Incorrect prefix/suffix names','    prefix/suffix: '+self.prefix+', '+str(self.source_suffix),'    Candidate Files: '+d_files]) 
            else:                      prefix_files = suffix_files 


        if len(prefix_files) == 1: 
            bridge_sumstats_warning('Sumstats files for pop '+self.pop_name+' are not separated by chromosome, attempting to split file in '+new_prefix_path) 
            self.source_suffix = 'NA'
            self.split_sumstats(prefix_path+'/'+prefix_files[0])  
        else:
            if self.source_suffix:
                chromosomes = [pf.split(self.source_suffix)[0].split(prefix_name)[1] for pf in prefix_files] 
                try:       chromosomes = [int(c) for c in chromosomes]
                except:    
                    bridge_sumstats_error([NN_ERROR,'    prefix/suffix: '+self.prefix+', '+str(self.source_suffix),'    Directory: '+', '.join(prefix_files)]) 
                
            if not self.source_suffix:
                self.TESTS['INFER_SUFFIX'] = True 
                if self.debug_level == 0: 
                    bridge_sumstats_error('Increase Debug Level (--debug_level) or supply sumstats suffix for pop '+self.pop_name+' (on the Command Line (--sumstats_suffix) or in a config_file (SUMSTATS_SUFFIX=))')  
                try:
                    new_prefix, new_suffix = self.get_prefix_suffix(prefix_files) 
                    chromosomes = [pf.split(new_suffix)[0].split(new_prefix)[-1] for pf in prefix_files] 
                    chromosomes = [int(c) for c in chromosomes]
                    self.source_prefix        = prefix_path+'/'+new_prefix
                    self.source_suffix = new_suffix 
                except: 
                    suffix_cands, test_suffixes = [], ''
                    try: 
                        for p in prefix_files:
                            k, c_cand, NEED_CHR = 0, p.split(prefix_name)[-1], False 
                            while k < len(c_cand) and c_cand[k] in NUM_STRS: k+=1  
                            suffix_cands.append(c_cand[k::]) 
                        
                        test_suffixes = ','.join(list(set(suffix_cands))) 
                        if len(list(set(suffix_cands))) != 1:  
                            bridge_sumstats_error('Ambiguous Sumstats Suffixes ('+test_suffixes+'), please supply sumstats suffix for pop '+self.pop_name+' (on Command Line (--sumstats_suffix) or in a config_file (SUMSTATS_SUFFIX=))')  
                        self.source_suffix = suffix_cands[0] 
                        chromosomes = [pf.split(self.source_suffix)[0].split(prefix_name)[-1] for pf in prefix_files] 
                        chromosomes = [int(c) for c in chromosomes]
                    except:
                        error_list = ['Ambiguous Sumstats Names','    prefix/suffix: '+self.prefix+', '+str(self.source_suffix),] 
                        if len(test_suffixes.split(',')) > 1: 
                            bridge_sumstats_error(error_list+['    Please Supply sumstats suffix for pop '+self.pop_name+' (on the Command Line (--sumstats_suffix) or in a config_file (SUMSTATS_SUFFIX=))'])  
                        bridge_sumstats_error(['Nonnumerical Chromosomes or Ambiguous Sumstats Names','    prefix/suffix: '+self.prefix+', '+str(self.source_suffix),'    Directory: '+', '.join(prefix_files)]) 

            for c,p in zip(chromosomes, prefix_files): 
                if self.debug_level > 0: self.split_sumstats(prefix_path+'/'+p, c) 
                else:                         self.map[str(c)] = prefix_path+'/'+p 


        if self.TESTS['NOSNPS']: 
            try: 
                self.snp_handle.close() 
                self.TESTS['NOSNPS'] = False 
                #if self.args.debug_level > 0: self.TESTS['NOSNPS'] = False 
                #else:                         self.TESTS['NOSNPS'] = True
            except: pass 
        


        self.X_fields   = ['--sumstats.snpID',self.fields['SNPID'],'--sumstats.allele0ID',self.fields['REF'],'--sumstats.allele1ID',self.fields['ALT'],'--sumstats.betaID',self.fields['BETA']]
        #self.X_fields.extend(['--sumstats.frqID',self.fields['MAF'],'--sumstats.nID',self.fields['N'],'--sumstats.seID',self.fields['SE'], '--sumstats.snpID',self.fields['SNPID']])
        
        return





































































