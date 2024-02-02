import sys, os, gzip  
from collections import defaultdict as dd
from .BridgeProgress  import BridgeProgress
from math import log 

def bridge_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeDataError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('\nBridgeDataError: '+eString+'\n')
    sys.exit(2) 

            
def zip_open(fp, HEADER = False): 
    if fp.split('.')[-1] == 'gz': gf = gzip.open(fp, 'rt') 
    else:                         gf = open(fp, 'rt')  
    if HEADER: 
        HL = gf.readline().split() 
        gf.close() 
        return HL 
    else: return gf 
        




class BridgePop: 
    def __init__(self, args, paths, pop_name, pop_key, ld_ref, pop_type = 'TARGET'):
        self.args, self.paths, self.name, self.key, self.type  = args, paths, pop_name, pop_key, pop_type
        self.ref_pop = self.key['ldpop'] 
        
        self.bdata      = BridgeData(self.args, 'bdata',      paths, self.ref_pop)
        self.sumstats   = BridgeData(self.args, 'sumstats',   paths, self.ref_pop) 
        self.phenotypes = BridgeData(self.args, 'phenotypes', paths, pop_name, pop_type) 
        self.bdata.add_panel(ld_ref[self.ref_pop]) 
        
        
        if pop_key['sumstats_prefix']: 
            self.sumstats.add_stats(pop_key['sumstats_prefix'],pop_key['sumstats_suffix'], self.args.snp_file)
        if pop_key['genotype_prefix'] and pop_key['phenotype_file']: 
            self.phenotypes.ph_fill(pop_key['genotype_prefix'],pop_key['phenotype_file'], self.args.validation_file) 
        
        
        self.verify_chromosomes() 

    def verify_chromosomes(self): 

        chromosomes       = list(set([k for k in self.sumstats.map.keys()]+[k for k in self.bdata.map.keys()]))
        valid_chromosomes = [k for k in self.sumstats.map.keys() if k in self.bdata.map.keys()] 
        c_int, c_str = [], [] 
        for c in valid_chromosomes: 
            try:               c_int.append(int(c))  
            except ValueError: c_str.append(c) 
        self.chromosomes = [str(c) for c in sorted(c_int) + sorted(c_str)]
        missing_chrs = [c for c in chromosomes if c not in valid_chromosomes] 
        if len(missing_chrs) > 0: 
            b_missing = [c for c in self.bdata.map.keys() if c not in self.sumstats.map.keys()]
            s_missing = [c for c in self.bdata.map.keys() if c not in self.sumstats.map.keys()]
            bridge_error('Missing Chromosomes') 

        return 



    def validate(self, F, P, L):
        
        mInputs, mGen = [], [] 
        if self.args.module == 'check': return True
        if self.args.cmd in ['go','run','clump','beta'] and not self.sumstats.VALID: mInputs.append('Sumstats Data')  
        if self.args.cmd in ['go','run','predict','quantify'] and not self.phenotypes.VALID: mInputs.append('Genotype/Phenotype Data')  
        
        if self.args.cmd in ['beta']: 
            if 'clump' not in P: mGen.append('Clump Data (Hint: Run '+self.args.module+' clump)')
        
        if self.args.cmd in ['predict','quantify']: 
            if 'beta' not in P and self.args.module != 'prs-port': mGen.append('Weight Data (Hint: Run '+self.args.module+' beta)')
            if self.args.cmd == 'quantify' and 'predict' not in P: mGen.append('Pred Data (Hint: Run '+self.args.module+' predict') 

        if self.args.module in ['prs-port','prs-prior'] and 'model' not in F: mInputs.append('Base Model --model_file (Hint: Run build-model)') 

        if len(mInputs) + len(mGen) == 0: return True 

        if len(mInputs) > 0: bridge_error(['Missing Input Data:']+mInputs) 
        if len(mGen) > 0:    bridge_error(['Missing Run Data:']+mGen) 
        return True  
        
        
            





class BridgeData: 
    def __init__(self, args, datatype, paths, pop_name,pop_type = 'TARGET'):
        self.args, self.datatype, self.paths, self.pop_name, self.pop_type = args, datatype, paths, pop_name, pop_type
        self.VALID, self.TESTS, self.map = False, dd(bool), dd(list) 
    
        if self.datatype == 'sumstats':    
            self.fields, self.map, self.total = {}, dd(list), 0 
        elif self.datatype == 'phenotypes':   
            self.fields, self.map = {}, dd(list) 
        
        


    
    def add_panel(self, ld_quad): 
        self.VALID = True 
        self.id_file, self.prefix, self.BYCHR, b_key  = ld_quad
        self.ldpath = "/".join(self.prefix.split('/')[0:-1])
        for c in b_key: self.map[c] = b_key[c] 
        self.X_fields = ['--clump-field',vars(self.args)['ssf-p'], '--clump-snp-field',vars(self.args)['ssf-snpid']] 
        return self 



    

    ####################################   SUMSTATS   #######################################
    
        
        
    def add_stats(self,prefix,suffix, snp_file, thin_snps = '0'): 
        self.VALID, self.BYCHR = True, True 
        self.prefix, self.prefix_path, pName, pLen, self.snp_file, self.thin_snps  = prefix, '/'.join(prefix.split('/')[0:-1]), prefix.split('/')[-1], len(prefix.split('/')[-1]), snp_file, thin_snps 
        p_files, suffix_cands = [self.prefix_path+'/'+x for x in os.listdir(self.prefix_path) if x[0:pLen] == pName], dd(int) 
        for pn in p_files: 
            k,p_tail = 0, pn.split(self.prefix)[-1] 
            while p_tail[k] in ['1','2','3','4','5','6','7','8','9','0']: k+=1 
            suffix_cands[p_tail[k::].split('.gz')[0]] += 1  
        suffix_pairs, suffix_cands = [], [sc[0] for sc in sorted(suffix_cands.items(), key = lambda X: X[1], reverse=True)] 
        if not suffix:                               suffix_prefix = suffix_cands[0] 
        elif suffix.split('.gz')[0] in suffix_cands: suffix_prefix = suffix.split('.gz')[0] 
        else:                                        bridge_error(['Invalid suffix supplied, --sumstats_suffix '+suffix,'Does not match file prefix: '+self.prefix]) 
        suffix_pairs = [[pn.split(self.prefix)[-1].split(suffix_prefix)[0], pn] for pn in p_files if len(pn.split(self.prefix)[-1].split(suffix_prefix)) == 2] 
        if suffix == None: self.TESTS['INFER_SUFFIX'] = True 
        if len(suffix_pairs) == 1: suffix_pairs  = self.ss_split_suffix_matches(suffix_pairs[0][1], suffix_prefix)  
        f_cnt = dd(int) 
        for f_chr, fp in suffix_pairs:    
            for x in zip_open(fp, HEADER = True): f_cnt[x] += 1 
            if fp.split('.')[-1] != 'gz': 
                os.system('gzip -f '+fp)  
                fp += '.gz' 
            self.map[f_chr] = fp 
        self.suffix = suffix_prefix 
        if self.suffix.split('.')[-1] != 'gz': self.suffix+='.gz' 
        sum_stats_fields    =  ['ssf-alt', 'ssf-beta', 'ssf-maf', 'ssf-p', 'ssf-ref', 'ssf-se', 'ssf-snpid', 'ssf-n']
        self.fields         =  {ks.split('-')[-1].upper(): vars(self.args)[ks] for ks in sum_stats_fields}
        cands = [x for x,y in f_cnt.items() if y == len(self.map)] 
        cand_found, cand_errors = [y for x,y in self.fields.items() if y in cands], ['--ssf-'+x.lower()+' '+y for x,y in self.fields.items() if y not in cands] 
        if len(cand_errors) > 0:  bridge_error(['Invalid Sumstats Fields(s):']+cand_errors) 
        self.X_fields   = ['--sumstats.allele0ID',self.fields['REF'],'--sumstats.allele1ID',self.fields['ALT'],'--sumstats.betaID',self.fields['BETA']]
        self.X_fields.extend(['--sumstats.frqID',self.fields['MAF'],'--sumstats.nID',self.fields['N'],'--sumstats.seID',self.fields['SE'], '--sumstats.snpID',self.fields['SNPID']])
        self.ss_test_snpfile() 
        return self 

    
    
    def ss_split_suffix_matches(self, sumstat_file, suffix): 
        gf, new_path = zip_open(sumstat_file), self.paths['save']+'/'+self.prefix.split('/')[-1] 
        header      = gf.readline().strip() 
        header_list = [ih for ih,h in enumerate(header.split()) if 'CHR' in h.upper()]
        if len(header_list) != 1:  bridge_error('Cannot Split Sumstats, CHR field not found') 
        SP, W, chr_loc = [], {}, header_list[0] 
        for line in gf: 
            line, lc = line.strip(), line.split()[chr_loc] 
            if lc not in W:    
                W[lc] = open(new_path+lc+suffix,'w') 
                W[lc].write(header+'\n') 
            W[lc].write(line+'\n') 
        for lc in W: 
            W[lc].close() 
            SP.append([lc, new_path+lc+suffix]) 
        self.prefix = new_path 
        gf.close() 
        return SP 

        

    def ss_test_snpfile(self): 
        if self.snp_file is None: 
            self.TESTS['NOSNPS'] = True 
            snp_handle = self.paths['save']+'/snps.'+self.pop_name+'.txt' 
            w = open(snp_handle,'w') 
        
        self.total  = 0 
        for f_chr, fp in self.map.items():  
            gf = zip_open(fp) 
            snp_loc = {h: i for i,h in enumerate(gf.readline().split())}[self.fields['SNPID']]
            snps = [line.split()[snp_loc] for line in gf] 
            self.total += len(snps) 
            if self.TESTS['NOSNPS']: w.write('\n'.join(snps)) 
            gf.close() 
        if self.TESTS['NOSNPS']: 
            self.snp_file = snp_handle 
            w.close()
        return  



    ####################################   PHENOTYPES   #######################################
    

    def ph_get_genotypes(self, genotype_prefix): 
        self.genotype_prefix = genotype_prefix
        CK, g_path, g_prefix = dd(int), "/".join(self.genotype_prefix.split('/')[0:-1]), self.genotype_prefix.split('/')[-1] 
        for f in [x for x in os.listdir(g_path) if x[0:len(g_prefix)] == g_prefix]:
            if f.split('.')[-1] in ['bed','bim','fam']: CK[".".join(f.split('.')[0:-1])]+=1                                                                                                                                                                                   
        cands = [k for k,v in CK.items() if v == 3]
        if len(cands) == 0: bridge_error('Invalid genotype data prefix '+self.genotype_prefix) 
        if len(cands) == 1: self.BYCHR = False 
        else: 
            self.BYCHR, k  = True, 0 
            cands = list(set(cands))                                                                                                                                                                                                                                                        
            while True:                                                                                                                                                                                                                                                                     
                my_prefix = list(set([c[0:k] for c in cands]))[0]                                                                                                                                                                                                                              
                ck1 = list(set([c[0:k+1] for c in cands]))                                                                                                                                                                                                                                  
                if len(ck1) == 1: k+=1                                                                                                                                                                                                                                                      
                else: break                                                                                                                                                                                                                                                                 
        return                                                                                                                                                                                                                                                                                  


    def split_phenotype_file(self, fname):
        self.TESTS['NOVALID'] = True 
        with open(fname, "r") as f: 
            header = f.readline().strip() 
            data   = [ln.strip() for ln in f]
            h_len = int(len(data)/2)
        d1,d2 = data[0:h_len], data[h_len::]
        test_file, valid_file = self.paths['save']+'/'+self.pop_name+'.test_phenos.dat', self.paths['save']+'/'+self.pop_name+'.valid_phenos.dat'
        w1, w2 = open(test_file,'w'), open(valid_file,'w') 
        w1.write(header+'\n') 
        w1.write("\n".join(d1)+'\n') 
        w2.write(header+'\n') 
        w2.write("\n".join(d2)+'\n') 
        w1.close() 
        w2.close() 
        return [test_file, valid_file] 


    def ph_fill(self, genotype_prefix, phenotype_file, validation_file = None):
        self.VALID, self.type, self.file, col_data = True, 'binary', phenotype_file, dd(list) 
        self.ph_get_genotypes(genotype_prefix) 
        if self.pop_type == 'TARGET': 
            if validation_file is not None: self.files = [phenotype_file, validation_file]
            else:                           self.files = self.split_phenotype_file(phenotype_file) 
            self.X_fields = ['--test.data',self.files[0],'--valid.data',self.files[1]]
        else: 
            self.files = [phenotype_file]
            self.X_fields = ['--test.data',self.files[0],'--valid.data','0'] 
        if self.args.module != 'check' and self.args.phenotype is None: bridge_error('Phenotype field required (--phenotype)') 
        for i,fn in enumerate(self.files): 
            f = open(fn, 'rt') 
            lp = f.readline().split() 
            if i == 0: self.header, lz = [x for x in lp] , ",".join(lp) 
            elif ",".join(lp) != lz: bridge_error('Phenotype File Headers Do Not Match: '+lz+' AND '+",".join(lp))
            for k,line in enumerate(f):
                line = line.split() 
                for j,c in enumerate(self.header): col_data[c].append(line[j]) 
                if k > 100: break 
            f.close() 
        if self.args.phenotype is not None: 
            self.fields['NAME'] = self.args.phenotype 
            self.X_fields.extend(['--pheno.name',self.args.phenotype]) 
            if self.args.phenotype not in self.header: bridge_error('Invalid phenotype field name(s) supplied '+self.args.phenotype+', Available Fields: '+','.join(self.header)) 
            if len((list(set(col_data[self.args.phenotype])))) > 2:   self.type = 'continuous' 
            else:                                                     self.X_fields.extend(['--binary','1']) 
        if self.args.covariates is not None: 
            self.fields['COVARIATES'] = self.args.covariates
            self.X_fields.extend(['--cov.names',self.args.covariates]) 
            for c in self.args.covariates.split(','): 
                if c not in self.header:  bridge_error('Invalid covariate field name(s) supplied '+self.args.covariates+', Available Fields: '+','.join(self.header)) 
        return self 


