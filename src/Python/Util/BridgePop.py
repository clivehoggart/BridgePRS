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
    def __init__(self, args, paths, pop_name, pop_key, ld_ref): 
        self.args, self.paths, self.name, self.key  = args, paths, pop_name, pop_key 
        self.ref_pop = self.key['ldpop'] 
        
        
        self.bdata, self.sumstats, self.phenotypes = BridgeData(self.args, 'bdata', paths, self.ref_pop), BridgeData(self.args, 'sumstats', paths, self.ref_pop), BridgeData(self.args, 'phenotypes', paths, pop_name) 
        self.bdata.add_panel(ld_ref[self.ref_pop]) 
        self.sumstats.add_stats(pop_key['sumstats_prefix'], pop_key['sumstats_suffix'], pop_key['snp_file'])
        self.phenotypes.ph_fill(pop_key['genotype_prefix'],pop_key['phenotype_files']) 
        
        
        
        self.chromosomes = [k for k in self.sumstats.map.keys() if k in self.bdata.map.keys()] 


class BridgeData: 
    def __init__(self, args, datatype, paths, pop_name):
        

        self.args, self.datatype, self.paths, self.pop_name = args, datatype, paths, pop_name
        self.VALID, self.TESTS = False, dd(bool)  
        if self.datatype == 'bdata':         self.map = dd(list) 
        elif self.datatype == 'sumstats':    self.fields, self.map, self.total = {}, dd(list), 0 
        elif self.datatype == 'phenotypes':   self.fields, self.map = {}, dd(list) 
        
        


    
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
        #self.chromosomes = sorted([x for x in self.map.keys()]) 
        sum_stats_fields    =  ['ssf-alt', 'ssf-beta', 'ssf-maf', 'ssf-p', 'ssf-ref', 'ssf-se', 'ssf-snpid', 'ssf-ss']
        self.fields         =  {ks.split('-')[-1].upper(): vars(self.args)[ks] for ks in sum_stats_fields}
        cands = [x for x,y in f_cnt.items() if y == len(self.map)] 
        cand_found, cand_errors = [y for x,y in self.fields.items() if y in cands], ['--ssf-'+x.lower()+' '+y for x,y in self.fields.items() if y not in cands] 
        if len(cand_errors) > 0:  bridge_error(['Invalid Sumstats Fields(s):']+cand_errors) 

        self.X_fields   = ['--sumstats.allele0ID',self.fields['REF'],'--sumstats.allele1ID',self.fields['ALT'],'--sumstats.betaID',self.fields['BETA']]
        self.X_fields.extend(['--sumstats.frqID',self.fields['MAF'],'--sumstats.nID',self.fields['SS'],'--sumstats.seID',self.fields['SE'], '--sumstats.snpID',self.fields['SNPID']])
        #self.ss_analyze_snpfile(TYPE='TEST')  
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


    def ss_analyze_snpfile99(self): 
        
        if self.snp_file is not None: self.NOSNPS = False 
        else: 
            self.NOSNPS = True 
            snp_handle = self.paths['save']+'/snps.'+self.pop_name+'.txt' 
            w = open(snp_handle,'w') 
        
        self.total, self.effects = 0, {} 
        for f_chr, fp in self.map.items():  
            gf = zip_open(fp) 
            h_key = {h: i for i,h in enumerate(gf.readline().split())} 
            SPB = [h_key[self.fields['SNPID']], h_key[self.fields['P']], h_key[self.fields['BETA']]]
            null_len, effect_data = 0, [] 
            for line in gf: 
                line = line.split() 
                spb = [line[k] for k in SPB] 
                if self.NOSNPS: w.write(spb[0]+'\n') 
                if float(spb[1]) > 0.05: null_len += 1 
                else: 
                    pv, effect = float(spb[1]), float(spb[2]) 
                    self.total += null_len + 1 
                    if null_len > 0: effect_data.append((null_len))  
                    null_len = 0 
                    effect_data.append((round(-log(pv,10),1),round(effect,1))) 
            if null_len > 0: effect_data.append((null_len))             
            self.total += null_len
            self.effects[f_chr] = effect_data 
        
            gf.close()     
        if self.TESTS['NOSNPS']: 
            self.snp_file = snp_handle 
            w.close()
        return  


    ####################################   PHENOTYPES   #######################################
    

    def ph_get_genotypes(self, genotype_prefix): 
        self.genotype_prefix = genotype_prefix
        g_path, g_prefix = "/".join(self.genotype_prefix.split('/')[0:-1]), self.genotype_prefix.split('/')[-1] 
        
        CK = dd(int) 
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


    def ph_fill(self, genotype_prefix, phenotype_files):
        # cols  
        self.VALID, self.type, self.files, col_data = True, 'binary', phenotype_files, dd(list) 
        self.ph_get_genotypes(genotype_prefix) 
        self.X_fields = ['--test.data',phenotype_files[0]] 
        
        if phenotype_files[0] == phenotype_files[-1]: 
                self.files = [phenotype_files[0]] 
                self.X_fields.extend(['--valid.data','0']) 
        else:   self.X_fields.extend(['--valid.data',phenotype_files[1]]) 
        

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


