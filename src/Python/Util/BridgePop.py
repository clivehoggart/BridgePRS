import sys, os, gzip  
from collections import defaultdict as dd
from .BridgeProgress  import BridgeProgress
from math import log 

#Missing Chromosomes 


def bridge_debug_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('BridgeDebugError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('BridgeDebugError: '+eString+'\n')
    sys.exit(2) 



def bridge_error(eString):
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
        


REV_COMP =  {'A': 'T', 'C': 'G','G': 'C', 'T': 'A'}
NUM_STRS =  ['1','2','3','4','5','6','7','8','9','0'] 



class BridgePop: 
    def __init__(self, args, progress, paths, pop_name, pop_key, ld_ref, pop_type = 'TARGET', prevPop = None):
        
        self.args, self.progress, self.paths, self.name, self.key, self.type  = args, progress, paths, pop_name, pop_key, pop_type
        self.ref_pop = self.key['ldpop'] 
        
        self.bdata      = BridgeData(self.args, self.progress, 'bdata',      paths, pop_name)
        self.bdata.add_panel(ld_ref[self.ref_pop]) 
        

        self.genopheno = BridgeData(self.args, self.progress, 'phenotypes', paths, pop_name, pop_type) 
        if   not pop_key['genotype_prefix']: bridge_error('At least one Genotype Prefix is Required on the Command Line (--genotype_prefix) or in a config_file (GENOTYPE_PREFIX=)')  
        elif not pop_key['phenotype_file']:  bridge_error('At least one Phenotype File is Required on the Command Line (--phenotype_file) or in a config_file (PHENOTYPE_FILE=)')  
        else:                                self.genopheno.ph_fill(pop_key['genotype_prefix'],pop_key['phenotype_file'], self.args.validation_file)         
        

        self.sumstats   = BridgeData(self.args, self.progress, 'sumstats',   paths, pop_name, pop_type) 
        if not pop_key['sumstats_prefix']: bridge_error('A Sumstats Prefix is Required on the Command Line (--sumstats_prefix) or in a config file (SUMSTATS_PREFIX=) for pop '+pop_name)  
        elif not prevPop:                  self.sumstats.add_sumstats(pop_key, self.genopheno, None)
        else:                              self.sumstats.add_sumstats(pop_key, self.genopheno, prevPop.sumstats) 
        
        self.verify_chromosomes() 
        if self.args.debug_level == 2: 
            if prevPop:                             self.verify_ss_data(prevPop.sumstats, self.sumstats) 
            elif self.args.module != 'easyrun':     self.verify_ss_data(self.sumstats, None) 
        
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
        if len(mInputs) > 0: bridge_error(['Missing Input Data:']+mInputs) 
        if len(mGen) > 0:    bridge_error(['Missing Run Data:']+mGen) 
        return True  
        

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
            if len(b_missing) > 0: bridge_error(['Missing Chromosomes','The following chromosomes are found in the ld-panel but not in sumstats:',",".join(b_missing)])  
            s_missing = [c for c in self.bdata.map.keys() if c not in self.sumstats.map.keys()]
            if len(s_missing) > 0: bridge_error(['Missing Chromosomes','The following chromosomes are found in sumstats but not in ld-panel:',",".join(b_missing)])  
            bridge_error('Missing Chromosomes') 

        return 


    def verify_ss_data(self, ss1, ss2): 
        self.progress.printout('\nBridgeDebugData(Lvl'+str(self.args.debug_level)+'):')
        pTypes = [] 
        for si,s in enumerate([ss1, ss2]): 
            if s is not None: 
                pk, pn, pt = s.pop_type, s.pop_name, s.phenoTYPE, 
                eR = min(s.WT),max(s.WT)
                self.progress.printout('\nPOP'+str(si+1)+' '+pn+' ('+pk+')')
                try: 
                    self.progress.printout('Phenotype: '+self.args.phenotype+' ('+pt+')') 
                    eStr = 'Effect Size Range: '+str(eR[0])+', '+str(eR[1]) 
                    self.progress.printout(eStr) 
                    if pt == 'binary' and eR[0] > 0:  bridge_debug_error('Truncated Effect Size Range, log(odds) is required.') 
                except: 
                    self.progress.printout('Phenotype: '+str("Not Declared")) 
                
                

                g1, g2 = [], [] 
                for nk,n in [['MATCH','Matching'],['SWAPREF','Ref/Alt Swap'],['REVCOMP','Reverse Complement'],['INVALID','Invalid Bases']]: 
                    if s.CK[nk] > 0 :        g1.append(str(s.CK[nk])+' ('+n+')') 
                    if s.CK['GENO_'+nk] > 0: g2.append(str(s.CK['GENO_'+nk])+' ('+n+')') 

                
                if si == 0: 
                    self.progress.printout('Initial Genotype Variants: '+str(s.snp_list_len)+' (SNP LIST), '+str(s.genome_snps)+' (VALID GENOTYPED SNPS)') 
                    self.progress.printout('Shared GWAS Variants:      '+', '.join(g1))  
                else: 
                    self.progress.printout('Shared Genotype Variants: '+','.join(g2)) 
                    self.progress.printout('Shared GWAS Variants:     '+', '.join(g1)) 

        return 







class BridgeData: 
    def __init__(self, args, progress, datatype, paths, pop_name,pop_type = 'TARGET'):
        self.args, self.progress, self.datatype, self.paths, self.pop_name, self.pop_type = args, progress, datatype, paths, pop_name, pop_type
        self.VALID, self.TESTS, self.map = False, dd(bool), dd(list) 



        if self.datatype == 'sumstats':       self.fields, self.map, self.total = {}, dd(list), 0 
        elif self.datatype == 'phenotypes':   self.fields, self.map = {}, dd(list) 
        
        


    
    def add_panel(self, ld_quad): 
        self.VALID = True 
        self.id_file, self.prefix, self.BYCHR, b_key  = ld_quad
        self.ldpath = "/".join(self.prefix.split('/')[0:-1])
        for c in b_key: self.map[c] = b_key[c] 
        self.X_fields = ['--clump-field',vars(self.args)['ssf-p'], '--clump-snp-field',vars(self.args)['ssf-snpid']] 
        return self 



    

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
        FK, NK, IK = {b: a for a,b in self.fields.items()}, {'CHR': 'CHR'}, {} 
        p_handle = zip_open(pf) 
        p_init   = p_handle.readline()  
        p_header = p_init.split()   
        for i,h in enumerate(p_header): 
            if   'CHR' in h.upper(): lk = 'CHR' 
            elif h.upper() in FK:    lk = FK[h.upper()] 
            else:                    continue 
            if lk in IK: bridge_sumstats_error(['Repeated Column']) 
            IK[lk], NK[lk] = i, h     
        for r in required: 
            if r not in IK: 
                if r != 'CHR':                  bridge_sumstats_error(['Invalid Sumstats File: '+pf,'    Missing Header Field: --SSF-'+r+' '+self.fields[r]]) 
                elif self.args.debug_level < 2: bridge_sumstats_error(['Cannot Split Sumstats, CHR field not found','    Fields Found: '+",".join(p_header)]) 
                else:                           IK['CHR'] = 'NA'  
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
        if self.TESTS['NOSNPS']: self.snp_handle.write(snpID+'\n')
        return


    def split_sumstats(self, p_file, CHR = 'NA'):
        self.s_key = {} 
        p_handle, self.header, locs, IK, iC, iS = self.begin_sumstats_file(p_file) 
        if self.args.debug_level < 2: 
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
                try: rsChr, rsLoc, rsRef, rsAlt = self.rs_key[lp[iS]] 
                except ValueError: 
                    print(self.rs_key[lp[iS]]) 
                    bridge_debug_error('this is Z info:'+','.join(self.rs_key[lp[iS]])) 

                LD = [lp[j] if j != 'NA' else j for j in locs] 
                chr_cands = list(set([c for c in [str(LD[0]), str(CHR), str(rsChr)] if c != 'NA'])) 
                if len(chr_cands) > 1: bridge_sumstats_error('Ambiguous Chromosomes For '+lp[iS]+': '+str(LD[0])+','+str(rsChr)+','+str(CHR)+' (Sumstats, Genotype, Filename(s))') 
                
                try: fC, LD[0] = int(chr_cands[0]), chr_cands[0]
                except ValueError: bridge_sumstats_error(['Nonnumerical Chromosome ('+chr_cands[0]+')','    Sumstats File: '+p_file])  
                try:                lpRef, lpAlt, lpMaf, lpN, lpWt, lpSE, lpP = LD[2], LD[3], float(LD[4]), float(LD[5]), float(LD[6]), float(LD[7]), float(LD[8])
                except ValueError:  bridge_sumstats_error(['Sumstats File Error(s), Incorrect DataType:',header,'\t'.join(LD)]) 
                
                #lpRef = "AC" 
                relationship = self.base_comp([lpRef, lpAlt], [rsRef, rsAlt])
                self.CK[relationship] += 1 
                if relationship != 'INVALID': 
                    if len(self.WT) < 2 or lpWt < min(self.WT) or lpWt > max(self.WT):  self.WT.append(lpWt) 
                    self.add_sumstats_line(fC, LD, lp[iS]) 
                    self.rs_key[lp[iS]].append(True) 

                #print(relationship,'yo') 

                #self.CK[self.base_comp([lpRef, lpAlt], [rsRef, rsAlt])] += 1 
                

                #if len(self.WT) < 2 or lpWt < min(self.WT) or lpWt > max(self.WT):  self.WT.append(lpWt) 
                #self.add_sumstats_line(fC, LD, lp[iS]) 
                #self.rs_key[lp[iS]].append(True) 
                 
        
        p_handle.close() 
        for sc,sd in self.s_key.items(): 
            sd[1].close() 
            os.system('gzip -f '+sd[0]) 
            self.map[str(sc)] = sd[0]+'.gz' 
        return 



    def initialize_sumstats(self,genoPheno): 
        if self.args.thinned_snp_file is None: self.thin_snps = '0' 
        else:                                  self.thin_snps = self.args.thinned_snp_file 
        self.rs_key, self.snp_file = {}, self.args.snp_file 
        if self.args.snp_file is None: 
            self.TESTS['NOSNPS'], self.snp_file, self.snp_handle  = True, self.paths['save']+'/snps.txt', open(self.paths['save']+'/snps.txt','w') 
        
        if self.args.debug_level >= 2: 
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
                        except: bridge_error('Invalid genotype file: '+gf) 
                        if self.TESTS['NOSNPS']: self.rs_key[rs] = [chr_name, loc, ref, alt] 
                        elif rs in self.rs_key:  self.rs_key[rs] = [chr_name, loc, ref, alt]  
                        else:                    continue  
        
            self.snp_list_len = len(self.rs_key) 
            self.rs_key = {a:b for a,b in self.rs_key.items() if len(b) > 0} 
            self.genome_snps = len(self.rs_key)


    def continue_sumstats(self,genoPheno, prevPop): 
        self.rs_key, self.snp_file, self.thin_snps = {}, prevPop.snp_file, prevPop.thin_snps 
        if self.args.debug_level >= 2:
            self.snp_list_len = prevPop.snp_list_len 
            geno_path, geno_name = '/'.join(genoPheno.genotype_prefix.split('/')[0:-1]), genoPheno.genotype_prefix.split('/')[-1] 
            geno_files = [f for f in os.listdir(geno_path) if f[0:len(geno_name)] == geno_name and f.split('.')[-1]=='bim'] 
            for gf in geno_files: 
                with open(geno_path+'/'+gf) as f: 
                    for line in f: 
                        line = line.split() 
                        try: chr_name, rs, loc, ref, alt = int(line[0]), line[1], int(line[3]), line[4], line[5] 
                        except: bridge_error('Invalid genotype file: '+gf) 
                        
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
        self.phenoTYPE      = genoPheno.type
        sum_stats_fields    =  ['ssf-alt', 'ssf-beta', 'ssf-maf', 'ssf-p', 'ssf-ref', 'ssf-se', 'ssf-snpid', 'ssf-n']
        sum_stats_args      =  [vars(self.args)[ks] for ks in sum_stats_fields] 
        self.fields         =  {a.split('-')[-1].upper(): b for a,b in zip(sum_stats_fields, sum_stats_args)}
        

         
        self.CK, self.WT, self.VALID, self.BYCHR = dd(int), [], True, True
        if not prevPop: self.initialize_sumstats(genoPheno) 
        else:           self.continue_sumstats(genoPheno, prevPop) 

           


        self.source_prefix, self.source_suffix = pk['sumstats_prefix'], pk['sumstats_suffix'] 
        prefix_path, prefix_name = '/'.join(self.source_prefix.split('/')[0:-1]), self.source_prefix.split('/')[-1] 
        prefix_files = [f for f in os.listdir(prefix_path) if f[0:len(prefix_name)] == prefix_name] 
        
        if   len(prefix_files) == 0: bridge_sumstats_error('Invalid Sumstats Prefix: '+prefix) 
        elif self.args.debug_level == 0 and len(prefix_files) > 1: self.prefix, self.suffix = self.source_prefix, self.source_suffix 
        else: 
            new_prefix_path = self.paths['save']+'/sumstats'
            if not os.path.exists(new_prefix_path): os.makedirs(new_prefix_path) 
            self.prefix, self.suffix = new_prefix_path +'/ss.'+self.pop_name+'.', '.out.gz'


        if len(prefix_files) == 1: 
            bridge_sumstats_warning('Sumstats files for pop '+self.pop_name+' are not separated by chromosome, attempting to split file in '+new_prefix_path) 
            self.source_suffix = 'NA' 
            self.split_sumstats(prefix_path+'/'+prefix_files[0])  
        else:     
            if not self.source_suffix:
                if self.args.debug_level == 0: bridge_sumstats_error('Increase Debug Level (--debug_level) or supply sumstats suffix for pop '+self.pop_name+' (on the Command Line (--sumstats_suffix) or in a config_file (SUMSTATS_SUFFIX=))')  
                suffix_cands = [] 
                for p in prefix_files: 
                    k, c_cand, NEED_CHR = 0, p.split(prefix_name)[-1], False 
                    while k < len(c_cand) and c_cand[k] in NUM_STRS: k+=1  
                    suffix_cands.append(c_cand[k::]) 
                if len(list(set(suffix_cands))) != 1:  bridge_sumstats_error('Ambiguous Sumstats Suffixes ('+','.join(list(set(suffix_cands)))+'), please supply sumstats suffix for pop '+self.pop_name+' (on the Command Line (--sumstats_suffix) or in a config_file (SUMSTATS_SUFFIX=))')  
                self.TESTS['INFER_SUFFIX'] = True 
                self.source_suffix = suffix_cands[0] 
            for p in prefix_files: 
                try: chr_num = int(p.split(self.source_suffix)[0].split(prefix_name)[-1])
                except ValueError:  bridge_sumstats_error(['Nonnumerical Chromosome ('+c_cand.split(self.source_suffix)[0]+')','    Sumstats File: '+p,'    Consider Prefix, Suffix Combination: '+prefix_name+' '+self.source_suffix]) 
                if self.args.debug_level > 0: self.split_sumstats(prefix_path+'/'+p, chr_num) 
                else:                         self.map[str(chr_num)] = prefix_path+'/'+p 

        if self.TESTS['NOSNPS']: self.snp_handle.close() 
        self.X_fields   = ['--sumstats.allele0ID',self.fields['REF'],'--sumstats.allele1ID',self.fields['ALT'],'--sumstats.betaID',self.fields['BETA']]
        self.X_fields.extend(['--sumstats.frqID',self.fields['MAF'],'--sumstats.nID',self.fields['N'],'--sumstats.seID',self.fields['SE'], '--sumstats.snpID',self.fields['SNPID']])
        
        return















    ####################################   PHENOTYPES   #######################################
    

    def ph_get_genotypes(self, genotype_prefix): 
        self.genotype_prefix = genotype_prefix
        CK, g_path, g_prefix = dd(int), "/".join(self.genotype_prefix.split('/')[0:-1]), self.genotype_prefix.split('/')[-1] 
        for f in [x for x in os.listdir(g_path) if x[0:len(g_prefix)] == g_prefix]:
            if f.split('.')[-1] in ['bed','bim','fam']: CK[".".join(f.split('.')[0:-1])]+=1                                                                                                                                                                                   
        cands = [k for k,v in CK.items() if v == 3]
        if len(cands) == 0: bridge_error('Invalid genotype data prefix '+self.genotype_prefix) 
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
            with open(fn, 'rt') as f: 
                lp = f.readline().split() 
                if i == 0: self.header, lz = [x for x in lp] , ",".join(lp) 
                elif ",".join(lp) != lz: bridge_error('Phenotype File Headers Do Not Match: '+lz+' AND '+",".join(lp))
                for k,line in enumerate(f):
                    line = line.split() 
                    for j,c in enumerate(self.header): col_data[c].append(line[j]) 
                    if k > 100: break 

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


