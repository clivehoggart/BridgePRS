import sys, os, gzip 
from collections import defaultdict as dd
from collections import Counter as cc 
from . import PopTools as ptools    














class SumStats: 
    def __init__(self, pop, mps): 

        self.VALID, self.TESTS, self.map, self.total = False, dd(bool), dd(list), 0 
        self.mps, self.pop, self.fields, self.debug_level, self.savepath = mps, pop, pop.field_key, pop.args.debug_level, pop.paths['save']   
        self.source_prefix, self.source_suffix, self.pop_size = pop.sumstats_prefix, pop.sumstats_suffix, pop.sumstats_size
        self.snp_file, self.max_clump_size, self.thinned_snp_file = pop.snp_file, pop.max_clump_size, pop.thinned_snp_file      

        
        if pop.sumstats_size is None and pop.args.module in ['pipeline', 'build-model']: self.mps.error('SumstatsError: Missing Argument (SUMSTATS_SIZE=), Required For This Module') 

        if self.thinned_snp_file is None:     self.thin_snps = '0' 
        if self.max_clump_size   is None: self.max_clump_size = '0' 
 

    def load(self, genoPheno, prevPop): 
        self.load_rs_key(genoPheno, prevPop) 
        self.load_sumstats_files() 
        self.X_fields   = ['--sumstats.snpID',self.fields['ID'],'--sumstats.allele0ID',self.fields['REF'],'--sumstats.allele1ID',self.fields['ALT'],'--sumstats.betaID',self.fields['BETA']]
        return self 
        





    ####################################   RS_KEY   #######################################

    def load_rs_key(self, genoPheno, PP): 
        self.CK, self.rs_key, self.WT, self.VALID, self.BYCHR = dd(int), {}, [], True, True
        


        if PP: self.TESTS['NOSNPS'], self.snp_file, self.thin_snps = PP.sumstats.TESTS['NOSNPS'], PP.sumstats.snp_file, PP.sumstats.thin_snps         
        elif self.snp_file is None: 
            self.snp_file   = self.savepath+'/snps.'+self.pop.name.lower()+'_valid.txt'
            self.snp_handle, self.TESTS['NEWSNPS'], self.TESTS['NOSNPS'] = open(self.snp_file,'w'), True, True 
        if self.debug_level == 0: return 
        
        geno_path, geno_name = '/'.join(genoPheno.genotype_prefix.split('/')[0:-1]), genoPheno.genotype_prefix.split('/')[-1] 
        geno_files = [f for f in os.listdir(geno_path) if f[0:len(geno_name)] == geno_name and f.split('.')[-1]=='bim']         
        
        if not PP: self.begin_rs_key(geno_path, geno_files) 
        else:      self.continue_rs_key(geno_path, geno_files, PP) 

    def begin_rs_key(self, geno_path, geno_files): 
        if not self.TESTS['NOSNPS']: 
            with open(self.snp_file) as f: self.rs_key = {line.strip(): [] for line in f} 

        for gf in geno_files: 
            with open(geno_path+'/'+gf) as f: 
                for line in f: 
                    line = line.split() 
                    try: chr_name, rs, loc, ref, alt = int(line[0]), line[1], int(line[3]), line[4].upper(), line[5].upper()  
                    except: self.mps.error('Invalid genotype file: '+gf) 
                    if self.TESTS['NOSNPS']: self.rs_key[rs] = [chr_name, loc, ref, alt] 
                    elif rs in self.rs_key:  self.rs_key[rs] = [chr_name, loc, ref, alt]  
                    else:                    continue  
        self.snp_list_len = len(self.rs_key) 
        self.rs_key = {a:b for a,b in self.rs_key.items() if len(b) > 0} 
        self.genome_snps = len(self.rs_key)
        return 


    def continue_rs_key(self,geno_path, geno_files, PP): 
        self.snp_list_len = PP.sumstats.snp_list_len 
        for gf in geno_files: 
            with open(geno_path+'/'+gf) as f: 
                for line in f: 
                    line = line.split() 
                    try: chr_name, rs, loc, ref, alt = int(line[0]), line[1], int(line[3]), line[4], line[5] 
                    except: bridge_pop_error('Invalid genotype file: '+gf) 
                    if rs not in PP.sumstats.rs_key or len(PP.sumstats.rs_key[rs]) < 5: self.CK['GENO_EXTRA'] += 1 
                    else:
                        rC, rL, rRef, rAlt, rBool = PP.sumstats.rs_key[rs] 
                        if chr_name != rC or loc != rL: 
                            azd = str(rC)+':'+str(rL)+', '+str(chr_name)+':'+str(loc) 
                            if   chr_name != rC: self.mps.error('Invalid Genome Builds For Base/Target Genotype, '+rs+' On Multiple Chromosomes: '+azd) 
                            else:                self.mps.error('Invalid Genome Builds For Base/Target Genotype, '+rs+' At Multiple Locations: '+azd) 
                        else:
                            self.CK['GENO_'+ptools.compare_alleles([ref,alt],[rRef, rAlt])] += 1 
                            self.rs_key[rs] = [rC, rL, rRef, rAlt] 
        self.genome_snps = len(self.rs_key)+self.CK['GENO_EXTRA']
        return 


    ####################################   SUMSTATS  OVERHEAD  ####################################### 

    

    def find_sumstats_pairs(self): 
        if self.source_suffix == '__FILE__': return [[self.source_suffix, self.source_prefix]] 
        prefix_path, prefix_name = '/'.join(self.source_prefix.split('/')[0:-1]), self.source_prefix.split('/')[-1] 
        file_cands = [f for f in os.listdir(prefix_path) if f[0:len(prefix_name)] == prefix_name] 
        if self.source_suffix: file_cands = [pf for pf in file_cands if pf[-1*len(self.source_suffix)::] == self.source_suffix] 
        if len(file_cands) == 1: 
            self.source_suffix = '__FILE__' 
            return [[self.source_suffix, prefix_path+'/'+file_cands[0]]] 
        chr_cands = [fc.split(prefix_name)[-1] for fc in file_cands] 
        if self.source_suffix: 
            try: 
                chr_cands = [cx.split(self.source_suffix)[0] for cx in chr_cands] 
                chr_ints =  [int(cx) for cx in chr_cands]
                if len(list(set(chr_ints))) == len(file_cands) and len(file_cands) < 25: return [[cr,prefix_path+'/'+cf]  for cr,cf in zip(chr_ints, file_cands)] 
            except: pass
            p_warn = ['Supplied Sumstats Prefix/Suffix For Pop "'+self.pop.name+'" Does Not Produce Unique Numerical Chromsomes',' Attempting To Infer Correct Prefix/Suffix...'] 
        else: p_warn = ['Sumstats Suffix Not Supplied for Pop "'+self.pop.name+'"',' Attempting to Infer Correct Prefix/Suffix...'] 
        try: 
            new_prefix, new_suffix = ptools.get_prefix_suffix(file_cands)             
            chr_ints = [int(pf.split(new_suffix)[0].split(new_prefix)[-1]) for pf in file_cands] 
            self.source_prefix, self.source_suffix        = prefix_path+'/'+new_prefix, new_suffix 
            if len(list(set(chr_ints))) == len(file_cands) and len(file_cands) < 25: 
                p_warn[-1]+= '...........FOUND! Using Valid Prefix|Suffix Pair: '+new_prefix+'|'+new_suffix
                ptools.warn(p_warn) 
                return [[cr,prefix_path+'/'+cf]  for cr,cf in zip(chr_ints, file_cands)] 
        except: pass 
        
        ptools.warn([p_warn[0],p_warn[1]+'...FAILED']) 

        self.mps.error('ConfigFileError: Invalid Sumstats Prefix|Suffix Pair: '+prefix_name+'|'+str(self.source_suffix)) 

    def load_sumstats_files(self): 
        self.s_key  = {} 
        sumstats_pairs = self.find_sumstats_pairs() 
        
        new_path = self.savepath+'/sumstats' 
        if not os.path.exists(new_path): os.makedirs(new_path) 
        self.prefix, self.suffix = new_path+'/ss.'+self.pop.name+'.', '.out.gz' 

        
        if len(sumstats_pairs) == 1: 
            self.TESTS['INFER_CHR'] = True 
            ptools.warn('Sumstats files for pop '+self.pop.name+' are not separated by chromosome, attempting to split file in '+new_path) 
            p_handle = self.start_sumstats_file(sumstats_pairs[0][1]) 
            if self.TESTS['INFER_CHR'] and self.debug_level < 1:  self.mps.error('BridgeSumstatsError: No CHR label in sumstats file, Please Increase Debug Level (--debug_level 2) Or Supply Chromosome Info in Sumstats File')           
            elif self.debug_level == 0: self.split_sumstats_file(p_handle) 
            else:                       self.process_sumstats_file(p_handle) 
            p_handle.close() 
       
        else: 
            for c_name, f_name in sumstats_pairs: 
                if self.debug_level == 0: self.map[str(c_name)] = f_name 
                else:
                    p_handle = self.start_sumstats_file(f_name)
                    self.process_sumstats_file(p_handle) 
                    p_handle.close() 
       
        for sc,sd in self.s_key.items(): 
            sd[1].close() 
            os.system('gzip -f '+sd[0]) 
            self.map[str(sc)] = sd[0]+'.gz' 
        return 







    ##############################   PROCESSING SUMSTATS FILE  #######################################
    
    def start_sumstats_file(self, pf, REQ=['ID','REF','ALT','P','BETA']): 
        f_header, f_reverse, f_extra = ['CHR'] + [self.fields[x] for x in REQ], {self.fields[x]: x for x in REQ}, [] 
        p_handle, p_header = ptools.zip_open(pf, WITH_HEADER=True) 
        self.header, self.loc_key = '\t'.join(f_header), {} 
        for i,h in enumerate(p_header): 
            if self.TESTS['INFER_CHR'] and 'CHR' in h.upper(): 
                self.loc_key['CHR'] = i 
                self.TESTS['INFER_CHR'] = False 
            elif h in f_header[1::]: 
                try: self.loc_key[f_reverse.pop(h)] = i 
                except: self.mps.error('SumstatsError: Repeated Column '+h) 
            else: f_extra.append(h) 

        if len(f_reverse) > 0:  self.mps.error('Invalid Sumstats File: '+pf,'         Missing Fields: '+','.join([x for x in f_reverse.values()])+'\n         Available Choices: '+",".join([str(x) for x in f_extra]))
        
        self.locs = [self.loc_key[r] for r in REQ] 
        return p_handle
    
    
    def add_sumstats_line(self,fC,snpID,LD): 
        self.total += 1 
        if fC not in self.s_key: 
            self.s_key[fC] = [self.prefix+str(fC)+'.out', open(self.prefix+str(fC)+'.out','w')]
            self.s_key[fC][1].write(self.header+'\n')
        self.s_key[fC][1].write("\t".join(LD)+'\n')  
        if self.TESTS['NOSNPS']: self.snp_handle.write(snpID+'\n')
        return


    def split_sumstats_file(self, p_handle): 
        for li, line in enumerate(p_handle): 
            lp = line.split() 
            xC, xS, lData   = lp[self.loc_key['CHR']], lp[self.loc_key['ID']], [lp[l] for l in self.locs]
            LD = [xC] + lData
            if 'NA' in LD: self.CK['NA_LINE'] += 1 
            else:          self.add_sumstats_line(xC, xS, [xC]+lData) 


    def process_sumstats_file(self, p_handle): 
        for li,line in enumerate(p_handle): 
            lp = line.split() 
            LD = [lp[j] if j != 'NA' else j for j in self.locs]
            if 'NA' in LD: self.CK['NA_LINE'] += 1 
            elif LD[0] not in self.rs_key: self.CK['NO_GENOTYPE'] +=1 
            else: 
                rsData = self.rs_key[LD[0]] 
                rsChr, rsLoc, rsRef, rsAlt = rsData[0], rsData[1], rsData[2], rsData[3] 
                LD[1], LD[2], lpWt = LD[1].upper(), LD[2].upper(), float(LD[-1]) 
                relationship = ptools.compare_alleles([LD[1], LD[2]],[rsRef,rsAlt]) 
                self.CK[relationship] += 1 
                if relationship != 'INVALID': 
                    if len(self.WT) < 2 or lpWt < min(self.WT) or lpWt > max(self.WT):  self.WT.append(lpWt) 
                    self.add_sumstats_line(rsChr, LD[0], [str(rsChr)]+LD) 
                    self.rs_key[LD[0]].append(True) 
       
            cMatch, cInvalid, cSwap, cRev, cTotal = self.CK['MATCH'], self.CK['INVALID'], self.CK['SWAPREF'], self.CK['REVCOMP'], sum([kk for kk in self.CK.values()]) 
            if cInvalid == cTotal: bridge_sumstats_error(['Sumstats File Error(s), No Matching Ref/Alt Bases']) 
            elif cInvalid > cTotal/2.0: bridge_sumstats_error(['Sumstats File Error(s), Majority MisMatching Ref/Alt Bases - Please Check Builds']) 
        p_handle.close() 

     
