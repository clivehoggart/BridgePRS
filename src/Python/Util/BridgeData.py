import sys, os, gzip  
from collections import defaultdict as dd
from .BridgeProgress  import BridgeProgress


def bridge_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeDataError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                     '+es+'\n')
    else: sys.stderr.write('\nBridgeDataError: '+eString+'\n')
    sys.exit(2) 

def parse_error(eString, SAVE=False):
    if not SAVE: 
        sys.stderr.write('\nBridgeParseError: '+eString+'\n')
        sys.exit(2) 
    sys.stderr.write('BridgeParseError: '+eString+'\n')
    return  
            
def parse_warning(eString, eRule): 
    if eRule == 'inferred': sys.stderr.write('                         BridgeWarning: --'+eString+' not supplied by command line or config file (inferring argument from process of elimination)\n') 
    else:                   sys.stderr.write('                         BridgeWarning: '+eString+'\n') 
    
















class BridgeData: 
    def __init__(self, bp, name):
        self.check, self.paths, self.args, self.name = bp.CHECK, bp.paths, bp.args, name
        self.names, self.chromosomes, self.fields, self.map = 'NA', [], {}, dd(list) 
        self.VALID, self.BYCHR = False, False 
        
    def validate_and_collect_prefix_files(self, pf): 
        self.prefix, self.prefix_path , pName, pLen = pf, "/".join(pf.split('/')[0:-1]), pf.split('/')[-1], len(pf.split('/')[-1]) 
        self.pfiles = [self.prefix_path+'/'+x for x in os.listdir(self.prefix_path) if x[0:pLen] == pName] 


    def find_suffix_matches(self, suffix_list): 
        suffix_key  = dd(list) 
        for pn in self.pfiles: 
            p_tail = pn.split(self.prefix)[-1] 
            for sf in suffix_list: 
                if len(p_tail.split(sf)) == 2: 
                    p_chr, p_ext = [x if len(x) > 0 else 'NA' for x in p_tail.split(sf)] 
                    suffix_key[sf].append([p_chr, pn]) 
        return suffix_key 
    


    def stream_file_cols(self, fp, COLS = [0], HEADER = False): 
        s_data = [] 
        if fp.split('/')[-1] == 'gz': gf = gzip.open(fp, 'rt') 
        else:                         gf = open(fp, 'rt') 
        for line in gf: 
            line = line.split()
            s_data.append(",".join([line[c] for c in COLS]))
        gf.close() 
        return s_data 



    def get_file_info(self, fp, RULE = 'NA', KEY = 'NA', GZ=False): 
        if GZ:           gf = gzip.open(fp, 'rt') 
        else:            gf = open(fp, 'rt') 
        header = gf.readline().split() 
        if RULE.upper() == 'HEADER': RT = list(set(header)) 
        else: 
            K  = {i: [h] for i,h in enumerate(header)}  
            for line in gf: 
                for j,x in enumerate(line.split()): K[j].append(x) 
            if RULE.upper() == 'SUMSNPS': RT =  K[header.index(KEY)] 
            elif RULE.upper() ==   'BIM': RT = K[1] 
            elif RULE.upper() ==   'FAM': RT = K[0] 
            else:                         parse_error('fileGET') 
        gf.close() 
        return(RT) 
    


    
    def confirm_helper_file(self, fp, fd, fstr, RULE = 'NA'): 
        if fp is not None: return fp 
        w = open(self.paths['tmp']+'/'+fstr,'w') 
        for i,n in enumerate(fd): 
            if i % 2 == 0: w.write('%s %s\n' % (n,n)) 
        w.close() 
        return self.paths['tmp']+'/'+fstr 

                

    def bf_fill(self, prefix, id_file, CHRS = []):  
        self.prefix, self.id_file, self.VALID, suffix_list = prefix, id_file, True, ['.bed','.bim','.fam'] 
        self.validate_and_collect_prefix_files(prefix) 
        suffix_pairs     = self.find_suffix_matches(suffix_list) 
        for sf in suffix_list: 
            for f_chr, fp in suffix_pairs[sf]: self.map[f_chr].append(fp) 
         
        for k,ft in self.map.items(): 
            k_ext = [fx.split('.')[-1] for fx in ft] 
            k_path = list(set([".".join(fx.split('.')[0:-1]) for fx in ft]))
            if k_ext != ['bed','bim','fam']: parse_error('A bed, bim and fam file are required: '+ft[0].split('.')[0]) 
            elif len(k_path) != 1:           parse_error('A bed, bim and fam file are required: '+ft[0].split('.')[0]) 
            else:                            self.map[k] = k_path[0]  

        if self.id_file is None:                   self.id_file = self.bf_get_id_file() 
        if len([k for k in self.map.keys()]) == 1: self.map = {c: [v for v in self.map.values()][0] for c in CHRS} 
        else:                                      self.BYCHR = True 
        self.X_fields = ['--clump-field',vars(self.args)['ssf-p'], '--clump-snp-field',vars(self.args)['ssf-snpid']] 
        return self  




    def bf_get_id_file(self): 
        parse_warning('--id_file not supplied, attempting to use all ids','missing') 
        id_names = [] 
        for k,ft in self.map.items(): 
            id_names = list(set(id_names + self.stream_file_cols(ft+'.fam', [0,1]))) 
        id_names = sorted([id_names[i] for i in range(0,len(id_names),2)]) 
        id_file = self.paths['tmp']+'/ids.txt'
        w = open(id_file,'w') 
        for x in id_names:  w.write(x.split(',')[0]+' '+x.split(',')[1]+'\n') 
        w.close()
        return id_file 
        
        




    ####################################   SUMSTATS   #######################################

    # pfiles 

    def ss_fill(self, prefix, snp_file): 
        self.prefix, self.snp_file, self.VALID, self.BYCHR = prefix, snp_file, True, True 
        self.validate_and_collect_prefix_files(prefix)
        suffix_list     = self.ss_validate_suffix(self.args.sumstats_suffix) 
        suffix_pairs    = self.find_suffix_matches(suffix_list)[suffix_list[0]] 
        if len(suffix_pairs) < 2: suffix_pairs  = self.ss_split_suffix_matches(suffix_pairs[0][1])  
        self.map                                = self.ss_confirm_suffix_files(suffix_pairs) 
        self.suffix, self.chromosomes           = self.args.sumstats_suffix, sorted([k for k in self.map.keys()]) 
        sum_stats_fields    =  ['ssf-alt', 'ssf-beta', 'ssf-maf', 'ssf-p', 'ssf-ref', 'ssf-se', 'ssf-snpid', 'ssf-ss']
        self.fields         =  {ks.split('-')[-1].upper(): vars(self.args)[ks] for ks in sum_stats_fields}
        if self.args.verify or self.snp_file is None: self.ss_verify_fields() 
        if self.snp_file is None:                     self.snp_file = self.ss_get_snpfile() 
        self.X_fields   = ['--sumstats.allele0ID',self.fields['REF'],'--sumstats.allele1ID',self.fields['ALT'],'--sumstats.betaID',self.fields['BETA']]
        self.X_fields.extend(['--sumstats.frqID',self.fields['MAF'],'--sumstats.nID',self.fields['SS'],'--sumstats.seID',self.fields['SE'], '--sumstats.snpID',self.fields['SNPID']])
        return self 



    def ss_confirm_suffix_files(self, suffix_pairs): 
        M = {} 
        for f_chr, fp in suffix_pairs:
            if fp.split('.')[-1] != 'gz':
                os.system('gzip -f '+fp)  
                fp += '.gz' 
            M[f_chr] = fp  
        return M 


    def ss_get_snpfile(self): 
        parse_warning('--snp_file not supplied, attempting to use all snps','missing') 
        snp_names = [] 
        for f_chr, fp in self.map.items():  snp_names.extend(self.get_file_info(fp, RULE = 'SUMSNPS', KEY = self.fields['SNPID'], GZ = True)[1::])
        snp_file = self.paths['tmp']+'/snps.txt' 
        w = open(snp_file,'w') 
        for snp in snp_names:  w.write(snp+'\n') 
        w.close()
        return snp_file 
  



    #def split_suffix_matches(self, sumstat_file): 
    def ss_split_suffix_matches(self, sumstat_file): 
        prefix = self.prefix.split('/')[-1] 
        suffix = self.args.sumstats_suffix.split('.gz')[0]  
        new_path = self.paths['tmp']+'/'+prefix 
        self.prefix = new_path 
        if sumstat_file.split('.')[-1] == 'gz': gf = gzip.open(sumstat_file, 'rt') 
        else:                                   gf = open(sumstat_file, 'rt')  
        header_str  =  gf.readline().strip() 
        header_list = header_str.split() 
        k = 0 
        parse_warning('Single sumstats file supplied, splitting sumstats','Splitting') 
        while k < len(header_list): 
            if 'CHR' in header_list[k].upper(): break 
            k+=1 
        if k == len(header_list):  parse_error('Cannot Split File')  
        SP, W = [], {} 
        for line in gf: 
            line = line.strip() 
            lc  = line.split()[k] 
            if lc not in W: 
                W[lc] = open(new_path+lc+suffix,'w') 
                W[lc].write(header_str+'\n') 
            W[lc].write(line+'\n') 
        for lc in W:
            W[lc].close() 
            SP.append([lc, new_path+lc+suffix]) 
        return SP 
            
   













    def ss_validate_suffix(self, SUFFIX):
        suffix_cands = dd(int) 
        for pn in self.pfiles: 
            k, p_tail = 0, pn.split(self.prefix)[-1]  
            while p_tail[k] in ['1','2','3','4','5','6','7','8','9','0']: k+=1 
            suffix_cands[p_tail[k::].split('.gz')[0]] += 1  
        suffix_cands = [sc[0] for sc in sorted(suffix_cands.items(), key = lambda X: X[1], reverse=True)] 
        if SUFFIX is None: 
            parse_warning('sumstats_suffix','inferred') 
            self.args.sumstats_suffix = suffix_cands[0]+'.gz'   
            return suffix_cands[0:1] 
        elif SUFFIX.split('.gz')[0] in suffix_cands: 
            return [SUFFIX.split('.gz')[0]] 
        else: 
            SF = ['Invalid suffix supplied, --sumstats_suffix '+SUFFIX] 
            SF.append('Does not match file prefix: '+self.prefix) 
            if len(suffix_cands) > 0: SF.append('Consider a matching suffix: '+suffix_cands[0]) 
            bridge_error(SF) 
            return  
            
            
    def ss_verify_fields(self): 
        f_cnt = dd(int) 
        for f_chr, fp in self.map.items():  
            if fp.split('.')[-1] == 'gz': gf =      gzip.open(fp, 'rt') 
            else:                                   gf = open(fp, 'rt')  
            for x in gf.readline().split(): f_cnt[x] += 1
            gf.close()  
        cands, cand_errors, cand_found =  [x for x,y in f_cnt.items() if y == len(self.map)] , [], [] 
        for x,y in self.fields.items(): 
            if y in cands: cand_found.append(y) 
            else:          cand_errors.append([x,y]) 
        if len(cand_errors) == 0: return 
        
        rem_cands = [c for c in cands if c not in cand_found] 
        if len(rem_cands) == 0: rem_cands = cands 
        b_errors, b_missing = [], [] 
        for x,y in cand_errors: 
            if y is None: b_missing.append('--ssf-'+x.lower()) 
            else:         b_errors.append('--ssf-'+x.lower()+' '+y) 
        
        ss_error = [] 
        if len(b_errors) > 0: ss_error.extend(['Invalid Sumstats Field(s) Supplied:']+b_errors) 
        if len(b_missing) > 0: ss_error.extend(['Missing Sumstats Field(s):']+b_missing) 
        ss_error.append('Available Fields:  '+",".join(rem_cands)) 
        
        bridge_error(ss_error) 
        
        sys.exit() 
        for x,y in cand_errors: b_error.append('--ssf-'+x.lower()+' '+y) 
        b_error.append('Available Fields:  '+",".join(rem_cands)) 
        bridge_error(b_error) 
        return 



    
    
    ####################################   PHENOTYPES   #######################################
    
    
    
    
    def ph_fill(self, pheno_files):
        self.VALID = True
        self.files = pheno_files  
        self.X_fields = ['--test.data',pheno_files[0]] 
        if pheno_files[0] == pheno_files[-1]: 
            self.names = pheno_files[0] 
            self.X_fields.extend(['--valid.data','0']) 
        else:                                 
            self.names = pheno_files[0]+','+pheno_files[1] 
            self.X_fields.extend(['--valid.data',pheno_files[1]]) 
       
        for v in vars(self.args): 
            if v.split('-')[0] in ['pf']: 
                nv, kV = v.split('-')[-1].upper(), vars(self.args)[v] 
                if kV is not None: 
                    self.fields[nv] = kV 
                    if nv == 'NAME':         self.X_fields.extend(['--pheno.name',kV])
                    elif nv == 'COVARIATES': self.X_fields.extend(['--cov.names',kV]) 
        
        if 'NAME' not in self.fields: bridge_error('Phenotype field required (--pf-name)') 
        self.ph_verify_fields() 
        if self.type == 'binary': self.X_fields.extend(['--binary','1']) 
        return self 
        #if self.args.verify:      self.ph_verify_fields() 
            

    def ph_verify_fields(self): 
        p_types = []
        for i,fname in enumerate(self.files): 
            f = open(fname,'rt') 
            lp = f.readline().split() 
            if i == 0: cands = [x for x in lp] 
            else:      cands = [x for x in lp if x in cands]
            for j,c in enumerate(lp):
                if c == self.fields['NAME']: 
                    nf = j 
                    break 
            pts = [] 
            for k,line in enumerate(f): 
                if k < 100: pts.append(line.split()[nf])
            pts = list(set(pts)) 
            if len(pts) == 2: p_types.append('binary') 
            else:             p_types.append('cont') 
            f.close() 
        
        p_types = list(set(p_types)) 
        if len(p_types) == 1: self.type = p_types[0] 
        else:                 bridge_error('ambiguous phenotype') 
        b_error = ['Invalid phenotype field name(s) supplied'] 
        for x,y in self.fields.items(): 
            if x == 'NAME' and y not in cands: bridge_error(b_error+['--pf-name '+y,'Available Fields: '+",".join(cands)]) 
            elif x[0:3] == 'COV': 
                print(y) 
                my_cov    = [c.strip() for c in y.split(',')] 
                fail_covs = [c for c in my_cov if c not in cands] 
                if len(fail_covs) > 0: bridge_error(b_error+['--pf-covariates '+y,'Available Fields: '+",".join(cands)]) 
        return             











































































