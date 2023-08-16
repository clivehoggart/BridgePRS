import sys, os 
from collections import defaultdict as dd

def bridge_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeSettingsError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeSettingsError: '+eString+'\n')
    sys.exit(2) 

def parse_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeParseError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeParseError: '+eString+'\n')
    sys.exit(2) 






class BridgeSettings: 
    def __init__(self,io):
        self.args, self.io, self.module, self.cmd = io.args, io, io.args.module, io.args.cmd
        self.config   = LOAD_CONFIGS(self.args.pop_config+self.args.spec, io.bridgedir)

    def verify(self): 
        self.files, self.prefixes, self.maps, self.fields, self.errors, self.warnings = BridgeParse(self).load_input(self.io.pipeline.input_key, self.config)
        if 'sumstats' in self.maps: self.chromosomes = sorted([k for k in self.maps['sumstats'].keys()])
        else:                       self.chromosomes = [] 
        if len(self.errors) > 0: 
            self.report_missing(self.errors['missing']) 
            self.report_invalid(self.errors['invalid'])  
            sys.stderr.write('\nFor help with this module type: ./bridgePRS '+self.args.module+'\n')            
            sys.exit()
        if len(self.warnings) > 0: self.report_warnings(self.warnings['invalid'], self.warnings['recovered'], self.warnings['missing']) 
        return self 


    def report_warnings(self, invalids, recovered, missing): 
        for k,k_type,k_found in invalids:
            if k_type == 'prefix':   sys.stderr.write('BridgeWarning: ---'+k+' '+k_found+' is an invalid file prefix (Prefix not required for this job and will be ignored)\n') 
            else:                    sys.stderr.write('BridgeWarning: ---'+k+' '+k_found+' is an invalid file name   (File not required for this job and will be ignored)\n') 
        for k,k_type,k_found in recovered:
            sys.stderr.write('BridgeWarning: --'+k+'_'+k_type+' not supplied by command line or config file (using file generated from previous run: '+self.io.pipeline.progress_file+')\n') 
            
        for k,k_type,k_found in missing:
            sys.stderr.write('BridgeWarning: --'+k+' not supplied by command line or config file (inferring argument from process of elimination)\n') 
        return 


    def report_invalid(self, invalids): 
        for k,k_type,k_found in invalids:
            if k_type == 'prefix':   sys.stderr.write('BridgeIOError: ---'+k+' '+k_found+' is an invalid file prefix\n') 
            else:                    sys.stderr.write('BridgeIOError: ---'+k+' '+k_found+' is an invalid file name\n') 
        return  

    def report_missing(self, missing): 
        for k,k_type,k_found in missing:
            hints = [] 
            if k == 'clump_prefix': 
                c_dir = self.io.paths['run']+'/clump'
                hints = ['clump files should be located at: '+c_dir+' (maybe call ./bridgePRS single clump?)']
                try:    c_cand = os.listdir(c_dir)[0].split('/')[-1].split('_')[0].split('-')[0].split('.')[0] 
                except: c_cand = None
                if c_cand == self.args.popname: hints = ['Try: --clump_prefix '+c_dir+'/'+c_cand] 
            elif k == 'model_prefix': 
                if self.module == 'port': c_dir = self.io.paths['run']+'/bridge'
                else:                      c_dir = self.io.paths['run']+'/model'
                hints = ['model files should be located at: '+c_dir+' (maybe call ./bridgePRS single model?)']
                try:    c_cand = os.listdir(c_dir)[0].split('/')[-1].split('-')[0].split('.')[0] 
                except: c_cand = None
                if self.args.popname in c_cand:  hints = ['Try: --model_prefix '+c_dir+'/'+c_cand] 
            elif k == 'eval_file': 
                c_dir = self.io.paths['run']+'/eval'
                hints = ['eval files should be located at: '+c_dir+' (maybe call ./bridgePRS single eval?)'] 
                try:    found = [f for f in os.listdir(c_dir) if self.args.popname in f and 'best_model_params' in f][0] 
                except: found = 'NA' 
                if found != 'NA': hints = ['Try: --eval_file '+c_dir+'/'+found] 
            sys.stderr.write('BridgeIOError: '+k+' is missing\n') 
            for h in hints: sys.stderr.write('   BridgeHint: '+h+'\n') 
        return  
    
        
        


class BridgeParse:
    def __init__(self,settings):
        self.args = settings.args 
        self.required = LOAD_REQUIREMENTS(settings.module, settings.cmd) 
        self.f_key, self.map, self.errors, self.warnings = dd(lambda: dd(lambda: None)), dd(bool), dd(list), dd(list) 
        self.fields = dd(lambda: dd(bool))  
    
    
    def load_input(self, pipeline_key, config_key):
        self.parse_filetypes_and_arguments(pipeline_key, config_key) 
        self.verify_prefixes(config_key) 
        self.verify_files(config_key)      
        return self.f_key['file'], self.f_key['prefix'], self.map, self.fields, self.errors, self.warnings 
        

        
    
    def parse_filetypes_and_arguments(self, pipeline_key, config_key): 
        for v in vars(self.args): 
            v_name, v_type = v.split('_')[0], v.split('_')[-1] 
            kV = vars(self.args)[v] 
            
            if v_type in ['file','prefix']: 
                if   kV in [None,'0'] and config_key[v]: kV = config_key[v] 
                elif kV in [None,'0'] and pipeline_key[v_type][v_name]: 
                    kV = pipeline_key[v_type][v_name] 
                    if v_name in self.required:  self.warnings['recovered'].append([v_name, v_type, kV])
                if kV in [None,'0']:  self.f_key[v_type][v_name] = kV 
                else: 
                    kp,kname = "/".join(kV.split('/')[0:-1]), kV.split('/')[-1] 
                    if os.path.exists(kp): self.f_key[v_type][v_name] = os.path.abspath(kp)+'/'+kname 
                    else: bridge_error('Invalid path: '+v+' '+kV) 
            elif kV is None and v in config_key:  
                vars(self.args)[v] = config_key[v] 
        return
    

    def verify_prefixes(self, config_key):
        for k,k_pre in self.f_key['prefix'].items(): 
            k_name = k+'_prefix' 
            if k_pre is None or k_pre == '0':  
                if k in self.required: self.errors['missing'].append([k_name,'prefix',k_pre]) 
                continue  
            
            
            
            pK, pF, p_dir, p_start = {}, {}, "/".join(k_pre.split('/')[0:-1]), k_pre.split('/')[-1]  
            p_files = [p_dir+'/'+f_name for f_name in os.listdir(p_dir) if p_start in f_name]
            if len(p_files) == 0: 
                    if k in self.required: self.errors['invalid'].append([k_name,'prefix',k_pre]) 
                    else:                  self.warnings['invalid'].append([k_name,'prefix', k_pre]) 
            if k == 'sumstats':  
                suffix_key = dd(list)  
                for pn in p_files: 
                    p_tail = pn.split(k_pre)[-1] 
                    p_chr = p_tail.split('.')[0]
                    p_suffix = '.'+".".join(p_tail.split('.')[1::]) 
                    if p_suffix.split('.')[-1] not in ['log']: suffix_key[p_suffix].append([p_chr, pn]) 
                suffix_cands = [kpp for kpp in suffix_key.keys()] 
                if self.args.sumstats_suffix is not None: 
                    if self.args.sumstats_suffix not in suffix_cands: parse_error('Invalid Sumstats Suffix: '+self.args.sumstats_suffix)
                    else:                                             my_suffix  = self.args.sumstats_suffix 
                else: 
                    if len(suffix_cands) > 1:                         parse_error('--sumstats_suffix required to resolve ambiguous sumstats extensions: '+",".join(suffix_cands)) 
                    else:                                             my_suffix = suffix_cands[0]    
                
                if self.args.sumstats_suffix is None: 
                    self.args.sumstats_suffix = my_suffix 
                    self.warnings['missing'].append(['sumstats_suffix', 'suffix', 'inferred']) 



                for cr, fp in suffix_key[my_suffix]:  pK[cr] = fp 
                self.map[k] = pK 
                sum_stats_fields, sum_stats_errors = ['ssf-alt', 'ssf-beta', 'ssf-maf', 'ssf-p', 'ssf-ref', 'ssf-se', 'ssf-snpid', 'ssf-ss'], [] 

                for ks in sum_stats_fields: 
                    kArg = vars(self.args)[ks] 
                    if kArg is None: sum_stats_errors.append(ks) 
                    else:            pF[ks.split('-')[-1].upper()] = kArg 
                if len(sum_stats_errors) > 0: parse_error('Sumstats field identifiers are missing: '+",".join(sum_stats_errors))  
                self.fields[k] = pF                 
        return 


    def verify_files(self, config_key): 


        for k,f_path in self.f_key['file'].items(): 
            k_name = k+'_file' 
            if f_path is None or f_path == '0':  
                if k in self.required: self.errors['missing'].append([k_name,'file',f_path]) 
                continue  
            
            if not os.path.isfile(f_path): 
                if k in self.required: self.errors['invalid'].append([k_name,'file',f_path]) 
                else:                  self.warnings['invalid'].append([k_name, 'file', f_path]) 
            pF = {} 
            if k in ['pheno']: 
                pheno_fields, pheno_errors = ['pf-name', 'pf-covariates'] , [] 
                for ks in pheno_fields: 
                    kArg = vars(self.args)[ks] 
                    if kArg is None: pheno_errors.append(ks) 
                    else:            pF[ks.split('-')[-1].upper()] = kArg 
                if len(pheno_errors) > 0: parse_error('Phenotype file field identifiers are missing: '+",".join(pheno_errors)) 
                self.fields[k] = pF 


def LOAD_REQUIREMENTS(module, cmd): 
    if module == 'prs-single': 
        if   cmd == 'clump':    return ['bfile','sumstats','id','snp'] 
        elif cmd == 'eval':    return  ['bfile','sumstats','id','clump'] 
        elif cmd == 'predict':  return ['bfile','eval','pheno','validation']
        elif cmd == 'quantify': return ['bfile','eval','predict','pheno','validation']
        else:                   return ['bfile','sumstats','id','snp','pheno','validation'] 
    
    elif module == 'prs-port': 
        if cmd == 'predict':  return ['bfile','eval','pheno','validation','model']
        elif cmd == 'quantify': return ['bfile','eval','predict','pheno','validation','model']
        else:                   return ['bfile','sumstats','id','snp','pheno','validation','model'] 
    
    elif module == 'prs-prior': 
        if   cmd == 'clump':    return ['bfile','sumstats','id','snp'] 
        elif cmd == 'eval':    return  ['bfile','sumstats','id','clump','model'] 
        elif cmd == 'predict':  return ['bfile','eval','pheno','validation','model']
        elif cmd == 'quantify': return ['bfile','eval','predict','pheno','validation','model']
        else:                   return ['bfile','sumstats','id','snp','pheno','validation','model'] 
    
    elif module == 'build-model': 
        if   cmd == 'clump':     return  ['bfile','sumstats','id','snp','clump-field','clump-snp-field'] 
        elif cmd == 'eval':      return   ['bfile','sumstats','id','clump'] 
        elif cmd == 'optimize':  return ['bfile','eval','pheno']
        elif cmd == 'prior':     return     ['bfile','eval','optimize','pheno']
        else:                   return    ['bfile','sumstats','id','snp','pheno','clump-field','clump-snp-field'] 
    return [] 
    



def LOAD_CONFIGS(configs, bd, X='NA'):
    K, MODULE = dd(bool) , X 
    if len(configs) == 0: return K 
    for config_file in configs: 
        local_path = "/".join(config_file.split('/')[0:-1]) 
        config_handle = open(config_file) 
        for line in config_handle: 
            if len(line) < 2 or line[0] == '#': continue 
            line=line.split('#')[0].strip() 
            if line[0] == '@': MODULE = line.split('@')[-1].lower() 
            elif len(line.split('=')) == 2: 
                kname = line.split('=')[0].lower() 
                kval = line.split('=')[1] 
                k1, k2 = kname.split('_')[0], kname.split('_')[-1] 
                if kname in K:  bridge_error('REPEATED OPTION IN CONFIGURATION FILES: '+kname) 
                if k2 not in ['file','prefix']: 
                    K[kname] = kval 
                else:
                    kp, kn = "/".join(kval.split('/')[0:-1]), kval.split('/')[-1] 
                    if len(kp) > 0: 
                        if os.path.exists(kp): full_name = os.path.abspath(kp)+'/'+kn 
                        elif os.path.exists(bd+'/'+kp): full_name = os.path.abspath(bd+'/'+kp)+'/'+kn 
                        elif os.path.exists(bd+'/'+local_path+'/'+kp): full_name = os.path.abspath(bd+'/'+local_path+'/'+kp)+'/'+kn 
                        elif (os.path.exists(local_path+'/'+kp)): full_name = os.path.abspath(local_path+'/'+kp)+'/'+kn  
                        else:                                     bridge_error('INVALID PATH SUPPLIED IN CONFIG FILE: '+kname+', '+kval) 
                    else:
                        if os.path.exists(bd+'/'+local_path+'/'+kp): full_name = os.path.abspath(bd+'/'+local_path)+'/'+kn 
                        else:                                        full_name = os.path.abspath(local_path)+'/'+kn 
                    
                    if k2 == 'file' and os.path.isfile(full_name): K[kname] = full_name 
                    elif k2 == 'file':                             bridge_error('INVALID FILE PATH SUPPLIED IN CONFIG FILE: '+kname+', '+kval) 
                    else: 
                        pp, pn = "/".join(full_name.split('/')[0:-1]), full_name.split('/')[-1] 
                        if len([x for x in os.listdir(pp) if x[0:len(pn)] == pn]) > 0: K[kname] = full_name 
                        else: bridge_error('INVALID FILE PREFIX SUPPLIED IN CONFIG FILE: '+kname+', '+kval) 
        config_handle.close() 
    return K 
        

