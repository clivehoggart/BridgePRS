import sys, os 
from collections import defaultdict as dd
from .BridgeProgress  import BridgeProgress
from .BridgeData      import BridgeData







def bridge_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeSettingsError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                     '+es+'\n')
    else: sys.stderr.write('\nBridgeSettingsError: '+eString+'\n')
    sys.exit(2) 

def parse_error(eString, SAVE=False):
    if not SAVE: 
        sys.stderr.write('\nBridgeParseError: '+eString+'\n')
        sys.exit(2) 
    sys.stderr.write('BridgeParseError: '+eString+'\n')
    return  
            
def parse_warning(eString, eRule): 
    if eRule == 'inferred': sys.stderr.write('BridgeWarning: --'+eString+' not supplied by command line or config file (inferring argument from process of elimination)\n') 
    else:                   sys.stderr.write('BridgeWarning: '+eString+'\n') 
    










class BridgeSettings: 
    def __init__(self,io):
        self.args, self.io, self.module, self.cmd = io.args, io, io.args.module, io.args.cmd
        self.config = self.load_configs(self.args.pop_config, self.io.bridgedir) 
        if len(self.args.pop) == 0: self.args.pop = self.config['name'] 
         

    def load_configs(self, configs, bp):
        K = dd(lambda: 'NA')  
        K['name'] = [] 
        if len(configs) == 0: return K
        for config_file in configs: 
            lp = "/".join(config_file.split('/')[0:-1]) 
            try:     f_handle = open(config_file) 
            except:  bridge_error('Invalid filepath supplied: '+config_file) 
            for ln in f_handle: 
                ln = ln.split('#')[0].strip().split('=') 
                if len(ln) != 2: continue 
                k_opt, k1, k2, k_val, k_data = ln[0].lower(), ln[0].lower().split('_')[0], ln[0].lower().split('_')[-1], ln[1], [] 
                if k_opt == 'name': K[k_opt].append(k_val) 
                elif k2 not in ['file','files','prefix']: K[k_opt] = k_val 
                else:
                    for kf in k_val.split(','): 
                        for kd in kf.split(): 
                            kp, kn = "/".join(kd.split('/')[0:-1]), kd.split('/')[-1] 
                            if os.path.exists(kp):          fn = os.path.abspath(kp)+'/'+kn 
                            elif os.path.exists(lp+'/'+kp): fn = os.path.abspath(lp+'/'+kp)+'/'+kn 
                            elif os.path.exists(bp+'/'+kp): fn = os.path.abspath(bp+'/'+kp)+'/'+kn 
                            else:                           bridge_error('INVALID PATH SUPPLIED IN CONFIG FILE: '+k_opt+': '+kd)  
                            k_data.append(fn) 
                    if k2 == 'files':       K[k_opt] = k_data 
                    elif len(k_data) == 1:  K[k_opt] = k_data[0] 
                    else:                   bridge_error('TWO MANY FILES SUPPLIED (MAX 1): '+k_opt+': '+k_val)  
            f_handle.close() 
        return K 
            



    def check(self):
        self.verify(dd(lambda: dd(bool)), True)  
        ss, bd, pf = self.input.sumstats, self.input.bdata, self.input.phenotypes
        if len(self.args.pop) == 0: bridge_error('A population name is required --pop') 
        pop_name = self.args.pop[0] 
        pf_name = pf.fields['NAME'] 
        config_name = self.io.paths['tmp']+'/'+pop_name+'.'+pf_name+'.config.txt'
        w = open(config_name, 'w') 
        w.write('NAME='+pop_name+'\n')
        w.write('SUMSTATS_PREFIX='+ss.prefix+'\n') 
        w.write('SUMSTATS_SUFFIX='+ss.suffix+'\n') 
        w.write('SNP_FILE='+ss.snp_file+'\n') 
        for n,v in ss.fields.items(): w.write('SSF-'+n.upper()+'='+v+'\n')
        w.write('BFILE_PREFIX='+bd.prefix+'\n') 
        w.write('ID_FILE='+bd.id_file+'\n') 
        pheno_files = pf.names 
        w.write('PHENO_FILES='+pf.names+'\n') 
        for n,v in pf.fields.items(): 
            if n == 'NAME': w.write('PF-NAME='+v+'\n')
            elif n == 'COVARIATES': w.write('PF_COVARIATES='+v+'\n') 
        w.close() 
        self.io.progress.finish_data_validation(config_name) 
        return 


    def verify(self, pipeline_key, CHECK = False):
        self.input = BridgeParser(self, CHECK).load_inputs(self.config, pipeline_key) 
        if len(self.input.warnings) > 0: self.report_warnings(self.input.warnings['invalid'], self.input.warnings['recovered'], self.input.warnings['missing']) 
        self.files, self.prefixes = self.input.f_key['file'], self.input.f_key['prefix']
        self.chromosomes = self.input.sumstats.chromosomes 
        return self 
        
    def report_warnings(self, invalids, recovered, missing): 
        for k,k_type,k_found in invalids:
            if k_type == 'prefix':   sys.stderr.write('BridgeWarning: ---'+k+' '+k_found+' is an invalid file prefix (Prefix not required for this job and will be ignored)\n') 
            else:                    sys.stderr.write('BridgeWarning: ---'+k+' '+k_found+' is an invalid file name   (File not required for this job and will be ignored)\n') 
        for k,k_type,k_found in recovered:  sys.stderr.write('BridgeWarning: --'+k+'_'+k_type+' not supplied by command line or config file (using file generated from previous run: '+self.io.pipeline.progress_file+')\n') 
        for k,k_type,k_found in missing:    sys.stderr.write('BridgeWarning: --'+k+' not supplied by command line or config file (inferring argument from process of elimination)\n') 
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
            elif k == 'beta_file': 
                c_dir = self.io.paths['run']+'/beta'
                hints = ['beta files should be located at: '+c_dir+' (maybe call ./bridgePRS single beta?)'] 
                try:    found = [f for f in os.listdir(c_dir) if self.args.popname in f and 'best_model_params' in f][0] 
                except: found = 'NA' 
                if found != 'NA': hints = ['Try: --beta_file '+c_dir+'/'+found] 
            sys.stderr.write('BridgeIOError: '+k+' is missing\n') 
            for h in hints: sys.stderr.write('   BridgeHint: '+h+'\n') 
        return  
    
        



class BridgeParser:
    def __init__(self, settings, CHECK):
        self.args, self.module, self.cmd, self.paths, self.CHECK = settings.args, settings.module, settings.cmd, settings.io.paths, CHECK
        self.defaults, self.recommended, self.required = self.load_requirements() 
        #self.required = LOAD_REQUIREMENTS(settings.module, settings.cmd) 
        #self.defaults, self.recommended, self.required = LOAD_REQUIREMENTS(settings.module, settings.cmd) 
        self.f_key = dd(lambda: dd(lambda: None))  
        self.errors, self.warnings = dd(list), dd(list) 


    def load_requirements(self): 
        DEFAULTS, REC, REQ = dd(lambda: None), [], [] 
        DEFAULTS['fst'] = 0.15 
        m1, m2 = self.module.split('-')[0], self.module.split('-')[-1] 
        if m1 in ['prs','build','easyrun']: 
            if m1 == 'prs': REQ =   ['bfile','sumstats','id','snp','pheno','validation'] 
            else:           REQ =   ['bfile','sumstats','id','snp','pheno'] 
            if m2 in ['port','prior']: REQ.append('model') 
            if self.cmd in ['predict','quantify']: REQ.append('beta') 
            if self.cmd in ['quantify']:           REQ.append('predict') 
        else: 
            if m1 == 'check' and self.cmd[0:3] != 'req': REQ = ['bfile', 'sumstats','pheno'] 
        return DEFAULTS, REC, REQ 
        
    
    def load_inputs(self, config_key, pipeline_key):
        
        FK = self.parse_triple(config_key, pipeline_key)
        self.sumstats, self.bdata, self.phenotypes = BridgeData(self,'sumstats'), BridgeData(self, 'bdata'), BridgeData(self, 'phenotypes') 


        if   'sumstats' in self.required and FK['prefix']['sumstats'] is not None: self.sumstats.ss_fill(FK['prefix'].pop('sumstats'), FK['file'].pop('snp'))
        elif 'sumstats' in self.required:                                          bridge_error('Sumstats Prefix Required --sumstats_prefix (or in config file)')  
        if   'bfile' in self.required and FK['prefix']['bfile'] is not None:       self.bdata.bf_fill(self.f_key['prefix'].pop('bfile'), self.f_key['file'].pop('id'), CHRS = self.sumstats.chromosomes) 
        elif 'sumstats' in self.required:                                          bridge_error('Bfile Prefix Required --bfile_prefix (or in config file)')  
        if   'pheno' in self.required and FK['files']['pheno'] is not None and len(FK['files']['pheno']) > 0: self.phenotypes.ph_fill(self.f_key['files'].pop('pheno'))
        elif 'sumstats' in self.required:                                                                     bridge_error('Phenotype File(s) Required --pheno_files (or in config file)')  
        return self  
        
        
        
        

    def parse_triple(self, config_key, pipeline_key): 
        for v, v_name, v_type, iVal, cVal in [[v, v.split('_')[0], v.split('_')[-1], vars(self.args)[v],config_key[v]] for v in vars(self.args)]: 
            if v in ['module','cmd','spec','popname','popnames','pop_config','pop_configs','platform','outpath','rpath','plinkpath']: continue 
            if v_type not in ['file','files','prefix']: 
                if cVal != 'NA':
                    if type(iVal) == bool: 
                        if cVal in ['True','False']: vars(self.args)[v] = bool(cVal) 
                        else:                        parse_error('Invalid Option Supplied in Config/Spec File: --'+v+' '+cVal+' (True or False required)')                    
                    elif iVal in [None, [], '0', 0, self.defaults[v]]: vars(self.args)[v] = cVal 
                    continue 
            else: 
                kV, kList = iVal, [] 
                if kV in [None,[]] and cVal != 'NA': kV = cVal 
                if kV in [None,[]] and pipeline_key[v_type][v_name]:     
                    kV = pipeline_key[v_type][v_name] 
                    if v_name in self.required:    self.warnings['recovered'].append([v_name, v_type, kV])
                if kV is None:
                    self.f_key[v_type][v_name] = kV 
                else: 
                    if v_type not in ['files','lists']: kV = [kV] 
                    for kv in kV:         
                        kp,kname = "/".join(kv.split('/')[0:-1]), kv.split('/')[-1] 
                        if os.path.exists(kp): kList.append(os.path.abspath(kp)+'/'+kname) 
                        else:                  bridge_error('Invalid path: '+v+' '+kv) 
                    if v_type not in ['files','lists']: self.f_key[v_type][v_name] = kList[0] 
                    else:                               self.f_key[v_type][v_name] = kList 
        return self.f_key  
                
    
    


        






















def LOAD_REQUIREMENTS99(module, cmd):
    DEFAULTS = dd(lambda: None) 
    DEFAULTS['fst'] = 0.15 
    REC, REQ = [], [] 
    if module == 'prs-single': 
        if   cmd == 'clump':    REQ =   ['bfile','sumstats','id','snp'] 
        elif cmd == 'beta':     REQ =   ['bfile','sumstats','id','clump'] 
        elif cmd == 'predict':  REQ =   ['bfile','beta','pheno','validation']
        elif cmd == 'quantify': REQ =   ['bfile','beta','predict','pheno','validation']
        else:                   REQ =   ['bfile','sumstats','id','snp','pheno','validation'] 
    elif module == 'prs-port': 
        if cmd == 'predict':    REQ =    ['bfile','beta','pheno','validation','model']
        elif cmd == 'quantify': REQ =   ['bfile','beta','predict','pheno','validation','model']
        else:                   REQ =   ['bfile','sumstats','id','snp','pheno','validation','model'] 
    elif module == 'prs-prior': 
        if   cmd == 'clump':    REQ =   ['bfile','sumstats','id','snp'] 
        elif cmd == 'beta':     REQ =  ['bfile','sumstats','id','clump','model'] 
        elif cmd == 'predict':  REQ =   ['bfile','beta','pheno','validation','model']
        elif cmd == 'quantify': REQ =   ['bfile','beta','predict','pheno','validation','model']
        else:                   REQ =   ['bfile','sumstats','id','snp','pheno','validation','model'] 
    
    elif module == 'build-model': 
        if   cmd == 'clump':     REQ =    ['bfile','sumstats','id','snp','clump-field','clump-snp-field'] 
        elif cmd == 'beta':      REQ =    ['bfile','sumstats','id','clump'] 
        elif cmd == 'optimize':  REQ =    ['bfile','beta','pheno']
        elif cmd == 'prior':     REQ =    ['bfile','beta','optimize','pheno']
        else:                    REQ =   ['bfile','sumstats','id','snp','pheno','clump-field','clump-snp-field'] 
    return DEFAULTS, REC, REQ 
    
    











def LOAD_CONFIGS99(configs, bp, X='NA'):
    K = dd(lambda: 'NA')  
    K['name'] = [] 
    if len(configs) == 0: return K
    for config_file in configs: 
        lp = "/".join(config_file.split('/')[0:-1]) 
        try:     f_handle = open(config_file) 
        except:  bridge_error('Invalid filepath supplied: '+config_file) 
        for ln in f_handle: 
            ln = ln.split('#')[0].strip().split('=') 
            if len(ln) != 2: continue 
            k_opt, k1, k2, k_val, k_data = ln[0].lower(), ln[0].lower().split('_')[0], ln[0].lower().split('_')[-1], ln[1], [] 
            if k_opt == 'name': K[k_opt].append(k_val) 
            elif k2 not in ['file','files','prefix']: K[k_opt] = k_val 
            else:
                for kf in k_val.split(','): 
                    for kd in kf.split(): 
                        kp, kn = "/".join(kd.split('/')[0:-1]), kd.split('/')[-1] 
                        if os.path.exists(kp):          fn = os.path.abspath(kp)+'/'+kn 
                        elif os.path.exists(lp+'/'+kp): fn = os.path.abspath(lp+'/'+kp)+'/'+kn 
                        elif os.path.exists(bp+'/'+kp): fn = os.path.abspath(bp+'/'+kp)+'/'+kn 
                        else:                           bridge_error('INVALID PATH SUPPLIED IN CONFIG FILE: '+k_opt+': '+kd)  
                        k_data.append(fn) 
                if k2 == 'files':       K[k_opt] = k_data 
                elif len(k_data) == 1:  K[k_opt] = k_data[0] 
                else:                   bridge_error('TWO MANY FILES SUPPLIED (MAX 1): '+k_opt+': '+k_val)  
        f_handle.close() 
    return K 
        

