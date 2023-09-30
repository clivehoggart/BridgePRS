import sys, os 
from collections import defaultdict as dd
from .BridgeProgress  import BridgeProgress
from .BridgePop      import BridgePop

def bridge_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeConfigError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                     '+es+'\n')
    else: sys.stderr.write('\nBridgeConfigError: '+eString+'\n')
    sys.exit(2) 

def validate_paths(v, kV): 
    valid_paths, vType = [], v.split('_')[-1] 
    if len(kV) == 0: return kV
    for kv in kV: 
        if type(kv) == list: valid_paths.append(validate_paths(v,kv)) 
        else:  
            kpath, kname = "/".join(kv.split('/')[0:-1]), kv.split('/')[-1] 
            if os.path.exists(kpath):  mp = os.path.abspath(kpath) 
            else:                      bridge_error('Invalid path supplied on command line: --'+v+': '+kv)       
            if vType != 'prefix' and kname not in os.listdir(mp): bridge_error('Invalid Filename supplied on command line: --'+v+': '+kv) 
            elif vType == 'prefix' and len([x for x in os.listdir(mp) if x[0:len(kname)] == kname]) == 0:  bridge_error('Invalid prefix supplied on command line: --'+v+': '+kv)       
            else: valid_paths.append(mp+'/'+kname)
    return valid_paths


def validate_requirements(module, cmd): 
    REQ, m1, m2 = [], module.split('-')[0], module.split('-')[-1] 
    if m1 == 'check' and cmd[0:3] != 'req': return ['ldpop','sumstats_prefix','genotype_prefix','phenotype_files'] 
    if m1 in ['prs','build','easyrun']: 
        if m1 == 'easyrun':                    REQ = ['ldpop','sumstats_prefix','genotype_prefix','phenotype_files'] 
        elif cmd in ['clump','beta']:       REQ = ['ldpop','sumstats_prefix']
        elif cmd in ['pred','quantify']:    REQ = ['genotype_prefix','phenotype_files','beta_prefix'] 
        if m2 in ['port','priot']:          REQ.append('model_file') 
        if  cmd  in ['quantify']:           REQ.append('predict_prefix') 
        
        #REQ = ['ldpop','sumstats_prefix','genotype_prefix','phenotype_files'] 
        #if m2 in ['port','prior']: REQ.append('model_file') 
        #if cmd in ['predict','quantify']: REQ.append('beta_prefix') 
        #if cmd in ['quantify']:           REQ.append('predict_prefix') 
    return REQ 


class BridgeConfig: 
    def __init__(self,settings):
        self.args, self.settings = settings.args, settings 
        self.load_defaults() 
        
        
    def load_defaults(self): 
        self.pop_specific = ['pop','ldpop','ldpath','snp_file',  'sumstats_prefix', 'sumstats_suffix', 'genotype_prefix', 'phenotype_files']
        self.command_args = ['pop_config','phenotype','verbose','plinkpath','rpath','platform','cmd','module','outpath','total_cores','cores','restart','silent','noPlots'] 
        self.DEFAULTS = {} 
        for a,b in [['fst',0.15]]: self.DEFAULTS[a] = b 
        for a,b in [["p","P"],["snpid","ID"],["se","SE"],["ss","OBS_CT"],["beta","BETA"],["ref","REF"],["alt","A1"],["maf","A1_FREQ"]]: self.DEFAULTS['ssf-'+a] = b 
        #max_clump_size 0




    def load_file(self, configs, bp):
        
        self.K = dd(list)
        if len(configs) == 0: return self
        for c_idx, config_file in enumerate(configs): 
            lp = "/".join(config_file.split('/')[0:-1]) 
            try:     f_handle = open(config_file) 
            except:  bridge_error('Invalid filepath supplied: '+config_file) 
            for ln in f_handle: 
                ln = ln.split('#')[0].strip().split('=') 
                if len(ln) != 2: continue 
                k_opt, k1, k2, k_val, k_data = ln[0].lower(), ln[0].lower().split('_')[0], ln[0].lower().split('_')[-1], ln[1], [] 
                if k_opt not in self.args:       
                    if k_opt in self.DEFAULTS: continue 
                    else:       bridge_error('Invalid Argument in Configuration File: '+k_opt) 
                elif k_opt in self.command_args: bridge_error('Invalid Argument in Configuration File: '+k_opt+' (system arguments can only be passed on command line --'+k_opt+')')
                elif k2 not in ['file','files','prefix']: self.K[k_opt].append(k_val)
                else:
                    for kf in k_val.split(','): 
                        for kd in kf.split(): 
                            kp, kn = "/".join(kd.split('/')[0:-1]), kd.split('/')[-1] 
                            if os.path.exists(kp):          fn = os.path.abspath(kp)+'/'+kn 
                            elif os.path.exists(lp+'/'+kp): fn = os.path.abspath(lp+'/'+kp)+'/'+kn 
                            elif os.path.exists(bp+'/'+kp): fn = os.path.abspath(bp+'/'+kp)+'/'+kn 
                            else:                           bridge_error('INVALID PATH SUPPLIED IN CONFIG FILE: '+k_opt+': '+kd)  
                            k_data.append(fn) 
                    if c_idx > 0 and len(self.K[k_opt]) == 0: self.K[k_opt] = ['NA'] 
                    if k2 == 'files':       self.K[k_opt].append(k_data)  
                    elif len(k_data) == 1:  self.K[k_opt].append(k_data[0]) 
                    else:                   bridge_error('TWO MANY FILES SUPPLIED (MAX 1): '+k_opt+': '+k_val)  
            f_handle.close()  
        return self 
            
   



        
    def get_pop_lists(self): 
        POP_LIST, POP_STRS = {} , {}  
        #self.ld_ref, self.pop_lists, d_strs  =  {}, {}, {} 
        #self.pop_specific = ['snp_file', 'ldpop', 'sumstats_prefix', 'sumstats_suffix', 'genotype_prefix', 'phenotype_files']

        for v in ['pop','ldpop','ldpath','sumstats_suffix']: 
            kV  = vars(self.args)[v] 
            if v ==    'pop':       vars(self.args)[v] = list(kV + self.K[v])[0:2] 
            elif v == 'ldpop':      vars(self.args)[v] = list(list(kV + self.K[v][0:2]) + self.args.pop)[0:2] 
            else:                   vars(self.args)[v] = list(kV + self.K[v]) 
        for v in ['snp_file','sumstats_prefix','genotype_prefix','phenotype_files']: vars(self.args)[v] = list(validate_paths(v,kV) + self.K[v])
        
        
        
        #if len(self.args.pop) == 0: bridge_error('A target population name is required --pop or using a config file [--pop_config]') 
        #if len(self.args.pop) >  2: bridge_error('At most two population name(s) are allowed --pop '+",".join(self.args.pop)) 
        #self.load_ld_ref() 
        for p in self.pop_specific: 
            POP_LIST[p] = vars(self.args)[p] 
            POP_STRS[p] = ", ".join([','.join(d) if type(d) == list else d for d in POP_LIST[p]])
        
        
        return POP_LIST, POP_STRS
        
        #if len(self.pop_lists['genotype_prefix']) != len(self.pop_lists['phenotype_files']): bridge_error(['Corresponding Phenotype And Genotype Data Required',d_strs['genotype_prefix'],d_strs['phenotype_files']])
        #for p,D in self.pop_lists.items(): 
        #    if   len(D) == 0 and p not in self.required:  self.pop_lists[p] = [None]  
        #    elif len(D) == 0 and p in self.required:      bridge_error('Missing Required Population Data: '+p)                      
        #    elif len(D) > 1 and len(self.args.pop) == 1:  bridge_error(['Too Much Population Data For: '+p,'At Most 1 File is Supported For 1 Population',d_strs[p]])
        #    elif len(D) > 2 and len(self.args.pop) == 2:  bridge_error(['Too Much Population Data For: '+p,'At Most 2 File(s) Supported For 2 Populations',d_strs[p]]) 
        #self.pop_data = [BridgePop(self.args, self.io.paths, self.args.pop[0], {p: D[0] for p,D in self.pop_lists.items()}, self.ld_ref)] 
        #if len(self.args.pop) > 1:  self.pop_data.append(BridgePop(self.args, self.io.paths, self.args.pop[1], {p: D[-1] for p,D in self.pop_lists.items()},self.ld_ref)) 
        #self.pop = self.pop_data[0] 
        
    


    def parse_arg_data(self): 
        for v in vars(self.args): 
            kV  = vars(self.args)[v] 
            if len(self.K[v]) == 0 or v in self.pop_specific + self.command_args: continue  
            if len(list(set(self.K[v]))) > 1: bridge_error('Incompatible Arugments In Configuration Files: '+v+': '+",".join(self.K[v]))
            cV = list(set(self.K[v]))[0] 
            if v.split('-')[0] == 'ssf' and kV == self.DEFAULTS[v]: vars(self.args)[v] = cV 
            elif v in ['fst']: 
                try: cF = float(cV)
                except: continue 
                if kV == self.DEFAULTS[v] and cF > 0 and cF < 1: vars(self.args)[v] = cF 
            continue               
        
        return 
        

    

