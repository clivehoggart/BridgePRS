import sys, os 
from collections import defaultdict as dd
from .BridgeProgress  import BridgeProgress
from .BridgePop      import BridgePop
from .BridgeConfig      import BridgeConfig

def bridge_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeSettingsError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                     '+es+'\n')
    else: sys.stderr.write('\nBridgeSettingsError: '+eString+'\n')
    sys.exit(2) 


def validate_requirements(module, cmd): 
    REQ, m1, m2 = [], module.split('-')[0], module.split('-')[-1] 
    if m1 == 'check' and cmd[0:3] != 'req': return ['ldpop','sumstats_prefix','genotype_prefix','phenotype_files'] 
    if m1 in ['prs','build','easyrun']: 
        if module == 'easyrun':             REQ = ['ldpop','sumstats_prefix','genotype_prefix','phenotype_files'] 
        elif cmd in ['clump','beta']:       REQ = ['ldpop','sumstats_prefix']
        elif cmd in ['pred','quantify']:    REQ = ['genotype_prefix','phenotype_files','beta_prefix'] 
        if m2 in ['port','priot']:          REQ.append('model_file') 
        if  cmd  in ['quantify']:           REQ.append('predict_prefix') 
        
    return REQ 



class BridgeSettings: 
    def __init__(self,io):
        self.args, self.io, self.module, self.cmd = io.args, io, io.args.module, io.args.cmd
        self.command_args = ['pop_config','phenotype','verbose','plinkpath','rpath','platform','cmd','module','outpath','total_cores','cores','restart','silent','noPlots'] 
        self.config = BridgeConfig(self).load_file(self.args.pop_config, self.io.bridgedir) 
        self.lists, self.files, self.prefixes = {}, {}, {} 

    
    def update_inputs(self, pipeline_key): 
        for v in vars(self.args): 
            kV = vars(self.args)[v] 
            if v in self.config.pop_specific or kV in [None, []]: continue 
            if v.split('_')[-1] not in ['file','prefix','files']: continue 
            if v.split('_')[-1] == 'file':     self.files[v.split('_')[0]] = kV 
            elif v.split('_')[-1] == 'prefix': self.prefixes[v.split('_')[0]] = kV 
            elif v.split('_')[-1] == 'files':  self.lists[v.split('_')[0]] = kV 

        for k,x in pipeline_key['prefix'].items(): self.prefixes[k] = x 
        for k,x in pipeline_key['file'].items():   self.files[k] = x 
        return
        
    
    
    
    
    
    def check_analysis_data(self):  
        return  
    
    
    
    
    
    
    
    #################   POP CHECK ############################
    
    def check_pop_data(self):
        self.required = validate_requirements(self.module, self.cmd) 
        self.disambiguate_pop_data()         
        self.config.parse_arg_data()        #### FIX THIS  
        self.io.progress.show_pop_data(self.pop_data) 
        return self 

       
    def disambiguate_pop_data(self): 
        #self.ld_ref, self.pop_lists, d_strs  =  {}, {}, {} 
        
        
        self.pop_lists, self.pop_strs = self.config.get_pop_lists() 
        self.load_ld_ref() 
        #if len(self.args.pop) == 0: bridge_error('A target population name is required --pop or using a config file [--pop_config]') 
        #if len(self.args.pop) >  2: bridge_error('At most two population name(s) are allowed --pop '+",".join(self.args.pop)) 
        
        if len(self.pop_lists['genotype_prefix']) != len(self.pop_lists['phenotype_files']): bridge_error(['Corresponding Phenotype And Genotype Data Required',d_strs['genotype_prefix'],d_strs['phenotype_files']])
        for p,D in self.pop_lists.items(): 
            if   len(D) == 0 and p not in self.required:  self.pop_lists[p] = [None]  
            elif len(D) == 0 and p in self.required:      bridge_error('Missing Required Population Data: '+p)                      
            elif len(D) > 1 and len(self.args.pop) == 1:  bridge_error(['Too Much Population Data For: '+p,'At Most 1 File is Supported For 1 Population',self.pop_strs[p]])
            elif len(D) > 2 and len(self.args.pop) == 2:  bridge_error(['Too Much Population Data For: '+p,'At Most 2 File(s) Supported For 2 Populations',self.pop_strs[p]]) 
        
        self.pop_data = [BridgePop(self.args, self.io.paths, self.args.pop[0], {p: D[0] for p,D in self.pop_lists.items()}, self.ld_ref)] 
        if len(self.args.pop) > 1:  self.pop_data.append(BridgePop(self.args, self.io.paths, self.args.pop[1], {p: D[-1] for p,D in self.pop_lists.items()},self.ld_ref)) 
        self.pop = self.pop_data[0] 
       
    def load_ld_ref(self): 
        self.ld_ref = {} 
        if '1000G_ref' in os.listdir(self.io.bridgedir+'/data'): self.load_ld(self.io.bridgedir+'/data/1000G_ref')                                                                                                                                                                       
        elif '1000G_sample' in os.listdir(self.io.bridgedir+'/data'): self.load_ld(self.io.bridgedir+'/data/1000G_sample')                                                                                                                                                               
        for p in self.args.ldpath:                                                                                                                                                                                                                                                 
            if os.path.exists(p):                      mp = os.path.abspath(p)                                                                                                                                                                                                                         
            elif os.path.exists(self.bridgedir+'/'+p): mp = os.path.abspath(self.bridgedir+'/'+p)                                                                                                                                                                                  
            else:                                      bridge_error(['--ldPath '+p+' Does not exist'])                                                                                                                                                                             
            self.load_ld(mp)

    def load_ld(self, ld_path):                                                                                                                                                                                                                                                     
        ids, cands, key, k = [], [], {}, 0                                                                                                                                                                                                                                          
        for f in os.listdir(ld_path):                                                                                                                                                                                                                                               
            if f.split('.')[-1] in ['bed','bim','fam']: cands.append(".".join(f.split('.')[0:-1]))                                                                                                                                                                                  
            elif 'ids' in f:                            ids.append([f.split('.')[0], ld_path+'/'+f])                                                                                                                                                                                
        cands = list(set(cands))                                                                                                                                                                                                                                                    
        while True:                                                                                                                                                                                                                                                                 
            my_prefix = list(set([c[0:k] for c in cands]))[0]                                                                                                                                                                                                                          
            ck1 = list(set([c[0:k+1] for c in cands]))                                                                                                                                                                                                                              
            if len(ck1) == 1: k+=1                                                                                                                                                                                                                                                  
            else: break                                                                                                                                                                                                                                                             
        for cand in cands:                                                                                                                                                                                                                                                          
            chr_cand = cand.split(my_prefix)[-1]                                                                                                                                                                                                                                       
            key[chr_cand] = ld_path+'/'+cand                                                                                                                                                                                                                                        
        if len(key) > 1: BYCHR = True 
        else:            BYCHR = False 
        for pop, pop_path in ids: self.ld_ref[pop.upper()] = [pop_path, ld_path+'/'+my_prefix, BYCHR, key]                                                                                                                                                                                    
        return                                                                                                                                                                                                                                                                      
    
