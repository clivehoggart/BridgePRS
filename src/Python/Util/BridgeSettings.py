import sys, os 
from collections import defaultdict as dd
from .BridgeProgress  import BridgeProgress
from .BridgePop      import BridgePop

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


DEFAULT = dd(lambda: None) 
POP_ARGS = ['pop','ldpop','sumstats_prefix','sumstats_suffix','genotype_prefix','phenotype_file']
for a,b in [['fst',0.15]]: DEFAULT[a] = b 
for a,b in [["p","P"],["snpid","ID"],["se","SE"],["n","OBS_CT"],["beta","BETA"],["ref","REF"],["alt","A1"],["maf","A1_FREQ"]]: DEFAULT['ssf-'+a] = b 


class BridgeSettings: 
    def __init__(self,io):
        self.args, self.io, self.module, self.cmd, self.pop = io.args, io, io.args.module, io.args.cmd, None
        self.files, self.prefixes, self.lists = {}, {}, {} 
        

    def update_inputs(self, pipeline_key):
        for v in vars(self.args): 
            kV = vars(self.args)[v] 
            if v in ['config','snp_file','validation_file'] + POP_ARGS or kV in [None, []]: continue 
            if v.split('_')[-1] not in ['file','prefix','files']: continue 
            if v.split('_')[-1] == 'file':     self.files[v.split('_')[0]] = kV 
            elif v.split('_')[-1] == 'prefix': self.prefixes[v.split('_')[0]] = kV 
            elif v.split('_')[-1] == 'files':  self.lists[v.split('_')[0]] = kV 
        
        for k,x in pipeline_key['prefix'].items(): self.prefixes[k] = x 
        for k,x in pipeline_key['file'].items():   self.files[k] = x 

        #if self.pop is not None: self.pop.validate(self.args.module, self.args.cmd, self.files, self.prefixes, self.lists) 
        # bro  
        if self.pop is not None: self.pop.validate(self.files, self.prefixes, self.lists) 
        return
        
    
    def check_analysis_data(self):  
        return  
    
    
    #################   POP CHECK ############################
    def check_pop_data(self):
        self.load_ld_ref()
        
        pop_key = self.resolve_pop_args()
        

        if len(pop_key['pop']) != self.popnum:   bridge_error(['Insufficient Number of Population Names ('+str(len(pop_key['pop']))+') For Subprogram: '+self.args.module+' '+self.args.cmd])         
        
        if self.module.split('-')[0] == 'build': self.pop_data = [BridgePop(self.args, self.io.paths, pop_key['pop'][0], {p: pop_key[p][0] for p in POP_ARGS}, self.ld_ref, 'BASE')] 
        else:                                    self.pop_data = [BridgePop(self.args, self.io.paths, pop_key['pop'][0], {p: pop_key[p][0] for p in POP_ARGS}, self.ld_ref, 'TARGET')] 

        if len(pop_key['pop']) > 1:              
            self.pop_data.append(BridgePop(self.args, self.io.paths, pop_key['pop'][1], {p: pop_key[p][-1] for p in POP_ARGS}, self.ld_ref, 'BASE', self.pop_data[0]))


        self.io.progress.show_pop_data(self.pop_data) 
        
        self.pop = self.pop_data[0] 
        return self 
    
    def resolve_pop_args(self): 
        KL, POP_KEY, ARG_KEY = dd(list), dd(list), {v: vars(self.args)[v] for v in vars(self.args)} 
        if self.args.cmd in ['go','pops']: self.popnum = 2 
        else:                              self.popnum = 1 
        

        for i,K in enumerate(self.args.config):
            for k,v in K.items(): 
                if i == 0 or k in KL:              KL[k].append(v)
                elif k not in ['genotype_prefix','phenotype_file','validation_file']: KL[k].append(v) 
                else:                              continue  
        
        for k,kl in KL.items():
            if k in ARG_KEY and k in POP_ARGS: vars(self.args)[k] = (ARG_KEY[k]+kl)[0:self.popnum] 
            elif k in ARG_KEY: 
                kA,kC = ARG_KEY[k], list(set(kl))[0] 
                if len(list(set(kl))) > 1: bridge_error('Incompatible Arguments In Configuration Files: '+k+': '+",".join(kl))
                elif kA ==  None:          vars(self.args)[k] = kC  
                elif kA == DEFAULT[k]:     vars(self.args)[k] = kC 
                else: 
                    print(kA, kC,'yo') 
                    sys.exit() 
        
        if len(self.args.ldpop) == 0: self.args.ldpop = self.args.pop 
        for lp in self.args.ldpop:
            if lp not in self.ld_ref: bridge_error('No LD-Reference Supplied For: '+lp) 
        for v in vars(self.args): 
            kV = vars(self.args)[v] 
            if v in POP_ARGS: 
                if kV not in [None, []]:        POP_KEY[v] = kV 
                elif v not in ['pop']:          POP_KEY[v] = [None] 
                else:                           continue 
            #elif v.split('_')[-1] == 'file':   self.files[v.split('_')[0]] = kV 
            #elif v.split('_')[-1] == 'prefix': self.prefixes[v.split('_')[0]] = kV 
        return POP_KEY  
        


    def load_ld_ref(self): 
        self.ld_ref = {} 
        if '1000G_ref' in os.listdir(self.io.bridgedir+'/data'): self.load_ld(self.io.bridgedir+'/data/1000G_ref')                                                                                                                                                                       
        elif '1000G_sample' in os.listdir(self.io.bridgedir+'/data'): self.load_ld(self.io.bridgedir+'/data/1000G_sample')                                                                                                                                                               
        for p in self.args.ld_path:     self.load_ld(p)                                                                                                                                                                                                                                          

    def load_ld(self, ld_path):                                                                                                                                                                                                                                                     
        ids, cands, key, k = [], [], {}, 0                                                                                                                                                                                                                                          
        for f in os.listdir(ld_path):
            if f[0] in ['.','_']: continue 
            if f.split('.')[-1] in ['bed','bim','fam']: cands.append(".".join(f.split('.')[0:-1]))                                                                                                                                                                                  
            elif 'IDS' in [z.upper() for z in f.split('.')]:                            ids.append([f.split('.')[0], ld_path+'/'+f])                                                                                                                                                                                
        cands = list(set(cands)) 
        

        if    len(cands) == 0: bridge_error('No plink files (bed, bim, fam) found in LD-Reference Path: '+ld_path) 
        elif  len(cands) == 1: 
            bridge_error('LD-Reference files must be split by chromosome: '+ld_path) 
            for pop, pop_path in ids: 
                #print(pop, cands, pop_path, ld_path+'/'+cands[0], False, {cands[0]: ld_path+'/'+cands[0]}) 
                self.ld_ref[pop.upper()] = [pop_path, ld_path+'/'+cands[0], False, {cands[0]: ld_path+'/'+cands[0]}]                                                                                                                                                                                    
            return  
        else: 
            BYCHR = True 
            while True:                                                                                                                                                                                                                                                                 
                my_prefix = list(set([c[0:k] for c in cands]))[0]                                                                                                                                                                                                                          
                ck1 = list(set([c[0:k+1] for c in cands]))                                                                                                                                                                                                                              
                if len(ck1) == 1: k+=1                                                                                                                                                                                                                                                  
                else: break                                                                                                                                                                                                                                                             
            for cand in cands:                                                                                                                                                                                                                                                          
                if len(my_prefix) > 0: chr_cand = cand.split(my_prefix)[-1]                                                                                                                                                                                                                                       
                else:                  chr_cand = cand 
                key[chr_cand] = ld_path+'/'+cand                                                                                                                                                                                                                                        
            
            for pop, pop_path in ids: self.ld_ref[pop.upper()] = [pop_path, ld_path+'/'+my_prefix, BYCHR, key]                                                                                                                                                                                    
            return  
        return 

        
   

    
    
