import sys, os 
from collections import defaultdict as dd
from collections import Counter as cc 
from .BridgeProgress  import BridgeProgress
from .BridgePop      import BridgePop

# yo

def bridge_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeSettingsError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                     '+es+'\n')
    else: sys.stderr.write('\nBridgeSettingsError: '+eString+'\n')
    sys.exit(2) 

def bridge_ld_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeLDPanelError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                     '+es+'\n')
    else: sys.stderr.write('\nBridgeLDPanelError: '+eString+'\n')
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
#POP_ARGS = ['pop','ldpop','sumstats_file','sumstats_prefix','sumstats_suffix','genotype_prefix','phenotype_file']
POP_ARGS = ['pop','ldpop','ld_path','ldref','sumstats_prefix','sumstats_suffix','genotype_prefix','phenotype_file']
for a,b in [['fst',0.15]]: DEFAULT[a] = b 
for a,b in [["p","P"],["snpid","ID"],["se","SE"],["n","OBS_CT"],["beta","BETA"],["ref","REF"],["alt","A1"],["maf","A1_FREQ"]]: DEFAULT['ssf-'+a] = b 


class BridgeSettings: 
    def __init__(self,io):
        self.args, self.io, self.module, self.cmd, self.pop = io.args, io, io.args.module, io.args.cmd, None
        self.files, self.prefixes, self.lists, self.paths = {}, {}, {}, {} 
        

    def update_inputs(self, pipeline_key):
        for v in vars(self.args): 
            kV = vars(self.args)[v] 
            if v in ['config','snp_file','validation_file'] + POP_ARGS or kV in [None, []]: continue 
            if v.split('_')[-1] not in ['file','prefix','files']: continue 
            if v.split('_')[-1] == 'file':     self.files[v.split('_')[0]] = kV 
            elif v.split('_')[-1] == 'prefix': self.prefixes[v.split('_')[0]] = kV 
            elif v.split('_')[-1] == 'files':  self.lists[v.split('_')[0]] = kV 
            elif v.split('_')[-1] == 'path':  self.paths[v.split('_')[0]] = kV 
       

        for k,x in pipeline_key['prefix'].items(): self.prefixes[k] = x 
        for k,x in pipeline_key['file'].items():   self.files[k] = x 
        if self.pop is not None: self.pop.validate(self.files, self.prefixes, self.lists) 
        return
        
    
    def check_analysis_data(self):  
        return  
    
    
    #################   POP CHECK ############################
    def check_pop_data(self):


        pop_key = self.resolve_pop_args()
        
        
        
        if self.module.split('-')[0] == 'build': self.pop_data = [BridgePop(self.args, self.io.progress, self.io.paths, pop_key['pop'][0], {p: pop_key[p][0] for p in POP_ARGS}, 'BASE')] 
        else:                                    self.pop_data = [BridgePop(self.args, self.io.progress, self.io.paths, pop_key['pop'][0], {p: pop_key[p][0] for p in POP_ARGS},  'TARGET')] 
        if len(pop_key['pop']) > 1:    
            self.pop_data.append(BridgePop(self.args, self.io.progress, self.io.paths, pop_key['pop'][1], {p: pop_key[p][-1] for p in POP_ARGS}, 'BASE', self.pop_data[0]))

        self.io.progress.show_pop_data(self.pop_data) 
        self.pop = self.pop_data[0] 
        return self 
    
    def resolve_pop_args(self): 
        KL, POP_KEY, ARG_KEY = dd(list), dd(list), {v: vars(self.args)[v] for v in vars(self.args)} 
        if self.args.cmd in ['go','pops']: self.popnum = 2 
        else:                              self.popnum = 1 
        for i,K in enumerate(self.args.config):
            k_keys = K.keys()  
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
                elif kA == []:             bridge_error('Unrecognized Argument In Configuration Files: '+k+': '+kC) 
                else:                      bridge_error('Universal Argument '+k+' Cannot Be Included In Configuration File (must be passed on command line --'+k+')')
        
        if len(self.args.pop) != len(self.args.ldpop): 
            if len(self.args.ldpop) == 0: self.args.ldpop = [pn for pn in self.args.pop]
            else:                         bridge_error('Ambiguous pop/ldpop assignments: '+",".join(self.args.pop)+" and "+",".join(self.args.ldpop)) 
        
        if len(self.args.pop) != len(self.args.ld_path): 
            if len(self.args.ld_path) == 1: self.args.ld_path.append(self.args.ld_path[0]) 
            else:                           bridge_error('LD Reference Path is Required')         

            

        for v in vars(self.args): 
            kV = vars(self.args)[v] 
            if v in POP_ARGS: 
                if kV not in [None, []]:        POP_KEY[v] = kV 
                elif v not in ['pop']:          POP_KEY[v] = [None] 
                else:                           continue 
        
        if len(POP_KEY['pop']) != self.popnum:   bridge_error(['Insufficient Number of Population Names ('+str(len(POP_KEY['pop']))+') For Subprogram: '+self.args.module+' '+self.args.cmd])         
        elif len(POP_KEY['pop']) > 1 and POP_KEY['pop'][0] == POP_KEY['pop'][-1]: bridge_error(['Insufficient Number of Unique Population Names ('+','.join(POP_KEY['pop'])+') For Subprogram: '+self.args.module+' '+self.args.cmd])         
        POP_KEY['ldref'] = [self.load_ld_panel(ldpop,ld_path) for ldpop,ld_path in zip(POP_KEY['ldpop'],POP_KEY['ld_path'])] 
        return POP_KEY  
        
    
    def load_ld_panel(self, ldpop, ld_path):
        


        ref_panel = ld_path.split('/')[-1] 

        key, cands, id_files = {}, [], [] 
        for f in os.listdir(ld_path): 
            fp = f.split('.')[0].upper().split('_') 
            if f[0] in ['.','_']: continue 
            if f.split('.')[-1] in ['bed','bim','fam']: cands.append(".".join(f.split('.')[0:-1]))                                                                                                                                                                                  
            elif 'IDS' in fp and ldpop.upper() in fp: id_files.append(f)
            else: continue 

        cands = [cn for cn in cc(cands) if cc(cands)[cn] == 3] 
        if    len(cands) == 0: bridge_ld_error('No plink triples (bed, bim, fam) found in LD-Reference Path: '+ld_path) 
        elif  len(cands) == 1: bridge_ld_error('LD-Reference files must be split by chromosome: '+ld_path) 
        elif len(id_files) == 0: bridge_ld_error('No Corresponding ID file found for population '+ldpop+' in: '+ld_path) 
        elif len(id_files) > 1: bridge_ld_error('Multiple ID files found for pop '+ldpop+' in: '+ld_path+'  ('+",".join(id_files)+')') 
        BYCHR = True 
        my_prefix, my_suffix = self.get_prefix_suffix(cands)
        for cand in cands:        
            chr_cand = cand.split(my_prefix)[-1].split(my_suffix)[0] 
            key[chr_cand] = ld_path+'/'+cand                                                                                                                                                                                                                                        


        return [ldpop, ld_path+'/'+id_files[0], ld_path+'/'+my_prefix, key] 


    

    def get_prefix_suffix(self, cands): 
        my_prefix,my_suffix, k = ' ', ' ', 0 
        while True:                                        
            my_prefixes = list(set([c[0:k] for c in cands]))                                                                                                                                                                                                               
            if len(my_prefixes) == 1: 
                my_prefix = my_prefixes[0] 
                k+=1 
            else: 
                break
        k = 0     
        while True:                                        
            my_suffixes = list(set([c[len(c)-k::] for c in cands]))                                                                                                                                                                                                               
            if len(my_suffixes) == 1: 
                my_suffix = my_suffixes[0] 
                k+=1 
            else: 
                break

        if len(my_prefix) == 0: my_prefix = ' '
        if len(my_suffix) == 0: my_suffix = ' '
        return my_prefix, my_suffix  


   

    
    
