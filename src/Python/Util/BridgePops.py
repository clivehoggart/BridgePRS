import sys, os, gzip, shutil
from collections import defaultdict as dd
from collections import Counter as cc 
from .Bridge_Pop.BData import BData 
from .Bridge_Pop.GenoPheno import GenoPheno
from .Bridge_Pop.GWAS import SumStats 



def bridge_debug_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('BridgeDebugError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('BridgeDebugError: '+eString+'\n')
    sys.exit(2) 


# Ambiguous debug_level Ref/Alt Missing snps.txt check printout

def bridge_pop_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgePopError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('\nBridgePopError: '+eString+'\n')
    sys.exit(2) 


class BridgePops: 
    def __init__(self,mps):
        self.paths = {} 
        self.mps   = mps 
    


    def system_tests(self,args): 
        if args.cmd == None: self.mps.error('Module: '+args.module+' needs a command')  
        if args.clean and args.restart: self.mps.error('ParseError: Collision. Only one of --clean (delete all) and --restart (restart jobs) are allowed\n') 
        if args.platform not in ['mac','linux']:  self.mps.error('Unrecognized platform '+args.platform+' (--platform linux/mac are supported)')                                                                                      
        if args.cores == 0: args.cores = args.total_cores - 1
        elif (args.cores < 1) or (args.cores > args.total_cores): args.cores = 1 
        return 

    
    def setup_paths(self, args): 
        if args.outpath is None: 
            if args.cmd == 'check-requirements': args.outpath = 'out' 
            else:                                self.mps.error('IOError: -o, [--out, an output path is required]')                                                                                                                                    
        args.outpath = os.path.abspath(args.outpath) 
        if os.path.isdir(args.outpath) and args.clean: shutil.rmtree(args.outpath)      
        if not os.path.isdir(args.outpath): os.makedirs(args.outpath)                                                                                                                                                           
        if os.path.isdir(args.outpath+'/tmp'): shutil.rmtree(args.outpath+'/tmp') 
        for d in ['logs','tmp','save']: 
            if not os.path.isdir(args.outpath+'/'+d): os.makedirs(args.outpath+'/'+d)    
            self.paths[d] = args.outpath+'/'+d 
        return 
    




    def parse(self): 
        args = self.mps.parse_args()
        argnames = [v for v in vars(args)] 
        self.system_tests(args) 
        self.setup_paths(args) 
        
        self.target, self.base, self.names = None, None, [None, None] 

        if 'config' in argnames and len(vars(args)['config']) > 0: self.parse_pop_configs(args) 
        return args, self 

        if args.module == 'pipeline' or args.module.split('-')[0] in ['prs','build']: self.parse_pop_configs(args)         


                #if args.module == 'pipeline':                      self.parse_two_pops(args) 
                #else:                                              self.parse_one_pop(args) 
        elif args.module == 'analyze': return args, self 


        else: 
            
            print(args.module) 

            print('hmmm') 
            sys.exit() 
        #self.arglen = str(len(vars(args)))
        return args, self 


    def parse_pop_configs(self, args): 
        if args.module == 'build-model': self.base   = BridgePop(self.mps, args.config[0], args, self.paths, 'base') 
        else:                            self.target = BridgePop(self.mps, args.config[0], args, self.paths, 'target') 
        if len(args.config) == 2:        self.base = BridgePop(self.mps, args.config[1], args, self.paths, 'base', self.target) 
        self.names = [x.name for x in [self.target, self.base] if x is not None] 
        return 
        
        







    

    
    def restart_directories(self,args, t_name, b_name): 
        for f in os.listdir(args.outpath): 
            if os.path.isdir(args.outpath+'/'+f): 
                names, fp, fn = f.split('_')[-1].split('-'), f.split('-')[0], f.split('_')[0] 
                if args.cmd == 'go':  
                    if fp == 'prs': 
                        if len(names) == 1 and t_name == names[0]:         shutil.rmtree(args.outpath+'/'+f) 
                        elif len(names) == 2 and f.split('_')[-1] == t_name+'-'+b_name: shutil.rmtree(args.outpath+'/'+f) 
                        else: continue 
                    elif b_name == names[0]:            shutil.rmtree(args.outpath+'/'+f) 
                    else: continue 
        return
                        

    def create_directory_structure2(self,args): 
        f_prev = dd(list) 
        for f in [f for f in os.listdir(args.outpath) if os.path.isdir(args.outpath+'/'+f)]: 
            ns, fp, fn = f.split('_')[-1].split('-'), f.split('-')[0], f.split('_')[0] 
            if args.restart: 
                if args.cmd == 'go': 
                    if (fn == 'prs-single' and ns[0] == self.names[0]) or (fn=='build-model' and ns[0] == self.names[1]): shutil.rmtree(args.outpath+'/'+f) 
                    elif fp == 'prs' and ns[0] == self.names[0] and ns[1] == self.names[1]:                               shutil.rmtree(args.outpath+'/'+f)    
                    else: continue 
                else:
                    print('oh shit') 
                    sys.exit() 
            elif f == 'save': 
                for x in os.listdir(args.outpath+'/'+f): 
                    if x.split('.')[-1] in ['config','source']: 
                        xn = x.split('.')[0]+'-'+x.split('.')[1]  
                        f_prev[xn].append(args.outpath+'/'+f+'/'+x) 
        return f_prev









    def __str__(self):                                                                                                                                                                                                      
        rpr = '<bridgeParseObject> Args='+self.arglen+'; PopData: '
        if self.target is not None: rpr += 'Target='+self.target.name+' '+self.target.valid_str+'; ' 
        if self.base is not None: rpr += 'Base='+self.base.name+' '+self.base.valid_str+';' 
        return rpr 

                                                                                                                                                                                                                            
    def __repr2__(self):                                                                                                                                                                                                     
        rpr = '<bridgeParseObject> Args='+self.arglen+'; PopData: '
        if self.target is not None: rpr += 'Target='+self.target.name+' '+self.target.valid_str+'; ' 
        if self.base is not None: rpr += 'Base='+self.base.name+' '+self.base.valid_str+';' 
        return rpr 






class BridgePop:                                                                                                                                                                                                            
    def __init__(self, mps, P, args, paths, pop_type, prevPop = None): 
        self.args, self.paths, self.type = args, paths, pop_type
        self.sumstats_fields = [pf.strip() for pf in P['sumstats_fields'].split(',')] 
        self.field_key       = {a: A for a,A in zip(['ID','REF','ALT','P','BETA'],self.sumstats_fields)}
        for k,v in P.items(): 
            if k == 'pop':           self.name = v 
            elif k == 'ldpop':       self.ref_pop = v 
            elif k == 'config_name': self.file = v 
            elif v is not None or prevPop is None:       vars(self)[k] = v 
            elif k.split('_')[0] != 'sumstats': vars(self)[k] = vars(prevPop)[k] 
            else:                             vars(self)[k] = v 
    


        self.bdata = BData(self, mps).load_panel()
        

        self.genopheno = GenoPheno(self, mps, args.phenotype).load()    


        self.sumstats  = SumStats(self, mps).load(self.genopheno, prevPop) 
        self.verify_chromosomes() 
        





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
        

        if self.args.module in ['prs-port','prs-prior']: 
            if 'model' not in F: mInputs.append('Base Model --model_file (Hint: Run build-model)') 
            else: 
                model_pairs = [line.strip().split('=') for line in open(F['model'],'rt') if line.strip().split('=')[0] == 'SUMSTATS_SIZE']      
                try:    self.sumstats.model_size = int(model_pairs[0][-1]) 
                except: self.sumstats.model_size = 10000 
        if len(mInputs) + len(mGen) == 0: return True 
        if len(mInputs) > 0: bridge_pop_error(['Missing Input Data:']+mInputs) 
        if len(mGen) > 0:    bridge_pop_error(['Missing Run Data:']+mGen) 
        return True  
        

    def get_chr_strs(self, itbl): 
        c_int, c_str = [], [] 
        for c in itbl: 
            try:               c_int.append(int(c))  
            except ValueError: c_str.append(c) 
        return [str(c) for c in sorted(c_int) + sorted(c_str)]



    def verify_chromosomes(self): 
        chromosomes       = list(set([k for k in self.sumstats.map.keys()]+[k for k in self.bdata.map.keys()]))
        ld_chrs = self.get_chr_strs(self.bdata.map.keys()) 
        ss_chrs = self.get_chr_strs(self.sumstats.map.keys())
        valid_chromosomes = [k for k in self.sumstats.map.keys() if k in self.bdata.map.keys()] 
        self.chromosomes       = self.get_chr_strs(list(set([k for k in self.sumstats.map.keys()]+[k for k in self.bdata.map.keys()])))
        if len(self.chromosomes) == len(ld_chrs) and len(self.chromosomes) == len(ss_chrs) and len(self.chromosomes) == len(valid_chromosomes): return
        my_error = ['Missing Chromosomes'] 
        my_error.append('The Following Chromosomes Are Found in the ld-panel: '+",".join(ld_chrs)) 
        my_error.append('The Following Chromosomes Are Found in the sumstats: '+",".join(ss_chrs))
        if len(valid_chromosomes) == 0: my_error.append('The Following Chromosomes Are Found in both sources: None') 
        else:                           my_error.append('The Following Chromosomes Are Found in both sources: '+",".join(valid_chromosomes)) 
        my_error.append('Note: A one-to-one mapping is required, consider changing sumstats or ld-panel chromosome names')  
        bridge_pop_error(my_error) 
        return 


    def verify_ss_data(self, ss1, ss2): 
        self.progress.quikout('\nBridgeDebugData(Lvl'+str(self.args.debug_level)+'):')
        pTypes = [] 
        for si,s in enumerate([ss1, ss2]): 
            if s is not None: 
                pk, pn, pt = s.pop_type, s.pop_name, s.phenoTYPE, 
                eR = min(s.WT),max(s.WT)
                self.progress.quikout('\nPOP'+str(si+1)+' '+pn+' ('+pk+')')
                try: 
                    self.progress.quikout('Phenotype: '+self.args.phenotype+' ('+pt+')') 
                    eStr = 'Effect Size Range: '+str(eR[0])+', '+str(eR[1]) 
                    self.progress.quikout(eStr) 
                    if pt == 'binary' and eR[0] > 0:  self.mps.error('Truncated Effect Size Range, log(odds) is required.') 
                except: 
                    self.progress.quikout('Phenotype: '+str("Not Declared")) 
                
                g1, g2 = [], [] 
                for nk,n in [['MATCH','Matching'],['SWAPREF','Ref/Alt Swap'],['REVCOMP','Reverse Complement'],['INVALID','Invalid Bases']]: 
                    if s.CK[nk] > 0 :        g1.append(str(s.CK[nk])+' ('+n+')') 
                    if s.CK['GENO_'+nk] > 0: g2.append(str(s.CK['GENO_'+nk])+' ('+n+')') 
                if si == 0: 
                    self.progress.quikout('Initial Genotype Variants: '+str(s.snp_list_len)+' (SNP LIST), '+str(s.genome_snps)+' (VALID GENOTYPED SNPS)') 
                    self.progress.quikout('Shared GWAS Variants:      '+', '.join(g1))  
                else: 
                    self.progress.quikout('Shared Genotype Variants: '+','.join(g2)) 
                    self.progress.quikout('Shared GWAS Variants:     '+', '.join(g1)) 
        return 




















