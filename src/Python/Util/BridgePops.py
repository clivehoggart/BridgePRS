import sys, os, gzip, shutil
from collections import defaultdict as dd
from collections import Counter as cc 
from .Bridge_Pop.BData import BData 
from .Bridge_Pop.GenoPheno import GenoPheno
from .Bridge_Pop.GWAS import SumStats 
from .Bridge_Pop import PopTools as ptools  

# multi_line_error

class BridgePops: 
    def __init__(self,mps):
        self.mps, self.paths, self.target, self.base, self.progress, self.names = mps, {}, None, None, None, [] 
   
    def set_progress(self, progress): self.progress = progress

    def parse(self): 
        self.args = self.mps.parse_args()
        self.argnames = [v for v in vars(self.args)] 
        self.system_tests(self.args) 
        self.setup_paths(self.args) 
        self.arglen, self.target, self.base, self.names = len(self.argnames), None, None, [None, None] 
        if 'config' in self.argnames and len(vars(self.args)['config']) > 0: self.names = [P['pop'] for P in self.args.config] 
        elif 'pop' in self.args:                                             self.names = [x for x in self.args.pop if x is not None] 
        else:                                                                self.names = [] 
        return self 



    def load_from_configs(self): 
        if 'config' in self.argnames and len(vars(self.args)['config']) > 0: 
            #if self.args.module == 'build-model': self.base   = BridgePop(self.mps, self.args.config[0], self.args, self.paths, 'base') 
            if self.args.module == 'build-model': self.base   = BridgePop(self, self.args.config[0], 'base') 
            else:                                 self.target = BridgePop(self, self.args.config[0], 'target') 
            #if len(self.args.config) == 2:        self.base = BridgePop(self.mps, self.args.config[1], self.args, self.paths, 'base', self.target) 
            if len(self.args.config) == 2:        self.base = BridgePop(self, self.args.config[1], 'base', self.target) 
            if self.args.restart: self.restart_directories(self.args)  
        return


    def load_from_args(self,args): 
        P1, P2 = {'MISSING': dd(bool)}, {'MISSING':  dd(bool)} 
        



        for v in vars(args):
            v1, v2, k = v.split('_')[0], v.split('_')[-1], vars(args)[v] 
            if type(k) != list or len(k) == 0: continue 
            if v1 != 'sumstats' and v2 not in  ['file','prefix'] and v not in ['pop','ldpop','ld_path','max_clump_size','covariates','clump_value']: continue
            P1[v], P2[v] = k[0], k[-1] 
        if 'pop' not in P1: self.mps.pop_error('Population Name(s) Are Required, use --pop')  
        if 'ldpop' not in P1: P1['ldpop'],P2['ldpop'] = P1['pop'], P2['pop'] 
        if 'ld_path' not in P1: self.mps.pop_error('LD Path(s) Are Required, use --ld_path')  
        if 'genotype_prefix' not in P1: self.mps.pop_error('A Genotype Prefix is Required, use --genotype_prefix')  
        if 'phenotype_file' not in P1: self.mps.pop_error('A Phenotype File is Required, use --phenotype_file')  
        for k in ['validation_file','thinned_snp_file','max_clump_size','covariates','sumstats_suffix','sumstats_size','snp_file','clump_value']: 
            if k not in P1: 
                P1[k], P1['MISSING'][k], P2[k], P2['MISSING'][k] = None, True, None, True
        

        self.target = BridgePop(self, P1, 'target') 
        if self.target.name != P2['pop']: self.base = BridgePop(self, P2, 'base', self.target) 
        self.names = [x.name for x in [self.target, self.base] if x is not None] 
        if args.restart: self.restart_directories(args)  
        return args, self 




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
        self.paths['home'] = args.outpath 
        for d in ['logs','tmp','save']: 
            if not os.path.isdir(args.outpath+'/'+d): os.makedirs(args.outpath+'/'+d)    
            self.paths[d] = args.outpath+'/'+d 
        return 
    
    def restart_directories(self, args): 
        for f in os.listdir(args.outpath): 
            if os.path.isdir(args.outpath+'/'+f) and len(f.split('_')) > 1:
                f_job, fn = f.split('_')[0], f.split('_')[-1].split('-') 
                if args.module == 'pipeline': 
                    if len(fn) == 2 and (fn[0] == self.names[0] and fn[1] == self.names[1]): shutil.rmtree(args.outpath+'/'+f) 
                    elif len(fn) == 1 and ((f_job == 'prs-single' and fn[0] == self.names[0]) or (f_job == 'build-model' and fn[0] == self.names[1])): shutil.rmtree(args.outpath+'/'+f)  
                    else: continue 
                elif args.cmd == 'run':                     
                    if f_job == args.module: 
                        if len(fn) == 1 and fn[0] == self.names[0]: shutil.rmtree(args.outpath+'/'+f) 
                        elif len(fn) == 2: 
                            print(f_job, 'hmmm') 
                            sys.exit() 
                        else: continue 
                else: 
                    print('wtf') 
                    continue 

    def __str__(self):                                                                                                                                                                                                      
        rpr = '<bridgeParseObject> Args='+str(self.arglen)+'; PopData: '
        if self.target is not None: rpr += 'Target='+self.target.name
        if self.base is not None: rpr += 'Base='+self.base.name
        return rpr 


    #else:                                 self.target = BridgePop(self, self.args.config[0], 'target') 
    #        if len(self.args.config) == 2:        self.base = BridgePop(self.mps, self.args.config[1], self.args, self.paths, 'base', self.target) 


class BridgePop:                                                                                                                                                                                                            
    def __init__(self, pd, P, pop_type, prevPop = None): 
        self.mps, self.args, self.paths, self.progress, self.type = pd.mps, pd.args, pd.paths, pd.progress, pop_type

        self.sumstats_fields = [pf.strip() for pf in P['sumstats_fields'].split(',')] 
        self.field_key       = {a: A for a,A in zip(['ID','REF','ALT','P','BETA'],self.sumstats_fields)}
        self.gen = {} 
        self.name, self.ref_pop = P.pop('pop'), P.pop('ldpop') 
        for k,v in P.items(): vars(self)[k] = v
        for k,v in P.items(): 
            if k == 'pop':           self.name = v 
            elif k == 'ldpop':       self.ref_pop = v 
            elif k == 'config_name': self.file = v 
            elif v is not None or prevPop is None:       vars(self)[k] = v 
            elif k.split('_')[0] != 'sumstats': vars(self)[k] = vars(prevPop)[k] 
            else:                             vars(self)[k] = v  
        self.bdata = BData(self).load_panel()
        self.genopheno = GenoPheno(self, self.args.phenotype).load()    
        self.sumstats  = SumStats(self).load(self.genopheno, prevPop) 
        




        self.verify_chromosomes() 
        self.validate_args(self.args, self.args.module, self.args.cmd, pd.mps) 

    def validate_args(self, args, module, cmd, progress): 
        
        self.gen, mInputs, mGen = {}, [], [] 
        self.module, self.cmd = module, cmd 
        for v in vars(args):
            kV = vars(args)[v] 
            if kV and v.split('_')[-1] in ['prefix','file']: self.gen[v.split('_')[0]] = kV 
            else: continue 
        
        if module in ['tools']: return True 
        if cmd in  ['beta']    and 'clump' not in self.gen: mGen.append('Clump Data (Hint: Run '+module+' clump)')
        if cmd in  ['predict'] and 'beta' not in self.gen and module not in ['prs-port']: mGen.append('Weight Data (Hint: Run '+self.args.module+' beta)')
        if cmd in ['quantify'] and 'predict' not in self.gen: mGen.append('Pred Data (Hint: Run '+module+' predict') 
        if module in ['prs-port','prs-prior']: 
            if 'model' not in self.gen: mGen.append('Base Model --model_file (Hint: Run build-model)') 
            else:
                model_pairs = [line.strip().split('=') for line in open(self.gen['model'],'rt') if line.strip().split('=')[0].split('_')[-1].upper() in ['SIZE','FIN']]  #'SUMSTATS_SIZE']      
                if len(model_pairs) != 5: mGen.append('Invalid Base Model --model_file (Hint: ReRun build-model)') 
        if len(mGen) > 0: progress.fail(['Missing/Invalid Input Data:']+mGen)
        return True


    def verify_chromosomes(self): 
        chromosomes       = list(set([k for k in self.sumstats.map.keys()]+[k for k in self.bdata.map.keys()]))
        ld_chrs = ptools.get_chr_strs(self.bdata.map.keys()) 
        ss_chrs = ptools.get_chr_strs(self.sumstats.map.keys())
        valid_chromosomes = [k for k in self.sumstats.map.keys() if k in self.bdata.map.keys()] 
        self.chromosomes       = ptools.get_chr_strs(list(set([k for k in self.sumstats.map.keys()]+[k for k in self.bdata.map.keys()])))
        if len(self.chromosomes) == len(ld_chrs) and len(self.chromosomes) == len(ss_chrs) and len(self.chromosomes) == len(valid_chromosomes): return
        my_error = ['Missing Chromosomes'] 
        my_error.append('The Following Chromosomes Are Found in the ld-panel: '+",".join(ld_chrs)) 
        my_error.append('The Following Chromosomes Are Found in the sumstats: '+",".join(ss_chrs))
        if len(valid_chromosomes) == 0: my_error.append('The Following Chromosomes Are Found in both sources: None') 
        else:                           my_error.append('The Following Chromosomes Are Found in both sources: '+",".join(valid_chromosomes)) 
        my_error.append('Note: A one-to-one mapping is required, consider changing sumstats or ld-panel chromosome names')  
        self.mps.fail(my_error) 
        return 














