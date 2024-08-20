import sys, os, gzip, shutil
from collections import defaultdict as dd
from collections import Counter as cc 
from .Bridge_Pop.BData import BData 
from .Bridge_Pop.GenoPheno import GenoPheno
from .Bridge_Pop.GWAS import SumStats 
from .Bridge_Pop import PopTools as ptools  


class BridgePops: 
    def __init__(self,mps):
        self.mps, self.paths, self.target, self.base, self.names = mps, {}, None, None, [] 
    
    def parse(self): 
        args = self.mps.parse_args()
        argnames = [v for v in vars(args)] 
        self.system_tests(args) 
        self.setup_paths(args) 
        self.arglen, self.target, self.base, self.names = len(vars(args)), None, None, [None, None] 
        if 'config' in argnames and len(vars(args)['config']) > 0: 
            if args.module == 'build-model': self.base   = BridgePop(self.mps, args.config[0], args, self.paths, 'base') 
            else:                            self.target = BridgePop(self.mps, args.config[0], args, self.paths, 'target') 
            if len(args.config) == 2:        self.base = BridgePop(self.mps, args.config[1], args, self.paths, 'base', self.target) 
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


                    continue 

    def __str__(self):                                                                                                                                                                                                      
        rpr = '<bridgeParseObject> Args='+str(self.arglen)+'; PopData: '
        if self.target is not None: rpr += 'Target='+self.target.name
        if self.base is not None: rpr += 'Base='+self.base.name
        return rpr 




class BridgePop:                                                                                                                                                                                                            
    def __init__(self, mps, P, args, paths, pop_type, prevPop = None): 
        self.mps, self.args, self.paths, self.type = mps, args, paths, pop_type
        self.sumstats_fields = [pf.strip() for pf in P['sumstats_fields'].split(',')] 
        self.field_key       = {a: A for a,A in zip(['ID','REF','ALT','P','BETA'],self.sumstats_fields)}
        self.gen = {} 
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
        self.mps.multi_line_error(my_error) 
        return 














