import sys,os 
from .BridgeProgress  import BridgeProgress
from .BridgeSettings  import BridgeSettings
from .BridgePipelines import BridgePipelines
from collections import defaultdict as dd 



def bridge_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeIOError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeIOError: '+eString+'\n')
    sys.exit(2) 


class BridgeIO:
        def __init__(self,args,bridgedir,rundir, command_line):
            self.args, self.bridgedir, self.rundir = args, bridgedir, rundir 
            self.homepath = os.path.abspath(self.args.outpath) 
            self.paths = {'home': self.homepath, 'logs': self.homepath+'/logs', 'tmp': self.homepath+'/tmp'} 
            self.progress = BridgeProgress(args, command_line).initialize(self.paths['home'])  
            self.set_programs_and_defaults() 
            
        
        def set_programs_and_defaults(self):

            try:    import matplotlib, matplotlib.pyplot
            except: self.args.noPlots = True 
            self.ld_ref, self.programs = {}, {} 
            for i,(p,pn) in enumerate([[self.args.rpath,'--rPath'],[self.args.plinkpath,'--plinkPath']]): #[self.args.ldpath,'--ldPath']]): 
                if   os.path.exists(p):                      mp = os.path.abspath(p) 
                elif os.path.exists(self.bridgedir+'/'+p):   mp = os.path.abspath(self.bridgedir+'/'+p) 
                else:                                        bridge_error([pn+' '+p+' Does not exist']) 
                for f in os.listdir(mp): self.programs[f.split('.')[0]] = mp+'/'+f 
                #if pn == '--tgPath': self.args.tgpath = mp   
                #else: 
            
            if '1000G_ref' in os.listdir(self.bridgedir+'/data'): self.load_ld(self.bridgedir+'/data/1000G_ref') 
            elif '1000G_sample' in os.listdir(self.bridgedir+'/data'): self.load_ld(self.bridgedir+'/data/1000G_sample') 
            for p in self.args.ldpath: 
                if os.path.exists(p):  mp = os.path.abspath(p) 
                elif os.path.exists(self.bridgedir+'/'+p): mp = os.path.abspath(self.bridgedir+'/'+p) 
                else:                                      bridge_error(['--ldPath '+p+' Does not exist']) 
                self.load_ld(mp) 
            
            if self.args.platform == 'mac' and self.args.plinkpath.split('/')[-1] == 'Xtra': self.programs['plink'] = self.programs['plink_mac'] 
            return 



        def load_ld(self, ld_path): 

            ids, cands, key, k = [], [], {}, 0                                                                                                                                                                                                                                                                
            for f in os.listdir(ld_path): 
                if f.split('.')[-1] in ['bed','bim','fam']: cands.append(".".join(f.split('.')[0:-1]))                                                                                                                                                                                      
                elif 'ids' in f:                            ids.append([f.split('.')[0], ld_path+'/'+f]) 
            cands = list(set(cands))    
            while True:                                                                                                                                                                                                                                                                     
            	prefix = list(set([c[0:k] for c in cands]))[0]                                                                                                                                                                                                                              
            	ck1 = list(set([c[0:k+1] for c in cands]))                                                                                                                                                                                                                                  
            	if len(ck1) == 1: k+=1                                                                                                                                                                                                                                                      
            	else: break          
            for cand in cands:                                                                                                                                                                                                                                                              
                chr_cand = cand.split(prefix)[-1]                                                                                                                                                                                                                                           
                key[chr_cand] = ld_path+'/'+cand       
            for pop, pop_path in ids: self.ld_ref[pop.upper()] = [pop_path, ld_path+'/'+prefix, key] 
            return 

        
        def start_progress(self): 
            self.progress.begin_module(self.module, self.cmd, self.paths['home']) 
            self.progress.start_major(self.cmd)  
            self.progress.show_settings(self.settings)
            return
        


        def initialize(self, module, cmd): 
            self.module, self.cmd = module, cmd 
            self.settings = BridgeSettings(self)
            if self.module != 'check':  
                if self.module == 'easyrun' and len(self.args.pop) < 2:  bridge_error('Two population names (--pop) or config files are required for easyrun (command line') 
                elif self.module != 'easyrun' and (len(self.args.pop) != 1 or len(self.args.pop_config) > 1):                            bridge_error('Ati most one population is required on command line (--pop) or in a config file (--pop_config)') 
                self.pipeline = BridgePipelines(self).verify() 
                if self.module != 'easyrun': self.settings.verify(self.pipeline.input_key) 
            else: 
                self.start_progress()
                self.check_requirements() 
                if self.cmd[0:3] in ['req']: self.progress.finish() 
                if self.cmd[0:3] != 'req': 
                    if len(self.args.pop) == 0:                            bridge_error('A population name is required on command line (--pop) or in a config file (--pop_config)') 
                    if len(self.args.pop) >  1:                            bridge_error('Module check requires data from at most one population') 
                    self.progress.start_minor('Checking Population Data', FIN = False) 
                    self.settings.verify() 
                    self.settings.check() 
                self.progress.finish() 
                sys.exit() 
            return self 
        

        def check_requirements(self): 
            self.progress.start_minor('Checking Requirements',SKIP=True) 
            check_log = self.paths['tmp']+'/bridge.validation.out' 
            check_err = self.paths['tmp']+'/bridge.validation.err' 
            pp = self.programs['check_availability'] 
            os.system('echo -n \"plink_loc \"  > '+check_log+' 2> '+check_err) 
            os.system('which plinki >> '+check_log+' 2>> '+check_err) 
            os.system('echo \"\" >> '+check_log+' 2>> '+check_err) 
            os.system('echo -n \"python_loc \"  >> '+check_log+' 2>> '+check_err) 
            os.system('which python3 >> '+check_log+' 2>> '+check_err) 
            os.system('echo \"\" >> '+check_log+' 2>> '+check_err) 
            os.system('echo -n \"R_loc \"  >> '+check_log+' 2>> '+check_err) 
            os.system('which R >> '+check_log+' 2>> '+check_err) 
            os.system('echo \"*****\" >> '+check_log+' 2>> '+check_err) 
            os.system('R --version >> '+check_log+' 2>> '+check_err) 
            os.system('echo \"*****\" >> '+check_log+' 2>> '+check_err) 
            os.system('Rscript --vanilla '+pp+' CHECK '+check_log+' >> '+check_err+' 2>> '+check_err) 
            RK = self.read_requirement_log(check_log)       
            ans = self.progress.show_requirements(RK) 
            if ans.lower()[0] == 'y': 
                self.progress.start_minor('Installing R Packages', FIN=False) 
                os.system('Rscript --vanilla '+pp+' INSTALL '+check_log) 
             
            
        
        def read_requirement_log(self, check_log): 
            MAKEPLOT = True 
            try:    import matplotlib, matplotlib.pyplot
            except: MAKEPLOT = False 
            K, M, stars, f = dd(bool), [], 0, open(check_log) 
            K['matplotlib'] = MAKEPLOT 
            for line in f: 
                if len(line) < 2: continue 
                lp, i, LOC = line.split(), 0, 'NA'
                if lp[0] == '*****': 
                    stars +=1 
                    continue
                if stars == 0:
                    if len(lp) == 2: K[lp[0]] = lp[1] 
                    else:            K[lp[0]] = 'NA' 
                elif stars == 1: 
                    if len(lp) > 2 and lp[0] == 'R' and lp[1] == 'version': K['R_version'] = lp[2] 
                else: 
                    M.append(lp[0]) 
            K['R_missing'] = M 
            return K 
        
        

