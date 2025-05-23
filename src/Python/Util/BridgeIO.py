import sys,os, multiprocessing 
from .BridgePops import BridgePops 
from .Bridge_IO.BridgeProgress  import BridgeProgress
from .Bridge_IO.BridgePipelines import BridgePipelines
from collections import defaultdict as dd 

class BridgeIO:
        def __init__(self, mps, command_line):
            self.pop_data, self.pop = BridgePops(mps).parse(), None 
            self.args, self.paths = self.pop_data.args, self.pop_data.paths  
            self.progress = BridgeProgress(self.args, command_line).initialize(self.paths['home'], [n for n in self.pop_data.names if n != None]) 
            self.pop_data.set_progress(self.progress) 
            self.set_programs_and_defaults()

        def set_programs_and_defaults(self):
            self.ROUNDS, self.FOUND, self.LOC, self.programs= 0, dd(bool), dd(lambda: 'NA'), {} 
            try:    
                import matplotlib, matplotlib.pyplot
                self.FOUND['matplotlib'] = True 
            except: self.args.noplots = True
            for f in os.listdir(self.args.rpath):     self.programs[f.split('.')[0]] = self.args.rpath+'/'+f 
            for f in os.listdir(self.args.plinkpath): self.programs[f.split('.')[0]] = self.args.plinkpath+'/'+f 
            self.find_which('plink') 
            if self.FOUND['plink']:                                                            self.programs['plink'], self.args.plinkpath = self.LOC['plink'], self.LOC['plink']
            elif self.args.platform == 'mac' and self.args.plinkpath.split('/')[-1] == 'Xtra': self.programs['plink'] = self.programs['plink_mac']
            self.LOC['plink'] = self.programs['plink']  
            return 

        def find_which(self, name): 
            p_log, p_err = self.paths['tmp']+'/tmp.'+name+'.out', self.paths['tmp']+'/tmp.'+name+'.err' 
            os.system('which '+name+' > '+p_log+' 2> '+p_err) 
            f = open(p_log, 'rt') 
            p_path = f.readline().split() 
            f.close() 
            if len(p_path) == 1 and p_path[0].split('/')[-1] == name: self.FOUND[name], self.LOC[name] = True, p_path[0] 
            return             
            
        def update(self, module, cmd): 
            self.module, self.cmd = module, cmd 
            self.progress.start_module(self.module, self.cmd, self.paths['home']) 
            if self.pop is not None: self.pop.validate_args(self.args, module, cmd, self.progress)  
            self.pipeline = BridgePipelines(self).verify_pipeline(self.pop)   

        def initialize(self, module, cmd): 
            self.module, self.cmd = module, cmd 
            self.check_requirements()
            if module in ['tools','analyze']: return self 
            
            self.pop_data.load_from_configs() 
            self.progress.show_pop_data([self.pop_data.target,self.pop_data.base]) 
            self.progress.show_settings()     
            if self.module == 'pipeline': return self 
            elif self.module == 'build-model': self.pop = self.pop_data.base 
            else:                              self.pop = self.pop_data.target 
            self.update(self.module, self.cmd) 
            return self


        def check_requirements(self): 
            R_data, pp = [], self.programs['check_availability'] 
            self.find_which('python3') 
            self.find_which('R') 
            p_cmd = self.programs['plink']+' > '+self.paths['tmp']+'/tmp.pk.out  2> '+self.paths['tmp']+'/tmp.pk.out'
            if self.FOUND['R']: 
                RV_log, RV_err = self.paths['tmp']+'/tmp.RV.out', self.paths['tmp']+'/tmp.RV.err' 
                RP_log, RP_err = self.paths['tmp']+'/tmp.RP.out', self.paths['tmp']+'/tmp.RP.err' 
                os.system('R --version > '+RV_log+' 2> '+RV_err) 
                os.system('Rscript --vanilla '+pp+' CHECK '+RP_log+' > '+RP_log+' 2> '+RP_err) 
                f = open(RV_log)
                lp = f.readline().split() 
                if len(lp) > 2 and lp[0] == 'R' and lp[1] == 'version': R_data = [lp[2]] 
                else:                                                   R_data = ['unknown'] 
                f.close() 
                f = open(RP_log)
                R_data.append([x.strip() for x in f.readlines() if len(x) > 1]) 
                f.close() 
            self.progress.show_requirements(self.FOUND, self.LOC, R_data, p_cmd) 
            return





