import sys,os, multiprocessing 
from .BridgeProgress  import BridgeProgress
from .BridgeSettings  import BridgeSettings
from .BridgeTools     import BridgeTools
from .BridgePipelines import BridgePipelines
from collections import defaultdict as dd 



class BridgeIO:
        def __init__(self,args,bridgedir,rundir, command_line):
            self.args, self.bridgedir, self.rundir, self.homepath = args, bridgedir, rundir, os.path.abspath(args.outpath) 
            self.paths = {'home': self.homepath, 'save': self.homepath+'/save', 'logs': self.homepath+'/logs', 'tmp': self.homepath+'/tmp'} 
            self.progress = BridgeProgress(args, command_line).initialize(self.paths['home'])  
            self.set_programs_and_defaults() 
             
        def set_programs_and_defaults(self):
            self.ROUNDS, self.FOUND, self.LOC, self.programs= 0, dd(bool), dd(lambda: 'NA'), {} 
            try:    
                import matplotlib, matplotlib.pyplot
                self.FOUND['matplotlib'] = True 
            except: self.args.noPlots = True 
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
            self.pipeline = BridgePipelines(self).verify_pipeline()  
            self.settings.update_inputs(self.pipeline.input_key)

        def initialize(self, module, cmd): 
            self.module, self.cmd = module, cmd 
            self.settings = BridgeSettings(self)
            self.check_requirements()
            
            if module in 'tools': BridgeTools(self).apply() 
            if module in 'tools' or self.cmd[0:3] in ['req']: self.progress.finish('Complete',FIN=True) 
            elif module == 'analyze':                         self.settings.check_analysis_data()  
            else:                                             self.settings.check_pop_data() 
            if module == 'check': self.progress.finish('Complete',FIN=True) 

            #if module == 'tools' or self.cmd[0:3] in ['req']: self.progress.finish('Complete',FIN=True) 
            #if module == 'analyze': self.settings.check_analysis_data()  
            #elif module == 'tools': BridgeTools(self).apply() 
            #else:                   self.settings.check_pop_data() 
            #if module == 'check': self.progress.finish('Complete',FIN=True) 
            self.progress.show_settings(self.settings) 
            self.progress.start_module(self.module, self.cmd, self.paths['home'])# .show_settings(self.settings) 
            self.pipeline = BridgePipelines(self).verify_pipeline() 
            self.settings.update_inputs(self.pipeline.input_key) 
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


