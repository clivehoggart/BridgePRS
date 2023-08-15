#!/usr/bin/env python3

import sys
from src.Util.BridgeIO       import BridgeIO
from src.Run.BridgeBase      import BridgeBase
from src.Run.BridgeMore      import BridgeMore
from src.Run.BridgeJobs      import BridgeJobs
from src.Util.BridgeProgress import BridgeProgress



def bridge_error(eString):
    if type(eString) in [list,tuple]:
        sys.stderr.write('\nBridgeError: '+eString[0]+'\n')                     
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeError: '+eString+'\n')                           
    sys.exit()





class BridgePRS:
        def __init__(self,args,bridgedir,rundir,command_line):
            self.args = args 
            self.io       =   BridgeIO(args, command_line, bridgedir, rundir) 
            
            self.io.initialize(self.args.module, self.args.cmd) 
            self.io.start_progress() 
            
            if self.args.module == 'easyrun':   self.easyrun() 
            elif self.args.module == 'analyze': self.analyze(self.args.cmd, self.args.results, PATH = self.io.paths['home']) 
            else:                               self.execute(self.io.pipeline) 
         



        def execute(self,pl):
            self.jobs, self.base = BridgeJobs(self), BridgeBase(self) 
            for i,command in enumerate(pl.commands): 
                self.io.progress.start_minor(pl.command_strings[i], REVEAL=[self.io.settings]) 
                if pl.FIN[command.upper()] and not self.args.repeatSteps: 
                    self.io.progress.update_minor('SKIPPING-JOB')
                    continue 
                elif command == 'clump':     self.jobs.run(self.base.run_clump, [[['chromosome'],[k]] for k in self.io.settings.chromosomes])
                elif command == 'eval':      self.jobs.run(self.base.run_eval, [[[],[]]]) 
                elif command == 'predict':   self.jobs.run(self.base.run_predict,  [[[],[]]]) 
                elif command == 'quantify':  self.jobs.run(self.base.run_quantify,  [[[],[]]])       
                elif command == 'optimize':  self.jobs.run(self.base.run_optimize,  [[[],[]]]) 
                elif command == 'prior':     self.jobs.run(self.base.run_prior,  [[[],[]]]) 
                pl.update(command)  
                self.io.progress.end() 
            self.collate(pl.commands[-1]) 
            self.io.progress.finish() 
            return




        def collate(self, cmd): 
            self.base.close_all() 
            if cmd == 'quantify' and not self.args.skipAnalysis: 
                self.analyze('result',[self.io.pipeline.progress_file], PATH = self.io.paths['run'])
            return 


        def analyze(self, cmd, prs_results, PATH):
            self.more     =   BridgeMore(self) 
            if cmd == 'result':
                self.io.progress.start_minor('Plotting Results') 
                self.more.run(cmd, prs_results, PATH) 
                self.io.progress.end() 
            else:
                self.io.progress.start_minor('Combining Results') 
                self.more.run(cmd, prs_results, PATH) 
                self.io.progress.end() 
            return 

        
        def easyrun(self): 
            modules = ['prs-single','build-model','prs-port','prs-prior'] 
            self.config = self.args.config 
            pop1, pop2 = self.args.pop_names 
            res_files = [] 
            
            for i,m in enumerate(modules): 
                self.args.module, self.args.cmd =  m, 'run' 
                if i == 1: self.args.popname, self.args.config = pop2, [self.args.pop_configs[1]] + self.config
                else:      self.args.popname, self.args.config = pop1, [self.args.pop_configs[0]] + self.config
                self.io.initialize(self.args.module, self.args.cmd) 
                self.io.start_progress() 
                self.execute(self.io.pipeline) 
                if i == 1:  self.args.model_file = str(self.io.pipeline.progress_file) 
                else:       res_files.append(str(self.io.pipeline.progress_file)) 
            
            self.args.results = res_files 
            self.args.module, self.args.cmd = 'analyze','combine' 
            self.io.initialize(self.args.module, self.args.cmd) 
            self.io.start_progress() 
            self.analyze(self.args.cmd, self.args.results, PATH = self.io.paths['home']) 


