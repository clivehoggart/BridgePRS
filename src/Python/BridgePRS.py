import sys
from .Util.BridgeIO       import BridgeIO
from .Util.BridgeRun      import BridgeRun 
from .Util.Bridge_Run.BridgeJobs import BridgeJobs 
from .Util.Bridge_Run.BridgeBase import BridgeBase

def bridge_error(eString):
    if type(eString) in [list,tuple]:
        sys.stderr.write('\nBridgeError: '+eString[0]+'\n')                     
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeError: '+eString+'\n')                           
    sys.exit()



class BridgePRS:
        def __init__(self, mps, command_line): 
            self.io   = BridgeIO(mps, command_line) 
            self.args = self.io.args 
            self.io.initialize(self.args.module, self.args.cmd)                         
            #try: 
            #    import matplotlib
            #except: 
            #    if not self.args.noplots: 
            #        bridge_error(['Matplotlib Not Found','*please install matplotlib or run program with --noPlots']) 
            if   self.args.module == 'tools':     BridgeRun(self).apply() 
            elif self.args.module == 'analyze':   self.analyze(self.args.cmd, self.args.result_files, PATH = self.io.paths['home']) 
            else: 
                if self.args.phenotype is None: bridge_error('Phenotype field required for this option (see above)') 
                if self.args.module == 'pipeline':  self.easyrun() 
                else:                               self.execute(self.io.pipeline) 
         
        def easyrun(self): 
            if self.args.module == 'pipeline' and not self.args.port: modules = ['prs-single','build-model','prs-prior'] 
            else:                                                     modules = ['prs-single','build-model','prs-port','prs-prior'] 
            res_files = [] 
            for i,m in enumerate(modules): 
                self.args.module, self.args.cmd =  m, 'run'
                if i == 1: self.io.pop = self.io.pop_data.base 
                else:      self.io.pop = self.io.pop_data.target 
                self.io.update(self.args.module, self.args.cmd)
                self.execute(self.io.pipeline) 
                if i == 1:  self.args.model_file = str(self.io.pipeline.progress_file) 
                else:       res_files.append(str(self.io.pipeline.progress_file)) 
            self.args.result_files = res_files 
            self.args.module, self.args.cmd = 'analyze','combine' 
            self.io.update(self.args.module, self.args.cmd)
            self.analyze(self.args.cmd, self.args.result_files, PATH = self.io.paths['home']) 

        def execute(self,pl):
            self.jobs, self.base = BridgeJobs(self), BridgeBase(self) 
            for i,command in enumerate(pl.commands): 
                self.io.progress.start_minor(pl.command_strings[i], RD=self.io)
                if len(pl.commands) > 1 and pl.FIN[command.upper()]: 
                    self.io.progress.write('SKIPPING-JOB')
                    continue 
                else:
                    if command == 'clump':     self.jobs.run(self.base.run_clump, [[['chromosome'],[k]] for k in self.io.pop.chromosomes])
                    elif command == 'beta':      self.jobs.run(self.base.run_beta, [[[],[]]]) 
                    elif command == 'predict':   self.jobs.run(self.base.run_predict,  [[[],[]]]) 
                    elif command == 'quantify':  self.jobs.run(self.base.run_quantify,  [[[],[]]])       
                    elif command == 'prior':     self.jobs.run(self.base.run_prior,  [[[],[]]]) 
                    pl.log_result(command)  
                    self.io.progress.end(RD=self.io) 
            if command == 'quantify' and not self.args.noplots: self.analyze('result',[self.io.pipeline.progress_file], PATH = self.io.paths['run'])
            self.io.progress.finish() 
            return 


        def analyze(self, cmd, prs_results, PATH):
            self.run  = BridgeRun(self).analyze(cmd, prs_results, PATH) 
            self.io.progress.end() 
            return 

        
