#!/usr/bin/env python2.7

import sys,os 
import os
from .BridgeProgress  import BridgeProgress
from .BridgeSettings  import BridgeSettings
from .BridgePipelines import BridgePipelines



def bridge_error(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeIOError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgeIOError: '+eString+'\n')
    sys.exit(2) 

class BridgeIO:
        def __init__(self,args,bridgedir,rundir, command_line):
            
            self.args, self.bridgedir, self.rundir = args, bridgedir, rundir 
            self.progress = BridgeProgress(args, command_line) 
            self.set_programs() 
            self.paths = {'home': os.path.abspath(self.args.outpath), 'logs': os.path.abspath(self.args.outpath)+'/logs'} 

            
            
        
        def set_programs(self): 
            self.programs = {} 
            for i,(p,pn) in enumerate([[self.args.rpath,'--rPath'],[self.args.plinkpath,'--plinkPath']]): 
                if os.path.exists(p):                      mp = os.path.abspath(p) 
                elif os.path.exists(self.bridgedir+'/'+p): mp = os.path.abspath(self.bridgedir+'/'+p) 
                else:                                      bridge_error([pn+' '+p+' Does not exist']) 
                for f in os.listdir(mp): self.programs[f.split('.')[0]] = mp+'/'+f 

            if self.args.platform == 'mac' and self.args.plinkpath.split('/')[-1] == 'Xtra': self.programs['plink'] = self.programs['plink_mac'] 
            return 


        def initialize(self, module, cmd): 
            self.module, self.cmd = module, cmd 
            self.pipeline = BridgePipelines(self).verify() 
            self.settings = BridgeSettings(self).verify()
            return self 
        
        def start_progress(self): 
            self.progress.begin_module(self.module, self.cmd, self.paths['logs']) 
            self.progress.start_major(self.cmd)  
            self.progress.show_settings(self.settings)
            return 
