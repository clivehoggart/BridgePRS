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
        def __init__(self,args,command_line, bridgedir,rundir):
            self.args, self.bdir, self.sdir = args, bridgedir, bridgedir + '/src/Scripts'
            self.progress = BridgeProgress(args, command_line)
            
            self.pdir = self.bdir+'/../Rscripts'  
            self.programs = {f.split('.')[0]: self.pdir+'/'+f for f in os.listdir(self.pdir)} 
            #self.programs = {f.split('.')[0]: self.sdir+'/'+f for f in os.listdir(self.sdir)} 
            if self.args.platform == 'mac': self.programs['plink'] = self.programs['plink_mac'] 
            self.paths = {'home': os.path.abspath(self.args.outpath), 'logs': os.path.abspath(self.args.outpath)+'/logs'} 
            

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
