import sys,os
from collections import defaultdict as dd





def bridge_error(eString):
    if type(eString) in [list,tuple]:
        sys.stderr.write('\nBridgePipeLineError: '+eString[0]+'\n')                     
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: sys.stderr.write('\nBridgePipeLineError: '+eString+'\n')                           
    sys.exit()





class BridgePipelines:
    def __init__(self,io): 
        self.args, self.io, self.module, self.cmd = io.args, io, io.module, io.cmd
        self.input_key = dd(lambda: dd(bool)) 
        if self.module != 'analyze': self.create() 


    def create(self): 
        self.commands, self.pop = [self.cmd], self.args.popname 
        
        
        if self.module != 'easyrun': 
        
            if self.module.split('-')[0] == 'prs': 
                self.io.paths['run'] = self.io.paths['home']+'/'+self.module+'_'+self.pop 
                if self.cmd == 'run': 
                    if self.module.split('-')[-1] == 'port':      self.commands = ['predict','quantify'] 
                    else:                                         self.commands = ['clump','eval','predict','quantify'] 
                else:                                             self.commands = [self.cmd] 
            elif self.module == 'build-model': 
                self.io.paths['run'] = self.io.paths['home']+'/model_'+self.pop 
                if self.cmd == 'run': self.commands = ['clump','eval','optimize','prior'] 
                else:                 self.commands = [self.cmd] 
            else: bridge_error('Unsupported Module') 
        
            self.progress_file = self.io.paths['run']+'/bridge.'+self.module+'.result'
            self.add_dirs(self.io.paths['run'], self.commands) 
            if not os.path.isfile(self.progress_file): 
                w = open(self.progress_file,'w')  
                w.write('POP_NAME='+self.pop+'\n') 
                w.write('MODULE_NAME='+self.module+'\n')
                w.close()
            return  
        else: 
            self.pop1, self.pop2 = self.args.pop_names 
            for i,x in enumerate(['prs-single','prs-prior','prs-port','prs-combined']): 
                dir1 = x+'_'+self.pop1 
                if i < 2: self.add_dirs(self.io.paths['home']+'/'+dir1, ['clump','eval','predict','quantify'])
                if i < 3: self.add_dirs(self.io.paths['home']+'/'+dir1, ['predict','quantify'])
                else:     self.add_dirs(self.io.paths['home']+'/'+dir1, []) 
                progress_file = self.io.paths['home']+'/'+dir1+'/bridge.'+x+'.result'
                if not os.path.isfile(progress_file): 
                    w = open(progress_file,'w')  
                    w.write('POP_NAME='+self.pop1+'\n') 
                    w.write('MODULE_NAME='+x+'\n')
                    w.close()
            
            self.add_dirs(self.io.paths['home']+'/model_'+self.pop2, ['clump','eval','optimize','prior']) 
            progress_file = self.io.paths['home']+'/model_'+self.pop2+'/bridge.build-model.result'
            if not os.path.isfile(progress_file): 
                w = open(progress_file,'w')  
                w.write('POP_NAME='+self.pop2+'\n') 
                w.write('MODULE_NAME=build-model\n')
                w.close()
            return 
       

    def add_dirs(self,parent,children,grandchildren=[]):
        if not os.path.exists(parent): os.makedirs(parent)
        for c in children:
            if not os.path.exists(parent+'/'+c): os.makedirs(parent+'/'+c)
            self.io.paths[c] = parent+'/'+c 
            if c == 'prior': 
                if not os.path.exists(parent+'/'+c+'/lambda'): os.makedirs(parent+'/'+c+'/lambda') 
            for g in grandchildren:
                if not os.path.exists(parent+'/'+c+'/'+g): os.makedirs(parent+'/'+c+'/'+g)
                self.io.paths[g] = parent+'/'+c+'/'+g 
    


    def verify(self):
        if self.module in ['easyrun','analyze']: return self  
        self.command_strings = [] 
        for i,c in enumerate(self.commands):
            if c == 'clump': self.command_strings.append('Clumping Population Data') 
            elif c == 'eval': self.command_strings.append('Evaluating SNP Data') 
            elif c == 'predict': self.command_strings.append('Generating Polygenic Predictions') 
            elif c == 'quantify': self.command_strings.append('Quantifying PRS Result') 
            elif c == 'optimize': self.command_strings.append('Optimizing Prior PRS Model') 
            elif c == 'prior': self.command_strings.append('Saving SNP Priors') 
            else:            self.command_strings.append(c) 
       
        self.FIN = dd(bool)
        for X1,X2 in self.progress_pair(): 
            try: x_pop, x_job, x_val = X1.split('_') 
            except: continue
            if x_pop == self.pop: # and x_job.lower() != self.commands[-1]: 
                if x_val == 'FIN': self.FIN[x_job] = True  
                elif x_val == 'PREFIX': self.input_key['prefix'][x_job.lower()] = X2 
                elif x_val == 'FILE': self.input_key['file'][x_job.lower()] = X2 
                else:                bridge_error('Bad Key')  
        
        return self

    
    def progress_pair(self): 
        p_handle = open(self.progress_file) 
        progress_pair = [x.strip().split('=') for x in p_handle.readlines()]  
        p_handle.close() 
        return progress_pair 



    def update(self, d): 
        f_pass,f_obs = dd(bool), [] 
        np, fp, D = self.args.popname+'_'+d, self.io.paths['run']+'/'+d, d.upper() 
        f_pairs = [[X1,X2] for X1,X2 in self.progress_pair() if X1.split('_')[0] != self.pop or X1.split('_')[1] != D] 
        
        for x in ['model','pheno','validation']: 
            if x in self.io.settings.files and self.io.settings.files[x] is not None and x.upper()+"_FILE" not in [x[0] for x in f_pairs]: f_pairs.append([x.upper()+"_FILE",self.io.settings.files[x]]) 
            if x == 'pheno': 
                for kk,vv in self.io.settings.fields['pheno'].items(): 
                    a,b = 'FIELD_PHENO-'+kk, vv 
                    if a not in [x[0] for x in f_pairs]: f_pairs.append([a,b]) 

        self.validate_path(np,fp, D) 
        self.io.settings.prefixes[d] = fp+'/'+np  
        f_pairs.extend([[self.pop+'_'+D+'_PREFIX',fp+'/'+np],[self.pop+'_'+D+'_FIN','TRUE']])     
        


        #if not os.path.isdir(fp) or len(os.listdir(fp)) < 4: self.io.progress.fail('EXITED UNSUCCESSFULLY') #,[self.args.module,self.args.cmd,D])     
        #d_files =   [fn for fn in os.listdir(fp) if fn.split('.')[-1] != 'log'] 
        #self.check_error_output([fn for fn in d_files if fn.split('.')[-1] == 'stderr'])
        #res_files = [fn for fn in d_files if np in fn]
        #self.check_error_output(fp, err_files, D) 
        self.io.settings.prefixes[d] = fp+'/'+np  
        f_pairs.extend([[self.pop+'_'+D+'_PREFIX',fp+'/'+np],[self.pop+'_'+D+'_FIN','TRUE']])     
        w = open(self.progress_file,'w') 
        for a,b in f_pairs: 
            if a not in f_obs: w.write(a+'='+b+'\n') 
            f_obs.append(a) 
        w.close() 
        return           



    def validate_path(self, np, fp, D): 
        if not os.path.isdir(fp): self.io.progress.fail('EXITED UNSUCCESSFULLY') #,[self.args.module,self.args.cmd,D])     
        d_files =   [fn for fn in os.listdir(fp) if fn.split('.')[-1] != 'log'] 
        err_files, res_files = [fn for fn in d_files if fn.split('.')[-1] == 'stderr'], [fn for fn in d_files if np in fn]
        self.check_error_output(fp, err_files, D) 
        if len(res_files) < 1: self.io.progress.fail('EXITED UNSUCCESSFULLY; NO OUTPUT CREATED') #,[self.args.module,self.args.cmd,D])     
        return 



    def check_error_output(self, fp, err_files, D): 
        
        if len(err_files) != 1: self.io.progress.fail('EXITED UNSUCCESSFULLY') #,[self.args.module,self.args.cmd,D])     
        
        fname = fp+'/'+err_files[0] 
        f_handle = open(fname)  
        f_lines = [lp.strip() for lp in f_handle] 
        f_handle.close() 
        f_errors = [] 

        for lp in f_lines: 
            ls = [x.lower() for x in lp.split()] 
            if len(ls) < 1: continue 
            elif ls[0][0:4] in ['warn','extr','load','usag']: continue 
            else: 
                if D == 'CLUMP' and " ".join(ls[-3::]) == 'see log file.': continue  
                elif D == 'QUANTIFY': continue                  
                if len(ls[0].split('/'))>1: ls[0] = ls[0].split('/')[-1] 
                f_errors.append(' '.join(ls)) 
        
        if len(f_errors) == 0: return 
        self.io.progress.fail('EXITED UNSUCCESSFULLY',f_errors) #,[self.args.module,self.args.cmd,D])     







