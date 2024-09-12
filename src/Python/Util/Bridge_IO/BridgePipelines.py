import sys,os,gzip
from collections import defaultdict as dd





class BridgePipelines:
    def __init__(self,io): 
        self.args, self.io, self.module, self.cmd = io.args, io, io.module, io.cmd
        self.eType, self.wType, self.FIN = 'BridgePipelineError:', 'BridgePipelineWarning:', dd(bool)
        self.pop = self.io.pop.name 
        self.create() 



    def create(self): 
        if    self.module in ['prs-single','build-model']: pn = self.pop  
        elif  self.module in ['prs-port','prs-prior']:     pn = self.pop + '-' + self.args.model_file.split('.')[1].upper() 
        else: return 
        self.io.paths['run'] = self.io.paths['home']+'/'+self.module+'_'+pn
        self.progress_file   = self.io.paths['run']+'/bridge.'+pn.lower()+'.'+self.module+'.result'
        if   self.cmd != 'run':                    self.commands = [self.cmd] 
        elif self.module == 'prs-single':          self.commands = ['clump','beta','predict','quantify']                                                                                                     
        elif self.module == 'build-model':         self.commands = ['clump','beta','predict','prior']                                                                                                                                    
        elif self.module == 'prs-port':            self.commands = ['predict','quantify'] 
        elif self.module == 'prs-prior':           self.commands = ['beta','predict','quantify'] 
        else:                                      self.commands = [] 
        self.add_dirs(self.io.paths['run'], self.commands) 
        if not os.path.isfile(self.progress_file): 
            w = open(self.progress_file,'w')  
            w.write('POP='+self.pop+'\n') 
            w.write('MODULE_NAME='+self.module+'\n')
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
   


    def verify_pipeline(self,pop=None):

        if self.module in ['easyrun','pipeline','analyze']: return self  
        self.command_strings = []
        for i,c in enumerate(self.commands):
            JN = ' ('+self.module+' '+c+')'
            if c == 'clump': self.command_strings.append('Clumping Pop Data'+JN) 
            elif c == 'beta': self.command_strings.append('Calculating SNP Weights'+JN) 
            elif c == 'predict': self.command_strings.append('Generating Polygenic Predictions'+JN) 
            elif c == 'quantify': self.command_strings.append('Quantifying PRS Result'+JN) 
            elif c == 'prior': self.command_strings.append('Saving SNP Priors'+JN) 
            else:            self.command_strings.append(c) 
        prelims = [] 
        for X1,X2 in self.progress_pair(): 
            xsp = X1.split('_') 
            if len(xsp) == 3 and xsp[0] == self.pop: 
                x_pop, x_job, x_val = xsp 
                if x_val == 'FIN':      self.FIN[x_job] = True
                elif x_val not in ['PREFIX','FILE']: self.io.progress.fail('Bad Key',ETYPE=self.eType) 
                else: prelims.append([x_pop, x_job, x_val, X2]) 
        
        if pop is not None: 
            for pname,job,val,res in prelims: 
                if pop.name == pname and self.FIN[job] and val == 'PREFIX': pop.gen[job.lower()] = res 
        return self

    
    def progress_pair(self):
        p_handle = open(self.progress_file) 
        progress_pair = [x.strip().split('=') for x in p_handle.readlines()]  
        p_handle.close() 
        return progress_pair 



    def log_result(self, d): 
        pop = self.io.pop
        prefix, path = self.pop+'_'+d, self.io.paths['run']+'/'+d 

        f_pairs = [x for x in self.progress_pair()] + [line.strip().split('=') for line in open(pop.config,'rt')]
        if pop.sumstats.VALID: f_pairs.extend([['SSFIELD_'+k,v] for k,v in pop.sumstats.fields.items()]) 
        if pop.genopheno.VALID: 
            f_pairs.append(['PHENOTYPE_FILES',",".join(pop.genopheno.files)])
            for x,y in pop.genopheno.fields.items():  f_pairs.append(['PHENOTYPE_FIELD-'+x,y]) 
            f_pairs.append(['PHENOTYPE_TYPE',pop.genopheno.type]) 
        self.validate_result(prefix, path, d.upper()) 
        pop.gen[d] = path+'/'+prefix 
        for k,v in pop.gen.items(): 
            if k in ['model']: f_pairs.append([k.upper()+'_FILE',v]) 
            else:              f_pairs.extend([[self.pop+'_'+k.upper()+'_PREFIX',v],[self.pop+'_'+k.upper()+'_FIN','TRUE']])             
        f_obs, w = [], open(self.progress_file,'w') 
        for a,b in f_pairs: 
            if a not in f_obs: w.write(a+'='+b+'\n') 
            f_obs.append(a)
        w.close() 
        return 



    def validate_result(self, prefix, path, D): 
        PK, RK = dd(list), dd(list) 
        if not os.path.isdir(path): self.io.progress.broke('No Working Directory: '+D) 
        for f in os.listdir(path):
            fe = f.split('.')[-1] 
            if prefix in f: PK[fe].append(f) 
            elif fe == 'stderr': RK['err'].append(f) 
            elif fe == 'stdout': RK['out'].append(f) 
            else: continue 

        self.check_error_output(path, RK['err'], RK['out'], PK, D) 
        if D == 'CLUMP':
            fins,logs = [f.split('.')[0] for f in PK['gz']], [f.split('.')[0] for f in PK['log']]
            if len(PK['gz']) < len(PK['log']): 
                missing = ",".join([str(x) for x in sorted([int(f.split('_')[-1]) for f in logs if f not in fins])])   
                self.io.progress.broke('no results created for chromosomes: '+missing, ['     Note: Consider running with custom clump_value']) 
            for f in PK['gz']:
                cI = f.split('.')[0].split('_')[-1] 
                with gzip.open(path+'/'+f,'rt') as fh:
                    fh.readline() 
                    for i,line in enumerate(fh): 
                        line = line.split() 
                        if i > 2 or len(line) == 0: break  
                if i < 2:  self.io.progress.broke('Insufficient Clumping On Chromosome: '+cI, ['     Note: Consider running with custom clump_value']) 
        return

        




    def check_error_output(self, fp, err_files, out_files, run_files, D): 
        if len(run_files) == 0: self.io.progress.broke('No results created for job: '+D) 
        if len(out_files) != 1: self.io.progress.broke('No stdout created for job: '+D) 
        if len(err_files) != 1: self.io.progress.broke('No stderr created for job: '+D) 
        f_errors, f_handle = [], open(fp+'/'+err_files[0]) 
        f_lines = [lp.strip() for lp in f_handle] 
        f_handle.close()
        for lp in f_lines: 
            ls = [x.lower() for x in lp.split()] 
            if len(ls) < 1: continue 
            elif ls[0][0:4] in ['warn','extr','load','usag','lade']: continue  
            else: 
                if D == 'CLUMP': 
                    if " ".join(ls[-3::]) == 'see log file.': continue  
                    if len(ls) > 3 and ls[3] == 'ignored,':   continue 
                if len(ls[0].split('/'))>1: ls[0] = ls[0].split('/')[-1] 
                f_errors.append(' '.join(ls)) 
        if len(f_errors) == 0: return
        else:                  self.io.progress.warn(['UNKNOWN R-OUTPUT IN STDERR (Program may have failed):']+f_errors, WTYPE=self.wType) 
        return
        







