import sys, os 
from collections import defaultdict as dd
import time 

LOCALTIME = time.asctime( time.localtime(time.time()) )


def req_error(eString): 
    if type(eString) in [list,tuple]:
        sys.stderr.write('\n\n      RequirementError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('    '+es+'\n')
    else: 
        sys.stderr.write('\n\n      RequirementError: '+eString+'\n')
    sys.exit(2)




class BridgeProgress:
    def __init__(self,args, command_line): 
        self.obs_input = [] 
        self.INIT, self.active, self.loud = True, True, False 
        if args.verbose:  self.active, self.loud = True, True 
        elif args.silent: self.active, self.loud = False, False 
        else:             self.active, self.loud = True, False 
        self.args, self.command_line = args, command_line  
        self.out = sys.stderr
        self.jobs = ['init'] 
        

    def initialize(self, runpath):
        self.runpath, self.INIT = runpath, False
        self.logfile = runpath + '/logs/bridgePRS.'+self.args.module+'.'+self.args.cmd+'.log' 
        self.loghandle = open(self.logfile,'w') 
        self.FILE, self.log = True, [[]] 
        self.write('BridgePRS Begins at '+LOCALTIME+' \n') 
        if self.command_line != None:  self.write('Bridge Command-Line: '+' '.join(self.command_line)+'\n')
        return self 
    





    def begin_module(self, module, cmd, runpath):
        self.module, self.cmd, self.runpath = module, cmd, runpath
        if self.INIT: 
            self.logfile = runpath + '/logs/bridgePRS.'+self.args.module+'.'+self.args.cmd+'.log' 
            self.loghandle = open(self.logfile,'w') 
            self.FILE, self.log = True, [[]] 
            self.write('BridgePRS Begins at '+LOCALTIME+'\n') 
            if self.command_line != None:  self.write('Bridge Command-Line: '+' '.join(self.command_line)+'\n')
            self.INIT = False  
        self.update_indent, self.update_len = '', 0 
        self.start_heading = 'Module: '
        self.job_heading   = 'Command: ' 
        self.write('\n'+self.start_heading+self.module+'\n') 
        self.h_blank = ''.join([' ' for x in range(len(self.start_heading))])
        self.sub_blank = ' '
        self.topic, self.subtopic, self.status, self.topics, self.events = None, None, None, 0, {} 

    def write(self, outstring): 
        if self.active: 
            self.out.write(outstring) 
            self.out.flush()  
        if self.FILE: self.loghandle.write(outstring) 

    def whisper(self, outstring): 
        if self.loud: 
            self.out.write(outstring) 
            self.out.flush() 
        if self.FILE: self.loghandle.write(outstring) 



    def missing_str_list(self, X): 
        if len(X) == 0: return '' 
        ki, chr_keep = 1, [[X[0]]] 
        while ki < len(X): 
                try:
                    while ki < len(X) and int(X[ki]) == int(chr_keep[-1][-1])+1:
                        chr_keep[-1].append(X[ki]) 
                        ki+=1 
                except: pass
                if ki < len(X): 
                    chr_keep.append([X[ki]])
                    ki+=1 
        chr_str = [] 
        for x in chr_keep: 
            if len(x) == 1: chr_str.append(x[0]) 
            else:           chr_str.append(x[0]+'-'+x[-1]) 
        return ('['+",".join(chr_str)+']') 




    def show_settings(self,settings):
        blank2 = self.h_blank+'  ' 
        kFiles, kPrefixes, kSkip = [], [], ['config','prefix','file','module','cmd','dataset','ssf','pf','sumstats_suffix','phenotype_files'] 
        if 'files' in vars(settings):
            kFiles     = [[a,b] for a,b in settings.files.items() if b is not None and a not in ['clumpz']] 
            kPrefixes  = [[a,b] for a,b in settings.prefixes.items() if b is not None and a not in ['clumpz']]
            if len(kFiles) > 0: 
                for i,(k,fn) in enumerate(sorted(kFiles)):
                    self.obs_input.append(k) 
                    if i == 0:  self.write(blank2+'       Files:  '+k+'_file='+fn+'\n') 
                    else:       self.write(blank2+'               '+k+'_file='+fn+'\n') 
            if len(kPrefixes) > 0: 
                for i,(k,fn) in enumerate(sorted(kPrefixes)): 
                    self.obs_input.append(k) 
                    if i == 0:  self.write(blank2+'    Prefixes:  '+k+'_prefix='+fn+'\n') 
                    else:       self.write(blank2+'               '+k+'_prefix='+fn+'\n') 
        elif settings.args.module == 'check':
            FSKIP = ['module','cmd','dataset'] 
        
        elif 'pop_config' in settings.args:
            pn = vars(settings.args)['pop'] 
            pc = vars(settings.args)['pop_config'] 
            self.write(blank2+' Populations:  '+", ".join(pn)+'\n')  
            self.write(blank2+'Config Files:  '+", ".join(pc)+'\n')  
            
        kFalse, kTrue, kNone, kPlatform, kOpts, kInputs, kObs = [] , [], [] , [], [], [], [] 
        for k in vars(settings.args):
            if k.split('_')[-1] in kSkip or k.split('-')[0] in kSkip or k in kSkip: continue 
            v = vars(settings.args)[k]
            kObs.append(k) 
            if k in ['cores','platform','plinkpath','rpath']: kPlatform.append([k,v])  
            elif str(v) == 'True':     kTrue.append([k,v]) 
            elif str(v) == 'False':    kFalse.append([k,v]) 
            elif v in [None,0,'0']:    kNone.append([k,v]) 
            elif type(v) != list:      kOpts.append([k,v])              
            else:                      kInputs.append([k,v])        
        
        kOpts.sort() 
        if len(kPlatform) > 0: self.write(blank2+'      System:  '+',  '.join([k+'='+str(fn) for k,fn in sorted(kPlatform)])+'\n')
        if len(kTrue) > 0:     self.write(blank2+'       Flags:  '+',  '.join([k+'='+str(fn) for k,fn in sorted(kTrue)+sorted(kFalse)])+'\n')
        if len(kOpts) > 0:     self.write(blank2+'     Options:  '+',  '.join([k+'='+str(fn) for k,fn in kOpts[0:6]])+'\n')
        if len(kOpts) > 6:     self.write(blank2+'               '+',  '.join([k+'='+str(fn) for k,fn in kOpts[6::]])+'\n')
        if 'input' in vars(settings): self.show_inputs(settings)  




    def show_inputs(self,settings):
        ss, bf, ph = settings.input.sumstats, settings.input.bdata, settings.input.phenotypes
        blank2 = self.h_blank+'  ' 
        if ss.VALID:
            chr_str = self.missing_str_list(sorted(ss.chromosomes)) 
            f_data = ss.prefix+chr_str+ss.suffix 
            self.write(blank2+'    Sumstats:\n')
            self.write(blank2+'             Input Files:   '+f_data+'\n')
            self.write(blank2+'             SNP-QC File:   '+ss.snp_file+'\n')
            self.write(blank2+'             Input Fields:  '+", ".join([x+'='+y for x,y in ss.fields.items()])+'\n') 
        
        if bf.VALID: 
            chr_str = self.missing_str_list(sorted(bf.chromosomes)) 
            f_data = bf.prefix+chr_str+'.'+'[bed,bim.fam]'
            self.write(blank2+'      Bfiles:\n')
            self.write(blank2+'             Input Files:  '+f_data+'\n')
            self.write(blank2+'             LD_ID File:   '+bf.id_file+'\n')
            
        if ph.VALID: 
            self.write(blank2+'  Phenotypes:\n')
            
            pi, (p1,p2) = 0, [fx.split('/') for fx in [ph.files[0],ph.files[-1]]]
            while pi < min(len(p1),len(p2)) and p1[pi] == p2[pi]: pi+=1 
            if len(ph.files) == 1:   self.write(blank2+'             Test/Validation File: '+ph.files[0]+', None\n') 
            elif pi < 2:             self.write(blank2+'             Test/Validation File: '+ph.files[0]+', '+ph.files[1]+'\n') 
            else:                    self.write(blank2+'             Test/Validation File: '+ph.files[0]+', '+"/".join(p2[pi-1::])+'\n')
            if 2 < 5:                self.write(blank2+'             Input Fields:         '+", ".join([x+'='+y for x,y in ph.fields.items()])+'\n') 
        self.write('\n') 
        return

    def reveal_new_settings(self, settings): 
        kFiles, kPrefixes = [[a,b] for a,b in settings.files.items() if b is not None and a not in self.obs_input], [[a,b] for a,b in settings.prefixes.items() if b is not None and a not in self.obs_input]
        if len(kFiles) + len(kPrefixes) > 0: 
            self.write(self.h_blank+'               Data Generated:       \n') 
            if len(kFiles) > 0: 
                for i,(k,fn) in enumerate(sorted(kFiles)):
                    self.obs_input.append(k) 
                    if i == 0:  self.write(self.h_blank+'                        Files:      '+k+'_file='+fn+'\n') 
                    else:       self.write(self.h_blank+'                     '+k+'_file='+fn+'\n') 
            if len(kPrefixes) > 0: 
                for i,(k,fn) in enumerate(sorted(kPrefixes)): 
                    self.obs_input.append(k) 
                    if i == 0:  self.write(self.h_blank+'                     Prefixes:      '+k+'_prefix='+fn+'\n') 
                    else:       self.write(self.h_blank+'                     '+k+'_prefix='+fn+'\n') 
            self.write('\n') 
        return            






    def update_minor(self,update):
        if 'SKIP' in update.upper(): self.write('..'+update+'....')  
        else:                        self.write('....'+update) 







    def localize_paths(self, jn, rPaths): 
        rundir = self.runpath.split('/')[-1] 
        FI, FZ, LP = [], {}, {} 
        for a,b in rPaths: 
            n, bs = a[2::], b.split(self.runpath) 
            if n == 'bfile':                   FI    = [a,n,b]  
            elif len(bs) > 1:                  LP[n] = [a,rundir+bs[-1]]
            else:                              FZ[n] = [a,b] 
        K = [k for k in FZ.keys()] + [k for k in LP.keys()] 
        if len(FI) == 0: 
            for n,(a,b) in FZ.items(): LP[n] = [a,b] 
        else:  
            fA, fN, fB = FI 
            K = [fN] + K 
            f_fields = fB.split('/')
            LP[fN] = [fA, fB] 
            for n,(a,b) in FZ.items(): 
                bi, b_fields = 0, b.split('/') 
                bLen = min(len(b_fields), len(f_fields)) 
                while bi < bLen and b_fields[bi] == f_fields[bi]: bi+= 1  
                if bi > 2:  b = "*/"+"/".join(b_fields[bi-1::])
                LP[n] = [a,b] 
        return K, LP  



    def show_local(self, K, LP): 
        my_len = 0
        for i,k in enumerate(K):
            a,b = LP[k]
            if b[-1] == '1': b = b.strip('1')+'*' 
            if i == 0:   self.whisper(self.new_bk+a+' '+b+'\n') 
            elif i == 1:     
                self.whisper(self.new_bk+a+' '+b) 
                my_len = len(a)+len(b) 
            else: 
                if my_len < 90: 
                    self.whisper('  '+a+' '+b) 
                    my_len += len(a)+len(b) 
                else:           
                    self.whisper('\n'+self.new_bk+a+' '+b) 
                    my_len = len(a)+len(b) 
        self.whisper('\n') 
        return  


    def finish_data_validation(self, CONFIG_FILE): 

        self.new_bk = self.sub_blank+'      ' 
        self.write('\n'+self.new_bk+'Complete: Population Data Successfully Validated'+'\n') 
        self.write(self.new_bk+'Temporary Config File Location: '+CONFIG_FILE+'\n') 
        sys.exit() 



    def show_requirements(self, RK):
        self.new_bk = self.sub_blank+'      ' 
        self.write('..................................................\n')  
        #if not RK['python_loc']: req_error('Python not found, Please Install Python')  
        
        if RK['python_loc']: self.write(self.new_bk+'Python Found: '+RK['python_loc']+'\n') 
        else:                self.write(self.new_bk+'Python: Not Found, (Install Python to use wrapper)\n') 


        if RK['matplotlib']: self.write(self.new_bk+'Matplotlib: Found\n') 
        else:                self.write(self.new_bk+'Matplotlib: Not Found, (Install Matplotlib to produce plots)\n') 
        if not RK['R_loc']: req_error('R not found, Please Install R')  
        else:               
            self.write(self.new_bk+'R Found: '+RK['R_loc']+'\n') 
            self.write(self.new_bk+'R Version: '+RK['R_version']+'\n') 
        ANS = 'n' 
        if len(RK['R_missing']) == 0 : self.write(self.new_bk+'R packages: Up To Date\n') 
        else:                          
            self.write(self.new_bk+'R Packages -Missing: '+",".join(RK['R_missing'])+'\n') 
            self.write(self.new_bk+'Would you like to attempt to install these packages (y/n)?: ') 
            ANS = input() 
        return ANS            






    def start_rJob(self, RJ, job_name): 
        rPaths, rCols, rSave, rOne, rVal = [], [], [], [], [] 
        self.new_bk = self.sub_blank+'         ' 
        if job_name == 'clump':    rInit,rJ = RJ[0:1], RJ[1::] 
        elif RJ[1] != '--vanilla': rInit,rJ = RJ[0:2], RJ[2::] 
        else:                      rInit,rJ = RJ[0:3], RJ[3::] 
        for i in range(1,len(rJ), 2): 
            a,n,b = rJ[i-1], rJ[i-1][2::], rJ[i] 
            if   n in ['n.cores','fpath']: continue 
            elif job_name == 'clump': 
                if n in  ['clump-field','clump-snp-field']: continue 
                elif n in ['clump','out','bfile','extract','keep']: rPaths.append([a,b])          
                else:                                               rSave.append([a,b])                                                
            else: 
                if len(a.split('.')) > 1 and a.split('.')[0] == '--sumstats': rCols.append([a.split('.')[-1], b])
                elif   len(rJ[i].split('/'))   > 1:                           rPaths.append([a,b]) 
                elif   rJ[i] == '1':                                          rOne.extend([rJ[i-1], rJ[i]]) 
                else:                                                         rSave.append([rJ[i-1],rJ[i]])
        K, LP = self.localize_paths(job_name, rPaths) 
        if job_name == 'clump': 
            self.whisper('\n'+self.sub_blank+'Running Plink:\n') 
            K = ['bfile', 'clump','extract', 'keep', 'out']
            self.show_local(K, LP) 
            if len(rSave) > 0: 
                self.whisper(self.new_bk+rSave[0][0]+' '+rSave[0][1]) 
                for rs in rSave[1::]: self.whisper(' '+rs[0]+' '+rs[1]) 
                self.whisper('\n') 
        else:    
            if RJ[1] == '--vanilla': rInit,rJ = RJ[0:3], RJ[3::] 
            else:                    rInit,rJ = RJ[0:2], RJ[2::] 
            self.whisper('..................................................\n'+self.sub_blank+'Running Rscript:\n') 
            self.whisper(self.new_bk+" ".join(rInit)+'\n') 
            self.show_local(K, LP) 
            if len(rOne) > 0: self.whisper(self.new_bk+" ".join(rOne)+'\n') 
            if len(rSave) > 0: 
                self.whisper(self.new_bk+rSave[0][0]+' '+rSave[0][1]) 
                for rs in rSave[1::]: self.whisper(' '+rs[0]+' '+rs[1]) 
                self.whisper('\n') 
        self.whisper(self.new_bk+'........................................................................') 
        return  




    def start_major(self,topic,msg=None,subtopics=[]): 
        if self.status == 'minor': self.end() 
        self.status = 'major'
        self.write(self.h_blank+self.job_heading+topic+'\n') 
        self.blank = self.h_blank+''.join([' ' for x in range(len(self.job_heading))])
        self.sub_blank = self.blank 
        if msg is not None: 
            self.write(self.blank+msg+' '+', '.join([str(s) for s in subtopics])+'\n') 
            self.sub_blank = ''.join([' ' for x in range(len(msg))])
        self.jblank = self.sub_blank 
        self.log.append([topic])
        self.topic, self.topics, self.subtopic, self.events[topic] = topic, 0, None, dd(list) 
        self.update_indent, self.update_len = self.h_blank, 0  
        


    def start_minor(self,topic,block_len=10,ENUMERATE=False, REVEAL=[], SKIP=False, FIN=True):
        if self.status == 'minor' and FIN: self.end() 
        if len(REVEAL) > 0 and self.topics > 0:  self.reveal_new_settings(REVEAL[0]) 
        self.topics += 1 
        self.counter, self.status, self.subtopic = 0, 'minor', topic 
        self.log[-1].append(topic)
        #self.ENUMERATE=ENUMERATE
        self.block_len,self.mp,self.dots,self.numbers = block_len,block_len,0,0 
        if block_len > 100: self.mp = 100 
        if SKIP: self.write('\n') 
        self.write(self.sub_blank+'JOB'+str(self.topics)+': '+topic+'...')
        if not FIN: self.write('........................................\n') 
        self.jblank = self.sub_blank + ''.join([' ' for x in range(6)]) 
        self.update_indent,self.update_len = self.sub_blank,0  
     


    def mark(self,dots=0):
        self.counter +=1 
        if self.active: 
            if dots > 0: self.write('.'.join(['' for x in range(dots)]))
            elif dots == 0: 
                if self.counter % self.mp == 0: 
                    self.write('.') 
                    self.dots +=1
                    if self.dots > 100: self.mp = self.counter - 1 
                    elif self.dots > 80: self.mp = self.counter / 2 
                    elif self.dots > 60: self.mp = self.counter / 4 
                    elif self.dots > 50: self.mp = self.counter / 8
                    elif self.dots > 40: self.mp = self.counter / 16
                    elif self.dots > 30: self.mp = self.counter / 8
                    elif self.dots > 20: self.mp = self.counter / 4
                    elif self.dots > 10: self.mp = self.counter / 2  
            else:
                if self.counter % self.mp == 0: 
                    self.counter *= 5 
                    self.write('.') 
                    self.dots +=1
                    if self.dots >   20: self.mp = self.counter - 1
                    elif self.dots > 10: self.mp = self.counter /  2   
                    else:                self.mp = self.counter * 2 

                        
    

    def fail(self, MSG, notes = []): 
        self.active = True
        self.write('ERROR!\n') 
        self.write(self.jblank+MSG+'\n') 
        for n in notes:
            self.write(self.jblank+n+'\n') 
        sys.exit(2) 

                    

    def update99(self, UPDATE): 
        self.update_len += len(UPDATE) + 2 
        if self.update_len > 40: 
            self.out.write('\n'+self.update_indent+'....') 
            self.update_len = 0 
        self.out.write(UPDATE+'..') 


    def end(self,NOTE=None): 
        if NOTE != None: self.write(NOTE+'\n')
        elif self.active and self.subtopic != None: self.write('Complete\n') 
        try:                        self.events[self.topic][self.subtopic], self.subtopic = [self.counter], None 
        except AttributeError:      self.subtopic = None 
        except KeyError:            self.subtopic, self.update_indent,self.update_len = None, '',0  



    def finish(self,NOTE=None):
        self.end(NOTE) 
        self.write('\n') 
