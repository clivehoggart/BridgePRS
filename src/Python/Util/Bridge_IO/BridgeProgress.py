import sys, os 
from collections import defaultdict as dd
import multiprocessing, time 
LOCALTIME = time.asctime( time.localtime(time.time()) )


# JOB optimize  fail Flags 
# Generated yooooo printout 

class BridgeProgress:
    def __init__(self,args, command_line): 
        self.obs_input = [] 
        self.topics = 0  
        self.status = None  
        self.INIT, self.active, self.loud, self.HOMEDIR = True, True, False, True 
        self.sub_blank = '     '
        self.bk1, self.bk2 = ' ', ' ' 
        if   args.verbose:  self.active, self.loud = True, True 
        else:               self.active, self.loud = True, False 
        self.args, self.command_line, self.out, = args, command_line, sys.stderr 



    def initialize(self, runpath, poplist):
        self.runpath, self.prepath = runpath, "/".join(runpath.split('/')[0:-1])+'/'  
        try: self.homedir = os.path.expanduser('~') 
        except: self.HOMEDIR = False 
        self.last_line, self.run_len, self.dot_loc, self.line_loc = '', 100, 0, 0 
        if len(poplist) == 0:   self.logfile = runpath + '/logs/bridgePRS.'+self.args.module+'.'+self.args.cmd+'.log' 
        else:                   self.logfile = runpath + '/logs/bridgePRS.'+'-'.join([p.lower() for p in poplist])+'.'+self.args.module+'.'+self.args.cmd+'.log'  
        self.loghandle = open(self.logfile,'w') 
        self.FILE  = True 
        self.write('BridgePRS Begins at '+LOCALTIME+' \n') 
        if self.command_line != None:  self.write('Bridge Command-Line:  '+' '.join(self.command_line)+'\n')
        self.JOB_RANK = dd(int)  
        for jb,ji in [['clump',1],['beta',2],['predict',3],['quantify',4],['prior',4],['result',5]]: self.JOB_RANK[jb] = ji 
        return self 

    # SHOWING SHOWING SHOWING SHOWING SHOWING SHOWING SHOWING SHOWING SHOWING SHOWING SHOWING SHOWING SHOWING  # 

    def show_requirements(self, RK, LOC, R_DATA, plink_cmd):
        self.new_bk = '            ' 
        self.write('Checking Requirements:\n') 
        cA, cX = multiprocessing.cpu_count(), self.args.cores
        tC, mC, dC = str(self.args.total_cores), str(self.args.cores), self.args.total_cores - self.args.cores 
        fs1 = '%22s  %-15s  %-25s\n' 
        fs2 = '%22s  %-15s  %-25s  %15s\n' 
        fs3 = '%22s  %-15s  %-25s  %15s       %-65s\n' 
        if self.args.cores > 1 or dC < 2: self.say(fs2,('System:','platform='+self.args.platform+',', 'cores(available)='+tC+',', 'cores(used)='+mC)) 
        else:                             self.say(fs3,('System:','platform='+self.args.platform+',', 'cores(available)='+tC+',', 'cores(used)='+mC, '(TIP: Using More Than One Core Will Improve Performace (e.g. --cores '+str(dC)+'))')) 
        if RK['plink']: self.say(fs1, ('Plink:','found=true,','path='+LOC['plink'])) 
        else:           self.say(fs3, ('Plink:','found=false,','path=NA','', '(Using included version: '+LOC['plink']+')')) 
        os.system(plink_cmd)  
        if RK['R'] and len(R_DATA[1]) == 0: 
            self.say(fs3, ('R:','found=true,','path='+LOC['R']+',','version='+R_DATA[0], '(packages=up to date)')) 
        else: 
            if not RK['R']: 
                self.say(fs1, ('R:','found=false','')) 
                self.say('\n'+fs2+'\n', ('RequirementError:', 'R Not Found (please install from r-project.org)','','')) 
            else: 
                rc = 'install.packages(c('+",".join(['\''+x+'\'' for x in R_DATA[1]])+', repos =\"http://cran.us.r-project.org\"))' 
                self.say(fs3, ('R:','found=true,','path='+LOC['R']+',','version='+R_DATA[0],'(packages=incomplete)'))
                self.say('\n'+fs2, ('RequirementError:', 'Missing R Packages',",".join(R_DATA[1]),'')) 
                self.write('            To Install, Type The Command in R: '+rc+'\n\n')  
            sys.exit() 
        if RK['python3'] and RK['matplotlib']: self.say(fs2, ('Python3:','found=true,','path='+LOC['python3']+',','matplotlib=true')) 
        elif RK['python3']: self.say(fs3, ('Python3:','found=true,','path='+LOC['python3']+',','matplotlib=false','(NOTE: Please Install Matplotlib To Enable Plotting)'))
        else:               self.say('\n'+fs1+'\n', ('RequirementError:', 'Python3 Not Found,','Please Install Python3 To Use Wrapper'))
        return 


    def record_me(self, n,v, INIT=False, KEEP=True, NOTE=False): 
        if INIT: 
            if NOTE: self.say('%s\n', (n+'='+str(v)+'    '+NOTE)) 
            else:    self.say('%s\n', (n+'='+str(v))) 
        else: 
            if NOTE: self.say('%25s %s\n',(' ',n+'='+str(v)+'     '+NOTE)) 
            else:    self.say('%25s %s\n',(' ',n+'='+str(v))) 
        if KEEP: self.REC.write(n+'='+str(v)+'\n') 
        else:    self.SEC.write(n+'='+str(v)+'\n') 
        return


    def show_pop_data(self, pop_data): 
        fs0, fs1, fs2 = '%30s  %-40s\n', '%30s  %-75s  %-25s\n', '%30s  %-22s  %-50s  %15s\n' 
        self.write('\nReading Population Data:\n') 
        for i, p in enumerate(pop_data):
            if p is None: continue 
            my_notes = [False,False,False]
            
            bd, ss, pt = p.bdata, p.sumstats, p.genopheno 
            self.tFile = self.runpath + '/save/'+p.name+'.'+p.type+'.config'
            self.sFile = self.runpath + '/save/'+p.name+'.'+p.type+'.source'
            self.REC, self.SEC = open(self.tFile,'w'), open(self.sFile, 'w') 
            self.say('%25s\n',(p.type.capitalize()+' Data:')) 
            self.say('%25s ',('Names:')) 
            self.record_me('POP',p.name, INIT=True) 
            self.record_me('LDPOP',p.ref_pop) 
            self.record_me('LD_PATH',bd.ld_path) 
            if p.clump_value is not None: self.record_me('CLUMP_VALUE',p.clump_value) 
            if ss.TESTS['INFER_SUFFIX']: my_notes[0] = '(WARNING: Not Given - Inferred From Directory)' 
            if ss.TESTS['NEWSNPS']:      my_notes[1] = '(WARNING: Not Given - Created Using All Snps)'
            if pt.TESTS['NOVALID']:      my_notes[2] = '(WARNING: Not Given - Created by Splitting Phenotype File)'
            self.say('%25s ',('Sumstats:')) 
            self.record_me('SOURCE_PREFIX',ss.source_prefix,INIT=True,KEEP=False) 
            self.record_me('SOURCE_SUFFIX',ss.source_suffix,KEEP=False,NOTE=my_notes[0]) 
            self.record_me('SUMSTATS_PREFIX',ss.prefix) 
            self.record_me('SUMSTATS_SUFFIX',ss.suffix) 
            self.record_me('SUMSTATS_FIELDS',",".join([p.field_key[k] for k in ['ID','REF','ALT','P','BETA']]))
            self.record_me('SUMSTATS_SIZE',ss.pop_size) 
            self.say('%25s ',('QC-Snps:')) 
            self.record_me('SNP_FILE',ss.snp_file,NOTE=my_notes[1],INIT=True) 
            self.say('%25s ',('Genotype/Phenotype:')) 
            if not pt.VALID:
                self.record_me('GENOTYPE_PREFIX','NONE',INIT=True,KEEP=False) 
                self.record_me('PHENOTYPE_FILE','NONE',KEEP=False) 
                self.record_me('PHENOTYPE_VARIABLES','NONE',KEEP=False) 
            else: 
                self.record_me('GENOTYPE_PREFIX',pt.genotype_prefix,INIT=True)  
                self.record_me('PHENOTYPE_FILE',pt.files[0])  
                if i == 0: self.record_me('VALIDATION_FILE',pt.files[1],NOTE=my_notes[2])
                self.record_me('CHROMOSOMES',','.join(p.chromosomes),KEEP=False)  
                self.record_me('PHENOTYPE_VARIABLES',",".join(pt.header[2::]),KEEP=False) 
            
            self.record_me('DECLARED_PHENOTYPE',str(self.args.phenotype),KEEP=False) 
            if self.args.debug_level > 0: 
                pk, pn = ss.pop.type, ss.pop.name
                g1, g2 = [], [] 
                if self.args.phenotype is not None: 
                    self.record_me('PHENOTYPE_CLASS',str(pt.type),KEEP=False) 
                    eR = min(ss.WT),max(ss.WT)                                                                                                                                                                                    
                    eStr = 'Effect Size Range: '+str(eR[0])+', '+str(eR[1])                                                                                                                                                 
                    self.record_me('EFFECT_SIZE_RANGE',str(round(eR[0],3))+','+str(round(eR[1],3)), KEEP=False) 
                    if pt == 'binary' and eR[0] > 0:  self.fail('Truncated Effect Size Range Found (odds-ratio). Please use log(odds) for binary traits') 
            
                    for nk,n in [['MATCH','Matching'],['SWAPREF','Ref/Alt Swap'],['REVCOMP','Reverse Complement'],['INVALID','Invalid Bases']]:                                                                                 
                        if ss.CK[nk] > 0 :        g1.append(str(ss.CK[nk])+' ('+n+')')                                                                                                                                            
                        if ss.CK['GENO_'+nk] > 0: g2.append(str(ss.CK['GENO_'+nk])+' ('+n+')')                                                                                                                                    
                    if i == 0:                                                                                                                                                                                                 
                        g0 = str(min(ss.snp_list_len, ss.genome_snps))+' (Total),'
                        self.record_me('GWAS_Variants','['+g0+' '+", ".join(g1)+']', KEEP=False) 
                    else:                                                                                                                                                                                                       
                        self.record_me('Shared_Variants','['+", ".join(g1)+']', KEEP=False) 
            self.say('%25s\n\n','    **'+p.type.capitalize()+' Config Made:'+' '+self.tFile) 
            self.REC.close()
            self.SEC.close() 
            p.config = self.tFile  
            p.info   = self.sFile  
        return                            

    
    def show_settings(self):
        kSkip = ['pop','ldpop','files','outpath','platform','cores', 'config','prefix','file','path','module','cmd','dataset','ssf','pf','sumstats_suffix','phenotype_files'] 
        kParams, kTrue, kFalse = {}, [], [] 
        self.write('Setting Program Parameters:\n') 
        for k in vars(self.args): 
            kv = vars(self.args)[k]  
            if k.split('_')[-1] in kSkip or k.split('-')[0] in kSkip or k in kSkip or k[-4::] in kSkip: continue 
            if k in ['verbose','silent','restart','noPlots','debug'] and  kv:      kTrue.append(k.upper()) 
            elif k in ['verbose','silent','restart','noPlots','debug'] and not kv: kFalse.append(k.upper()) 
            elif k == 'max_clump_size' and int(kv) == 0:             kParams[k] = 'NO_LIMIT' 
            else:                                                    kParams[k] = str(kv) 
        TSTR = [] 
        if len(kTrue) > 0:  TSTR.append('ON='+','.join(kTrue))
        if len(kFalse) > 0: TSTR.append('OFF='+','.join(kFalse)) 
        self.write('         Flags: '+", ".join(TSTR)+'\n')  
        try: 
            MSTR =  ', '.join(['Phenotype='+kParams['phenotype'], 'Covariates='+kParams['covariates'], 'Max_clump_size='+kParams['max_clump_size']]) 
            self.write('         Model: '+MSTR+'\n')  
        except: pass 
        self.write('\n') 
        return 
   



    def show_new_data(self, io=None): 
        rData = [] 
        if io is None: return [] 
        elif io.pop is None:
            return [] 
            print(io.args) 
            print() 
            print('reveal_new_data in progress') 
            print('here') 
            print('yes')
            sys.exit() 
        else:
            for k,P in io.pop.gen.items(): 
                if P in self.obs_input: continue 
                self.obs_input.append(P) 
                if k in ['clump','beta','predict','quantify','prior']: rData.append([self.JOB_RANK[k], k+'_prefix',P.split(self.prepath)[-1]]) 
                elif k in ['model','result']:                          rData.append([0, k+'_file',P.split(self.prepath)[-1]]) 
                else: 
                    print()
                    print('cant reveal') 
                    print(k,P) 

        if len(rData) == 0: return rData 
        rData.sort() 
        rN, rD = [rd[1] for rd in rData], self.condense_paths([rd[2] for rd in rData]) 
        return [a+'='+b for a,b in zip(rN, rD)]


       




    # MODULE MODULE MODULE MODULE MODULE MODULE # MODULE MODULE MODULE MODULE MODULE MODULE # 
    # MODULE MODULE MODULE MODULE MODULE MODULE # MODULE MODULE MODULE MODULE MODULE MODULE # wtf 
    
    
    
    
    def start_module(self, module, cmd, runpath):
        self.bk1 = '  ' 
        self.start_heading = 'Begin Module: ' 
        if module in ['prs-single','stage1']:   self.write(self.bk1+self.start_heading+' Stage1: '+module+'\n') 
        elif module in ['prs-port','stage1.5']: self.write(self.bk1+self.start_heading+' Stage1.5: '+module+'\n') 
        elif module in ['prs-prior','stage2']:  self.write(self.bk1+self.start_heading+' Stage2: '+module+'\n') 
        else:                                       self.write(self.bk1+self.start_heading+' '+module+'\n') 
        self.topic, self.status, self.topics  = 'major', 'major', 0 
        return self
    
    def start_minor(self,topic,RD = None, block_len = 10):
        if self.status == 'minor': self.end(RD) 
        else: 
            rJ = self.show_new_data(RD)
            if len(rJ) > 0: 
                if self.status is None and self.topics == 0: self.write(self.bk2+'Previously Generated: '+", ".join(rJ)+'\n') 
                elif self.status == 'major':                 self.write(self.bk1+'Previously Generated: '+", ".join(rJ)+'\n') 
        self.status = 'minor' 
        self.topics += 1 
        self.write(self.sub_blank+'JOB'+str(self.topics)+': '+topic+'...')
        return self 

    
    def start_warning(self, m1, m2=None, WTYPE = 'BridgeWarning: '):
        ws = "".join([' ' for x in range(len(WTYPE)+1)]) 
        if m2 is not None: 
            self.write('\n'+WTYPE+' '+m1+'\n'+ws+m2+'...')                                                                                                                                                                                                                     
            self.line_loc = len(ws) + len(m2) 
        else:              
            self.write('\n'+WTYPE+' '+m1+'...')                                                                                                                                                                                                                     
            self.line_loc =  len(WTYPE) + len(m1)
        if self.line_loc > 100: self.line_loc = 50
        return 

    def end_warning(self,m1): 
        self.write(m1+'\n') 
        return



























    def condense_paths(self, X): 
        X = [self.homeshrink(x) for x in X] 
        if len(X) < 2: return X 
        cL, cList, xLocs = X[0].split('/') , [X[0]],  [x.split('/') for x in X[1::]] 
        for xL in xLocs: 
            for i in range(min(len(xL),len(cL))): 
                if xL[i] != cL[i]: break  
            cList.append("/"+"/".join(xL[i::]))
        return cList    
    
   

    def homeshrink(self, f_name): 
        try:    f_tmp = f_name.split(self.prepath)[-1] 
        except: f_tmp = f_name 
        if not self.HOMEDIR: return f_tmp 
        if len(f_tmp.split('/')) < 3: return f_tmp
        if self.homedir not in f_tmp: return f_tmp 
        try: f_bit = '~'+f_tmp.split(self.homedir)[-1] 
        except: f_bit = t_fmp 
        return f_bit 
    

    
    def end(self,RD=None): 
        rJ = self.show_new_data(RD)
        if self.status == 'minor':
            dl = max(self.run_len - self.line_loc, 1) 
            cl = '.'.join(['' for x in range(dl)])+'Complete' 
            if len(rJ) > 0: self.write(cl+'\n     [File(s) Generated: '+", ".join(rJ)+']\n') 
            else:           self.write(cl+'\n') 
            self.status = None  
        return 
        
    def finish(self,MSG=None, FIN=False):
        if MSG is not None:  self.write(MSG+'\n') 
        if self.status == 'minor': self.write('...Complete\n') 
        if FIN: sys.exit() 

    # ERRORS ERRORS ERRORS ERRORS ERRORS ERRORS ERRORS ERRORS ERRORS # 

    


    def warn(self, MSG, WTYPE = 'BridgeWarning:'): 
        try: 
            e1 = self.sub_blank 
            e2 = self.sub_blank + ' '.join(['' for x in range(len(WTYPE))])+'  ' 
        except: 
            e1 = '' 
            e2 = ' '.join(['' for x in range(len(WTYPE))])+'  ' 
        if type(MSG) not in [list,tuple]: self.write('\n'+e1+WTYPE+' '+MSG) 
        else: 
            self.write('\n'+e1+WTYPE+' '+MSG[0]+'\n')                                                                                                                                                                                                                     
            for es in MSG[1::]:   self.write(e2+es+'\n')                                                                                                                                                                                                                        





    def fail(self, MSG, ETYPE = 'BridgeError:'): 
        EB = ' '.join(['' for x in range(len(ETYPE))])+'  '
        if type(MSG) not in [list,tuple]: self.write('\n'+ETYPE+' '+MSG) 
        else: 
            self.write('\n'+ETYPE+' '+MSG[0]+'\n')                                                                                                                                                                                                                     
            for es in MSG[1::]:   self.write(EB+es+'\n')                                                                                                                                                                                                                        
        self.write('\n') 
        sys.exit(2) 



    def broke(self, MSG, extra = [], WTYPE='BridgePipelineError'): 
        try: 
            e1 = self.sub_blank 
            e2 = self.sub_blank + ' '.join(['' for x in range(len(WTYPE))])+'  ' 
        except: 
            e1 = '' 
            e2 = ' '.join(['' for x in range(len(WTYPE))])+'  ' 
        if type(MSG) not in [list,tuple]: self.write('\n'+e1+WTYPE+' '+MSG) 
        else: 
            self.write('\n'+e1+WTYPE+' '+MSG[0]+'\n')                                                                                                                                                                                                                     
            for es in MSG[1::]:   self.write(e2+es+'\n')                                                                                                                                                                                                                        
        sys.exit()         





















    
    # WRITING WRITING WRITING WRITING WRITING # WRITING WRITING WRITING WRITING WRITING # 
    def write(self, outstring, DOTS = False): 
        if self.active: 
            self.out.write(outstring) 
            self.out.flush()  
        if self.FILE: self.loghandle.write(outstring) 
        if DOTS: self.dot_loc += len(outstring) 
        else:    self.dot_loc = 0 
        if '\n' in outstring: self.line_loc = len(outstring.split('\n')[-1]) 
        else:                 self.line_loc += len(outstring) 
        return 

    def whisper(self, outstring): 
        if self.loud: 
            self.out.write(outstring) 
            self.out.flush() 
        if self.FILE: self.loghandle.write(outstring) 
        return 

    def say(self, outformat, outtuple): 
        if self.active: 
            self.out.write(outformat % outtuple) 
            self.out.flush()  
        if self.FILE: self.loghandle.write(outformat % outtuple)

    
    

    def mark(self,dots=1):
        if self.line_loc >   110:    return 
        elif self.line_loc >  80:    dl = 1
        elif self.line_loc >  65:    dl = dots 
        elif self.line_loc >  50:    dl = 2*dots 
        else:                        dl = 3*dots 
        mark_string = '.'.join(['' for x in range(dl+1)])  
        
        self.write(mark_string, DOTS = True)  





    # RCODE RCODE RCODE RCODE # RCODE RCODE RCODE RCODE # RCODE RCODE RCODE RCODE 
    # RCODE RCODE RCODE RCODE # RCODE RCODE RCODE RCODE # RCODE RCODE RCODE RCODE 
    def start_rJob(self, RJ, job_name): 
        rPaths, rCols, rSave, rOne, rVal = [], [], [], [], [] 
        self.new_bk = self.sub_blank+'           ' 
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
        dl = max(self.run_len - self.line_loc, 5) 
        cl = '.'.join(['' for x in range(dl)])
        if job_name == 'clump': 
            self.whisper(cl+'\n'+self.sub_blank+'Running Plink:\n') 
            K = ['bfile', 'clump','extract', 'keep', 'out']
            self.show_local(K, LP) 
            if len(rSave) > 0: 
                self.whisper(self.new_bk+rSave[0][0]+' '+rSave[0][1]) 
                for rs in rSave[1::]: self.whisper(' '+rs[0]+' '+rs[1]) 
                self.whisper('\n') 
        else:    
            if RJ[1] == '--vanilla': rInit,rJ = RJ[0:3], RJ[3::] 
            else:                    rInit,rJ = RJ[0:2], RJ[2::] 
            self.whisper(cl+'\n'+self.sub_blank+'Running Rscript:\n') 
            self.whisper(self.new_bk+" ".join(rInit)+'\n') 
            self.show_local(K, LP) 
            if len(rOne) > 0: self.whisper(self.new_bk+" ".join(rOne)+'\n') 
            if len(rSave) > 0: 
                self.whisper(self.new_bk+rSave[0][0]+' '+rSave[0][1]) 
                for rs in rSave[1::]: self.whisper(' '+rs[0]+' '+rs[1]) 
                self.whisper('\n') 
        nl = self.new_bk+'.......................................'
        self.whisper(nl) 
        if self.args.verbose: self.line_loc, self.last_line = len(nl), nl 
        return  

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
            if k not in LP: continue 
            a,b = LP[k]

            if b[-1] == '1': b = b.strip('1')+'*' 
            if i == 0:   self.whisper(self.new_bk+a+' '+b+'\n') 
            elif i == 1:     
                self.whisper(self.new_bk+a+' '+b) 
                my_len = len(a)+len(b) 
            else: 
                if my_len < 70: 
                    self.whisper('  '+a+' '+b) 
                    my_len += len(a)+len(b) 
                else:           
                    self.whisper('\n'+self.new_bk+a+' '+b) 
                    my_len = len(a)+len(b) 
        self.whisper('\n') 
        return  

