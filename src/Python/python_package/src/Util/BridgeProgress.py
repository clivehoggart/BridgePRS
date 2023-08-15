import sys, os 
from collections import defaultdict as dd
import time 

LOCALTIME = time.asctime( time.localtime(time.time()) )


class BridgeProgress:
    #def __init__(self,VERBOSE,msg=' New_Module ',PARENT=False):
    #def __init__(self,args,COMMAND_LINE=None, OUT_FILE=None):
    def __init__(self,args, command_line): 


        self.obs_input = [] 
        self.INIT, self.active, self.loud = True, True, False 

        if args.verbose:  self.active, self.loud = True, True 
        elif args.silent: self.active, self.loud = False, False 
        else:             self.active, self.loud = True, False 

        self.args, self.command_line = args, command_line  
        self.out = sys.stderr
        
        
        
        
    def begin_module(self, module, cmd, runpath):
        self.module, self.cmd = module, cmd 

        if self.INIT: 
            self.logfile = runpath + '/bridgePRS.'+self.args.module+'.'+self.args.cmd+'.log' 
            self.loghandle = open(self.logfile,'w') 
            self.FILE, self.log = True, [[]] 
            self.write('BridgePRS Begins at '+LOCALTIME+' \n\n') 
            if self.command_line != None:  self.write('Bridge Command-Line: '+' '.join(self.command_line)+'\n\n')
            self.INIT = False  

        self.update_indent, self.update_len = '', 0 
        self.start_heading = 'Module: '
        self.job_heading   = 'Command: ' 
        self.write('\n'+self.start_heading+self.module+'\n') 
        self.h_blank = ''.join([' ' for x in range(len(self.start_heading))])
        self.sub_blank = ' '
        self.topic, self.subtopic, self.status, self.topics, self.events = None, None, None, [], {} 

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

    def show_settings(self,settings): 
        self.write(self.h_blank+'Settings:\n')  
        kFiles, kPrefixes = [[a,b] for a,b in settings.files.items() if b is not None], [[a,b] for a,b in settings.prefixes.items() if b is not None]
         

        if len(kFiles) > 0: 
            for i,(k,fn) in enumerate(sorted(kFiles)):
                self.obs_input.append(k) 
                if i == 0:  self.write(self.h_blank+'         Files:      '+k+'_file='+fn+'\n') 
                else:       self.write(self.h_blank+'                     '+k+'_file='+fn+'\n') 
            self.write('\n') 
        if len(kPrefixes) > 0: 
            for i,(k,fn) in enumerate(sorted(kPrefixes)): 
                self.obs_input.append(k) 
                if i == 0:  self.write(self.h_blank+'         Prefixes:   '+k+'_prefix='+fn+'\n') 
                else:       self.write(self.h_blank+'                     '+k+'_prefix='+fn+'\n') 
            self.write('\n') 
            
        kFlags, kOpts, kInputs, kObs = [] , [] , [], [] 
        for k in vars(settings.args):
            if k.split('_')[-1] in ['config','prefix','file','module','cmd','dataset']: continue 
            v = vars(settings.args)[k] 
            if v in [None,0]: continue
            kObs.append(k) 
            if str(v) in ['True','False']: kFlags.append([k,v]) 
            elif type(v) != list:          kOpts.append([k,v])              
            else:                          kInputs.append([k,v])  

        for k,v in settings.config.items(): 
            if k.split('_')[-1] in ['config','prefix','file','module','cmd','dataset']: continue 
            if k in kObs: continue
            if str(v) in ['True','False']: kFlags.append([k,v]) 
            else:                          kOpts.append([k,v])              
        
        if len(kOpts) > 0:
            kOpts.sort() 
            self.write(self.h_blank+'         Options:   '+',  '.join([k+'='+str(fn) for k,fn in kOpts[0:6]])+'\n')
            if len(kOpts) > 6: 
                for xp in range(6,len(kOpts),6):   self.write(self.h_blank+'                    '+',  '.join([k+'='+str(fn) for k,fn in kOpts[xp:xp+6]])+'\n')


            #self.write(self.h_blank+'         Options:   '+',  '.join([k+'='+str(fn) for k,fn in sorted(kOpts)])+'\n')
        if len(kFlags) > 0: 
            self.write(self.h_blank+'         Flags:     '+',  '.join([k+'='+str(fn) for k,fn in sorted(kFlags)])+'\n')
        
        if len(kInputs) > 0:
            for k,v in kInputs: 
                self.write(self.h_blank+'         '+k+':     '+',  '.join(v)+'\n') # [k+'='+str(fn) for k,fn in sorted(kFlags)])+'\n')
        self.write('\n') 

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


    def start_rJob(self, RJ): 

        if RJ[1] == '--vanilla': rInit,rJ = RJ[0:3], RJ[3::] 
        else:                    rInit,rJ = RJ[0:2], RJ[2::] 
        self.new_bk = self.sub_blank+'         ' 
        if self.status.strip() == 'minor':
            
            self.whisper('\n'+self.sub_blank+'Running Rscript:\n') 
            self.whisper(self.new_bk+" ".join(rInit)+'\n') 
            for i in range(1,len(rJ),2):
                #print(i, rJ[i-1], rJ[i]) 
                self.whisper(self.new_bk+rJ[i-1]+' '+rJ[i]+'\n') 
                
        
        self.whisper(self.new_bk+'.....................................') 
        return 




    def start_major(self,topic,msg=None,subtopics=[]): 
        
        if self.status == 'minor': self.end() 
        self.status = 'major'

        self.write(self.h_blank+self.job_heading+topic+'\n') 

        self.blank = self.h_blank+''.join([' ' for x in range(len(self.job_heading))])
        
        #print('yooo',msg)
        #sys.exit() 

        if msg: 
                self.write(self.blank+msg+' '+', '.join([str(s) for s in subtopics])+'\n') 
                self.sub_blank = ''.join([' ' for x in range(len(msg))])
        else:
                    self.sub_blank = self.blank    
        self.jblank = self.sub_blank 
        self.log.append([topic])
        self.topic, self.subtopic = topic, None 
        self.topics.append(topic) 
        self.events[self.topic] = dd(list) 
        self.topics = 0
        self.update_indent,self.update_len = self.h_blank,0  
        


    def start_minor(self,topic,block_len=10,ENUMERATE=False, REVEAL=[]):
        
        if self.status == 'minor': self.end() 

        if len(REVEAL) > 0 and self.topics > 0:  self.reveal_new_settings(REVEAL[0]) 
        self.topics += 1 
        self.status = 'minor'
        self.subtopic = topic
        self.log[-1].append(topic)
        self.counter = 0 
        self.ENUMERATE=ENUMERATE
        self.block_len,self.mp,self.dots,self.numbers = block_len,block_len,0,0 
        if block_len > 100: self.mp = 100 
        self.write(self.sub_blank+'JOB'+str(self.topics)+': '+topic+'...')
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

                    

    def update(self, UPDATE): 
        self.update_len += len(UPDATE) + 2 
        if self.update_len > 40: 
            self.out.write('\n'+self.update_indent+'....') 
            self.update_len = 0 

        self.out.write(UPDATE+'..') 


    def end(self,NOTE=None): 
        if NOTE != None: 
            self.write(NOTE+'\n')
        elif self.active and self.subtopic != None: self.write('Complete\n') 
        try:    
            self.events[self.topic][self.subtopic] = [self.counter]
            self.subtopic = None 
        except AttributeError:
            self.subtopic = None 
        except KeyError:
            self.subtopic = None
            self.update_indent,self.update_len = '',0  



    def finish(self,NOTE=None):
        self.end(NOTE) 
        self.write('\n') 
