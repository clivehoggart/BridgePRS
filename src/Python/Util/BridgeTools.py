import sys, os 
import gzip
from collections import defaultdict as dd


def bridge_error(eString):
    if type(eString) in [list,tuple]:
        sys.stderr.write('\n     BridgeToolBoxError: '+eString[0]+'\n')                     
        for es in eString[1::]: sys.stderr.write('                         '+es+'\n')
    else: sys.stderr.write('\n    BridgeToolBoxError: '+eString+'\n')                           
    sys.exit()

                                                                                                                
def bridge_warning(eString):                                                                           
    if type(eString) in [list,tuple]:                                                                           
        sys.stderr.write('\n     BridgeToolBoxWarning: '+eString[0]+'\n')                                           
        for es in eString[1::]: sys.stderr.write('                           '+es+'\n')                                   
    else: sys.stderr.write('\n     BridgeToolBoxWarning: '+eString+'\n')                                            


def zip_open(fp, HEADER = False):                                                                               
    if fp.split('.')[-1] == 'gz': gf = gzip.open(fp, 'rt')                                                      
    else:                         gf = open(fp, 'rt')                                                           
    if HEADER:                                                                                                  
        HL = gf.readline().split()                                                                              
        gf.close()                                                                                              
        return HL                                                                                               
    else: return gf  



class BridgeTools:
    def __init__(self, io):
        self.args, self.io = io.args, io 
                
        self.GT = False 
        self.io.progress.start_module(self.io.module, self.io.cmd, self.io.paths['home'])# .show_settings(self.settings) 
    
        if self.io.cmd == 'reformat_sumstats': self.apply = self.reformat_sumstats 


    def read_genotype_prefix(self, gt_prefix): 
        g_path, g_name = "/".join(gt_prefix.split('/')[0:-1]), gt_prefix.split('/')[-1] 
        self.gt_lookup = dd(bool) 
        for f in os.listdir(g_path): 
            if f.split('.')[-1] == 'bim' and f[0:len(g_name)] == g_name: 
                with open(g_path+'/'+f) as filehandle:
                    for line in filehandle: self.gt_lookup[line.split()[1]] = line.split()[0]  
        self.GT = True 


    def reformat_sumstats(self): 
        if len(self.args.sumstats_prefix) != 1:    bridge_error('A Single Sumstats Prefix Is Required')
        if len(self.args.genotype_prefix) != 1:    bridge_error('A Single Genotype Prefix Is Required') 
        self.read_genotype_prefix(self.args.genotype_prefix[0]) 

        outprefix = self.io.paths['home']+'/'+".".join(self.args.sumstats_prefix[0].split('/')[-1].split('.')[0:-1])
        self.io.progress.start_minor('Verifying Sumstats File') 
        s_args = ['chr','snpid','ref','alt','p','beta','se','maf','n'] 
        s_names = ['#CHR','ID','REF','A1','P','BETA','SE','A1_FREQ','OBS_CT'] 
        s_placeholders  = {'n': '1000', 'maf': '0.1'}
        s_found, s_available = {}, [] 
        s_heads = {vars(self.args)[v]: v.split('-')[-1] for v in vars(self.args) if v.split('-')[0] == 'ssf'}



        s_handle = zip_open(self.args.sumstats_prefix[0]) 
        s_header = s_handle.readline().split()        
        for i,k in enumerate(s_header): 
            if 'CHR' in k.upper(): s_found['chr'] = i 
            elif k in s_heads:     s_found[s_heads[k]] = i 
            else:                  s_available.append(k) 

        s_missing = [s for s in s_args[1::] if s not in s_found]
        for s in s_missing: 
            if s not in s_placeholders:  bridge_error(['Missing Header Field (Required): '+s,'Assign using --ssf-'+s,'Available Values: '+", ".join(s_available)])  
            else:                        bridge_warning(['Missing Header Field: '+s,'Assign using --ssf-'+s+' (Available Values: '+", ".join(s_available)+')','Otherwise placeholder will be assigned: '+str(s_placeholders[s])]) 
        
        my_locs, my_defs   = [], []  
        for s in s_args[1::]: 
            if s in s_found: my_locs.append(s_found[s]) 
            else:            my_defs.append(s_placeholders[s]) 

        if len(s_missing) > 0: self.io.progress.write('     ') 
        self.io.progress.start_minor('Writing Formatted Sumstats File') 
        w = open(outprefix+'.reformat','w') 
        my_format = '%-6s %16s %5s %5s %16s %22s %16s %16s %16s\n'
        w.write(my_format % tuple(s_names)) 
        snp_loc = s_found['snpid'] 
        for total_snps,line in enumerate(s_handle): 
            line = line.split()
            s_id = line[snp_loc] 
            if s_id not in self.gt_lookup: continue 
            if s_id[0:2] != 'rs':          continue 
            ld = tuple([self.gt_lookup[s_id]]+[line[loc] for loc in my_locs]+my_defs)
            w.write(my_format % ld) 
            #if total_snps > 10: break 



        w.close()  
        #if 'chr' not in s_found: 
        #    self.io.progress.start_minor('Sorting Initial Sumstats File') 
        #    os.system('sort -k1,1 '+outprefix+'.tmp > '+outprefix+'.srt') 
        #    self.io.progress.start_minor('Merging Sumstats with Chromosome Key') 
        #    os.system('join '+self.args.snpdb+' '+outprefix+'.srt > '+outprefix+'.reformat') 
        #    CK = {} 
        #    with open(self.args.snpdb) as d_handle: 
        #        for line in d_handle: 
        #            break
        #            line = line.split() 
        #            CK[line[0]] = line[1] 

        self.io.progress.start_minor('Zipping Sumstats File') 
        os.system('gzip -rf '+outprefix+'.reformat')

        
        self.io.progress.end() 
        self.io.progress.write('     Reformatted Sumstats File Created Successfully: '+outprefix+'.reformat.gz\n')             
        sys.exit() 







        












def bridge_offer_summary(): 
    print(r""" 
    
    
    This Manual demonstrates how to bridgePRS on the included test data, 
        hit any key, or choose a specfic chapter (1-5), or (q) to quit: 
      
      
      Manual Chapters: 1: Run prs step-by-step
                         2: Run prs using a multistep pipeline 
                         3: Build a prior model (typically european data) 
                         4: PRS Portability (prs-port) 
                         5: Bridge PRS      (prs-bridge) 

------------------------------------------------------------------------------------ 
            """) 


def bridge_manual_image(): 
    print(r"""                        WELCOME TO THE BRIDGE MANUAL
                                                                        
         ^^       ..                                       ..       ^^
                  []                                       []   ^^
      ^^        .:[]:_          ^^             ^^        ,:[]:.     ^^
              .: :[]: :-.                             ,-: :[]: :.
            .: : :[]: : :`._                       ,.': : :[]: : :.
          .: : : :[]: : : : :-._               _,-: : : : :[]: : : :.
      _..: : : : :[]: : : : : : :-._________.-: : : : : : :[]: : : : :-._
      _:_:_:_:_:_:[]:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:[]:_:_:_:_:_:_
      !!!!!!!!!!!![]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!![]!!!!!!!!!!!!!
      ^^^^^^^^^^^^[]^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^[]^^^^^^^^^^^^^
                  []                                       []
                  []                                       []
                  []                                       []
       ~~^-~^_~^~/  \~^-~^~_~^-~_^~-^~_^~~-^~_~^~-~_~-^~_^/  \~^-~_~^-~~-
      ~ _~~- ~^-^~-^~~- ^~_^-^~~_ -~^_ -~_-~~^- _~~_~-^_ ~^-^~~-_^-~ ~^
        ~ ^- _~~_-  ~~ _ ~  ^~  - ~~^ _ -  ^~-  ~ _  ~~^  - ~_   - ~^_~
          ~-  ^_  ~^ -  ^~ _ - ~^~ _   _~^~-  _ ~~^ - _ ~ - _ ~~^ -
             ~^ -_ ~^^ -_ ~ _ - _ ~^~-  _~ -_   ~- _ ~^ _ -  ~ ^-
                 ~^~ - _ ^ - ~~~ _ - _ ~-^ ~ __- ~_ - ~  ~^_-
                     ~ ~- ^~ -  ~^ -  ~ ^~ - ~~  ^~ - ~
    """)


def bridge_splash_image(): 
    print(r"""                            WELCOME TO
   ^^      ..               BRIDGE PRS                        ..
            []                                       []
          .:[]:_          ^^                       ,:[]:.
        .: :[]: :-.                             ,-: :[]: :.
      .: : :[]: : :`._                       ,.': : :[]: : :.
    .: : : :[]: : : : :-._               _,-: : : : :[]: : : :.
_..: : : : :[]: : : : : : :-._________.-: : : : : : :[]: : : : :-._
_:_:_:_:_:_:[]:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:[]:_:_:_:_:_:_
!!!!!!!!!!!![]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!![]!!!!!!!!!!!!!
^^^^^^^^^^^^[]^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^[]^^^^^^^^^^^^^
            []                                       []
            []                                       []
            []                                       []
 ~~^-~^_~^~/  \~^-~^~_~^-~_^~-^~_^~~-^~_~^~-~_~-^~_^/  \~^-~_~^-~~-
~ _~~- ~^-^~-^~~- ^~_^-^~~_ -~^_ -~_-~~^- _~~_~-^_ ~^-^~~-_^-~ ~^
   ~ ^- _~~_-  ~~ _ ~  ^~  - ~~^ _ -  ^~-  ~ _  ~~^  - ~_   - ~^_~
     ~-  ^_  ~^ -  ^~ _ - ~^~ _   _~^~-  _ ~~^ - _ ~ - _ ~~^ -
        ~^ -_ ~^^ -_ ~ _ - _ ~^~-  _~ -_   ~- _ ~^ _ -  ~ ^-
            ~^~ - _ ^ - ~~~ _ - _ ~-^ ~ __- ~_ - ~  ~^_-
                ~ ~- ^~ -  ~^ -  ~ ^~ - ~~  ^~ - ~
    """)


def bridge_tutorial_image(): 
    print(r"""                                 BRIDGE TUTORIAL  
                                                                        
         ^^       ..                                       ..       ^^
                  []                                       []   ^^
      ^^        .:[]:_          ^^             ^^        ,:[]:.     ^^
              .: :[]: :-.                             ,-: :[]: :.
            .: : :[]: : :`._                       ,.': : :[]: : :.
          .: : : :[]: : : : :-._               _,-: : : : :[]: : : :.
      _..: : : : :[]: : : : : : :-._________.-: : : : : : :[]: : : : :-._
      _:_:_:_:_:_:[]:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:[]:_:_:_:_:_:_
      !!!!!!!!!!!![]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!![]!!!!!!!!!!!!!
      ^^^^^^^^^^^^[]^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^[]^^^^^^^^^^^^^
                  []                                       []
                  []                                       []
                  []                                       []
       ~~^-~^_~^~/  \~^-~^~_~^-~_^~-^~_^~~-^~_~^~-~_~-^~_^/  \~^-~_~^-~~-
      ~ _~~- ~^-^~-^~~- ^~_^-^~~_ -~^_ -~_-~~^- _~~_~-^_ ~^-^~~-_^-~ ~^
        ~ ^- _~~_-  ~~ _ ~  ^~  - ~~^ _ -  ^~-  ~ _  ~~^  - ~_   - ~^_~
          ~-  ^_  ~^ -  ^~ _ - ~^~ _   _~^~-  _ ~~^ - _ ~ - _ ~~^ -
             ~^ -_ ~^^ -_ ~ _ - _ ~^~-  _~ -_   ~- _ ~^ _ -  ~ ^-
                 ~^~ - _ ^ - ~~~ _ - _ ~-^ ~ __- ~_ - ~  ~^_-
                     ~ ~- ^~ -  ~^ -  ~ ^~ - ~~  ^~ - ~
    """)













def bridge_start_chapter(): 
    print(r"""                        
    
                """)                                                     


def bridge_make_chapter(): 
    print(r"""                        
                                                                     
-------------------------------------------------------------------------------------    
                                                                                   
                                                                                    
                                                                    """)


















class BridgeManual:
    def __init__(self, dc = None):
        self.dc = dc
        
    
    def say(self, msg): 
        sys.stderr.write(msg+'\n') 

    def begin(self): 
        bridge_manual_image() 
        bridge_offer_summary() 
        chapter = input('Make Selection:') 
        if chapter in ['1','2','3','4','5']: chapter = int(chapter) 
        elif chapter not in ['q','quit','Q']: chapter = 1
        else: self.quit() 
        self.run_chapter(chapter) 
        
    
    def quit(self): 
        self.say('Quitting Tutorial')
        sys.exit() 


    def run_chapter(self,chapter):  
        if chapter == 1: self.chapter_1() 
        elif chapter == 2: self.chapter_2() 
        elif chapter == 3: self.chapter_3() 
        elif chapter == 4: self.chapter_4() 
        elif chapter == 5: self.chapter_5() 
        else:              self.finish() 
    
    
    def complete_chapter(self,chap): 
        self.say('---------------------------------------------------------------------------------------') 
        self.say('Chapter Complete, Select Chapter (1-5), any key for next chapter, or q to quit:') 
        chapter = input('') 
        if chapter in ['1','2','3','4','5']: chapter = int(chapter) 
        elif chapter not in ['q','quit','Q']: chapter = chap+1 
        else:                                self.quit() 
        self.run_chapter(chapter) 
        return 

    def chapter_1(self): 
        bridge_start_chapter() 
        self.say('*********************************************************************************************************************')
        self.say('**********************************   Chapter 1: Running PRS step-by-step  *******************************************')
        self.say('*********************************************************************************************************************\n')
        self.say('The prs module includes 5 steps: clump, eval, predict, quantify\n(hit any key to scroll through the steps)\n') 
        self.say('') 
        foo = input()
        self.say('Clump Step:')  
        self.say('./bridgePRS prs clump -o out --p afr --bfile_prefix test_data/african/bfiles/bfile3 --sumstat_prefix test_data/african/sumstats/AFR --snp_file test_data/african/snps.txt --id_file test_data/african/ids.txt') 
        self.say('\n1) This command generates clump data in out/prs_AFR/clump') 
        self.say('2) To avoid many command line arguments, a config file can be passed instead:\n./bridgePRS prs clump -o out --popname afr --config test_data/african/afr.config') 
        self.say('--------------------------------------------------------------------------------------------------------------------------------------------------------------') 
        foo = input()
        self.say('\nEval Step:')  
        self.say('./bridgePRS prs eval -o out --p afr --config test_data/african/afr.config --clump_prefix out/prs_AFR/clump/AFR_clump') 
        self.say('\n1) Uses the clump data produced previously') 
        self.say('2) Generates evaluation data in out/prs_AFR/eval') 
        self.say('--------------------------------------------------------------------------------------------------------------------------------------------------------------') 
        foo = input()
        self.say('\nPredict Step:')  
        self.say('./bridgePRS prs predict -o out --p afr --config test_data/african/afr.config --eval_prefix out/prs_AFR/eval/AFR_eval') 
        self.say('\n1) Uses the eval data produced previously') 
        self.say('2) Generates prediction data in out/prs_AFR/predict') 
        self.say('--------------------------------------------------------------------------------------------------------------------------------------------------------------') 
        foo = input()
        self.say('\nQuantify Step:')  
        self.say('./bridgePRS prs quantify -o out --p afr --config test_data/african/afr.config --predict_prefix out/prs_AFR/predic/AFR_predict') 
        self.say('\n1) Uses the prediction data produced previously') 
        self.say('2) Generates quantification data in out/prs_AFR/quantify') 
        self.say('--------------------------------------------------------------------------------------------------------------------------------------------------------------') 
        foo = input()
        self.complete_chapter(1)  
        return 


    def chapter_2(self): 
        bridge_start_chapter() 
        self.say('*********************************************************************************************************************')
        self.say('**********************************  Chapter 2: Running PRS w one command  *******************************************')
        self.say('*********************************************************************************************************************\n')
        self.say('The prs module ran described in Chapter1 can be ran using a single command:') 
        self.say('\n./bridgePRS prs run -o out --p afr --config test_data/african/afr.config') 
        self.say('\n\n1) The pipeline will skip already completed steps') 
        self.say('') 
        self.complete_chapter(2)  
        foo = input()

    def chapter_3(self): 
        bridge_start_chapter() 
        self.say('*********************************************************************************************************************')
        self.say('**********************************  Chapter 3: Building Population Model ********************************************')
        self.say('*********************************************************************************************************************\n')
        
        self.say('Building A high coverage population model (usually european) can be built to improve results') 
        self.say('\n./bridgePRS build-model run -o out --p eur --config test_data/african/eur.config') 
        self.say('\n\n1) The pipeline will skip already completed steps') 
        foo = input()
        self.complete_chapter(3)  
    
    def chapter_4(self): 
        bridge_start_chapter() 
        self.say('*********************************************************************************************************************')
        self.say('**********************************  Chapter 4: Running PRS-PORT     *************************************************')
        self.say('*********************************************************************************************************************\n')
        
        self.say('PRS values from the European Model can be ported to the African Data') 
        self.say('\n./bridgePRS prs-port run -o out --p afr --config test_data/african/afr.config --model_file out/model_EUR/bridge.pipeline.result') 
        self.say('\n\n1) The low scores demonstrate the "portability problem.')  
        self.say('') 
        foo = input()
        self.complete_chapter(4)  

    def chapter_5(self): 
        bridge_start_chapter() 
        self.say('*********************************************************************************************************************')
        self.say('**********************************  Chapter 5: Running PRS-BRIDGE  **************************************************')
        self.say('*********************************************************************************************************************\n')
        self.say('The European Model Can be Used to Generate a Prior SNP Distribution') 
        self.say('./bridgePRS prs-bridge run -o out --p afr --config test_data/african/afr.config --model_file out/model_EUR/bridge.pipeline.result') 
        self.say('') 
        foo = input()
        self.complete_chapter(5)  


    def finish(self): 
        self.say('Tutorial Finished')
        sys.exit() 




















def bridge_start_chapter(): 
    print(r"""                                     """)                                                     


def bridge_make_chapter(): 
    print(r"""                        
                                                                     
-------------------------------------------------------------------------------------    
                                                                                   
                                                                                    
                                                                    """)








class BridgeTutorial:
    def __init__(self, dc = None):
        self.dc = dc
        
    
    def say(self, msg): 
        sys.stderr.write(msg+'\n') 

    def begin(self): 
        bridge_tutorial_image() 
        #bridge_manual_summary() 
        self.list_files() 
        self.list_prefixes() 
        self.list_options() 
        
    
    def quit(self): 
        self.say('Quitting Tutorial')
        sys.exit() 


    def run_chapter(self,chapter):  
        if chapter == 1: self.chapter_1() 
        elif chapter == 2: self.chapter_2() 
        elif chapter == 3: self.chapter_3() 
        elif chapter == 4: self.chapter_4() 
        elif chapter == 5: self.chapter_5() 
        else:              self.finish() 
    
    
    def complete_chapter(self,chap): 
        self.say('---------------------------------------------------------------------------------------') 
        self.say('Chapter Complete, Select Chapter (1-5), any key for next chapter, or q to quit:') 
        chapter = input('') 
        if chapter in ['1','2','3','4','5']: chapter = int(chapter) 
        elif chapter not in ['q','quit','Q']: chapter = chap+1 
        else:                                self.quit() 
        self.run_chapter(chapter) 
        return 

    def list_files(self): 
        bridge_start_chapter() 
        self.say('*****************************   QUICKSTART COMMANDS     *******************************')
        self.say('***************************************************************************************\n')
        self.say('From the directory "tests", type the following five command:\n(listed in tutorialCommands.txt)\n')  
        self.say('../bridgePRS prs run          -o out -p afr --config ../test_data/afr.pop.config') 
        self.say('../bridgePRS build-model run  -o out -p eur --config ../test_data/eur.pop.config')     # build a large european model to improve prs prediction in the african population
        self.say('../bridgePRS prs-port run     -o out -p afr --config ../test_data/afr.pop.config') #--model_file out/model_EUR/bridge.build-model.result  # run prs porting snp weights from europeans
        self.say('../bridgePRS prs-bridge run   -o out -p afr --config ../test_data/afr.pop.config') #--model_file out/model_EUR/bridge.build-model.result  # run prs using prior snps distributions dervied from europeans
        self.say('../bridgePRS analyze result -o out -p afr --results out/prs_AFR/bridge.prs.result out/prs-port_AFR/bridge.prs-port.result out/prs-bridge_AFR/bridge.prs-bridge.result') 
        return 
    
    def list_prefixes(self): 
        bridge_start_chapter() 
        self.say('*******************************   DESCRIPTION     ***************************************')
        self.say('*****************************************************************************************\n')

        self.say('../bridgePRS prs run         ->  runs prs using a single small african populaton in directory: out/prs_AFR') 
        self.say('../bridgePRS build-model run ->  builds a large european model to improve prs prediction in the african population (out/model_EUR)') 
        self.say('../bridgePRS prs-port run    ->  uses snp weights to port the european result to the african populaton (out/prs-port_AFR') 
        self.say('../bridgePRS prs-bridge run  ->  uses snp priors to bridge the european result to the african population') 
        return    



        self.say('') 
        return 

    def list_options(self): 
        bridge_start_chapter() 
        self.say('*******************************     OUTPUT     ****************************************')
        self.say('****************************************************************************************\n')
        self.say('../bridgePRS analyze result   ->  analyzes the three prs methods (single-pop, port, bridge)') 
        self.say('                              ->  a plot is produced in: /out/plot_AFR_analysis.png')
        self.say('') 
        self.say('Notice the different results from different prs methods') 
        return 



    def finish(self): 
        self.say('Tutorial Finished')
        sys.exit() 



def bridge_requirement_image(): 
    print(r"""                                 BRIDGE REQUIREMENTS  
                                                                        
         ^^       ..                                       ..       ^^
                  []                                       []   ^^
      ^^        .:[]:_          ^^             ^^        ,:[]:.     ^^
              .: :[]: :-.                             ,-: :[]: :.
            .: : :[]: : :`._                       ,.': : :[]: : :.
          .: : : :[]: : : : :-._               _,-: : : : :[]: : : :.
      _..: : : : :[]: : : : : : :-._________.-: : : : : : :[]: : : : :-._
      _:_:_:_:_:_:[]:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:[]:_:_:_:_:_:_
      !!!!!!!!!!!![]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!![]!!!!!!!!!!!!!
      ^^^^^^^^^^^^[]^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^[]^^^^^^^^^^^^^
                  []                                       []
                  []                                       []
                  []                                       []
       ~~^-~^_~^~/  \~^-~^~_~^-~_^~-^~_^~~-^~_~^~-~_~-^~_^/  \~^-~_~^-~~-
      ~ _~~- ~^-^~-^~~- ^~_^-^~~_ -~^_ -~_-~~^- _~~_~-^_ ~^-^~~-_^-~ ~^
        ~ ^- _~~_-  ~~ _ ~  ^~  - ~~^ _ -  ^~-  ~ _  ~~^  - ~_   - ~^_~
          ~-  ^_  ~^ -  ^~ _ - ~^~ _   _~^~-  _ ~~^ - _ ~ - _ ~~^ -
             ~^ -_ ~^^ -_ ~ _ - _ ~^~-  _~ -_   ~- _ ~^ _ -  ~ ^-
                 ~^~ - _ ^ - ~~~ _ - _ ~-^ ~ __- ~_ - ~  ~^_-
                     ~ ~- ^~ -  ~^ -  ~ ^~ - ~~  ^~ - ~
    """)













class BridgeRequirements:
    def __init__(self, dc = None):
        self.dc = dc
        
    
    def say(self, msg): 
        sys.stderr.write(msg+'\n') 

    def begin(self): 
        bridge_requirement_image() 
        #bridge_manual_summary() 
        self.list_basic() 
        self.list_packages() 
        self.list_options() 
        
    
    def quit(self): 
        self.say('Quitting Tutorial')
        sys.exit() 


    

    def list_basic(self): 
        bridge_start_chapter() 
        self.say('***************************************************************************************\n')
        self.say('bridgePRS requires python3, plink, and R.\n') 
        self.say('Plink is included in the src directory,\nPython can be downloaded from python.org,\nR can be downloaded from www.r-project.org\n') 
        return 
    
    def list_packages(self): 
        bridge_start_chapter() 
        self.say('*******************************   R Packages     ***************************************')
        self.say('*****************************************************************************************\n')
        self.say('bridgePRS requires that the following R-packages be installed:\n')
        self.say('BEDMatrix, R.utils, boot, data.table, doMC, glmnet, MASS, optparse, and parallel') 
        self.say('') 
        self.say('The command install.packages() can be used to install them') 
        #self.say('> install.packages(c('BEDMatrix', R.util, boot, data.table, doMC, glmnet, MASS, optparse, and parallel') 
        return    




    def list_options(self): 
        bridge_start_chapter() 
        self.say('*******************************   Python Libraries  ************************************')
        self.say('*****************************************************************************************\n')
        self.say('For plotting (optional) bridgePRS relies on matplotlib\n') 
        return 



    def finish(self): 
        self.say('Tutorial Finished')
        sys.exit() 




























#!/usr/bin/env python3






















HELP_FORMAT="%-50s %80s"
# Description Steps

def create_help_str(module):
    if module == 'prs':           return (HELP_FORMAT % ('Single Population PRS (Ridge Regression)','./bridgePRS prs run -o out --popname africa --config afro_pop.config')) 
    if module == 'build':         return (HELP_FORMAT % ('Build large pop model to allow (port/bridge)','./bridgePRS build-model run -o out --popname euro --config euro_pop.config'))
    if module == 'prs-port':      return (HELP_FORMAT % ('Port PRS snp-weights from large to small pop','./bridgePRS prs-port -o out --popname afr --model_file euro.model')) 
    if module == 'prs-bridge':    return (HELP_FORMAT % ('Bridge snp-priors from large to small pop','./bridgePRS prs-bridge -o out --popname afr --model_file euro.model'))
    if module == 'analyze':       return (HELP_FORMAT % ('Analyze One or More Results','./bridgePRS analyze single -o out --results out/afr.result')) 
    if module == 'easyrun':  return (HELP_FORMAT  % ('run','Run All Five Modules Consecutively (for advanced users)')) 





SUB_HELP_FORMAT="%-50s %100s"
SUB_HELP_FORMAT="%-80s %10s"
def create_sub_help_str(cmd, module=None):
    HS = ('INFO','INFO') 
    if cmd == 'run': HS = ('Run The Following Module Commands:','') 
    if cmd == 'clump': HS = ('Load Sumstats Data, Produce SNP Clumps','') 
    if cmd == 'eval': HS =  ('Evaluate SNP Clumps, Produce SNP Weights','') 
    if cmd == 'predict': HS =  ('Apply SNP Weights, Produce Polygenic Predictions','') 
    if cmd == 'quantify': HS =  ('Test Predictions, Produce Quantifiable Results','') 
    if cmd == 'prior': HS = ('Test Predictions,Produce SNP Prior Proability Matrix','') 
    if cmd == 'result': HS = ('Analyze PRS Results','') 
    if cmd == 'combine': HS = ('Combine And Analyze PRS Results','') 
    return(SUB_HELP_FORMAT % HS) 


def evaluate_detail(X,lp):
    if len(lp) == 1: return 'SPLASH','NA'
    if X in ['too few arguments','the following arguments are required: cmd']: return 'NO_COMMAND',lp[1] 
    if lp[1].upper() in ['REQUIREMENTS','MANUAL','TUTORIAL']: return lp[1].upper(), 'FULL' 
    return 'NA','NA'  


########################################################################################################################################################################

if __name__ == '__main__':

    import sys, os, argparse, platform, multiprocessing             
    from   src.Util.BridgeHelp import BridgeHelp
    bridgeHelp = BridgeHelp()

    class MyParser2(argparse.ArgumentParser):
        def error(self, detail=None): #values=None,choice=None): #,settings = None):
             
            if bridgeHelp.evaluate_error(detail, sys.argv): 
                sys.stderr.write(detail+'\n')  
                self.print_help() 
            sys.exit() 
            
            dt,dc = self.evaluate_detail(detail,sys.argv) 
            if dt == 'SPLASH': 
                bridge_draw_image() 
                bridge_offer_tutorial() 
                self.print_help()
            elif dt == 'NO_COMMAND': 
                self.print_help() 
            elif dt == 'TUTORIAL': 
                bridge_tutorial(dc) 
            elif dt == 'MANUAL': 
                bridge_manual(dc) 
            elif dt == 'REQUIREMENTS': 
                bridge_requirements(dc) 
            else:
                sys.stderr.write(detail+'\n')  
                self.print_help() 
            sys.exit(2)
    
        def evaluate_detail(self,X,lp):
            if len(lp) == 1: return 'SPLASH','NA'
            if X in ['too few arguments','the following arguments are required: cmd']: return 'NO_COMMAND',lp[1] 
            if lp[1].upper() in ['REQUIREMENTS','MANUAL','TUTORIAL']: return lp[1].upper(), 'FULL' 
            return 'NA','NA'  


        def dir_path(self,path):
            if not os.path.isdir(path): os.makedirs(path) 
            if not os.path.isdir(path+'/logs'): os.makedirs(path+'/logs') 
            return path 
    
        
        def find_platform(self, p): 
            if p == None: P = platform.system().upper() 
            else:         P = p.upper() 
            if P[0:3] == 'MAC':                     return 'mac' 
            elif P[0:3] in ['LIN','UNI','SER']:     return 'linux' 
            else:                                   self.error('Unrecognized --platform '+p+' (Linux and Mac are supported)') 
            
        def force_caps(self,x): 
            return x.upper() 
