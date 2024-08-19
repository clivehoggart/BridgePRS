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
        self.args, self.io, self.pop_data = io.args, io, io.pop_data
                
        self.GT = False 
    


    def apply(self): 
        
        if self.io.cmd == 'check-requirements':  self.io.progress.finish('Complete',FIN=True) 
        if self.io.cmd == 'reformat-sumstats': self.reformat_sumstats() 
        sys.exit() 


    def read_genotype_prefix(self, gt_prefix): 
        g_path, g_name = "/".join(gt_prefix.split('/')[0:-1]), gt_prefix.split('/')[-1] 
        self.gt_lookup = dd(bool) 
        for f in os.listdir(g_path): 
            if f.split('.')[-1] == 'bim' and f[0:len(g_name)] == g_name: 
                with open(g_path+'/'+f) as filehandle:
                    for line in filehandle: self.gt_lookup[line.split()[1]] = line.split()[0]  
        self.GT = True 


    def reformat_sumstats(self): 
        self.io.progress.start_module(self.io.module, self.io.cmd, self.io.paths['home'])# .show_settings(self.settings) 
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







        







