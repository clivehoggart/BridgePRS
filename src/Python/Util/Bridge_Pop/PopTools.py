import sys, os, gzip, shutil
from collections import defaultdict as dd
from collections import Counter as cc 
#import BridgeTools 



#from .BridgePop      import BridgePop

# model_file
REV_COMP =  {'A': 'T', 'C': 'G','G': 'C', 'T': 'A'}
NUM_STRS =  ['1','2','3','4','5','6','7','8','9','0'] 


def compare_alleles(g1, g2):                                                                                                                                                                                            
	try:                                                                                                                                                                                                                
		if g1[0] == g2[0] and g1[1] == g2[1]: return 'MATCH'                                                                                                                                                            
		if g1[0] == g2[1] and g1[1] == g2[0]: return 'SWAPREF'                                                                                                                                                          
		if REV_COMP[g1[0]] == g2[0] and REV_COMP[g1[1]] == g2[1]: return 'REVCOMP'                                                                                                                                      
	except: pass                                                                                                                                                                                                        
	return 'INVALID'         



def warn(eString): 
    if type(eString) in [list,tuple]:  
        sys.stderr.write('BridgePopWarning: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('BridgePopWarning: '+eString+'\n')
    return 


def flat_warn(eString): 
    if type(eString) in [list,tuple]:  
        sys.stderr.write('BridgePopWarning: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('BridgePopWarning: '+eString+'\n')
    return 


def get_chr_strs(itbl):                                                                                                                                                                                           
    c_int, c_str = [], []                                                                                                                                                                                               
    for c in itbl:                                                                                                                                                                                                      
        try:               c_int.append(int(c))                                                                                                                                                                         
        except ValueError: c_str.append(c)                                                                                                                                                                              
    return [str(c) for c in sorted(c_int) + sorted(c_str)]                                                                                                                                                              
                                                                










def bridge_debug_error2(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('BridgeDebugError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('BridgeDebugError: '+eString+'\n')
    sys.exit(2) 


# Ambiguous debug_level Ref/Alt Missing snps.txt check printout

def bridge_pop_error2(eString):
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgePopError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('\nBridgePopError: '+eString+'\n')
    sys.exit(2) 


def bridge_sumstats_error2(eString): 
    if type(eString) in [list,tuple]:  
        sys.stderr.write('\nBridgeSumstatsError: '+eString[0]+'\n')
        for es in eString[1::]: sys.stderr.write('                 '+es+'\n')
    else: sys.stderr.write('\nBridgeSumstatsError: '+eString+'\n')
    sys.exit(2) 




            
def zip_open(fp, WITH_HEADER = False, ONLY_HEADER = False): 
    if fp.split('.')[-1] == 'gz': gf = gzip.open(fp, 'rt') 
    else:                         gf = open(fp, 'rt')  
    if WITH_HEADER: 
        hl = gf.readline().split() 
        return gf, hl 
    if ONLY_HEADER: 
        HL = gf.readline().split() 
        gf.close() 
        return HL 
    else: return gf 
        








def get_prefix_suffix(cands): 
    my_prefix,my_suffix, k = ' ', ' ', 0                                                                                                                                                                                
    while True:                                                                                                                                                                                                         
        my_prefixes = list(set([c[0:k] for c in cands]))                                                                                                                                                                
        if len(my_prefixes) == 1:                                                                                                                                                                                       
            my_prefix = my_prefixes[0]                                                                                                                                                                                  
            k+=1                                                                                                                                                                                                        
        else:                                                                                                                                                                                                           
            break                                                                                                                                                                                                       
    k = 0                                                                                                                                                                                                               
    while True:                                                                                                                                                                                                         
        my_suffixes = list(set([c[len(c)-k::] for c in cands]))                                                                                                                                                         
        if len(my_suffixes) == 1:                                                                                                                                                                                       
            my_suffix = my_suffixes[0]                                                                                                                                                                                  
            k+=1                                                                                                                                                                                                        
        else:                                                                                                                                                                                                           
            break                                                                                                                                                                                                       
    
    if len(my_prefix) == 0: my_prefix = ' '                                                                                                                                                                             
    if len(my_suffix) == 0: my_suffix = ' '  
    return my_prefix, my_suffix    






