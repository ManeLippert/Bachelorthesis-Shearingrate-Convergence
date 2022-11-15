import sys
import time  
import math  

def progressbar(current_value, required_value, barsize=40,
                prefix='', 
                progress_fill='=', progress_fill_top='', progress_fill_bot='',
                progress_unfill='.', 
                progress_bracket=['[',']']):
    
    x = int(barsize*current_value/required_value)
    percent = int(100*current_value/required_value)
    
    try:
        percent_format = '  (' + (int(math.log10(100)) - int(math.log10(percent)))*' ' + '{}%)'
    except ValueError:
        percent_format = '  (' + int(math.log10(100))*' ' + '{}%)'
        
    try:
        ratio_format = '  ' + (int(math.log10(required_value)) - int(math.log10(current_value)))*' ' + '{}/{}'
    except ValueError:
        ratio_format = '  ' + int(math.log10(required_value))*' ' + '{}/{}'

        
    bar_format =  '{}' + progress_bracket[0] + progress_fill_bot + '{}' + progress_fill_top + '{}' + progress_bracket[1] + percent_format + ratio_format
    bar = bar_format.format(prefix, progress_fill*x, progress_unfill*(barsize-x), percent, current_value, required_value)
        
    return bar

i = 0
timestep_req = 70000

while i <= timestep_req:
    print(progressbar(i, timestep_req, barsize=100, 
                      progress_fill='=', progress_unfill=' ', progress_bracket=['','']), flush=True, file=sys.stdout, end='\r')
    time.sleep(0.001)
    i += 10

print('\n')

