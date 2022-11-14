import sys
import time    

def progressbar(current_value, required_value, barsize=60,
                prefix='', 
                progress_fill='=', progress_fill_top='', progress_fill_bot='',
                progress_unfill='.', 
                progress_bracket=['[',']']):
    
    x = int(barsize*current_value/required_value)
    
    percent = int(100*current_value/required_value)
    
    percent_format = ' ' + int(log10(percent))*' ' + '{}% {}/{}'
    ratio_format = ' ' + int(log10(current_value))*' ' + '{}/{}'
        
    bar_format =  '{}' + progress_bracket[0] + progress_fill_bot + '{}' + progress_fill_top + '{}' + progress_bracket[1] + percent_format
    bar = bar_format.format(prefix, progress_fill*x, progress_unfill*(barsize-x), percent, current_value, required_value)
        
    return bar

i = 0
timestep_req = 70000

while i <= timestep_req:
    print(progressbar(i, timestep_req, progress_fill='─', progress_unfill=' ', progress_bracket=['[',']']))
    time.sleep(0.1)
    i += 1000
