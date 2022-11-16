def get_time_in_seconds(time):
    
    # Format D-HH:MM:SS or HH:MM:SS
    
    try:
        day_split = time.split('-')
        time_split = day_split[1].split(':')

        day_sec = int(day_split[0])*24*60*60
        hour_sec = int(time_split[0])*60*60
        min_sec = int(time_split[1])*60
        sec = int(time_split[2])

        time_sec = day_sec + hour_sec + min_sec + sec
        
    except IndexError:
        time_split = time.split(':')
        
        hour_sec = int(time_split[0])*60*60
        min_sec = int(time_split[1])*60
        sec = int(time_split[2])

        time_sec = hour_sec + min_sec + sec

    return time_sec
    
print(get_time_in_seconds('12:23:49'))

print(get_time_in_seconds('1-12:23:49'))