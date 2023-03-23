table_inner_width = 76

jobStatusHeader = ["                     NAME     USER    STATUS TASK NODES TIME (W:DD:HH:MM:SS)"]

jobStatusInfo_cols   = [table_inner_width - 59, 8, 9, 10, 5, 6, 20, 1]
jobStatusInfo_format = "".join(["{:>" + str(col) + "}" for col in jobStatusInfo_cols])
        
jobStatusInfo =  ["", "6.4/3x3", "bt712347", "STARTING", "32", "3", "0:03:15:52:16", ""]
jobStatusInfo = [jobStatusInfo_format.format(*jobStatusInfo)]

print(jobStatusHeader)
print(jobStatusInfo)