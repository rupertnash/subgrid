"""Get an integer task ID number, either from a GridEngine array job
or a number given on the command line.

"""

import sys
import os

def getTaskID():
    try:
        taskID = int(os.environ['SGE_TASK_ID'])
    except KeyError:
        taskStr = sys.argv[1]
        taskID = int(taskStr)
        pass
    
    return taskID
