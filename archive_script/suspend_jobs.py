import subprocess
import sys

def suspend_jobs():
    # run squeue to get all jobs
    # ensure full job names and job IDs is printed
    output = subprocess.check_output('squeue -u $USER -o "%.18i %.9P %.100j %.8u %.2t %.10M %.6D %R"', shell=True).decode('utf-8')
    # split the output into lines
    lines = output.split('\n')
    num_lines = len(lines)
    # get the job IDs and names
    job_ids = [line.split()[0] for line in lines[1:num_lines-1]]
    job_names = [line.split()[2] for line in lines[1:num_lines-1]]
    job_states = [line.split()[4] for line in lines[1:num_lines-1]]
    # suspend all jobs found
    for job_id, job_name, job_state in zip(job_ids, job_names, job_states):
       subprocess.call('scontrol suspend ' + job_id, shell=True)

def resume_jobs():
    # run squeue to get all jobs that are suspended
    # ensure full job names and job IDs is printed
    output = subprocess.check_output('squeue -u $USER -o "%.18i %.9P %.100j %.8u %.2t %.10M %.6D %R"', shell=True).decode('utf-8')
    # split the output into lines
    lines = output.split('\n')
    num_lines = len(lines)
    # get the job IDs and names
    job_ids = [line.split()[0] for line in lines[1:num_lines-1]]
    job_names = [line.split()[2] for line in lines[1:num_lines-1]]
    job_states = [line.split()[4] for line in lines[1:num_lines-1]]
    # resume all jobs found in suspended state
    for job_id, job_name, job_state in zip(job_ids, job_names, job_states):
       if job_state == 'S':
           subprocess.call('scontrol resume ' + job_id, shell=True)

# get input arguments
args = sys.argv[1:]
if len(args) == 0:
    print('Please provide the action to perform: suspend or resume')
    sys.exit()
else:
    action = args[0]
    if action == 'suspend':
        suspend_jobs()
    elif action == 'resume':
        resume_jobs()
    else:
        print('Please provide the action to perform: suspend or resume')
        sys.exit()

